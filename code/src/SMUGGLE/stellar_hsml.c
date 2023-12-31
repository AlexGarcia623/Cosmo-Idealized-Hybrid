/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/SMUGGLE/stellar_hsml.c
 * \date        MM/YYYY
 * \author      Paul Torrey, Alex Qi, Federico Marinacci, and Laura Sales
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"


#ifdef SMUGGLE_STELLAR_EVOLUTION

static int find_hsml_from_cells(int target, int mode, int thread_id);


#ifdef DETAILEDTIMINGS
static double tstart;
static int current_timebin;
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];
  in->Hsml = STP(StarParticle[i].index).Hsml;

  in->Firstnode = firstnode;
}

 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  int exact_NumNgb;
  MyFloat NumNgb;
  MyFloat NormSph;
  MyFloat Dhsmlrho;
} data_out;

static data_out *DataResult, *DataOut;

 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      StarParticle[i].exact_NumNgb = out->exact_NumNgb;
      StarParticle[i].NumNgb = out->NumNgb;
      StarParticle[i].NormSph = out->NormSph;
      StarParticle[i].Dhsmlrho = out->Dhsmlrho;
    }
  else                          /* combine */
    {
      StarParticle[i].exact_NumNgb += out->exact_NumNgb;
      StarParticle[i].NumNgb += out->NumNgb;
      StarParticle[i].NormSph += out->NormSph;
      StarParticle[i].Dhsmlrho += out->Dhsmlrho;
    }
}

#include "../generic_comm_helpers2.h"

static int Npart;
static MyFloat *Left, *Right;
static unsigned char *Todo;

static void kernel_local(void)
{
#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

  /* do local particles */
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Npart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Npart)
          break;

        if(Todo[i] > 0)         /* do we already have hsml for this star? */
          {
            int p = StarParticle[i].index;

            if(P[p].Ti_Current != All.Ti_Current)
              {
                terminate("we should not get here");
              }
            find_hsml_from_cells(i, MODE_LOCAL_PARTICLES, threadid);
          }
      }
  }


#ifdef DETAILEDTIMINGS
  double t1 = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
          MODE_LOCAL_PARTICLES, timediff(tstart, t0), timediff(tstart, t1));
#endif
}

static void kernel_imported(void)
{
#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

  /* now do the particles that were sent to us */
  int i, cnt = 0;

#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        find_hsml_from_cells(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

#ifdef DETAILEDTIMINGS
  double t1 = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
          MODE_IMPORTED_PARTICLES, timediff(tstart, t0), timediff(tstart, t1));
#endif
}


void smuggle_calculate_stellar_hsml(int npart)
{
  int i, npleft, iter = 0;
  long long ntot, npartall;

  Npart = npart;
  sumup_large_ints(1, &npart, &npartall);
  if(npartall == 0)
    return;

  CPU_Step[CPU_MISC] += measure_time();
  double t0 = second();

  Left = (MyFloat *) mymalloc("Left", npart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", npart * sizeof(MyFloat));
  Todo = (unsigned char *) mymalloc("Todo", npart * sizeof(unsigned char));

  for(i = 0; i < npart; i++)
    {
      Left[i] = Right[i] = 0;
      Todo[i] = 1;
    }

  generic_set_MaxNexport();

#ifdef DETAILEDTIMINGS
  tstart = second();
  current_timebin = All.HighestActiveTimeBin;
#endif

  /* we will repeat the whole thing for those stars where we didn't find enough neighbours */
  do
    {
      double tA = second();

      generic_comm_pattern(Npart, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < npart; i++)
        {
          if(Todo[i] > 0)
            {
              /* now check whether we had enough neighbours */
              if( (StarParticle[i].NumNgb < (All.DesNumNgbEnrichment - All.MaxNumNgbDeviationEnrichment)) || 
                  (StarParticle[i].NumNgb > (All.DesNumNgbEnrichment + All.MaxNumNgbDeviationEnrichment))  )
                {
                  /* need to redo this particle */
                  npleft++;

                  if(StarParticle[i].NumNgb > 0)
                    {
                      StarParticle[i].Dhsmlrho *= STP(StarParticle[i].index).Hsml / (NUMDIMS * StarParticle[i].NumNgb / (NORM_COEFF * pow(STP(StarParticle[i].index).Hsml, 3)));

                      if(StarParticle[i].Dhsmlrho > -0.9)       /* note: this would be -1 if only a single particle at zero lag is found */
                        StarParticle[i].Dhsmlrho = 1 / (1 + StarParticle[i].Dhsmlrho);
                      else
                        StarParticle[i].Dhsmlrho = 1;
                    }
                  else
                    StarParticle[i].Dhsmlrho = 1;

#ifndef GFM_EXACT_NUMNGB
                  if(Left[i] > 0 && Right[i] > 0 && StarParticle[i].NumNgb >  All.DesNumNgbEnrichment/2.0 )
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        Todo[i] = 0;    /* done */
                        continue;
                      }
#endif

                  if(StarParticle[i].NumNgb < (All.DesNumNgbEnrichment - All.MaxNumNgbDeviationEnrichment))
                    Left[i] = fmax(STP(StarParticle[i].index).Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(STP(StarParticle[i].index).Hsml < Right[i])
                            Right[i] = STP(StarParticle[i].index).Hsml;
                        }
                      else
                        Right[i] = STP(StarParticle[i].index).Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g  Right-Left=%g\n   pos=(%g|%g|%g)\n",
                             i, ThisTask, (int) P[StarParticle[i].index].ID, STP(StarParticle[i].index).Hsml, Left[i],
                             Right[i], (float) StarParticle[i].NumNgb, Right[i] - Left[i], P[StarParticle[i].index].Pos[0], P[StarParticle[i].index].Pos[1], P[StarParticle[i].index].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    STP(StarParticle[i].index).Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur");  /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          double fac = 1.26;

                          if(fabs(StarParticle[i].NumNgb - All.DesNumNgbEnrichment) < 0.5 * All.DesNumNgbEnrichment)
                            {
                              fac = 1 - (StarParticle[i].NumNgb - All.DesNumNgbEnrichment) / (NUMDIMS * StarParticle[i].NumNgb) * StarParticle[i].Dhsmlrho;

                              if(fac > 1.26)
                                fac = 1.26;
                            }

                          STP(StarParticle[i].index).Hsml *= fac;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          double fac = 1 / 1.26;

                          if(fabs(StarParticle[i].NumNgb - All.DesNumNgbEnrichment) < 0.5 * All.DesNumNgbEnrichment)
                            {
                              fac = 1 - (StarParticle[i].NumNgb - All.DesNumNgbEnrichment) / (NUMDIMS * StarParticle[i].NumNgb) * StarParticle[i].Dhsmlrho;

                              if(fac < 1 / 1.26)
                                fac = 1 / 1.26;
                            }
                          STP(StarParticle[i].index).Hsml *= fac;
                        }
                    }
                }
              else
                Todo[i] = 0;    /* done */
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      double tB = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("SMUGGLE_STELLAR_EVOLUTION: star ngb iteration %3d: need to repeat for %12lld particles. (previous iteration took %g sec)\n", iter, ntot, timediff(tA, tB));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in smuggle_calculate_stellar_hsml()\n");
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

#ifdef DETAILEDTIMINGS
  double tend = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
           MODE_FINISHED, timediff(tstart, tend), timediff(tstart, tend));
  fflush(FdDetailed);
#endif

  CPU_Step[CPU_GFM_ENRICH] += measure_time();
  double t1 = second();

  mpi_printf("SMUGGLE_STELLAR_EVOLUTION: active particles %lld, stellar density iterations took = %g sec\n", npartall, timediff(t0, t1));
}


/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int find_hsml_from_cells(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  double wk, dwk;

  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  MyDouble *pos = target_data->Pos;
  double h = target_data->Hsml;
  double h2 = h*h;

  double hinv = 1.0 / h;
#ifndef  TWODIMS
  double hinv3 = hinv * hinv * hinv;
#else
  double hinv3 = hinv * hinv / boxSize_Z;
#endif

  double h3 = 1.0 / hinv3;

  double hinv4 = hinv3 * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  double normsph = 0;
  double weighted_numngb = 0;
  double dhsmlrho = 0;
  int numngb = 0;

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[thread_id].R2list[n];
          double r = sqrt(r2);
          if(r2 > h2)
            continue;
              
          double u = r * hinv;

          if(u < 0.5)
            {
              wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
            }
          else
            {
              wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
              dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
            }

          weighted_numngb += NORM_COEFF * wk * h3;      /* 4.0/3 * PI = 4.188790204786, swallowed cells are not counted here, because they have mass_j=0 */
          dhsmlrho += (-(NUMDIMS * hinv * wk + u * dwk));
          numngb++;

#ifndef SMUGGLE_OMEGA_WEIGHT_SN
#ifndef GFM_TOPHAT_KERNEL
          normsph += SphP[j].Volume * wk;
#else
          normsph += SphP[j].Volume;
#endif
#else
          double cell_radius = get_cell_radius(j);
          normsph += 0.5 * (1.0 - sqrt(r2 / (r2 + cell_radius * cell_radius )));
#endif
        }
    }

#ifdef GFM_EXACT_NUMNGB
  weighted_numngb = numngb;
#endif


  out.exact_NumNgb = numngb;
  out.NumNgb = weighted_numngb;
  out.NormSph = normsph;
  out.Dhsmlrho = dhsmlrho;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#endif
