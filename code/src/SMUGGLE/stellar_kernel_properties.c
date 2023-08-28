/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_density.c
 * \date        MM/YYYY
 * \author     
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


#ifdef GFM_STELLAR_EVOLUTION

#if defined(GFM_SINGLE_CELL_INJECTION) && !defined(GFM_EXACT_NUMNGB)
#error "FAIL"
#endif

static int calc_kernel_properties_evaluate(int target, int mode, int thread_id);


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
  MyFloat ClosestNeighbourDistance;
#ifdef SMUGGLE_VAR_SN_EFF
  MyFloat AvgMetalNgb;
#endif
#ifdef SMUGGLE_MASS_WEIGHT_SN
  MyFloat TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyFloat TotSolidAngle;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
  MyFloat RadFeed_MinGasDist;
  MyFloat GasColumnDensity;
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
  MyFloat LocISMdens; 
  MyFloat LocISMdensH; 
  MyFloat LocISMmet;
#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyFloat FaceAreas[7];
#endif 
#endif
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      StarParticle[i].ClosestNeighbourDistance = out->ClosestNeighbourDistance;
#ifdef SMUGGLE_VAR_SN_EFF
      StarParticle[i].AvgMetalNgb = out->AvgMetalNgb;
#endif
#ifdef SMUGGLE_MASS_WEIGHT_SN
      StarParticle[i].TotNgbMass = out->TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
      StarParticle[i].TotSolidAngle = out->TotSolidAngle;
#endif
#ifdef  SMUGGLE_RADIATION_FEEDBACK
      StarParticle[i].RadFeed_MinGasDist = out->RadFeed_MinGasDist;
      StarParticle[i].GasColumnDensity = out->GasColumnDensity;
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
      StarParticle[i].LocISMdens  = out->LocISMdens;
      StarParticle[i].LocISMdensH  = out->LocISMdensH;
      StarParticle[i].LocISMmet = out->LocISMmet;
#ifdef SMUGGLE_FACE_AREA_BALANCE
      for(int k=0; k < 7; k++) { StarParticle[i].FaceAreas[k] = out->FaceAreas[k];}
#endif
#endif
    }
  else                          /* combine */
    {
      if(out->ClosestNeighbourDistance < StarParticle[i].ClosestNeighbourDistance)
        {
          StarParticle[i].ClosestNeighbourDistance = out->ClosestNeighbourDistance;
        }

#ifdef SMUGGLE_VAR_SN_EFF
      StarParticle[i].AvgMetalNgb += out->AvgMetalNgb;
#endif
#ifdef SMUGGLE_MASS_WEIGHT_SN
      StarParticle[i].TotNgbMass += out->TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
      StarParticle[i].TotSolidAngle += out->TotSolidAngle;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
      StarParticle[i].RadFeed_MinGasDist = fmin(out->RadFeed_MinGasDist, StarParticle[i].RadFeed_MinGasDist);
      StarParticle[i].GasColumnDensity += out->GasColumnDensity;
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
      StarParticle[i].LocISMdens  += out->LocISMdens;
      StarParticle[i].LocISMdensH  += out->LocISMdensH;
      StarParticle[i].LocISMmet += out->LocISMmet;

#ifdef SMUGGLE_FACE_AREA_BALANCE
      for(int k=0; k < 7; k++) { StarParticle[i].FaceAreas[k] += out->FaceAreas[k];}
#endif
#endif
    }
}

#include "../generic_comm_helpers2.h"


static int Npart;

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

        if(1)         /* do we already have hsml for this star? */
          {
            int p = StarParticle[i].index;

            if(P[p].Ti_Current != All.Ti_Current)
              {
                terminate("we should not get here");
              }

            calc_kernel_properties_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        calc_kernel_properties_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

#ifdef DETAILEDTIMINGS
  double t1 = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
          MODE_IMPORTED_PARTICLES, timediff(tstart, t0), timediff(tstart, t1));
#endif
}


void smuggle_caluclate_kernel_properties(int npart)
{

  int i;
  long long npartall;

  Npart = npart;
  sumup_large_ints(1, &npart, &npartall);
  if(npartall == 0)
    return;

  CPU_Step[CPU_MISC] += measure_time();
  double t0 = second();

  generic_set_MaxNexport();

#ifdef DETAILEDTIMINGS
  tstart = second();
  current_timebin = All.HighestActiveTimeBin;
#endif

  generic_comm_pattern(Npart, kernel_local, kernel_imported);


#if defined(SMUGGLE_STAR_FEEDBACK) || defined(SMUGGLE_RADIATION_FEEDBACK)
  for(i = 0; i < npart; i++)
    {
      if(StarParticle[i].NumNgb <= 0)
        terminate("No neighbours found within star Hsml radius");

#ifdef SMUGGLE_VAR_SN_EFF
#ifdef SMUGGLE_MASS_WEIGHT_SN
      StarParticle[i].AvgMetalNgb /= StarParticle[i].TotNgbMass;
#else
      StarParticle[i].AvgMetalNgb /= StarParticle[i].NumNgb;
#endif
      STP(StarParticle[i].index).AvgMetalNgb = StarParticle[i].AvgMetalNgb;
#endif

#ifdef SMUGGLE_MASS_WEIGHT_SN
      STP(StarParticle[i].index).TotNgbMass = StarParticle[i].TotNgbMass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
      STP(StarParticle[i].index).TotSolidAngle = StarParticle[i].TotSolidAngle;
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
      STP(StarParticle[i].index).LocISMdens = StarParticle[i].LocISMdens;
      StarParticle[i].LocISMmet /= StarParticle[i].LocISMdens;
      STP(StarParticle[i].index).LocISMmet= StarParticle[i].LocISMmet;
         
      StarParticle[i].FeedbackRadiusLimiter = get_feedback_radius_limiter(StarParticle[i].LocISMdensH, STP(StarParticle[i].index).Metallicity, STP(StarParticle[i].index).InitialMass);
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
      double hsml = fmin(STP(StarParticle[i].index).Hsml, StarParticle[i].FeedbackRadiusLimiter);
      StarParticle[i].GasColumnDensity = StarParticle[i].LocISMdens * (hsml + StarParticle[i].LocISMdens / StarParticle[i].GasColumnDensity); /* comoving units */
      StarParticle[i].GasColumnDensity /= All.cf_atime * All.cf_atime; /* physical units */
      STP(StarParticle[i].index).GasColumnDensity = StarParticle[i].GasColumnDensity;
#endif
    }
#endif

#ifdef DETAILEDTIMINGS
  double tend = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
           MODE_FINISHED, timediff(tstart, tend), timediff(tstart, tend));
  fflush(FdDetailed);
#endif

  CPU_Step[CPU_GFM_ENRICH] += measure_time();
  double t1 = second();

  mpi_printf("SMUGGLE_STELLAR_EVOLUTION: active particles %lld, stellar kernel property evaluation took = %g sec\n", npartall, timediff(t0, t1));
}


/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int calc_kernel_properties_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  double wk, dwk;

#ifdef SMUGGLE_VAR_SN_EFF
  double sum_ngb_metallicity = 0;
#endif

#ifdef SMUGGLE_MASS_WEIGHT_SN
  double tot_ngb_mass = 0;
#endif

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  double tot_solid_angle = 0;
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
#endif

#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
  double loc_ism_dens = 0;
  double loc_ism_densH = 0;
  double loc_ism_z_dens = 0;
  double s_lenght = 0;
#endif

  double minDistGas = MAX_DOUBLE_NUMBER;
  double tol = 1.0e-2;

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

  double hinv = 1.0 / h;
#ifndef  TWODIMS
  double hinv3 = hinv * hinv * hinv;
#else
  double hinv3 = hinv * hinv / boxSize_Z;
#endif

  double h3 = 1.0 / hinv3;

  double hinv4 = hinv3 * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

#ifdef SMUGGLE_FACE_AREA_BALANCE
  for(int k=0; k < 7; k++) {out.FaceAreas[k] = 0.0; } 
#endif

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[thread_id].R2list[n];
          double r = sqrt(r2);
//#ifdef SMUGGLE_STAR_FEEDBACK
//          if(P[j].Mass < 0.3 * All.TargetGasMass)   /* this is the gas particle mass ... */
//            continue;
//#endif

          if(r < h)
          {
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

//#ifdef SMUGGLE_OMEGA_WEIGHT_SN
            double cell_radius = get_cell_radius(j);
            double cell_area = cell_radius * cell_radius;
//#endif

#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION)
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
            tot_solid_angle += 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
            loc_ism_dens   += NORM_COEFF * P[j].Mass * wk;
            loc_ism_densH  += NORM_COEFF * SphP[j].MetalsFraction[element_index("Hydrogen")] * P[j].Mass * wk;
            loc_ism_z_dens += NORM_COEFF * P[j].Mass * wk * SphP[j].Metallicity;
            s_lenght += NORM_COEFF * dwk * P[j].Mass;

            /* this works because of the SMUGGLE_OMEGA_WEIGHT_SN definion above, but would break otherwise.  should be cleaned up */
 //           cell_area   = cell_radius * cell_radius;
            double omega       = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));

#ifdef SMUGGLE_FACE_AREA_BALANCE
            double dr[3];

            dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
            dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
            dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

            out.FaceAreas[0] += omega;

            if(dr[0] > 0)
                out.FaceAreas[1] += omega * dr[0] / r;
            else
                out.FaceAreas[2] -= omega * dr[0] / r;

            if(dr[1] > 0)
                out.FaceAreas[3] += omega * dr[1] / r;
            else
                out.FaceAreas[4] -= omega * dr[1] / r;

            if(dr[2] > 0)
                out.FaceAreas[5] += omega * dr[2] / r;
            else
                out.FaceAreas[6] -= omega * dr[2] / r;

#endif

#endif

#ifdef SMUGGLE_MASS_WEIGHT_SN
            tot_ngb_mass += NORM_COEFF * P[j].Mass * wk * h3;
#endif

#ifdef SMUGGLE_VAR_SN_EFF
#ifdef SMUGGLE_MASS_WEIGHT_SN
            sum_ngb_metallicity += SphP[j].Metallicity * NORM_COEFF * P[j].Mass * wk * h3;
#else
            sum_ngb_metallicity += SphP[j].Metallicity * NORM_COEFF * wk * h3;
#endif
#endif

#endif

            if(r < minDistGas)
              {
                minDistGas = r;
              }
          }
        }
    }

#ifdef SMUGGLE_VAR_SN_EFF
  out.AvgMetalNgb = sum_ngb_metallicity;
#endif
#ifdef SMUGGLE_MASS_WEIGHT_SN
  out.TotNgbMass = tot_ngb_mass;
#endif
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  out.TotSolidAngle = tot_solid_angle;
#endif
#ifdef SMUGGLE_RADIATION_FEEDBACK
  out.RadFeed_MinGasDist = minDistGas;
  out.GasColumnDensity = s_lenght;
#endif
#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_RADIATION_FEEDBACK)
  out.LocISMdens  = loc_ism_dens; 
  out.LocISMdensH = loc_ism_densH; 
  out.LocISMmet = loc_ism_z_dens; 
#endif
  out.ClosestNeighbourDistance = (1.0 + tol) * minDistGas;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#endif
