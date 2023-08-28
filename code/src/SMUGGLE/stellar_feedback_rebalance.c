/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/stellar_feedback_rebalance.c
 * \date        05/2020
 * \author      Paul Torrey, Alex Qi, Federico Marinacci, Laura Sales
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SMUGGLE_STELLAR_EVOLUTION

static int feedback_rebalance_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyFloat FaceAreas[7];
#endif
  MyFloat FeedbackRadiusLimiter;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  int iel;

  in->Pos[0]                 = P[StarParticle[i].index].Pos[0];
  in->Pos[1]                 = P[StarParticle[i].index].Pos[1];
  in->Pos[2]                 = P[StarParticle[i].index].Pos[2];
  in->Vel[0]                 = P[StarParticle[i].index].Vel[0];
  in->Vel[1]                 = P[StarParticle[i].index].Vel[1];
  in->Vel[2]                 = P[StarParticle[i].index].Vel[2];
  in->Hsml                   = STP(StarParticle[i].index).Hsml;
  in->NormSph                = StarParticle[i].NormSph;

#ifdef SMUGGLE_FACE_AREA_BALANCE
  for (int k = 0; k < 7; k++) 
    in->FaceAreas[k] = StarParticle[i].FaceAreas[k];
#endif
  in->FeedbackRadiusLimiter = StarParticle[i].FeedbackRadiusLimiter;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{  
#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyDouble TotFaceNorm;
#endif
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
#ifdef SMUGGLE_FACE_AREA_BALANCE
      StarParticle[i].TotFaceNorm = out->TotFaceNorm;
#endif
    }
  else /* merge */
    {
#ifdef SMUGGLE_FACE_AREA_BALANCE
      StarParticle[i].TotFaceNorm += out->TotFaceNorm;
#endif
    }
}

#include "../generic_comm_helpers2.h"

static int Ncount;

static void kernel_local(void)
{
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
              if(generic_polling_primary(count, Ncount))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ncount)
          break;

        feedback_rebalance_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
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

        feedback_rebalance_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void smuggle_rebalance_face_areas(void)
{
  Ncount = Nstar;

  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

  double t1 = second();
  mpi_printf("SMUGGLE_STELLAR_EVOLUTION: rebalancing stellar face areas took %g sec\n", timediff(t0, t1));
}


static int feedback_rebalance_evaluate(int target, int mode, int threadid)
{
  int j, n, iel;
  int numnodes, *firstnode;
  double h;
#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  double wk, u, r, hinv, hinv3;
#endif
  double weight_fac;
  MyDouble *pos;
  MyFloat *vel;
  MyFloat normsph;

  MyDouble rlim;
#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyDouble *FaceAreas;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
#endif

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos                    = in->Pos;
  vel                    = in->Vel;
  h                      = in->Hsml;
  normsph                = in->NormSph;

  double h2 = h * h ;
#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
  FaceAreas = in->FaceAreas;
  out.TotFaceNorm = 0.0;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
//#ifdef SMUGGLE_STAR_FEEDBACK
//          if(P[j].Mass < 0.3 * All.TargetGasMass)
//            continue;
//#endif

          rlim = in->FeedbackRadiusLimiter;
          double r2 = Thread[threadid].R2list[n];
          double r = sqrt(r2);
          //if( (r > h) || (r > rlim) )
          if(r > h)
            continue;

#ifdef SMUGGLE_OMEGA_WEIGHT_SN				/* Weight by cell opening angle */
          double cell_radius = get_cell_radius(j);
          double cell_area = cell_radius * cell_radius;
          double omega = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_area)));
          weight_fac = omega / normsph;
#elif !defined(GFM_TOPHAT_KERNEL)                       /* Weight by normal spline kernel */
          u = r * hinv;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          weight_fac = SphP[j].Volume * wk / normsph;
#else                                                   /* Top hat kernel */
          weight_fac = SphP[j].Volume / normsph;
#endif

          if(r > rlim)
            continue;

#ifdef SMUGGLE_FACE_AREA_BALANCE
          double dr[3];

          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          MyFloat x_weight_fac, y_weight_fac, z_weight_fac;

          if(dr[0] > 0)
            x_weight_fac = (omega * dr[0] / r) / FaceAreas[1];
          else
            x_weight_fac = (omega * dr[0] / r) / FaceAreas[2];

          if(dr[1] > 0)
            y_weight_fac = (omega * dr[1] / r) / FaceAreas[3];
          else
            y_weight_fac = (omega * dr[1] / r) / FaceAreas[4];

          if(dr[2] > 0)
            z_weight_fac = (omega * dr[2] / r) / FaceAreas[5];
          else
            z_weight_fac = (omega * dr[2] / r) / FaceAreas[6];

          out.TotFaceNorm += sqrt(x_weight_fac * x_weight_fac + y_weight_fac * y_weight_fac + z_weight_fac * z_weight_fac);

#endif
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

#ifdef TRACER_MC
  myfree(prob_tracers);
#endif

  return 0;
}

#endif
