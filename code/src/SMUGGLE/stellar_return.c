/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/stellar_return.c
 * \date        05/2020
 * \author      Paul Torrey, Alex Qi, Federico Marinacci, & Laura Sales
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

#ifdef GFM_STELLAR_EVOLUTION

#ifdef TRACER_MC
#define MAXLEN_TRACERTARGET_LIST 6
typedef struct
{
  int attach_to_task;
  int attach_to_index;
#ifdef TRACER_MC_CHECKS
  MyIDType attach_to_ID;
#endif
} tracertarget_data;
#endif

static int return_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat MetalsReleased[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat DustReleased[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  MyFloat n_SNII;
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MetalReleasedChemTags[GFM_N_CHEM_TAGS];
#endif
#ifdef GFM_RPROCESS_CHANNELS
  MyFloat MassReleasedRProcess[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  MyIDType ID;
  MyFloat SNIaEnergy;
#endif
  MyFloat TotalMetalMassReleased;
  MyFloat TotalMassReleased;
  MyFloat TotalMassReleasedSNII;
  MyFloat TotalMassReleasedSNIa;
  MyFloat TotalMassReleasedAGB;

  MyFloat n_SNII;
  MyFloat n_SNIa;

#ifdef TRACER_MC
  MyFloat Mass;
  int NumberOfTracers;
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyFloat FaceAreas[7];
#endif

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  MyFloat TotSolidAngle;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  MyFloat SNIaEnergyReleased;
  MyFloat AGBMomentumReleased;
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyFloat TotFaceNorm;
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
  MyFloat TotalEnergyReleased;
  MyFloat LocISMdensH;
  MyFloat LocISMmet;
  MyFloat FeedbackRadiusLimiter;

#ifdef SMUGGLE_AGB_WINDS
  MyFloat AGBWindSpeed;
#endif
#endif 

  MyDouble RadiationMomentumReleased;
  MyFloat StromgrenRadius;
  MyFloat RadCoolShutoffTime;
  MyFloat NormSphRadFeedback_cold;
  MyFloat StromgrenMass;

  MyFloat Lum;

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
  in->TotalMassReleased      = StarParticle[i].TotalMassReleased;	/* total mass released from all channels */
  in->TotalMassReleasedSNII  = StarParticle[i].TotalMassReleasedSNII;
  in->TotalMassReleasedSNIa  = StarParticle[i].TotalMassReleasedSNIa;
  in->TotalMassReleasedAGB   = StarParticle[i].TotalMassReleasedAGB;   

  in->n_SNII = StarParticle[i].NumSNII;
  in->n_SNIa = StarParticle[i].NumSNIa;

  in->TotalMetalMassReleased = StarParticle[i].TotalMetalMassReleased;  /* total metal mass released from all channels */

  for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
    {
      in->MetalsReleased[iel] = StarParticle[i].MetalMassReleased[iel];
    }
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          in->DustReleased[chan][iel] = StarParticle[i].DustMassReleased[chan][iel];
        }
    }
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  in->ID = P[StarParticle[i].index].ID;
  in->SNIaEnergy = 1e51 / All.UnitEnergy_in_cgs * StarParticle[i].NumSNIa;
#endif
#ifdef GFM_CHEMTAGS
  for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
      in->MetalReleasedChemTags[iel] = StarParticle[i].MetalMassReleasedChemTags[iel];
#endif
#ifdef GFM_RPROCESS_CHANNELS
  for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
      in->MassReleasedRProcess[iel] = StarParticle[i].MassReleasedRProcess[iel];
#endif
#ifdef TRACER_MC
  in->Mass            = P[StarParticle[i].index].Mass;
  in->NumberOfTracers = P[StarParticle[i].index].NumberOfTracers;
#endif

#ifdef SMUGGLE_OMEGA_WEIGHT_SN
  in->TotSolidAngle = StarParticle[i].TotSolidAngle;
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
  for(int k = 0; k < 7; k++)
    in->FaceAreas[k] = StarParticle[i].FaceAreas[k];
  in->TotFaceNorm = StarParticle[i].TotFaceNorm;
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
  in->TotalEnergyReleased = StarParticle[i].TotalEnergyReleased;
  in->LocISMdensH  = StarParticle[i].LocISMdensH;      /* local ISM H density (code units) */
  in->LocISMmet  = StarParticle[i].LocISMmet;     /* local ISM metallicity          */
  in->FeedbackRadiusLimiter = StarParticle[i].FeedbackRadiusLimiter;
#ifdef SMUGGLE_AGB_WINDS
  in->AGBWindSpeed = StarParticle[i].AGBWindSpeed;
#endif
#endif 

  in->RadiationMomentumReleased = StarParticle[i].RadiationMomentumReleased;
  in->StromgrenRadius           = StarParticle[i].StromgrenRadius;
  in->RadCoolShutoffTime        = StarParticle[i].RadCoolShutoffTime;

#ifdef SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION
  in->StromgrenMass = StarParticle[i].StromgrenMass;
  in->Lum           = StarParticle[i].Lum;
  if(StarParticle[i].GasColumnDensity > 0)
    {
      StarParticle[i].RadFeedTau =
          StarParticle[i].GasColumnDensity * All.DustOpacityRadiationFeedback * StarParticle[i].LocISMmet / GFM_SOLAR_METALLICITY;
      StarParticle[i].RadFeedTau =  fmin(StarParticle[i].RadFeedTau, 100.); // cap max TauIR
    }
  in->RadiationMomentumReleased *= (1. + StarParticle[i].RadFeedTau);
#endif

  in->NormSphRadFeedback_cold = StarParticle[i].NormSphRadFeedback_cold;


  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{  
#ifdef GFM_SNIA_ENERGY_INJECTION
  int CountSNIa;
#endif
#ifdef TRACER_MC
  int ntracertargets;
  tracertarget_data tracertargets[MAXLEN_TRACERTARGET_LIST];
#else
  char dummy;
#endif

#ifdef SMUGGLE_FACE_AREA_BALANCE
  MyDouble TotFaceNorm;
#endif

} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
#ifdef GFM_SNIA_ENERGY_INJECTION
  STP(StarParticle[i].index).NumSNIa += out->CountSNIa;
#endif
  
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
#ifdef TRACER_MC
      if(out->ntracertargets)
        {
          for(int k = 0; k < out->ntracertargets; k++)
            {
#ifdef TRACER_MC_CHECKS
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index,
                              out->tracertargets[k].attach_to_ID);
#else
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index, 0);
#endif
            }
        }
#endif
    }
  else /* merge */
    {
#ifdef TRACER_MC
      if(out->ntracertargets)
        {
          for(int k = 0; k < out->ntracertargets; k++)
            {
#ifdef TRACER_MC_CHECKS
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index,
                              out->tracertargets[k].attach_to_ID);
#else
              move_one_tracer(StarParticle[i].index, out->tracertargets[k].attach_to_task, out->tracertargets[k].attach_to_index, 0);
#endif
            }
        }
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

        return_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        return_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void smuggle_execute_stellar_return(void)
{
  Ncount = Nstar;

  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  generic_set_MaxNexport();

  double t0 = second();

#ifdef TRACER_MC
  start_MC_tracer(N_tracer); /* allocate buffer for tracer exchange */
#endif

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

#ifdef TRACER_MC
  finish_MC_tracer();
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
  for(int i = 0; i < Nstar; i++)
    {
      STP(StarParticle[i].index).FeedbackEnergy = StarParticle[i].TotalEnergyInjected / All.cf_atime / All.cf_atime;
      STP(StarParticle[i].index).FeedbackMomentum = StarParticle[i].TotalMomentumInjected / All.cf_atime;
      STP(StarParticle[i].index).FeedbackMomentumAGB = StarParticle[i].TotalMomentumInjectedAGB / All.cf_atime;
      STP(StarParticle[i].index).Cum_FeedbackEnergy += StarParticle[i].TotalEnergyInjected / All.cf_atime / All.cf_atime;
      STP(StarParticle[i].index).Cum_FeedbackMomentum += StarParticle[i].TotalMomentumReleased / All.cf_atime;
      STP(StarParticle[i].index).Cum_InjFeedbackMomentum += StarParticle[i].TotalMomentumInjected / All.cf_atime;
      STP(StarParticle[i].index).Cum_InjFeedbackMomentumAGB += StarParticle[i].TotalMomentumInjectedAGB / All.cf_atime;
    }
#endif
#endif

  double t1 = second();
  mpi_printf("SMUGGLE_STELLAR_EVOLUTION: returning mass, momentum, and energy from stars took %g sec\n", timediff(t0, t1));
}

#ifdef TRACER_MC
static int sort_probs_kernel(const void *a, const void *b)
{
  if(*((double *)a) < *((double *)b))
    return -1;

  if(*((double *)a) > *((double *)b))
    return +1;

  return 0;
}
#endif


static int return_evaluate(int target, int mode, int threadid)
{
  int j, n, iel, numnodes, *firstnode;
  double h, weight_fac, inj_energy;
  double r, r2, dr[3], FaceAreas[7], dpv[3], bf[2];

#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  double wk, u, r, hinv, hinv3;
#endif
  MyDouble *pos;
  MyFloat *vel;
  MyFloat dm_total, dm_metals, dm_metal[GFM_N_CHEM_ELEMENTS], dm_sn, dm_snii, dm_snia, dm_agb;
  MyFloat TotalMassReleased, TotalMetalMassReleased, MetalMassReleased[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat DustMassReleased[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS], dm_metal_dust[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  MyFloat n_SNII, dNumSNII;
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MetalMassReleasedChemTags[GFM_N_CHEM_TAGS], dm_metal_tags[GFM_N_CHEM_TAGS];
#endif
#ifdef GFM_RPROCESS_CHANNELS
  MyFloat MassReleasedRProcess[GFM_RPROCESS_CHANNELS], dm_metal_rprocess[GFM_RPROCESS_CHANNELS];
#endif
#ifdef GFM_SNIA_ENERGY_INJECTION
  MyFloat SNIaEnergy;
#endif
  MyFloat normsph;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  double momc[3], dpSN, msq, cdot, massr, de_feedback, Ekin, dp;

/* adding here new variable definitions that needed to be imported compared against the old routine */
  double n_SNII, n_SNIa;
  double SNIa_velocity = 0.0, SNII_velocity = 0.0;
  double du;

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

//#ifdef GFM_DISCRETE_ENRICHMENT
//  /* if this particle has no mass to distribute, quit before we actually do the neighbor search */
//  if(in->TotalMassReleased == 0.0)
//    return 0;
//#endif

  pos                    = in->Pos;
  vel                    = in->Vel;
  h                      = in->Hsml;
  normsph                = in->NormSph;
  TotalMassReleased      = in->TotalMassReleased;
  TotalMetalMassReleased = in->TotalMetalMassReleased;

  double RadiationMomentumReleased = in->RadiationMomentumReleased * All.cf_atime; /* now in comoving code units; */
  double RadCoolShutoffTime = in->RadCoolShutoffTime;                              /* time in physical code units */
  double lum = in->Lum;


#ifdef SMUGGLE_AGB_WINDS
  MyDouble TotalMassReleasedAGB, TotalEnergyReleasedAGB, TotalMomentumReleasedAGB;
  MyDouble TotEnenergyReleasedAGBPhys;

  TotalMassReleasedAGB = in->TotalMassReleasedAGB;
  if(TotalMassReleasedAGB < 0.0)
    terminate("Negative mass released by AGB stars %g\n", TotalMassReleasedAGB);
//  TotalEnergyReleasedAGB   = 0.5 * TotalMassReleasedAGB * in->AGBWindSpeed * in->AGBWindSpeed * All.cf_atime *
//                           All.cf_atime;                                             /* comoving total energy released */
  TotEnenergyReleasedAGBPhys = 0.5 * TotalMassReleasedAGB * in->AGBWindSpeed * in->AGBWindSpeed;
  TotalEnergyReleasedAGB = TotEnenergyReleasedAGBPhys * All.cf_atime * All.cf_atime; /* comoving total energy released */
  TotEnenergyReleasedAGBPhys *= (All.FeedbackEfficiency / All.one_SNe_energy); /* energy released in units of 10^51 erg */

  TotalMomentumReleasedAGB = TotalMassReleasedAGB * in->AGBWindSpeed * All.cf_atime; /* comoving total momentum released */
#endif


  /* Aggregate properties for SNII here */
  n_SNII = in->n_SNII;
  n_SNIa = in->n_SNIa;

  if(in->TotalMassReleasedSNII > 0.0)
    SNII_velocity = sqrt(2.0 * n_SNII * All.one_SNe_energy / in->TotalMassReleasedSNII); /* physical velocity in code units */
  if(in->TotalMassReleasedSNIa > 0.0)  
    SNIa_velocity = sqrt(2.0 * n_SNIa * All.one_SNe_energy / in->TotalMassReleasedSNIa); /* physical velocity in code units */

  double TotalMomentumReleased = (in->TotalMassReleasedSNII * SNII_velocity + in->TotalMassReleasedSNIa * SNIa_velocity) * All.cf_atime; /* comoving total momentum released */

  double TotalMassReleasedSN     = in->TotalMassReleasedSNII + in->TotalMassReleasedSNIa;
  double TotalSNMomentumReleased = (in->TotalMassReleasedSNII * SNII_velocity + in->TotalMassReleasedSNIa * SNIa_velocity) * All.cf_atime; /* comoving total momentum released */
  double TotalSNEnergyReleased   = (n_SNII + n_SNIa) * All.one_SNe_energy * All.cf_atime * All.cf_atime; /* comoving total energy released */


  /* Expression for terminal momentum in Hopkins+ 2017, this is the maximum boost momentum can have */
  MyFloat terminal_mom = 4.8e10 * SOLAR_MASS * All.HubbleParam / (All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s); /* physical code units  */
// federico
  /* Fire adopts the same treatment for OB/AGB winds as well (can this be done better?) */
  MyFloat terminal_mom_AGB =
      4.8e10 * SOLAR_MASS * All.HubbleParam / (All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s); /* physical code units  */

  terminal_mom *= pow((n_SNII + n_SNIa) * All.FeedbackEfficiency, 13. / 14.) * All.cf_atime;      /* comoving code units  */
  terminal_mom_AGB *= pow(TotEnenergyReleasedAGBPhys, 13. / 14.) * All.cf_atime;                  /* comoving code units  */

  MyFloat pi_kernel_mass = in->NormSphRadFeedback_cold;

  /* properties of the local ISM for cooling radius calc */
  /* to avoid zero values: facZ saturates at 2 below about 0.01 Zsun anyway */
  double z0   = fmax(in->LocISMmet / GFM_SOLAR_METALLICITY, 1e-5);
  double facZ = pow(fmin(2.0, pow(z0, -0.14)), 1.5);
  double n0   = in->LocISMdensH * All.cf_a3inv / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* physical density in cgs */
  double facn0 = pow(n0, -1. / 7.);
  double rcool = 28.4 * PARSEC * facZ * pow(n0, -3. / 7.) * (n_SNII + n_SNIa) * pow(All.FeedbackEfficiency, 2. / 7.) / All.UnitLength_in_cm * All.HubbleParam / All.cf_atime;

  terminal_mom *= (facZ * facn0);
  terminal_mom_AGB *= (facZ * facn0);

  double m_swept = 4.0/3.0 *3.14159 * rcool * rcool * rcool * in->LocISMdensH;

  double rlim = in->FeedbackRadiusLimiter;

  double momentum_norm=1.0;
#ifdef SMUGGLE_FACE_AREA_BALANCE
  momentum_norm = in->TotFaceNorm;
  for(int k = 0; k < 7; k++) FaceAreas[k] = in->FaceAreas[k];
  if(momentum_norm != 0.0){
       momentum_norm = 1.0 / momentum_norm;
  } else {
//    printf("SMUGGLE -- stellar return -- momentum_norm is zero\n");
//    printf("           The momentum_norm (set during face balance renormalization) is zero.\n");
//    printf("           Not returning mass can result in breaking mass conservation... \n");
//    printf("           normsph = %16.8e \n", normsph);
    return 0;
  }
#endif

  for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
    {
      MetalMassReleased[iel] = in->MetalsReleased[iel];
#ifdef GFM_CHEMTAGS
      MetalMassReleasedChemTags[iel] = in->MetalReleasedChemTags[iel];
#endif
    }

#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          DustMassReleased[chan][iel] = in->DustReleased[chan][iel];
        }
    }
#endif

/* this will be more complicated once feedback is included */
#if defined(GFM_DUST) || defined(DUST_LIVE)
  n_SNII = in->n_SNII; 
#endif
#ifdef GFM_CHEMTAGS
  for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
    {
      MetalMassReleasedChemTags[iel] = in->MetalReleasedChemTags[iel];
    }
#endif
#ifdef GFM_RPROCESS_CHANNELS
  for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
    {
      MassReleasedRProcess[iel] = in->MassReleasedRProcess[iel];
    }
#endif

/* this will be more complicated once feedback is included */
#ifdef GFM_SNIA_ENERGY_INJECTION
  SNIaEnergy = in->SNIaEnergy;
  out.CountSNIa = 0;
#endif

#if !defined(GFM_TOPHAT_KERNEL) && !defined(SMUGGLE_OMEGA_WEIGHT_SN)
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef TRACER_MC
  out.ntracertargets = 0;
#else
  out.dummy = 0;
#endif

#ifdef TRACER_MC
  double *prob_tracers = (double *)mymalloc("prob_tracers", in->NumberOfTracers * sizeof(double));
  for(int k = 0; k < in->NumberOfTracers; k++)
    prob_tracers[k] = get_random_number();
  double pcum1 = 0, pcum2 = 0;

  mysort(prob_tracers, in->NumberOfTracers, sizeof(double), sort_probs_kernel);

  int kpoint = 0;
#endif




  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);


  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {

#ifdef SMUGGLE_STAR_FEEDBACK
          if(P[j].Mass < 0.3 * All.TargetGasMass)
            continue;
#endif

          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);
          r2    = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
          r = sqrt(r2);

          /* If outside the kernel, this particle should not recieve anything */
          if( r > h )
            continue;


/* I think we should make this into a 'smuggle_weight_fac' function that is then called in the different places. */
          double unnormalized_weight_fac = 0.0;
#ifdef SMUGGLE_OMEGA_WEIGHT_SN
          double cell_radius = get_cell_radius(j);
          unnormalized_weight_fac = 0.5 * (1.0 - sqrt(r2 / (r2 + cell_radius*cell_radius)));
//          weight_fac = omega / normsph;
#elif !defined(GFM_TOPHAT_KERNEL)
          r2 = Thread[threadid].R2list[n];

          r = sqrt(r2);
          u = r * hinv;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          unnormalized_weight_fac = SphP[j].Volume * wk;
//          weight_fac = SphP[j].Volume * wk / normsph;
#else
          unnormalized_weight_fac = SphP[j].Volume;
//          weight_fac = SphP[j].Volume / normsph;
#endif


          MyDouble x_weight_fac, y_weight_fac, z_weight_fac;
#ifdef SMUGGLE_FACE_AREA_BALANCE

          if(dr[0] > 0)  x_weight_fac = (unnormalized_weight_fac * dr[0] / r) / FaceAreas[1];
          else           x_weight_fac = (unnormalized_weight_fac * dr[0] / r) / FaceAreas[2];

          if(dr[1] > 0)  y_weight_fac = (unnormalized_weight_fac * dr[1] / r) / FaceAreas[3];
          else           y_weight_fac = (unnormalized_weight_fac * dr[1] / r) / FaceAreas[4];

          if(dr[2] > 0)  z_weight_fac = (unnormalized_weight_fac * dr[2] / r) / FaceAreas[5];
          else           z_weight_fac = (unnormalized_weight_fac * dr[2] / r) / FaceAreas[6];
#else
          x_weight_fac = unnormalized_weight_fac / normsph * dr[0] / r;
          y_weight_fac = unnormalized_weight_fac / normsph * dr[1] / r;
          z_weight_fac = unnormalized_weight_fac / normsph * dr[2] / r;
          momentum_norm = 1.0;
#endif
          weight_fac = momentum_norm * sqrt( x_weight_fac*x_weight_fac + y_weight_fac*y_weight_fac + z_weight_fac*z_weight_fac);

          dm_total  = weight_fac * TotalMassReleased;               /* total mass returned to this cell       */
          dm_sn     = weight_fac * TotalMassReleasedSN;               /* total mass returned to this cell       */
          dm_agb    = weight_fac * TotalMassReleasedAGB;               /* total mass returned to this cell       */
          dm_metals = weight_fac * TotalMetalMassReleased;          /* total metal mass released to this cell */
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              dm_metal[iel] = weight_fac * MetalMassReleased[iel];  /* mass released for each element to this cell */
            }
#ifdef GFM_DUST
          for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
                {
                  dm_metal_dust[chan][iel] = weight_fac * DustMassReleased[chan][iel];
                }
            }
#endif

/* this will likely change when feedback is being evaluated here also */
#if defined(GFM_DUST) || defined(DUST_LIVE)
          dNumSNII = weight_fac * n_SNII;
#endif
#ifdef GFM_CHEMTAGS
          for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
            {
              dm_metal_tags[iel] = weight_fac * MetalMassReleasedChemTags[iel];
            }
#endif
#ifdef GFM_RPROCESS_CHANNELS
          for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
            {
              dm_metal_rprocess[iel] = weight_fac * MassReleasedRProcess[iel];
            }
#endif


          /* Now begin updating the SphP properties by evaluating mass/metal/momentum/energy return */
          double initial_mass = P[j].Mass;

          double initial_lab_total_sph_energy = SphP[j].Energy;
          double initial_lab_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
          double initial_lab_thermal_sph_energy = initial_lab_total_sph_energy - initial_lab_kinetic_sph_energy;


          /* boost into the stars frame */
          for(int k=0; k<3; k++) SphP[j].Momentum[k] -= P[j].Mass * vel[k];  /* boost to star's frame */

          double initial_star_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
          double initial_star_thermal_sph_energy = initial_lab_thermal_sph_energy; /* same, regardless of frame */
          double initial_star_total_sph_energy = initial_star_kinetic_sph_energy + initial_star_thermal_sph_energy;

          for(int k=0; k<3; k++)  SphP[j].Momentum[k] += P[j].Mass * vel[k];  /* boost back to lab frame */

          /* First carry out the mass return  */
          P[j].Mass += dm_total;

#if defined(REFINEMENT_HIGH_RES_GAS) || defined(REFINEMENT_CGM)
          double mass_fac = P[j].Mass / initial_mass;

#ifdef REFINEMENT_HIGH_RES_GAS
          /* mass scale factor to new total mass */
          SphP[j].HighResMass *= mass_fac;
#endif
#ifdef REFINEMENT_CGM
          SphP[j].HighResMassCGM *= mass_fac;
#endif

#endif

          /* Second, carry out the metal return (add metals; primitive variables updated below)  */
          SphP[j].MassMetallicity += dm_metals;
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              SphP[j].MassMetals[iel] += dm_metal[iel];
            }
#ifdef GFM_DUST
          for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
                {
                  SphP[j].MassMetalsDust[chan][iel] += dm_metal_dust[chan][iel];
                }
            }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
          SphP[j].NumSNII += dNumSNII;
#endif
#ifdef GFM_CHEMTAGS
          for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
            {
              SphP[j].MassMetalsChemTags[iel] += dm_metal_tags[iel];
            }
#endif
#ifdef GFM_RPROCESS_CHANNELS
          for(iel = 0; iel < GFM_RPROCESS_CHANNELS; iel++)
            {
              SphP[j].MassRProcess[iel] += dm_metal_rprocess[iel];
            }
#endif


          /* No momentum or energy is deposited beyond the feedback radius limiter */ 
          if( r > rlim )
            continue;

            Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] +
                          SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                          SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                            initial_mass;

            double Utherm = (SphP[j].Energy - Ekin) / (All.cf_atime * All.cf_atime * initial_mass);

              /* identify gas that could be photoionized, but make sure we can
                 "reselect" gas that has just returned to non-photoionized state */
            if((Utherm < 1.1 * All.PhotoionizationEgySpec) && (SphP[j].GasRadCoolShutoffTime <= 0.0))
              {
                //double alpha_rec = 2.6e-13;
                //double s = weight_fac * lum / (All.RadiationFeedbackAvgPhotonEnergyineV * ELECTRONVOLT_IN_ERGS); /* phot/sec */
                double s = lum / (All.RadiationFeedbackAvgPhotonEnergyineV * ELECTRONVOLT_IN_ERGS); /* phot/sec */
                double p_ionize = s / pi_kernel_mass;

                //double gas_dens = SphP[j].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs * All.cf_a3inv;
                //double nh       = gas_dens / PROTONMASS * SphP[j].MetalsFraction[element_index("Hydrogen")];
                //double rec = SphP[j].MetalsFraction[element_index("Hydrogen")] * nh * alpha_rec * P[j].Mass * All.UnitMass_in_g /
                //               All.HubbleParam / PROTONMASS;
                //double p_ionize = s / rec;

                if(get_random_number() < p_ionize)
                  {
                    du = fmax(All.PhotoionizationEgySpec - Utherm, 0.0);         /* physical */
                    de_feedback = All.cf_atime * All.cf_atime * du * initial_mass; /* add *comoving* change in temperature to energy */

                    SphP[j].GasRadCoolShutoffTime = RadCoolShutoffTime;
                    //out.PhotoionizationEvents += 1;
                    SphP[j].Energy += de_feedback; /* de_feedback has a momentum and heat component and in comoving code units*/
                  }
              }


          /* First input momentum if applicable */
          if(RadiationMomentumReleased > 0)
            {
              Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                            SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                            SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                initial_mass;
              dp = momentum_norm * sqrt( x_weight_fac * x_weight_fac + y_weight_fac * y_weight_fac + z_weight_fac * z_weight_fac ) * RadiationMomentumReleased;

              /* injected feedback momentum radially away */
              dpv[0] = -1.0 * momentum_norm * RadiationMomentumReleased * x_weight_fac;
              dpv[1] = -1.0 * momentum_norm * RadiationMomentumReleased * y_weight_fac;
              dpv[2] = -1.0 * momentum_norm * RadiationMomentumReleased * z_weight_fac;

              for(int k=0; k<3; k++) SphP[j].Momentum[k] += dpv[k];

              de_feedback = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                                   SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                   SphP[j].Momentum[2] * SphP[j].Momentum[2]) /
                                initial_mass - Ekin;
              // dp_tot += dp; /* LVS cumulative on each stellar part */

              SphP[j].Energy += de_feedback; /* de_feedback has a momentum and heat
                                                    component and in comoving code units*/
            }


          /* Before doing SN/AGB feedback, update the SphP 'initial' properties after radiation pressure & photoionization */ // initial_mass = P[j].Mass;

	  initial_lab_total_sph_energy = SphP[j].Energy;
          initial_lab_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / initial_mass;
          initial_lab_thermal_sph_energy = initial_lab_total_sph_energy - initial_lab_kinetic_sph_energy;

          /* boost into the stars frame */
          for(int k=0; k<3; k++) SphP[j].Momentum[k] -= initial_mass * vel[k];  /* boost to star's frame */

          initial_star_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / initial_mass;
          initial_star_thermal_sph_energy = initial_lab_thermal_sph_energy; /* same, regardless of frame */
          initial_star_total_sph_energy = initial_star_kinetic_sph_energy + initial_star_thermal_sph_energy;

          for(int k=0; k<3; k++)  SphP[j].Momentum[k] += initial_mass * vel[k];  /* boost back to lab frame */

          /* Third, carry out the momentum and energy injection and update for SN */
          if(TotalSNMomentumReleased > 0)
            {
              dpv[0] = -1.0 * momentum_norm * TotalSNMomentumReleased  * x_weight_fac; /* dp of a single cell in comoving units */
              dpv[1] = -1.0 * momentum_norm * TotalSNMomentumReleased  * y_weight_fac; /* dp of a single cell in comoving units */
              dpv[2] = -1.0 * momentum_norm * TotalSNMomentumReleased  * z_weight_fac; /* dp of a single cell in comoving units */

              /* Boost into the rest frame of the acting star for momentum injection  */
              for(int k=0; k<3; k++) momc[k] = SphP[j].Momentum[k] - vel[k] * initial_mass;     

              dpSN = dpv[0] * dpv[0] + dpv[1] * dpv[1] + dpv[2] * dpv[2];   /* Delta P squared, in the frame of the star */
              msq   = momc[0] * momc[0] + momc[1] * momc[1] + momc[2] * momc[2];  /* gas part p squared, in the star frame */
              cdot  = momc[0] * dpv[0] + momc[1] * dpv[1] + momc[2] * dpv[2];/* gas part p dotted into Delta p */
              massr = initial_mass / dm_sn;                              /* mass ratio for the returned mass */

              /* For now, assume cooling radius is unresolved and terminal momentum is injected. */
              double boost_fac = terminal_mom / TotalSNMomentumReleased;

              double sol, boost_fac_max;
              if(cdot > 0.0)
                {
                   sol = cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));
                   boost_fac_max = (1. + massr + msq / (dpSN * massr)) / sol;
                }
              else
                boost_fac_max = -cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));

              if(boost_fac > 0.9*boost_fac_max) boost_fac = 0.9*boost_fac_max;

              dpSN = sqrt(dpSN);
              for(int k=0; k<3; k++) dpv[k] *= boost_fac;
              double dp = sqrt(dpv[0] * dpv[0] + dpv[1] * dpv[1] + dpv[2] * dpv[2]);

              /* boost to star frame */
              for(int k=0; k<3; k++) SphP[j].Momentum[k] -= initial_mass * vel[k];

              /* actually inject the momentum, as calcualted in the rest frame of the star */
              for(int k=0; k<3; k++) SphP[j].Momentum[k] += dpv[k];

              double new_star_kinetic_energy = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] +                      /* new KE includes mass injection */
                                                      SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                                      SphP[j].Momentum[2] * SphP[j].Momentum[2]) / (initial_mass+dm_sn);

              /* boost back to lab/sim frame  */
              for(int k=0; k<3; k++) SphP[j].Momentum[k] += (initial_mass+dm_sn) * vel[k];

              double new_lab_kinetic_energy = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 	/* in the frame of the star */
                                                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / ( initial_mass +dm_sn ); 

              double delta_kinetic_energy = new_star_kinetic_energy - initial_star_kinetic_sph_energy;

              inj_energy =  0.5 * dpSN * dpSN / dm_sn;
              if(inj_energy > 0.0) 
		{
                  du = inj_energy * (1. - delta_kinetic_energy / inj_energy);  /* alternatively, inj_energy - delta_kinetic_energy */
		  du = 10.0 * fmax(0.0, du); //20.0 JONAH_EDIT
	        }
              else
                  terminate("SMUGGLE: SN Feedback total energy injection should never be negative");

              if(TotalSNEnergyReleased>0 && weight_fac > 0) 
                {
                  if( abs(TotalSNEnergyReleased*weight_fac - inj_energy)/inj_energy > 0.01 )
                      printf("Bad Sanity Check:  TotalSNEnergyReleased = %16.8e  weight_fac = %16.8e boost_fac = %16.8e TotalSNEnergyReleased*weight_fac = %16.8e  inj_energy = %16.8e du = %16.8e \n", 
                             TotalSNEnergyReleased, weight_fac, boost_fac, TotalSNEnergyReleased*weight_fac, inj_energy, du);
                }

              if(du < 0)
                  warn(  "negative thermal energy variation du=%g, boost %g, inj_mass %g, "
                         "inj_energy %g, dpv %g|%g|%g\n",
                              du, boost_fac, dm_sn, inj_energy, dpv[0], dpv[1], dpv[2]);

              double new_thermal_energy =  initial_lab_thermal_sph_energy + du;

              SphP[j].Energy = new_lab_kinetic_energy + new_thermal_energy;
            }



#ifdef SMUGGLE_AGB_WINDS

          /* Fourth, carry out the momentum and energy injection and update for AGBs  */
          if(TotalMomentumReleasedAGB > 0)
            {
              /* update the 'initial' SphP properties */
              initial_lab_total_sph_energy = SphP[j].Energy;
              initial_lab_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / initial_mass;
              initial_lab_thermal_sph_energy = initial_lab_total_sph_energy - initial_lab_kinetic_sph_energy;

              /* boost into the stars frame */
              for(int k=0; k<3; k++) SphP[j].Momentum[k] -= initial_mass * vel[k];  /* boost to star's frame */

              initial_star_kinetic_sph_energy =  0.5*(SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / initial_mass;
              initial_star_thermal_sph_energy = initial_lab_thermal_sph_energy; /* same, regardless of frame */
              initial_star_total_sph_energy = initial_star_kinetic_sph_energy + initial_star_thermal_sph_energy;

              for(int k=0; k<3; k++)  SphP[j].Momentum[k] += initial_mass * vel[k];  /* boost back to lab frame */


              dpv[0] = -1.0 * momentum_norm * TotalMomentumReleasedAGB  * x_weight_fac; /* dp of a single cell in comoving units */
              dpv[1] = -1.0 * momentum_norm * TotalMomentumReleasedAGB  * y_weight_fac; /* dp of a single cell in comoving units */
              dpv[2] = -1.0 * momentum_norm * TotalMomentumReleasedAGB  * z_weight_fac; /* dp of a single cell in comoving units */

              /* Boost into the rest frame of the acting star for momentum injection  */
              for(int k=0; k<3; k++) momc[k] = SphP[j].Momentum[k] - vel[k] * initial_mass; 

              dpSN  = dpv[0] * dpv[0] + dpv[1] * dpv[1] + dpv[2] * dpv[2];   /* Delta P squared, in the frame of the star */
              msq   = momc[0] * momc[0] + momc[1] * momc[1] + momc[2] * momc[2];  /* gas part p squared, in the star frame */
              cdot  = momc[0] * dpv[0] + momc[1] * dpv[1] + momc[2] * dpv[2];/* gas part p dotted into Delta p */
              massr = initial_mass / dm_agb;                              /* mass ratio for the returned mass */

              /* For now, assume boost factor is 1.0 */  // cooling radius is unresolved and terminal momentum is injected. */
              double boost_fac = fmin(sqrt(1. + massr), terminal_mom_AGB / TotalMomentumReleasedAGB);

              //double boost_fac = 1.0;	//terminal_mom / TotalSNMomentumReleased;

              double sol, boost_fac_max;
              if(cdot > 0.0)
                {
                   sol = cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));
                   boost_fac_max = (1. + massr + msq / (dpSN * massr)) / sol;
                }
              else
                boost_fac_max = -cdot / dpSN + sqrt(1. + massr + msq / (dpSN * massr) + cdot * cdot / (dpSN * dpSN));
              //if(boost_fac > 0.90*boost_fac_max) boost_fac = 0.90*boost_fac_max;

              boost_fac = fmax(fmin(boost_fac, 0.99 * boost_fac_max), 1.0);


              dpSN = sqrt(dpSN);
              for(int k=0; k<3; k++) dpv[k] *= boost_fac;
              double dp = sqrt(dpv[0] * dpv[0] + dpv[1] * dpv[1] + dpv[2] * dpv[2]);

              for(int k=0; k<3; k++) SphP[j].Momentum[k] -= initial_mass * vel[k];    /* boost to star frame */
              for(int k=0; k<3; k++) SphP[j].Momentum[k] += dpv[k];                   /* inject the momentum */

              double new_star_kinetic_energy = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                                                      SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                                      SphP[j].Momentum[2] * SphP[j].Momentum[2]) / (initial_mass + dm_agb);  /* The mass has been added now... */

              for(int k=0; k<3; k++) SphP[j].Momentum[k] += (initial_mass + dm_agb) * vel[k];    /* boost back, but with the new mass... */

              double new_lab_kinetic_energy = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
                                                     SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                                                     SphP[j].Momentum[2] * SphP[j].Momentum[2]) / (initial_mass + dm_agb);  /* The mass has been added now... */ 

              double delta_kinetic_energy = new_star_kinetic_energy - initial_star_kinetic_sph_energy;	/* change of kinetic energy in the star's frame */

              inj_energy =  0.5 * dpSN * dpSN / (TotalMassReleasedAGB * weight_fac);
              du = inj_energy - delta_kinetic_energy;  

              if(TotalEnergyReleasedAGB>0 && weight_fac > 0) 
                {
                  if( abs(TotalEnergyReleasedAGB*weight_fac - inj_energy)/inj_energy > 0.01 )
                      printf("Bad Sanity Check:  TotalAGBEnergyReleased = %16.8e  weight_fac = %16.8e boost_fac = %16.8e TotalAGBEnergyReleased*weight_fac = %16.8e  inj_energy = %16.8e du = %16.8e \n", 
                             TotalEnergyReleasedAGB, weight_fac, boost_fac, TotalEnergyReleasedAGB*weight_fac, inj_energy, du);
                }

              if(du < 0)
                {
                //  warn(  "negative thermal energy variation in AGB calculation du=%16.8e, boost %g, inj_mass %g, "
                //         "inj_energy %g, p %g|%g|%g,  dpv %g|%g|%g\n",
                //              du, boost_fac, dm_agb, inj_energy, SphP[j].Momentum[0], SphP[j].Momentum[1], SphP[j].Momentum[2], dpv[0], dpv[1], dpv[2]);
                  du = 0.0;
                }
           
              double new_thermal_energy =  initial_lab_thermal_sph_energy + du;

              SphP[j].Energy = new_lab_kinetic_energy + new_thermal_energy;

            }
#endif


 



#ifdef TRACER_MC
          if(in->Mass > 0)
            {
              double prob = dm_total / in->Mass;
              pcum2 += prob;

              while(kpoint < in->NumberOfTracers && (prob_tracers[kpoint] < pcum1))
                kpoint++;

              while(kpoint < in->NumberOfTracers && (pcum1 < prob_tracers[kpoint] && prob_tracers[kpoint] < pcum2))
                {
                  if(out.ntracertargets < MAXLEN_TRACERTARGET_LIST)
                    {
                      int n                                = out.ntracertargets++;
                      out.tracertargets[n].attach_to_task  = ThisTask;
                      out.tracertargets[n].attach_to_index = j;
#ifdef TRACER_MC_CHECKS
                      out.tracertargets[n].attach_to_ID = P[j].ID;
#endif
                    }
                  else
                    warn("reached MAXLEN_TRACERTARGET_LIST");

                  kpoint++;
                }

              pcum1 = pcum2;
            }
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
