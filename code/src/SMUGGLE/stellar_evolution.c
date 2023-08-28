/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/stellar_evolution.c
 * \date        05/2020
 * \author      Alex Qi, Paul Torrey, Federico Marinacci, and Laura Sales
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef SMUGGLE_STELLAR_EVOLUTION

#ifndef GFM_STELLAR_EVOLUTION
#error "SMUGGLE_STELLAR_EVOLUTION requires GFM_STELLAR_EVOLUTION"
#endif

void do_smuggle_stellar_evolution()
{
  start_enrichment();           /* same routine from GFM stellar evolution to allocate and initialize StarParticle structure */

  smuggle_calculate_stellar_hsml(Nstar);        /* first non-local loop to calculate star particle HSML values */
  smuggle_caluclate_kernel_properties(Nstar);   /* second non-local loop to calculate kernel properties for each star particle */
  smuggle_rebalance_face_areas();               /* third non-local loop to rebalance face areas prior to mass/metal/momentum/energy return */

  evolve_active_stars();        /* same routine from GFM stellar evolution to determine mass and metals returned during timestep */

  find_radiation_feedback_cells();
  smuggle_execute_stellar_return();        /* fourth (and final) non-local loop to actually return the mass/metals/momentum/energy */

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  create_dust_particles();
#endif
#ifdef GFM_DUST
  dust_growth_and_destruction();
#endif

//#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
//  output_stellar_feedback_statistics();
//#endif

  /* logs only for highest time bin */
//  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
//    output_stellar_evolution_statistics();

  end_enrichment();
}

#endif

