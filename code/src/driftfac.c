/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/driftfac.c
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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/** table for the cosmological drift factors */
static double DriftTable[DRIFT_TABLE_LENGTH];

/** table for the cosmological kick factor for gravitational forces */
static double GravKickTable[DRIFT_TABLE_LENGTH];

/** table for the cosmological kick factor for hydrodynmical forces */
static double HydroKickTable[DRIFT_TABLE_LENGTH];

static double logTimeBegin;
static double logTimeMax;

/*! \brief Integrand for drift factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scalefactor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for drift factor calculation.
 */
double drift_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * a * a * a);
}

/*! \brief Integrand for gravitational kick factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scalefactor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for gravitational kick factor calculation.
 */
double gravkick_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * a * a);
}

/*! \brief Integrand for hydrodynamics kick factor calculation.
 *
 *  For cosmological simulations.
 *
 *  \param[in] a Scalefactor.
 *  \param[in] param (unused)
 *
 *  \return Integrand for hydrodynamics kick factor calculation.
 */
double hydrokick_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
}

double growthfactor_integ(double a, void *param)
{
  double s;

  s = hubble_function(a) / All.Hubble * sqrt(a * a * a);

  return pow(sqrt(a) / s, 3);
}

/*! \brief Initializes lookup table for cosmological pre-factors for a drift.
 *
 *  Numerical integrals using the integrand functions defined above.
 *
 *  \return void
 */
void init_drift_table(void)
{
#define WORKSIZE 100000
  int i;
  double result, abserr;

  /*---------- DEBUG ! ------------*/
#ifdef DARKENERGY_DEBUG
  FILE *FdDKfac;
  char buf[MAXLEN_PATH];

  file_path_sprintf(buf, "%s/driftkickfac.txt", All.OutputDir);
  FdDKfac = fopen(buf, "w");
  fprintf(FdDKfac, "i a drift GravKick HydroKick\n");
  /*---------- DEBUG ! ------------*/
#endif

  gsl_function F;
  gsl_integration_workspace *workspace;

  logTimeBegin = log(All.TimeBegin);
  logTimeMax   = log(All.TimeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      F.function = &drift_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      DriftTable[i] = result;

      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;

      F.function = &hydrokick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      HydroKickTable[i] = result;

#ifdef DARKENERGY_DEBUG
      /*---------- DEBUG ! ------------*/
      fprintf(FdDKfac, "%d %e %e %e %e \n", i, exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)),
              DriftTable[i], GravKickTable[i], HydroKickTable[i]);
      /*---------- DEBUG ! ------------*/
#endif
    }

  gsl_integration_workspace_free(workspace);

#ifdef DARKENERGY_DEBUG
  /*---------- DEBUG ! ------------*/
  fclose(FdDKfac);
  /*---------- DEBUG ! ------------*/
#endif
}

/*! \brief This function integrates the cosmological prefactor for a drift
 *         step between time0 and time1. A lookup-table is used for reasons
 *         of speed.
 *
 *  \param[in] time0 Start time.
 *  \param[in] time1 End time.
 *
 *   \return \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a)} \f].
 */
double get_drift_factor(integertime time0, integertime time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);

  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * DriftTable[0];
  else
    df2 = DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

/*! \brief This function integrates the cosmological prefactor for a
 *         gravitational kick between time0 and time1. A lookup-table is used
 *         for reasons of speed.
 *
 *  \param[in] time0 Start time.
 *  \param[in] time1 End time.
 *
 *   \return Gravkick factor.
 */
double get_gravkick_factor(integertime time0, integertime time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);

  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

/*! \brief This function integrates the cosmological prefactor for a
 *         hydrodynamical kick between time0 and time1. A lookup-table is
 *         used for reasons of speed.
 *
 *  \param[in] time0 Start time
 *  \param[in] time1 End time
 *
 *   \return Hydro kick factor.
 */
double get_hydrokick_factor(integertime time0, integertime time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * HydroKickTable[0];
  else
    df1 = HydroKickTable[i1 - 1] + (HydroKickTable[i1] - HydroKickTable[i1 - 1]) * (u1 - i1);

  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * HydroKickTable[0];
  else
    df2 = HydroKickTable[i2 - 1] + (HydroKickTable[i2] - HydroKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}
