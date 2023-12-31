#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        NTYPES=6\n"
"        PERIODIC\n"
"        COOLING\n"
"        UVB_SELF_SHIELDING\n"
"        USE_SFR\n"
"        VORONOI\n"
"        NODEREFINE_BACKGROUND_GRID\n"
"        REGULARIZE_MESH_CM_DRIFT\n"
"        REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED\n"
"        REGULARIZE_MESH_FACE_ANGLE\n"
"        TREE_BASED_TIMESTEPS\n"
"        REFINEMENT_MERGE_CELLS\n"
"        ENLARGE_DYNAMIC_RANGE_IN_TIME\n"
"        NTYPES_ICS=6\n"
"        SELFGRAVITY\n"
"        EVALPOTENTIAL\n"
"        ADAPTIVE_HYDRO_SOFTENING\n"
"        OUTPUT_PRESSURE\n"
"        OUTPUTPOTENTIAL\n"
"        GFM\n"
"        GFM_STELLAR_EVOLUTION=0\n"
"        GFM_METALLICITY_IN_ICS\n"
"        GFM_COOLING_METAL\n"
"        GFM_OUTPUT_MASK=1+4+8+16+32+64\n"
"        GFM_NORMALIZED_METAL_ADVECTION\n"
"        SMUGGLE_SFR\n"
"        SMUGGLE_STAR_FEEDBACK\n"
"        SMUGGLE_STAR_FEEDBACK_TIME_LIMITER\n"
"        SMUGGLE_VARIABLE_EFFICIENCY\n"
"        SMUGGLE_OUTPUT_STELLAR_FEEDBACK\n"
"        SMUGGLE_COMPUTE_SFR_FROM_H2\n"
"        SMUGGLE_OUTPUT_VIRIAL_PARAM\n"
"        SMUGGLE_RADIATION_FEEDBACK\n"
"        SMUGGLE_RADIATION_FEEDBACK_DEBUG\n"
"        SMUGGLE_OMEGA_WEIGHT_SN\n"
"        SMUGGLE_MOLEC_COOLING\n"
"        SMUGGLE_COSMIC_RAY_HEATING\n"
"        SMUGGLE_PHOTOELECTRIC_HEATING\n"
"        SMUGGLE_SN_COOLING_RADIUS_BOOST\n"
"        SMUGGLE_DISCRETE_SN\n"
"        SMUGGLE_AGB_WINDS\n"
"        SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION=0\n"
"        SMUGGLE_FACE_AREA_BALANCE\n"
"        SMUGGLE_STELLAR_EVOLUTION\n"
"        NSOFTTYPES=6\n"
"        MULTIPLE_NODE_SOFTENING\n"
"        INDIVIDUAL_GRAVITY_SOFTENING=4+8+16+32\n"
"        RCUT=5.5\n"
"        CHUNKING\n"
"        DOUBLEPRECISION=1\n"
"        DOUBLEPRECISION_FFTW\n"
"        NGB_TREE_DOUBLEPRECISION\n"
"        COMBINETYPES\n"
"        PROCESS_TIMES_OF_OUTPUTLIST\n"
"        VORONOI_DYNAMIC_UPDATE\n"
"        NO_MPI_IN_PLACE\n"
"        NO_ISEND_IRECV_IN_DOMAIN\n"
"        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG\n"
"        LONGIDS\n"
"        OFFSET_FOR_NON_CONTIGUOUS_IDS\n"
"        SPLIT_PARTICLE_TYPE=2+4+8\n"
"        HAVE_HDF5\n"
"        OUTPUT_SOFTENINGS\n"
"\n");
}
