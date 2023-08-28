#include <stdio.h>
#include "arepoconfig.h"
#ifdef HAVE_HDF5
#include <hdf5.h>
#include "hdf5_util.h"

void write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES");
my_H5Aclose(hdf5_attribute, "NTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PERIODIC", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PERIODIC");
my_H5Aclose(hdf5_attribute, "PERIODIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "COOLING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "COOLING");
my_H5Aclose(hdf5_attribute, "COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "UVB_SELF_SHIELDING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "UVB_SELF_SHIELDING");
my_H5Aclose(hdf5_attribute, "UVB_SELF_SHIELDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "USE_SFR", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "USE_SFR");
my_H5Aclose(hdf5_attribute, "USE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI");
my_H5Aclose(hdf5_attribute, "VORONOI");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NODEREFINE_BACKGROUND_GRID", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NODEREFINE_BACKGROUND_GRID");
my_H5Aclose(hdf5_attribute, "NODEREFINE_BACKGROUND_GRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_FACE_ANGLE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_FACE_ANGLE");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_FACE_ANGLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "TREE_BASED_TIMESTEPS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "TREE_BASED_TIMESTEPS");
my_H5Aclose(hdf5_attribute, "TREE_BASED_TIMESTEPS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_MERGE_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_MERGE_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_MERGE_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENLARGE_DYNAMIC_RANGE_IN_TIME", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Aclose(hdf5_attribute, "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NTYPES_ICS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES_ICS");
my_H5Aclose(hdf5_attribute, "NTYPES_ICS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVALPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVALPOTENTIAL");
my_H5Aclose(hdf5_attribute, "EVALPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ADAPTIVE_HYDRO_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADAPTIVE_HYDRO_SOFTENING");
my_H5Aclose(hdf5_attribute, "ADAPTIVE_HYDRO_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_PRESSURE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_PRESSURE");
my_H5Aclose(hdf5_attribute, "OUTPUT_PRESSURE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUTPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUTPOTENTIAL");
my_H5Aclose(hdf5_attribute, "OUTPUTPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM");
my_H5Aclose(hdf5_attribute, "GFM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_STELLAR_EVOLUTION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_STELLAR_EVOLUTION");
my_H5Aclose(hdf5_attribute, "GFM_STELLAR_EVOLUTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_METALLICITY_IN_ICS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_METALLICITY_IN_ICS");
my_H5Aclose(hdf5_attribute, "GFM_METALLICITY_IN_ICS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_COOLING_METAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_COOLING_METAL");
my_H5Aclose(hdf5_attribute, "GFM_COOLING_METAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_MASK", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+4+8+16+32+64;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_OUTPUT_MASK");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_MASK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_NORMALIZED_METAL_ADVECTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Aclose(hdf5_attribute, "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_SFR", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_SFR");
my_H5Aclose(hdf5_attribute, "SMUGGLE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STAR_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_STAR_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STAR_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_VARIABLE_EFFICIENCY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_VARIABLE_EFFICIENCY");
my_H5Aclose(hdf5_attribute, "SMUGGLE_VARIABLE_EFFICIENCY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OUTPUT_STELLAR_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OUTPUT_STELLAR_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OUTPUT_STELLAR_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_COMPUTE_SFR_FROM_H2", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_COMPUTE_SFR_FROM_H2");
my_H5Aclose(hdf5_attribute, "SMUGGLE_COMPUTE_SFR_FROM_H2");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OUTPUT_VIRIAL_PARAM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OUTPUT_VIRIAL_PARAM");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OUTPUT_VIRIAL_PARAM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_RADIATION_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_RADIATION_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_RADIATION_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_RADIATION_FEEDBACK_DEBUG", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_RADIATION_FEEDBACK_DEBUG");
my_H5Aclose(hdf5_attribute, "SMUGGLE_RADIATION_FEEDBACK_DEBUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OMEGA_WEIGHT_SN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OMEGA_WEIGHT_SN");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OMEGA_WEIGHT_SN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_MOLEC_COOLING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_MOLEC_COOLING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_MOLEC_COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_COSMIC_RAY_HEATING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_COSMIC_RAY_HEATING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_COSMIC_RAY_HEATING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_PHOTOELECTRIC_HEATING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_PHOTOELECTRIC_HEATING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_PHOTOELECTRIC_HEATING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_SN_COOLING_RADIUS_BOOST", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_SN_COOLING_RADIUS_BOOST");
my_H5Aclose(hdf5_attribute, "SMUGGLE_SN_COOLING_RADIUS_BOOST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_DISCRETE_SN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_DISCRETE_SN");
my_H5Aclose(hdf5_attribute, "SMUGGLE_DISCRETE_SN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_AGB_WINDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_AGB_WINDS");
my_H5Aclose(hdf5_attribute, "SMUGGLE_AGB_WINDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_FACE_AREA_BALANCE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_FACE_AREA_BALANCE");
my_H5Aclose(hdf5_attribute, "SMUGGLE_FACE_AREA_BALANCE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STELLAR_EVOLUTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_STELLAR_EVOLUTION");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STELLAR_EVOLUTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NSOFTTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NSOFTTYPES");
my_H5Aclose(hdf5_attribute, "NSOFTTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MULTIPLE_NODE_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MULTIPLE_NODE_SOFTENING");
my_H5Aclose(hdf5_attribute, "MULTIPLE_NODE_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "INDIVIDUAL_GRAVITY_SOFTENING", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 4+8+16+32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Aclose(hdf5_attribute, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RCUT", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 5.5;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "RCUT");
my_H5Aclose(hdf5_attribute, "RCUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CHUNKING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CHUNKING");
my_H5Aclose(hdf5_attribute, "CHUNKING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION_FFTW", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DOUBLEPRECISION_FFTW");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION_FFTW");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NGB_TREE_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NGB_TREE_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "NGB_TREE_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "COMBINETYPES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "COMBINETYPES");
my_H5Aclose(hdf5_attribute, "COMBINETYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROCESS_TIMES_OF_OUTPUTLIST", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Aclose(hdf5_attribute, "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI_DYNAMIC_UPDATE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI_DYNAMIC_UPDATE");
my_H5Aclose(hdf5_attribute, "VORONOI_DYNAMIC_UPDATE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NO_MPI_IN_PLACE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NO_MPI_IN_PLACE");
my_H5Aclose(hdf5_attribute, "NO_MPI_IN_PLACE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NO_ISEND_IRECV_IN_DOMAIN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NO_ISEND_IRECV_IN_DOMAIN");
my_H5Aclose(hdf5_attribute, "NO_ISEND_IRECV_IN_DOMAIN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG");
my_H5Aclose(hdf5_attribute, "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "LONGIDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "LONGIDS");
my_H5Aclose(hdf5_attribute, "LONGIDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OFFSET_FOR_NON_CONTIGUOUS_IDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OFFSET_FOR_NON_CONTIGUOUS_IDS");
my_H5Aclose(hdf5_attribute, "OFFSET_FOR_NON_CONTIGUOUS_IDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SPLIT_PARTICLE_TYPE", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2+4+8;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "SPLIT_PARTICLE_TYPE");
my_H5Aclose(hdf5_attribute, "SPLIT_PARTICLE_TYPE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_SOFTENINGS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_SOFTENINGS");
my_H5Aclose(hdf5_attribute, "OUTPUT_SOFTENINGS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif
