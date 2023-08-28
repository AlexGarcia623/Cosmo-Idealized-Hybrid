#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code
NTYPES=6                       # number of particle types
PERIODIC
#TWODIMS
#AXISYMMETRY                    # This is for axisymmetry in cylindrical coordinates (requires TWODIMS and a stationary mesh)
#ONEDIMS
#ONEDIMS_PARALLEL
#LONG_X=10.0
#LONG_Y=2.0
#LONG_Z=10.0
#REFLECTIVE_X=1 #=2            # if set to 2, the boundary is inflow/outflow
#REFLECTIVE_Y=1 #=2
#REFLECTIVE_Z=1 #=2

COOLING
UVB_SELF_SHIELDING            # gas is self-shielded from the cosmic background based on its density
USE_SFR
#QUICK_LYALPHA                # turns dense and cold gas to stars immediately
#QUICK_LYALPHA_LATETIMEONLY   # cooling and star formation only after a certain time
#SFR_KEEP_CELLS
#GAMMA=1.4
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#PREHEATING
#NOHYDRO

#----------------------------------------MPI/Threading Hybrid
#NUM_THREADS=4                           # use OpenMP, with the given number of threads per MPI task
# IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE
#GENERIC_ASYNC                           # enables asynchronous communication scheme

#--------------------------------------- Mesh Type
#AMR
VORONOI

#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # for VORONOI_STATIC_MESH force domain decomposition if there exist non-gas particles
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
#REGULARIZE_MESH_LLOYD
#OUTPUT_MESH_FACE_ANGLE
#STICKY_POINTS_ON_REFLECTIVE_SURFACE     # if reflective boundaries are used, allows points to move only tangentially at boundary

#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global timestep
TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)
#DECOUPLE_TIMESTEPS       # allows different timebins for gravity and hydro. use only WITHOUT FORCE_EQUAL_TIMESTEPS
#MUSCL_HANCOCK           # original (now depreciated) time integration scheme, only first order
#RUNGE_KUTTA_FULL_UPDATE
#PM_TIMESTEP_BASED_ON_TYPES=2+4      # select particle types that should be considered in setting the PM timestep
#NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP  # if this is on, PM force is not included in short-range timestep criterion

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
#REFINEMENT_SPLIT_MOST_DISTANCE_NEIGHBOUR
#REFINEMENT_MERGE_PAIRS
#REFINEMENT_VOLUME_LIMIT
#REFINEMENT_HIGH_RES_GAS
#REFINEMENT_CGM
#REFINEMENT_CGM_USE_R200M
#REFINEMENT_AROUND_BH=0                    # spatial refinement scheme near BHs (0: default, 1: ignore cell shape constraints and always refine)
#DEREFINE_ONLY_DENSE_GAS
NODEREFINE_BACKGROUND_GRID
#DEREFINE_GENTLY
#OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT       # deletes the mesh structures not needed for refinement/derefinemet to lower the peak memory consumption
#REFINEMENT_AROUND_DM                      # refine around DM particles according to their softening length (useful for binary systems)
#GMC_REFINEMENT
#JEANS_DEREFINEMENT_DENSITY_THRESHOLD
#NO_TARGET_MASS_CONDITION
#DISC_REFINE_ONLY
#REFINE_ONLY_WITH_TRACER
#ROTATING_HIGHRES_REGION
#TRACK_ROTATING_HIGHRES_REGION
#RAMP_REFINE
#SNE_RAMP_REFINE
#MHD_REFINE_ON_DIVB_FACTOR
#REFINE_ABOVE_WNM_DENSITY
#BH_BASED_CGM_ZOOM
#REFINE_MCTR

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16
#AMR_REMAP

#--------------------------------------- Gravity treatment
SELFGRAVITY                   # switch on for self-gravity
#HIERARCHICAL_GRAVITY         # use hierarchical splitting of the time integration of the gravity
#CELL_CENTER_GRAVITY          # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
#NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree
GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically
#GRAVITY_TALLBOX               # special switch for making treating gravity in z-extended box, with x/y periodic, and z nonperiodic. LONG_Z may be used but must be an integer.
#ALLOW_DIRECT_SUMMATION
#DIRECT_SUMMATION_THRESHOLD=1000
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4 #N-squared fashion gravity for a small number of particles of the given type
#NO_SELFGRAVITY_TYPE=1         # exclude particle type from self-gravity (can be used with exact gravity)
#NO_GRAVITY_TYPE=1             # disable computation of gravity on particle type
#EXACT_GRAVITY_REACTION        # include reaction to other particle types when using exact gravity
#EXTERNALGRAVITY               # switch on for external potential
#EXTERNALGY=0.0
#EXTERNALDISKPOTENTIAL
#EXTERNALSHEARBOX
#EXTERNALSHEARBOX_KSRATE_RANDOM
#EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
#ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
# ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
EVALPOTENTIAL                 # computes gravitational potential
#EXTERNALSHEETY
#COMPUTE_POTENTIAL_ENERGY
#ACCRETE_ONTO_CENTRAL_POTENTIAL # Allow mass to be accreted onto the central potential (needs CENTRAL_MASS_POTENTIAL)


#--------------------------------------- Gravity softening
#NSOFTTYPES=4                  # Number of different softening values to which particle types can be mapped.
#MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
#INDIVIDUAL_GRAVITY_SOFTENING=2+4  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
#ADAPTIVE_HYDRO_SOFTENING
#NSOFTTYPES_HYDRO=64           # this is only relevant for ADAPTIVE_HYDRO_SOFTENING can can be set to override default value of 64


#--------------------------------------- TreePM Options
# PMGRID=512
#ASMTH=1.25
#RCUT=6.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2
#ONLY_PM

#FFT_COLUMN_BASED
#PM_ZOOM_OPTIMIZED

#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC                # Enables automatic ENDIAN swapping for reading ICs
# CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done


#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
#OUTPUT_COORDINATES_IN_DOUBLEPRECISION    # will always output coordinates in double precision
#NGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension


#-------------------------------------------- Things for special behaviour
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE            # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE
#IDS_OFFSET=1                   # offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES                   # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE          # optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS
#PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH          # This extends the ghost search to the full 3x3 domain instead of the principal domain
#DOUBLE_STENCIL                 # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells
#TETRA_INDEX_IN_FACE            # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge
VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
#VORONOI_MESH_KEEP_DT_AND_DTC    # keeps DTC and DT in memory, i.e. for anisotropic transport solvers
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV       # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#TIMESTEP_OUTPUT_LIMIT          # Limit timesteps to write snaps on time for output lists with huge range
#DO_NOT_CREATE_STAR_PARTICLES
#DMPIC                          # enable special image code for dark matter simulations
#ALLOWEXTRAPARAMS
#RADIATIVE_RATES                # used in non-equilibrium chemistry model
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                  # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#VEL_POWERSPEC_BOX
#ADJ_BOX_POWERSPEC              # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#DISABLE_OPTIMIZE_DOMAIN_MAPPING
#RECOMPUTE_POTENTIAL_IN_SNAPSHOT   # needed for postprocess option 18 that can be used to calculate potential values for a snapshot
#ACTIVATE_MINIMUM_OPENING_ANGLE   # this does not open tree nodes under the relative opening criterion any more if their opening angle has dropped below a minimum angle
#USE_DIRECT_IO_FOR_RESTARTS     # Try to use O_DIRECT for low-level read/write operations of restart files to circumvent the linux kernel page caching
#PERTURB_VELOCITIES             # continuously perturb velocities when running simulation
#UVB_OFF                        #No UVB
#UVB_START                      #UVB switched on after a redshift supplied in parameterfile

#CUDA                       # enables CUDA support in Arepo
#CUDA_INSTRUMENT            # This enables instrumentation support for the nvidia profiler
#USE_DSDE                   # try to use a dynamic sparse data exchange paradigm to get rid off sparse MPI_Alltoall patterns on large partitions

#HUGEPAGES                  # use huge pages for memory allocation, through hugetlbfs library
#DETAILEDTIMINGS            # creates individual timings entries for primary/secondary kernels to diagnose work-load balancing

#PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL
#BITS_PER_DIMENSION=42      # Peano-Hilbert order
#OVERRIDE_PEANOGRID_WARNING


#--------------------------------------- Output/Input options
#READ_IN_ALL_IC_FIELDS      # Read in all fields that are in the initial condition files
#UPDATE_GRADIENTS_FOR_OUTPUT
#REDUCE_FLUSH
#OUTPUT_REFBHCOUNTER
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
#OUTPUT_CPU_CSV
#OUTPUT_TASK
#OUTPUT_TIMEBIN_HYDRO
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_VOLUME
#OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS            # output particle softenings
#OUTPUTGRAVINTERACTIONS       # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                     # needed when HDF5 I/O support is desired
#HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                   #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_HE_IONIZATION_STATE    # outputs the helium ionization state in detail (HeI, HeII, HeIII)
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN
#MEASURE_DISSIPATION_RATE      # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                # output maximum mach number of a cell
#OUTPUT_TASK
#OUTPUT_ENTROPY
#OUTPUT_CSND

#--------------------------------------- Testing and Debugging options
#DEBUG                         # enables core-dumps
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#RESTART_DEBUG
#VERBOSE                       # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo
#VTUNE_INSTRUMENT
#FORCETEST=0.001               # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_TESTFORCELAW=1      # this enables a special test to measure the effective force law of the code, can be set to 1 or 2

#--------------------------------------- Static Disk Potential
#DISK_POTENTIAL
#DISK_MASS_M0=1.0
#DISK_SCALE_R0=1.0

#--------------------------------------- Static NFW Potential
#STATICNFW
#NFW_C=12
#NFW_M200=100.0
#NFW_Eps=0.01
#NFW_DARKFRACTION=0.87
#NFW_h=0.7

#--------------------------------------- Static Isothermal Sphere Potential
#STATICISO
#ISO_M200=100.0
#ISO_R200=160.0
#ISO_Eps=0.1
#ISO_FRACTION=0.9

#--------------------------------------- Static Hernquist Potential
#STATICHQ
#HQ_M200=112.233
#HQ_C=12.0
#HQ_A=10.0
#HQ_DARKFRACTION=0.9554374

#--------------------------------------- Growing Disk Potential
#GROWING_DISK_POTENTIAL

#--------------------------------------- Dark energy
#DARKENERGY # Enables Dark Energy
#TIMEDEPDE  # read w(z) from a DE file
#RESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#EXTERNALHUBBLE # reads the hubble function from the DE file
#TIMEDEPGRAV # resacles H and G according to DE model
#DARKENERGY_DEBUG # enable writing of drift/kick table

#--------------------------------------- Glass making/ 2nd-order initial conditions / Initial conditions options
#SECOND_ORDER_ICS
#LONGIDS
#OFFSET_FOR_NON_CONTIGUOUS_IDS
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=4+8
#NTYPES_ICS=6 # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)

#-------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128

#AB_TURB
#AB_TURB_DECAYING

#READ_LEGACY_ICS

#-------------------------------------- GFM - Galaxy Formation Module
GFM                                    #master switch
GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine
#GFM_STELLAR_EVOLUTION_NO_ELEMENTS
#GFM_CONST_IMF=1                        #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                     #0 for a pure power-law that depends on DM-veldisp
#GFM_PREENRICH                          #pre enrich gas at given redshift
GFM_SET_METALLICITY                    #set the metallicity of gas in solar metallicity units
#GFM_NO_METAL_ENRICHMENT                #disable metal production by SNII and AGB stars
#GFM_EXACT_NUMNGB                       #use direct neighbor count instead of kernel weighted neighbor count
#GFM_WINDS                              #decoupled ISM winds
#GFM_WINDS_VARIABLE=0                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
#GFM_WINDS_VARIABLE_HUBBLE              #add an additional H(z)^(-1/3) factor to the wind scaling, such that it scales with halo mass not halo velocity dispersion
#GFM_WINDS_HUBBLESCALING                #scale the wind energy fraction with the Hubble rate, limit the maximum to 1
#GFM_WINDS_MASSSCALING                  #scale the wind energy mass loading with halo mass (equivalent to scaling the wind energy fraction with halo virial radius)
#GFM_WIND_ENERGY_METAL_DEPENDENCE       #this can be used to decrease the wind energy for high metallicity (mimicking higher cooling losses)
#GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH  #this selects an alternative functional form for the transition, requires GFM_WIND_ENERGY_METAL_DEPENDENCE
#GFM_WINDS_STRIPPING                    #wind metal stripping
#GFM_WINDS_THERMAL                      #not only give the wind kinetic energy but also thermal energy
#GFM_WINDS_THERMAL_NEWDEF               #with this switch, the thermal energy is specified as a fraction of the total energy
#GFM_BIPOLAR_WINDS=1                    #decoupled ISM winds: bipolar winds: 0->default, 1->relative to motion of FOF group, 3->parallel to spin of star-forming gas in halo
#GFM_WINDS_LOCAL                        #energy-driven decoupled local sigma winds
#GFM_STELLAR_FEEDBACK                   #local SNIa and AGB energy and momentum feedback
#GFM_PRIMORDIAL_RATES                   #updated coefficients for primordial chemistry and cooling
GFM_COOLING_METAL                      #metal line cooling
#GFM_UVB_CORRECTIONS                    #reionization energy corrections
#GFM_AGN_RADIATION                      #cooling suppression/heating due to AGN radiation field (proximity effect)
#GFM_STELLAR_PHOTOMETRICS               #calculate stellar magnitudes for different filters based on GALAXEV/BC03
GFM_OUTPUT_MASK=1+2+4+8+16+32+64+128   #which fields to output (see io_fields.c)
#GFM_CHECKS                             #this checks the consistency of the AuxDataID/PID indices of stars and black holes every timestep
#GFM_DISCARD_ENRICHMENT_GRADIENTS       #this disables the gradient extrapolation of the passively advected metallicity scalar variables
GFM_NORMALIZED_METAL_ADVECTION         #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
#GFM_OUTPUT_BIRTH_POS                   #output BirthPos and BirthVel for all star particles
#GFM_CHEMTAGS                           #see documentation/modules_GFM_chemtags
#GFM_WINDS_SAVE_PARTTYPE=2              #save wind particles as separate particle type instead of mixed with 4 (stars)
GFM_DISCRETE_ENRICHMENT                #allow stars to enrich nearby gas from stellar evolution only above some delta mass fraction threshold
#GFM_SPLITFE                            #see documentation/modules_GFM_chemtags
#GFM_SPLITFE_ADDINAGB                   #add in the AGB iron half-half on the two iron SNIa/SNII tags such that the sum of them should be equal to the total iron
#GFM_RPROCESS                           #see documentation/modules_GFM_chemtags, must have GFM_SPLITFE toggled as well
#GFM_LAMBDA                             #output all cooling rates
#GFM_RPROCESS_CHANNELS=10               #alternative to GFM_PROCESS, use many different channels, number is number of independent r-process channels
#GFM_RPROCESS_CHANNELS_NS_KICKS         #include neutron star kicks for NSNS mergers
#GFM_RPROCESS_NSNS=7                    #the number of channels of GFM_RPROCESS_CHANNELS that are NSNS channels, the rest are SN channels
#GFM_SPROCESS                           #adds s-process elements, need to exist in yield tables
#GFM_SNIA_ENERGY_INJECTION              #add thermal energy if Ia's
GFM_NO_NEGATIVE_ELEMENT_MASS_RELEASED  #do not allow that negative yields for each element consume more mass than contained in ejecta elemental composition

#-------------------------------------- SMUGGLE - Star formation and feedback module

SMUGGLE_AGB_WINDS                          #returns momentum and energy to the ISM accounting OB and AGB stellar winds

#SMUGGLE_SFR_THRESH		# this should probably be entered back in...
SMUGGLE_STELLAR_EVOLUTION


SMUGGLE_SFR                                #turns on star formation (needs USE_SFR)
SMUGGLE_STAR_FEEDBACK                      #turns on stellar feedback
###SMUGGLE_STAR_FEEDBACK_TIME_LIMITER         #turns on time step limiter for stellar evolution
#SMUGGLE_VARIABLE_EFFICIENCY                #allows for variation of star formation efficiency based on virial parameter
#SMUGGLE_OUTPUT_SF_PROBABILITY              #enables output of the probability of transforming gas cell into star particle
#SMUGGLE_TEST_SFR                           #only calls the SF initialization and saves Kennicutt law (and the gas effective EOS if available)
#SMUGGLE_USE_POLYTROPIC_EQSTATE             #imposes a minimum temperature to star forming gas (through a polytropic equation of state)
SMUGGLE_COMPUTE_SFR_FROM_H2                #links the SFR to the H2 gas fraction
###SMUGGLE_OUTPUT_STELLAR_FEEDBACK            #outputs SNII number, feedback energy and mass released for stellar particles and log files for feedback (requires GFM_STELLAR_EVOLUTION)
SMUGGLE_OUTPUT_MOLECULAR_FRACTION          #outputs the H2 gas fraction (requires SMUGGLE_COMPUTE_SFR_FROM_H2 switched on)
#SMUGGLE_OUTPUT_OPTICAL_DEPTH               #outputs the gas optical depth (requires SMUGGLE_COMPUTE_SFR_FROM_H2 switched on)
SMUGGLE_OUTPUT_VIRIAL_PARAM                #outputs the gas cell virial parameter
#SMUGGLE_RADPRESS_OPT_THIN                  #adds radiative pressure in optically thin approximation. If GFM active only young stars are considered. Needs OTVET
#SMUGGLE_RADPRESS_OPT_THIN_LUMPERMASS       #source emits at a rate proportional to mass (IonizingLumPerSolarMass in parameterfile). Otherwise constant given by IonizingLumPerSolarMass
#SMUGGLE_RADPRESS_OPT_THICK                 #adds radiation pressure feedback using radiative transfer. Needs OTVET active
SMUGGLE_RADIATION_FEEDBACK                 #inputs momentum to gas particles within stromgren radius, keep cells at 10^4 K and prevents star formation in them.
# SMUGGLE_RADIATION_FEEDBACK_DEBUG           #extra output fields for SMUGGLE_RADIATION_FEEDBACK
#SMUGGLE_MASS_WEIGHT_SN                     #feedback energy weighted by mass instead of volume
SMUGGLE_OMEGA_WEIGHT_SN                    #feedback energy weighted by solid angle instead of volume
#SMUGGLE_VAR_SN_EFF	                    #SN efficiency scales with neighboring gas metallicity
SMUGGLE_MOLEC_COOLING                      #approx extra molecular cooling contribution addition based on fit to CLOUDY cooling curves
SMUGGLE_DUST_HEATING_COOLING               #approx extra gas-dust collisional heating cooling (Meijerink & Spaans 2005)
SMUGGLE_COSMIC_RAY_HEATING                 #approx cosmic rate heating based on formula of Guo & Oh (2008)
SMUGGLE_PHOTOELECTRIC_HEATING              #approx photoelectric heating based on formula of Wolfire (2003)
SMUGGLE_SN_COOLING_RADIUS_BOOST            #returns momentum and energy to the ISM accounting for an unresolved energy conserving early ST blast wave phase
SMUGGLE_DISCRETE_SN                        #SN feedback is done in discrete SN rather then as a continuous injection
SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION=0   #stochastic photoionization based on ionizing photon budget (0=each time step; 1=one ionization event per star particle)
#SMUGGLE_SUPERBUBBLE_LIMITER                #sets feedback coupling radius according to estimates of superbubble sizes
SMUGGLE_FACE_AREA_BALANCE
