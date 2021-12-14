project('eonclient', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++17'])

host_system = host_machine.system()
eon_extra_args = []

add_languages('c', required: true)

cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

# Add conditionals
if host_system == 'darwin'
  eon_extra_args += ['-DOSX=TRUE']
  # Workaround for
  # ../MinModeSaddleSearch.cpp:173:25: error: aligned allocation function of
  # type 'void *(unsigned long, enum std::align_val_t)' is only available on macOS 10.13 or newer
  add_global_arguments(['-faligned-allocation'], language: 'cpp')
endif

# Dependencies
eigen_dep = dependency('Eigen3', method : 'cmake', modules : ['Eigen3::Eigen'], required : true)
eon_deps = [ eigen_dep ]

# Optional, only for AMS
absl_dep = dependency('absl', method : 'cmake',
  modules : [
    'absl::base',
    'absl::strings',
    'absl::numeric',
    'absl::utility',
  ],
  required : false,
)

boost_ams_dep = dependency('boost',
  modules : [
    'headers',
    'iostreams',
  ],
  required : false,
)

fmt_dep = dependency('fmt', method : 'cmake',  modules : ['fmt::fmt'],  required : false)

ams_dep = [ boost_ams_dep, absl_dep, fmt_dep  ]

# Potentials
subdir('potentials/EAM')
subdir('potentials/EMT')
subdir('potentials/ExtPot')
subdir('potentials/IMD')
subdir('potentials/LJ')
subdir('potentials/LJCluster')
subdir('potentials/Morse')
subdir('potentials/NewPot')

potentials = [
               eam,
               emt,
               extpot,
               imd,
               lennard_jones,
               lennard_jones_cluster,
               morse,
               newpot
             ]

eon_pot = declare_dependency(link_with : potentials)
  # include_directories : my_inc)

eonclient_sources = [
           'Optimizer.cpp',
           'PrefactorJob.cpp',
           'LBFGS.cpp',
           'Parameters.cpp',
           'ReplicaExchangeJob.cpp',
           'LowestEigenmode.cpp',
           'BondBoost.cpp',
           'Job.cpp',
           'INIFile.cpp',
           'GlobalOptimization.cpp',
           'Log.cpp',
           'MinModeSaddleSearch.cpp',
           'StructureComparisonJob.cpp',
           'SteepestDescent.cpp',
           'ImprovedDimer.cpp',
           'PointJob.cpp',
           'Prefactor.cpp',
           'ConjugateGradients.cpp',
           'Matter.cpp',
           'FiniteDifferenceJob.cpp',
           'Lanczos.cpp',
           'HessianJob.cpp',
           'TADJob.cpp',
           'ProcessSearchJob.cpp',
           'Bundling.cpp',
           'NudgedElasticBand.cpp',
           'MonteCarloJob.cpp',
           'DynamicsJob.cpp',
           'MonteCarlo.cpp',
           'Hessian.cpp',
           'NudgedElasticBandJob.cpp',
           'DynamicsSaddleSearch.cpp',
           'HelperFunctions.cpp',
           'Dimer.cpp',
           'Dynamics.cpp',
           'GlobalOptimizationJob.cpp',
           'BiasedGradientSquaredDescent.cpp',
           'SafeHyperJob.cpp',
           'MinimizationJob.cpp',
           'Quickmin.cpp',
           'ParallelReplicaJob.cpp',
           'Potential.cpp',
           'ClientEON.cpp',
           'BasinHoppingJob.cpp',
           'CommandLine.cpp',
           'FIRE.cpp',
           'EpiCenters.cpp',
           'SaddleSearchJob.cpp',
           'BasinHoppingSaddleSearch.cpp',
]

# ------------------- Options

if get_option('with_gprd')
  libgprd_proj = subproject('gprdimer')
  gprd_dep = libgprd_proj.get_variable('libgprd_dep')
  subdir('potentials/GPRPotential')
  potentials += [ gprpot ]
  eon_extra_args += ['-DWITH_GPRD']
  eonclient_sources += [ 'AtomicGPDimer.cpp', 'GPRHelpers.cpp' ]
  eon_deps += [ gprd_dep ]
endif

if get_option('with_vasp')
  if host_system != 'windows'
    subdir('potentials/VASP')
    potentials += [ vasp ]
    eon_extra_args += ['-DWITH_VASP']
  endif
endif

if get_option('with_water')
  subdir('potentials/Water')
  subdir('potentials/Water_H')
  subdir('potentials/Water_Pt')
  potentials += [ water, water_h, water_pt ]
  eon_extra_args += ['-DWITH_WATER']
endif

if get_option('with_ams')
  subdir('potentials/AMS')
  subdir('potentials/AMS_IO')
  potentials += [ams, ams_io]
  eon_extra_args += ['-DWITH_AMS']
endif

if get_option('with_lammps')
  subdir('potentials/LAMMPS')
  potentials += [ lammps_pot ]
  eon_extra_args += ['-DLAMMPS_POT']
endif

if get_option('with_python')
  subdir('potentials/PyAMFF')
  subdir('potentials/QSC')
  potentials += [ pyamff, qsc ]
  eon_extra_args += ['-DWITH_PYTHON']
endif

if get_option('with_mpi')
  subdir('potentials/MPIPot')
  potentials += [ mpipot ]
  eon_extra_args += ['-DEONMPI']
  mpi_dep = dependency('mpi')
  eon_deps += [ mpi_dep ]
endif

if get_option('with_fortran')
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
  subdir('potentials/Aluminum')
  subdir('potentials/EDIP')
  subdir('potentials/FeHe')
  subdir('potentials/Lenosky')
  subdir('potentials/SW')
  subdir('potentials/Tersoff')
  subdir('potentials/Water_H')
  potentials += [aluminum, edip, fehe, lenosky, sw, tersoff, water_h]
  eon_extra_args += ['-DWITH_FORTRAN']
endif

# ---------------------- Executable

ec = executable('eonclient',
                eonclient_sources,
                dependencies : eon_deps,
                link_with : potentials,
                cpp_args : eon_extra_args,
                install : true)