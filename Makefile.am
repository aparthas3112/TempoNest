.PHONY: temponest temponest-install
ACLOCAL_AMFLAGS = -I autoconf

# Source files used by plugins that are not part of temponest
xtra_dist_sources=triplot.py

# Specify to include the T2runtime in distributions
EXTRA_DIST = $(xtra_dist_sources)

plugindir=@TEMPO2_PLUG_INSTALL@

# Optimization flags
AM_CXXFLAGS = -O3 -march=native -mtune=native -ffast-math $(OPENMP_CXXFLAGS)
AM_CFLAGS = $(AM_CXXFLAGS)

AM_CXXFLAGS += @MPI_CXXFLAGS@
AM_CFLAGS += @MPI_CFLAGS@

# Warning flags
AM_CXXFLAGS += -Wall -Wextra -Wpedantic
AM_CFLAGS += -Wall -Wextra -Wpedantic

# Debug flags (uncomment for debug builds)
# AM_CXXFLAGS += -g -DDEBUG
# AM_CFLAGS += -g -DDEBUG

if HAVE_MULTINEST

plugin_LTLIBRARIES = temponest_@T2ARCH@_plug.la

common_sources = tempo2pred_int.h T2accel.h tempo2pred.h tempo2Util.h  \
                 TKmatrix.h TKfit.h T2toolkit.h tempo2.h config.h \
                 TempoNest.h pulsar_utils.h json_loader.h json_loader.cpp \
                 TempoNestLikeFuncs.cpp   multinest.h \
                 TempoNestTextOutput.cpp pulsar_utils.cpp\
                 TempoNestPlugin.cpp TempoNestUtilities.cpp
                 

common_sources += model_elements/timing_model.cpp model_elements/efac.cpp model_elements/equad.cpp

common_sources += types/model.cpp types/parameter.cpp
common_sources += namespaces/sampler.cpp namespaces/settings.cpp



temponest_@T2ARCH@_plug_la_SOURCES = $(common_sources)
temponest_@T2ARCH@_plug_la_CXXFLAGS = $(AM_CXXFLAGS)
temponest_@T2ARCH@_plug_la_CFLAGS = $(AM_CFLAGS)

AM_CXXFLAGS += $(MPI_CXXFLAGS)
AM_CFLAGS += $(MPI_CXXFLAGS)

temponest_@T2ARCH@_plug_la_LIBADD = -lgfortran -ltempo2 -lsofa @MULTINEST_LIBS@ @GSL_LIBS@ @BLAS_LIBS@ @LIBS_LAPACK@ @TEMPO2_LIBS@ $(MPI_CXXLIBS)

temponest_@T2ARCH@_plug_la_CPPFLAGS = @MULTINEST_CFLAGS@ @GSL_CFLAGS@ @TEMPO2_CFLAGS@ @EIGEN_CPPFLAGS@ $(MPI_CXXFLAGS)

temponest_@T2ARCH@_plug_la_LDFLAGS = -fopenmp  -shared -module -avoid-version -shrext .t2 $(MPI_CXXLDFLAGS) 

temponest: temponest_@T2ARCH@_plug.la

temponest-install: temponest_@T2ARCH@_plug.la
	$(INSTALL) -d $(plugindir)
	$(INSTALL) .libs/temponest_@T2ARCH@_plug.t2 $(plugindir)

endif  # HAVE_MULTINEST

# Add this to define a custom clean-local target
clean-local:
	-find . -type d -name .deps -exec rm -rf {} +