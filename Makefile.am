.PHONY: temponest temponest-install test
ACLOCAL_AMFLAGS = -I autoconf

# Source files used by plugins that are not part of temponest
xtra_dist_sources=triplot.py

# Specify to include the T2runtime in distributions
EXTRA_DIST = $(xtra_dist_sources)

plugindir=@TEMPO2_PLUG_INSTALL@

# Optimization flags
AM_CXXFLAGS = -O3 -g -march=native -mtune=native -ffast-math $(OPENMP_CXXFLAGS)
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

common_sources = TempoNestTextOutput.cpp pulsar_utils.cpp\
                 TempoNestPlugin.cpp \
				 utils.cpp logger.cpp
                 

common_sources += model_elements/timing_model.cpp model_elements/efac.cpp model_elements/equad.cpp model_elements/power_law_red_noise.cpp \
					model_elements/power_law_dm_noise.cpp
common_sources += types/parameter.cpp
common_sources += namespaces/settings.cpp
common_sources += json/json_loader.cpp json/json_node.cpp
common_sources += model/model_space.cpp
common_sources += likelihoods/temponest_v1.cpp
common_sources += likelihoods/gpu_functions.cpp
common_sources += samplers/sampler.cpp samplers/multinest.cpp

# Tests 
common_sources += tests/tests.cpp tests/test_likelihood.cpp

temponest_@T2ARCH@_plug_la_SOURCES = $(common_sources)
temponest_@T2ARCH@_plug_la_CXXFLAGS = $(AM_CXXFLAGS)
temponest_@T2ARCH@_plug_la_CFLAGS = $(AM_CFLAGS)

AM_CXXFLAGS += $(MPI_CXXFLAGS)
AM_CFLAGS += $(MPI_CXXFLAGS)

temponest_@T2ARCH@_plug_la_LIBADD = -lgfortran -ltempo2 -lsofa @MULTINEST_LIBS@ @GSL_LIBS@ @BLAS_LIBS@ @LIBS_LAPACK@ @TEMPO2_LIBS@ $(MPI_CXXLIBS)

temponest_@T2ARCH@_plug_la_CPPFLAGS = @MULTINEST_CFLAGS@ @GSL_CFLAGS@ @TEMPO2_CFLAGS@ @EIGEN_CPPFLAGS@ $(MPI_CXXFLAGS)

temponest_@T2ARCH@_plug_la_LDFLAGS = -fopenmp  -shared -module -avoid-version -shrext .t2 $(MPI_CXXLDFLAGS) 

# ArrayFire support
if HAVE_ARRAYFIRE
temponest_@T2ARCH@_plug_la_CPPFLAGS += @ARRAYFIRE_CPPFLAGS@ -DHAVE_ARRAYFIRE
temponest_@T2ARCH@_plug_la_LIBADD += @ARRAYFIRE_LIBS@
temponest_@T2ARCH@_plug_la_LDFLAGS += @ARRAYFIRE_LDFLAGS@
endif

temponest: temponest_@T2ARCH@_plug.la

temponest-install: temponest_@T2ARCH@_plug.la
	$(INSTALL) -d $(plugindir)
	$(INSTALL) .libs/temponest_@T2ARCH@_plug.t2 $(plugindir)

endif  # HAVE_MULTINEST

# Add this to define a custom clean-local target
clean-local:
	-find . -type d -name .deps -exec rm -rf {} +
	-find . -type f -name .dirstamp -delete


run-tests: temponest
	@echo "Running TempoNest tests..."
	@if [ -z "$$TEMPO2" ]; then \
		echo "ERROR: TEMPO2 environment variable is not set"; \
		exit 1; \
	fi
	@if [ ! -x "$$TEMPO2/bin/tempo2" ]; then \
		echo "ERROR: $$TEMPO2/bin/tempo2 does not exist or is not executable"; \
		exit 1; \
	fi
	@$$TEMPO2/bin/tempo2 -gr temponest -f tests/test_data/test.par tests/test_data/test.tim -cfile tests/test_data/test.json || (echo "Tests failed"; exit 1)

test: temponest temponest-install run-tests
