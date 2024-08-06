.PHONY: temponest temponest-install
ACLOCAL_AMFLAGS = -I autoconf

# Source files used by plugins that are not part of temponest
xtra_dist_sources=triplot.py

# Specify to include the T2runtime in distributions
EXTRA_DIST = $(xtra_dist_sources)

# This clears out any CVS files that might get in our distribution
#dist-hook:
#	rm -rf `find $(distdir)/T2runtime -name CVS`
#	rm -rf $(distdir)/T2runtime/plugins/*
#	rm -rf $(distdir)/unsupported_plugins/built_plugins/*
#	rm -rf `find $(distdir)/unsupported_plugins/ -name CVS`
#	rm -rf `find $(distdir)/unsupported_plugins/ -name '*.o'`



plugindir=@TEMPO2_PLUG_INSTALL@

# Optimization flags
AM_CXXFLAGS = -O3 -march=native -mtune=native -ffast-math
AM_CFLAGS = $(AM_CXXFLAGS)

# Warning flags
AM_CXXFLAGS += -Wall -Wextra -Wpedantic
AM_CFLAGS += -Wall -Wextra -Wpedantic

# Debug flags (uncomment for debug builds)
# AM_CXXFLAGS += -g -DDEBUG
# AM_CFLAGS += -g -DDEBUG

if HAVE_MULTINEST

plugin_LTLIBRARIES = temponest_@T2ARCH@_plug.la

common_sources = tempo2pred_int.h T2accel.h tempo2pred.h tempo2Util.h constraints.h \
                 TKmatrix.h TKfit.h T2toolkit.h tempo2.h config.h configfile.h configfile.cpp \
                 TempoNest.h  \
                 TempoNestLikeFuncs.cpp   multinest.h \
                 TempoNestTextOutput.cpp  MultiNestParams.cpp\
                 TempoNestPlugin.cpp TempoNestParams.cpp TempoNestUtilities.cpp \
                 TNpriors.cpp 

temponest_@T2ARCH@_plug_la_SOURCES = $(common_sources)
temponest_@T2ARCH@_plug_la_CXXFLAGS = $(AM_CXXFLAGS)
temponest_@T2ARCH@_plug_la_CFLAGS = $(AM_CFLAGS)

temponest_@T2ARCH@_plug_la_LIBADD = -lgfortran -ltempo2 -lsofa @MULTINEST_LIBS@ @GSL_LIBS@ @BLAS_LIBS@ @LIBS_LAPACK@ @TEMPO2_LIBS@
temponest_@T2ARCH@_plug_la_CPPFLAGS = @MULTINEST_CFLAGS@ @GSL_CFLAGS@ @TEMPO2_CFLAGS@ @EIGEN_CPPFLAGS@

if HAVE_CULA
temponest_@T2ARCH@_plug_la_SOURCES += TempoNestGPUFuncs.cu TempoNestGPULikeFuncs.cpp
temponest_@T2ARCH@_plug_la_LIBADD += -lTNGPU @CULA_LIBS@
temponest_@T2ARCH@_plug_la_CPPFLAGS += @CULA_CFLAGS@
endif

# This sets the options that are needed to build nice modules with libtool
temponest_@T2ARCH@_plug_la_LDFLAGS = -shared -module -avoid-version -shrext .t2

temponest: temponest_@T2ARCH@_plug.la

temponest-install: temponest_@T2ARCH@_plug.la
	$(INSTALL) -d $(plugindir)
	$(INSTALL) .libs/temponest_@T2ARCH@_plug.t2 $(plugindir)

endif  # HAVE_MULTINEST

#.PHONY: plugins-install plugins plugins-clean \
#   unsupported unsupported-install unsupported-clean \
#   complete complete-install complete-clean \
#   core core-install core-clean \
#   temponest temponest-install temponest-clean
#
#install-data-hook:
#	@echo ""
#	@echo ""
#	@echo "*** NOTICE ***"
#	@echo "Plugins are no longer installed by default"
#	@echo "Use: 'make plugins-install' to install them"
#	@echo ""
#	@echo ""
