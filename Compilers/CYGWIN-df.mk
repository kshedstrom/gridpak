#
# Include file for Compaq Visual Fortran compiler on Cygwin
# -----------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#

              BIN := $(BIN).exe

               FC := df
           FFLAGS := /stand:f95
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN
            CLEAN := Bin/cpp_clean
               LD := ncargf90
          LDFLAGS := /nodefaultlib:libcmt /stack:67108864 
               AR := ar
          ARFLAGS := r
	    MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
	     TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
# These are be specified in Unix form and will be converted as
# necessary to Windows form for Windows-native commands. The default
# values below assume that Cygwin mounts have been defined pointing to
# the NETCDF and ARPACK library locations.
#

#   NETCDF_INCDIR ?= /netcdf-win32/include
#   NETCDF_LIBDIR ?= /netcdf-win32/lib

    NETCDF_INCDIR ?= 'C:\\include\\'
    NETCDF_LIBDIR ?= 'C:\\lib\\'

         CPPFLAGS += -I$(NETCDF_INCDIR)
       NETCDF_LIB := $(NETCDF_LIBDIR)/netcdfs.lib

ifdef ARPACK
    ARPACK_LIBDIR ?= /arpack-win32/lib
       ARPACK_LIB := $(ARPACK_LIBDIR)/arpack.lib
endif

#
# Compiler flags
#

ifdef DEBUG
           FFLAGS += /debug:full /traceback
else
           FFLAGS += /fast
endif

#
# For a Windows compiler, create variables pointing to the Windows
# file names needed when linking. Use of the "=" sign means that
# variables will be evaluated only when needed.
#

         BIN_WIN32 = "$$(cygpath --windows $(BIN))"
        LIBS_WIN32 = "$$(cygpath --windows $(NETCDF_LIB))"
ifdef ARPACK
        LIBS_WIN32 += "$$(cygpath --windows $(ARPACK_LIB))"
endif

#
# For a Windows compiler, override the compilation rule
#

%.o: %.f90
	$(FC) $(FFLAGS) /compile $< /object:$@


