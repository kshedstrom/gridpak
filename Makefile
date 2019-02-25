#:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Kate Hedstrom :::
#                                                                       :::
#  ROMS/TOMS Gridpak Code Master Makefile                              :::
#                                                                       :::
#  This Makefile is designed to work only with GNU Make version 3.77 or :::
#  higher. It can be used in any architecture provided that there is a  :::
#  machine/compiler rules file in the  "Compilers"  subdirectory.  You  :::
#  may need to modify the rules file to specify the  correct path  for  :::
#  the NetCDF and ARPACK libraries. The ARPACK library is only used in  :::
#  the Generalized Stability Theory analysis.                           :::
#                                                                       :::
#  If appropriate,  the USER needs to modify the  macro definitions in  :::
#  in user-defined section below.  To activate an option set the macro  :::
#  to "on". For example, if you want to compile with debugging options  :::
#  set:                                                                 :::
#                                                                       :::
#      DEBUG := on                                                      :::
#                                                                       :::
#  Otherwise, leave macro definition blank.                             :::
#                                                                       :::
#  The USER needs to provide a value for the  macro FORT.  Choose  the  :::
#  appropriate value from the list below.                               :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::


ifneq (3.80,$(firstword $(sort $(MAKE_VERSION) 3.80)))
 $(error This Makefile requires GNU make version 3.80 or higher. \
  Your current version is: $(MAKE_VERSION))
endif

#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  sources    :=
  libraries  :=
  includes   := Include
# includes   := Include_Is
# includes   := Include_Ben
# includes   := Include_S_Africa
# includes   := Include_supercritical


#==========================================================================
#  Start of user-defined options. Modify macro variables: on is TRUE while
#  blank is FALSE.
#==========================================================================
#
#  Activate debugging compiler options:

       DEBUG := on

#--------------------------------------------------------------------------
#  We are going to include a file with all the settings that depend on
#  the system and the compiler. We are going to build up the name of the
#  include file using information on both. Set your compiler here from
#  the following list:
#
#  Operating System        Compiler(s)
#
#     CYGWIN:                 g95, df
#     Darwin:                 f90
#     Linux:                  ifort, pgi, path, g95, mpif90
#
#  Feel free to send us additional rule files to include! Also, be sure
#  to check the appropriate file to make sure it has the right paths to
#  NetCDF and so on.
#--------------------------------------------------------------------------

        FORT ?= gfortran

#--------------------------------------------------------------------------
#  Set directory for executable.
#--------------------------------------------------------------------------

      BINDIR := .

#==========================================================================
#  End of user-defined options. See also the machine-dependent include
#  file being used above.
#==========================================================================

#--------------------------------------------------------------------------
#  Set directory for temporary objects.
#  #--------------------------------------------------------------------------

SCRATCH_DIR ?= Build
 clean_list := core *.ipo $(SCRATCH_DIR)

ifeq "$(strip $(SCRATCH_DIR))" "."
  clean_list := core *.o *.oo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif
ifeq "$(strip $(SCRATCH_DIR))" "./"
  clean_list := core *.o *.oo *.ipo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif

#--------------------------------------------------------------------------
#  Make functions for putting the temporary files in $(SCRATCH_DIR)
#  #  DO NOT modify this section; spaces and blank lines are needed.
#--------------------------------------------------------------------------

# $(call source-dir-to-binary-dir, directory-list)
source-dir-to-binary-dir = $(addprefix $(SCRATCH_DIR)/, $(notdir $1))

# $(call source-to-object, source-file-list)
source-to-object = $(call source-dir-to-binary-dir,      \
                   $(subst .F,.o,$1))

# $(call source-to-object, source-file-list)
c-source-to-object = $(call source-dir-to-binary-dir,      \
                     $(subst .c,.o,$(filter %.c,$1))        \
                     $(subst .cc,.o,$(filter %.cc,$1)))

$(call make-library, library-name, source-file-list)
define make-library
   libraries += $(SCRATCH_DIR)/$1
   sources   += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

$(call make-c-library, library-name, source-file-list)
define make-c-library
   libraries += $(SCRATCH_DIR)/$1
   c_sources += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .c,.o,$(filter %.c,$2)) \
                      $(subst .cc,.o,$(filter %.cc,$2)))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call f90-source, source-file-list)
f90-source = $(call source-dir-to-binary-dir,     \
                   $(subst .F,.f90,$1))

# $(compile-rules)
define compile-rules
  $(foreach f, $(local_src),       \
      $(call one-compile-rule,$(call source-to-object,$f), \
      $(call f90-source,$f),$f))
endef

# $(c-compile-rules)
define c-compile-rules
  $(foreach f, $(local_c_src),       \
    $(call one-c-compile-rule,$(call c-source-to-object,$f), $f))
endef

# $(call one-compile-rule, binary-file, f90-file, source-file)
define one-compile-rule
  $1: $2 $3
	cd $$(SCRATCH_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

# $(call one-c-compile-rule, binary-file, source-file)
define one-c-compile-rule
  $1: $2
	cd $$(SCRATCH_DIR); $$(CXX) -c $$(CXXFLAGS) $$<

endef

#--------------------------------------------------------------------------
#  Set executable file names.
#--------------------------------------------------------------------------

COAST := $(BINDIR)/coast
GRID := $(BINDIR)/grid
#TOLAT := $(BINDIR)/tolat
BATHTUB := $(BINDIR)/bathtub
BATHSUDS := $(BINDIR)/bathsuds
BATHSOAP := $(BINDIR)/bathsoap
SPHERE := $(BINDIR)/sphere
#ifdef DEBUG
#   COAST := $(BINDIR)/coastG
#   GRID := $(BINDIR)/gridG
#   TOLAT := $(BINDIR)/tolatG
#   BATHTUB := $(BINDIR)/bathtubG
#   BATHSUDS := $(BINDIR)/bathsudsG
#   BATHSOAP := $(BINDIR)/bathsoapG
#   SPHERE := $(BINDIR)/sphereG
#endif

#--------------------------------------------------------------------------
#  Set name of module files for netCDF F90 interface. On some platforms
#  these will need to be overridden in the machine-dependent include file.
#--------------------------------------------------------------------------

   NETCDF_MODFILE := netcdf.mod
TYPESIZES_MODFILE := typesizes.mod

#--------------------------------------------------------------------------
#  "uname -s" should return the OS or kernel name and "uname -m" should
#  return the CPU or hardware name. In practice the results can be pretty
#  flaky. Run the results through sed to convert "/" and " " to "-",
#  then apply platform-specific conversions.
#--------------------------------------------------------------------------

OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))

CPU := $(shell uname -m | sed 's/[\/ ]/-/g')

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
  include Compilers/$(OS)-$(strip $(FORT)).mk
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(SCRATCH_DIR) $(SCRATCH_DIR)/MakeDepend $(COAST) $(GRID) $(SQGRID) \
	$(TOLAT) $(BATHTUB) $(BATHSUDS) $(BATHSOAP) $(SPHERE)

modules  := Utility Drivers

vpath %.F $(modules)
vpath %.h $(includes)
vpath %.f90 $(SCRATCH_DIR)
vpath %.o $(SCRATCH_DIR)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes)) --silent --moddir=$(SCRATCH_DIR)

CPPFLAGS += $(patsubst %,-I%,$(includes))

$(SCRATCH_DIR):
	$(shell $(TEST) -d $(SCRATCH_DIR) || $(MKDIR) $(SCRATCH_DIR) )

#--------------------------------------------------------------------------
.PHONY: libraries

libraries: $(libraries)

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS dependecies.
#--------------------------------------------------------------------------

$(SCRATCH_DIR)/$(NETCDF_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/$(TYPESIZES_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/MakeDepend: Makefile \
                           $(SCRATCH_DIR)/$(NETCDF_MODFILE) \
                           $(SCRATCH_DIR)/$(TYPESIZES_MODFILE) \
                           | $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend

.PHONY: depend

SFMAKEDEPEND := ./Bin/sfmakedepend

depend: $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(SCRATCH_DIR)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar -cvf gridpak.tar *

.PHONY: zipfile

zipfile:
		zip -r gridpak.zip *


#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

.PHONY: clean

clean:
	$(RM) -r $(clean_list)

#--------------------------------------------------------------------------
#  A handy debugging target. This will allow to print the value of any
#  Makefile defined macro (see http://tinyurl.com/8ax3j). For example,
#  to find the value of CPPFLAGS execute:
#
#        gmake print-CPPFLAGS
#  or
#        make print-CPPFLAGS
#--------------------------------------------------------------------------

.PHONY: print-%

print-%:
	@echo $* = $($*)
