local_sub  := Drivers

local_src  := $(wildcard $(local_sub)/*.F)
local_objs := $(subst .F,.o,$(local_src))
local_objs := $(addprefix $(SCRATCH_DIR)/, $(notdir $(local_objs)))

sources    += $(local_src)

ifeq ($(OS)-$(strip $(FORT)),CYGWIN-df)
$(COAST):	$(libraries) $(SCRATCH_DIR)/coast.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/coast.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(GRID):	$(libraries) $(SCRATCH_DIR)/grid.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/grid.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(TOLAT):	$(libraries) $(SCRATCH_DIR)/tolat.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/tolat.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(BATHTUB):	$(libraries) $(SCRATCH_DIR)/bathtub.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/bathtub.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(BATHSUDS):	$(libraries) $(SCRATCH_DIR)/bathsuds.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/bathsuds.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(BATHSOAP):	$(libraries) $(SCRATCH_DIR)/bathsoap.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/bathsoap.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
$(SPHERE):	$(libraries) $(SCRATCH_DIR)/sphere.o
	$(LD) $(FFLAGS) $(SCRATCH_DIR)/sphere.o $(libraries) \
	    $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
else
$(COAST):	$(libraries) $(SCRATCH_DIR)/coast.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/coast.o -o $@ \
	    $(libraries) $(LIBS)
$(GRID):	$(libraries) $(SCRATCH_DIR)/grid.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/grid.o -o $@ \
	    $(libraries) $(LIBS)
$(TOLAT):	$(libraries) $(SCRATCH_DIR)/tolat.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/tolat.o -o $@ \
	    $(libraries) $(LIBS)
$(BATHTUB):	$(libraries) $(SCRATCH_DIR)/bathtub.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/bathtub.o -o $@ \
	    $(libraries) $(LIBS)
$(BATHSUDS):	$(libraries) $(SCRATCH_DIR)/bathsuds.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/bathsuds.o -o $@ \
	    $(libraries) $(LIBS)
$(BATHSOAP):	$(libraries) $(SCRATCH_DIR)/bathsoap.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/bathsoap.o -o $@ \
	    $(libraries) $(LIBS)
$(SPHERE):	$(libraries) $(SCRATCH_DIR)/sphere.o
	$(LD) $(FFLAGS) $(LDFLAGS) $(SCRATCH_DIR)/sphere.o -o $@ \
	    $(libraries) $(LIBS)
endif

$(eval $(compile-rules))
