#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS := $(shell basename `pwd`)
LIB  = lib$(THIS).a

BIN_SRCS = letkf.f90 obsop.f90 oceanobs_nc2bin.f90
BIN      = oceanda.x oceanobs_nc2bin.x oceanobsop.x

esma_install install: $(LIB) $(BIN) 
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p $(BIN)      $(ESMABIN)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := common.f90 common_mpi.f90 common_mtx.f90 mod_sfmt.f90 \
        params_letkf.f90 params_model.f90 params_obs.f90 vars_model.f90 \
        vars_obs.f90 common_letkf.f90 common_mom4.f90 common_mpi_mom4.f90 \
        common_obs_mom4.f90 letkf_obs.f90 letkf_local.f90 \
        letkf_tools.f90 gsw_pot_to_insitu.f90 netlib2.f 

OBJS := $(addsuffix .o, $(basename $(SRCS)))

ALL_SRCS = $(SRCS) $(BIN_SRCS)
DEPS := $(addsuffix .d, $(basename $(ALL_SRCS)))

FOPT = $(FOPT3) $(TRACEBACK)
#FOPT = $(FOPTG) $(TRACEBACK)
FREAL = $(FREAL4)



MOD_DIRS = . $(INC_NETCDF) 
USER_FDEFS = $(D)_LAPACK_ $(D)gmao_intf $(D)HAVE_ESMF
USER_FFLAGS = $(BIG_ENDIAN) $(BYTERECLEN)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINC   = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

oceanda.x:  $(LIB) letkf.o
	$(LD) $(LDFLAGS) -o $@ letkf.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

oceanobsop.x: $(LIB) obsop.o
	$(LD) $(LDFLAGS) -o $@ obsop.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

oceanobs_nc2bin.x: $(LIB) oceanobs_nc2bin.o
	$(LD) $(LDFLAGS) -o $@ oceanobs_nc2bin.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
