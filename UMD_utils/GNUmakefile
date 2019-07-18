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
       ESMADIR := $(PWD)/../../../../..
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

BIN_SRCS = ocean_recenter.f90 anaice2rst.f90 ocean_moments.f90 ocean_iau.f90
BIN      = ocean_recenter.x anaice2rst.x ocean_moments.x ocean_iau.x read_merra2_bcs.so

PYTHONFILES = $(wildcard *.py)
PYTHONSCRIPT = $(wildcard $(PYTHONFILES))

pythonscript = $(PYTHONSCRIPT)
SCRIPTS  = $(PYTHONSCRIPT)

scripts  = $(SCRIPTS)

esma_install install: $(LIB) $(BIN) 
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p $(BIN) $(PYTHONSCRIPT)     $(ESMABIN)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *.so

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

pythonscript:
	@echo ""
	@for scr in $($@); do  \
           ( $(SED) -e "$(SEDCMD)"< $$scr > $(ESMABIN)/$$scr; \
             echo "perl script to \$$(ESMABIN): $$scr";  \
             chmod 755 $(ESMABIN)/$$scr ) \
        done
	@echo ""
	@echo "\$$(ESMABIN) = $(ESMABIN)"
	@echo ""

scripts:
	@echo ""
	@for scr in $($@); do  \
	   ( $(SED) -e "$(SEDCMD)"< $$scr > $(ESMABIN)/$$scr; \
	     echo "script to \$$(ESMABIN): $$scr";  \
	     chmod 755 $(ESMABIN)/$$scr ) \
	done
	@echo ""
	@echo "\$$(ESMABIN) = $(ESMABIN)"
	@echo ""


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

SRCS := ocean_recenter_utils.f90 ocnice_utils.f90
#SRCS := ocean_recenter_utils.f90 ocnice_utils.f90 discrete_operators.f90 interp2d.f
OBJS := $(addsuffix .o, $(basename $(SRCS)))

ALL_SRCS = $(SRCS) $(BIN_SRCS)
DEPS := $(addsuffix .d, $(basename $(ALL_SRCS)))

FOPT = $(FOPT3) $(TRACEBACK)
FREAL = $(FREAL4)

MOD_DIRS = . $(INC_NETCDF) $(INC_MPI)
USER_FDEFS = $(D)_LAPACK_ $(D)gmao_intf $(D)HAVE_ESMF
USER_FFLAGS = $(BIG_ENDIAN) $(BYTERECLEN)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINC   = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

ocean_recenter.x: $(LIB) ocean_recenter.o
	$(LD) $(LDFLAGS) -o $@ ocean_recenter.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

anaice2rst.x: $(LIB) anaice2rst.o
	$(LD) $(LDFLAGS) -o $@ anaice2rst.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

ocean_moments.x: $(LIB) ocean_moments.o
	$(LD) $(LDFLAGS) -o $@ ocean_moments.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

ocean_iau.x: $(LIB) ocean_iau.o
	$(LD) $(LDFLAGS) -o $@ ocean_iau.o $(LIB) \
              $(LIB_NETCDF) \
              $(LIB_MPI) $(LIB_SYS) $(LIB_SCI)

read_merra2_bcs.so: read_bin.f90
	f2py -c -m read_merra2_bcs read_bin.f90

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
