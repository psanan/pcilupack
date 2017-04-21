PCILUPACK_DIR  := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PCILUPACK_ARCH := $(if $(PETSC_ARCH),$(PETSC_ARCH),build)

PCILUPACK_ILUPACK_DIR ?= ilupack
PCILUPACK_ILUPACK_PLATFORM ?= GNU64

.SECONDEXPANSION:   # to expand $$(@D)/.DIR
.SUFFIXES:          # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:   # Delete likely-corrupt target file if rule fails

OBJDIR ?= $(PCILUPACK_ARCH)/obj
LIBDIR ?= $(PCILUPACK_ARCH)/lib
INCDIR ?= $(PCILUPACK_ARCH)/include

libpcilupack.c := pcilupack.c
libpcilupack.o := $(patsubst %.c,$(OBJDIR)/%.o,$(libpcilupack.c))
libpcilupack.h := pcilupack.h

libpcildl.c := pcildl.c
libpcildl.o := $(patsubst %.c,$(OBJDIR)/%.o,$(libpcildl.c))
libpcildl.h := pcildl.h

include $(PETSC_DIR)/lib/petsc/conf/variables

# $(call SONAME_FUNCTION,libfoo,abiversion)
SONAME_FUNCTION ?= $(1).$(SL_LINKER_SUFFIX).$(2)
# $(call SL_LINKER_FUNCTION,libfoo,abiversion,libversion)
SL_LINKER_FUNCTION ?= -shared  -Wl,-soname,$(call SONAME_FUNCTION,$(notdir $(1)),$(2)) 

PCILUPACK_VERSION_MAJOR := 0
PCILUPACK_VERSION_MINOR := 1
PCILUPACK_VERSION_PATCH := 0

libpcilupack_abi_version := $(PCILUPACK_VERSION_MAJOR).$(PCILUPACK_VERSION_MINOR)
libpcilupack_lib_version := $(libpcilupack_abi_version).$(PCILUPACK_VERSION_PATCH)
soname_function  = $(call SONAME_FUNCTION,$(1),$(libpcilupack_abi_version))
libname_function = $(call SONAME_FUNCTION,$(1),$(libpcilupack_lib_version))
basename_all     = $(basename $(basename $(basename $(basename $(1)))))
sl_linker_args   = $(call SL_LINKER_FUNCTION,$(call basename_all,$@),$(libpcilupack_abi_version),$(libpcilupack_lib_version))

libpcilupack_shared  := $(LIBDIR)/libpcilupack.$(SL_LINKER_SUFFIX) 
libpcilupack_soname  := $(call soname_function,$(LIBDIR)/libpcilupack)
libpcilupack_libname := $(call libname_function,$(LIBDIR)/libpcilupack)
libpcilupack_static  := $(LIBDIR)/libpcilupack.$(AR_LIB_SUFFIX)
libpcilupack         := $(if $(filter-out no,$(BUILDSHAREDLIB)),$(libpcilupack_shared) $(libpcilupack_soname) $(libpcilupack_libname),$(libpcilupack_static))

libpcildl_shared  := $(LIBDIR)/libpcildl.$(SL_LINKER_SUFFIX) 
libpcildl_soname  := $(call soname_function,$(LIBDIR)/libpcildl)
libpcildl_libname := $(call libname_function,$(LIBDIR)/libpcildl)
libpcildl_static  := $(LIBDIR)/libpcildl.$(AR_LIB_SUFFIX)
libpcildl         := $(if $(filter-out no,$(BUILDSHAREDLIB)),$(libpcildl_shared) $(libpcildl_soname) $(libpcildl_libname),$(libpcildl_static))

# ILUPACK headers
PCILUPACK_ILUPACK_INCLUDE := -I${PCILUPACK_ILUPACK_DIR}/include
PCILUPACK_ILUPACK_LIBS := -L${PCILUPACK_ILUPACK_DIR}/lib/${PCILUPACK_ILUPACK_PLATFORM} -lilupack_mc64 -lcamd -lamd -lsuitesparseconfig -lmetisomp -lmetis -lmetisomp -lsparspak -llapack -lblaslike -lblas ${PCILUPACK_ILUPACK_DIR}/notdistributed/MC64D.o ${PCILUPACK_ILUPACK_DIR}/notdistributed/MC21D.o ${PCILUPACK_ILUPACK_DIR}/notdistributed/MC64S.o ${PCILUPACK_ILUPACK_DIR}/notdistributed/MC21S.o
PCILUPACK_OPENMP_FLAG=-fopenmp

all: libpcilupack libpcildl

libpcilupack : $(libpcilupack)
libpcildl : $(libpcildl)

### Rules

PCILUPACK_LIB_DIR       = $(PCILUPACK_DIR)/$(PCILUPACK_ARCH)/lib
PCILUPACK_C_SH_LIB_PATH = $(CC_LINKER_SLFLAG)$(PCILUPACK_LIB_DIR)
PCILUPACK_LIB           = $(PCILUPACK_C_SH_LIB_PATH) -L$(PCILUPACK_LIB_DIR) -lpcilupack $(PETSC_KSP_LIB)

# compile an object from a generated c file
$(OBJDIR)/%.o: %.c | $$(@D)/.DIR
	$(PETSC_COMPILE) $(PCILUPACK_ILUPACK_INCLUDE) $(C_DEPFLAGS) $< -o $@ 


$(libpcilupack_static): $(libpcilupack.o) | $$(@D)/.DIR
	$(AR) $(AR_FLAGS) $@ $^ $(PCILUPACK_ILUPACK_LIBS)
	$(RANLIB) $@

# shared library linking
$(libpcilupack_libname): $(libpcilupack.o) | $$(@D)/.DIR
	$(CLINKER) $(PCILUPACK_OPENMP_FLAG) $(sl_linker_args) -o $@ $^ $(UISCE_EXTERNAL_LIB_BASIC) $(PCILUPACK_ILUPACK_LIBS)
ifneq ($(DSYMUTIL),true)
	$(DSYMUTIL) $@
endif

$(libpcilupack_shared): $(libpcilupack_libname)
	@ln -sf $(notdir $<) $@

$(libpcilupack_soname): $(libpcilupack_libname)
	@ln -sf $(notdir $<) $@



$(libpcildl_static): $(libpcildl.o) | $$(@D)/.DIR
	$(AR) $(AR_FLAGS) $@ $^ $(PCILUPACK_ILUPACK_LIBS)
	$(RANLIB) $@

# shared library linking
$(libpcildl_libname): $(libpcildl.o) | $$(@D)/.DIR
	$(CLINKER) $(PCILUPACK_OPENMP_FLAG) $(sl_linker_args) -o $@ $^ $(UISCE_EXTERNAL_LIB_BASIC) $(PCILUPACK_ILUPACK_LIBS)
ifneq ($(DSYMUTIL),true)
	$(DSYMUTIL) $@
endif

$(libpcildl_shared): $(libpcildl_libname)
	@ln -sf $(notdir $<) $@

$(libpcildl_soname): $(libpcildl_libname)
	@ln -sf $(notdir $<) $@

# make print VAR=the-variable
print:
	@echo $($(VAR))

# make directories that we need
%/.DIR :
	@mkdir -p $(@D)
	@touch $@

clean:
	rm -rf $(OBJDIR) $(LIBDIR) ex50

install:
	./install.py --prefix=$(DESTDIR)

test: test_pcildl test_pcilupack

test_pcildl : ex50
	-@${MPIEXEC} -n 1 ./ex50 -da_grid_x 67 -da_rid_y 79 -dll_append ${PETSC_ARCH}/lib/libpcildl.so -pc_type ildl -ksp_view -ksp_monitor_short > ex50_pcildl.tmp 2>&1; \
	   if (${DIFF} testref/ex50_pcildl.ref ex50_pcildl.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with ex50_pcildl, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex50_pcildl.tmp

test_pcilupack : ex50
	-@${MPIEXEC} -n 1 ./ex50 -da_grid_x 67 -da_rid_y 79 -dll_append ${PETSC_ARCH}/lib/libpcilupack.so -pc_type ilupack -ksp_view -ksp_monitor_short > ex50_pcilupack.tmp 2>&1; \
	   if (${DIFF} testref/ex50_pcilupack.ref ex50_pcilupack.tmp) then true; \
	   else printf "${PWD}\nPossible problem with with ex50_pcilupack, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex50_pcilupack.tmp

ex50.c :
	cp ${PETSC_DIR}/src/ksp/ksp/examples/tutorials/ex50.c .

ex50: ex50.o
	-${CLINKER} -o ex50 ex50.o ${PETSC_KSP_LIB}
	${RM} ex50.o ex50.d

ex50.o : ex50.c
	$(PETSC_COMPILE) $(C_DEPFLAGS) $< -o $@ 

.PHONY: all clean print libpcilupack libpcildl install test test_pcildl test_pcilupack

.PRECIOUS: %/.DIR

libpcilupack.d := $(libpcilupack.o:%.o=%.d)

$(libpcilupack.d) : ;

-include $(libpcilupack.d)

libpcildl.d := $(libpcildl.o:%.o=%.d)

$(libpcildl.d) : ;

-include $(libpcildl.d)
