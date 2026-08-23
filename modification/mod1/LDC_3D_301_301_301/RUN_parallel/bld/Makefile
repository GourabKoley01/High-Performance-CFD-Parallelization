PLATFORM = "roddam"
EXENAME = VYOM3D_mod1.exe
#-----------
# Compilers
#-----------
ifeq ($(PLATFORM),"roddam")
   FC = mpif90
   CC =
   FFLAGS =-byteswapio -r8 -march=native -acc -gpu=cc80 -Minfo=all
#   FFLAGS =-O3 -march=native -fconvert=big-endian -acc -ta=tesla:cc80 -Minfo=all
#   FFLAGS = -traceback -g -ftrapuv \
#       -r8 -convert big_endian -fp-model strict \
#	-check all -debug all -warn all
#	 FFLAGS = -O2 -traceback -g -r8 -convert big_endian \
						-fpe0
endif

#---------------------
# Program directories
#---------------------
#MPILDC = /work/ranjan.25/RR_LDC/FDL3D_NEWMAKE
MPILDC = ./
SRCDIR = $(MPILDC)/src
BLDDIR = .
EXEDIR = ../

#----------------------------------
# Rules for compiling source files
#----------------------------------
#
# Load dependency search path
#
VPATH := $(shell cat $(BLDDIR)/Filepath)
#
# Build the model
#
all: $(EXEDIR)/$(EXENAME)
#
# Source files are listed in $(BLDDIR)/Srcfiles
#
SRCS := $(shell cat $(BLDDIR)/Srcfiles)
#
# Execute mkDepends to create $(DEPS) (to ensure that new dependencies
# are taken into account)
#
DEPS = depends # Contains list of dependencies
$(BLDDIR)/$(DEPS): $(BLDDIR)/Srcfiles $(BLDDIR)/Filepath
	$(BLDDIR)/mkDepends Filepath Srcfiles > $@
#
# Execute mkSrcfiles to create 'Srcfiles' (to ensure that entire list
# of source files are listed in 'Srcfiles')
#
$(BLDDIR)/Srcfiles: $(BLDDIR)/Filepath
	$(BLDDIR)/mkSrcfiles > $@

OBJS := $(addsuffix .o, $(basename $(SRCS)))

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) \
        -L/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-8.5.0/nvhpc-24.11-bphfxrlohbqt7qd5igbuvb7va2ve6znl/Linux_x86_64/24.11/cuda/11.8/targets/x86_64-linux/lib \
        -L/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-8.5.0/nvhpc-24.11-bphfxrlohbqt7qd5igbuvb7va2ve6znl/Linux_x86_64/24.11/compilers/lib \
        -lnvToolsExt -lnvhpcwrapnvtx \
        -acc -gpu=cc80
#	$(FC) -o $@ $(OBJS) -L /opt/nvidia/hpc_sdk/Linux_x86_64/24.9/cuda/12.6/lib64 -lnvToolsExt -lnvhpcwrapnvtx  -acc    #gpu02 
#	$(FC) -o $@ $(OBJS) -L /opt/hpc-sdk-24.3/Linux_x86_64/24.3/cuda/12.3/lib64 -lnvToolsExt -lnvhpcwrapnvtx -acc     #airawat

.SUFFIXES:
.SUFFIXES: .f90 .c .o

.f90.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

#----------
# Clean up
#----------
.PHONY: clean

clean:
	rm -f Srcfiles $(DEPS) \
	rm -f $(BLDDIR)/*.mod $(BLDDIR)/*.o \
	rm -f *genmod* \ 
#	rm -f $(EXEDIR)/$(EXENAME)

include $(BLDDIR)/$(DEPS)
