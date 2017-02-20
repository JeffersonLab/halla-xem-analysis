## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## ARGONNE DEFAULT SETUP FLAGS:
#simcdir = .
#Csoft = /disk1/users/reinhold/Csoft/05Dec1996
#CERN_ROOT = /disk1/lib/cern/97a

## CEBAF DEFAULT SETUP FLAGS: ALSO NEED TO UNCOMMENT CERN_ROOT DEFNS FOR HP/SUN/...
simcdir = .

## THE REST SHOULD BE OK WITHOUT MODIFICATION.

RM        = rm -f 
SHELL     = /bin/sh

SPEC      = ./hrsr/
SPECB     = ./hrsl/

SHARE	  = ./shared/

my_objs = $(SPEC)mc_hrsr.o $(SPEC)mc_hrsr_hut.o $(SPEC)mc_hrsr_recon.o \
	  $(SPECB)mc_hrsl.o $(SPECB)mc_hrsl_hut.o $(SPECB)mc_hrsl_recon.o \
	  $(SHARE)project.o $(SHARE)transp.o $(SHARE)rotate_haxis.o \
	  $(SHARE)locforunt.o $(SHARE)musc.o $(SHARE)musc_ext.o \
	  $(SHARE)stringlib.o mt19937.o gauss1.o loren.o enerloss_new.o

sub_objs = $(SPEC)mc_hrsr.o $(SPEC)mc_hrsr_hut.o $(SPEC)mc_hrsr_recon_sub.o \
	  $(SPECB)mc_hrsl.o $(SPECB)mc_hrsl_hut.o $(SPECB)mc_hrsl_recon_sub.o \
	  $(SHARE)project.o $(SHARE)transp_sub.o $(SHARE)rotate_haxis.o \
	  $(SHARE)locforunt.o $(SHARE)musc.o $(SHARE)musc_ext.o \
	  $(SHARE)stringlib.o mt19937.o gauss1.o loren.o enerloss_new.o

MYOS := $(subst -,,$(shell uname))
#mess $(MYOS)

#CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
CERNLIBS = -Wl,-static -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib -Wl,-dy
ifeq ($(MYOS),HPUX)
  ifneq (,$(findstring 09,$(shell uname -r)))
    HPUXVERSION := 09
  else
    HPUXVERSION := 10
  endif
  LIBROOT = $(Csoft)/../$(MYOS)$(HPUXVERSION)/lib
else
  LIBROOT = $(Csoft)/../$(MYOS)/lib
endif

ifeq ($(MYOS),HPUX)
  CERN_ROOT = /site/cernlib/hp700_ux90/96a
  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = \
	-L$(CERN_ROOT)/lib $(CERNLIBS) \
	-L/usr/lib/X11R5 -lX11 -lm
endif

#For use with gfortran compiler
# -fno-automatic - all program storage treated as static
ifeq ($(MYOS),Linux)
#  CERN_ROOT = /apps/cernlib/i386_fc8/2005
  FFLAGSA=-O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib64
  FC  := gfortran
  F77 := gfortran
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS)
endif

ifeq ($(MYOS),SunOS)
  CERN_ROOT = /site/cernlib/sun4_solaris2/97a
  FFLAGS=-e -O -I$(Csoft)/SRC/INCLUDE -ftrap=underflow
  ifeq ($(OSTYPE),SunOS4)
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lX11
  else
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lsocket -lX11
  endif
endif

ifeq ($(MYOS),AIX)
  F77=f77
  FFLAGS=-g -qfixed=132 -qextname -O -I$(Csoft)/SRC/INCLUDE
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) -lX11
endif

ifeq ($(MYOS),OSF1)
  F77=f77
  CERN_ROOT = /disk1/lib/cern/new
  LIBROOT = $(Csoft)/OSF1/lib
  FFLAGS= -r8 -extend_source -Wl,-taso -I -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused
  LDFLAGS= 
  OTHERLIBS = -Wl,-L$(CERN_ROOT)/lib \
        -lpacklib $(CERNLIBS) -Wl,-L/usr/lib/X11R5 -lX11 -lm 
endif

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

none: mc_hrs_single

all: mc_hrs_single

mc_hrs_single: mc_hrs_single.o $(my_objs) Makefile
	$(F77) -o $@ $(FFLAGS) $(my_objs) mc_hrs_single.o $(OTHERLIBS)

clean:
	rm -f *.o $(SPEC)*.o $(SPECB)*.o $(SHARE)*.o mc_hrs_single
