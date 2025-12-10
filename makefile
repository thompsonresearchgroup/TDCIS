mqcroot = $(HOME)/mqcPack
blalibs = mkl
FC = gfortran

ifeq ($(blalibs),blalap)
  ifeq ($(FC),gfortran)
	FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
	MQCOBJS = -I$(mqcroot)/GNU/mod
	LIBS = $(mqcroot)/GNU/lib/libmqc.a -llapack -lblas
  else ifeq ($(FC),pgfortran)
	FCFLAGS = -Mallocatable=03 -r8 -i8 -mp
	MQCOBJS = -module $(mqcroot)/PGI/mod
	LIBS = $(mqcroot)/PGI/lib/libmqc.a -llapack -lblas
  endif
else ifeq ($(blalibs),mkl)
  ifeq ($(FC),gfortran)
	FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
	MQCOBJS = -I$(mqcroot)/GNU/mod
	LIBS = $(mqcroot)/GNU/lib/libmqc.a -L/opt/intel/oneapi/mkl/2022.2.0/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl
  else ifeq ($(FC),pgfortran)
	FCFLAGS = -Mallocatable=03 -r8 -i8 -mp
	MQCOBJS = -module $(mqcroot)/PGI/mod
	LIBS = $(mqcroot)/PGI/lib/libmqc.a -L/opt/intel/oneapi/mkl/2022.2.0/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl
  endif
endif

all: compile

compile: tdcis.f03
	$(FC) $(FCFLAGS) $(MQCOBJS) -o tdcis.exe tdcis.f03 $(LIBS)

clean:
	rm -f tdcis.exe tdcis.o