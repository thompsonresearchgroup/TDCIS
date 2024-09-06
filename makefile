mqcroot = "$(HOME)/.mqc"
blalibs = blalap
FC = gfortran
ifeq ($(blalibs),blalap)
  ifeq ($(FC),gfortran)
	FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
	MQCOBJS = -I$(mqcroot)/GNU/mod
	LIBS = -llapack -lblas $(mqcroot)/GNU/lib/libmqc.a
  else ifeq ($(FC),pgfortran)
	FCFLAGS = -Mallocatable=03 -r8 -i8 -mp
	MQCOBJS = -module $(mqcroot)/PGI/mod
	LIBS = -llapack -lblas $(mqcroot)/PGI/lib/libmqc.a
  endif
else ifeq ($(blalibs),mkl)
  ifeq ($(FC),gfortran)
	FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl
	MQCOBJS = -I$(mqcroot)/GNU/mod
	LIBS =  $(mqcroot)/GNU/lib/libmqc.a -L/opt/intel/oneapi/mkl/2022.2.0/lib/intel64
  else ifeq ($(FC),pgfortran)
	FCFLAGS = -Mallocatable=03 -r8 -i8 -mp -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl
	MQCOBJS = -module $(mqcroot)/PGI/mod
	LIBS =  $(mqcroot)/GNU/lib/libmqc.a -L/opt/intel/oneapi/mkl/2022.2.0/lib/intel64
  endif
endif

all: compile

compile: tdcis.f03
	$(FC) $(FCFLAGS) $(MQCOBJS) -o tdcis.exe tdcis.f03 $(LIBS) 

clean:
	rm -f tdcis.exe tdcis.o
