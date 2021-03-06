CMP = ifort
CMPFLAGS = -extend-source -O4
DEBUG   = -check bounds
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE = -I/usr/local/opt/lapack/include
LAPACK = -Wl,--start-group  -larpack_Intel $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lm -parallel
ARPACK =  -L/opt/ARPACK/
OBJS = besselnew.o matrix_stuff.o zgensub.o DipoleDipole.o logder.o

logder.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o logder.x

logder.o: logder.f90
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c logder.f90

matrix_stuff.o: matrix_stuff.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -extend-source 132 -c besselnew.f

zgensub.o: zgensub.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c zgensub.f

DipoleDipole.mod: DipoleDipole.o
	${CMP} ${FORCEDP} DipoleDipole.o

DipoleDipole.o: DipoleDipole.f90
	${CMP} ${FORCEDP} ${CMPFLAGS} -c DipoleDipole.f90

clean:
	rm -f *.mod *.o *.x
