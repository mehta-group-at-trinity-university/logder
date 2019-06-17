CMP = ifort
CMPFLAGS = -extend-source -O4
DEBUG   = -check bounds
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE = -I/usr/local/opt/lapack/include
LAPACK = -Wl,--start-group  -larpack_Intel $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lm -parallel
ARPACK =  -L/opt/ARPACK/
OBJS = DataStructures.o besselnew.o matrix_stuff.o zgensub.o  DipoleDipole.o logderdipole.o

logderdipole.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o logderdipole.x

logderdipole.o: logderdipole.f90
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c logderdipole.f90

matrix_stuff.o: matrix_stuff.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

zgensub.o: zgensub.f
	${CMP} ${FORCEDP} -c zgensub.f

DipoleDipole.mod: DipoleDipole.o
	${CMP} ${FORCEDP} DipoleDipole.o

DipoleDipole.o: DipoleDipole.f90
	${CMP} ${FORCEDP} -c DipoleDipole.f90

DataStructures.mod: DataStructures.o
	${CMP} ${FORCEDP} DataStructures.o

DataStructures.o: DataStructures.f90
	${CMP} ${FORCEDP} -c DataStructures.f90

clean:
	rm -f *.mod *.o *.x


