CMP = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3 -fdec
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE = -I/usr/local/opt/lapack/include
LAPACK =  -framework accelerate
ARPACK = -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX
OBJS = DataStructures.o besselnew.o matrix_stuff.o zgensub.o MorsePotential.o logder.o

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
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

zgensub.o: zgensub.f
	${CMP} ${FORCEDP} -c zgensub.f

MorsePotential.mod: MorsePotential.o
	${CMP} ${FORCEDP} MorsePotential.o

MorsePotential.o: MorsePotential.f90
	${CMP} ${FORCEDP} -c MorsePotential.f90

DataStructures.mod: DataStructures.o
	${CMP} ${FORCEDP} DataStructures.o

DataStructures.o: DataStructures.f90
	${CMP} ${FORCEDP} -c DataStructures.f90

clean:
	rm -f *.mod *.o *.x
