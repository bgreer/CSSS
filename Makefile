
MPIF90 = mpif90

FLAGS = -O3

seismology : seismology.f90 kernels.f90
	${MPIF90} ${FLAGS} -c kernels.f90
	${MPIF90} ${FLAGS} seismology.f90 kernels.f90 -o seismology

clean :
	rm -f *.o *.mod seismology

