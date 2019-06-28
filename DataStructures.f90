!****************************************************************************************************
MODULE DataStructures
  IMPLICIT NONE

  !****************************************************************************************************
  TYPE ScatData

    DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), sigma(:,:)
    complex*16, ALLOCATABLE :: f(:,:), S(:,:), T(:,:)
    double precision delta, tandel, sindel,sin2del
  END TYPE ScatData

CONTAINS
  !****************************************************************************************************
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%R(N,N),SD%sigma(N,N))
    ALLOCATE(SD%f(N,N),SD%S(N,N),SD%T(N,N))
    SD%K=0d0
    SD%R=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
  END SUBROUTINE AllocateScat
  !****************************************************************************************************
  SUBROUTINE DeAllocateScat(SD)
    IMPLICIT NONE
    TYPE(ScatData) SD

    DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T)

  END SUBROUTINE DeAllocateScat

END MODULE DataStructures
