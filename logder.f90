!****************************************************************************************************
MODULE GlobalVars
  IMPLICIT NONE
  INTEGER NumParticles, NumChannels,  NumAllChan, NumE
  INTEGER PointsPerBox,NumBoxes,lmax
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION reducedmass, xStart, xEnd, energy
  DOUBLE PRECISION SpatialDim, EffDim, Emin, Emax
  DOUBLE PRECISION, ALLOCATABLE :: mass(:)
  CHARACTER*64 InputFile
  COMPLEX*16 II
  PARAMETER(II=(0.0d0,1.0d0))
  !****************************************************************************************************
 TYPE ScatData
   DOUBLE PRECISION, ALLOCATABLE :: K(:,:), sigma(:,:), sigmatot(:)
   complex*16, ALLOCATABLE :: f(:,:), S(:,:), T(:,:)
 END TYPE ScatData

CONTAINS
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%sigma(N,N))
    ALLOCATE(SD%f(N,N),SD%S(N,N),SD%T(N,N))
    ALLOCATE(SD%sigmatot(0:N))  !for each value of m_l=0,1...l
    SD%K=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
    SD%sigmatot=0d0
  END SUBROUTINE AllocateScat

  SUBROUTINE ReadGlobal
    IMPLICIT NONE
    INTEGER n
    ! Be sure to match the format of the input file to the format of the read statements below
    OPEN(unit=7,file=InputFile(1:INDEX(InputFile,' ')-1),action='read')
    READ(7,*)
    READ(7,*) NumParticles, NumChannels, SpatialDim, NumAllChan
    ALLOCATE(mass(NumParticles))
    READ(7,*)
    READ(7,*)
    READ(7,*) (mass(n), n=1,NumParticles)
    READ(7,*)
    READ(7,*)
    READ(7,*) xStart, xEnd, xNumPoints
    READ(7,*)
    READ(7,*)
    READ(7,*) NumBoxes
    READ(7,*)
    READ(7,*)
    READ(7,*) Emin,  Emax,  NumE
    PointsPerBox=9  ! careful! if you change this you need to go change the weights array.
    lmax = 2*NumChannels
    CLOSE(unit=7)
    EffDim = NumParticles*SpatialDim - SpatialDim
    !AlphaFactor = 0d0
    AlphaFactor = (EffDim-1d0)/2d0

    IF (NumParticles.EQ.2) THEN
       reducedmass = mass(1)*mass(2)/(mass(1)+mass(2))
    ELSE
       WRITE(6,*) "Reduced mass not set. Must set reduced mass"
       STOP
    END IF
  END SUBROUTINE ReadGlobal
    !****************************************************************************************************
END MODULE GlobalVars

module logderprop
  use GlobalVars
  implicit none
  double precision h ! step size
  double precision, allocatable :: weights(:) ! weights for pointesperbox points for each box

  double precision, allocatable :: identity(:,:) !identity matrix
  double precision, allocatable :: ycurrent(:,:), yprevious(:,:) ! y_(n} and y_{n-1}
  double precision, allocatable :: u(:,:,:) !Johnson's u matrix for indices i,j,n
  double precision, allocatable :: VV(:,:,:) ! Johnson's curly V matrix.  indices for i, j, n

contains
  subroutine allocateprop
    implicit none


    allocate(VV(NumChannels,NumChannels,PointsPerBox))
    allocate(u(NumChannels,NumChannels,PointsPerBox))
    allocate(identity(NumChannels,NumChannels))
    allocate(ycurrent(NumChannels,NumChannels))
    allocate(yprevious(NumChannels,NumChannels))
    allocate(weights(PointsPerBox))
  end subroutine allocateprop
  subroutine initprop
    implicit none
    integer i,j

    call allocateprop

    weights(1)=1d0
    weights(2)=4d0
    weights(3)=2d0
    weights(4)=4d0
    weights(5)=2d0
    weights(6)=4d0
    weights(7)=2d0
    weights(8)=4d0
    weights(9)=1d0
    identity = 0d0
    yprevious = 0d0
    ycurrent = 0d0
    do i=1,NumChannels
      identity(i,i)=1d0
      yprevious(i,i)=1d20
    enddo

  end subroutine initprop


  subroutine boxstep(x,yi,yf)
    double precision x(pointsperpot)
    integer i,j,step
    double precision y1(NumChannels,NumChannels),yf(NumChannels,NumChannels)





  end subroutine boxstep

end module logderprop

subroutine MakePot(Pot,x,DP,BoxGrid)
  use GlobalVars
  use DipoleDipole
  implicit none
  double precision Pot(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax),BoxGrid(NumBoxes+1)
  double precision x(PointsPerBox,NumBoxes),xm2(PointsPerBox,NumBoxes),xm3(PointsPerBox,NumBoxes)
  integer mch,nch,ml,l,lp,step,ibox

  do ibox=1,NumBoxes
    do step=1,PointsPerBox
      call GridMaker(x(:,iBox),BoxGrid(iBox),BoxGrid(iBox+1),"linear")
      xm2(:,iBox)=1d0/x(:,iBox)**2
      xm3(:,iBox)=1d0/x(:,iBox)**3
    enddo
  enddo
  do ml=0,lmax
    DO mch = 1,NumChannels
      DO nch = 1, mch
        IF(DP%even) THEN
          l=2*(mch-1) ! l = 0, 2, 4, ...
          lp=2*(nch-1) ! l' = 0, 2, 4, ...
        ELSE
          l=2*(mch-1)+1 ! l = 1, 3, 5, ...
          lp=2*(nch-1)+1 ! l' = 1, 3, 5, ...
        ENDIF
        do iBox=1,NumBoxes
          do step=1,PointsPerBox
            Pot(mch,nch,step,iBox,ml) = -2d0*DP%Cllp(l,lp,DP%ml)*xm3(step,iBox)**3
            IF(mch.EQ.nch) Pot(mch,nch,step,iBox,ml) = Pot(mch,nch,step,iBox,ml) + 0.5d0*l*(l+1)*xm2(step,iBox)
            Pot(nch,mch,step,iBox,ml) = Pot(mch,nch,step,iBox,ml)
          enddo
        enddo
    ENDDO
  ENDDO

end subroutine MakePot

program main
  use GlobalVars
  implicit none
  type(DPData) DP
  double precision PotOdd(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax)
  double precision PotEven(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax)
  double precision BoxGrid(NumBoxes+1)

  InputFile = 'logder.inp'
  call ReadGlobal()
  call AllocateDP(DP,NumChannels)
  call MakeDipoleDipoleCouplingMatrix(DP)
  DP%even = .true.
  call MakePot(PotEven,DP)
  DP%even = .false.
  call MakePot(PotOdd,DP)


end program
SUBROUTINE GridMaker(grid,numpts,E1,E2,scale)
  implicit none
  DOUBLE PRECISION grid(numpts)
  DOUBLE PRECISION E1,E2,LE1,LE2,DE,LDE
  INTEGER numpts, iE
  CHARACTER(LEN=*), INTENT(IN) :: scale
  !--------------------------------------------
  ! Linear grid:
  !--------------------------------------------
  grid(1)=E1
  IF((scale.EQ."linear").and.(numpts.gt.1)) THEN
     DE=(E2-E1)/DBLE(numpts-1)
     DO iE=1,numpts
        grid(iE) = E1 + (iE-1)*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! Log grid:
  !--------------------------------------------
  IF((scale.EQ."log").and.(numpts.gt.1)) THEN
     LE1=dlog(E1)
     LE2=dlog(E2)

     LDE=(LE2-LE1)/DBLE(numpts-1d0)
     DO iE=1,numpts
        grid(iE) = dexp(LE1 + (iE-1)*LDE)
!        write(6,*) LE1, LE2, LDE, grid(iE)
     ENDDO
  ENDIF
  !--------------------------------------------
  ! Quadratic grid:
  !--------------------------------------------
  IF((scale.EQ."quadratic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**2*DE
     ENDDO
  ENDIF
END SUBROUTINE GridMaker
