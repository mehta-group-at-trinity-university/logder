!****************************************************************************************************
MODULE GlobalVars
  IMPLICIT NONE
  INTEGER NumParticles, NumChannels,  NumAllChan, NumE
  INTEGER PointsPerBox,NumBoxes,lmax
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION mu, xStart, xEnd, energy,onethird
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
    READ(7,*) xStart, xEnd, PointsPerBox
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
       mu = mass(1)*mass(2)/(mass(1)+mass(2))
    ELSE
       WRITE(6,*) "Reduced mass not set. Must set reduced mass"
       STOP
    END IF

    onethird = 0.3333333333333333333333333d0
  END SUBROUTINE ReadGlobal
    !****************************************************************************************************
END MODULE GlobalVars

module logderprop
  use GlobalVars
  implicit none
  double precision h ! step size
  double precision, allocatable :: weights(:) ! weights for pointesperbox points for each box
  double precision, allocatable :: identity(:,:) !identity matrix
  double precision, allocatable :: ycurrent(:,:), yprevious(:,:),ystart(:,:) ! y_(n} and y_{n-1}
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
    allocate(ystart(NumChannels,NumChannels))
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
      ystart(i,i)=1d20
    enddo

  end subroutine initprop

  subroutine boxstep(xx,yi,yf,Pot)
    implicit none
    integer i,j,step
    double precision h
    double precision xx(PointsPerBox)
    double precision yi(NumChannels,NumChannels),yf(NumChannels,NumChannels)
    double precision Pot(NumChannels,NumChannels,PointsPerBox) !make sure Pot includes the threshold offsets
    double precision temp(NumChannels,NumChannels)

    h=xx(2)-xx(1)

    yprevious = yi
    do step = 2, PointsPerBox-1,2
        VV(:,:,step)=2d0*mu*(identity*Energy-Pot(:,:,step))
        u(:,:,step)=VV(:,:,step)
        VV(:,:,step+1)=2d0*mu*(identity*Energy-Pot(:,:,step+1))
        temp = identity + h**2/6d0 * VV(:,:,step+1)
        call sqrmatinv(temp,NumChannels)
        u(:,:,step+1)=MATMUL(temp,VV(:,:,step+1))
    enddo

    do step = 2, PointsPerBox
      temp = identity + h*yprevious
      call sqrmatinv(temp,NumChannels)
      ycurrent = MATMUL(temp,yprevious) - onethird*h*weights(step)*u(:,:,step)
      yprevious = ycurrent
    enddo
    yf = ycurrent

  end subroutine boxstep
end module logderprop
!============================================================================================
module scattering
  !****************************************************************************************************
 TYPE ScatData

   DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), sigma(:,:), sigmatot(:)
   complex*16, ALLOCATABLE :: f(:,:), S(:,:), T(:,:)
 END TYPE ScatData

CONTAINS
 !****************************************************************************************************
 SUBROUTINE AllocateScat(SD,N)
   IMPLICIT NONE
   TYPE(ScatData) SD
   INTEGER N
   ALLOCATE(SD%K(N,N),SD%R(N,N),SD%sigma(N,N))
   ALLOCATE(SD%f(N,N),SD%S(N,N),SD%T(N,N))
   ALLOCATE(SD%sigmatot(-N:N))
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
   DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T,SD%sigmatot)
 END SUBROUTINE DeAllocateScat
  !****************************************************************************************************
 SUBROUTINE CalcK(B,BPD,SD,mu,d,alpha,EE,Eth)
   IMPLICIT NONE
   TYPE(ScatData) :: SD
   DOUBLE PRECISION mu, EE,rm,d,alpha
   DOUBLE PRECISION, ALLOCATABLE :: s(:),c(:),sp(:),cp(:)
   DOUBLE PRECISION, ALLOCATABLE :: Imat(:,:),Jmat(:,:)
   DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp,Pi
   DOUBLE PRECISION k(BPD%NumChannels),Eth(BPD%NumChannels)
   complex*16, allocatable :: tmp(:,:),Identity(:,:)
   complex*16  II
   INTEGER i,j,no,nw,nc,beta

   II=(0d0,1d0)
   Pi=dacos(-1d0)
   rm=BPD%xr
   no=0
   nw=0
   nc=0

   DO i = 1,BPD%NumChannels
      IF (EE.GE.Eth(i)) THEN
         k(i) = dsqrt(2d0*mu*(EE-Eth(i))) ! k is real
         no=no+1
      ELSE
         k(i) = dsqrt(2d0*mu*(Eth(i)-EE)) ! k->kappa, kappa is real
         IF( (k(i)*rm).LT.10d0) nw = nw+1
         IF( (k(i)*rm).GE.10d0) nc = nc+1
      ENDIF
   ENDDO
   !      write(6,*) "no = ", no
   IF((no+nw+nc).NE.BPD%NumChannels) THEN
      WRITE(6,*) "Channel miscount in calcK"
      STOP
   ENDIF

   ALLOCATE(s(no),c(no),Imat(no,no),Jmat(no,no),tmp(no,no))
   ALLOCATE(sp(no),cp(no))
   allocate(Identity(no,no))
   Identity = 0d0;

   DO i = 1,no
     Identity(i,i) = 1d0
     CALL hyperrjry(INT(d),alpha,BPD%lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
     s(i) = rhypj/dsqrt(Pi*k(i))
     c(i) = -rhypy/dsqrt(Pi*k(i))
     sp(i) = dsqrt(k(i)/Pi)*rhypjp
     cp(i) = -dsqrt(k(i)/Pi)*rhypyp
   ENDDO
   Imat=0d0
   Jmat=0d0
   DO i=1,no
      DO beta=1,no
         Imat(i,beta) = (B%Zf(i,beta)*cp(i) - B%Zfp(i,beta)*c(i))/(s(i)*cp(i)-c(i)*sp(i))
         Jmat(i,beta) = (B%Zf(i,beta)*sp(i) - B%Zfp(i,beta)*s(i))/(c(i)*sp(i)-s(i)*cp(i))
      ENDDO
   ENDDO

   CALL SqrMatInv(Imat, no)
   SD%K = MATMUL(Jmat,Imat)

   SD%R=0.0d0
   DO i=1,B%NumOpenR
      DO j=1,B%NumOpenR
         DO beta = 1, B%NumOpenR
            SD%R(i,j) = SD%R(i,j) - B%Zf(i,beta)*B%Zf(j,beta)/B%bf(beta) !
         ENDDO
      ENDDO
   ENDDO

   tmp = Identity - II*SD%K
   !write(6,*) "Identity matrix:"
   !call printmatrix(Identity,no,no,6)
   !write(6,*) "K matrix:"
   !call printmatrix(SD%K,no,no,6)
   !write(6,*) "real part of 1 - II*K"
   !call printmatrix(realpart(tmp),no,no,6)
   !write(6,*) "Imaginary part of 1 - II*K"
   !call printmatrix(aimag(tmp),no,no,6)
   SD%S = Identity + II*SD%K
   call CompSqrMatInv(tmp,no)
   SD%S = MATMUL(SD%S,tmp)
   SD%T = -II*0.5d0*(SD%S-Identity)

   SD%sigma = conjg(SD%T)*SD%T*Pi/(2d0*mu*EE)

   DEALLOCATE(s,c,Imat,Jmat,sp,cp)
 END SUBROUTINE CalcK
end module scattering
!=========================================================================================
program main
  use GlobalVars
  use DipoleDipole
  use logderprop
  implicit none
  type(DPData) DP
  double precision, allocatable :: PotOdd(:,:,:,:,:), PotEven(:,:,:,:,:), BoxGrid(:)
  double precision, allocatable :: x(:,:),yin(:,:),yout(:,:)
  double precision time1, time2
  integer iBox,ml

  InputFile = 'logder.inp'
  call ReadGlobal()
  DP%lmax=lmax
  allocate(PotOdd(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax))
  allocate(PotEven(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax))
  allocate(BoxGrid(NumBoxes+1))
  allocate(yin(NumChannels,NumChannels),yout(NumChannels,NumChannels))
  allocate(x(PointsPerBox,NumBoxes))
  call AllocateDP(DP,NumChannels)
  call MakeDipoleDipoleCouplingMatrix(DP)
  write(6,*) "cllp, m=0"
  call printmatrix(DP%cllp(:,:,0),NumChannels,NumChannels,6)
  call GridMaker(BoxGrid,NumBoxes+1,xStart,xEnd,"quadratic")
  call printmatrix(BoxGrid,NumBoxes+1,1,6)

  call cpu_time(time1)
  DP%even = .true.
  call MakePot(PotEven,x,DP,BoxGrid)
  DP%even = .false.
  call MakePot(PotOdd,x,DP,BoxGrid)
  call cpu_time(time2)
  write(6,*) "time for potential calculation:", time2-time1

  call initprop

  do ml=0,lmax
    DP%ml=ml
    DP%even = .true.
    yin = ystart
    do iBox=1,NumBoxes
      call boxstep(x(:,iBox),yin,yout,PotEven(:,:,:,iBox,ml))
      yin = yout
      write(6,*)"ml, iBox=", ml, iBox
    enddo
  enddo

end program
!=========================================================================================
!=========================================================================================
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


subroutine MakePot(Pot,x,DP,BoxGrid)
  use GlobalVars
  use DipoleDipole
  implicit none
  double precision Pot(NumChannels,NumChannels,PointsPerBox,NumBoxes,0:lmax),BoxGrid(NumBoxes+1)
  double precision x(PointsPerBox,NumBoxes),xm2(PointsPerBox,NumBoxes),xm3(PointsPerBox,NumBoxes)
  integer mch,nch,ml,l,lp,step,ibox
  type(DPData) DP

  do ibox=1,NumBoxes
    do step=1,PointsPerBox
      call GridMaker(x(:,iBox),PointsPerBox,BoxGrid(iBox),BoxGrid(iBox+1),"linear")
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
            Pot(mch,nch,step,iBox,ml) = -2d0*DP%Cllp(l,lp,ml)*xm3(step,iBox)**3
            IF(mch.EQ.nch) Pot(mch,nch,step,iBox,ml) = Pot(mch,nch,step,iBox,ml) + 0.5d0*l*(l+1)*xm2(step,iBox)
            Pot(nch,mch,step,iBox,ml) = Pot(mch,nch,step,iBox,ml)
          enddo
        enddo
      ENDDO
    ENDDO
  enddo
end subroutine MakePot
!=========================================================================================
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)

  DO j = 1,nr
     WRITE(file,20) (M(j,k), k = 1,nc)
  ENDDO

20 FORMAT(1P,100D20.12)
30 FORMAT(100F12.6)
END SUBROUTINE printmatrix
!=========================================================================================
SUBROUTINE checkpot(Pot,x,file)
  use GlobalVars
  IMPLICIT NONE
  INTEGER file,mch,nch,step,iBox
  double precision x(PointsPerBox,NumBoxes),Pot(NumChannels,NumChannels,PointsPerBox,NumBoxes)

  DO mch = 1,NumChannels
     DO nch = 1,mch
        do iBox=1,NumBoxes
           DO step = 1,PointsPerBox
              WRITE(file,*) x(step,iBox), Pot(mch,nch,step,iBox)
           ENDDO
        ENDDO
        WRITE(file,*)
     ENDDO
     WRITE(file,*)
  ENDDO

END SUBROUTINE checkpot
