!****************************************************************************************************

MODULE GlobalVars
  use datastructures
  IMPLICIT NONE
  INTEGER NumParticles, NumAllChan, NumE, NumChannels
  INTEGER PointsPerBox,NumBoxes,lmax
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION mu, xStart, xEnd, energy,onethird,onesixth
  DOUBLE PRECISION SpatialDim, EffDim, Emin, Emax
  DOUBLE PRECISION, ALLOCATABLE :: mass(:)
  CHARACTER*64 InputFile
  COMPLEX*16 II
  PARAMETER(II=(0.0d0,1.0d0))
  !****************************************************************************************************

CONTAINS

  SUBROUTINE ReadGlobal
    IMPLICIT NONE
    INTEGER n
    ! Be sure to match the format of the input file to the format of the read statements below
    OPEN(unit=7,file=InputFile(1:INDEX(InputFile,' ')-1),action='read')
    READ(7,*)
    READ(7,*) NumParticles, SpatialDim, NumAllChan,lmax
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
    !PointsPerBox=9  ! careful! if you change this you need to go change the weights array.
    !lmax = 2*NumChannels
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
    onesixth = onethird*0.5d0
    
  END SUBROUTINE ReadGlobal

  SUBROUTINE  CalcNumChannels(m,l)
    IMPLICIT NONE
    INTEGER, INTENT(IN):: l,m
    IF(MOD(m,2).EQ.1)THEN
       NumChannels = (lmax - m + 1)/2
    ELSE
       IF(MOD(l,2).EQ.0)THEN
          NumChannels = (lmax -m)/2 + 1
       ELSE
          NumChannels = (lmax - m)/2
       END IF
    END IF
  END SUBROUTINE CalcNumChannels
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

    !    weights(0)=1d0
    !write(6,*) "check that PointsPerBox is even"
    !if (mod(PointsPerBox,2).ne.0) then
       !write(6,*) "PointsPerBox is not even, adding one point"
       !PointsPerBox = PointsPerBox + 1
       !stop
    !endif
    do i = 1, PointsPerBox-1
       if (mod(i,2).eq.1)  weights(i)=4d0
       if (mod(i,2).eq.0)  weights(i)=2d0
    enddo
    weights(PointsPerBox)=1d0

    identity = 0d0
    yprevious = 0d0
    ycurrent = 0d0
    ystart = 0d0
    do i=1,NumChannels
      identity(i,i)=1d0
      ystart(i,i)=1d20
    enddo

  end subroutine initprop

  subroutine boxstep(xx,yi,yf,Pot,iBox,NumBoxes)
    implicit none
    integer i,j,step,iBox,NumBoxes
    double precision h
    double precision xx(0:PointsPerBox)
    double precision yi(NumChannels,NumChannels),yf(NumChannels,NumChannels)
    double precision Pot(NumChannels,NumChannels,PointsPerBox) !make sure Pot includes the threshold offsets
    double precision tempy(NumChannels,NumChannels)
    double precision vtemp1(NumChannels,NumChannels), vtemp2(NumChannels,NumChannels), un(NumChannels,NumChannels)

    h=xx(2)-xx(1)

    yprevious = yi
    weights(PointsPerBox)=2d0
    if(iBox.eq.NumBoxes) weights(PointsPerBox)=1d0
    
    do step = 1, PointsPerBox
       !write(6,*) weights(step)
       vtemp1 = 2d0*mu*(identity*Energy-Pot(:,:,step))
       if (mod(step,2).eq.0) then
          un = vtemp1
       else
          vtemp2 = identity + h*h*onesixth*vtemp1
          call sqrmatinv(vtemp2,NumChannels)
          un = matmul(vtemp2,vtemp1)
       endif
       
       tempy = identity + h*yprevious
       call sqrmatinv(tempy,NumChannels)
       ycurrent = MATMUL(tempy,yprevious) - onethird*h*weights(step)*un
       yprevious = ycurrent

    enddo
    !stop
    yf = ycurrent
    
  end subroutine boxstep
 
  SUBROUTINE DeallocateProp
    IMPLICIT NONE
    DEALLOCATE(VV,u,ycurrent,yprevious,ystart,identity,weights)
  END SUBROUTINE DeallocateProp
    
  !****************************************************************************************************

end module logderprop
!============================================================================================

module scattering
  use datastructures
  !****************************************************************************************************

CONTAINS

 SUBROUTINE CalcK(Y,rm,SD,mu,d,alpha,EE,Eth,NumChannels,NumOpen,lam)
   use DipoleDipole
   IMPLICIT NONE
   TYPE(ScatData) :: SD

   DOUBLE PRECISION mu, EE, rm, d, alpha,Y(NumChannels,NumChannels)
   DOUBLE PRECISION, ALLOCATABLE :: JJ(:),NN(:),JJp(:),NNp(:)
   double precision, allocatable :: Ktemp1(:,:),Ktemp2(:,:)
   DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp,Pi,rhypi,rhypk,rhypip,rhypkp,ldrhk,ldrhi
   DOUBLE PRECISION k(NumChannels),Eth(NumChannels),lam(NumChannels)
   double precision rj, drj, ry, dry
   double precision ri, dri, rk, drk, ldi, ldk
   complex*16, allocatable :: tmp(:,:),Identity(:,:)
   complex*16  II
   INTEGER i,j,NumOpen,Numchannels,no, nw, nc, beta,itest

   
   II=(0d0,1d0)
   Pi=dacos(-1d0)

   no=0
   nw=0
   nc=0
   

   DO i = 1,NumChannels
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
   IF((no+nw+nc).NE.NumChannels) THEN
      WRITE(6,*) "Channel miscount in calcK"
      STOP
   ENDIF
!   write(6,*) "no = ", no
   deallocate(SD%S,SD%T,SD%sigma)
   
   allocate(SD%S(no,no),SD%T(no,no),SD%sigma(no,no))
   ALLOCATE(JJ(NumChannels),NN(NumChannels),tmp(no,no))
   ALLOCATE(JJp(NumChannels),NNp(NumChannels))
   allocate(Ktemp1(NumChannels,NumChannels))
   allocate(Ktemp2(NumChannels,NumChannels))
   allocate(Identity(NumChannels,NumChannels))
   Identity = 0d0;

   DO i = 1,no
      Identity(i,i) = 1d0
      !write(6,*) k(i), rm
!      CALL hyperrjry(INT(d),alpha,lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
      call fdfgdg(lam(i),k(i),rm,rj,drj,ry,dry)
      JJ(i) = rj
      NN(i) = ry
      JJP(i) = drj
      NNP(i) = dry
!!$      JJ(i) = rhypj/dsqrt(Pi*k(i))
!!$      NN(i) = -rhypy/dsqrt(Pi*k(i))
!!$      JJp(i) = dsqrt(k(i)/Pi)*rhypjp
!!$      NNp(i) = -dsqrt(k(i)/Pi)*rhypyp
   ENDDO
   do i=no+1,NumChannels  !These should never be getting called for this dipole dipole calculation anyhow.  Be careful with other calculations.
      Identity(i,i) = 1d0
!!$      CALL hyperrirk(INT(d),alpha,lam(i),k(i)*rm,rhypi,rhypk,rhypip,rhypkp,ldrhi,ldrhk)
!!$      JJ(i) = 1d0
!!$      NN(i) = -1d0
!!$      JJp(i) = ldrhi
!!$      NNp(i) = ldrhk
      call bfdfgdg(lam(i),k(i),rm,k(i),ri,dri,rk,drk,ldi,ldk) ! really x*sphbes
      JJ(i) = 1d0
      NN(i) = 1d0
      JJP(i) = ldi
      NNP(i) = ldk
   ENDDO

   Ktemp1=0d0
   Ktemp2=0d0
   SD%K=0d0
   do i=1,NumChannels
      Ktemp1(i,i) = -NNp(i)
      Ktemp2(i,i) = JJp(i)
      do j=1,NumChannels
         Ktemp1(i,j) = Ktemp1(i,j) + Y(i,j)*NN(j)
         Ktemp2(i,j) = Ktemp2(i,j) - Y(i,j)*JJ(j)
      enddo
   enddo
   call sqrmatinv(Ktemp1,NumChannels)
   SD%K = MATMUL(Ktemp1,Ktemp2)

   tmp = Identity(1:no,1:no) - II*SD%K(1:no,1:no)
 
   SD%S = Identity(1:no,1:no) + II*SD%K(1:no,1:no)
   call CompSqrMatInv(tmp,no)
   SD%S = MATMUL(SD%S,tmp)
   SD%T = -II*0.5d0*(SD%S-Identity(1:no,1:no))
   SD%sigma = conjg(SD%T)*SD%T*Pi/(2d0*mu*EE)

   SD%tandel = SD%K(1,1)
   SD%delta = atan(SD%tandel)
   SD%sindel = sin(SD%delta)
   SD%sin2del = SD%sindel**2

   DEALLOCATE(JJ, NN, JJp, NNp, Ktemp1, Ktemp2, tmp, Identity)

 END SUBROUTINE CalcK

end module scattering
!=========================================================================================
program main

  use GlobalVars
  use DipoleDipole
  use logderprop
  use scattering
  implicit none
  type(DPData) DP
  type(ScatData) SD
  double precision, allocatable :: VPot(:,:,:),sigmagrandtotal(:),sigmatot(:,:)
  double precision, allocatable :: BoxGrid(:)
  double precision, allocatable :: x(:),yin(:,:),yout(:,:),Egrid(:)
  double precision time1, time2,rj,drj,ry,dry,Pi
  integer iBox,ml,iE,m,itest
  Pi = dacos(-1d0)
  call Setupam
  InputFile = 'logder.inp'
  call ReadGlobal()
  NumChannels = 0
  
!!$  m=0
!!$  call CalcNumChannels(0,1)
!!$  DP%lmax = lmax
!!$  call AllocateDP(DP,NumChannels)
!!$  call MakeDipoleDipoleCouplingMatrix(DP) 
!!$
!!$!stop 
!!$  allocate(VPot(NumChannels,NumChannels,0:PointsPerBox))
!!$  allocate(x(0:PointsPerBox))
!!$  allocate(BoxGrid(NumBoxes+1))
!!$  call GridMaker(BoxGrid,NumBoxes+1,xStart,xEnd,"quadratic")
!!$ 
!!$  DO iBox=1,NumBoxes
!!$     VPot=0d0
!!$     x=0d0
!!$     call SetDipoleDipolePot(VPot,DP,PointsPerBox,x,BoxGrid(iBox),BoxGrid(iBox+1),NumChannels,m,lmax)
!!$     !call PlotPot(Vpot,x,NumChannels,PointsPerBox,2,2,15)
!!$  END DO
!!$ 
!!$ stop 

!  call testbessik
!  stop
 allocate(sigmatot(0:lmax,NumE))
 allocate(Egrid(NumE))
 allocate(sigmagrandtotal(NumE))
 allocate(x(0:PointsPerBox))
 allocate(BoxGrid(NumBoxes+1))
 call GridMaker(Egrid,NumE,Emin,Emax,"log")
 call GridMaker(BoxGrid,NumBoxes+1,xStart,xEnd,"linear")
! call GridMaker(BoxGrid,NumBoxes+1,xStart,xEnd,"quadratic")
! call printmatrix(BoxGrid,NumBoxes+1,1,6)
! stop
 write(6,*) "lmax = ", lmax
 Do m=0,lmax-1
    !m = 0
    write(6,*) "doing calculation for m = ",m
      call CalcNumChannels(m,1)
      !write(6,*) NumChannels, m
      DP%lmax = lmax 
      call AllocateDP(DP,NumChannels)
      call MakeDipoleDipoleCouplingMatrix(DP)
      allocate(VPot(NumChannels,NumChannels,0:PointsPerBox))
      allocate(yin(NumChannels,NumChannels),yout(NumChannels,NumChannels))
      call AllocateScat(SD,NumChannels)

     
      call initprop ! sets the weights and the initial Y matrix.
      ! call cpu_time(time1)
      !  write(6,"(3A15)") "energy","sigma","time")
      DO iE = 1,NumE
         Energy = Egrid(iE)
         yin = ystart
         DO iBox=1,NumBoxes
            !write(6,*) "calling set morse"
            VPot = 0d0
            x=0d0
            
            call SetDipoleDipolePot(VPot,DP,PointsPerBox,x,BoxGrid(iBox),BoxGrid(iBox+1),NumChannels,m,lmax)
            
            call boxstep(x,yin,yout,VPot(:,:,1:PointsPerBox),iBox,NumBoxes)
            yin = yout
         END DO
         
         call CalcK(yout,BoxGrid(NumBoxes+1),SD,mu,EffDim,AlphaFactor,Energy,DP%Eth,NumChannels,NumChannels,DP%lam)

         write(6,*) "K:"
         call printmatrix(SD%K,NumChannels,NumChannels,6)
         write(6,*) "Y:"
         call printmatrix(yout,NumChannels,NumChannels,6)
         itest=1
         write(6,*) INT(EFFDIM),ALPHAFACTOR,DP%lam(itest),dsqrt(2d0*mu*Energy)*BoxGrid(NumBoxes+1)
         write(6,*)
         !CALL hyperrjry(INT(EFFDIM),ALPHAFACTOR,DP%lam(itest),dsqrt(2d0*mu*Energy)*BoxGrid(NumBoxes+1),rj,ry,drj,dry)
         !rj = rj/dsqrt(Pi*dsqrt(2d0*mu*Energy))
         !drj = drj*dsqrt(dsqrt(2d0*mu*Energy)/Pi)
         call fdfgdg(INT(DP%lam(itest)),dsqrt(2d0*mu*Energy),BoxGrid(NumBoxes+1),rj,drj,ry,dry)
         write(6,*) "looking at diagonal element:",itest
         write(6,*) "Energy, bessel functions, log-der:", Energy, BoxGrid(NumBoxes+1), drj/rj
         write(6,*) "error in log-der of channel", drj/rj - yout(itest,itest)
         write(6,*) "... corresponding T matrix: ", SD%T(itest,itest)
         
         stop
!!$           call printmatrix(SD%T,NumChannels,NumChannels,6)
         sigmatot(m,iE) = sum(SD%sigma)
!!$           write(10,*) Energy, sigmatot(m,iE)
!!$           WRITE(6,*)  Energy, sigmatot(m,iE)
      END DO
      
      deallocate(VPot)
      deallocate(yin)
      deallocate(yout)
      call DeallocateDP(DP)
      call DeAllocateScat(SD)
      Call DeallocateProp
      
      
   END DO
   
   sigmagrandtotal = 0
   DO iE=1,NumE
      Energy = Egrid(iE) 
      DO m = 1,lmax-1
         sigmagrandtotal(iE)=sigmagrandtotal(iE)+sigmatot(m,iE)
      END DO
    sigmagrandtotal(iE)=2*sigmagrandtotal(iE)
    sigmagrandtotal(iE)=sigmagrandtotal(iE)+sigmatot(0,iE)
   Write(6,*) Energy, sigmagrandtotal(iE)
   WRITE(20,*) Energy, sigmagrandtotal(iE)
 END DO

  
  !call cpu_time(time2)
  !write(6,*) "total time for calculation = ", time2-time1
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
!=========================================================================================
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)

  DO j = 1,nr
     WRITE(file,30) (M(j,k), k = 1,nc)
  ENDDO

20 FORMAT(1P,100D20.12)
30 FORMAT(1p,100d15.7)
END SUBROUTINE printmatrix
!=========================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fdfgdg(L,k,r,f,df,g,dg) ! returns the oscillatory Riccati functions and their radial derivatives.
!  use modb
  implicit none
  integer L
  double precision k, r, f, df, g, dg, nu, x, factor
  double precision J, dJ, RTPIO2, Y, dY
  PARAMETER (RTPIO2=1.25331413731550d0)

  x = k*r
  if(L.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesjy'

  nu=L+0.5d0
!  write(6,*) "k, r:", k, r
  call bessjy(x,nu,J,Y,dJ,dY)
  factor=RTPIO2*sqrt(x)
  f=factor*J
  g=factor*Y
  df=k*factor*(dJ+J/(2.d0*x))
  dg=k*factor*(dY+Y/(2.d0*x))
  
end subroutine fdfgdg
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine bfdfgdg(L,k,r,xscale,f,df,g,dg,ldi,ldk) ! returns the exponential-type Riccati functions and their radial derivatives.
!  use modb
  implicit none
  integer L
  double precision k, r, f, df, g, dg, nu,xscale
  double precision RTPIO2, RT2OPI, factor, x
  double precision alpha, beta, alphap, betap, ldi, ldk
  PARAMETER (RTPIO2=1.25331413731550d0) 
  PARAMETER (RT2OPI=0.7978845608028654d0)
      x=k*r
      if(L.lt.0.d0.or.x.le.0.d0) then
         write(6,*) 'bad arguments in sphbesik, (L, x) = ', L, x
      endif
  nu = L + 0.5d0
  call MyScaledBessIK(x, nu, alpha, beta, alphap, betap, ldi,ldk)
  factor=RT2OPI*sqrt(x)

  f = factor*exp(x-xscale)*alpha
  df = k*factor*exp(x-xscale)*( alpha*(1.d0 + 1.d0/(2d0*x)) + alphap )

  g = factor*exp(xscale-x)*beta
  dg = k*factor*exp(xscale-x)*( beta*(-1.d0 + 1.d0/(2.d0*x)) + betap )
  !  ldi = (alphap/alpha + 1.d0 + 1.d0/(2.d0*x))*k
  !  ldk = (betap/beta - 1.d0 + 1.d0/(2.d0*x))*k

  ldi = df/f 
  ldk = dg/g
  
end subroutine bfdfgdg
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine testbessik
  implicit none
  real*8 ri,rk,dri,drk, x, dx, ldi, ldk,xscale,nu0, r, k
  real*8 si, sk, sip, skp, alpha, beta, alphap, betap
  real*8 v, vm,pi,xinit,xfinal,f,df,g,dg,rj,ry,rjp,ryp
  complex*16, allocatable :: YCI(:), YCK(:), YDI(:), YDK(:)
!  complex*16 z
  integer ixtot, L, ix
  pi = acos(-1.d0)
  !write(6,*) 'pi = ',pi

  L=1
  nu0 = 0.5d0
  xinit=0.001d0
  xfinal=100.d0
  ixtot=1000
  dx=(xfinal-xinit)/dble(ixtot-1)

  xscale = 0.d0
  
!  write(6,*) dx, ixtot
!  read(*,*)
  !  x=0.0d0
  k = 1.d0

  do ix = 1, ixtot
     x = xinit + (ix-1)*dx
     r = x/k
!     CALL hyperrjry(3,1,1,k*r,rj,ry,rjp,ryp)
!     xscale = x
     !     call Mysphbesik(L,x,xscale,si,sk,sip,skp,ldk,ldi) ! really x*sphbes
!     call bfdfgdg(L,k,r,xscale,ri,dri,rk,drk,ldi,ldk) ! really x*sphbes
     call fdfgdg(L,k,r,f,df,g,dg) ! returns the oscillatory Riccati functions and their radial derivatives.
!     call sphbesik(L,x,si,sk,sip,skp) ! really x*sphbes
     !     write(6,10) x, ri, dri, rk, drk, ldi, ldk, f, df, g, dg
     write(30,10) r, f, df, g, dg
!     write(30,10) r, rj,ry,rjp,ryp
     
  end do

10 format(100e24.14)
end subroutine testbessik
