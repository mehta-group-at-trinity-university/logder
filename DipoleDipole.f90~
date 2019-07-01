MODULE DipoleDipole 
  TYPE DPData
     INTEGER lmax,ml
     DOUBLE PRECISION, ALLOCATABLE :: cllp(:,:,:)
     DOUBLE PRECISION, ALLOCATABLE :: Eth(:)
     DOUBLE PRECISION, ALLOCATABLE :: lam(:)
     LOGICAL even
  END TYPE DPData

CONTAINS
  SUBROUTINE AllocateDP(DP,N)
    IMPLICIT NONE
    integer N ! Numchannels
    TYPE(DPData) DP
    ALLOCATE(DP%cllp(0:DP%lmax,0:DP%lmax,0:DP%lmax))
    ALLOCATE(DP%Eth(N))
    ALLOCATE(DP%lam(N))

  END SUBROUTINE AllocateDP

  SUBROUTINE MakeDipoleDipoleCouplingMatrix(DP)
    IMPLICIT NONE
    TYPE(DPData) :: DP
    INTEGER l,lp,m
    DOUBLE PRECISION, EXTERNAL :: THRJ
    DOUBLE PRECISION tj1,tj2,phase,prefact

    DP%Eth=0d0
    DP%cllp=0d0
    DO l=0,DP%lmax
       DO lp=Max(l-2,0),Min(DP%lmax,l+2)
          prefact=DBLE((2*l+1)*(2*lp+1))
          prefact=dsqrt(prefact)
          DO m = 0,l
             phase=(-1)**m
             tj1 = THRJ(2*l,2*2,2*lp,-2*m,0,2*m)
             tj2 = THRJ(2*l,2*2,2*lp,0,0,0)
             DP%cllp(l,lp,m)=prefact*phase*tj1*tj2
             DP%cllp(lp,l,m)=DP%cllp(l,lp,m)
             
            !Write(1000,*) l, lp, m,  DP%cllp(l,lp,m)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE MakeDipoleDipoleCouplingMatrix
    

  SUBROUTINE DeallocateDP(DP)
    IMPLICIT NONE
    TYPE(DPData) DP
    DEALLOCATE (DP%Eth,DP%cllp,DP%lam)
  END SUBROUTINE DeallocateDP


  SUBROUTINE SetDipoleDipolePot(VPot,DP,NPTS,x,x1,x2,Nchan,m,lmax)
  
    IMPLICIT NONE
    TYPE(DPData) DP
    INTEGER, INTENT(IN) :: m,lmax, Nchan,NPTS
    INTEGER lp,l1
    INTEGER kx, mch, nch,j
    DOUBLE PRECISION x(0:NPTS), x1, x2, dx, VPot(Nchan,Nchan,0:NPTS)

    lp=0
    l1=0
    VPot=0d0
    x=0d0
    dx = (x2-x1)/DBLE(NPTS)
    If(mod(m,2).eq.0)then 
       Do kx=0,NPTS
          x(kx) = x1 + kx*dx
          mch = 1
          Do l1 = m+1,lmax-1, 2
           nch = 1
           DP%lam(mch) = l1
           Do lp = m+1,lmax-1, 2
             If(mch.EQ.nch)THEN
                VPot(mch,nch,kx)=-2d0*DP%cllp(l1,lp,m)/x(kx)**3 + 0.50*l1*(l1+1)/x(kx)**2
             Else 
                VPot(mch, nch, kx) = -2d0*DP%cllp(l1,lp,m)/x(kx)**3
                VPot(nch,mch,kx) = VPot(mch,nch,kx)
             End if
            ! Write(25,*)  l1,lp,m, Vpot(mch,nch,kx)
             nch = nch + 1
          End do
          mch = mch + 1
       End do
    End do
 Else
    Do kx=0,NPTS
          x(kx) = x1 + kx*dx
          mch = 1
          Do l1 = m,lmax-1, 2
           nch = 1
           DP%lam(mch) = l1
           Do lp = m,lmax-1, 2
             If(mch.EQ.nch)THEN
                VPot(mch,nch,kx)=-2d0*DP%cllp(l1,lp,m)/x(kx)**3 + 0.50*l1*(l1+1)/x(kx)**2
             Else 
                VPot(mch, nch, kx) = -2d0*DP%cllp(l1,lp,m)/x(kx)**3
                VPot(nch,mch,kx) = VPot(mch,nch,kx)
             End if
            ! Write(25,*)  l1,lp,m, Vpot(mch,nch,kx)
             nch = nch + 1
          End do
          mch = mch + 1
       End do
    End do
 end if

!!$    lp=0
!!$    l1=0
!!$    VPot=0d0
!!$    x=0d0
!!$    dx = (x2-x1)/DBLE(NPTS)
!!$    Do kx=0,NPTS
!!$       x(kx) = x1 + kx*dx
!!$       mch = 1
!!$       Do l1 = m,lmax, 2
!!$           nch = 1
!!$           Do lp = m,lmax, 2
!!$             If(mch.EQ.nch)THEN
!!$                VPot(mch,nch,kx)=-2d0*DP%cllp(l1,lp,m)/x(kx)**3 + 0.50*l1*(l1+1)/x(kx)**2
!!$             Else 
!!$                VPot(mch, nch, kx) = -2d0*DP%cllp(l1,lp,m)/x(kx)**3
!!$                VPot(nch,mch,kx) = VPot(mch,nch,kx)
!!$             End if
!!$             Write(25,*)  l1,lp,m, Vpot(mch,nch,kx)
!!$             nch = nch + 1
!!$          End do
!!$          mch = mch + 1
!!$       End do
!!$    End do



  END SUBROUTINE SetDipoleDipolePot

  SUBROUTINE PlotPot(VPot,x,Nchan,NPTS,i,j,file)
    IMPLICIT NONE
    INTEGER i,j,NPTS,Nchan,file,xi
    DOUBLE PRECISION VPot(Nchan,Nchan,0:NPTS),x(0:NPTS)
    
    Do xi = 0, NPTS
       Write(file,*)  x(xi), VPot(i,j,xi)
    End Do 

  END SUBROUTINE PlotPot
       

    


END MODULE DipoleDipole
