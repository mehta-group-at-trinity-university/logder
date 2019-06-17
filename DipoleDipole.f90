MODULE DipoleDipole
  TYPE DPData
     DOUBLE PRECISION, ALLOCATABLE :: Eth(:)
  END TYPE DPData

CONTAINS
  SUBROUTINE AllocateDP(DP,N)
    IMPLICIT NONE
    integer N ! Numchannels
    TYPE(DPData) DP
    ALLOCATE(DP%Eth(N))
  END SUBROUTINE AllocateDP

  FUNCTION  Coupling(m,l,lp)
    IMPLICIT NONE
    INTEGER l,lp,m 
    DOUBLE PRECISION, EXTERNAL :: THRJ
    DOUBLE PRECISION tj1,tj2,phase,prefact
    DOUBLE PRECISION Coupling
    
    phase= (-1)**m
    prefact= DBLE((2*l+1)*(2*lp+1))
    prefact= dsqrt(prefact)
    tj1 = THRJ(2*l,2*2,2*lp,-2*m,0,2*m)
    tj2 = THRJ(2*l,2*2,2*lp,0,0,0)
    Coupling = phase*prefact*tj1*tj2

  END FUNCTION Coupling
  
  SUBROUTINE SetDipoleDipolePot(VPot,DP,NPTS,x,x1,x2,Nchan,m,lmax)
  
    IMPLICIT NONE
    TYPE(DPData) DP
    INTEGER, INTENT(IN) :: m,lmax,Nchan,NPTS
    INTEGER lp,l1
    INTEGER kx, mch, nch
    DOUBLE PRECISION x(0:NPTS), x1, x2, dx, VPot(Nchan,Nchan,0:NPTS)

    lp=0
    l1=0
    VPot=0d0
    x=0d0
    dx = (x2-x1)/DBLE(NPTS)
    
    Do kx=0,NPTS
       x(kx) = x1 + kx*dx
       mch = 1
       Do l1 = m, lmax, 2
          nch = 1
          Do lp = m, lmax, 2
             If(mch.EQ.nch)THEN
                VPot(mch,nch,kx)=-2d0*Coupling(m,l1,lp)/x(kx)**3 + 0.50*l1*(l1+1)/x(kx)**2
             Else 
             VPot(mch, nch, kx) = -2d0*Coupling(m,l1,lp)/x(kx)**3
             VPot(nch,mch,kx )= VPot(mch,nch,kx)
             End if
             nch = nch + 1
          End do
          mch = mch + 1
       End do
    End do

  END SUBROUTINE SetDipoleDipolePot

  SUBROUTINE PlotPot(VPot,x,Nchan,NPTS,i,j,file)
    IMPLICIT NONE
    INTEGER i,j,NPTS,Nchan,file,xi
    DOUBLE PRECISION VPot(Nchan,Nchan,0:NPTS),x(0:NPTS)
    
    Do xi = 0, NPTS
       Do i = 1,Nchan
          Do j = 1,Nchan
             Write(file,*) VPot(i,j,x(xi))
          End Do
        End Do
    End Do 

  END SUBROUTINE PlotPot
       

    


END MODULE DipoleDipole
