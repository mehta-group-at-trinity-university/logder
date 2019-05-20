!****************************************************************************************************
    MODULE MorsePotential
    USE DataStructures
!    USE Quadrature
    IMPLICIT NONE
    TYPE Morse
       DOUBLE PRECISION D(3),V(3,3),rc(3,3),a(3,3)
       DOUBLE PRECISION Eth(3)
       INTEGER NumMorseOpen,l
    END TYPE Morse

  CONTAINS
    SUBROUTINE InitMorse(M)
      TYPE(Morse) :: M
      M%D(1)=1.d0
      M%D(2)=4.d0
      M%D(3)=6.d0
      M%V=0d0
      M%V(1,2)=0.2d0
      M%V(2,1)=M%V(1,2)
      M%V(2,3)=0.2d0
      M%V(3,2)=M%V(2,3)
      M%a=1d0
      M%a(1,2)=0.5d0
      M%a(2,1)=M%a(1,2)
      M%a(2,3)=0.5d0
      M%a(3,2)=M%a(2,3)
      M%rc=1d0                !initialize coupling ranges to 1
      M%rc(1,2) = 0.5d0       !set any coupling ranges different from 1d0
      M%rc(2,1) = M%rc(1,2)
      M%l=0
      M%Eth(1) = 0d0
      M%Eth(2) = 3d0
      M%Eth(3) = 5d0

    END SUBROUTINE InitMorse

    FUNCTION Morse1(a,D,re,r)
      DOUBLE PRECISION Morse1, a, D, re, r
      Morse1 = D*(1d0 - dexp(-(r - re)/a))**2 - D
      RETURN
    END FUNCTION Morse1

    SUBROUTINE SetMorsePotential(VPot,M,NPTS,x,x1,x2)

      IMPLICIT NONE
      TYPE(Morse), INTENT(in) :: M
      INTEGER kx,mch,nch,NChan,NPTS
      DOUBLE PRECISION pot(3,3),VPot(3,3,0:NPTS)
      DOUBLE PRECISION x1,x2,dx,x(0:NPTS)
!      double precision, dimension(:), allocatable :: x(:)
!      double precision, dimension(:), allocatable :: VPot(:,:,:)

      Nchan=3
      VPot(:,:,:) = 0d0
      x=0d0
      dx = (x2-x1)/dble(NPTS)
      
      DO kx = 0,NPTS
         x(kx) = x1 + kx*dx
         
         DO mch = 1,NChan
            VPot(mch,mch,kx) = Morse1(M%a(mch,mch),M%D(mch),M%rc(mch,mch),x(kx)) + M%Eth(mch)  ! Diagonal Morse potential
            DO nch = 1, mch-1
               VPot(mch,nch,kx) = M%V(mch,nch)*dexp(-(x(kx)-M%rc(mch,nch))**2/M%a(mch,nch)**2)  ! Gaussian off-diagonal couplings
               VPot(nch,mch,kx) = VPot(mch,nch,kx)
            ENDDO
         ENDDO
      ENDDO
      
   
    END SUBROUTINE SetMorsePotential
  END MODULE MorsePotential
