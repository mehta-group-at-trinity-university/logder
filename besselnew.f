c!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! this subroutine calculates the hyerpspehrical besselfunctions multiplied by x**alpha
c!     rj = (x**alpha*hypj)
c!     rjp = alpha*x**(alpha-1)*hypj + x**alpha*hypjp
c!
      subroutine hyperrjry(d,alpha,lam,x,rj,ry,rjp,ryp)
      implicit none
      integer d,df
      double precision, external :: mygamma
      double precision lam, x, j, y, jp, yp,gam,halfd,alpha
      double precision hypj,hypy,hypjp,hypyp
      double precision rj,rjp, ry, ryp
      double precision order,prefact
      halfd=0.5d0*dble(d)
      order = halfd + lam - 1d0
      call doubfact(d-4,df)
!write(6,*) "order=",order, halfd, lam, d

      if (x.le.10d0*order**2) then
         call bessjy(x,order,j,y,jp,yp)
      else
         call BesselAsymNew(order,x,j,y,jp,yp)
      endif
      
      !call bessjy(x,order,j,y,jp,yp)
      prefact = mygamma(halfd-1.d0)*2**(halfd-2d0)/df
      hypj = prefact*j*x**(-halfd+1d0)
      hypy = prefact*y*x**(-halfd+1d0)
      hypjp = prefact*x**(-halfd+1d0)*(jp - (halfd - 1d0)*j/x)
      hypyp = prefact*x**(-halfd+1d0)*(yp - (halfd - 1d0)*y/x)

      rj = x**alpha*hypj
      rjp = x**alpha*(alpha*hypj/x + hypjp)
      ry = x**alpha*hypy
      ryp = x**alpha*(alpha*hypy/x + hypyp)

      end
c!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! This is a subroutine that returns the hyperspherical bessel functions as defined in Avery
      subroutine hyperjy(d,lam,x,hypj,hypy,hypjp,hypyp)
      implicit none
      integer d,df
      double precision, external :: mygamma
      double precision lam, x, j, y, jp, yp,gam,halfd
      double precision hypj,hypy,hypjp,hypyp
      double precision order,prefact
      halfd=0.5d0*d
      order = halfd + lam - 1d0
      call doubfact(d-4,df)
      call bessjy(x,order,j,y,jp,yp)
      prefact = mygamma(halfd-1.d0)*2**(halfd-2d0)/df
      hypj = prefact*j*x**(-halfd+1d0)
      hypy = prefact*y*x**(-halfd+1d0)
      hypjp = prefact*x**(-halfd+1d0)*(jp - (halfd - 1)*j/x)
      hypyp = prefact*x**(-halfd+1d0)*(yp - (halfd - 1)*y/x)
      end
c!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! This is a subroutine that returns the hyperspherical bessel functions as defined in Avery
      subroutine hyperrirk(d,alpha,lam,x,rhypi,rhypk,rhypip,rhypkp,ldrhi,ldrhk)
      implicit none
      integer d,df
      double precision, external :: mygamma
      double precision lam, x, i, k, ip, kp,gam,halfd,alpha
      double precision hypi,hypk,hypip,hypkp
      double precision rhypi,rhypk,rhypip,rhypkp,ldrhi,ldrhk
      double precision order,prefact
      double precision ai,bk,aip,bkp,ldi,ldk
      integer dfact(12)
      dfact = (/1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080/)
      halfd=0.5d0*d
      order = halfd + lam - 1d0
      df=1
      if((d-4).gt.1) df = dfact(d-4)

c!     Inu(x) = exp(x) * alpha(x)
c!     Knu(x) = exp(-x) * beta(x)
c!     This routine returns alpha, beta, alpha', beta', I'/I, and K'/K
      call MyScaledBessIK(x, order, ai, bk, aip, bkp, ldi,ldk)

      prefact = mygamma(halfd-1.d0)*2**(halfd-2d0)/df
      hypi = prefact*dexp(x)*ai*x**(-halfd+1d0)
      hypk = prefact*dexp(-x)*bk*x**(-halfd+1d0)
      hypip = prefact*0.5d0*dexp(x)*x**(-halfd)*(2d0*x*aip + ai*(-d+2*x+2))
      hypkp = prefact*0.5d0*dexp(-x)*x**(-halfd)*(2d0*x*bkp - bk*(d+2*x-2))

      rhypi = x**alpha*hypi
      rhypip = x**alpha*(alpha*hypi/x + hypip)
      rhypk = x**alpha*hypk
      rhypkp = x**alpha*(alpha*hypk/x + hypkp)

      ldrhi=dble(6-d)/(2d0*x) + ldi
      ldrhk=dble(6-d)/(2d0*x) + ldk
      end
c!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c! This is a subroutine that returns the modified hyperspherical bessel functions.
      subroutine hyperik(d,lam,x,hypi,hypk,hypip,hypkp)
      implicit none
      integer d,df
      double precision, external :: mygamma
      double precision lam, x, i, k, ip, kp,gam,halfd
      double precision hypi,hypk,hypip,hypkp
      double precision order,prefact
      double precision ai,bk,aip,bkp,ldi,ldk
      integer dfact(12)
      dfact = (/1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080/)
      halfd=0.5d0*d
      order = halfd + lam - 1d0 ! same as l+1/2 for d=3
      df = 1
      if((d-4).gt.1) df = dfact(d-4)
   
c!     Inu(x) = exp(x) * alpha(x)
c!     Knu(x) = exp(-x) * beta(x)
c!     This routine returns alpha, beta, alpha', beta', I'/I, and K'/K
      call MyScaledBessIK(x, order, ai, bk, aip, bkp, ldi,ldk)

      prefact = mygamma(halfd-1.d0)*2**(halfd-2d0)/df
      hypi = prefact * dexp(x)*ai * x**(-halfd+1d0)
      hypk = prefact * dexp(-x)*bk * x**(-halfd+1d0)
      hypip = prefact*0.5d0*dexp(x)*x**(-halfd)* (2d0*x*aip + ai*(-d+2*x+2))
      hypkp = prefact*0.5d0*dexp(-x)*x**(-halfd)*(2d0*x*bkp - bk*(d+2*x-2))

      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine doubfact(n,df)
      implicit none
      integer n,m,df

      df=1
      m=n

      do while (m.gt.1)
         df=df*m
         write(6,*) m, df
         m=m-2
      end do

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sphbesjy(n,x,sj,sy,sjp,syp) ! change norm
!NB:  This is NOT the spherical bessel function jn(x), instead it is x*jn(x), the Riccati function
!leaving a factor of sqrt(x) in the numerator of the prefactor
      INTEGER n
      DOUBLE PRECISION sj,sjp,sy,syp,x
!U    USES bessjy
      DOUBLE PRECISION factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.25331413731550d0)
      if(n.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesjy'
      order=n+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2*sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp+sj/(2.d0*x)
      syp=factor*ryp+sy/(2.d0*x)
      return
      END SUBROUTINE sphbesjy
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sphbesik(n,x,si,sk,sip,skp) ! change norm
!NB:  This is NOT the spherical bessel function kn(x), instead it is x*kn(x),
!     leaving a factor of sqrt(x) in the numerator of the prefactor
      INTEGER n
      DOUBLE PRECISION si,sip,sk,skp,x,ri,rk,rip,rkp,ldi,ldk
      DOUBLE PRECISION factor,order,RTPIO2
      PARAMETER (RTPIO2=1.25331413731550d0)
      PARAMETER (RT2OPI=0.7978845608028654d0)
      if(n.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesik'
      order=n+0.5d0
      call bessik(x,order,ri,rk,rip,rkp)
      factor=RT2OPI*sqrt(x)
      si=factor*ri
      sk=factor*rk
      sip=factor*rip+si/(2.d0*x)
      skp=factor*rkp+sk/(2.d0*x)
      return
      END SUBROUTINE sphbesik
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Mysphbesik(n,x,xscale,si,sk,sip,skp,ldk,ldi) ! change norm
      implicit none
!     NB:  This is NOT the spherical bessel function kn(x), instead it is x*kn(x),
!     leaving a factor of sqrt(x) in the numerator of the prefactor to the cylindrical Kn(x).
!     inu and knu are multiplied by exp(-xscale) and exp(xscale), respectively.
!     xscale is a constant number that must be chosen at the beginning of a scattering calcualtion.
      INTEGER n
      DOUBLE PRECISION si,sip,sk,skp,x,ri,rk,rip,rkp,ldi,ldk,xscale
      DOUBLE PRECISION factor,order,RTPIO2
      double precision RT2OPI,bigk, bigi, bigkp, bigip
      PARAMETER (RTPIO2=1.25331413731550d0)
      PARAMETER (RT2OPI=0.7978845608028654d0)
      if(n.lt.0.d0.or.x.le.0.d0) write(6,*) 'bad arguments in sphbesik, (n, x) = ', n, x
      order=n+0.5d0

c$$$  call bessik(x,order,ri,rk,rip,rkp)
c$$$  factor=RT2OPI*sqrt(x)
c$$$  si=factor*ri
c$$$  sk=factor*rk
c$$$  sjp=factor*rip+si/(2.d0*x)
c$$$  syp=factor*rkp+sk/(2.d0*x)

      call MyScaledBessIK(x, order, ri, rk, rip, rkp, ldi,ldk)
!     Inu(x) = exp(x) * alpha(x)
!     Knu(x) = exp(-x) * beta(x)
!     This routine returns alpha, beta, alpha', beta', I'/I, and K'/K
      factor=RT2OPI*sqrt(x)
      si = ri*factor*exp(x-xscale)
      sk = rk*factor*exp(xscale-x)
      bigkp = exp(xscale-x)*(rkp-rk)
      bigip = exp(x-xscale)*(rip+ri)
      sip = factor*bigip + si/(2.d0*x)
      skp = factor*bigkp + sk/(2.d0*x)
      ldi = sip/si
      ldk = skp/sk
      return
      END SUBROUTINE Mysphbesik
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,PI=3.141592653589793d0)
!CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e
      DOUBLE PRECISION f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q
      DOUBLE PRECISION r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1
      DOUBLE PRECISION temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0.d0.or.xnu.lt.0.d0) write(6,*) 'bad arguments in bessjy', x, xnu
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
 11   continue
      write(6,*)'x too large in bessjy; try asymptotic expansion', x,xnu
      read(*,*)
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        write(6,*) 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        write(6,*)'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END SUBROUTINE bessjy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine BesselAsymNew(v,x,sj,sy,sjp,syp)
      integer k,b,kstop
      double precision v,x,sj,sy,sjp,syp,u
      double precision P,Q,R,S,pi,factor,chi
      double precision denomP,denomQ,denomP2,denomQ2
      double precision Pfunction,Qfunction,sign
      double precision coefR,coefS,factrl

c     asymptotic expansions from Abramowitz and Stegun, pg. 364
c     following assumes v = n+1/2 where n=0,1,2,......

      pi = dacos(-1.0d0)
      factor = dsqrt(2.0d0/(pi*x))
      chi = x - pi*(0.5d0*v + 0.25d0)
      u = 4.0d0*v*v

      P = 0.0d0
      R = 0.0d0
      do k=0,int(v/2.d0-0.25d0)
         sign = (-1.0d0)**k
         denomP = (2.0d0*x)**(2*k)
         denomP2 = factrl(2*k)
         coefR = (u-1.0d0+16.0d0*dfloat(k*k))/
     .       (u-(4.0d0*dfloat(k)-1.0d0)*(4.0d0*dfloat(k)-1.0d0))
         Pfunction = 1.0d0
         do b=0,4*k-1
          Pfunction = Pfunction*(v-0.5d0+float(2*k-b))
         enddo
      P = P + sign*Pfunction/(denomP*denomP2)
      R = R + sign*coefR*Pfunction/(denomP*denomP2)
      enddo

      Q = 0.0d0
      S = 0.0d0
      do k=0,int(v/2.d0-0.75d0)
         sign = (-1.0d0)**k
         denomQ = denomP*(2.0d0*x)
         denomQ2 = factrl(2*k+1)
         coefS = (u+4.0d0*(2.0d0*dfloat(k)+1.0d0)*(2.0d0*dfloat(k)+1.0d0)-
     .           1.0d0)/(u-(4.0d0*dfloat(k)+1.0d0)*(4.0d0*dfloat(k)+1.0d0))
         Qfunction = 1.0d0
         do b=0,4*k+1
         Qfunction = Qfunction*(v+0.5d0+float(2*k-b))
         enddo       
      Q = Q + sign*Qfunction/(denomQ*denomQ2)
      S = S + sign*coefS*Qfunction/(denomQ*denomQ2)
      enddo


      if (v.eq.0.5d0) then
         R = 1.d0
         S = 1.d0/(2.d0*x)
      endif

      sj = factor*(P*dcos(chi) - Q*dsin(chi)) 
      sy = factor*(P*dsin(chi) + Q*dcos(chi))
      sjp=-factor*(R*dsin(chi) - S*dcos(chi))
      syp= factor*(R*dcos(chi) - S*dsin(chi))

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
!CU    USES chebev
      DOUBLE PRECISION xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END SUBROUTINE beschb
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      DOUBLE PRECISION chebev,a,b,x,c(m)
      INTEGER j
      DOUBLE PRECISION d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) write(6,*)'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END FUNCTION chebev
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     energy normalized!!!
      Subroutine sphbessel(rr,energy,LL,xmass,b)
!      use modb
      implicit none
      double precision b(4),prefac,ek,xmass,energy,xx,rr
      double precision sj,sy,sjp,syp,pi
      Integer LL
      pi=dacos(-1d0)
      ek = Dsqrt(2.d0*xmass*energy)
!      prefac = Dsqrt(2.d0*xmass/(ek*pi))
      prefac = Dsqrt(2.d0/(ek*pi))
      xx = ek*rr
      call sphbesjy(LL,xx,sj,sy,sjp,syp)

      b(1)=prefac*xx*sj
      b(2)=prefac*ek*(sj+xx*sjp)
      b(3)= prefac*xx*sy
      b(4)=prefac*ek*(sy+xx*syp)
      write(6,*)b

      End
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MyScaledBessIK(x, xnu, alpha, beta, alphap, betap, ldi,ldk)
!     returns the modified bessel functions Inu and Knu multiplied by exp(-x) and exp(x) respectively.  Uses asymptotic forms for large x>15.0
!     log derivatives of the bessel functions themselves returned in ldk and ldi
!
!     Inu(x) = exp(x) * alpha(x)
!     Knu(x) = exp(-x) * beta(x)
!
!     This routine returns alpha, beta, alpha', beta', I'/I, and K'/K
!      use modb
      implicit none
      double precision mu, x, xnu, alpha, beta, alphap, betap,ldi,ldk
      double precision ri, rip, rk, rkp, xswitch,pi
      pi = dacos(-1d0)
      mu = 4.d0*xnu**2
      xswitch = 15.d0!*(mu-1)/8.d0
      if(x.gt.xswitch) then
!         print*, 'returning asymptotic form'
         alpha =  1.d0/sqrt(2.d0*pi*x)*(1.d0
     >        - (mu - 1.d0)/(8.d0*x)
     >        + (mu - 1.d0)*(mu - 9.d0)/(2.d0*(8.d0*x)**2)
     >        - (mu - 1.d0)*(mu - 9.d0)*(mu - 25.d0)/(6.d0*(8.d0*x)**3)
     >        + (mu - 1.d0)*(mu - 9.d0)*(mu - 25.d0)*(mu - 49.d0)/(24.d0*(8.d0*x)**4))

         alphap = ((1.d0/x)**5.5d0*(-99225.d0 + 464976.d0*xnu**2 -
     -        32.d0*(3.d0*x*(525.d0 + 8.d0*x*(45.d0 + 16.d0*x*(3.d0 + 8.d0*x))) -
     -        4.d0*x*(1813.d0 + 48.d0*x*(25.d0 + 24.d0*x))*xnu**2 +
     -        (8883.d0 + 80.d0*x*(49.d0 + 24.d0*x))*xnu**4 -
     -        56.d0*(27.d0 + 8.d0*x)*xnu**6.d0 + 72.d0*xnu**8)))/
     -        (196608.d0*Sqrt(2.d0*Pi))

         beta =  sqrt(pi/(2.d0*x))*(1.d0 + (mu - 1.d0)/(8.d0*x)
     >        + (mu - 1.d0)*(mu - 9.d0)/(2.d0*(8.d0*x)**2)
     >        + (mu - 1.d0)*(mu - 9.d0)*(mu - 25.d0)/(6.d0*(8.d0*x)**3)
     >        + (mu - 1.d0)*(mu - 9.d0)*(mu - 25.d0)*(mu - 49.d0)/(24.d0*(8.d0*x)**4))

         betap = (Sqrt(Pi/2.d0)*(1/x)**5.5d0*
     -        (-99225 + 464976*xnu**2 -
     -        32*(3*x*(-525 + 8*x*(45 + 16*x*(-3 + 8*x))) +
     -        4*x*(1813 + 48*x*(-25 + 24*x))*xnu**2 +
     -        (8883 + 80*x*(-49 + 24*x))*xnu**4 +
     -        56*(-27 + 8*x)*xnu**6 + 72*xnu**8)))/196608.d0

      else

         call bessik(x,xnu,ri,rk,rip,rkp)
         alphap = dexp(-x)*(rip-ri)
         alpha = dexp(-x)*ri
         betap = dexp(x)*(rkp+rk)
         beta = dexp(x)*rk
      endif
      ldi = -1.d0 +  alphap/alpha
      ldk = +1.d0 +  betap/beta

      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      INTEGER MAXIT
      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-10,FPMIN=1.d-30,MAXIT=10000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     *gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,
     *ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) then
         read(*,*)
         write(6,*)  'bad arguments in bessik'
      endif
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN) h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
         b=b+xi2
         d=1.d0/(b+d)
         c=b+1.d0/c
         del=c*d
         h=del*h
         if(abs(del-1.d0).lt.EPS) goto 1
 11   continue
      write(6,*)  'x too large in bessik; try asymptotic expansion: x = ', x
      read(*,*)
 1    continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
         ritemp=fact*ril+ripl
         fact=fact-xi
         ripl=fact*ritemp+ril
         ril=ritemp
 12   continue
      f=ripl/ril
      if(x.lt.XMIN) then
         x2=.5d0*x
         pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
 13    continue
       write(6,*)  'bessk series failed to converge'
        read(*,*)
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
 14    continue
       write(6,*)  'bessik: failure to converge in cf2'
       read(*,*)
3       continue
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do 15 i=1,nl
        rktemp=(xmu+i)*xi2*rk1+rkmu
        rkmu=rk1
        rk1=rktemp
15    continue
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.


      double precision FUNCTION gammln(xx)
      double precision xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.


      double precision function mygamma(xx)
      double precision gammln, xx
      mygamma = dexp(gammln(xx))
      end
      
      double precision FUNCTION factrl(n)
      INTEGER n
CU    USES gammln
      INTEGER j,ntop
      double precision a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(float(n)+1.0d0))
      endif
      return
      END
