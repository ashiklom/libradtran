c     Response function factor = 0.34804
      subroutine avhrr15(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(2,mm)
      real wght(2)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *ux(mq),uo3(mm),tau(6,mm),tauo3(6,mm),delw,x(mm,8),y(mm,8),
     *sx(6,mm),sy(6,mm),f(6),fac(mm),tsx(mq),tsy(mq),fo3(6)
      data xi /.01986,.10167,.23723,.40828,
     *         .59172,.76277,.89833,.98014/
      data wi /.05062,.11119,.15685,.18134,
     *         .18134,.15685,.11119,.05062/
      tb=0
      nlyr = nlev-1
      if(nlyr.GT.mm)then
        write (*,*) 'error: increase mm to at least', nlyr
        write (*,*) '       in avhrr.inc'
        return
      end if
      if(nlev.GT.mq)then
        write (*,*) 'error: increase mq to at least', nlev
        write (*,*) '       in avhrr.inc'
        return
      end if
      omega=17550.0
      delw=500.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_15(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_15(back,omega,t0(1))
      call flux_15(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_15(ni,mlv,u,delw,f,p,t,tau)
      call cko3_15(no3,mlv,uo3,delw,fo3,p,t,tauo3)
c     **fo3=1.0000**
      do 10 m=1,nlyr
        do i=1,ni
          tau(i,m)=tau(i,m)+tauo3(1,m)
c          tau(i,m)=tauo3(1,m)
        end do
 10   continue
      do 90 i=1,ni
        do 30 m=nlyr,1,-1
          do 20 j=1,8
           if(m.eq.nlyr)then
c           y(m,j)=bnu(m)*(1.0-exp(-tau(i,m)/xi(j)))
            y(m,j)=bsun*0.6*exp(-tau(i,m)/xi(j))+
     *             bnu(m)*(1.0-exp(-tau(i,m)/xi(j)))
           else
            y(m,j)=y(m+1,j)*exp(-tau(i,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tau(i,m)/xi(j)))
           end if
 20       continue
 30     continue
        do 50 m=1,nlyr
          sx(i,m)=0.0
          sy(i,m)=0.0
          do 40 j=1,8
           if(m.eq.1)then
            x(m,j)=(0.10*y(m,j)+0.90*back)*exp(-tau(i,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tau(i,m)/xi(j)))
           else
            x(m,j)=x(m-1,j)*exp(-tau(i,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tau(i,m)/xi(j)))
           end if
 40       continue
 50     continue
        do 70 m=1,nlyr
          do 60 j=1,8
            sx(i,m)=sx(i,m)+2.0*x(m,j)*wi(j)*xi(j)
            sy(i,m)=sy(i,m)+2.0*y(m,j)*wi(j)*xi(j)
 60       continue
 70     continue
        do 80 m=1,nlyr
           tsx(m+1)=tsx(m+1)+sx(i,m)*delw*f(i)
           tsy(m)=tsy(m)+sy(i,m)*delw*f(i)
 80     continue
c        tb=tb+(0.10*tsy(1))+0.90*back*delw*f(i)
        tb=tb+back*delw*f(i)
c        tsx(1)=tsx(1)+0.10*tsy(1)
c        tsy(nlev)=tsy(nlev)+bsun*0.60*delw*f(i)
 90   continue
      tb=0.10*tsy(1)+0.90*tb
      tsx(1)=tb
      tsy(nlev)=tsy(nlev)+bsun*0.60*delw
      do 95 m=1,nlev
cbm         write(6,100)m-1,tsx(m),tb-tsx(m),tsy(m)
 95   continue
      do i1=1,2
       do m=1,nlyr
        od(1+(i1-1),m) =
     *tau(i1,m)
        wght(1+(i1-1)) =
     *f(i1)
       enddo
      enddo
      ntau=2
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_15(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_15(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.991211
      f(2)=0.008789
      k(1)=0.00093679579
      do i=2,ni
        k(i)=81.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.5361,0.5361,0.5361,0.5361,0.5361,0.5361,0.5361,0.5361,0.5361,
     *0.5469,0.5628,0.5876,0.6230,0.6749,0.7442,0.8327,0.9398,1.0548,
     *1.1672,
     *5.886e-04,5.886e-04,5.886e-04,5.886e-04,5.886e-04,5.886e-04,
     *5.886e-04,5.886e-04,5.886e-04,5.257e-04,4.337e-04,2.940e-04,
     *9.450e-05,-1.580e-04,-4.641e-04,-7.900e-04,-1.068e-03,-1.252e-03,
     *-1.410e-03,
     *-3.597e-06,-3.597e-06,-3.597e-06,-3.597e-06,-3.597e-06,-3.597e-06,
     *-3.597e-06,-3.597e-06,-3.597e-06,-3.531e-06,-2.994e-06,-2.963e-06,
     *-1.737e-06,-1.287e-06,-7.094e-07,1.875e-07,-1.406e-07,-6.969e-07,
     *-1.625e-07/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.6283,1.6283,1.6283,1.6283,1.6283,1.6283,1.6283,1.6283,1.6283,
     *1.6152,1.5919,1.5593,1.5081,1.4383,1.3422,1.2188,1.0756,0.9140,
     *0.7518,
     *-3.567e-03,-3.567e-03,-3.567e-03,-3.567e-03,-3.567e-03,-3.567e-03,
     *-3.567e-03,-3.567e-03,-3.567e-03,-3.484e-03,-3.346e-03,-3.153e-03,
     *-2.889e-03,-2.526e-03,-2.104e-03,-1.644e-03,-1.258e-03,-1.007e-03,
     *-8.077e-04,
     *7.300e-06,7.300e-06,7.300e-06,7.300e-06,7.300e-06,7.300e-06,
     *7.300e-06,7.300e-06,7.300e-06,6.225e-06,6.137e-06,4.259e-06,
     *4.769e-06,3.050e-06,2.063e-06,2.041e-06,2.287e-06,2.459e-06,
     *2.419e-06/
	data stp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     *	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     *	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0/
      do i=1,ni
      do m=1,mlv
        ml=1
        pmb(m)=p(m)*1013.25
        if(pmb(m).lt.stp(1))then
         x1=coefk(i,1,1)+coefk(i,2,1)*(t(m)-250.0)
     *   +coefk(i,3,1)*(t(m)-250.0)**2
         fkg(i,m)=x1*pmb(m)/stp(1)
        else if (pmb(m).gt.stp(19)) then
         x1=coefk(i,1,18)+coefk(i,2,18)*(t(m)-250.0)
     *   +coefk(i,3,18)*(t(m)-250.0)**2
         x2=coefk(i,1,19)+coefk(i,2,19)*(t(m)-250.0)
     *   +coefk(i,3,19)*(t(m)-250.0)**2
         fkg(i,m)=x1+(x2-x1)/(stp(19)-stp(18))
     *  *(pmb(m)-stp(18))
        else
         do while(pmb(m).ge.stp(ml))
           ml=ml+1
         end do
         x1=coefk(i,1,ml-1)+coefk(i,2,ml-1)*(t(m)-250.0)
     *   + coefk(i,3,ml-1)*(t(m)-250.0)**2
         x2=coefk(i,1,ml)+coefk(i,2,ml)*(t(m)-250.0)
     *   + coefk(i,3,ml)*(t(m)-250.0)**2
         fkg(i,m)=x1+(x2-x1)/(stp(ml)-stp(ml-1))
     *   *(pmb(m)-stp(ml-1))
        end if
      end do
      end do
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)*fkg(i,m)
c         write(6,*)m,i,u(m),k(i),fkg(i,m),tau(i,m)
        end do
      end do
      return
      end
c *********************************************************************
      subroutine cko3_15(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=4.43e-21*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
