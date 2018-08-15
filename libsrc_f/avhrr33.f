c     Response function factor = 0.96178
      subroutine avhrr33(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(4,mm)
      real wght(4)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *uch4(mm),tau(6,mm),tauch4(6,mm),tauf(6,6,mm),x(mm,8),y(mm,8),
     *sx(6,6,mm),sy(6,6,mm),f(6),fch4(6),fac(mm),tsx(mq),tsy(mq)
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
c0      omega=2660.5
c0      delw=360.0
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
c1      omega=2525.0
c1      delw=70.0
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
c2      omega=2595.0
c2      delw=70.0
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
      omega=2667.5
      delw=75.0
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
c4      omega=2740.0
c4      delw=70.0
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
c5      omega=2817.5
c5      delw=85.0
      omegas=omega
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_33(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
c       u(m)=fac(m)*(u0(m)+u0(m+1))/2.0*6.023e+23/2.687e+19/18.01534
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*ch4*1e-6*fac(m)*
     *       1.0e+05/6.023e+23*16.04303
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_33(back,omega,t0(1))
      call flux_33(bsun,omegas,tsun)
      bsun=bsun/fsun*sconst
      call ck_33(ni,mlv,u,delw,f,p,t,tau)
      call ckch4_33(nch4,mlv,uch4,delw,fch4,p,t,tauch4)
      do 10 m=1,nlyr
        do i=1,ni
          do ia=1,nch4
            tauf(i,ia,m)=tau(i,m)+tauch4(ia,m)
c           tauf(i,ia,m)=tauch4(ia,m)
c           tauf(i,ia,m)=tau(i,m)
          end do
        end do
 10   continue
      do 90 i=1,ni
       do 85 ia=1,nch4
        do 30 m=nlyr,1,-1
          do 20 j=1,8
           if(m.eq.nlyr)then
c           y(m,j)=bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
            y(m,j)=bsun*0.6*exp(-tauf(i,ia,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
           else
            y(m,j)=y(m+1,j)*exp(-tauf(i,ia,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
           end if
 20       continue
 30     continue
        do 50 m=1,nlyr
          sx(i,ia,m)=0.0
          sy(i,ia,m)=0.0
          do 40 j=1,8
           if(m.eq.1)then
c            x(m,j)=back*exp(-tauf(i,ia,m)/xi(j))+
c     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
            x(m,j)=(0.10*y(m,j)+0.90*back)*exp(-tauf(i,ia,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
           else
            x(m,j)=x(m-1,j)*exp(-tauf(i,ia,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
           end if
 40       continue
 50     continue
        do 70 m=1,nlyr
          do 60 j=1,8
            sx(i,ia,m)=sx(i,ia,m)+2.0*x(m,j)*wi(j)*xi(j)
            sy(i,ia,m)=sy(i,ia,m)+2.0*y(m,j)*wi(j)*xi(j)
 60       continue
 70     continue
        do 80 m=1,nlyr
           tsx(m+1)=tsx(m+1)+sx(i,ia,m)*delw*f(i)*fch4(ia)
           tsy(m)=tsy(m)+sy(i,ia,m)*delw*f(i)*fch4(ia)
 80     continue
        tb=tb+back*delw*f(i)*fch4(ia)
        tsx(1)=tb
 85    continue
 90   continue
      tb=0.10*tsy(1)+0.90*tb
      tsx(1)=tb
      tsy(nlev)=tsy(nlev)+bsun*0.60*delw
      do 95 m=1,nlev
cbm         write(6,100)m-1,tsx(m),tb-tsx(m),tsy(m)
 95   continue
      do i1=1,2
       do i2=1,2
        do m=1,nlyr
         od(1+(i2-1)+(i1-1)*2,m) =
     *tauf(i1,i2,m)
         wght(1+(i2-1)+(i1-1)*2) =
     *f(i1) * fch4(i2)
        enddo
       enddo
      enddo
      ntau=4
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_33(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_33(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.933205
      f(2)=0.066795
      k(1)=0.01953061
      do i=2,ni
        k(i)=48.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0386,0.0386,0.0386,0.0386,0.0386,0.0386,0.0386,0.0386,0.0386,
     *0.0599,0.0914,0.1366,0.2010,0.2930,0.4182,0.5861,0.8354,1.1911,
     *1.6361,
     *6.587e-05,6.587e-05,6.587e-05,6.587e-05,6.587e-05,6.587e-05,
     *6.587e-05,6.587e-05,6.587e-05,5.450e-05,3.075e-05,-1.113e-05,
     *-9.463e-05,-2.285e-04,-4.003e-04,-7.846e-04,-1.541e-03,-2.457e-03,
     *-3.263e-03,
     *1.531e-07,1.531e-07,1.531e-07,1.531e-07,1.531e-07,1.531e-07,
     *1.531e-07,1.531e-07,1.531e-07,2.875e-07,3.250e-07,5.469e-07,
     *8.031e-07,7.500e-07,7.875e-07,2.466e-06,3.372e-06,3.375e-06,
     *4.438e-07/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.2684,1.2684,1.2684,1.2684,1.2684,1.2684,1.2684,1.2684,1.2684,
     *1.2613,1.2526,1.2388,1.2195,1.1925,1.1565,1.1065,1.0439,0.9380,
     *0.7896,
     *-1.517e-03,-1.517e-03,-1.517e-03,-1.517e-03,-1.517e-03,-1.517e-03,
     *-1.517e-03,-1.517e-03,-1.517e-03,-1.535e-03,-1.530e-03,-1.501e-03,
     *-1.475e-03,-1.434e-03,-1.397e-03,-1.272e-03,-1.019e-03,-5.909e-04,
     *-3.600e-04,
     *-4.331e-06,-4.331e-06,-4.331e-06,-4.331e-06,-4.331e-06,-4.331e-06,
     *-4.331e-06,-4.331e-06,-4.331e-06,-3.694e-06,-3.491e-06,-3.741e-06,
     *-3.363e-06,-3.316e-06,-3.931e-06,-3.337e-06,-5.700e-06,-5.041e-06,
     *-2.831e-06/
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
c**********************************************************************
      subroutine ckch4_33(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.833250
      f(2)=0.166750
      k(1)=2.7367128
      do i=2,ni
        k(i)=24.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0433,0.0433,0.0433,0.0433,0.0433,0.0433,0.0433,0.0433,0.0433,
     *0.0685,0.1059,0.1605,0.2364,0.3391,0.4776,0.6578,0.8761,1.1501,
     *1.4554,
     *7.837e-05,7.837e-05,7.837e-05,7.837e-05,7.837e-05,7.837e-05,
     *7.837e-05,7.837e-05,7.837e-05,1.093e-04,1.633e-04,2.420e-04,
     *3.453e-04,4.586e-04,6.061e-04,8.044e-04,9.617e-04,1.254e-03,
     *1.535e-03,
     *-2.906e-07,-2.906e-07,-2.906e-07,-2.906e-07,-2.906e-07,-2.906e-07,
     *-2.906e-07,-2.906e-07,-2.906e-07,-6.937e-07,-1.050e-06,-1.369e-06,
     *-1.800e-06,-1.809e-06,-2.691e-06,-3.972e-06,-6.725e-06,-8.666e-06,
     *-9.725e-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.2130,1.2130,1.2130,1.2130,1.2130,1.2130,1.2130,1.2130,1.2130,
     *1.2074,1.2007,1.1892,1.1737,1.1521,1.1232,1.0850,1.0295,0.9803,
     *0.9205,
     *2.866e-03,2.866e-03,2.866e-03,2.866e-03,2.866e-03,2.866e-03,
     *2.866e-03,2.866e-03,2.866e-03,2.848e-03,2.842e-03,2.821e-03,
     *2.809e-03,2.786e-03,2.761e-03,2.716e-03,2.601e-03,2.579e-03,
     *2.537e-03,
     *1.194e-06,1.194e-06,1.194e-06,1.194e-06,1.194e-06,1.194e-06,
     *1.194e-06,1.194e-06,1.194e-06,1.822e-06,1.103e-06,1.694e-06,
     *1.194e-06,1.356e-06,1.622e-06,1.713e-06,1.700e-06,2.503e-06,
     *2.784e-06/
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
