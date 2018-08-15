c     Response function factor = 0.96461
      subroutine avhrr12(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(9,mm)
      real wght(9)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *ux(mq),uo3(mm),uo2(mm),tau(6,mm),tauo3(6,mm),tauo2(6,mm),
     *tauf(6,6,mm),delw,x(mm,8),y(mm,8),sx(6,6,mm),sy(6,6,mm),
     *f(6),fo2(6),fac(mm),tsx(mq),tsy(mq),fo3(6)
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
      omega=15250.0
      delw=1500.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_12(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
        uo2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.20954*fac(m)*
     *         1.0e+05/6.023e+23*31.9988
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_12(back,omega,t0(1))
      call flux_12(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_12(ni,mlv,u,delw,f,p,t,tau)
      call cko3_12(no3,mlv,uo3,delw,fo3,p,t,tauo3)
      call cko2_12(no2,mlv,uo2,delw,fo2,p,t,tauo2)
c     **fo3=1.0000**
      do 10 m=1,nlyr
        do i=1,ni
          do ia=1,no2
            tauf(i,ia,m)=tau(i,m)+tauo3(1,m)+tauo2(ia,m)
c            tauf(i,ia,m)=tau(i,m)
c            tauf(i,ia,m)=tauo2(ia,m)
c            tauf(i,ia,m)=tauo3(1,m)
          end do
        end do
 10   continue
      do 90 i=1,ni
       do 85 ia=1,no2
        do 30 m=nlyr,1,-1
          do 20 j=1,8
           if(m.eq.nlyr)then
c           y(m,j)=bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
            y(m,j)=bsun*0.6*exp(-tauf(i,ia,m)/xi(j))+
     *             bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
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
           tsx(m+1)=tsx(m+1)+sx(i,ia,m)*delw*f(i)*fo2(ia)
           tsy(m)=tsy(m)+sy(i,ia,m)*delw*f(i)*fo2(ia)
 80     continue
c        tb=tb+(0.10*tsy(1))+0.90*back*delw*f(i)*fo2(ia)
        tb=tb+back*delw*f(i)*fo2(ia)
c        tsx(1)=tsx(1)+0.10*tsy(1)
c        tsy(nlev)=tsy(nlev)+bsun*0.60*delw*f(i)*fo2(ia)
 85    continue
 90   continue
      tb=0.10*tsy(1)+0.90*tb
      tsx(1)=tb
      tsy(nlev)=tsy(nlev)+bsun*0.60*delw
      do 95 m=1,nlev
cbm         write(6,100)m-1,tsx(m),tb-tsx(m),tsy(m)
 95   continue
      do i1=1,3
       do i2=1,3
        do m=1,nlyr
         od(1+(i2-1)+(i1-1)*3,m) =
     *tauf(i1,i2,m)
         wght(1+(i2-1)+(i1-1)*3) =
     *f(i1) * fo2(i2)
        enddo
       enddo
      enddo
      ntau=9
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_12(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_12(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
      f(1)=0.957728
      f(2)=0.030620
      f(3)=0.011652
      k(1)=0.0009989559
      do i=2,ni
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0939,0.0939,0.0939,0.0939,0.0939,0.0939,0.0939,0.0939,0.0939,
     *0.1177,0.1531,0.2052,0.2804,0.3847,0.5232,0.6964,0.8954,1.1242,
     *1.3656,
     *5.264e-04,5.264e-04,5.264e-04,5.264e-04,5.264e-04,5.264e-04,
     *5.264e-04,5.264e-04,5.264e-04,4.783e-04,4.101e-04,3.150e-04,
     *1.790e-04,-1.063e-05,-2.395e-04,-4.475e-04,-6.051e-04,-7.571e-04,
     *-8.071e-04,
     *-3.127e-09,-3.127e-09,-3.127e-09,-3.127e-09,-3.127e-09,-3.127e-09,
     *-3.127e-09,-3.127e-09,-3.127e-09,-1.937e-07,-4.844e-07,-7.313e-07,
     *-1.350e-06,-1.653e-06,-1.744e-06,-1.369e-06,-1.719e-07,4.344e-07,
     *4.781e-07/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.7227,0.7227,0.7227,0.7227,0.7227,0.7227,0.7227,0.7227,0.7227,
     *0.7305,0.7425,0.7605,0.7878,0.8241,0.8690,0.9250,0.9783,1.0327,
     *1.0655,
     *1.500e-03,1.500e-03,1.500e-03,1.500e-03,1.500e-03,1.500e-03,
     *1.500e-03,1.500e-03,1.500e-03,1.401e-03,1.230e-03,9.871e-04,
     *7.241e-04,4.560e-04,2.803e-04,2.220e-04,2.360e-04,1.479e-04,
     *-4.874e-06,
     *-4.425e-06,-4.425e-06,-4.425e-06,-4.425e-06,-4.425e-06,-4.425e-06,
     *-4.425e-06,-4.425e-06,-4.425e-06,-1.937e-07,-4.844e-07,-7.313e-07,
     *-1.350e-06,-1.653e-06,-1.744e-06,-1.369e-06,-1.719e-07,4.344e-07,
     *4.781e-07/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.5979,1.5979,1.5979,1.5979,1.5979,1.5979,1.5979,1.5979,1.5979,
     *1.5834,1.5605,1.5258,1.4782,1.4107,1.3227,1.2123,1.0757,0.9328,
     *0.7860,
     *-1.363e-03,-1.363e-03,-1.363e-03,-1.363e-03,-1.363e-03,-1.363e-03,
     *-1.363e-03,-1.363e-03,-1.363e-03,-1.290e-03,-1.221e-03,-1.131e-03,
     *-9.840e-04,-8.320e-04,-6.782e-04,-5.007e-04,-3.427e-04,-2.757e-04,
     *-2.515e-04,
     *2.597e-06,2.597e-06,2.597e-06,2.597e-06,2.597e-06,2.597e-06,
     *2.597e-06,2.597e-06,2.597e-06,2.328e-06,1.963e-06,2.856e-06,
     *2.131e-06,2.137e-06,1.638e-06,1.069e-06,8.125e-07,9.250e-07,
     *7.125e-07/
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
      subroutine cko3_12(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=2.21e-21*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
c *********************************************************************
      subroutine cko2_12(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
      f(1)=0.810000
      f(2)=0.183997
      f(3)=0.006003
      k(1)=0.0
      k(2)=0.000052595402
      do i=3,ni
        k(i)=280.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,
     *1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,
     *1.0000,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0448,0.0448,0.0448,0.0448,0.0448,0.0448,0.0448,0.0448,0.0448,
     *0.0676,0.1000,0.1474,0.2162,0.3186,0.4648,0.6549,0.8667,1.2724,
     *1.8227,
     *1.585e-04,1.585e-04,1.585e-04,1.585e-04,1.585e-04,1.585e-04,
     *1.585e-04,1.585e-04,1.585e-04,1.396e-04,1.013e-04,3.900e-05,
     *-4.975e-05,-1.995e-04,-4.345e-04,-7.880e-04,-1.337e-03,-2.225e-03,
     *-3.379e-03,
     *2.750e-07,2.750e-07,2.750e-07,2.750e-07,2.750e-07,2.750e-07,
     *2.750e-07,2.750e-07,2.750e-07,2.281e-07,2.625e-07,4.500e-07,
     *9.813e-07,1.344e-06,2.437e-06,3.406e-06,4.337e-06,6.703e-06,
     *7.969e-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.3076,1.3076,1.3076,1.3076,1.3076,1.3076,1.3076,1.3076,1.3076,
     *1.3041,1.3002,1.2940,1.2833,1.2650,1.2363,1.1659,1.0447,1.0536,
     *1.0073,
     *-1.145e-03,-1.145e-03,-1.145e-03,-1.145e-03,-1.145e-03,-1.145e-03,
     *-1.145e-03,-1.145e-03,-1.145e-03,-1.142e-03,-1.136e-03,-1.142e-03,
     *-1.099e-03,-1.073e-03,-1.015e-03,-8.623e-04,-8.010e-04,-7.271e-04,
     *-5.045e-04,
     *2.550e-06,2.550e-06,2.550e-06,2.550e-06,2.550e-06,2.550e-06,
     *2.550e-06,2.550e-06,2.550e-06,3.281e-06,2.541e-06,2.894e-06,
     *2.775e-06,2.722e-06,2.906e-06,2.606e-06,1.775e-06,2.347e-06,
     *1.025e-06/
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
