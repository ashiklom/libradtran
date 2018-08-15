c     Response function factor = 0.24758 (channel 1)
c     Response function factor = 0.44552 (channel 2)
      subroutine avhrr11(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(12,mm)
      real wght(12)
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
      omega=13900.0
      delw=1200.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_11(bnu(m),omega,t(m))
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
      call flux_11(back,omega,t0(1))
      call flux_11(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_11(ni,mlv,u,delw,f,p,t,tau)
      call cko3_11(no3,mlv,uo3,delw,fo3,p,t,tauo3)
      call cko2_11(no2,mlv,uo2,delw,fo2,p,t,tauo2)
c     **fo3=1.0000**
      do 10 m=1,nlyr
        do i=1,ni
          do ia=1,no2
            tauf(i,ia,m)=tau(i,m)+tauo3(1,m)+tauo2(ia,m)
c            tauf(i,ia,m)=tau(i,m)
c            tauf(i,ia,m)=tauo3(1,m)
c            tauf(i,ia,m)=tauo2(ia,m)
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
      do i1=1,4
       do i2=1,3
        do m=1,nlyr
         od(1+(i2-1)+(i1-1)*3,m) =
     *tauf(i1,i2,m)
         wght(1+(i2-1)+(i1-1)*3) =
     *f(i1) * fo2(i2)
        enddo
       enddo
      enddo
      ntau=12
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_11(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_11(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=4
      f(1)=0.709180
      f(2)=0.215184
      f(3)=0.067047
      f(4)=0.008589
      k(1)=0.003165849
      do i=2,ni
        k(i)=9.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0086,0.0086,0.0086,0.0086,0.0086,0.0086,0.0086,0.0086,0.0086,
     *0.0216,0.0446,0.0831,0.1439,0.2391,0.3809,0.5818,0.8410,1.1868,
     *1.6047,
     *-1.100e-05,-1.100e-05,-1.100e-05,-1.100e-05,-1.100e-05,-1.100e-05,
     *-1.100e-05,-1.100e-05,-1.100e-05,-7.750e-06,-1.413e-05,-1.475e-05,
     *-3.400e-05,-6.000e-05,-9.475e-05,-8.675e-05,-5.000e-05,1.500e-04,
     *6.385e-04,
     *-6.250e-09,-6.250e-09,-6.250e-09,-6.250e-09,-6.250e-09,-6.250e-09,
     *-6.250e-09,-6.250e-09,-6.250e-09,-1.250e-08,-5.937e-08,-2.563e-07,
     *-3.750e-07,-7.937e-07,-1.150e-06,-1.856e-06,-1.619e-06,-2.200e-06,
     *-3.950e-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.1847,0.1847,0.1847,0.1847,0.1847,0.1847,0.1847,0.1847,0.1847,
     *0.2122,0.2504,0.3029,0.3729,0.4659,0.5859,0.7341,0.9012,1.1114,
     *1.3361,
     *6.537e-04,6.537e-04,6.537e-04,6.537e-04,6.537e-04,6.537e-04,
     *6.537e-04,6.537e-04,6.537e-04,5.916e-04,5.249e-04,4.499e-04,
     *3.880e-04,3.387e-04,3.236e-04,3.783e-04,5.040e-04,6.245e-04,
     *5.881e-04,
     *-8.625e-07,-8.625e-07,-8.625e-07,-8.625e-07,-8.625e-07,-8.625e-07,
     *-8.625e-07,-8.625e-07,-8.625e-07,-6.594e-07,-4.531e-07,-2.906e-07,
     *2.688e-07,2.437e-07,8.438e-08,5.438e-07,4.313e-07,5.437e-07,
     *1.209e-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.6392,0.6392,0.6392,0.6392,0.6392,0.6392,0.6392,0.6392,0.6392,
     *0.6485,0.6640,0.6887,0.7247,0.7747,0.8375,0.9102,0.9693,1.0319,
     *1.0717,
     *2.020e-03,2.020e-03,2.020e-03,2.020e-03,2.020e-03,2.020e-03,
     *2.020e-03,2.020e-03,2.020e-03,1.933e-03,1.815e-03,1.635e-03,
     *1.411e-03,1.125e-03,8.254e-04,5.374e-04,2.741e-04,2.975e-05,
     *-1.946e-04,
     *-1.137e-06,-1.137e-06,-1.137e-06,-1.137e-06,-1.137e-06,-1.137e-06,
     *-1.137e-06,-1.137e-06,-1.137e-06,-7.156e-07,-3.625e-07,-2.500e-07,
     *5.625e-07,7.875e-07,1.391e-06,6.781e-07,7.969e-07,2.251e-07,
     *-3.656e-07/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.7033,1.7033,1.7033,1.7033,1.7033,1.7033,1.7033,1.7033,1.7033,
     *1.6847,1.6574,1.6165,1.5555,1.4734,1.3659,1.2350,1.0747,0.9284,
     *0.7825,
     *-2.582e-03,-2.582e-03,-2.582e-03,-2.582e-03,-2.582e-03,-2.582e-03,
     *-2.582e-03,-2.582e-03,-2.582e-03,-2.509e-03,-2.373e-03,-2.208e-03,
     *-1.974e-03,-1.720e-03,-1.463e-03,-1.166e-03,-9.144e-04,-7.999e-04,
     *-7.045e-04,
     *2.887e-06,2.887e-06,2.887e-06,2.887e-06,2.887e-06,2.887e-06,
     *2.887e-06,2.887e-06,2.887e-06,3.194e-06,2.803e-06,1.919e-06,
     *2.519e-06,1.637e-06,1.953e-06,1.647e-06,1.566e-06,1.478e-06,
     *1.913e-06/
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
      subroutine cko3_11(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=6.71e-22*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
c *********************************************************************
      subroutine cko2_11(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
c     Regression coefficients for response function = 1.000000
      f(1)=0.842000
      f(2)=0.154773
      f(3)=0.003227
c     Regression coefficients for AVHRR channel 1 response function
c      f(1)=0.412289
c      f(2)=0.575708
c      f(3)=0.012003
c     Regression coefficients for AVHRR channel 2 response function
c     Note: For this case the spectral signature of O2 may be neglected
c      f(1)=1.000000
c      f(2)=0.000000
c      f(3)=0.000000
      k(1)=0.0
      k(2)=0.000031155735
      do i=3,ni
        k(i)=400.0*k(i-1)
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
     *0.0633,0.0633,0.0633,0.0633,0.0633,0.0633,0.0633,0.0633,0.0633,
     *0.0846,0.1146,0.1572,0.2166,0.3032,0.4326,0.6217,0.8459,1.2979,
     *1.9226,
     *2.597e-04,2.597e-04,2.597e-04,2.597e-04,2.597e-04,2.597e-04,
     *2.597e-04,2.597e-04,2.597e-04,2.865e-04,3.353e-04,4.073e-04,
     *5.233e-04,7.623e-04,1.173e-03,1.769e-03,2.607e-03,3.833e-03,
     *5.477e-03,
     *-5.375e-07,-5.375e-07,-5.375e-07,-5.375e-07,-5.375e-07,-5.375e-07,
     *-5.375e-07,-5.375e-07,-5.375e-07,-6.063e-07,-6.563e-07,-9.500e-07,
     *-1.006e-06,-8.875e-07,-1.606e-06,-2.909e-06,-4.331e-06,-7.578e-06,
     *-1.254e-05/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.2941,1.2941,1.2941,1.2941,1.2941,1.2941,1.2941,1.2941,1.2941,
     *1.2917,1.2878,1.2829,1.2741,1.2595,1.2366,1.1847,1.0536,1.0351,
     *0.9694,
     *3.098e-03,3.098e-03,3.098e-03,3.098e-03,3.098e-03,3.098e-03,
     *3.098e-03,3.098e-03,3.098e-03,3.095e-03,3.090e-03,3.083e-03,
     *3.066e-03,3.011e-03,2.937e-03,2.759e-03,2.573e-03,2.348e-03,
     *2.040e-03,
     *-8.803e-06,-8.803e-06,-8.803e-06,-8.803e-06,-8.803e-06,-8.803e-06,
     *-8.803e-06,-8.803e-06,-8.803e-06,-8.800e-06,-8.778e-06,-8.731e-06,
     *-8.716e-06,-8.391e-06,-8.250e-06,-7.956e-06,-7.366e-06,-7.184e-06,
     *-6.956e-06/
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
