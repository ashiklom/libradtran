c     Response function factor = 0.98912
      subroutine avhrr24(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(5,mm)
      real wght(5)
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
      omega=13100.0
      delw=400.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_24(bnu(m),omega,t(m))
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
      call flux_24(back,omega,t0(1))
      call flux_24(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_24(ni,mlv,u,delw,f,p,t,tau)
      call cko3_24(no3,mlv,uo3,delw,fo3,p,t,tauo3)
      call cko2_24(no2,mlv,uo2,delw,fo2,p,t,tauo2)
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
      do i1=1,1
       do i2=1,5
        do m=1,nlyr
         od(1+(i2-1)+(i1-1)*5,m) =
     *tauf(i1,i2,m)
         wght(1+(i2-1)+(i1-1)*5) =
     *f(i1) * fo2(i2)
        enddo
       enddo
      enddo
      ntau=5
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_24(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_24(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=0.0
      do i=1,ni
        do m=1,mlv
          tau(i,m)=0.0
        end do
      end do
      return
      end
c *********************************************************************
      subroutine cko3_24(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=2.98e-22*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
c *********************************************************************
      subroutine cko2_24(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=5
      f(1)=0.737335
      f(2)=0.076314
      f(3)=0.137479
      f(4)=0.037464
      f(5)=0.011408
      k(1)=0.00002766977
      do i=2,ni
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     *0.0000,0.0000,0.00186,0.00879,0.0344,0.1410,0.3661,0.7652,1.2757,
     *2.0981,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.300e-05,
     *8.813e-05,4.689e-04,1.435e-03,3.610e-03,7.085e-03,1.092e-02,
     *1.536e-02,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.750e-08,
     *2.156e-07,1.641e-06,1.309e-06,5.587e-06,1.910e-05,2.782e-05,
     *1.949e-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     *0.0000,0.00064,0.00328,0.0313,0.1347,0.2266,0.4068,0.7821,1.2653,
     *2.0214,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,5.637e-06,1.250e-05,
     *9.050e-05,6.830e-04,8.766e-04,8.690e-04,7.786e-04,7.881e-04,
     *8.647e-04,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,1.781e-08,3.750e-08,
     *-6.313e-07,-1.144e-06,2.272e-06,1.013e-06,-8.541e-06,-1.850e-05,
     *-3.014e-05/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.00277,0.00277,0.00277,0.00277,0.00277,0.00277,0.00277,0.00277,
     *0.00277,0.0140,0.0390,0.0773,0.1297,0.1938,0.2822,0.4552,0.7939,
     *1.2534,1.9295,
     *4.050e-05,4.050e-05,4.050e-05,4.050e-05,4.050e-05,4.050e-05,
     *4.050e-05,4.050e-05,4.050e-05,1.009e-04,1.087e-04,1.105e-04,
     *1.171e-04,8.887e-05,1.253e-04,3.450e-05,-3.951e-04,-8.293e-04,
     *-1.453e-03,
     *2.375e-07,2.375e-07,2.375e-07,2.375e-07,2.375e-07,2.375e-07,
     *2.375e-07,2.375e-07,2.375e-07,4.031e-07,7.313e-07,7.000e-07,
     *5.156e-07,4.781e-07,-1.250e-08,-3.751e-08,-4.031e-07,-7.500e-07,
     *1.006e-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0400,0.0400,0.0400,0.0400,0.0400,0.0400,0.0400,0.0400,0.0400,
     *0.0400,0.0400,0.1000,0.2416,0.3333,0.4540,0.6500,0.8391,1.1662,
     *1.4855,
     *2.479e-04,2.479e-04,2.479e-04,2.479e-04,2.479e-04,2.479e-04,
     *2.479e-04,2.479e-04,2.479e-04,2.347e-04,2.466e-04,2.667e-04,
     *3.233e-04,3.931e-04,4.757e-04,4.771e-04,1.279e-04,-1.650e-04,
     *-4.517e-04,
     *2.078e-06,2.078e-06,2.078e-06,2.078e-06,2.078e-06,2.078e-06,
     *2.078e-06,2.078e-06,2.078e-06,1.894e-06,1.516e-06,1.144e-06,
     *1.181e-06,1.178e-06,1.575e-06,1.378e-06,3.281e-07,-9.812e-07,
     *-2.125e-06/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.3805,0.9000,0.8000,0.7000,0.6000,0.5000,0.3000,0.3000,0.2000,
     *0.1900,0.2500,0.5000,1.2463,1.1416,1.0260,1.0139,1.0499,0.9286,
     *0.7629,
     *-1.838e-05,-1.838e-05,-1.838e-05,-1.838e-05,-1.838e-05,-1.838e-05,
     *-1.838e-05,-1.838e-05,-1.838e-05,3.525e-05,1.211e-04,1.869e-04,
     *3.837e-04,5.924e-04,6.476e-04,4.061e-04,-1.177e-04,-1.697e-04,
     *-1.691e-04,
     *-8.469e-07,-8.469e-07,-8.469e-07,-8.469e-07,-8.469e-07,-8.469e-07,
     *-8.469e-07,-8.469e-07,-8.469e-07,-4.000e-07,-7.594e-07,-6.532e-07,
     *-1.594e-06,-1.191e-06,-9.094e-07,-9.686e-08,-1.300e-06,-1.719e-06,
     *-5.344e-07/
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
