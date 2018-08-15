c     Response function factor = 0.64054
      subroutine avhrr22(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(5,mm)
      real wght(5)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *ux(mq),tau(6,mm),delw,x(mm,8),y(mm,8),
     *sx(6,mm),sy(6,mm),f(6),fac(mm),tsx(mq),tsy(mq)
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
      omega=10600.0
      delw=800.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_22(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_22(back,omega,t0(1))
      call flux_22(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_22(ni,mlv,u,delw,f,p,t,tau)
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
      do i1=1,5
       do m=1,nlyr
        od(1+(i1-1),m) =
     *tau(i1,m)
        wght(1+(i1-1)) =
     *f(i1)
       enddo
      enddo
      ntau=5
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_22(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_22(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=5
      f(1)=0.082056
      f(2)=0.540035
      f(3)=0.298713
      f(4)=0.072936
      f(5)=0.006260
      k(1)=0.003088904
      do i=2,ni
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     *0.0000,0.00089,0.00762,0.0323,0.0974,0.2229,0.4383,0.7740,1.2958,
     *2.0899,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,-4.750e-06,-4.500e-06,
     *1.475e-05,6.225e-05,1.215e-04,2.625e-04,5.704e-04,9.451e-04,
     *1.794e-03,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,-4.375e-08,-1.063e-07,
     *-4.250e-07,-1.031e-06,-1.731e-06,-3.775e-06,-4.028e-06,-8.328e-06,
     *-1.568e-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0048,0.0048,0.0048,0.0048,0.0048,0.0048,0.0048,0.0048,0.0048,
     *0.0144,0.0330,0.0652,0.1187,0.2040,0.3324,0.5251,0.8187,1.2202,
     *1.7427,
     *1.425e-05,1.425e-05,1.425e-05,1.425e-05,1.425e-05,1.425e-05,
     *1.425e-05,1.425e-05,1.425e-05,3.850e-05,7.850e-05,1.431e-04,
     *2.311e-04,3.545e-04,5.670e-04,8.243e-04,1.211e-03,1.816e-03,
     *2.560e-03,
     *-6.875e-08,-6.875e-08,-6.875e-08,-6.875e-08,-6.875e-08,-6.875e-08,
     *-6.875e-08,-6.875e-08,-6.875e-08,-1.937e-07,-3.375e-07,-6.344e-07,
     *-1.028e-06,-1.494e-06,-2.100e-06,-3.050e-06,-4.481e-06,-6.525e-06,
     *-5.297e-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0853,0.0853,0.0853,0.0853,0.0853,0.0853,0.0853,0.0853,0.0853,
     *0.1094,0.1427,0.1906,0.2579,0.3513,0.4792,0.6491,0.8781,1.1362,
     *1.4139,
     *5.881e-04,5.881e-04,5.881e-04,5.881e-04,5.881e-04,5.881e-04,
     *5.881e-04,5.881e-04,5.881e-04,6.406e-04,6.406e-04,6.811e-04,
     *7.330e-04,7.964e-04,8.890e-04,9.834e-04,1.115e-03,1.235e-03,
     *1.395e-03,
     *1.084e-06,1.084e-06,1.084e-06,1.084e-06,1.084e-06,1.084e-06,
     *1.084e-06,1.084e-06,1.084e-06,7.906e-07,7.031e-07,3.719e-07,
     *4.375e-08,-1.844e-07,-3.500e-07,-8.281e-07,-1.813e-06,-2.475e-06,
     *-3.419e-06/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.5405,0.5405,0.5405,0.5405,0.5405,0.5405,0.5405,0.5405,0.5405,
     *0.5545,0.5764,0.6096,0.6567,0.7180,0.7933,0.8747,0.9673,1.0333,
     *1.0743,
     *2.195e-03,2.195e-03,2.195e-03,2.195e-03,2.195e-03,2.195e-03,
     *2.195e-03,2.195e-03,2.195e-03,2.093e-03,1.947e-03,1.751e-03,
     *1.512e-03,1.249e-03,1.008e-03,7.746e-04,5.860e-04,3.907e-04,
     *1.921e-04,
     *-1.600e-06,-1.600e-06,-1.600e-06,-1.600e-06,-1.600e-06,-1.600e-06,
     *-1.600e-06,-1.600e-06,-1.600e-06,-1.481e-06,-7.656e-07,-1.875e-07,
     *1.688e-07,8.250e-07,2.469e-07,5.312e-08,6.227e-05,-1.125e-07,
     *1.066e-06/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.7739,1.7739,1.7739,1.7739,1.7739,1.7739,1.7739,1.7739,1.7739,
     *1.7518,1.7172,1.6684,1.5989,1.5045,1.3830,1.2377,1.0877,0.9158,
     *0.7495,
     *-1.725e-03,-1.725e-03,-1.725e-03,-1.725e-03,-1.725e-03,-1.725e-03,
     *-1.725e-03,-1.725e-03,-1.725e-03,-1.637e-03,-1.530e-03,-1.364e-03,
     *-1.147e-03,-9.054e-04,-6.973e-04,-5.326e-04,-3.981e-04,-2.995e-04,
     *-2.284e-04,
     *2.200e-06,2.200e-06,2.200e-06,2.200e-06,2.200e-06,2.200e-06,
     *2.200e-06,2.200e-06,2.200e-06,1.269e-06,2.163e-06,1.566e-06,
     *1.259e-06,5.156e-07,1.206e-06,1.072e-06,1.047e-06,1.519e-06,
     *2.084e-06/
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
