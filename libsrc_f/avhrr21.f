c     Response function factor = 0.20102
      subroutine avhrr21(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(4,mm)
      real wght(4)
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
      omega=9900.0
      delw=600.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_21(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_21(back,omega,t0(1))
      call flux_21(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_21(ni,mlv,u,delw,f,p,t,tau)
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
      do i1=1,4
       do m=1,nlyr
        od(1+(i1-1),m) =
     *tau(i1,m)
        wght(1+(i1-1)) =
     *f(i1)
       enddo
      enddo
      ntau=4
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_21(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_21(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=4
      f(1)=0.931585
      f(2)=0.049960
      f(3)=0.017682
      f(4)=0.000773
      k(1)=0.001329955
      do i=2,ni
        k(i)=12.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0176,0.0176,0.0176,0.0176,0.0176,0.0176,0.0176,0.0176,0.0176,
     *0.0340,0.0601,0.1004,0.1612,0.2523,0.3873,0.5786,0.8480,1.1580,
     *1.5010,
     *9.500e-05,9.500e-05,9.500e-05,9.500e-05,9.500e-05,9.500e-05,
     *9.500e-05,9.500e-05,9.500e-05,1.801e-04,3.204e-04,5.326e-04,
     *8.487e-04,1.307e-03,1.962e-03,2.852e-03,4.076e-03,5.406e-03,
     *6.773e-03,
     *-3.125e-08,-3.125e-08,-3.125e-08,-3.125e-08,-3.125e-08,-3.125e-08,
     *-3.125e-08,-3.125e-08,-3.125e-08,-8.438e-08,-1.656e-07,-3.281e-07,
     *-3.688e-07,-3.250e-07,-6.687e-07,-1.175e-06,-1.434e-06,-9.531e-07,
     *1.250e-07/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.2393,0.2393,0.2393,0.2393,0.2393,0.2393,0.2393,0.2393,0.2393,
     *0.2680,0.3091,0.3671,0.4470,0.5526,0.6788,0.8132,0.9488,1.0321,
     *1.0617,
     *2.811e-03,2.811e-03,2.811e-03,2.811e-03,2.811e-03,2.811e-03,
     *2.811e-03,2.811e-03,2.811e-03,2.888e-03,3.009e-03,3.181e-03,
     *3.433e-03,3.765e-03,4.125e-03,4.445e-03,4.832e-03,5.108e-03,
     *5.351e-03,
     *7.441e-06,7.441e-06,7.441e-06,7.441e-06,7.441e-06,7.441e-06,
     *7.441e-06,7.441e-06,7.441e-06,7.034e-06,6.747e-06,6.469e-06,
     *5.903e-06,4.328e-06,3.472e-06,2.878e-06,5.319e-06,4.981e-06,
     *7.497e-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.9871,0.9871,0.9871,0.9871,0.9871,0.9871,0.9871,0.9871,0.9871,
     *0.9877,0.9893,0.9936,0.9991,1.0060,1.0116,1.0118,1.0143,0.9725,
     *0.9195,
     *5.866e-03,5.866e-03,5.866e-03,5.866e-03,5.866e-03,5.866e-03,
     *5.866e-03,5.866e-03,5.866e-03,5.859e-03,5.857e-03,5.859e-03,
     *5.909e-03,6.005e-03,6.207e-03,6.434e-03,6.736e-03,6.882e-03,
     *6.983e-03,
     *1.675e-05,1.675e-05,1.675e-05,1.675e-05,1.675e-05,1.675e-05,
     *1.675e-05,1.675e-05,1.675e-05,1.746e-05,1.794e-05,1.764e-05,
     *1.720e-05,1.609e-05,1.516e-05,1.347e-05,1.266e-05,1.142e-05,
     *1.040e-05/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.9989,1.9989,1.9989,1.9989,1.9989,1.9989,1.9989,1.9989,1.9989,
     *1.9737,1.9330,1.8709,1.7817,1.6594,1.4962,1.3003,1.1075,0.8862,
     *0.6960,
     *1.254e-02,1.254e-02,1.254e-02,1.254e-02,1.254e-02,1.254e-02,
     *1.254e-02,1.254e-02,1.254e-02,1.247e-02,1.235e-02,1.210e-02,
     *1.168e-02,1.102e-02,1.008e-02,8.826e-03,7.467e-03,5.807e-03,
     *4.388e-03,
     *-1.931e-05,-1.931e-05,-1.931e-05,-1.931e-05,-1.931e-05,-1.931e-05,
     *-1.931e-05,-1.931e-05,-1.931e-05,-1.905e-05,-1.843e-05,-1.796e-05,
     *-1.698e-05,-1.551e-05,-1.215e-05,-9.609e-06,-7.909e-06,-4.919e-06,
     *-5.197e-06/
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
