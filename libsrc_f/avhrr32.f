      subroutine avhrr32(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(12,mm)
      real wght(12)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *uch4(mm),un2o(mm),tau(6,mm),tauch4(6,mm),taun2o(6,mm),
     *tauf(6,6,6,mm),x(mm,8),y(mm,8),sx(6,6,6,mm),sy(6,6,6,mm),f(6),
     *fch4(6),fn2o(6),fac(mm),tsx(mq),tsy(mq)
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
      omega=2595.0
      delw=70.0
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
c3      omega=2667.5
c3      delw=75.0
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
        call flux_32(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
c       u(m)=fac(m)*(u0(m)+u0(m+1))/2.0*6.023e+23/2.687e+19/18.01534
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*ch4*1e-6*fac(m)*
     *       1.0e+05/6.023e+23*16.04303
        un2o(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*n2o*1e-6*fac(m)*
     *       1.0e+05/6.023e+23*44.0128
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_32(back,omega,t0(1))
      call flux_32(bsun,omegas,tsun)
      bsun=bsun/fsun*sconst
      call ck_32(ni,mlv,u,delw,f,p,t,tau)
      call ckch4_32(nch4,mlv,uch4,delw,fch4,p,t,tauch4)
      call ckn2o_32(nn2o,mlv,un2o,delw,fn2o,p,t,taun2o)
      do 10 m=1,nlyr
        do i=1,ni
         do ia=1,nch4
          do ib=1,nn2o
c           tauf(i,ia,ib,m)=tau(i,m)
c           tauf(i,ia,ib,m)=tauch4(ia,m)
c           tauf(i,ia,ib,m)=taun2o(ib,m)
            tauf(i,ia,ib,m)=tau(i,m)+tauch4(ia,m)+taun2o(ib,m)
          end do
         end do
        end do
 10   continue
      do 90 i=1,ni
       do 85 ia=1,nch4
       do 83 ib=1,nn2o
        do 30 m=nlyr,1,-1
          do 20 j=1,8
           if(m.eq.nlyr)then
c           y(m,j)=bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
            y(m,j)=bsun*0.6*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           else
            y(m,j)=y(m+1,j)*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           end if
 20       continue
 30     continue
        do 50 m=1,nlyr
          sx(i,ia,ib,m)=0.0
          sy(i,ia,ib,m)=0.0
          do 40 j=1,8
           if(m.eq.1)then
c            x(m,j)=back*exp(-tauf(i,ia,ib,m)/xi(j))+
c     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
            x(m,j)=(0.10*y(m,j)+0.90*back)*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           else
            x(m,j)=x(m-1,j)*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           end if
 40       continue
 50     continue
        do 70 m=1,nlyr
          do 60 j=1,8
            sx(i,ia,ib,m)=sx(i,ia,ib,m)+2.0*x(m,j)*wi(j)*xi(j)
            sy(i,ia,ib,m)=sy(i,ia,ib,m)+2.0*y(m,j)*wi(j)*xi(j)
 60       continue
 70     continue
        do 80 m=1,nlyr
           tsx(m+1)=tsx(m+1)+sx(i,ia,ib,m)*delw*f(i)*fch4(ia)*fn2o(ib)
           tsy(m)=tsy(m)+sy(i,ia,ib,m)*delw*f(i)*fch4(ia)*fn2o(ib)
 80     continue
        tb=tb+back*delw*f(i)*fch4(ia)*fn2o(ib)
        tsx(1)=tb
 83    continue
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
        do i3=1,3
         do m=1,nlyr
          od(1+(i3-1)+(i2-1)*3+(i1-1)*3*2,m) =
     *tauf(i1,i2,i3,m)
          wght(1+(i3-1)+(i2-1)*3+(i1-1)*3*2) =
     *f(i1) * fch4(i2) * fn2o(i3)
         enddo
        enddo
       enddo
      enddo
      ntau=12
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_32(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_32(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.832661
      f(2)=0.167339
      k(1)=0.003144265
      do i=2,ni
        k(i)=48.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,0.0209,
     *0.0413,0.0720,0.1181,0.1857,0.2860,0.4310,0.6411,0.8655,1.3127,
     *1.9164,
     *1.466e-04,1.466e-04,1.466e-04,1.466e-04,1.466e-04,1.466e-04,
     *1.466e-04,1.466e-04,1.466e-04,2.931e-04,5.130e-04,8.469e-04,
     *1.339e-03,2.050e-03,3.059e-03,4.495e-03,5.951e-03,9.188e-03,
     *1.396e-02,
     *-3.031e-07,-3.031e-07,-3.031e-07,-3.031e-07,-3.031e-07,-3.031e-07,
     *-3.031e-07,-3.031e-07,-3.031e-07,-5.094e-07,-9.437e-07,-1.303e-06,
     *-1.619e-06,-2.250e-06,-3.119e-06,-5.266e-06,-8.319e-06,-1.471e-05,
     *-2.123e-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.2408,1.2408,1.2408,1.2408,1.2408,1.2408,1.2408,1.2408,1.2408,
     *1.2385,1.2563,1.2300,1.2219,1.2113,1.1955,1.1732,1.0584,1.0669,
     *1.0368,
     *8.653e-03,8.653e-03,8.653e-03,8.653e-03,8.653e-03,8.653e-03,
     *8.653e-03,8.653e-03,8.653e-03,8.640e-03,8.624e-03,8.596e-03,
     *8.552e-03,8.479e-03,8.367e-03,8.219e-03,7.444e-03,7.612e-03,
     *7.355e-03,
     *-1.087e-05,-1.087e-05,-1.087e-05,-1.087e-05,-1.087e-05,-1.087e-05,
     *-1.087e-05,-1.087e-05,-1.087e-05,-1.078e-05,-1.065e-05,-1.063e-05,
     *-1.005e-05,-1.000e-05,-9.950e-06,-9.622e-06,-8.787e-06,-7.806e-06,
     *-7.594e-06/
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
      subroutine ckch4_32(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.854106
      f(2)=0.145894
      k(1)=2.0090607
      do i=2,ni
        k(i)=24.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0349,0.0349,0.0349,0.0349,0.0349,0.0349,0.0349,0.0349,0.0349,
     *0.0576,0.0951,0.1442,0.2190,0.3227,0.4629,0.6455,0.8713,1.1482,
     *1.4573,
     *-5.525e-05,-5.525e-05,-5.525e-05,-5.525e-05,-5.525e-05,-5.525e-05,
     *-5.525e-05,-5.525e-05,-5.525e-05,-8.713e-05,-1.256e-04,-1.724e-04,
     *-2.186e-04,-3.213e-04,-4.741e-04,-7.211e-04,-1.007e-03,-1.513e-03,
     *-2.032e-03,
     *1.875e-08,1.875e-08,1.875e-08,1.875e-08,1.875e-08,1.875e-08,
     *1.875e-08,1.875e-08,1.875e-08,-1.562e-08,-1.344e-07,-4.219e-06,
     *3.437e-08,1.575e-06,3.247e-06,5.997e-06,9.216e-06,1.511e-05,
     *1.507e-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.2321,1.2321,1.2321,1.2321,1.2321,1.2321,1.2321,1.2321,1.2321,
     *1.2273,1.2186,1.2075,1.1900,1.1651,1.1323,1.0885,1.0302,0.9709,
     *0.9005,
     *-1.090e-03,-1.090e-03,-1.090e-03,-1.090e-03,-1.090e-03,-1.090e-03,
     *-1.090e-03,-1.090e-03,-1.090e-03,-1.101e-03,-1.076e-03,-1.088e-03,
     *-1.057e-03,-1.035e-03,-1.004e-03,-9.739e-04,-8.859e-04,-7.886e-04,
     *-6.649e-04,
     *5.572e-06,5.572e-06,5.572e-06,5.572e-06,5.572e-06,5.572e-06,
     *5.572e-06,5.572e-06,5.572e-06,5.206e-06,5.787e-06,5.563e-06,
     *5.694e-06,5.263e-06,4.628e-06,4.459e-06,3.122e-06,2.509e-06,
     *2.503e-06/
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
      subroutine ckn2o_32(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
      f(1)=0.668178
      f(2)=0.248540
      f(3)=0.083282
      k(1)=17.0622393
      do i=2,ni
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.01703,0.01703,0.01703,0.01703,0.01703,0.01703,
     *0.01703,0.01703,0.01703,0.03731,0.06718,0.10397,
     *0.16562,0.25877,0.39737,0.59672,0.86234,1.16586,
     *1.47071,
     *1.965E-04,1.965E-04,1.965E-04,1.965E-04,1.965E-04,1.965E-04,
     *1.965E-04,1.965E-04,1.965E-04,2.875E-04,4.894E-04,7.641E-04,
     *1.162E-03,1.781E-03,2.691E-03,4.006E-03,5.979E-03,8.487E-03,
     *1.088E-02,
     *8.437E-07,8.437E-07,8.437E-07,8.437E-07,8.437E-07,8.437E-07,
     *8.437E-07,8.437E-07,8.437E-07,1.225E-06,8.656E-07,1.834E-06,
     *1.631E-06,2.028E-06,3.047E-06,2.178E-06,1.250E-06,8.769E-06,
     *1.780E-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.05480,0.05480,0.05480,0.05480,0.05480,0.05480,
     *0.05480,0.05480,0.05480,0.07811,0.11178,0.16211,
     *0.23725,0.34539,0.48627,0.65651,0.87918,1.15664,
     *1.45212,
     *3.943E-04,3.943E-04,3.943E-04,3.943E-04,3.943E-04,3.943E-04,
     *3.943E-04,3.943E-04,3.943E-04,5.045E-04,6.142E-04,7.469E-04,
     *9.075E-04,1.192E-03,1.556E-03,1.566E-03,9.789E-04,1.514E-04,
     *-4.506E-04,
     *1.363E-06,1.363E-06,1.363E-06,1.363E-06,1.363E-06,1.363E-06,
     *1.363E-06,1.363E-06,1.363E-06,1.944E-06,3.287E-06,4.072E-06,
     *3.737E-06,3.259E-06,3.912E-06,1.136E-05,1.465E-05,1.280E-05,
     *9.872E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.38653,1.38653,1.38653,1.38653,1.38653,1.38653,
     *1.38653,1.38653,1.38653,1.35160,1.32052,1.29430,
     *1.26967,1.24171,1.20387,1.15143,1.06761,0.94384,
     *0.79703,
     *1.688E-03,1.688E-03,1.688E-03,1.688E-03,1.688E-03,1.688E-03,
     *1.688E-03,1.688E-03,1.688E-03,1.615E-03,1.477E-03,1.339E-03,
     *1.156E-03,9.571E-04,6.251E-04,4.246E-04,3.944E-04,3.654E-04,
     *2.709E-04,
     *6.000E-06,6.000E-06,6.000E-06,6.000E-06,6.000E-06,6.000E-06,
     *6.000E-06,6.000E-06,6.000E-06,5.522E-06,6.284E-06,6.441E-06,
     *6.112E-06,5.891E-06,5.609E-06,3.409E-06,1.603E-06,1.141E-06,
     *9.219E-07/
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
