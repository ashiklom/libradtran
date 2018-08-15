c     Continuum calculations according to Clough et al (1994)
      subroutine avhrr41(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      implicit real*8 (v)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(20,mm)
      real wght(20)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *uco2(mm),ucfc(mm),tau(6,mm),tauco2(6,mm),taucfc(6,mm),
     *tauf(6,6,6,mm),x(mm,8),y(mm,8),sx(6,6,6,mm),sy(6,6,6,mm),
     *f(6),fco2(6),fcfc(6),fac(mm),tsx(mq),tsy(mq)
      real dvabs,absrb(2030),pres(mm),war(mm),wn2(mm),wo2(mm),wco2(mm),
     *w(mm)
      data absrb /2030*0.0/
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
      omega=925.0
      delw=90.0
      v1abs=880.0
      v2abs=970.0
      dvabs=5.0
      nptabs=19
      mlv=nlyr
      h=5.0
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_41(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
c       u(m)=fac(m)*(u0(m)+u0(m+1))/2.0*6.023e+23/2.687e+19/18.01534
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        w(m)=u(m)*6.023e+23/18.01534
        uco2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*co2*1e-6*fac(m)*
     *          1.0e+05/6.023e+23*44.00995
        wco2(m)=uco2(m)*6.023e+23/44.00995
        ucfc(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*f12*1e-6*fac(m)*
     *          1.0e+05/6.023e+23*120.9139
        if(z0(m).ge.24.999)
     *ucfc(m)=ucfc(m)*exp((-7.5-(z0(m)-20.))/h)
        if(z0(m).lt.24.999.and.z0(m).ge.14.999)
     *ucfc(m)=ucfc(m)*exp((14.5-z0(m))/h)
        wn2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.78084*fac(m)*
     *         1.0e+05
        wo2(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.20948*fac(m)*
     *         1.0e+05
        war(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*0.00934*fac(m)*
     *         1.0e+05
        pres(m)=p(m)*1013.25
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_41(back,omega,t0(1))
      call ck_41(ni,mlv,u,delw,f,p,t,tau)
      call ckco2_41(nco2,mlv,uco2,delw,fco2,p,t,tauco2)
      call ckcfc_41(ncfc,mlv,ucfc,delw,fcfc,p,t,taucfc)
      do 10 m=1,nlyr
        do k=1,nptabs
          absrb(k)=0.0
        end do
        call clough_41(pres(m),t(m),war(m),wn2(m),wo2(m),wco2(m),w(m),
     *       v1abs,v2abs,dvabs,nptabs,absrb)
        twin=0.5*(absrb(1)+absrb(19))
        do k=2,nptabs-1
          twin=twin+absrb(k)
        end do
        twin=twin/18.0
        do i=1,ni
          do ia=1,nco2
            do ib=1,ncfc
c             tauf(i,ia,ib,m)=twin
c             tauf(i,ia,ib,m)=tau(i,m)+twin
c             tauf(i,ia,ib,m)=tauco2(ia,m)
c             tauf(i,ia,ib,m)=taucfc(ib,m)
              tauf(i,ia,ib,m)=tau(i,m)+twin+tauco2(ia,m)+taucfc(ib,m)
            end do
          end do
        end do
 10   continue
      do 90 i=1,ni
       do 85 ia=1,nco2
       do 82 ib=1,ncfc
        do 30 m=1,nlyr
          sx(i,ia,ib,m)=0.0
          sy(i,ia,ib,m)=0.0
          do 20 j=1,8
           if(m.eq.1)then
            x(m,j)=back*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           else
            x(m,j)=x(m-1,j)*exp(-tauf(i,ia,ib,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           end if
 20       continue
 30     continue
        do 50 m=nlyr,1,-1
          do 40 j=1,8
           if(m.eq.nlyr)then
            y(m,j)=bnu(m)*(1.0-exp(-tauf(i,ia,ib,m)/xi(j)))
           else
            y(m,j)=y(m+1,j)*exp(-tauf(i,ia,ib,m)/xi(j))+
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
           tsx(m+1)=tsx(m+1)+sx(i,ia,ib,m)*delw*f(i)*fco2(ia)*fcfc(ib)
           tsy(m)=tsy(m)+sy(i,ia,ib,m)*delw*f(i)*fco2(ia)*fcfc(ib)
 80     continue
        tb=tb+back*delw*f(i)*fco2(ia)*fcfc(ib)
        tsx(1)=tb
 82    continue
 85    continue
 90   continue
      do 95 m=1,nlev
cbm         write(6,100)m-1,tsx(m),tb-tsx(m),tsy(m)
 95   continue
      do i1=1,5
       do i2=1,2
        do i3=1,2
         do m=1,nlyr
          od(1+(i3-1)+(i2-1)*2+(i1-1)*2*2,m) =
     *tauf(i1,i2,i3,m)
          wght(1+(i3-1)+(i2-1)*2+(i1-1)*2*2) =
     *f(i1) * fco2(i2) * fcfc(i3)
         enddo
        enddo
       enddo
      enddo
      ntau=20
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_41(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_41(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=5
      f(1)=0.765304
      f(2)=0.128953
      f(3)=0.082958
      f(4)=0.017931
      f(5)=0.004854
      k(1)=0.000598265
      do i=2,ni
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0000167,
     *0.00102,0.00685,0.0248,0.0679,0.1501,0.2830,0.4910,0.7879,
     *1.2714,1.9667,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,2.801e-07,4.075e-05,2.989e-04,1.088e-03,
     *2.712e-03,5.415e-03,9.694e-03,1.664e-02,2.641e-02,4.239e-02,
     *6.551e-02,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,0.000e+00,0.000e+00,4.937e-07,3.916e-06,1.440e-05,
     *3.386e-05,6.345e-05,1.105e-04,1.881e-04,2.961e-04,4.738e-04,
     *7.334e-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.00000806,0.00397,
     *0.0143,0.0384,0.0829,0.1471,0.2377,0.3721,0.5681,0.8267,1.2174,
     *1.7459,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *1.196e-07,2.385e-05,1.710e-04,6.129e-04,1.523e-03,2.969e-03,
     *4.899e-03,7.784e-03,1.218e-02,1.857e-02,2.708e-02,3.962e-02,
     *5.674e-02,
     *0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,
     *0.000e+00,2.931e-07,2.225e-06,7.972e-06,1.888e-05,3.427e-05,
     *5.417e-05,8.588e-05,1.346e-04,2.064e-04,3.014e-04,4.369e-04,
     *6.262e-04/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0,0.0,0.0,0.0000476,0.000577,0.00250,0.00735,0.0181,0.0365,
     *0.0639,0.1024,0.1554,0.2314,0.3385,0.4868,0.6717,0.8775,1.1431,
     *1.4355,
     *0.000e+00,0.000e+00,0.000e+00,1.431e-06,2.650e-05,1.176e-04,
     *3.227e-04,7.110e-04,1.317e-03,2.183e-03,3.414e-03,5.135e-03,
     *7.596e-03,1.102e-02,1.576e-02,2.188e-02,2.833e-02,3.693e-02,
     *4.499e-02,
     *0.000e+00,0.000e+00,0.000e+00,1.071e-08,3.581e-07,1.634e-06,
     *4.294e-06,8.875e-06,1.552e-05,2.488e-05,3.848e-05,5.792e-05,
     *8.533e-05,1.227e-04,1.738e-04,2.451e-04,3.156e-04,4.193e-04,
     *5.074e-04/
      data ( (coefk(4,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0111,0.0156,0.0150,0.0224,0.0313,0.0444,0.0626,
     *0.0857,0.1164,0.1603,0.2226,0.3035,0.4011,0.5187,0.6469,0.7962,
     *0.9257,1.0795,1.1977,
     *7.722e-04,8.094e-04,8.878e-04,1.079e-03,1.387e-03,1.785e-03,
     *2.331e-03,3.056e-03,4.094e-03,5.540e-03,7.560e-03,1.019e-02,
     *1.334e-02,1.712e-02,2.119e-02,2.454e-02,2.643e-02,2.820e-02,
     *2.940e-02,
     *1.298e-05,1.115e-05,1.401e-05,1.511e-05,1.882e-05,2.287e-05,
     *2.837e-05,3.642e-05,4.893e-05,6.557e-05,8.772e-05,1.168e-04,
     *1.525e-04,1.964e-04,2.466e-04,2.751e-04,2.839e-04,2.764e-04,
     *2.602e-04/
      data ( (coefk(5,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.8026,1.8013,1.7999,1.7975,1.7946,1.7903,1.7857,1.7770,
     *1.7657,1.7467,1.7146,1.6669,1.6006,1.5096,1.3988,1.2596,1.0787,
     *0.9212,0.7461,
     *4.801e-02,4.797e-02,4.790e-02,4.780e-02,4.765e-02,4.750e-02,
     *4.719e-02,4.676e-02,4.619e-02,4.535e-02,4.416e-02,4.252e-02,
     *4.027e-02,3.732e-02,3.368e-02,2.981e-02,2.543e-02,2.172e-02,
     *1.756e-02,
     *4.678e-04,4.674e-04,4.665e-04,4.657e-04,4.642e-04,4.633e-04,
     *4.589e-04,4.541e-04,4.467e-04,4.365e-04,4.236e-04,4.051e-04,
     *3.789e-04,3.449e-04,3.010e-04,2.613e-04,2.183e-04,1.884e-04,
     *1.577e-04/
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
      subroutine ckco2_41(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.982902
      f(2)=0.017098
      k(1)=0.009177011
      k(2)=64.0*k(1)
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.00221,0.00345,0.00554,0.00875,0.0132,0.0197,0.0294,0.0432,
     *0.0635,0.0931,0.1347,0.1888,0.2534,0.3404,0.4690,0.6445,0.8630,
     *1.1031,1.3334,
     *1.647e-04,2.326e-04,3.267e-04,4.613e-04,6.520e-04,9.253e-04,
     *1.312e-03,1.841e-03,2.609e-03,3.744e-03,5.349e-03,7.415e-03,
     *1.002e-02,1.305e-02,1.719e-02,2.286e-02,3.001e-02,3.796e-02,
     *4.559e-02,
     *2.856e-06,3.866e-06,5.106e-06,6.800e-06,9.319e-06,1.290e-05,
     *1.776e-05,2.413e-05,3.322e-05,4.682e-05,6.605e-05,9.045e-05,
     *1.239e-04,1.585e-04,2.009e-04,2.590e-04,3.330e-04,4.163e-04,
     *4.957e-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *2.1075,2.1008,2.0902,2.0741,2.0506,2.0144,1.9763,1.9195,1.8559,
     *1.7967,1.7389,1.6784,1.6199,1.5448,1.4355,1.2796,1.0863,0.8725,
     *0.6639,
     *6.946e-02,6.918e-02,6.890e-02,6.847e-02,6.757e-02,6.640e-02,
     *6.555e-02,6.383e-02,6.182e-02,6.025e-02,5.854e-02,5.593e-02,
     *5.372e-02,5.123e-02,4.806e-02,4.330e-02,3.716e-02,3.012e-02,
     *2.307e-02,
     *7.279e-04,7.238e-04,7.214e-04,7.182e-04,7.057e-04,6.930e-04,
     *6.878e-04,6.704e-04,6.495e-04,6.378e-04,6.230e-04,5.883e-04,
     *5.617e-04,5.347e-04,5.068e-04,4.618e-04,4.005e-04,3.275e-04,
     *2.523e-04/
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
      subroutine ckcfc_41(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=2
c     f(1)=ratio of spectral interval used vs. cfc spectral interval
c      cor=1.000000
      cor=0.998434
      f(1)=0.77778*cor
      f(2)=1.0-f(1)
      do m=1,mlv
        k(1)=3814.461*(1.0000-2.3034e-04*(t(m)-250.0)+
     *       5.1926e-07*(t(m)-250.0)**2)
        tau(1,m)=k(1)*u(m)
        tau(2,m)=0.0
c       write(6,*)m,u(m),k(1),tau(1,m)
      end do
      return
      end
c**********************************************************************
      subroutine clough_41(pres,temp,wa,wn2,wo2,wco2,wh2o,v1abs,v2abs,
     *dvabs,nptabs,absrb)
      INCLUDE "avhrr.inc"
      implicit real*8 (v)
      dimension absrb(*)
      dimension wk(60)
      data xlosmt/2.68675e+19/
      data wk/60*0/
c
c
c   this program calculates the continuum optical depth
c         for an homogeneous layer
c
c   the following quantities must be specified:
c
c          pressure                   pave (mb)
c
c          temperature                tave ( k)
c
c          column amount
c            nitrogen                 wn2    (molec/cm**2)
c            oxygen                   wk(7)  (molec/cm**2)
c            carbon dioxide           wk(2)  (molec/cm**2)
c            water vapor              wk(1)  (molec/cm**2)
c
c          number of molecules        nmol
c
c          beginning wavenumber       v1abs (cm-1)
c
c          ending wavenumber          v2abs (cm-1)
c
c          sampling interval          dvabs (cm-1)
c
c          number of values           nptabs
c
c
c   the results are in array absorb
c
c   note that for an atmospheric layer:
c
c            wtot   = xlosmt * (pave/1013) * (273/tave) * (path length)
c
c            wbroad = the column amount for all species not explicitly provided
c
c            wk(m)  = (volume mixing ratio) * (column of dry air)
c
c
c
c      note that continua are included for
c
c             m=1      water vapor       air & self         0 - 20,000 cm-1
c             wn2      nitrogen          air             2020 -   2800 cm-1
c
c   the following is an example for a one km path (see cntnout for results)
c
      pave = pres
      tave = temp 
      wk(7) = wo2
      wk(2) = wco2
      wk(1) =wh2o
      wbroad=wn2+wa
      nmol = 7
c
      jrad=1
c
      call contnm_41(jrad,pave,tave,wk,wbroad,nmol,v1abs,v2abs,dvabs,
     *nptabs,absrb)
c
      return
      end
c**********************************************************************
      subroutine xint_41(v1a,v2a,dva,a,afact,vft,dvr3,r3,n1r3,n2r3)
      INCLUDE "avhrr.inc"
c
      implicit double precision (v)
c
c     this subroutine interpolates the a array stored
c     from v1a to v2a in increments of dva using a multiplicative
c     factor afact, into the r3 array from location n1r3 to n2r3 in
c     increments of dvr3
c
      dimension a(*),r3(*)
      data onepl/1.001/, onemi/0.999/, argmin/34./
c
      recdva = 1./dva
      ilo = (v1a+dva-vft)/dvr3+1.+onemi
      ilo = max(ilo,n1r3)
      ihi = (v2a-dva-vft)/dvr3+onemi
      ihi = min(ihi,n2r3)
c
      do 10 i = ilo, ihi
         vi = vft+dvr3*float(i-1)
         j = (vi-v1a)*recdva+onepl
         vj = v1a+dva*float(j-1)
         p = recdva*(vi-vj)
         c = (3.-2.*p)*p*p
         b = 0.5*p*(1.-p)
         b1 = b*(1.-p)
         b2 = b*p
         conti = -a(j-1)*b1+a(j)*(1.-c+b2)+a(j+1)*(c+b1)-a(j+2)*b2
         r3(i) = r3(i)+conti*afact
   10 continue
c
      return
c
      end
      function radfn41 (vi,xkt)
c
      implicit double precision (v)
c
c     function radfn41 calculates the radiation term for the line shape
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c               last modification:    12 august 1991
c
c                  implementation:    r.d. worsham
c
c             algorithm revisions:    s.a. clough
c                                     r.d. worsham
c                                     j.l. moncet
c
c
c                     atmospheric and environmental research inc.
c                     840 memorial drive,  cambridge, ma   02139
c
c----------------------------------------------------------------------
c
c               work supported by:    the arm program
c                                     office of energy research
c                                     department of energy
c
c
c      source of original routine:    afgl line-by-line model
c
c                                             fascod3
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      data onepl/1.001/, onemi/0.999/, argmin/34./
c
c      in the small xviokt region 0.5 is required
c
      xvi = vi
c
      if (xkt.gt.0.0) then
c
         xviokt = xvi/xkt
c
         if (xviokt.le.0.01) then
            radfn41 = 0.5*xviokt*xvi
c
         elseif (xviokt.le.10.0) then
            expvkt = exp(-xviokt)
            radfn41 = xvi*(1.-expvkt)/(1.+expvkt)
c
         else
            radfn41 = xvi
         endif
c
      else
         radfn41 = xvi
      endif
c
      return
c
      end
c**********************************************************************
      subroutine contnm_41(jrad,pave,tave,wk,wbroad,nmol,v1abs,v2abs,
     *dvabs,nptabs,absrb)
      INCLUDE "avhrr.inc"
c
      implicit double precision (v)
c
c     subroutine contnm contains the continuum data
c     which is interpolated into the array absrb
c
c
      dimension absrb(*)
      dimension wk(60)
      dimension ch(2050),csh2o(2050),cfh2o(2050)
      dimension c(2050),c0(2050),c1(2050),c2(2050)
      dimension xfac(0:50)
c
      data pi /3.1415927/
      data planck / 6.626176e-27 /,boltz / 1.380662e-16 /,
     *     clight / 2.99792458e10 /,avog / 6.022045e23 /
      data p0 / 1013. /,t0 / 296. /
      data xlosmt / 2.68675e+19 /
c
c     these are self-continuum modification factors from 700-1200 cm-1
c
      data (xfac(i),i=0,50)/
     1    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     2    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     3    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     4    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     5    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     6    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     7    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     8    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     9    1.00000,1.00000,1.00000/
c
      radcn1=2.*planck*clight*clight*1.e-07
      radcn2=planck*clight/boltz
c
c
c     assign sccs version number to module
c
cbm      hvrcnt = '3.3'
c
      rhoave = (pave/p0)*(t0/tave)
      xkt = tave/radcn2
      wtot = wbroad
      do 10 m = 1, nmol
         wtot = wtot+wk(m)
   10 continue
      wtot = wtot*1.e-20
      patm = pave/p0
c
c     ********    water    ********
c
      call sl296_41(v1c,v2c,dvc,nptc,c0,v1abs,v2abs,dvabs,nptabs)
      call sl260_41(v1c,v2c,dvc,nptc,c1,v1abs,v2abs,dvabs,nptabs)
      call frn296_41(v1c,v2c,dvc,nptc,c2,v1abs,v2abs,dvabs,nptabs)
c
      w1 = wk(1)*1.0e-20
c
c     the factor of 1.e-20 is handled this way to avoid underflows
c
      ph2o = patm*(w1/wtot)
      rh2o = ph2o*(t0/tave)
      pfrgn = patm-ph2o
      rfrgn = pfrgn*(t0/tave)
      xkt = tave/radcn2
      tfac = (tave-t0)/(260.-t0)
c
c--------------------------------------------------------------------
c                             self
 
      alpha2 = 200.**2
c
      alphs2= 120.**2
      betas = 5.e-06
      v0s=1310.
      factrs= 0.15
c
c--------------------------------------------------------------------
c                             foreign
      hwsqf= 330.**2
      betaf = 8.  e-11
      v0f =1130.
      factrf = 0.97
c
      v0f2 =1900.
      hwsqf2 = 150.**2
      beta2 = 3.e-06
c
c--------------------------------------------------------------------
c
      v1h=v1c
      dvh=dvc
      npth=nptc
c
      do 20 j = 1, nptc
         vj = v1c+dvc*float(j-1)
         vs2 = (vj-v0s)**2
         sh2o = 0.
         if (c0(j).gt.0.) then
            sh2o = c0(j)*(c1(j)/c0(j))**tfac
c
         sfac = 1.
         if (vj.ge.700. .and.  vj.le.1200.) then
            jfac = (vj-700.)/10. + 0.00001
            sfac = xfac(jfac)
         endif
c
c     correction to self continuum (1 sept 85); factor of 0.78 at 1000
c                             and  .......
c
      sh2o = sfac * sh2o*(1.-0.2333*(alpha2/((vj-1050.)**2+alpha2))) *
     *                  (1.-factrs*(alphs2/(vs2+(betas*vs2**2)+alphs2)))
         endif
c
c     correction to foreign continuum
c
        vf2 = (vj-v0f)**2
        vf6 = vf2 * vf2 * vf2
        fscal  = (1.-factrf*(hwsqf/(vf2+(betaf*vf6)+hwsqf)))
        vf2 = (vj-v0f2)**2
        vf4 = vf2*vf2
        fscal = fscal* (1.- 0.6*(hwsqf2/(vf2 + beta2*vf4 + hwsqf2)))
c
        c2(j)=c2(j)*fscal
c
         c(j) = w1*(sh2o*rh2o+c2(j)*rfrgn)
c
      csh2o(j)=1.e-20 * sh2o
      cfh2o(j)=1.e-20 * c2(j)
c
c     radiation field
c
         if (jrad.eq.1) c(j) = c(j)*radfn41(vj,xkt)
   20 continue
      call xint_41(v1c,v2c,dvc,c,1.0,v1abs,dvabs,absrb,1,nptabs)
c
      return
c
      end
c**********************************************************************
      subroutine sl296_41(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
      INCLUDE "avhrr.inc"
c
      implicit double precision (v)
c
      dimension s(305)
      dimension c(*)
c
c               06/28/82
c               units of (cm**3/mol) * 1.e-20
c
      data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/
c
      data (s(i),i=1,2)/
     *     1.1109e-01 ,1.0573e-01/
      data (s(i),i=3,52)/
     *     1.0162e-01, 1.0573e-01, 1.1109e-01, 1.2574e-01, 1.3499e-01,
     *     1.4327e-01, 1.5065e-01, 1.5164e-01, 1.5022e-01, 1.3677e-01,
     *     1.3115e-01, 1.2253e-01, 1.1271e-01, 1.0070e-01, 8.7495e-02,
     *     8.0118e-02, 6.9940e-02, 6.2034e-02, 5.6051e-02, 4.7663e-02,
     *     4.2450e-02, 3.6690e-02, 3.3441e-02, 3.0711e-02, 2.5205e-02,
     *     2.2113e-02, 1.8880e-02, 1.6653e-02, 1.4626e-02, 1.2065e-02,
     *     1.0709e-02, 9.1783e-03, 7.7274e-03, 6.7302e-03, 5.6164e-03,
     *     4.9089e-03, 4.1497e-03, 3.5823e-03, 3.1124e-03, 2.6414e-03,
     *     2.3167e-03, 2.0156e-03, 1.7829e-03, 1.5666e-03, 1.3928e-03,
     *     1.2338e-03, 1.0932e-03, 9.7939e-04, 8.8241e-04, 7.9173e-04/
      data (s(i),i=53,102)/
     *     7.1296e-04, 6.4179e-04, 5.8031e-04, 5.2647e-04, 4.7762e-04,
     *     4.3349e-04, 3.9355e-04, 3.5887e-04, 3.2723e-04, 2.9919e-04,
     *     2.7363e-04, 2.5013e-04, 2.2876e-04, 2.0924e-04, 1.9193e-04,
     *     1.7618e-04, 1.6188e-04, 1.4891e-04, 1.3717e-04, 1.2647e-04,
     *     1.1671e-04, 1.0786e-04, 9.9785e-05, 9.2350e-05, 8.5539e-05,
     *     7.9377e-05, 7.3781e-05, 6.8677e-05, 6.3993e-05, 5.9705e-05,
     *     5.5788e-05, 5.2196e-05, 4.8899e-05, 4.5865e-05, 4.3079e-05,
     *     4.0526e-05, 3.8182e-05, 3.6025e-05, 3.4038e-05, 3.2203e-05,
     *     3.0511e-05, 2.8949e-05, 2.7505e-05, 2.6170e-05, 2.4933e-05,
     *     2.3786e-05, 2.2722e-05, 2.1736e-05, 2.0819e-05, 1.9968e-05/
      data (s(i),i=103,152)/
     *     1.9178e-05, 1.8442e-05, 1.7760e-05, 1.7127e-05, 1.6541e-05,
     *     1.5997e-05, 1.5495e-05, 1.5034e-05, 1.4614e-05, 1.4230e-05,
     *     1.3883e-05, 1.3578e-05, 1.3304e-05, 1.3069e-05, 1.2876e-05,
     *     1.2732e-05, 1.2626e-05, 1.2556e-05, 1.2544e-05, 1.2604e-05,
     *     1.2719e-05, 1.2883e-05, 1.3164e-05, 1.3581e-05, 1.4187e-05,
     *     1.4866e-05, 1.5669e-05, 1.6717e-05, 1.8148e-05, 2.0268e-05,
     *     2.2456e-05, 2.5582e-05, 2.9183e-05, 3.3612e-05, 3.9996e-05,
     *     4.6829e-05, 5.5055e-05, 6.5897e-05, 7.5360e-05, 8.7213e-05,
     *     1.0046e-04, 1.1496e-04, 1.2943e-04, 1.5049e-04, 1.6973e-04,
     *     1.8711e-04, 2.0286e-04, 2.2823e-04, 2.6780e-04, 2.8766e-04/
      data (s(i),i=153,202)/
     *     3.1164e-04, 3.3640e-04, 3.6884e-04, 3.9159e-04, 3.8712e-04,
     *     3.7433e-04, 3.4503e-04, 3.1003e-04, 2.8027e-04, 2.5253e-04,
     *     2.3408e-04, 2.2836e-04, 2.4442e-04, 2.7521e-04, 2.9048e-04,
     *     3.0489e-04, 3.2646e-04, 3.3880e-04, 3.3492e-04, 3.0987e-04,
     *     2.9482e-04, 2.8711e-04, 2.6068e-04, 2.2683e-04, 1.9996e-04,
     *     1.7788e-04, 1.6101e-04, 1.3911e-04, 1.2013e-04, 1.0544e-04,
     *     9.4224e-05, 8.1256e-05, 7.3667e-05, 6.2233e-05, 5.5906e-05,
     *     5.1619e-05, 4.5140e-05, 4.0273e-05, 3.3268e-05, 3.0258e-05,
     *     2.6440e-05, 2.3103e-05, 2.0749e-05, 1.8258e-05, 1.6459e-05,
     *     1.4097e-05, 1.2052e-05, 1.0759e-05, 9.1400e-06, 8.1432e-06/
      data (s(i),i=203,252)/
     *     7.1460e-06, 6.4006e-06, 5.6995e-06, 4.9372e-06, 4.4455e-06,
     *     3.9033e-06, 3.4740e-06, 3.1269e-06, 2.8059e-06, 2.5558e-06,
     *     2.2919e-06, 2.0846e-06, 1.8983e-06, 1.7329e-06, 1.5929e-06,
     *     1.4631e-06, 1.3513e-06, 1.2461e-06, 1.1519e-06, 1.0682e-06,
     *     9.9256e-07, 9.2505e-07, 8.6367e-07, 8.0857e-07, 7.5674e-07,
     *     7.0934e-07, 6.6580e-07, 6.2580e-07, 5.8853e-07, 5.5333e-07,
     *     5.2143e-07, 4.9169e-07, 4.6431e-07, 4.3898e-07, 4.1564e-07,
     *     3.9405e-07, 3.7403e-07, 3.5544e-07, 3.3819e-07, 3.2212e-07,
     *     3.0714e-07, 2.9313e-07, 2.8003e-07, 2.6777e-07, 2.5628e-07,
     *     2.4551e-07, 2.3540e-07, 2.2591e-07, 2.1701e-07, 2.0866e-07/
      data (s(i),i=253,302)/
     *     2.0082e-07, 1.9349e-07, 1.8665e-07, 1.8027e-07, 1.7439e-07,
     *     1.6894e-07, 1.6400e-07, 1.5953e-07, 1.5557e-07, 1.5195e-07,
     *     1.4888e-07, 1.4603e-07, 1.4337e-07, 1.4093e-07, 1.3828e-07,
     *     1.3569e-07, 1.3270e-07, 1.2984e-07, 1.2714e-07, 1.2541e-07,
     *     1.2399e-07, 1.2102e-07, 1.1878e-07, 1.1728e-07, 1.1644e-07,
     *     1.1491e-07, 1.1305e-07, 1.1235e-07, 1.1228e-07, 1.1224e-07,
     *     1.1191e-07, 1.1151e-07, 1.1098e-07, 1.1068e-07, 1.1109e-07,
     *     1.1213e-07, 1.1431e-07, 1.1826e-07, 1.2322e-07, 1.3025e-07,
     *     1.4066e-07, 1.5657e-07, 1.7214e-07, 1.9449e-07, 2.2662e-07,
     *     2.6953e-07, 3.1723e-07, 3.7028e-07, 4.4482e-07, 5.3852e-07/
      data (s(i),i=303,305)/
     *     6.2639e-07, 7.2175e-07, 7.7626e-07/
c
      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc
c
      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1
c
      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)
      do 10 j = 1, nptc
         i = i1+j
         c(j) = 0.
         if ((i.lt.1).or.(i.gt.npts)) go to 10
         c(j) = s(i)
   10 continue
c
      return
c
      end
c**********************************************************************
      subroutine sl260_41(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
      INCLUDE "avhrr.inc"
c
      implicit double precision (v)
c
      dimension s(305)
      dimension c(*)
c
c               06/28/82
c               units of (cm**3/mol) * 1.e-20
c
       data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/
c
c
      data (s(i),i=1,2)/
     *     1.7750e-01, 1.7045e-01/
      data (s(i),i=3,52)/
     *     1.6457e-01, 1.7045e-01, 1.7750e-01, 2.0036e-01, 2.1347e-01,
     *     2.2454e-01, 2.3428e-01, 2.3399e-01, 2.3022e-01, 2.0724e-01,
     *     1.9712e-01, 1.8317e-01, 1.6724e-01, 1.4780e-01, 1.2757e-01,
     *     1.1626e-01, 1.0098e-01, 8.9033e-02, 7.9770e-02, 6.7416e-02,
     *     5.9588e-02, 5.1117e-02, 4.6218e-02, 4.2179e-02, 3.4372e-02,
     *     2.9863e-02, 2.5252e-02, 2.2075e-02, 1.9209e-02, 1.5816e-02,
     *     1.3932e-02, 1.1943e-02, 1.0079e-02, 8.7667e-03, 7.4094e-03,
     *     6.4967e-03, 5.5711e-03, 4.8444e-03, 4.2552e-03, 3.6953e-03,
     *     3.2824e-03, 2.9124e-03, 2.6102e-03, 2.3370e-03, 2.1100e-03,
     *     1.9008e-03, 1.7145e-03, 1.5573e-03, 1.4206e-03, 1.2931e-03/
      data (s(i),i=53,102)/
     *     1.1803e-03, 1.0774e-03, 9.8616e-04, 9.0496e-04, 8.3071e-04,
     *     7.6319e-04, 7.0149e-04, 6.4637e-04, 5.9566e-04, 5.4987e-04,
     *     5.0768e-04, 4.6880e-04, 4.3317e-04, 4.0037e-04, 3.7064e-04,
     *     3.4325e-04, 3.1809e-04, 2.9501e-04, 2.7382e-04, 2.5430e-04,
     *     2.3630e-04, 2.1977e-04, 2.0452e-04, 1.9042e-04, 1.7740e-04,
     *     1.6544e-04, 1.5442e-04, 1.4425e-04, 1.3486e-04, 1.2618e-04,
     *     1.1817e-04, 1.1076e-04, 1.0391e-04, 9.7563e-05, 9.1696e-05,
     *     8.6272e-05, 8.1253e-05, 7.6607e-05, 7.2302e-05, 6.8311e-05,
     *     6.4613e-05, 6.1183e-05, 5.8001e-05, 5.5048e-05, 5.2307e-05,
     *     4.9761e-05, 4.7395e-05, 4.5197e-05, 4.3155e-05, 4.1256e-05/
      data (s(i),i=103,152)/
     *     3.9491e-05, 3.7849e-05, 3.6324e-05, 3.4908e-05, 3.3594e-05,
     *     3.2374e-05, 3.1244e-05, 3.0201e-05, 2.9240e-05, 2.8356e-05,
     *     2.7547e-05, 2.6814e-05, 2.6147e-05, 2.5551e-05, 2.5029e-05,
     *     2.4582e-05, 2.4203e-05, 2.3891e-05, 2.3663e-05, 2.3531e-05,
     *     2.3483e-05, 2.3516e-05, 2.3694e-05, 2.4032e-05, 2.4579e-05,
     *     2.5234e-05, 2.6032e-05, 2.7119e-05, 2.8631e-05, 3.0848e-05,
     *     3.3262e-05, 3.6635e-05, 4.0732e-05, 4.5923e-05, 5.3373e-05,
     *     6.1875e-05, 7.2031e-05, 8.5980e-05, 9.8642e-05, 1.1469e-04,
     *     1.3327e-04, 1.5390e-04, 1.7513e-04, 2.0665e-04, 2.3609e-04,
     *     2.6220e-04, 2.8677e-04, 3.2590e-04, 3.8624e-04, 4.1570e-04/
      data (s(i),i=153,202)/
     *     4.5207e-04, 4.9336e-04, 5.4500e-04, 5.8258e-04, 5.8086e-04,
     *     5.6977e-04, 5.3085e-04, 4.8020e-04, 4.3915e-04, 4.0343e-04,
     *     3.7853e-04, 3.7025e-04, 3.9637e-04, 4.4675e-04, 4.7072e-04,
     *     4.9022e-04, 5.2076e-04, 5.3676e-04, 5.2755e-04, 4.8244e-04,
     *     4.5473e-04, 4.3952e-04, 3.9614e-04, 3.4086e-04, 2.9733e-04,
     *     2.6367e-04, 2.3767e-04, 2.0427e-04, 1.7595e-04, 1.5493e-04,
     *     1.3851e-04, 1.1874e-04, 1.0735e-04, 9.0490e-05, 8.1149e-05,
     *     7.4788e-05, 6.5438e-05, 5.8248e-05, 4.8076e-05, 4.3488e-05,
     *     3.7856e-05, 3.3034e-05, 2.9592e-05, 2.6088e-05, 2.3497e-05,
     *     2.0279e-05, 1.7526e-05, 1.5714e-05, 1.3553e-05, 1.2145e-05/
      data (s(i),i=203,252)/
     *     1.0802e-05, 9.7681e-06, 8.8196e-06, 7.8291e-06, 7.1335e-06,
     *     6.4234e-06, 5.8391e-06, 5.3532e-06, 4.9079e-06, 4.5378e-06,
     *     4.1716e-06, 3.8649e-06, 3.5893e-06, 3.3406e-06, 3.1199e-06,
     *     2.9172e-06, 2.7348e-06, 2.5644e-06, 2.4086e-06, 2.2664e-06,
     *     2.1359e-06, 2.0159e-06, 1.9051e-06, 1.8031e-06, 1.7074e-06,
     *     1.6185e-06, 1.5356e-06, 1.4584e-06, 1.3861e-06, 1.3179e-06,
     *     1.2545e-06, 1.1951e-06, 1.1395e-06, 1.0873e-06, 1.0384e-06,
     *     9.9250e-07, 9.4935e-07, 9.0873e-07, 8.7050e-07, 8.3446e-07,
     *     8.0046e-07, 7.6834e-07, 7.3800e-07, 7.0931e-07, 6.8217e-07,
     *     6.5648e-07, 6.3214e-07, 6.0909e-07, 5.8725e-07, 5.6655e-07/
      data (s(i),i=253,302)/
     *     5.4693e-07, 5.2835e-07, 5.1077e-07, 4.9416e-07, 4.7853e-07,
     *     4.6381e-07, 4.5007e-07, 4.3728e-07, 4.2550e-07, 4.1450e-07,
     *     4.0459e-07, 3.9532e-07, 3.8662e-07, 3.7855e-07, 3.7041e-07,
     *     3.6254e-07, 3.5420e-07, 3.4617e-07, 3.3838e-07, 3.3212e-07,
     *     3.2655e-07, 3.1865e-07, 3.1203e-07, 3.0670e-07, 3.0252e-07,
     *     2.9749e-07, 2.9184e-07, 2.8795e-07, 2.8501e-07, 2.8202e-07,
     *     2.7856e-07, 2.7509e-07, 2.7152e-07, 2.6844e-07, 2.6642e-07,
     *     2.6548e-07, 2.6617e-07, 2.6916e-07, 2.7372e-07, 2.8094e-07,
     *     2.9236e-07, 3.1035e-07, 3.2854e-07, 3.5481e-07, 3.9377e-07,
     *     4.4692e-07, 5.0761e-07, 5.7715e-07, 6.7725e-07, 8.0668e-07/
      data (s(i),i=303,305)/
     *     9.3716e-07, 1.0797e-06, 1.1689e-06/
c
      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc
c
      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1
c
      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)
      do 10 j = 1, nptc
         i = i1+j
         c(j) = 0.
         if ((i.lt.1).or.(i.gt.npts)) go to 10
         c(j) = s(i)
   10 continue
c
      return
c
      end
c**********************************************************************
      subroutine frn296_41(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
      INCLUDE "avhrr.inc"
c
      implicit double precision (v)
c
      dimension s(305)
      dimension c(*)
c
c               06/28/82
c               units of (cm**3/mol)*1.e-20
c
       data v1s,v2s,dvs,npts / -20.0, 3020.0, 10.0, 305/
c
      data (s(i),i=1,2)/
     *     1.2859e-02, 1.1715e-02/
      data (s(i),i=3,52)/
     *     1.1038e-02, 1.1715e-02, 1.2859e-02, 1.5326e-02, 1.6999e-02,
     *     1.8321e-02, 1.9402e-02, 1.9570e-02, 1.9432e-02, 1.7572e-02,
     *     1.6760e-02, 1.5480e-02, 1.3984e-02, 1.2266e-02, 1.0467e-02,
     *     9.4526e-03, 8.0485e-03, 6.9484e-03, 6.1416e-03, 5.0941e-03,
     *     4.4836e-03, 3.8133e-03, 3.4608e-03, 3.1487e-03, 2.4555e-03,
     *     2.0977e-03, 1.7266e-03, 1.4920e-03, 1.2709e-03, 9.8081e-04,
     *     8.5063e-04, 6.8822e-04, 5.3809e-04, 4.4679e-04, 3.3774e-04,
     *     2.7979e-04, 2.1047e-04, 1.6511e-04, 1.2993e-04, 9.3033e-05,
     *     7.4360e-05, 5.6428e-05, 4.5442e-05, 3.4575e-05, 2.7903e-05,
     *     2.1374e-05, 1.6075e-05, 1.3022e-05, 1.0962e-05, 8.5959e-06/
      data (s(i),i=53,102)/
     *     6.9125e-06, 5.3808e-06, 4.3586e-06, 3.6394e-06, 2.9552e-06,
     *     2.3547e-06, 1.8463e-06, 1.6036e-06, 1.3483e-06, 1.1968e-06,
     *     1.0333e-06, 8.4484e-07, 6.7195e-07, 5.0947e-07, 4.2343e-07,
     *     3.4453e-07, 2.7830e-07, 2.3063e-07, 1.9951e-07, 1.7087e-07,
     *     1.4393e-07, 1.2575e-07, 1.0750e-07, 8.2325e-08, 5.7524e-08,
     *     4.4482e-08, 3.8106e-08, 3.4315e-08, 2.9422e-08, 2.5069e-08,
     *     2.2402e-08, 1.9349e-08, 1.6152e-08, 1.2208e-08, 8.9660e-09,
     *     7.1322e-09, 6.1028e-09, 5.2938e-09, 4.5350e-09, 3.4977e-09,
     *     2.9511e-09, 2.4734e-09, 2.0508e-09, 1.8507e-09, 1.6373e-09,
     *     1.5171e-09, 1.3071e-09, 1.2462e-09, 1.2148e-09, 1.2590e-09/
      data (s(i),i=103,152)/
     *     1.3153e-09, 1.3301e-09, 1.4483e-09, 1.6944e-09, 2.0559e-09,
     *     2.2954e-09, 2.6221e-09, 3.2606e-09, 4.2392e-09, 5.2171e-09,
     *     6.2553e-09, 8.2548e-09, 9.5842e-09, 1.1280e-08, 1.3628e-08,
     *     1.7635e-08, 2.1576e-08, 2.4835e-08, 3.0014e-08, 3.8485e-08,
     *     4.7440e-08, 5.5202e-08, 7.0897e-08, 9.6578e-08, 1.3976e-07,
     *     1.8391e-07, 2.3207e-07, 2.9960e-07, 4.0408e-07, 5.9260e-07,
     *     7.8487e-07, 1.0947e-06, 1.4676e-06, 1.9325e-06, 2.6587e-06,
     *     3.4534e-06, 4.4376e-06, 5.8061e-06, 7.0141e-06, 8.4937e-06,
     *     1.0186e-05, 1.2034e-05, 1.3837e-05, 1.6595e-05, 1.9259e-05,
     *     2.1620e-05, 2.3681e-05, 2.7064e-05, 3.2510e-05, 3.5460e-05/
      data (s(i),i=153,202)/
     *     3.9109e-05, 4.2891e-05, 4.7757e-05, 5.0981e-05, 5.0527e-05,
     *     4.8618e-05, 4.4001e-05, 3.7982e-05, 3.2667e-05, 2.7794e-05,
     *     2.4910e-05, 2.4375e-05, 2.7316e-05, 3.2579e-05, 3.5499e-05,
     *     3.8010e-05, 4.1353e-05, 4.3323e-05, 4.3004e-05, 3.9790e-05,
     *     3.7718e-05, 3.6360e-05, 3.2386e-05, 2.7409e-05, 2.3626e-05,
     *     2.0631e-05, 1.8371e-05, 1.5445e-05, 1.2989e-05, 1.1098e-05,
     *     9.6552e-06, 8.0649e-06, 7.2365e-06, 5.9137e-06, 5.2759e-06,
     *     4.8860e-06, 4.1321e-06, 3.5918e-06, 2.7640e-06, 2.4892e-06,
     *     2.1018e-06, 1.7848e-06, 1.5855e-06, 1.3569e-06, 1.1986e-06,
     *     9.4693e-07, 7.4097e-07, 6.3443e-07, 4.8131e-07, 4.0942e-07/
      data (s(i),i=203,252)/
     *     3.3316e-07, 2.8488e-07, 2.3461e-07, 1.7397e-07, 1.4684e-07,
     *     1.0953e-07, 8.5396e-08, 6.9261e-08, 5.4001e-08, 4.5430e-08,
     *     3.2791e-08, 2.5995e-08, 2.0225e-08, 1.5710e-08, 1.3027e-08,
     *     1.0229e-08, 8.5277e-09, 6.5249e-09, 5.0117e-09, 3.9906e-09,
     *     3.2332e-09, 2.7847e-09, 2.4570e-09, 2.3359e-09, 2.0599e-09,
     *     1.8436e-09, 1.6559e-09, 1.4910e-09, 1.2794e-09, 9.8229e-10,
     *     8.0054e-10, 6.0769e-10, 4.5646e-10, 3.3111e-10, 2.4428e-10,
     *     1.8007e-10, 1.3291e-10, 9.7974e-11, 7.8271e-11, 6.3833e-11,
     *     5.4425e-11, 4.6471e-11, 4.0209e-11, 3.5227e-11, 3.1212e-11,
     *     2.8840e-11, 2.7762e-11, 2.7935e-11, 3.2012e-11, 3.9525e-11/
      data (s(i),i=253,302)/
     *     5.0303e-11, 6.8027e-11, 9.3954e-11, 1.2986e-10, 1.8478e-10,
     *     2.5331e-10, 3.4827e-10, 4.6968e-10, 6.2380e-10, 7.9106e-10,
     *     1.0026e-09, 1.2102e-09, 1.4146e-09, 1.6154e-09, 1.7510e-09,
     *     1.8575e-09, 1.8742e-09, 1.8700e-09, 1.8582e-09, 1.9657e-09,
     *     2.1204e-09, 2.0381e-09, 2.0122e-09, 2.0436e-09, 2.1213e-09,
     *     2.0742e-09, 1.9870e-09, 2.0465e-09, 2.1556e-09, 2.2222e-09,
     *     2.1977e-09, 2.1047e-09, 1.9334e-09, 1.7357e-09, 1.5754e-09,
     *     1.4398e-09, 1.4018e-09, 1.5459e-09, 1.7576e-09, 2.1645e-09,
     *     2.9480e-09, 4.4439e-09, 5.8341e-09, 8.0757e-09, 1.1658e-08,
     *     1.6793e-08, 2.2694e-08, 2.9468e-08, 3.9278e-08, 5.2145e-08/
      data (s(i),i=303,305)/
     *     6.4378e-08, 7.7947e-08, 8.5321e-08/
c
      dvc = dvs
      v1c = v1abs-dvc
      v2c = v2abs+dvc
c
      i1 = (v1c-v1s)/dvs
      if (v1c.lt.v1s) i1 = i1-1
c
      v1c = v1s+dvs*float(i1)
      i2 = (v2c-v1s)/dvs
      nptc = i2-i1+3
      v2c = v1c+dvs*float(nptc-1)
      do 10 j = 1, nptc
         i = i1+j
         c(j) = 0.
         if ((i.ge.1).and.(i.le.npts)) then
            c(j) = s(i)
         endif
   10 continue
c
      return
c
      end
