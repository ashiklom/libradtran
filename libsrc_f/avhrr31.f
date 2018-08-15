      subroutine avhrr31(z0, p0, t0, u0, ux, nlev,
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
      omega=2525.0
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
        call flux_31(bnu(m),omega,t(m))
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
      call flux_31(back,omega,t0(1))
      call flux_31(bsun,omegas,tsun)
      bsun=bsun/fsun*sconst
      call ck_31(ni,mlv,u,delw,f,p,t,tau)
      call ckch4_31(nch4,mlv,uch4,delw,fch4,p,t,tauch4)
      call ckn2o_31(nn2o,mlv,un2o,delw,fn2o,p,t,taun2o)
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
      subroutine flux_31(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_31(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.985956
      f(2)=0.014044
      k(1)=0.0010079803
      do i=2,ni
        k(i)=48.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.1841,0.1841,0.1841,0.1841,0.1841,0.1841,0.1841,0.1841,0.1841,
     *0.2224,0.2757,0.3466,0.4364,0.5444,0.6644,0.7943,0.9301,1.0653,
     *1.1844,
     *4.106e-03,4.106e-03,4.106e-03,4.106e-03,4.106e-03,4.106e-03,
     *4.106e-03,4.106e-03,4.106e-03,4.731e-03,5.649e-03,6.892e-03,
     *8.545e-03,1.057e-02,1.288e-02,1.544e-02,1.820e-02,2.092e-02,
     *2.336e-02,
     *3.184e-05,3.184e-05,3.184e-05,3.184e-05,3.184e-05,3.184e-05,
     *3.184e-05,3.184e-05,3.184e-05,3.496e-05,3.990e-05,4.726e-05,
     *5.824e-05,7.242e-05,8.933e-05,1.072e-04,1.262e-04,1.418e-04,
     *1.572e-04/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *2.1845,2.1845,2.1845,2.1845,2.1845,2.1845,2.1845,2.1845,2.1845,
     *2.1297,2.0531,1.9491,1.8168,1.6580,1.4840,1.2939,1.0983,0.8994,
     *0.7275,
     *4.325e-02,4.325e-02,4.325e-02,4.325e-02,4.325e-02,4.325e-02,
     *4.325e-02,4.325e-02,4.325e-02,4.235e-02,4.102e-02,3.918e-02,
     *3.682e-02,3.381e-02,3.042e-02,2.665e-02,2.263e-02,1.854e-02,
     *1.493e-02,
     *2.785e-04,2.785e-04,2.785e-04,2.785e-04,2.785e-04,2.785e-04,
     *2.785e-04,2.785e-04,2.785e-04,2.734e-04,2.656e-04,2.547e-04,
     *2.399e-04,2.182e-04,1.924e-04,1.662e-04,1.373e-04,1.128e-04,
     *8.950e-05/
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
      subroutine ckch4_31(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.970893
      f(2)=0.029107
      k(1)=3.05253987
      do i=2,ni
        k(i)=36.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.1450,0.1450,0.1450,0.1450,0.1450,0.1450,0.1450,0.1450,0.1450,
     *0.1801,0.2288,0.2943,0.3799,0.4855,0.6084,0.7492,0.9109,1.0927,
     *1.2683,
     *6.250e-04,6.250e-04,6.250e-04,6.250e-04,6.250e-04,6.250e-04,
     *6.250e-04,6.250e-04,6.250e-04,5.906e-04,5.633e-04,5.335e-04,
     *5.476e-04,5.810e-04,6.444e-04,6.384e-04,4.586e-04,2.983e-04,
     *3.597e-04,
     *4.375e-07,4.375e-07,4.375e-07,4.375e-07,4.375e-07,4.375e-07,
     *4.375e-07,4.375e-07,4.375e-07,-3.437e-08,-7.063e-07,-1.044e-06,
     *-1.959e-06,-3.237e-06,-3.853e-06,-3.622e-06,-2.397e-06,-2.225e-06,
     *-5.350e-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.7740,1.7740,1.7740,1.7740,1.7740,1.7740,1.7740,1.7740,1.7740,
     *1.7740,1.7281,1.6673,1.5870,1.4886,1.3722,1.2365,1.0837,0.9117,
     *0.7453,
     *1.727e-04,1.727e-04,1.727e-04,1.727e-04,1.727e-04,1.727e-04,
     *1.727e-04,1.727e-04,1.727e-04,1.727e-04,2.139e-04,2.355e-04,
     *2.547e-04,2.380e-04,1.720e-04,1.867e-04,3.383e-04,5.025e-04,
     *4.486e-04,
     *-9.919e-06,-9.919e-06,-9.919e-06,-9.919e-06,-9.919e-06,-9.919e-06,
     *-9.919e-06,-9.919e-06,-9.919e-06,-9.919e-06,-9.134e-06,-8.131e-06,
     *-7.256e-06,-7.169e-06,-6.625e-06,-6.300e-06,-7.937e-06,-7.394e-06,
     *-5.078e-06/
	data stp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     *	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     *	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
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
      subroutine ckn2o_31(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
      f(1)=0.696936
      f(2)=0.237465
      f(3)=0.065599
      k(1)=12.5009444
      do i=2,ni
        k(i)=9.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.02212,0.02212,0.02212,0.02212,0.02212,0.02212,
     *0.02212,0.02212,0.02212,0.04383,0.07673,0.12265,
     *0.18576,0.28034,0.41269,0.59252,0.84899,1.18822,
     * 1.59001,
     * 1.895E-04, 1.895E-04, 1.895E-04, 1.895E-04, 1.895E-04, 1.895E-04,
     * 1.895E-04, 1.895E-04, 1.895E-04, 2.763E-04, 3.986E-04, 6.550E-04,
     * 9.229E-04, 1.417E-03, 2.067E-03, 3.018E-03, 4.352E-03, 6.331E-03,
     * 8.956E-03,
     * 9.250E-07, 9.250E-07, 9.250E-07, 9.250E-07, 9.250E-07, 9.250E-07,
     * 9.250E-07, 9.250E-07, 9.250E-07, 1.269E-06, 1.572E-06, 1.744E-06,
     * 2.541E-06, 3.822E-06, 4.906E-06, 5.531E-06, 2.416E-06,-1.138E-06,
     *-3.788E-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.06808,0.06808,0.06808,0.06808,0.06808,0.06808,
     *0.06808,0.06808,0.06808,0.09098,0.12207,0.16351,
     *0.22360,0.31496,0.44511,0.62763,0.86851,1.15743,
     * 1.45730,
     * 4.320E-04, 4.320E-04, 4.320E-04, 4.320E-04, 4.320E-04, 4.320E-04,
     * 4.320E-04, 4.320E-04, 4.320E-04, 5.049E-04, 5.926E-04, 6.813E-04,
     * 8.382E-04, 1.047E-03, 1.273E-03, 1.334E-03, 1.198E-03, 9.530E-04,
     * 6.486E-04,
     * 1.656E-06, 1.656E-06, 1.656E-06, 1.656E-06, 1.656E-06, 1.656E-06,
     * 1.656E-06, 1.656E-06, 1.656E-06, 2.084E-06, 2.759E-06, 3.475E-06,
     * 3.919E-06, 3.750E-06, 5.825E-06, 6.909E-06, 6.906E-06, 6.981E-06,
     * 5.866E-06/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.42388,1.42388,1.42388,1.42388,1.42388,1.42388,
     *1.42388,1.42388,1.42388,1.38818,1.35694,1.33104,
     *1.30827,1.28260,1.24359,1.18016,1.07683,0.92628,
     * 0.75171,
     * 2.339E-03, 2.339E-03, 2.339E-03, 2.339E-03, 2.339E-03, 2.339E-03,
     * 2.339E-03, 2.339E-03, 2.339E-03, 2.242E-03, 2.110E-03, 1.967E-03,
     * 1.837E-03, 1.689E-03, 1.565E-03, 1.493E-03, 1.435E-03, 1.323E-03,
     * 1.103E-03,
     * 7.500E-07, 7.500E-07, 7.500E-07, 7.500E-07, 7.500E-07, 7.500E-07,
     * 7.500E-07, 7.500E-07, 7.500E-07, 9.687E-07, 9.687E-07, 4.594E-07,
     * 2.000E-07,-6.000E-07,-1.728E-06,-2.272E-06,-3.497E-06,-2.519E-06,
     *-1.641E-06/
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
