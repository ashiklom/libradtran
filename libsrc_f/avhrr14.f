c     Response function factor = 0.78896
      subroutine avhrr14(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(3,mm)
      real wght(3)
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
      omega=16950.0
      delw=700.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_14(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_14(back,omega,t0(1))
      call flux_14(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_14(ni,mlv,u,delw,f,p,t,tau)
      call cko3_14(no3,mlv,uo3,delw,fo3,p,t,tauo3)
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
      do i1=1,3
       do m=1,nlyr
        od(1+(i1-1),m) =
     *tau(i1,m)
        wght(1+(i1-1)) =
     *f(i1)
       enddo
      enddo
      ntau=3
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_14(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_14(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=3
      f(1)=0.914321
      f(2)=0.054775
      f(3)=0.030904
      k(1)=0.0024289518
      do i=2,ni
        k(i)=8.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.1013,0.1013,0.1013,0.1013,0.1013,0.1013,0.1013,0.1013,0.1013,
     *0.1235,0.1569,0.2061,0.2763,0.3743,0.5061,0.6784,0.8855,1.1379,
     *1.4178,
     *5.336e-04,5.336e-04,5.336e-04,5.336e-04,5.336e-04,5.336e-04,
     *5.336e-04,5.336e-04,5.336e-04,4.715e-04,3.877e-04,2.794e-04,
     *1.297e-04,-6.213e-05,-2.831e-04,-5.150e-04,-6.824e-04,-8.774e-04,
     *-1.020e-03,
     *-7.031e-07,-7.031e-07,-7.031e-07,-7.031e-07,-7.031e-07,-7.031e-07,
     *-7.031e-07,-7.031e-07,-7.031e-07,-8.500e-07,-1.156e-06,-1.422e-06,
     *-1.387e-06,-1.353e-06,-9.406e-07,-7.187e-07,3.781e-07,8.406e-07,
     *-1.594e-07/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.5929,0.5929,0.5929,0.5929,0.5929,0.5929,0.5929,0.5929,0.5929,
     *0.6013,0.6154,0.6377,0.6729,0.7252,0.7952,0.8817,0.9649,1.0435,
     *1.0898,
     *1.616e-03,1.616e-03,1.616e-03,1.616e-03,1.616e-03,1.616e-03,
     *1.616e-03,1.616e-03,1.616e-03,1.529e-03,1.392e-03,1.196e-03,
     *9.329e-04,6.607e-04,3.900e-04,1.035e-04,-9.275e-05,-2.953e-04,
     *-4.020e-04,
     *-1.747e-06,-1.747e-06,-1.747e-06,-1.747e-06,-1.747e-06,-1.747e-06,
     *-1.747e-06,-1.747e-06,-1.747e-06,-7.375e-07,1.281e-07,1.787e-06,
     *3.016e-06,3.037e-06,3.469e-06,2.931e-06,1.950e-06,7.625e-07,
     *1.252e-08/
      data ( (coefk(3,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.5250,1.5250,1.5250,1.5250,1.5250,1.5250,1.5250,1.5250,1.5250,
     *1.5128,1.4948,1.4677,1.4271,1.3719,1.2951,1.1977,1.0700,0.9381,
     *0.7974,
     *-1.931e-03,-1.931e-03,-1.931e-03,-1.931e-03,-1.931e-03,-1.931e-03,
     *-1.931e-03,-1.931e-03,-1.931e-03,-1.891e-03,-1.820e-03,-1.716e-03,
     *-1.613e-03,-1.465e-03,-1.284e-03,-1.130e-03,-9.167e-04,-8.091e-04,
     *-7.313e-04,
     *3.766e-06,3.766e-06,3.766e-06,3.766e-06,3.766e-06,3.766e-06,
     *3.766e-06,3.766e-06,3.766e-06,4.028e-06,3.250e-06,3.303e-06,
     *3.447e-06,2.297e-06,2.563e-06,1.969e-06,2.044e-06,2.366e-06,
     *2.619e-06/
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
      subroutine cko3_14(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=4.55e-21*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
