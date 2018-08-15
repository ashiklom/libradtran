c     Response function factor = 0.89060
      subroutine avhrr13(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(1,mm)
      real wght(1)
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
      omega=16300.0
      delw=600.0
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_13(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uo3(m)=fac(m)*(ux(m)+ux(m+1))/2.0
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_13(back,omega,t0(1))
      call flux_13(bsun,omega,tsun)
      bsun=bsun/fsun*sconst
      call ck_13(ni,mlv,u,delw,f,p,t,tau)
      call cko3_13(no3,mlv,uo3,delw,fo3,p,t,tauo3)
c     **fo3=1.0000**
      do 10 m=1,nlyr
        do i=1,ni
          tau(i,m)=tau(i,m)+tauo3(1,m)
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
      do i1=1,1
       do m=1,nlyr
        od(1+(i1-1),m) =
     *tau(i1,m)
        wght(1+(i1-1)) =
     *f(i1)
       enddo
      enddo
      ntau=1
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_13(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_13(ni,mlv,u,delw,f,p,t,tau)
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
      subroutine cko3_13(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm)
      ni=1
      f(1)=1.000000
      k(1)=4.35e-21*6.023e+23/47.9982
      do i=1,ni
        do m=1,mlv
          tau(i,m)=k(i)*u(m)
        end do
      end do
      return
      end
