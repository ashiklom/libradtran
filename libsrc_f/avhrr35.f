c     Response function factor = 0.62487
      subroutine avhrr35(z0, p0, t0, u0, ux, nlev,
     *co2, n2o, f11, f12, ch4,
     *od, wght, ntau)
      INCLUDE "avhrr.inc"
      real co2, n2o, ch4, f11, f12
      real z0(mm)
      real od(4,mm)
      real wght(4)
      integer ntau,nlev
      real wi(8),xi(8),p0(mq),p(mm),t0(mq),t(mm),bnu(mm),u0(mq),u(mm),
     *uch4(mm),tau(6,mm),tauch4(6,mm),tauf(6,6,mm),x(mm,8),y(mm,8),
     *sx(6,6,mm),sy(6,6,mm),f(6),fch4(6),fac(mm),tsx(mq),tsy(mq)
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
      omega=2817.5
      delw=85.0
      omegas=omega
      sigma=5.67e-05
      sconst=1.368e+06
      tsun=5710.0
      fsun=sigma*tsun**4
      mlv=nlyr
      do 5 m=1,nlyr
        t(m)=(t0(m)+t0(m+1))/2.0
        p(m)=(p0(m)+p0(m+1))/2.0
        call flux_35(bnu(m),omega,t(m))
        fac(m)=z0(m+1)-z0(m)
c       u(m)=fac(m)*(u0(m)+u0(m+1))/2.0*6.023e+23/2.687e+19/18.01534
        u(m)=fac(m)*(u0(m)+u0(m+1))/2.0
        uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*ch4*1e-6*fac(m)*
     *       1.0e+05/6.023e+23*16.04303
 5    continue
      do 7 m=1,nlev
        tsx(m)=0.0
        tsy(m)=0.0
 7    continue
      call flux_35(back,omega,t0(1))
      call flux_35(bsun,omegas,tsun)
      bsun=bsun/fsun*sconst
      call ck_35(ni,mlv,u,delw,f,p,t,tau)
      call ckch4_35(nch4,mlv,uch4,delw,fch4,p,t,tauch4)
      do 10 m=1,nlyr
        do i=1,ni
          do ia=1,nch4
            tauf(i,ia,m)=tau(i,m)+tauch4(ia,m)
c           tauf(i,ia,m)=tauch4(ia,m)
c           tauf(i,ia,m)=tau(i,m)
          end do
        end do
 10   continue
      do 90 i=1,ni
       do 85 ia=1,nch4
        do 30 m=nlyr,1,-1
          do 20 j=1,8
           if(m.eq.nlyr)then
c           y(m,j)=bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
            y(m,j)=bsun*0.6*exp(-tauf(i,ia,m)/xi(j))+
     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
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
c            x(m,j)=back*exp(-tauf(i,ia,m)/xi(j))+
c     *      bnu(m)*(1.0-exp(-tauf(i,ia,m)/xi(j)))
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
           tsx(m+1)=tsx(m+1)+sx(i,ia,m)*delw*f(i)*fch4(ia)
           tsy(m)=tsy(m)+sy(i,ia,m)*delw*f(i)*fch4(ia)
 80     continue
        tb=tb+back*delw*f(i)*fch4(ia)
        tsx(1)=tb
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
        do m=1,nlyr
         od(1+(i2-1)+(i1-1)*2,m) =
     *tauf(i1,i2,m)
         wght(1+(i2-1)+(i1-1)*2) =
     *f(i1) * fch4(i2)
        enddo
       enddo
      enddo
      ntau=4
 100  format(i3,3f15.7)
      end
c *********************************************************************
      subroutine flux_35(bnu,omega,t)
      INCLUDE "avhrr.inc"
      bnu=0.0
      z=3.14159
      bnu=z/1000.0*((2.*6.625e-27)*(2.997925e+10)**2
     *    *omega**3)/(exp(min(88.0,1.4388*omega/t))-1.0)
      return
      end
c *********************************************************************
      subroutine ck_35(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.935168
      f(2)=0.064832
      k(1)=0.03154426
      do i=2,ni
        k(i)=32.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0582,0.0582,0.0582,0.0582,0.0582,0.0582,0.0582,0.0582,0.0582,
     *0.0855,0.1251,0.1819,0.2593,0.3621,0.4956,0.6624,0.8760,1.1128,
     *1.3704,
     *1.243e-04,1.243e-04,1.243e-04,1.243e-04,1.243e-04,1.243e-04,
     *1.243e-04,1.243e-04,1.243e-04,1.189e-04,1.074e-04,8.775e-05,
     *5.737e-05,2.275e-05,4.413e-05,1.951e-04,4.390e-04,7.075e-04,
     *1.062e-03,
     *3.937e-07,3.937e-07,3.937e-07,3.937e-07,3.937e-07,3.937e-07,
     *3.937e-07,3.937e-07,3.937e-07,2.719e-07,1.969e-07,-2.813e-07,
     *-2.406e-07,-4.687e-07,-1.203e-06,-2.197e-06,-1.906e-06,-2.037e-06,
     *-2.606e-06/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.3796,1.3796,1.3796,1.3796,1.3796,1.3796,1.3796,1.3796,1.3796,
     *1.3689,1.3503,1.3252,1.2899,1.2433,1.1833,1.1077,1.0378,0.9315,
     *0.8047,
     *1.695e-04,1.695e-04,1.695e-04,1.695e-04,1.695e-04,1.695e-04,
     *1.695e-04,1.695e-04,1.695e-04,1.767e-04,1.914e-04,2.091e-04,
     *2.524e-04,2.587e-04,2.546e-04,1.834e-04,7.612e-05,-1.219e-04,
     *-3.417e-04,
     *1.737e-06,1.737e-06,1.737e-06,1.737e-06,1.737e-06,1.737e-06,
     *1.737e-06,1.737e-06,1.737e-06,7.437e-07,1.947e-06,1.478e-06,
     *1.447e-06,1.175e-06,1.616e-06,2.541e-06,1.841e-06,8.656e-07,
     *2.537e-06/
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
      subroutine ckch4_35(ni,mlv,u,delw,f,p,t,tau)
      INCLUDE "avhrr.inc"
      real f(6),k(mm),tau(6,mm),delw,p(mm),t(mm),u(mm),coefk(6,3,19),
     *     stp(19),fkg(6,mm),pmb(mm)
      ni=2
      f(1)=0.862979
      f(2)=0.137021
      k(1)=18.91527
      do i=2,ni
        k(i)=16.0*k(i-1)
      end do
      data ( (coefk(1,jt,jp), jp = 1, 19), jt = 1, 3)/
     *0.0927,0.0927,0.0927,0.0927,0.0927,0.0927,0.0927,0.0927,0.0927,
     *0.1244,0.1683,0.2295,0.3121,0.4199,0.5579,0.7275,0.9079,1.1418,
     *1.3839,
     *4.329e-04,4.329e-04,4.329e-04,4.329e-04,4.329e-04,4.329e-04,
     *4.329e-04,4.329e-04,4.329e-04,4.309e-04,4.271e-04,4.199e-04,
     *4.284e-04,4.320e-04,4.501e-04,4.654e-04,5.425e-04,5.360e-04,
     *5.393e-04,
     *1.853e-06,1.853e-06,1.853e-06,1.853e-06,1.853e-06,1.853e-06,
     *1.853e-06,1.853e-06,1.853e-06,1.884e-06,2.278e-06,2.497e-06,
     *2.597e-06,3.287e-06,4.484e-06,7.666e-06,1.052e-05,1.502e-05,
     *2.008e-05/
      data ( (coefk(2,jt,jp), jp = 1, 19), jt = 1, 3)/
     *1.4115,1.4115,1.4115,1.4115,1.4115,1.4115,1.4115,1.4115,1.4115,
     *1.3999,1.3820,1.3584,1.3265,1.2833,1.2280,1.1610,1.0553,0.9829,
     *0.8988,
     *-6.399e-04,-6.399e-04,-6.399e-04,-6.399e-04,-6.399e-04,-6.399e-04,
     *-6.399e-04,-6.399e-04,-6.399e-04,-6.407e-04,-6.421e-04,-6.420e-04,
     *-6.183e-04,-6.207e-04,-6.384e-04,-6.399e-04,-5.795e-04,-6.334e-04,
     *-6.507e-04,
     *1.528e-05,1.528e-05,1.528e-05,1.528e-05,1.528e-05,1.528e-05,
     *1.528e-05,1.528e-05,1.528e-05,1.435e-05,1.516e-05,1.512e-05,
     *1.449e-05,1.457e-05,1.412e-05,1.314e-05,1.200e-05,1.020e-05,
     *8.069e-06/
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
