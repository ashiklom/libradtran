c=======================================================================
c SBDART gaseous absorption module
c
c BM adopted code from SBDART version 1.21
c=======================================================================

      SUBROUTINE SBTAUGAS(nz,z,p,t,wh,wo,dtausc,wl,sza,
     $     xn2, xo2, xco2, xch4, xn2o, xno2, 
     $     o4abs, n2abs, coabs, so2abs, nh3abs, noabs, hno3abs,
     $     dtau, wght, nk)

      IMPLICIT NONE

      INCLUDE 'sbtaugas.inc'

c  input 
      integer nz
      real z(mxly), p(mxly), t(mxly), wh(mxly), wo(mxly)
      real dtausc(mxly)
      real wl, sza
      real xn2, xo2, xco2, xch4, xn2o, xno2
      integer o4abs, n2abs, coabs, so2abs, nh3abs, noabs, hno3abs
c  output 
      real dtau(mxly,3), wght(3)
      integer nk

      integer iz, kd

      double precision uu(mxq,mxly)
      double precision amu0, wt,  dwl

      double precision dz(mxly), dp(mxly), dt(mxly), 
     &     dwh(mxly), dwo(mxly), ddtausc(mxly)

      double precision gwk(mxkd), 
     &     dtauk(mxly,mxkd), dtaug(mxly), dtaugc(mxly)

      integer kdist
      double precision xxn2,xxo2,xxco2,xxch4,xxn2o,
     &     xxco,xxno2,xxso2,xxnh3,xxno,xxhno3,xxo4

      real tausum

      data kdist/3/
      
      data    xxco, xxso2, xxnh3, xxno, xxhno3
     &     /  zip,    zip,   zip,  zip,   zip /
      
c switch on/off O4 absorption
      if (o4abs .NE. 0) then
         xxo4 = 1.0
      else
         xxo4 = 0.0
      endif

c switch on/off N2 absorption

      if (n2abs .NE. 0) then
         xxn2 = xn2
      else
         xxn2 = 0.0
      endif

c switch on/off absorption by CO, SO2, NH3, NO, HNO3

      if (coabs .EQ. 0) then
         xxco = 0.0
      endif

      if (so2abs .EQ. 0) then
         xxso2 = 0.0
      endif

      if (nh3abs .EQ. 0) then
         xxnh3 = 0.0
      endif

      if (noabs .EQ. 0) then
         xxno = 0.0
      endif

      if (hno3abs .EQ. 0) then
         xxhno3 = 0.0
      endif

      xxo2  = xo2 
      xxco2 = xco2
      xxch4 = xch4
      xxn2o = xn2o
      xxno2 = xno2

      do iz=1,nz
         dz (iz) = z (iz)
         dp (iz) = p (iz)
         dt (iz) = t (iz)
         dwh(iz) = wh(iz)
         dwo(iz) = wo(iz)
         ddtausc(iz) = dtausc(iz)
      enddo

      dwl = wl
      amu0=cos(sza*3.1415926/180.0)

      call modmix(xxn2,xxo2,xxco2,xxch4,xxn2o,
     &     xxco,xxno2,xxso2,xxnh3,xxno,xxhno3,xxo4)

      call absint(uu,nz,dz,dp,dt,dwh,dwo,0)
      
c  gas absorption

      call gasset(kdist,dwl,uu,amu0,nz,dz,nk,gwk,dtauk,dtaugc,0)

c   gasset output:
c   gwk     k-distribution weights
c   dtauk   k-distribution optical depth increments
c   dtaugc  continuum

c  loop through k-distribution terms

      do kd=1,nk 
         
         call depthscl(kdist, kd, nk, nz, dwl, ddtausc, 
     &        gwk, dtauk, dtaugc,
     &        wt, dtaug)

         tausum=0.0
         do iz=1,nz
            dtau(iz,kd) = dtaug(iz)
            tausum = tausum+dtau(iz,kd)
         enddo
         
         wght(kd) = wt
      enddo

      RETURN
      END
