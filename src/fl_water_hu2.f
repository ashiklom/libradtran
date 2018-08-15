! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MAKE CHANGE IN RAD TO USE HU's WATER CLOUD OPTICS
!      do 20 ib = mbn, mb
!       call ice ( ib )
! !! BEGIN CHANGE
!     call water ( ib ) ! REPLACE WITH FOLLOWING IF THEN ELSE
! !!!
!     if ( ib <=6 ) then
!       call water_hu ( ib ) ! CALL Yong Hu's WATER CLOUD OPTICS For SW banib=1:6
!     else
!       call water ( ib )
!     endif
! !! END CHANGE
!       call rain ( ib )
!       call graup ( ib )
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine water_hu(ib)
      INCLUDE 'fl_radparams.inc'
 
      integer i, ib, j
      real pre, plwc, pde, piwc, re, fl, bz, wz, gz, dz, tw, ww, www
      real gg, x1, x2, x3, x4
 
      real xre, bext, coalb, asy, trc, trcasy
 
      common /clouds/pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /thick/dz(nvx)
      common /wat/tw(nvx), ww(nvx), www(nvx,4)
 
      do i = 1, nv
         if (plwc(i) .lt. 1.0E-5) then
            tw(i) = 0.0
            ww(i) = 0.0
            www(i,1) = 0.0
            www(i,2) = 0.0
            www(i,3) = 0.0
            www(i,4) = 0.0
         else
 
            xre = pre(i)
            if (xre .le. 20.) then
               call cldopt_l20 (ib, xre, bext, coalb, asy)
            else
               call cldopt_g20 (ib, xre, bext, coalb, asy, trc, trcasy)
            endif
            tw(i) = dz(i)*plwc(i)*bext           ! TAU from extinction
            ww(i) = 1.0000 - coalb               !SSA
            gg = asy                             !ASY
 
            x1 = gg
            x2 = x1*x1
            x3 = x2*x1
            x4 = x3*x1
            www(i,1) = 3.0*x1
            www(i,2) = 5.0*x2
            www(i,3) = 7.0*x3
            www(i,4) = 9.0*x4
! !      print'(I5,f10.2,f10.4,f10.6,f10.2)',ib,pre(i),gg,ww(i),bext
         endif
 
      end do
      return 
      end 
!============================================================================
c   Fu-Liou band water cloud optical properties:
c   (if you have questions, please call Yong Hu at 757-864-9824
c                                 or email yhu@adm.larc.nasa.gov )
c----------------------------------------------------------------------------
c
c   The cloud optical properties parameterization code for Fu-Liou bands
cc   PARTICLE SIZE LESS THAN 20 MICRON
      subroutine cldopt_l20(iband, xre, bext, coalb, asy)
c
c   input variables:
c     iband:   the fu-liou band ( 1 - 6 )
c                iband=1: 0.2-0.7; iband=2: 0.7-1.3, ...
! !NO     xlwp:   liquid water path  ( gram per square meter )
c     xre:    effective droplet size
c
c
      real xlwp, xre
      integer iband
 
c   output variables:
c     bext:  extinction coefficient
c     coalb:  co-albedo ( = 1. - s.s.a )
c     asy:   asymmetry factor
c
      real bext, coalb, asy
c
      real coefbt(4), coefca(4), coefas(4)
      real coefbt1(4), coefca1(4), coefas1(4)
      real coefbt2(4), coefca2(4), coefas2(4)
      real coefbt3(4), coefca3(4), coefas3(4)
      real coefbt4(4), coefca4(4), coefas4(4)
      real coefbt5(4), coefca5(4), coefas5(4)
      real coefbt6(4), coefca6(4), coefas6(4)
 
      integer i
      real xxre
      data coefbt1/ -0.50313E-03, 0.15272E+01, 0.45841E+00,  - 
     1   0.29436E+00/ 
      data coefca1/-0.36882E-06,0.30867E-06,-0.22233E-07,0.71604E-09/
      data coefas1/0.80079E+00,0.11321E-01,-0.64926E-03,0.12850E-04/
      data coefbt2/ -0.75121E-03, 0.15492E+01, 0.57192E+00,  - 
     1   0.53427E-01/ 
      data coefca2/ 0.48969E-05, 0.30593E-04, 0.37730E-07,  - 
     1   0.49663E-08/ 
      data coefas2/0.70709E+00,0.32338E-01,-0.22309E-02,0.50809E-04/
      data coefbt3/-0.22040E-02,0.16009E+01,0.50647E+00,0.86071E+00/
      data coefca3/-0.70109E-03,0.14401E-02,-0.28917E-04,0.60054E-06/
      data coefas3/0.70192E+00,0.30346E-01,-0.19202E-02,0.41345E-04/
      data coefbt4/-0.11612E-01,0.19447E+01,-0.27681E+01,0.11911E+02/
      data coefca4/-0.91441E-02,0.75884E-02,-0.24023E-03,0.46077E-05/
      data coefas4/ 0.85008E+00, -0.11273E-01, 0.15253E-02,  - 
     1   0.43008E-04/ 
      data coefbt5/ 0.48334E-02, 0.13934E+01, 0.28501E+01,  - 
     1   0.35198E+01/ 
      data coefca5/0.18698E+00,0.32091E-01,-0.19028E-02,0.39806E-04/
      data coefas5/0.85272E+00,0.16497E-01,-0.97105E-03,0.19844E-04/
      data coefbt6/ 0.52556E-02, 0.14316E+01, 0.24743E+01, 0.16562E+01/ 
      data coefca6/-0.28624E-01,0.21779E-01,-0.79026E-03,0.13885E-04/
      data coefas6/ 0.90047E+00, -0.41211E-01, 0.42051E-02,  - 
     1   0.10977E-03/ 
 
      if (iband .eq. 1) then
         do i = 1, 4
            coefbt(i) = coefbt1(i)
            coefca(i) = coefca1(i)
            coefas(i) = coefas1(i)
         end do
      endif
      if (iband .eq. 2) then
         do i = 1, 4
            coefbt(i) = coefbt2(i)
            coefca(i) = coefca2(i)
            coefas(i) = coefas2(i)
         end do
      endif
      if (iband .eq. 3) then
         do i = 1, 4
            coefbt(i) = coefbt3(i)
            coefca(i) = coefca3(i)
            coefas(i) = coefas3(i)
         end do
      endif
      if (iband .eq. 4) then
         do i = 1, 4
            coefbt(i) = coefbt4(i)
            coefca(i) = coefca4(i)
            coefas(i) = coefas4(i)
         end do
      endif
      if (iband .eq. 5) then
         do i = 1, 4
            coefbt(i) = coefbt5(i)
            coefca(i) = coefca5(i)
            coefas(i) = coefas5(i)
         end do
      endif
      if (iband .eq. 6) then
         do i = 1, 4
            coefbt(i) = coefbt6(i)
            coefca(i) = coefca6(i)
            coefas(i) = coefas6(i)
         end do
      endif
 
      xxre = 1./xre
      bext = 0.
      coalb = 0.
      asy = 0.
      do i = 1, 4
         bext = bext + coefbt(i)*xxre**(i - 1.)
         coalb = coalb + coefca(i)*xre**(i - 1.)
         asy = asy + coefas(i)*xre**(i - 1.)
      end do
      bext = bext*1000.                      !! (m2g-1) --> (km-1 m3g-1)
!TAU  bext=bext*xlwp
      return 
      end 
c
c----------------------------------------------------------------------------
c
c   The cloud optical properties parameterization code for Fu-Liou bands
c   PARTICLE SIZE GREATER THAN 20 MICRON
      subroutine cldopt_g20(iband, xre, bext, coalb, asy, trc, trcasy)
c
c   input variables:
c     iband:   the fu-liou band ( 1 - 6 )
c                iband=1: 0.2-0.7; iband=2: 0.7-1.3, ...
!  !NO    xlwp:   liquid water path  ( gram per square meter )
c     xre:    effective droplet size
c
c
      real xlwp, xre
      integer iband
 
c   output variables:
c     bext:  extinction coefficient
c     coalb:  co-albedo ( = 1. - s.s.a )
c     asy:   asymmetry factor
c
      real bext, coalb, asy, trc, trcasy
c
      real coefbt(4), coefca(4), coefas(4), coeftr(4), coefta(4)
      real coefbt1(4), coefca1(4), coefas1(4), coeftr1(4), coefta1(4)
      real coefbt2(4), coefca2(4), coefas2(4), coeftr2(4), coefta2(4)
      real coefbt3(4), coefca3(4), coefas3(4), coeftr3(4), coefta3(4)
      real coefbt4(4), coefca4(4), coefas4(4), coeftr4(4), coefta4(4)
      real coefbt5(4), coefca5(4), coefas5(4), coeftr5(4), coefta5(4)
      real coefbt6(4), coefca6(4), coefas6(4), coeftr6(4), coefta6(4)
 
      integer i
      real xxre
      data coefbt1/ 0.23418E-01, -0.32565E+00, 0.48340E+02,  - 
     1   0.41043E+03/ 
      data coefca1/-0.21786E-04,0.26767E-05,-0.10008E-06,0.12921E-08/
      data coefas1/ 0.84794E+00, 0.14857E-02, -0.67237E-05,  - 
     1   0.40686E-06/ 
      data coeftr1/-0.52739E+00,0.11323E+00,-0.42522E-02,0.53457E-04/
      data coefta1/ 0.90311E+00, -0.18028E-01, 0.77449E-03,  - 
     1   0.10696E-04/ 
      data coefbt2/ 0.67806E-02, 0.96850E+00, 0.15519E+02,  - 
     1   0.12923E+03/ 
      data coefca2/ 0.44711E-02, -0.51398E-03, 0.22043E-04,  - 
     1   0.29747E-06/ 
      data coefas2/0.67621E+00,0.20439E-01,-0.71104E-03,0.83388E-05/
      data coeftr2/ 0.60176E+00, -0.23758E-01, 0.11284E-02,  - 
     1   0.16278E-04/ 
      data coefta2/ 0.90592E+00, -0.16972E-01, 0.67994E-03,  - 
     1   0.88442E-05/ 
      data coefbt3/ 0.18332E-02, 0.13641E+01, 0.57843E+01,  - 
     1   0.42550E+02/ 
      data coefca3/-0.21214E-01,0.37580E-02,-0.11065E-03,0.14348E-05/
      data coefas3/0.74915E+00,0.12089E-01,-0.38639E-03,0.43741E-05/
      data coeftr3/ 0.32956E+00, 0.67449E-02, -0.49333E-04,  - 
     1   0.44926E-06/ 
      data coefta3/ 0.82010E+00, -0.60738E-02, 0.24663E-03,  - 
     1   0.29954E-05/ 
      data coefbt4/ 0.10629E-01, 0.68353E+00, 0.23881E+02,  - 
     1   0.19519E+03/ 
      data coefca4/-0.25940E-02,0.55061E-02,-0.74017E-04,0.69346E-06/
      data coefas4/ 0.81411E+00, 0.46326E-02, -0.38183E-04,  - 
     1   0.59648E-06/ 
      data coeftr4/0.19206E+00,0.20226E-01,-0.44379E-03,0.42292E-05/
      data coefta4/ 0.78841E+00, -0.16535E-02, 0.15572E-03,  - 
     1   0.28824E-05/ 
      data coefbt5/ 0.71208E-02, 0.97718E+00, 0.16073E+02,  - 
     1   0.12230E+03/ 
      data coefca5/ 0.32924E+00, 0.31633E-02, -0.31895E-05,  - 
     1   0.62926E-06/ 
      data coefas5/0.80288E+00,0.15341E-01,-0.51183E-03,0.57613E-05/
      data coeftr5/0.21655E+00,0.40081E-01,-0.10761E-02,0.10858E-04/
      data coefta5/0.78023E+00,0.10683E-01,-0.45877E-03,0.55770E-05/
      data coefbt6/ 0.96706E-02, 0.76972E+00, 0.22557E+02,  - 
     1   0.17571E+03/ 
      data coefca6/ 0.10667E+00, 0.20174E-02, 0.21095E-03,  - 
     1   0.36733E-05/ 
      data coefas6/0.58782E+00,0.28048E-01,-0.89508E-03,0.10214E-04/
      data coeftr6/-0.22692E+00,0.55017E-01,-0.14222E-02,0.14334E-04/
      data coefta6/0.60193E+00,0.20226E-01,-0.74250E-03,0.91387E-05/
 
      if (iband .eq. 1) then
         do i = 1, 4
            coefbt(i) = coefbt1(i)
            coefca(i) = coefca1(i)
            coefas(i) = coefas1(i)
            coeftr(i) = coeftr1(i)
            coefta(i) = coefta1(i)
         end do
      endif
      if (iband .eq. 2) then
         do i = 1, 4
            coefbt(i) = coefbt2(i)
            coefca(i) = coefca2(i)
            coefas(i) = coefas2(i)
            coeftr(i) = coeftr2(i)
            coefta(i) = coefta2(i)
         end do
      endif
      if (iband .eq. 3) then
         do i = 1, 4
            coefbt(i) = coefbt3(i)
            coefca(i) = coefca3(i)
            coefas(i) = coefas3(i)
            coeftr(i) = coeftr3(i)
            coefta(i) = coefta3(i)
         end do
      endif
      if (iband .eq. 4) then
         do i = 1, 4
            coefbt(i) = coefbt4(i)
            coefca(i) = coefca4(i)
            coefas(i) = coefas4(i)
            coeftr(i) = coeftr4(i)
            coefta(i) = coefta4(i)
         end do
      endif
      if (iband .eq. 5) then
         do i = 1, 4
            coefbt(i) = coefbt5(i)
            coefca(i) = coefca5(i)
            coefas(i) = coefas5(i)
            coeftr(i) = coeftr5(i)
            coefta(i) = coefta5(i)
         end do
      endif
      if (iband .eq. 6) then
         do i = 1, 4
            coefbt(i) = coefbt6(i)
            coefca(i) = coefca6(i)
            coefas(i) = coefas6(i)
            coeftr(i) = coeftr6(i)
            coefta(i) = coefta6(i)
         end do
      endif
 
      xxre = 1./xre
      bext = 0.
      coalb = 0.
      asy = 0.
      trc = 0.
      trcasy = 0.
      do i = 1, 4
         bext = bext + coefbt(i)*xxre**(i - 1.)
         coalb = coalb + coefca(i)*xre**(i - 1.)
         asy = asy + coefas(i)*xre**(i - 1.)
         trc = trc + coeftr(i)*xre**(i - 1.)
         trcasy = trcasy + coefta(i)*xre**(i - 1.)
      end do
      bext = bext*1000.                      !! (m2g-1) --> (km-1 m3g-1)
! ! TAU bext=bext*xlwp
      return 
      end 
 
