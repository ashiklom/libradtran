cbm      subroutine thicks
cbm      INCLUDE 'fl_radparams.inc'
cbm      integer i
cbm      real pp, pt, pa, ph, po, dz
cbmc##      include 'rad_0698.h'
cbm      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
cbm      common /thick/dz(nvx)
cbm 
cbm      do i = 1, nv
cbm         dz(i) = 0.0146337*(pt(i)+pt(i+1))*alog(pp(i+1)/pp(i))
cbm      end do
cbm 
cbm      return 
cbm      end 
 
      subroutine ice(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, ibr
      real pre, plwc, pde, piwc, ap, bp, cps, dps, cpir, dz, ti
      real wi,wwi,fw1,fw2,fw3,fd,wf1,wf2,wf3,wf4,gg,x1,x2,x3,x4
 
      common /clouds/pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /ic1/ap(3,mbx), bp(4,mbx), cps(4,4,mbsx), dps(4,mbsx), cpir
     1   (4,mbirx)
      common /thick/dz(nvx)
      common /ic/ti(nvx), wi(nvx), wwi(nvx,4)
 
c ******************************************************************
c The constant 1000.0 below is to consider the units of dz(i) is km.
c ******************************************************************
      do i = 1, nv
         if (piwc(i) .lt. 1.0E-5) then
            ti(i) = 0.0
            wi(i) = 0.0
            wwi(i,1) = 0.0
            wwi(i,2) = 0.0
            wwi(i,3) = 0.0
            wwi(i,4) = 0.0
         else
            fw1 = pde(i)
            fw2 = fw1*pde(i)
            fw3 = fw2*pde(i)
            ti(i) = dz(i)*1000.0*piwc(i)*(ap(1,ib)+ap(2,ib)/fw1+ap(3,ib)
     1         /fw2)
            wi(i) = 1.0 - (bp(1,ib)+bp(2,ib)*fw1+bp(3,ib)*fw2+bp(4,ib)*
     1         fw3)
 
            if (ib .le. mbs) then
 
               fd=dps(1,ib)+dps(2,ib)*fw1+dps(3,ib)*fw2+dps(4,ib)*fw3
               wf1 = cps(1,1,ib) + cps(2,1,ib)*fw1 + cps(3,1,ib)*fw2 + 
     1            cps(4,1,ib)*fw3
               wwi(i,1) = (1.0 - fd)*wf1 + 3.0*fd
               wf2 = cps(1,2,ib) + cps(2,2,ib)*fw1 + cps(3,2,ib)*fw2 + 
     1            cps(4,2,ib)*fw3
               wwi(i,2) = (1.0 - fd)*wf2 + 5.0*fd
               wf3 = cps(1,3,ib) + cps(2,3,ib)*fw1 + cps(3,3,ib)*fw2 + 
     1            cps(4,3,ib)*fw3
               wwi(i,3) = (1.0 - fd)*wf3 + 7.0*fd
               wf4 = cps(1,4,ib) + cps(2,4,ib)*fw1 + cps(3,4,ib)*fw2 + 
     1            cps(4,4,ib)*fw3
               wwi(i,4) = (1.0 - fd)*wf4 + 9.0*fd
 
            else
 
               ibr = ib - mbs
               gg = cpir(1,ibr) + cpir(2,ibr)*fw1 + cpir(3,ibr)*fw2 + 
     1            cpir(4,ibr)*fw3
               x1 = gg
               x2 = x1*gg
               x3 = x2*gg
               x4 = x3*gg
               wwi(i,1) = 3.0*x1
               wwi(i,2) = 5.0*x2
               wwi(i,3) = 7.0*x3
               wwi(i,4) = 9.0*x4
 
            endif
         endif
      end do
 
      return 
      end 
 
      subroutine water(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, j
      real pre, plwc, pde, piwc, re, fl, bz, wz, gz, dz, tw, ww, www
      real gg, x1, x2, x3, x4
      common /clouds/pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common/wat1/re(ncx),fl(ncx),bz(ncx,mbx),wz(ncx,mbx),gz(ncx,mbx)
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
 
            if (pre(i) .lt. re(1)) then
c            ----------------------------------------------------
c            A cloud with the effective radius smaller than 4.18
c            um is assumed to have an effective radius of 4.18 um
c            with respect to the single scattering properties.
c            ----------------------------------------------------
               tw(i) = dz(i)*plwc(i)*bz(1,ib)/fl(1)
               ww(i) = wz(1,ib)
               x1 = gz(1,ib)
               x2 = x1*gz(1,ib)
               x3 = x2*gz(1,ib)
               x4 = x3*gz(1,ib)
               www(i,1) = 3.0*x1
               www(i,2) = 5.0*x2
               www(i,3) = 7.0*x3
               www(i,4) = 9.0*x4
 
            else if (pre(i) .gt. re(nc)) then
c            ----------------------------------------------------
c            A cloud with the effective radius larger than 31.23
c            um is assumed to have an effective radius of 31.18 um
c            with respect to the single scattering properties.
c            ----------------------------------------------------
               tw(i) = dz(i)*plwc(i)*bz(nc,ib)/fl(nc)
               ww(i) = wz(nc,ib)
               x1 = gz(nc,ib)
               x2 = x1*gz(nc,ib)
               x3 = x2*gz(nc,ib)
               x4 = x3*gz(nc,ib)
               www(i,1) = 3.0*x1
               www(i,2) = 5.0*x2
               www(i,3) = 7.0*x3
               www(i,4) = 9.0*x4
 
            else
 
               j = 1
 1002          continue
               if (pre(i).lt.re(j) .or. pre(i).gt.re(j+1)) then
                  j = j + 1
                  go to 1002
               endif
 
               tw(i) = dz(i)*plwc(i)*(bz(j,ib)/fl(j)+(bz(j+1,ib)/fl(j+1)
     1            -bz(j,ib)/fl(j))/(1.0/re(j+1)-1.0/re(j))*(1.0/pre(i)-
     2            1.0/re(j)))
               ww(i) = wz(j,ib) + (wz(j+1,ib)-wz(j,ib))/(re(j+1)-re(j))*
     1            (pre(i)-re(j))
               gg = gz(j,ib) + (gz(j+1,ib)-gz(j,ib))/(re(j+1)-re(j))*(
     1            pre(i)-re(j))
               x1 = gg
               x2 = x1*gg
               x3 = x2*gg
               x4 = x3*gg
               www(i,1) = 3.0*x1
               www(i,2) = 5.0*x2
               www(i,3) = 7.0*x3
               www(i,4) = 9.0*x4
            endif
         endif
 
      end do
      return 
      end 
 
      subroutine rayle2
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, trp
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /ray2/trp(nvx)
 
      do i = 1, nv
         trp(i) = 14.6337*(pp(i)+pp(i+1))*alog(pp(i+1)/pp(i))
      end do
      return 
      end 
 
      subroutine rayle(ib, u0, ig, lray)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      logical lray
      integer i, ib, ig
      real u0, ri, trp, tr, wr, wwr, x, rix, riy
 
      common /ray1/ri(mbsx), rix(10), riy(10)
      common /ray2/trp(nvx)
      common /ray/tr(nvx), wr(nvx), wwr(nvx,4)
 
      if (ib .le. mbs) then
         if (ib .eq. 1) then
            if (lray) then
               x = riy(ig)                       !LOG & Lin
            else
               x = (-3.902860E-6*u0*u0) + 6.120070E-6*u0 + 4.177440E-6
!          x= rix(ig) !LIN
            endif
         else
            x = ri(ib)
         endif
! !!  x=1.0E-09 !!! ZERO RAYLEIGH
         do i = 1, nv
            tr(i) = trp(i)*x
            wr(i) = 1.0
            wwr(i,1) = 0.0
            wwr(i,2) = 0.5
            wwr(i,3) = 0.0
            wwr(i,4) = 0.0
         end do
      else
         do i = 1, nv
            tr(i) = 0.0
            wr(i) = 0.0
            wwr(i,1) = 0.0
            wwr(i,2) = 0.0
            wwr(i,3) = 0.0
            wwr(i,4) = 0.0
         end do
      endif
 
      return 
      end 
 
      subroutine rain(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib
      real prwc, rwc, brn, wrnf, grn, dz, trn, wrn, wwrn
      real x1, x2, x3, x4, y1, y2, y3, y4
      common /rains/prwc(nvx)
      common /rai1/rwc, brn(mbx), wrnf(mbx), grn(mbx)
      common /thick/dz(nvx)
      common /rai/trn(nvx), wrn(nvx), wwrn(nvx,4)
 
      x1 = grn(ib)
      x2 = x1*grn(ib)
      x3 = x2*grn(ib)
      x4 = x3*grn(ib)
      y1 = 3.0*x1
      y2 = 5.0*x2
      y3 = 7.0*x3
      y4 = 9.0*x4
 
      do i = 1, nv
         if (prwc(i) .lt. 1.0E-5) then
            trn(i) = 0.0
            wrn(i) = 0.0
            wwrn(i,1) = 0.0
            wwrn(i,2) = 0.0
            wwrn(i,3) = 0.0
            wwrn(i,4) = 0.0
         else
            trn(i) = dz(i)*prwc(i)*brn(ib)/rwc
            wrn(i) = wrnf(ib)
            wwrn(i,1) = y1
            wwrn(i,2) = y2
            wwrn(i,3) = y3
            wwrn(i,4) = y4
         endif
      end do
      return 
      end 
 
      subroutine graup(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib
      real pgwc,gwc,bg,wgf,gg,dz,tgr,wgr,wwgr,x1,x2,x3,x4,y1,y2,y3,y4
      common /graups/pgwc(nvx)
      common /gra1/gwc, bg(mbx), wgf(mbx), gg(mbx)
      common /thick/dz(nvx)
      common /gra/tgr(nvx), wgr(nvx), wwgr(nvx,4)
 
      x1 = gg(ib)
      x2 = x1*gg(ib)
      x3 = x2*gg(ib)
      x4 = x3*gg(ib)
      y1 = 3.0*x1
      y2 = 5.0*x2
      y3 = 7.0*x3
      y4 = 9.0*x4
 
      do i = 1, nv
         if (pgwc(i) .lt. 1.0E-5) then
            tgr(i) = 0.0
            wgr(i) = 0.0
            wwgr(i,1) = 0.0
            wwgr(i,2) = 0.0
            wwgr(i,3) = 0.0
            wwgr(i,4) = 0.0
         else
            tgr(i) = dz(i)*pgwc(i)*bg(ib)/gwc
            wgr(i) = wgf(ib)
            wwgr(i,1) = y1
            wwgr(i,2) = y2
            wwgr(i,3) = y3
            wwgr(i,4) = y4
         endif
      end do
 
      return 
      end 
 
 
 
c---------- 4/1/97 (5) -- PREVIOUS 465 LINES
 
      subroutine ql_rh(rh, tl, pl, ql)
 
c      rh (0-100)
c      tl (K)
c      pl (mb)
c      q (g/g)
 
      es = satvap(tl)
      ws = 0.622*es/(pl - es)
      rh = ql/ws*100.
      return 
      end 
 
      real function satvap (temp2)
      temp = temp2 - 273.155
      if (temp .lt. (-20.)) then                 !!!! ice saturation
         toot = 273.16/temp2
         toto = 1/toot
         eilog = (-9.09718*(toot - 1)) - 3.56654*(log(toot)/log(10.)) + 
     1      .876793*(1 - toto) + log(6.1071)/log(10.)
         satvap = 10**eilog
      else
         tsot = 373.16/temp2
         ewlog = (-7.90298*(tsot - 1)) + 5.02808*(log(tsot)/log(10.))
         ewlog2 = ewlog - 1.3816E-07*(10**(11.344*(1 - 1/tsot)) - 1)
         ewlog3 = ewlog2 + .0081328*(10**(-3.49149*(tsot - 1)) - 1)
         ewlog4 = ewlog3 + log(1013.246)/log(10.)
         satvap = 10**ewlog4
      endif
 
      return 
      end 
 
      subroutine gascon(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib
      real tgm, vv(18)
 
      common /con/tgm(nvx)
      data vv/10*0.0,1175.0,1040.0,890.0,735.0,605.0,470.0,340.0,0.0/
 
      if (ib.gt.10 .and. ib.lt.18) then
         call qopcon (vv(ib), tgm)
      else
         do i = 1, nv
            tgm(i) = 0.0
         end do
      endif
 
      return 
      end 
 
      subroutine gases(ib, ig, hk)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, ig
      real hk, pp, pt, pa, ph, po, umco2, umch4, umn2o, hk1, fk1o3
      real hk2, c2h2o, hk3, c3h2o, hk4, c4h2o, hk5, c5h2o, hk6, c6h2o
      real hk7,c7h2o,hk8,c8h2o,hk9,c9h2o,hk10,c10h2o,c10ch4,c10n2o
      real hk11,c11h2o,c11ch4,c11n2o,hk12,c12o3,c12h2o,hk13,c13h2o
      real hk14, c14hca, c14hcb, hk15, c15hca, c15hcb, hk16, c16h2o
      real hk17, c17h2o, hk18, c18h2o, tg, fk
      real fkg(nv1x), fkga(nv1x), fkgb(nv1x), pq(nv1x)
      real tg1(nvx), tg2(nvx), tg3(nvx)
      real sol_spect, fk1h2o
 
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
 
      common /band1/hk1(10), fk1o3(10), sol_spect(0:7), fk1h2o(10)
      common /band2/hk2(8), c2h2o(3,11,8)
      common /band3/hk3(12), c3h2o(3,11,12)
      common /band4/hk4(7), c4h2o(3,11,7)
      common /band5/hk5(12), c5h2o(3,11,12)
      common /band6/hk6(5), c6h2o(3,11,5)
      common /band7/hk7(2), c7h2o(3,19,2)
      common /band8/hk8(3), c8h2o(3,19,3)
      common /band9/hk9(4), c9h2o(3,19,4)
      common /band10/hk10(4), c10h2o(3,19,4), c10ch4(3,19), c10n2o(3,19)
      common /band11/hk11(3), c11h2o(3,19,3), c11ch4(3,19), c11n2o(3,19)
      common /band12/hk12(5), c12o3(3,19,5), c12h2o(3,19)
      common /band13/hk13(2), c13h2o(3,19,2)
      common /band14/hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
      common /band15/hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
      common /band16/hk16(7), c16h2o(3,19,7)
      common /band17/hk17(7), c17h2o(3,19,7)
      common /band18/hk18(8), c18h2o(3,19,8)
      common /gas/tg(nvx)
      integer idkfr
      common /dkfrwn/idkfr
      logical lchou, lband6a, lpar, lray
      common /chou/lchou, lband6a, lpar, lray
 
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) ib
      stop 
c  ---------------------------------------------------------------------
c  In this band ( 50000 - 14500 cm**-1 ), we have considered the nongray
c  gaseous absorption of O3.    619.618 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    1 continue
      fk = fk1o3(ig)
! !   fk=1.0E-09 !ZERO OZONE
      call qopo3s (fk, tg)
 
      if (lpar) call qoph2o_chou (fk1h2o(ig), tg)
!      hk = 619.618 * hk1(ig)
      hk = sol_spect(1)*hk1(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 14500 - 7700 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.  484.295 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    2 continue
      call qks (c2h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
!      hk = 484.295 * hk2(ig)
      hk = sol_spect(2)*hk2(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 7700 - 5250 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O. 149.845 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    3 continue
      call qks (c3h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
!      hk = 149.845 * hk3(ig)
      hk = sol_spect(3)*hk3(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 5250 - 4000 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O. 48.7302 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    4 continue
      call qks (c4h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
!      hk = 48.7302 * hk4(ig)
      hk = sol_spect(4)*hk4(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 4000 - 2850 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O. 31.6576 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    5 continue
      call qks (c5h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      if (lband6a) then
         hk = (sol_spect(5)+sol_spect(7))*hk5(ig)
      else
!      hk = 31.6576 * hk5(ig)
         hk = sol_spect(5)*hk5(ig)
      endif
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 2850 - 2500 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O. 5.79927 is the solar energy contained in
c  the band in units of Wm**-2.
c  ---------------------------------------------------------------------
    6 continue
      call qks (c6h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
!      hk = 5.79927 * hk6(ig)
      hk = sol_spect(6)*hk6(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 2200 - 1900 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
    7 continue
      call qki (c7h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk7(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 1900 - 1700 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
    8 continue
      call qki (c8h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk8(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 1700 - 1400 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
    9 continue
      call qki (c9h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk9(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 1400 - 1250 cm**-1 ), we have considered the
c  overlapping absorption of H2O, CH4, and N2O by approach one of
c  Fu(1991).
c  ---------------------------------------------------------------------
   10 continue
      call qki (c10h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg1)
      call qki (c10ch4, fkg)
      call qopch4 (fkg, tg2)
      call qki (c10n2o, fkg)
      call qopn2o (fkg, tg3)
      do i = 1, nv
         tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
      end do
      hk = hk10(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 1250 - 1100 cm**-1 ), we have considered the
c  overlapping absorption of H2O, CH4, and N2O by approach one of
c  Fu(1991).
c  ---------------------------------------------------------------------
   11 continue
      if (idkfr .eq. 0) then
 
         call qki (c11h2o(1,1,ig), fkg)
         call qoph2o (fkg, tg1)
 
         call qki (c11ch4, fkg)
         call qopch4 (fkg, tg2)
         call qki (c11n2o, fkg)
         call qopn2o (fkg, tg3)
         do i = 1, nv
            tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
         end do
 
         hk = hk11(ig)
 
      else if (idkfr .eq. 1) then
         call ck_dkfr_win_all (ib, ig, hk, tg)
      else if (idkfr .eq. 2) then
         call ck_dkfr_win_hyb (ib, ig, hk, tg)
      endif
 
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 1100 - 980 cm**-1 ), we have considered the overlapping
c  absorption of H2O and O3 by approach one of Fu(1991).
c  ---------------------------------------------------------------------
   12 continue
      if (idkfr .eq. 0) then
         call qkio3 (c12o3(1,1,ig), fkg)
         call qopo3i (fkg, tg1)
         call qki (c12h2o, fkg)
         call qoph2o (fkg, tg2)
 
         do i = 1, nv
            tg(i) = tg1(i) + tg2(i)
         end do
 
         hk = hk12(ig)
 
      else if (idkfr .eq. 1) then
         call ck_dkfr_win_all (ib, ig, hk, tg)
      else if (idkfr .eq. 2) then
         call ck_dkfr_win_hyb (ib, ig, hk, tg)
      endif
 
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 980 - 800 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
   13 continue
      if (idkfr .eq. 0) then
         call qki (c13h2o(1,1,ig), fkg)
         call qoph2o (fkg, tg)
         hk = hk13(ig)
      else if (idkfr .eq. 1) then
         call ck_dkfr_win_all (ib, ig, hk, tg)
      else if (idkfr .eq. 2) then
         call ck_dkfr_win_hyb (ib, ig, hk, tg)
      endif
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 800 - 670 cm**-1), we have considered the overlapping
c  absorption of H2O and CO2 by approach two of Fu(1991).
c  ---------------------------------------------------------------------
   14 continue
      do i = 1, nv1
         if (pp(i) .ge. 63.1) then
            pq(i) = ph(i)
         else
            pq(i) = 0.0
         endif
      end do
      call qki (c14hca(1,1,ig), fkga)
      call qki (c14hcb(1,1,ig), fkgb)
      do i = 1, nv1
         fkg(i) = fkga(i)/330.0*umco2 + pq(i)*fkgb(i)
      end do
      call qophc (fkg, tg)
      hk = hk14(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 670 - 540 cm**-1), we have considered the overlapping
c  absorption of H2O and CO2 by approach two of Fu(1991).
c  ---------------------------------------------------------------------
   15 continue
      do i = 1, nv1
         if (pp(i) .ge. 63.1) then
            pq(i) = ph(i)
         else
            pq(i) = 0.0
         endif
      end do
      call qki (c15hca(1,1,ig), fkga)
      call qki (c15hcb(1,1,ig), fkgb)
      do i = 1, nv1
         fkg(i) = fkga(i)/330.0*umco2 + pq(i)*fkgb(i)
      end do
      call qophc (fkg, tg)
      hk = hk15(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 540 - 400 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
   16 continue
      call qki (c16h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk16(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 400 - 280 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
   17 continue
      call qki (c17h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk17(ig)
      go to 20
 
c  ---------------------------------------------------------------------
c  In this band ( 280 - 000 cm**-1 ), we have considered the nongray
c  gaseous absorption of H2O.
c  ---------------------------------------------------------------------
   18 continue
      call qki (c18h2o(1,1,ig), fkg)
      call qoph2o (fkg, tg)
      hk = hk18(ig)
 
   20 continue
      return 
      end 

c 338kfix
      subroutine qks(coefks, fkg)
      INCLUDE 'fl_radparams.inc'
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      real coefks(3,11)
      real fkg(nv1x)  
c  in Fu code here dimension fkg(nv1x)
      real stanp(11)
      real tk(nv1x)
      data stanp/ 10.0, 15.8, 25.1, 39.8, 63.1, 100.0, 158.0, 251.0, 
     1   398.0, 631.0, 1000.0/ 
      do i = 1, nv1
         if (pt(i) .le. 180.0) then
            tk(i) = 180.0
         else if (pt(i) .ge. 320.0) then
            tk(i) = 320.0
         else
            tk(i) = pt(i)
         endif
      end do
      i1 = 1
      do i = 1, nv1
         if (pp(i) .lt. stanp(1)) then
            x1 = exp(coefks(1,1)+coefks(2,1)*(tk(i)-245.0)+coefks(3,1)*(
     1         tk(i)-245.0)**2)
            fkg(i) = x1*pp(i)/stanp(1)
         else if (pp(i) .ge. stanp(11)) then
            y1 = (tk(i)-245.0)*(tk(i)-245.0)
            x1 = exp(coefks(1,10)+coefks(2,10)*(tk(i)-245.0)+coefks(3,10
     1         )*y1)
            x2 = exp(coefks(1,11)+coefks(2,11)*(tk(i)-245.0)+coefks(3,11
     1         )*y1)
            fkg(i)=x1+(x2-x1)/(stanp(11)-stanp(10))*(pp(i)-stanp(10))
         else
   30       continue
            if (pp(i) .lt. stanp(i1)) then
               y1 = (tk(i)-245.0)*(tk(i)-245.0)
               x1 = exp(coefks(1,i1-1)+coefks(2,i1-1)*(tk(i)-245.0)+
     1            coefks(3,i1-1)*y1)
               x2 = exp(coefks(1,i1)+coefks(2,i1)*(tk(i)-245.0)+coefks(3
     1            ,i1)*y1)
               fkg(i) = x1 + (x2 - x1)/(stanp(i1)-stanp(i1-1))*(pp(i)-
     1            stanp(i1-1))
               go to 1002
            endif
            i1 = i1 + 1
            go to 30
         endif
 1002    continue
      end do
      return 
      end 
 
      subroutine qki(coefki, fkg)
      INCLUDE 'fl_radparams.inc'
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      real coefki(3,19)
      real fkg(nv1x)
      real stanp(19)
      real tk(nv1x)
      data stanp/ 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 3.98, 6.31, 
     1   10.0, 15.8, 25.1, 39.8, 63.1, 100.0, 158.0, 251.0, 398.0, 631.0
     2   , 1000.0/ 
      do i = 1, nv1
         if (pt(i) .le. 180.0) then
            tk(i) = 180.0
         else if (pt(i) .ge. 320.0) then
            tk(i) = 320.0
         else
            tk(i) = pt(i)
         endif
      end do
      i1 = 1
      do i = 1, nv1
         if (pp(i) .lt. stanp(1)) then
            x1 = exp(coefki(1,1)+coefki(2,1)*(tk(i)-245.0)+coefki(3,1)*(
     1         tk(i)-245.0)**2)
            fkg(i) = x1*pp(i)/stanp(1)
         else if (pp(i) .ge. stanp(19)) then
            y1 = (tk(i)-245.0)*(tk(i)-245.0)
            x1 = exp(coefki(1,18)+coefki(2,18)*(tk(i)-245.0)+coefki(3,18
     1         )*y1)
            x2 = exp(coefki(1,19)+coefki(2,19)*(tk(i)-245.0)+coefki(3,19
     1         )*y1)
            fkg(i)=x1+(x2-x1)/(stanp(19)-stanp(18))*(pp(i)-stanp(18))
         else
   30       continue
            if (pp(i) .lt. stanp(i1)) then
               y1 = (tk(i)-245.0)*(tk(i)-245.0)
               x1 = exp(coefki(1,i1-1)+coefki(2,i1-1)*(tk(i)-245.0)+
     1            coefki(3,i1-1)*y1)
               x2 = exp(coefki(1,i1)+coefki(2,i1)*(tk(i)-245.0)+coefki(3
     1            ,i1)*y1)
               fkg(i) = x1 + (x2 - x1)/(stanp(i1)-stanp(i1-1))*(pp(i)-
     1            stanp(i1-1))
               go to 1002
            endif
            i1 = i1 + 1
            go to 30
         endif
 1002    continue
      end do
      return 
      end 
 
      subroutine qkio3(coefki, fkg)
      INCLUDE 'fl_radparams.inc'
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      real coefki(3,19)
      real fkg(nv1x)
      real stanp(19)
      real tk(nv1x)
      data stanp/ 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 3.98, 6.31, 
     1   10.0, 15.8, 25.1, 39.8, 63.1, 100.0, 158.0, 251.0, 398.0, 631.0
     2   , 1000.0/ 
      do i = 1, nv1
         if (pt(i) .le. 180.0) then
            tk(i) = 180.0
         else if (pt(i) .ge. 320.0) then
            tk(i) = 320.0
         else
            tk(i) = pt(i)
         endif
      end do
      i1 = 1
      do i = 1, nv1
         if (pp(i) .lt. stanp(1)) then
            x1 = exp(coefki(1,1)+coefki(2,1)*(tk(i)-250.0)+coefki(3,1)*(
     1         tk(i)-250.0)**2)
            fkg(i) = x1*pp(i)/stanp(1)
         else if (pp(i) .ge. stanp(19)) then
            y1 = (tk(i)-250.0)*(tk(i)-250.0)
            x1 = exp(coefki(1,18)+coefki(2,18)*(tk(i)-250.0)+coefki(3,18
     1         )*y1)
            x2 = exp(coefki(1,19)+coefki(2,19)*(tk(i)-250.0)+coefki(3,19
     1         )*y1)
            fkg(i)=x1+(x2-x1)/(stanp(19)-stanp(18))*(pp(i)-stanp(18))
         else
   30       continue
            if (pp(i) .lt. stanp(i1)) then
               y1 = (tk(i)-250.0)*(tk(i)-250.0)
               x1 = exp(coefki(1,i1-1)+coefki(2,i1-1)*(tk(i)-250.0)+
     1            coefki(3,i1-1)*y1)
               x2 = exp(coefki(1,i1)+coefki(2,i1)*(tk(i)-250.0)+coefki(3
     1            ,i1)*y1)
               fkg(i) = x1 + (x2 - x1)/(stanp(i1)-stanp(i1-1))*(pp(i)-
     1            stanp(i1-1))
               go to 1002
            endif
            i1 = i1 + 1
            go to 30
         endif
 1002    continue
      end do
      return 
      end 
 
c 338kfix
 
      subroutine qopo3s(fk, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fk, tg(nvx), fq

      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
      fq = 238.08*fk
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (po(i)+po(i+1))*(pp(i+1)-pp(i))*fq
         tg(i) = (po(i)+po(i+1))*( pa(i)*0.098 )*fq
      end do
 
c      do 20 i = 1, nv
c         tg(i) = tg(i) * 476.16 * fk
c20      continue
c 476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for O3.
 
      return 
      end 
 
      subroutine qoph2o(fkg, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fkg(nv1x),tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (fkg(i)*ph(i)+fkg(i+1)*ph(i+1))*(pp(i+1)-pp(i))*
cbm  1      634.9205
         tg(i) = (fkg(i)*ph(i)+fkg(i+1)*ph(i+1))*( pa(i)*0.098 )*
     1      634.9205
      end do

c      do 20 i = 1, nv
c         tg(i) = tg(i) * 1269.841
c20      continue
c 1269.841 = 2.24e4 / M * 10.0 / 9.8, where M = 18 for H2O.
 
      return 
      end 
 
 
      subroutine qoph2o_chou(fk, tg)
      INCLUDE 'fl_radparams.inc'
      integer i
      real pp, pt, pa, ph, po, fk, tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
!        do  i = 1, nv
!       tg(i) = tg(i) +  fk* 1.02
!     & * (ph(i)+ph(i+1)) *0.5
!     & * ( pp(i+1)-pp(i) )
!     & * (.5*(pp(i)+pp(i+1))/300.)**.8
!     & * (1.+0.00135*(pt(i)-240.)) +1.e-11
!     enddo
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = tg(i) + fk*0.51*(ph(i)+ph(i+1))*(pp(i+1)-pp(i))
         tg(i) = tg(i) + fk*0.51*(ph(i)+ph(i+1))*( pa(i)*0.098 )
      end do
 
      return 
      end 
 
      subroutine qopch4(fkg, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fkg(nv1x),tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (fkg(i)+fkg(i+1))*(pp(i+1)-pp(i))*6.3119E-4
         tg(i) = (fkg(i)+fkg(i+1))*( pa(i)*0.098 )*6.3119E-4
      end do
 
c      do 20 i = 1, nv
c         tg(i) = tg(i) * 1.26238e-3
c20      continue
c 1.26238e-3 = 2.24e4 / M * 10.0 / 9.8 * 1.6e-6 * M / 28.97, where
c M = 16 for CH4.
 
      return 
      end 
 
      subroutine qopn2o(fkg, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fkg(nv1x),tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (fkg(i)+fkg(i+1))*(pp(i+1)-pp(i))*1.10459E-4
         tg(i) = (fkg(i)+fkg(i+1))*( pa(i)*0.098 )*1.10459E-4
      end do
 
c      do 20 i = 1, nv
c         tg(i) = tg(i) * 2.20918e-4
c20      continue
c 2.20918e-4 = 2.24e4 / M * 10.0 / 9.8 * 0.28e-6 * M / 28.97, where
c M = 44 for N2O.
 
      return 
      end 
 
      subroutine qopo3i(fkg, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fkg(nv1x),tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (fkg(i)*po(i)+fkg(i+1)*po(i+1))*(pp(i+1)-pp(i))*238.08
         tg(i) = (fkg(i)*po(i)+fkg(i+1)*po(i+1))*( pa(i)*0.098 )*238.08
      end do
 
c      do 20 i = 1, nv
c         tg(i) = tg(i) * 476.16
c20      continue
 
      return 
      end 
 
      subroutine qophc(fkg, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp, pt, pa, ph, po, fkg(nv1x), tg(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (fkg(i)+fkg(i+1))*(pp(i+1)-pp(i))*0.5
         tg(i) = (fkg(i)+fkg(i+1))*( pa(i)*0.098 )*0.5
      end do
c  ------------------------
c  See page 86 of Fu (1991).
c  ------------------------
 
      return 
      end 
 
      subroutine qopcon(vv, tg)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real pp,pt,ph,po,ff(nv1x),pe(nv1x),tg(nvx),x,y,z,r,s,vv,w
 
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
      x = 4.18
      y = 5577.8
      z = 0.00787
      r = 0.002
      s = (x + y*exp((-z*vv)))/1013.25
 
      do i = 1, nv1
         pe(i) = pp(i)*ph(i)/(0.622 + 0.378*ph(i))
         w = exp(1800.0/pt(i)-6.08108)
         ff(i) = s*(pe(i)+r*pp(i))*w
      end do
 
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
      do i = 1, nv
cbm      tg(i) = (ff(i)*ph(i)+ff(i+1)*ph(i+1))*(pp(i+1)-pp(i))*0.5098835
         tg(i) = (ff(i)*ph(i)+ff(i+1)*ph(i+1))*( pa(i)*0.098 )*0.5098835
      end do
 
c      do 7 i = 1, nv
c       tg(i) = tg(i) * 10.0 / 9.80616
c7      continue
 
      return 
      end 
 
      subroutine comscp
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
 
      integer i, j, iac
      real ti, wi, wwi, tw, ww, www, trn, wrn, wwrn, tgr, wgr, wwgr
      real tr, wr, wwr, tgm, tg, wc1, wc2, wc3, wc4, wc, tt, tc(nvx)
      real tis, tws, trns, tgrs, fw
 
c---------- 10/29/96 (4)
      real tae, wae, wwae, taes(0:4)
      common /aer/tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)
c---------- 10/29/96 (4)
 
      common /ic/ti(nvx), wi(nvx), wwi(nvx,4)
      common /wat/tw(nvx), ww(nvx), www(nvx,4)
      common /rai/trn(nvx), wrn(nvx), wwrn(nvx,4)
      common /gra/tgr(nvx), wgr(nvx), wwgr(nvx,4)
      common /ray/tr(nvx), wr(nvx), wwr(nvx,4)
      common /con/tgm(nvx)
      common /gas/tg(nvx)
      common/dfsin/wc1(nvx),wc2(nvx),wc3(nvx),wc4(nvx),wc(nvx),tt(nvx)
      integer j1
C *** Do construct levels ***
 
      do i = 1, nv
 
         tc(i)=ti(i)+tw(i)+tr(i)+tgm(i)+tg(i)+trn(i)+tgr(i)
 
         do iac = 1, nac
            tc(i) = tc(i) + tae(i,iac)
         end do
 
c---------- 10/29/96 (5)
 
         tis = ti(i)*wi(i)
         tws = tw(i)*ww(i)
         trns = trn(i)*wrn(i)
         tgrs = tgr(i)*wgr(i)
 
c---------- 10/29/96 (6)
         do j1 = 1, 5
            taes(j1-1) = 0.0
         end do
         do iac = 1, nac
            taes(0) = taes(0) + tae(i,iac)*wae(i,iac)
            do j = 1, 4
               taes(j) = taes(j) + tae(i,iac)*wae(i,iac)*wwae(i,j,iac)
            end do
         end do
 
         fw = tis + tws + tr(i) + trns + tgrs + taes(0)
c---------- 10/29/96 (6)
 
         wc(i) = fw/tc(i)
 
         if (fw .lt. 1.0E-20) then
            wc1(i) = 0.0
            wc2(i) = 0.0
            wc3(i) = 0.0
            wc4(i) = 0.0
         else
c---------- 10/29/96 (7)
 
!         wc1(i) = ( tis * wwi(i,1) + tws * www(i,1) +
!     1      tr(i) * wwr(i,1) + trns * wwrn(i,1) + tgrs * wwgr(i,1) +
!     1      taes * wwae(i,1) )/fw
 
!        wc2(i) = ( tis * wwi(i,2) + tws * www(i,2) +
!     1      tr(i) * wwr(i,2) + trns * wwrn(i,2) + tgrs * wwgr(i,2) +
!     1       taes * wwae(i,2) )/fw
 
!         wc3(i) = ( tis * wwi(i,3) + tws * www(i,3) +
!     1      tr(i) * wwr(i,3) + trns * wwrn(i,3) + tgrs * wwgr(i,3) +
!     1       taes * wwae(i,3) )/fw
 
!         wc4(i) = ( tis * wwi(i,4) + tws * www(i,4) +
!     1      tr(i) * wwr(i,4) + trns * wwrn(i,4) + tgrs * wwgr(i,4) +
!     1       taes * wwae(i,4) )/fw
!----------------------------------------------------------------------
            wc1(i) = (tis*wwi(i,1)+tws*www(i,1)+tr(i)*wwr(i,1)+trns*wwrn
     1         (i,1)+tgrs*wwgr(i,1)+taes(1))/fw
            wc2(i) = (tis*wwi(i,2)+tws*www(i,2)+tr(i)*wwr(i,2)+trns*wwrn
     1         (i,2)+tgrs*wwgr(i,2)+taes(2))/fw
            wc3(i) = (tis*wwi(i,3)+tws*www(i,3)+tr(i)*wwr(i,3)+trns*wwrn
     1         (i,3)+tgrs*wwgr(i,3)+taes(3))/fw
            wc4(i) = (tis*wwi(i,4)+tws*www(i,4)+tr(i)*wwr(i,4)+trns*wwrn
     1         (i,4)+tgrs*wwgr(i,4)+taes(4))/fw
         endif
      end do
 
 
 
      tt(1) = tc(1)
      do i = 2, nv
         tt(i) = tt(i-1) + tc(i)
      end do
      return 
      end 
 
      subroutine planck(ib, pts)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, nv11, ibr, j, nd(mbirx)
      real pts, pp, pt, pa, ph, po, bf, bs, ve(mbirx), bt(nv1x)
      real bts, v1, v2, w, fq1, fq2, x
 
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /planci/bf(nv1x), bs
      data ve/ 2200.0, 1900.0, 1700.0, 1400.0, 1250.0, 1100.0, 980.0, 
     1   800.0, 670.0, 540.0, 400.0, 280.001/ 
      data nd/ 30, 20, 30, 15, 15, 12, 18, 13, 13, 14, 12, 28/ 
 
      nv11 = nv1 + 1
      ibr = ib - mbs
      bts = 0.0
 
      do i = 1, nv1
         bt(i) = 0.0
      end do
      v1 = ve(ibr)
      do j = 1, nd(ibr)
         v2 = v1 - 10.0
         w = (v1 + v2)*0.5
         fq1 = 1.19107E-8*w*w*w
         fq2 = 1.43884*w
         do i = 1, nv11
            if (i .eq. nv11) then
               x = fq1/(exp(fq2/pts) - 1.0)
               bts = bts + x
            else
               x = fq1/(exp(fq2/pt(i))-1.0)
               bt(i) = bt(i) + x
            endif
         end do
         v1 = v2
      end do
      do i = 1, nv1
         bf(i) = bt(i)*10.0
      end do
      bs = bts*10.0
 
      return 
      end 
 
      subroutine coeff1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, j, ib
      real a,u,p0d,p1d,p2d,p3d,p11d,p22d,p33d,x,w0w,w1w,w2w,w3w,fw
      real w, w1, w2, w3, t0, t1, u0, f0, b, c(4,5), q1, q2, q3, fq
      common /dis/a(4)
      common /point/u(4)
      common /legen/p0d(4), p1d(4), p2d(4), p3d(4)
      common /legen1/p11d(4,4), p22d(4,4), p33d(4,4)
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/b(4,3)
 
      x = 0.5*w
      w0w = x
      w1w = x*w1
      w2w = x*w2
      w3w = x*w3
      if (ib .le. mbs) then
         fw = u0*u0
         q1 = -w1w*u0
         q2 = w2w*(1.5*fw - 0.5)
         q3 = -w3w*(2.5*fw - 1.5)*u0
      endif
      fq = 0.5*w0w
 
      do i = 3, 4
         do j = 1, 4
            c(i,j) = fq + w1w*p11d(i,j) + w2w*p22d(i,j) + w3w*p33d(i,j)
            if (i .eq. j) then
               c(i,j) = (c(i,j)-1.0)/u(i)
            else
               c(i,j) = c(i,j)/u(i)
            endif
         end do
      end do
 
      do i = 1, 4
         if (ib .le. mbs) then
            c(i,5) = w0w + q1*p1d(i) + q2*p2d(i) + q3*p3d(i)
         else
            c(i,5) = 1.0
         endif
         c(i,5) = c(i,5)/u(i)
      end do
 
      b(1,1) = c(4,4) - c(4,1)
      b(1,2) = c(4,4) + c(4,1)
      b(2,1) = c(4,3) - c(4,2)
      b(2,2) = c(4,3) + c(4,2)
      b(3,1) = c(3,4) - c(3,1)
      b(3,2) = c(3,4) + c(3,1)
      b(4,1) = c(3,3) - c(3,2)
      b(4,2) = c(3,3) + c(3,2)
      b(1,3) = c(4,5) - c(1,5)
      b(2,3) = c(3,5) - c(2,5)
      b(3,3) = c(3,5) + c(2,5)
      b(4,3) = c(4,5) + c(1,5)
      return 
      end 
 
      subroutine coeff2
c  *****************************************************************
c  coefficient calculations for second order differential equations.
c  *****************************************************************
      integer ib
      real w, w1, w2, w3, t0, t1, u0, f0, b, a, d, fw1, fw2, fw3, fw4
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/b(4,3)
      common /coedf2/a(2,2,2), d(4)
 
      fw1 = b(1,1)*b(1,2)
      fw2 = b(2,1)*b(3,2)
      fw3 = b(3,1)*b(2,2)
      fw4 = b(4,1)*b(4,2)
      a(2,2,1) = fw1 + fw2
      a(2,1,1) = b(1,1)*b(2,2) + b(2,1)*b(4,2)
      a(1,2,1) = b(3,1)*b(1,2) + b(4,1)*b(3,2)
      a(1,1,1) = fw3 + fw4
      a(2,2,2) = fw1 + fw3
      a(2,1,2) = b(1,2)*b(2,1) + b(2,2)*b(4,1)
      a(1,2,2) = b(3,2)*b(1,1) + b(4,2)*b(3,1)
      a(1,1,2) = fw2 + fw4
      d(1) = b(3,2)*b(4,3) + b(4,2)*b(3,3) + b(2,3)/u0
      d(2) = b(1,2)*b(4,3) + b(2,2)*b(3,3) + b(1,3)/u0
      d(3) = b(3,1)*b(1,3) + b(4,1)*b(2,3) + b(3,3)/u0
      d(4) = b(1,1)*b(1,3) + b(2,1)*b(2,3) + b(4,3)/u0
 
      return 
      end 
 
      subroutine coeff4
c  *****************************************************************
c  coefficient calculations for fourth-order differential equations.
c  *****************************************************************
      integer ib
      real w, w1, w2, w3, t0, t1, u0, f0, a, d, b1, c1, z, x
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf2/a(2,2,2), d(4)
      common /coedf4/b1, c1, z(4)
      x = u0*u0
      b1 = a(2,2,1) + a(1,1,1)
      c1 = a(2,1,1)*a(1,2,1) - a(1,1,1)*a(2,2,1)
      z(1) = a(2,1,1)*d(3) + d(4)/x - a(1,1,1)*d(4)
      z(2) = a(1,2,1)*d(4) - a(2,2,1)*d(3) + d(3)/x
      z(3) = a(2,1,2)*d(1) + d(2)/x - a(1,1,2)*d(2)
      z(4) = a(1,2,2)*d(2) - a(2,2,2)*d(1) + d(1)/x
      return 
      end 
 
      subroutine coeffl
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ib, i
      real w,w1,w2,w3,t0,t1,u0,f0,a,d,b,b1,c1,z,x,aa,zz,a1,z1,fk1,fk2
      real dt, fw, a2, b2, fw1, fw2, y, zx, fq0, fq1
 
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedf1/b(4,3)
      common /coedf2/a(2,2,2), d(4)
      common /coedf4/b1, c1, z(4)
      common /coedfl/aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
 
      dt = t1 - t0
      x = sqrt(b1*b1 + 4.0*c1)
      fk1 = sqrt((b1 + x)*0.5)
      fk2 = sqrt((b1 - x)*0.5)
      fw = u0*u0
      x = 1.0/(fw*fw) - b1/fw - c1
 
c---------- 4/2/97 (4)
      if (abs(x) .lt. 1.0E-16) then
         if (x .lt. 0.0) then
            x = -1.0E-6
         else
            x = 1.0E-6
         endif
      endif
c---------- 4/2/97 (4)
 
      fw = 0.5*f0/x
      z(1) = fw*z(1)
      z(2) = fw*z(2)
      z(3) = fw*z(3)
      z(4) = fw*z(4)
      z1(1) = 0.5*(z(1)+z(3))
      z1(2) = 0.5*(z(2)+z(4))
      z1(3) = 0.5*(z(2)-z(4))
      z1(4) = 0.5*(z(1)-z(3))
      a2 = (fk1*fk1 - a(2,2,1))/a(2,1,1)
      b2 = (fk2*fk2 - a(2,2,1))/a(2,1,1)
      x = b(1,1)*b(4,1) - b(3,1)*b(2,1)
      fw1 = fk1/x
      fw2 = fk2/x
      y = fw2*(b2*b(2,1)-b(4,1))
      zx = fw1*(a2*b(2,1)-b(4,1))
      a1(1,1) = 0.5*(1 - y)
      a1(1,2) = 0.5*(1 - zx)
      a1(1,3) = 0.5*(1 + zx)
      a1(1,4) = 0.5*(1 + y)
      y = fw2*(b(3,1)-b2*b(1,1))
      zx = fw1*(b(3,1)-a2*b(1,1))
      a1(2,1) = 0.5*(b2 - y)
      a1(2,2) = 0.5*(a2 - zx)
      a1(2,3) = 0.5*(a2 + zx)
      a1(2,4) = 0.5*(b2 + y)
      a1(3,1) = a1(2,4)
      a1(3,2) = a1(2,3)
      a1(3,3) = a1(2,2)
      a1(3,4) = a1(2,1)
      a1(4,1) = a1(1,4)
      a1(4,2) = a1(1,3)
      a1(4,3) = a1(1,2)
      a1(4,4) = a1(1,1)
 
      if (ib .le. mbs) then
         fq0 = exp((-t0/u0))
         fq1 = exp((-t1/u0))
      else
         fq0 = 1.0
         fq1 = exp((-dt/u0))
      endif
      x = exp((-fk1*dt))
      y = exp((-fk2*dt))
 
      do i = 1, 4
         zz(i,1) = z1(i)*fq0
         zz(i,2) = z1(i)*fq1
         aa(i,1,1) = a1(i,1)
         aa(i,2,1) = a1(i,2)
         aa(i,3,1) = a1(i,3)*x
         aa(i,4,1) = a1(i,4)*y
         aa(i,3,2) = a1(i,3)
         aa(i,4,2) = a1(i,4)
         aa(i,1,2) = a1(i,1)*y
         aa(i,2,2) = a1(i,2)*x
      end do
      return 
      end 
 
      subroutine coefft
c  ********************************************************************
c  See the paper by Liou, Fu and Ackerman (1988) for the formulation of
c  the delta-four-stream approximation in a homogeneous layer.
c  ********************************************************************
      call coeff1
      call coeff2
      call coeff4
      call coeffl
      return 
      end 
 
      subroutine coefft0
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, jj, j, k
      real u,w,w1,w2,w3,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2,y,fw,dt,x
 
      common /point/u(4)
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      common /coedfl/aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
      fk1 = 4.7320545
      fk2 = 1.2679491
      y = exp((-(t1 - t0)/u0))
      fw = 0.5*f0
 
      do i = 1, 4
         if (ib .le. mbs) then
            z1(i) = 0.0
            zz(i,1) = 0.0
            zz(i,2) = 0.0
         else
            jj = 5 - i
            z1(i) = fw/(1.0 + u(jj)/u0)
            zz(i,1) = z1(i)
            zz(i,2) = z1(i)*y
         endif
 
c  ***************************************************************
c  ******************  10/24/96 **********************************
c  ***************************************************************
c  multiple references to 10 continue replaced by loops 11 and 12.
c  ***************************************************************
c  ***************************************************************
c  ***************************************************************
         do j = 1, 4
            a1(i,j) = 0.0
            do k = 1, 2
               aa(i,j,k) = 0.0
            end do
         end do
      end do
 
      do i = 1, 4
         j = 5 - i
         a1(i,j) = 1.0
      end do
 
      dt = t1 - t0
      x = exp((-fk1*dt))
      y = exp((-fk2*dt))
      aa(1,4,1) = y
      aa(2,3,1) = x
      aa(3,2,1) = 1.0
      aa(4,1,1) = 1.0
      aa(1,4,2) = 1.0
      aa(2,3,2) = 1.0
      aa(3,2,2) = x
      aa(4,1,2) = y
 
      return 
      end 
 
      subroutine qccfe(ib, asbs, ee)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i,ib,j,k,n,n4,ibn,i8,kf,i1,i2,i3,j1,j2,j3,m1,m2,m18,m28
      real a, u, w1, w2, w3, w, t, u0, f0, wn, w1n, w2n, w3n, t0n, t1n
      real u0n,f0n,aa,zz,a1,z1,fk1t,fk2t,fk1,fk2,a4,z4,g4,ab,bx,xx
      real fu(4,4), wu(4), v1, v2, v3, asbs, ee, fw1, fw2
 
      common /dis/a(4)
      common /point/u(4)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0(ndfsx), f0(ndfsx)
      common /coedfi/ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
      common /coedfl/aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1t, fk2t
      common /qccfeo/fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), z4(4,ndfsx)
     1   , g4(4,ndfsx)
      common /qcfelc/ab(13,ndfs4x), bx(ndfs4x), xx(ndfs4x)
 
      n = ndfs
      n4 = ndfs4
 
      do i = 1, n4
         do j = 1, 13
            ab(j,i) = 0.0
         end do
      end do
 
      ibn = ib
      wn = w(1)
      w1n = w1(1)
      w2n = w2(1)
      w3n = w3(1)
      t0n = 0.0
      t1n = t(1)
      u0n = u0(1)
      f0n = f0(1)
      if (wn .ge. 0.999999) wn = 0.999999
      if (wn .le. 1.0E-4) then
         call coefft0
         fk1(1) = fk1t
         fk2(1) = fk2t
      else
         call coefft
         fk1(1) = fk1t
         fk2(1) = fk2t
      endif
 
      do i = 1, 4
         z4(i,1) = z1(i)
         do j = 1, 4
            a4(i,j,1) = a1(i,j)
         end do
      end do
 
      do i = 1, 2
         bx(i) = -zz(i+2,1)
         i8 = i + 8
         do j = 1, 4
            ab(i8-j,j) = aa(i+2,j,1)
         end do
      end do
 
      do i = 1, 4
         wu(i) = zz(i,2)
         do j = 1, 4
            fu(i,j) = aa(i,j,2)
         end do
      end do
 
      do k = 2, n
         wn = w(k)
         w1n = w1(k)
         w2n = w2(k)
         w3n = w3(k)
         t0n = t(k-1)
         t1n = t(k)
         u0n = u0(k)
         f0n = f0(k)
         if (wn .ge. 0.999999) wn = 0.999999
         if (wn .le. 1.0E-4) then
            call coefft0
            fk1(k) = fk1t
            fk2(k) = fk2t
         else
            call coefft
            fk1(k) = fk1t
            fk2(k) = fk2t
         endif
 
         do i = 1, 4
            z4(i,k) = z1(i)
            do j = 1, 4
               a4(i,j,k) = a1(i,j)
            end do
         end do
 
         kf = k + k + k + k
         i1 = kf - 5
         i2 = i1 + 3
         j1 = kf - 7
         j2 = j1 + 3
         i3 = 0
         do i = i1, i2
            i3 = i3 + 1
            bx(i) = (-wu(i3)) + zz(i3,1)
            j3 = 0
            i8 = i + 8
            do j = j1, j2
               j3 = j3 + 1
               ab(i8-j,j) = fu(i3,j3)
            end do
            j3 = 0
            do j = j2 + 1, j2 + 4
               j3 = j3 + 1
               ab(i8-j,j) = -aa(i3,j3,1)
            end do
         end do
 
         do i = 1, 4
            wu(i) = zz(i,2)
            do j = 1, 4
               fu(i,j) = aa(i,j,2)
            end do
         end do
 
      end do
 
      if (ib .le. mbs) then
         v1 = 0.2113247*asbs
         v2 = 0.7886753*asbs
         v3 = asbs*u0(1)*f0(1)*exp((-t(n)/u0(1)))
         m1 = n4 - 1
         m2 = n4
         m18 = m1 + 8
         m28 = m2 + 8
         fw1 = v1*wu(3)
         fw2 = v2*wu(4)
         bx(m1) = -(wu(1)-fw1-fw2-v3)
         bx(m2) = -(wu(2)-fw1-fw2-v3)
         do j = 1, 4
            j1 = n4 - 4 + j
            fw1 = v1*fu(3,j)
            fw2 = v2*fu(4,j)
            ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
            ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
         end do
      else
         v1 = 0.2113247*(1.0 - ee)
         v2 = 0.7886753*(1.0 - ee)
         v3 = asbs
         m1 = n4 - 1
         m2 = n4
         m18 = m1 + 8
         m28 = m2 + 8
         fw1 = v1*wu(3)
         fw2 = v2*wu(4)
         bx(m1) = -(wu(1)-fw1-fw2-v3)
         bx(m2) = -(wu(2)-fw1-fw2-v3)
         do j = 1, 4
            j1 = n4 - 4 + j
            fw1 = v1*fu(3,j)
            fw2 = v2*fu(4,j)
            ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
            ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
         end do
      endif
 
      call qcfel
      do k = 1, n
         j = k + k + k + k - 4
         do i = 1, 4
            j = j + 1
            g4(i,k) = xx(j)
         end do
      end do
      return 
      end 
 
      subroutine qcfel
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer n,n4,k,l,i,j,k44,m1,i0m1,i0,m18,i0f,n44,m,im1,n3,n2,n1
      integer m2, m3, m4, m1f, m28, m38, m48, ifq
      real ab, b, x, p, t, yy, xx
 
      common /qcfelc/ab(13,ndfs4x), b(ndfs4x), x(ndfs4x)
      n = ndfs
      n4 = ndfs4
 
      do k = 1, n - 1
         k44 = 4*k - 4
         do l = 1, 4
            m1 = k44 + l
            p = 0.0
 
            do i = 8, 14 - l
               if (abs(ab(i,m1)) .gt. abs(p)) then
                  p = ab(i,m1)
                  i0 = i
               endif
            end do
            i0m1 = i0 + m1
            m18 = m1 + 8
 
c       **********************************************************
c       ****************** 10/24/96 ******************************
c       **********************************************************
c       Replaced "if i0=8 goto 20" statement with if "i0<>8 then".
c       line "20 replaced by the endif.
c       **********************************************************
c       **********************************************************
c       **********************************************************
            if (i0 .ne. 8) then
               do j = m1, m1 + 8 - l
                  i0f = i0m1 - j
                  m1f = m18 - j
                  t = ab(i0f,j)
                  ab(i0f,j) = ab(m1f,j)
                  ab(m1f,j) = t
               end do
 
               i0f = i0m1 - 8
               t = b(i0f)
               b(i0f) = b(m1)
               b(m1) = t
            endif
 
            yy = ab(8,m1)
            ab(8,m1) = 1.0
            do j = m1 + 1, m1 + 8 - l
               m1f = m18 - j
               ab(m1f,j) = ab(m1f,j)/yy
            end do
            b(m1) = b(m1)/yy
 
            do i = 9, 14 - l
               xx = ab(i,m1)
               ab(i,m1) = 0.0
               im1 = i + m1
               do j = m1 + 1, m1 + 8 - l
                  ifq = im1 - j
                  m1f = m18 - j
                  ab(ifq,j) = ab(ifq,j) - ab(m1f,j)*xx
               end do
               ifq = im1 - 8
               b(ifq) = b(ifq) - b(m1)*xx
            end do
 
         end do
      end do
 
      n44 = n4 - 4
      do l = 1, 3
         m1 = n44 + l
         p = 0.0
         do i = 8, 12 - l
            if (abs(ab(i,m1)) .gt. abs(p)) then
               p = ab(i,m1)
               i0 = i
            endif
         end do
         i0m1 = i0 + m1
         m18 = m1 + 8
 
c      **********************************************************
c      ****************** 10/24/96 ******************************
c      **********************************************************
c      Replaced "if i0=8 goto 55" statement with if "i0<>8 then".
c      line "55 replaced by the endif.
c      **********************************************************
c      **********************************************************
c      **********************************************************
         if (i0 .ne. 8) then
            do j = m1, m1 + 4 - l
               i0f = i0m1 - j
               m1f = m18 - j
               t = ab(i0f,j)
               ab(i0f,j) = ab(m1f,j)
               ab(m1f,j) = t
            end do
            i0f = i0m1 - 8
            t = b(i0f)
            b(i0f) = b(m1)
            b(m1) = t
         endif
 
         yy = ab(8,m1)
         ab(8,m1) = 1.0
         do j = m1 + 1, m1 + 4 - l
            m1f = m18 - j
            ab(m1f,j) = ab(m1f,j)/yy
         end do
         b(m1) = b(m1)/yy
         do i = 9, 12 - l
            xx = ab(i,m1)
            ab(i,m1) = 0.0
            im1 = i + m1
            do j = m1 + 1, m1 + 4 - l
               ifq = im1 - j
               m1f = m18 - j
               ab(ifq,j) = ab(ifq,j) - ab(m1f,j)*xx
            end do
            ifq = im1 - 8
            b(ifq) = b(ifq) - b(m1)*xx
         end do
      end do
 
      yy = ab(8,n4)
      ab(8,n4) = 1.0
      b(n4) = b(n4)/yy
      n3 = n4 - 1
      n2 = n3 - 1
      n1 = n2 - 1
      x(n4) = b(n4)
      x(n3) = b(n3) - ab(7,n4)*x(n4)
      x(n2) = b(n2) - ab(7,n3)*x(n3) - ab(6,n4)*x(n4)
      x(n1) = b(n1) - ab(7,n2)*x(n2) - ab(6,n3)*x(n3) - ab(5,n4)*x(n4)
      do k = 1, n - 1
         m4 = 4*(n - k)
         m3 = m4 - 1
         m2 = m3 - 1
         m1 = m2 - 1
         m48 = m4 + 8
         m38 = m3 + 8
         m28 = m2 + 8
         m18 = m1 + 8
         x(m4) = b(m4)
         do m = m4 + 1, m4 + 4
            x(m4) = x(m4) - ab(m48-m,m)*x(m)
         end do
         x(m3) = b(m3)
         do m = m3 + 1, m3 + 5
            x(m3) = x(m3) - ab(m38-m,m)*x(m)
         end do
         x(m2) = b(m2)
         do m = m2 + 1, m2 + 6
            x(m2) = x(m2) - ab(m28-m,m)*x(m)
         end do
         x(m1) = b(m1)
         do m = m1 + 1, m1 + 7
            x(m1) = x(m1) - ab(m18-m,m)*x(m)
         end do
      end do
      return 
      end 
 
      subroutine adjust4
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, n
      real ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t, u0a, f0a
      real dtt(ndfsx), dt(ndfsx),tt0,f,fw
 
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0a(ndfsx), f0a(ndfsx)
 
      n = ndfs
      tt0 = 0.0
      do i = 1, n
         f = ww4(i)/9.0
         fw = 1.0 - f*ww(i)
         w1(i) = (ww1(i)-3.0*f)/(1.0 - f)
         w2(i) = (ww2(i)-5.0*f)/(1.0 - f)
         w3(i) = (ww3(i)-7.0*f)/(1.0 - f)
         w(i) = (1.0 - f)*ww(i)/fw
         dtt(i) = tt(i) - tt0
         tt0 = tt(i)
         dt(i) = dtt(i)*fw
      end do
 
      t(1) = dt(1)
      do i = 2, n
         t(i) = dt(i) + t(i-1)
      end do
 
      return 
      end 
 
      subroutine qfts(ib, as, u0, f0)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ib, n, m, i, k, jj, ii
      real a, u, ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t, u0a, f0a
      real fk1, fk2, a4, z4, g4, ffu, ffd, as, u0, f0, asbs
      real x(4), fi(4), ee, fw1, fw2, fw3, y, y1, fw4
c---------- 10/28/96 (5)
      real ffdr, ffdf, ydr
      common /dirdiff/ffdr(nv1x), ffdf(nv1x)
c---------- 10/28/96 (5)
 
      common /dis/a(4)
      common /point/u(4)
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0a(ndfsx), f0a(ndfsx)
      common /qccfeo/fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), z4(4,ndfsx)
     1   , g4(4,ndfsx)
      common /dfsout/ffu(mdfsx), ffd(mdfsx)
 
      n = ndfs
      m = mdfs
      ee = 0.0
      asbs = as
      call adjust4
      do i = 1, n
         u0a(i) = u0
         f0a(i) = f0
      end do
 
      call qccfe (ib, asbs, ee)
      fw1 = 0.6638961
      fw2 = 2.4776962
      fw3 = u0*3.14159*f0
      do i = 1, m
         if (i .eq. 1) then
            x(1) = 1.0
            x(2) = 1.0
            x(3) = exp((-fk1(1)*t(1)))
            x(4) = exp((-fk2(1)*t(1)))
            k = 1
            y = 1.0
c---------- 10/28/96 (6a)
            ydr = 1.0
c---------- 10/28/96 (6a)
         else if (i .eq. 2) then
            x(1) = exp((-fk2(1)*t(1)))
            x(2) = exp((-fk1(1)*t(1)))
            x(3) = 1.0
            x(4) = 1.0
            k = 1
            y = exp((-t(1)/u0))
c---------- 10/28/96 (6b)
            ydr = exp((-tt(1)/u0))
c---------- 10/28/96 (6b)
         else
            k = i - 1
            y1 = t(k) - t(k-1)
            x(1) = exp((-fk2(k)*y1))
            x(2) = exp((-fk1(k)*y1))
            x(3) = 1.0
            x(4) = 1.0
            y = exp((-t(k)/u0))
c---------- 10/28/96 (6c)
            ydr = exp((-tt(k)/u0))
c---------- 10/28/96 (6c)
         endif
         do jj = 1, 4
            fi(jj) = z4(jj,k)*y
         end do
 
         do ii = 1, 4
            fw4 = g4(ii,k)*x(ii)
            do jj = 1, 4
               fi(jj) = fi(jj) + a4(jj,ii,k)*fw4
            end do
         end do
 
         ffu(i) = fw1*fi(2) + fw2*fi(1)
         ffd(i) = fw1*fi(3) + fw2*fi(4) + fw3*y
c---------- 10/28/96 (7)
         ffdr(i) = fw3*ydr
         ffdf(i) = ffd(i) - ffdr(i)
c---------- 10/28/96 (7)
      end do
      return 
      end 
 
      subroutine qfti(ib, ee)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, ib, n, k, ii, jj, m
      real a,u,ww1,ww2,ww3,ww4,ww,tt,w1,w2,w3,w,t,u0,f0,fk1
      real fk2,a4,z4,g4,ffu,ffd,bf,bs,asbs,q1,q2,t0,fw1,fw2,xy
      real x(4), fi(4), fw3, ee, y1
      real timt0
      common /dis/a(4)
      common /point/u(4)
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0(ndfsx), f0(ndfsx)
      common /qccfeo/fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), z4(4,ndfsx)
     1   , g4(4,ndfsx)
      common /dfsout/ffu(mdfsx), ffd(mdfsx)
      common /planci/bf(nv1x), bs
 
      n = ndfs
      m = mdfs
      asbs = bs*ee
      call adjust4
      t0 = 0.0
 
      do i = 1, n
         q1 = alog(bf(i+1)/bf(i))
 
c Paul Stackhouse: July 5, 98
c          q2 = 1.0 / ( t(i) - t0 )
         timt0 = t(i) - t0
         timt0 = amax1(1.0E-12,timt0)
         q2 = 1.0/timt0
c Paul Stackhouse: July 5, 98
 
         f0(i) = 2.0*(1.0 - w(i))*bf(i)
         if (abs(q1) .le. 1.0E-10) then
            u0(i) = -1.0E+10/q2
         else
            u0(i) = -1.0/(q1*q2)
         endif
 
c---------- 4/2/97 (5)
         if (abs(u0(i)) .gt. 4.25E+09) then
            if (u0(i) .lt. 0.0) then
               u0(i) = -4.25E+09
            else
               u0(i) = 4.25E+09
            endif
         endif
c---------- 4/2/97 (5)
 
         t0 = t(i)
      end do
 
      call qccfe (ib, asbs, ee)
      fw1 = 0.6638958
      fw2 = 2.4776962
      do i = 1, m
         if (i .eq. 1) then
            x(1) = 1.0
            x(2) = 1.0
            x(3) = exp((-fk1(1)*t(1)))
            x(4) = exp((-fk2(1)*t(1)))
            k = 1
            xy = 1.0
         else if (i .eq. 2) then
            x(1) = exp((-fk2(1)*t(1)))
            x(2) = exp((-fk1(1)*t(1)))
            x(3) = 1.0
            x(4) = 1.0
            k = 1
            xy = exp((-t(1)/u0(1)))
         else
            k = i - 1
            y1 = t(k) - t(k-1)
            x(1) = exp((-fk2(k)*y1))
            x(2) = exp((-fk1(k)*y1))
            x(3) = 1.0
            x(4) = 1.0
            xy = exp((-y1/u0(k)))
         endif
         do jj = 1, 4
            fi(jj) = z4(jj,k)*xy
         end do
         do ii = 1, 4
            fw3 = g4(ii,k)*x(ii)
            do jj = 1, 4
               fi(jj) = fi(jj) + a4(jj,ii,k)*fw3
            end do
         end do
         ffu(i) = fw1*fi(2) + fw2*fi(1)
         ffd(i) = fw1*fi(3) + fw2*fi(4)
      end do
 
      return 
      end 
 
      subroutine cfgts0(gamma1, gamma2, gamma3, gamma4, ugts1)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ib
      real w, w1, w2, w3, t0, t1, u0, f0, x, y, gamma1, gamma2
      real gamma3, gamma4, z, ugts1
 
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
      if (edding) then
         x = 0.25*w1
         y = w*x
         gamma1 = 1.75 - w - y
         gamma2 = (-0.25) + w - y
         gamma3 = 0.0
         gamma4 = 0.0
         if (ib .le. mbs) then
            gamma3 = 0.5 - x*u0
            gamma4 = 1.0 - gamma3
         endif
         ugts1 = 0.5
      endif
 
      if (quadra) then
         x = 0.866*w
         y = 0.2887*w1
         z = y*w
         gamma1 = 1.732 - x - z
         gamma2 = x - z
         gamma3 = 0.0
         gamma4 = 0.0
         if (ib .le. mbs) then
            gamma3 = 0.5 - y*u0
            gamma4 = 1.0 - gamma3
         endif
         ugts1 = 0.57735
      endif
 
      if (hemisp) then
         x = w*w1/3.0
         gamma1 = 2.0 - w - x
         gamma2 = w - x
         gamma3 = 0.0
         gamma4 = 0.0
         ugts1 = 0.5
      endif
 
      if (mquadr) then
         y = 0.2767*w
         x = y + y + y
         z = y*w1
         gamma1 = 1.66 - x - z
         gamma2 = x - z
         gamma3 = 0.0
         gamma4 = 0.0
         ugts1 = 0.6024
      endif
 
      return 
      end 
 
      subroutine cfgts(lamda,gamma,cadd0,cadd1,cmin0,cmin1,g1g2,fkb)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ib
      real lamda, w, w1, w2, w3, t0, t1, u0, f0, gamma1, gamma2
      real gamma3, gamma4, ugts1, gamma, g1g2, alfa, beta, cadd0, cmin0
      real fw, cadd1, cmin1, fkb, x, z, fq
 
      common /coedfi/ib, w, w1, w2, w3, t0, t1, u0, f0
 
      call cfgts0 (gamma1, gamma2, gamma3, gamma4, ugts1)
      lamda = sqrt((gamma1 + gamma2)*(gamma1 - gamma2))
      gamma = gamma2/(gamma1 + lamda)
      g1g2 = gamma1 + gamma2
      fq = 1.0/u0
 
      x = exp((-fq*(t1 - t0)))
      z = lamda*lamda - fq*fq
 
c---------- 4/2/97 (6)
      if (z .eq. 0.0) z = 0.001
c---------- 4/2/97 (6)
 
      fkb = z
 
      if (ib .le. mbs) then
         alfa = gamma3
         beta = gamma4
         fw = 3.1415927*f0*w*exp((-fq*t0))
         cadd0 = fw*((gamma1 - fq)*alfa + beta*gamma2)/z
         cmin0 = fw*((gamma1 + fq)*beta + alfa*gamma2)/z
      else
         fw = 3.1415927*f0
         cadd0 = fw*(g1g2 - fq)/z
         cmin0 = fw*(g1g2 + fq)/z
      endif
 
      cadd1 = cadd0*x
      cmin1 = cmin0*x
      return 
      end 
 
c +++ PAUSE FIX      subroutine qccgts ( ib, asbs, ee )
      subroutine qccgts(ib, asbs, ee, itrierr)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ndfs2x, ib, k, ibn, k1, k2
      integer itrierr
      parameter (ndfs2x = ndfsx*2)
      real lamdan, w1, w2, w3, w, t, u0, f0, wn, w1n, w2n
      real w3n, t0n, t1n, u0n, f0n, gamman, caddn, cminn, caddn0, cminn0
      real aa, bb, expn, g1g2n, fkbn, ee, asbs, wm1, wm2, rsfc, ssfc
      real a(ndfs2x), b(ndfs2x), c(ndfs2x), r(ndfs2x), u(ndfs2x)
      real xn(ndfsx), yn(ndfsx), zn(ndfsx),gam(ndfs2x)
 
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0(ndfsx), f0(ndfsx)
      common /coedfi/ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
      common /gtscoe/lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx), cminn(
     1   ndfsx), caddn0(ndfsx), cminn0(ndfsx), aa(ndfsx), bb(ndfsx), 
     2   expn(ndfsx), g1g2n(ndfsx), fkbn(ndfsx)
 
      ibn = ib
      do k = 1, ndfs
         wn = w(k)
         w1n = w1(k)
         if (k .eq. 1) then
            t0n = 0.0
         else
            t0n = t(k-1)
         endif
         t1n = t(k)
         u0n = u0(k)
         f0n = f0(k)
         if (wn .ge. 0.999999) wn = 0.999999
         call cfgts (lamdan(k), gamman(k), caddn0(k), caddn(k), cminn0(k
     1      ), cminn(k), g1g2n(k), fkbn(k))
         expn(k) = exp((-lamdan(k)*(t1n-t0n)))
         xn(k) = gamman(k)*expn(k)
         yn(k) = (expn(k)-gamman(k))/(xn(k)-1.0)
         zn(k) = (expn(k)+gamman(k))/(xn(k)+1.0)
      end do
 
      a(1) = 0.0
      b(1) = xn(1) + 1.0
      c(1) = xn(1) - 1.0
      r(1) = -cminn0(1)
      do k = 1, ndfs - 1
         k1 = k + k
         k2 = k + k + 1
         a(k1) = 1.0 + xn(k) - yn(k+1)*(gamman(k)+expn(k))
         b(k1) = 1.0 - xn(k) - yn(k+1)*(gamman(k)-expn(k))
         c(k1) = yn(k+1)*(1.0 + xn(k+1)) - expn(k+1) - gamman(k+1)
         r(k1) = caddn0(k+1) - caddn(k) - yn(k+1)*(cminn0(k+1)-cminn(k))
         a(k2) = gamman(k) - expn(k) - zn(k)*(1.0 - xn(k))
         b(k2) = (-1.0) - xn(k+1) + zn(k)*(expn(k+1)+gamman(k+1))
         c(k2) = zn(k)*(expn(k+1)-gamman(k+1)) - xn(k+1) + 1.0
         r(k2) = cminn0(k+1) - cminn(k) - zn(k)*(caddn0(k+1)-caddn(k))
      end do
      if (ib .le. mbs) then
         rsfc = asbs
         ssfc = 3.1415927*u0(1)*exp((-t(ndfs)/u0(1)))*rsfc*f0(1)
      else
         rsfc = 1.0 - ee
         ssfc = 3.1415927*asbs
      endif
 
      wm1 = 1.0 - rsfc*gamman(ndfs)
      wm2 = xn(ndfs) - rsfc*expn(ndfs)
      a(ndfs2) = wm1 + wm2
      b(ndfs2) = wm1 - wm2
      c(ndfs2) = 0.0
      r(ndfs2) = rsfc*cminn(ndfs) - caddn(ndfs) + ssfc
c +++ PAUSE FIX      call tridag ( a, b, c, r, u, gam, ndfs2 )
      call tridag (a, b, c, r, u, gam, ndfs2, itrierr)
      if (itrierr .ne. 0) return 
      do k = 1, ndfs
         k1 = k + k - 1
         k2 = k + k
         aa(k) = u(k1) + u(k2)
         bb(k) = u(k1) - u(k2)
      end do
      return 
      end 
 
c +++ PAUSE FIX      subroutine tridag ( a, b, c, r, u, gam, n )
      SUBROUTINE TRIDAG(A, B, C, R, U, GAM, N, ITRIERR)
c *******************************************************************
c
c   | b1 c1 0  ...                |   | u1   |   | r1   |
c   | a2 b2 c2 ...                |   | u2   |   | r2   |
c   |          ...                | . | .    | = | .    |
c   |          ... an-1 bn-1 cn-1 |   | un-1 |   | rn-1 |
c   |              0    an   bn   |   | un   |   | rn   |
c
c This  subroutine solves for  a vector U of length N the tridiagonal
c linear set given by above equation. A, B, C and R are input vectors
c and are not modified (Numerical Recipes by Press et al. 1989).
c *******************************************************************
      INTEGER N, J
      INTEGER ITRIERR
      REAL GAM(N), A(N), B(N), C(N), R(N), U(N), BET
      ITRIERR = 0
c +++ PAUSE FIX      if ( b(1) .eq. 0. ) pause
      IF (B(1) .EQ. 0.) THEN
         ITRIERR = 1
         RETURN 
      ENDIF
c  *********************************************************
c  If this happens then you should rewrite your equations as
c  a set of order n-1, with u2 trivially eliminated.
c  *********************************************************
      BET = B(1)
      U(1) = R(1)/BET
 
c  ***************************************
c  Decomposition and forward substitution.
c  ***************************************
      DO J = 2, N
         GAM(J) = C(J-1)/BET
         BET = B(J) - A(J)*GAM(J)
c +++ PAUSE FIX       if ( bet .eq. 0. ) pause
         IF (BET .EQ. 0.) THEN
            ITRIERR = 2
            RETURN 
         ENDIF
c      ---------------------------------------
c      Algorithm fails; see Numerical Recipes.
c      ---------------------------------------
         U(J) = (R(J)-A(J)*U(J-1))/BET
      END DO
 
c  ****************
c  Backsubstitution
c  ****************
      DO J = N - 1, 1, -1
         U(J) = U(J) - GAM(J+1)*U(J+1)
      END DO
 
      RETURN 
      END 
 
c +++ PAUSE FIX      subroutine qftsts ( ib, as, u0, f0 )
      subroutine qftsts(ib, as, u0, f0, itrierr)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer itrierr
      integer ib, n, m, i, k
      real ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t, u0a, f0a
      real lamdan, gamman, caddn, cminn, caddn0, cminn0, aa, bb, expn
      real g1g2n, fkbn, ffu, ffd, as, u0, f0, ee, asbs, fw3, xx
      real yy(ndfsx)
 
c---------- 10/28/96 (8)
      real ffdr, ffdf, yydr(ndfsx)
      common /dirdiff/ffdr(nv1x), ffdf(nv1x)
c---------- 10/28/96 (8)
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0a(ndfsx), f0a(ndfsx)
      common /gtscoe/lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx), cminn(
     1   ndfsx), caddn0(ndfsx), cminn0(ndfsx), aa(ndfsx), bb(ndfsx), 
     2   expn(ndfsx), g1g2n(ndfsx), fkbn(ndfsx)
      common /dfsout/ffu(mdfsx), ffd(mdfsx)
 
      n = ndfs
      m = mdfs
      ee = 0.0
      asbs = as
      call adjust2
 
      do i = 1, n
         u0a(i) = u0
         f0a(i) = f0
      end do
 
c +++ PAUSE FIX      call qccgts ( ib, asbs, ee )
      call qccgts (ib, asbs, ee, itrierr)
      if (itrierr .ne. 0) return 
      fw3 = u0*3.1415927*f0
      do k = 1, ndfs
         yy(k) = exp((-t(k)/u0))
c---------- 10/28/96 (9)
         yydr(k) = exp((-tt(k)/u0))
c---------- 10/28/96 (9)
      end do
 
      xx = aa(1)*expn(1)
      ffu(1) = xx + gamman(1)*bb(1) + caddn0(1)
      ffd(1) = gamman(1)*xx + bb(1) + cminn0(1) + fw3
c---------- 10/28/96 (10)
      ffdr(1) = fw3
      ffdf(1) = ffd(1) - ffdr(1)
c---------- 10/28/96 (10)
      do i = 2, m
         k = i - 1
         xx = bb(k)*expn(k)
         ffu(i) = aa(k) + gamman(k)*xx + caddn(k)
         ffd(i) = gamman(k)*aa(k) + xx + cminn(k) + fw3*yy(k)
c---------- 10/28/96 (11)
         ffdr(i) = fw3*yydr(k)
         ffdf(i) = ffd(i) - ffdr(i)
c---------- 10/28/96 (11)
      end do
 
      return 
      end 
 
c +++ PAUSE FIX      subroutine qftits ( ib, ee )
      subroutine qftits(ib, ee, itrierr)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer itrierr
      integer ib, n, m, i, k
      real lamdan, ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t
      real gamman, caddn, cminn, caddn0, cminn0, aa, bb, expn
      real g1g2n,fkbn,ffu,ffd,u0,f0,ee,asbs,xx,q1,q2,bf,bs,t0
      real timt0
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0(ndfsx), f0(ndfsx)
      common /gtscoe/lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx), cminn(
     1   ndfsx), caddn0(ndfsx), cminn0(ndfsx), aa(ndfsx), bb(ndfsx), 
     2   expn(ndfsx), g1g2n(ndfsx), fkbn(ndfsx)
      common /dfsout/ffu(mdfsx), ffd(mdfsx)
      common /planci/bf(nv1x), bs
      n = ndfs
      m = mdfs
      asbs = bs*ee
      call adjust2
      t0 = 0.0
 
      do i = 1, n
         q1 = alog(bf(i+1)/bf(i))
 
c Paul Stackhouse: July 5, 98
c          q2 = 1.0 / ( t(i) - t0 )
         timt0 = t(i) - t0
         timt0 = amax1(1.0E-12,timt0)
         q2 = 1.0/timt0
c Paul Stackhouse: July 5, 98
 
         if (mquadr) then
            f0(i) = 1.66*(1.0 - w(i))*bf(i)
         else
            f0(i) = 2.0*(1.0 - w(i))*bf(i)
         endif
 
         if (abs(q1) .le. 1.0E-10) then
            u0(i) = -1.0E+10/q2
         else
            u0(i) = -1.0/(q1*q2)
         endif
         t0 = t(i)
      end do
 
c +++ PAUSE FIX      call qccgts ( ib, asbs, ee )
      call qccgts (ib, asbs, ee, itrierr)
      if (itrierr .ne. 0) return 
      xx = aa(1)*expn(1)
      ffu(1) = xx + gamman(1)*bb(1) + caddn0(1)
      ffd(1) = gamman(1)*xx + bb(1) + cminn0(1)
 
      do i = 2, m
         k = i - 1
         xx = bb(k)*expn(k)
         ffu(i) = aa(k) + gamman(k)*xx + caddn(k)
         ffd(i) = gamman(k)*aa(k) + xx + cminn(k)
      end do
      return 
      end 
 
c 6-24-98 (8)
c     subroutine qftisf ( ib, ee )
c +++ PAUSE FIX subroutine qftisf ( ib, ee, ur )
      subroutine qftisf(ib, ee, ur, itrierr)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer itrierr
      integer ib, nq, n, m, i, j, i1
      parameter (nq = 2)
      real ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t
      real lamdan, gamman, caddn, cminn, caddn0, cminn0, aa, bb, expn
      real g1g2n, fkbn, ffu, ffd, u0, f0, ee, asbs, xx, bf, bs, ugts1
      real t0, q1, q2, x, y, z, y1, yy
      real fg(ndfsx), fh(ndfsx), fj(ndfsx), fk(ndfsx)
      real alfa(ndfsx+1), beta(ndfsx)
      real fiu(mdfsx,nq), fid(mdfsx,nq),ub(ndfsx,nq)
      real fx(ndfsx,nq), fy(ndfsx), fz1(ndfsx,nq), fz2(ndfsx,nq)

      real fuq1(ndfsx), fuq2(ndfsx)
      real ug(nq), wg(nq), ugwg(nq)
 
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0(ndfsx), f0(ndfsx)
      common /gtscoe/lamdan(ndfsx), gamman(ndfsx), caddn(ndfsx), cminn(
     1   ndfsx), caddn0(ndfsx), cminn0(ndfsx), aa(ndfsx), bb(ndfsx), 
     2   expn(ndfsx), g1g2n(ndfsx), fkbn(ndfsx)
      common /dfsout/ffu(mdfsx), ffd(mdfsx)
c 6-24-98 (9/10a)
      real fiurt, fiurw, fiur, ur
      real timt0
      real  fxr(ndfsx), fz1r(ndfsx), fz2r(ndfsx), ubr(ndfsx)
      common /radiance/ fiurt(nv1x), fiurw(nv1x), fiur(nv1x)
c 6-24-98 (9/10a)
      common /planci/bf(nv1x), bs
      data ug/ 0.2113248, 0.7886752/ 
      data wg/ 0.5, 0.5/ 
      data ugwg/ 0.105662, 0.394338/ 
      if (mquadr) then
         ugts1 = 0.60241
      else
         ugts1 = 0.5
      endif
      n = ndfs
      m = mdfs
      asbs = bs*ee
      call adjust2
      t0 = 0.0
 
      do i = 1, n
         q1 = alog(bf(i+1)/bf(i))
 
c Paul Stackhouse: July 5, 98
c          q2 = 1.0 / ( t(i) - t0 )
         timt0 = t(i) - t0
         timt0 = amax1(1.0E-12,timt0)
         q2 = 1.0/timt0
c Paul Stackhouse: July 5, 98
 
         if (mquadr) then
            f0(i) = 1.66*(1.0 - w(i))*bf(i)
         else
            f0(i) = 2.0*(1.0 - w(i))*bf(i)
         endif
         if (abs(q1) .le. 1.0E-10) then
            u0(i) = -1.0E+10/q2
         else
            u0(i) = -1.0/(q1*q2)
         endif
         t0 = t(i)
         beta(i) = -1.0/u0(i)
      end do
 
c +++ PAUSE FIX      call qccgts ( ib, asbs, ee )
      call qccgts (ib, asbs, ee, itrierr)
      if (itrierr .ne. 0) return 
      do i = 1, n
         x = (1.0 - w(i))*w(i)/fkbn(i)/ugts1
         y1 = w1(i)/3.0
         y = g1g2n(i)
         z = -y1*beta(i)
         fuq1(i) = x*(y - z) + 1.0 - w(i)
         fuq2(i) = x*(y + z) + 1.0 - w(i)
      end do
 
      do i = 1, n + 1
         alfa(i) = 6.2832*bf(i)
      end do
 
      do i = 1, n
         x = lamdan(i)*ugts1
         y = gamman(i)*(1.0 + x)
         z = 1.0 - x
         fg(i) = aa(i)*(z + z)
         fh(i) = bb(i)*(y + y)
         fj(i) = aa(i)*(y + y)
         fk(i) = bb(i)*(z + z)
      end do
 
      do j = 1, 2
         fid(1,j) = 0.0
      end do
 
      do j = 1, 2
         t0 = 0.0
         do i = 2, mdfs
            i1 = i - 1
            fx(i1,j) = exp((-(t(i1)-t0)/ug(j)))
            fy(i1) = expn(i1)
            xx = lamdan(i1)*ug(j)
            fz1(i1,j) = (1.0 - fx(i1,j)*fy(i1))/(xx + 1.0)
            fz2(i1,j) = (fx(i1,j)-fy(i1))/(xx - 1.0)
            ub(i1,j) = ug(j)*beta(i1)
 
c---------- 4/2/97 (7)
            if (ub(i1,j) .eq. 1.0) ub(i1,j) = 1.001
c---------- 4/2/97 (7)
 
 
            fid(i,j) = fid(i1,j)*fx(i1,j) + fj(i1)*fz1(i1,j) + fk(i1)*
     1         fz2(i1,j) + fuq2(i1)/(ub(i1,j)+1.0)*(alfa(i)-alfa(i1)*fx(
     2         i1,j))
            t0 = t(i1)
         end do
      end do
 
      yy = 0.0
      do j = 1, 2
         yy = yy + ugwg(j)*fid(mdfs,j)
      end do
 
      xx = yy*(1.0 - ee)*2.0 + 6.2831854*ee*bs
      do j = 1, 2
         fiu(mdfs,j) = xx
      end do
c 6-24-98 (11)
      fiur(mdfs) = xx
c 6-24-98 (11)
      do j = 1, 2
         do i = mdfs - 1, 1, -1
            fiu(i,j) = fiu(i+1,j)*fx(i,j) + fg(i)*fz2(i,j) + fh(i)*fz1(i
     1         ,j) + fuq1(i)/(ub(i,j)-1.0)*(alfa(i+1)*fx(i,j)-alfa(i))
         end do
      end do
 
      do i = 1, mdfs
         ffu(i) = 0.0
         ffd(i) = 0.0
      end do
 
      do i = 1, mdfs
         do j = 1, 2
            ffu(i) = ffu(i) + ugwg(j)*fiu(i,j)
            ffd(i) = ffd(i) + ugwg(j)*fid(i,j)
         end do
      end do
      if (ur .eq. 0.5) ur = 0.500001
      t0 = 0.0
      do i = 2, mdfs
         i1 = i - 1
         fxr(i1) = exp((-(t(i1)-t0)/ur))
         xx = lamdan(i1)*ur
         fz1r(i1) = (1.0 - fxr(i1)*fy(i1))/(xx + 1.0)
         fz2r(i1) = (fxr(i1)-fy(i1))/(xx - 1.0)
         ubr(i1) = ur*beta(i1)
         t0 = t(i1)
      end do
      do i = mdfs - 1, 1, -1
         fiur(i) = fiur(i+1)*fxr(i) + fg(i)*fz2r(i) + fh(i)*fz1r(i) + 
     1      fuq1(i)/(ubr(i)-1.0)*(alfa(i+1)*fxr(i)-alfa(i))
      end do
!     else
!     do i = 1, mdfs
!           fiur(i) = 0.0
!     enddo
!     endif
c 6-24-98 (12)
      return 
      end 
 
      subroutine adjust2
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, n
      real ww1, ww2, ww3, ww4, ww, tt, w1, w2, w3, w, t, u0a
      real f0a, dtt(ndfsx), dt(ndfsx), tt0, f, fw
      common /dfsin/ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx), ww(
     1   ndfsx), tt(ndfsx)
      common /qccfei/w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), t(ndfsx)
     1   , u0a(ndfsx), f0a(ndfsx)
 
      n = ndfs
      tt0 = 0.0
      do i = 1, n
         f = ww2(i)/5.0
         fw = 1.0 - f*ww(i)
         w1(i) = (ww1(i)-3.0*f)/(1.0 - f)
         w(i) = (1.0 - f)*ww(i)/fw
         dtt(i) = tt(i) - tt0
         tt0 = tt(i)
         dt(i) = dtt(i)*fw
      end do
      t(1) = dt(1)
 
      do i = 2, n
         t(i) = dt(i) + t(i-1)
      end do
 
      return 
      end 
 
c---------- 4/2/97 (3) -- NEXT 801 LINES
      subroutine gascon_ckd_parm(ib, h2o_cont)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer ib, h2o_cont
      integer m 
      real log_average
      real pp, pt, pa, ph, po, tgm, dz
      real dp, amnt, patm, temp, parm_ckd
      integer iflb(18)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /con/tgm(nvx)
      common /thick/dz(nvx)
      data iflb/ 6*0, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1/ 
 
      do m = 1, nv
         tgm(m) = 0.0
      end do
 
      if (iflb(ib) .eq. 0) return 
 
      do m = 1, nv
cbm      dp = pp(m+1) - pp(m)
cbm      amnt = 1.02*dp*(ph(m)+ph(m+1))*0.5
         amnt = pa(m)*0.1*(ph(m)+ph(m+1))*0.5
cbm logarithmic interpolation of pressure instead of linear
cbm ??? this should be made consistent with uvspec ???
cbm      patm = ((pp(m)+pp(m+1))*0.5)/1013.25
         patm = log_average(pp(m), pp(m+1))/1013.25 
         temp = (pt(m)+pt(m+1))*0.5
cbm         dz = dz_quick(pp(m),pp(m+1),pt(m),pt(m+1),ph(m),ph(m+1))
         if (h2o_cont .eq. 1) then 
           tgm(m)=parm_ckd(iflb(ib),amnt,patm,temp,dz(m))
         elseif (h2o_cont .eq. 2) then 
           tgm(m)=parm_ckd24(iflb(ib),amnt,patm,temp,dz(m))
         endif
c         tgm(m)=parm_ckd(iflb(ib),amnt,patm,temp,dz(m))
      end do
 
      return 
      end 

c --------------------------------------------
cbm      real function dz_quick (p1, p2, t1, t2, h1, h2)
cbm      temp = (t1 + t2)*0.5
cbm      hum = (h1 + h2)*0.5
cbm      tv = temp*(1 + .61*hum)
cbm      dz_quick = 29.3*tv*log(p2/p1)*0.001
cbm      return 
cbm      end 
c --------------------------------------------
c*********************************************************************
      real function parm_ckd (iband, amnt, patm, temp, dz)
c Parameterization of CKD_2.1 continuum over Fu-Liou Bands
c Input:
c iband  =  integer (1-8) where
c       Band 1 ='  5:280cm-1'
c         Band 2 ='280:400cm-1'
c         Band 3 ='400:540cm-1'
c         Band 4 ='540:670cm-1'
c         Band 5 ='670:800cm-1'
c         Band 6 ='800:980cm-1'
c         Band 7 ='980:1100cm-1'
c         Band 8 ='1100:1250cm-1'
c         Band 9 ='1250:1400cm-1'
c         Band10 ='1400:1700cm-1'
c         Band11 ='1700:1900cm-1'
c         Band12 ='1900:2200cm-1'
c amnt = h2O ammount (g/cm**2)
c patm= pressure (atm)
c temp = temperature (k)
c dz = pathlength (Km)
c Output:
c clough_parm = parameterized CKD_2.1optical depth for band
c234567890123456789012345678901234567890123456789012345678901234567890
 
      parameter (ncoef = 7, nband = 12)
      real aa(ncoef,nband), aa1, aa2, aa3, aa4
      common /pcont_aa0/aa0(ncoef,nband)         ! ORIGINAL 280:1400
      common /pcont_aa1/aa1(ncoef,nband)         ! lin no plank
      common /pcont_aa2/aa2(ncoef,nband)         ! log no plank
      common /pcont_aa3/aa3(ncoef,nband)         ! lin Plank wgt
      common /pcont_aa4/aa4(ncoef,nband)         ! log Plank wgt
      common /cont_tas/iwtas
      integer j1, j2
      if (iwtas.lt.0 .or. iwtas.gt.4) stop 'iwtas'
      if (iwtas .eq. 0) then
         do j2 = 1, 12
            do j1 = 1, 7
               aa(j1,j2) = aa0(j1,j2)
            end do
         end do
      endif
      if (iwtas .eq. 1) then
         do j2 = 1, 12
            do j1 = 1, 7
               aa(j1,j2) = aa1(j1,j2)
            end do
         end do
      endif
      if (iwtas .eq. 2) then
         do j2 = 1, 12
            do j1 = 1, 7
               aa(j1,j2) = aa2(j1,j2)
            end do
         end do
      endif
      if (iwtas .eq. 3) then
         do j2 = 1, 12
            do j1 = 1, 7
               aa(j1,j2) = aa3(j1,j2)
            end do
         end do
      endif
      if (iwtas .eq. 4) then
         do j2 = 1, 12
            do j1 = 1, 7
               aa(j1,j2) = aa4(j1,j2)
            end do
         end do
      endif
 
      ph2o = amnt*(8.314D+07*temp)/(dz*1.0D+05*18.01534*1.01325D+06)
 
      if (iwtas .eq. 0) patmx = patm
      if (iwtas .gt. 0) patmx = log(patm)
 
 
      tau_log = aa(1,iband) + aa(2,iband)*log(amnt) + aa(3,iband)*temp
     1    + aa(4,iband)*patmx + aa(5,iband)*ph2o + aa(6,iband)*amnt + aa
     2   (7,iband)*log(ph2o)
 
 
      parm_ckd = exp(tau_log)
 
      return 
      end 
 
      subroutine gascon_off
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i
      real tgm
      common /con/tgm(nvx)
 
      do i = 1, nv
         tgm(i) = 0.0
      end do
 
      return 
      end 
 


c=====================================================================
        function parm_ckd24(iband,amnt,patm,temp,dz)
c Parameterization of CKD_2.4 continuum over Fu-Liou Bands
c Input:
c iband  =  integer (1-12) where
c	  Band 1 ='  5:280cm-1'
c         Band 2 ='280:400cm-1'
c         Band 3 ='400:540cm-1'
c         Band 4 ='540:670cm-1'
c         Band 5 ='670:800cm-1'
c         Band 6 ='800:980cm-1'
c         Band 7 ='980:1100cm-1'
c         Band 8 ='1100:1250cm-1'
c         Band 9 ='1250:1400cm-1'
c         Band10 ='1400:1700cm-1'
c         Band11 ='1700:1900cm-1'
c         Band12 ='1900:2200cm-1'
c amnt = h2O ammount (g/cm**2)
c patm = pressure (atm)
c temp = temperature (k)
c dz   = pathlength (Km)
c Output:
c parm_ckd24 = parameterized CKD_2.4optical depth for band
c234567890123456789012345678901234567890123456789012345678901234567890
	parameter (ncoef=7,nreg=2,nband=12)

	real  h2obnd(nband)
	real ck24_3(ncoef,nreg,nband)! lin Plank wgt

	data h2obnd /-5,-3.5,-2.0,-2,-1.,-4,-4,-4,-3,-3.5,-3,-2/

       data ck24_3/  !ckd24fu.fuliou.lin.plnk.out
! band        1
     x 1.667e+00, 9.421e-01,-7.358e-03, 1.355e+00,
     x 2.557e+03, 5.798e+01,-4.570e-01,
! band        1
     x 6.417e+00, 1.002e+00,-6.991e-03, 1.010e+00,
     x 1.203e+01, 4.501e-02,-2.428e-02,
! band        2
     x 2.390e+00, 9.528e-01,-6.058e-03, 1.071e+00,
     x 2.676e+02, 9.848e+00,-1.459e-01,
! band        2
     x 4.849e+00, 1.002e+00,-6.910e-03, 8.961e-01,
     x 1.635e+01, 2.115e-02, 7.243e-02,
! band        3
     x 2.326e+00, 9.720e-01,-6.551e-03, 8.739e-01,
     x 6.984e+01, 8.346e-01, 4.824e-02,
! band        3
     x 5.002e+00, 1.005e+00,-9.286e-03, 6.222e-01,
     x 1.168e+01, 3.611e-03, 3.148e-01,
! band        4
     x-4.865e+00, 8.455e-01,-6.911e-03, 1.475e+00,
     x 2.905e+02, 7.078e+00,-6.846e-01,
! band        4
     x 4.596e+00, 1.012e+00,-1.152e-02, 5.713e-01,
     x 1.270e+01,-1.395e-03, 3.447e-01,
! band        5
     x-5.396e+00, 8.596e-01,-8.479e-03, 1.619e+00,
     x 1.664e+02, 3.236e+00,-7.782e-01,
! band        5
     x 7.478e+00, 1.007e+00,-1.963e-02, 2.771e-01,
     x 6.021e+00,-4.489e-03, 6.709e-01,
! band        6
     x 1.262e+00, 2.347e-01,-2.360e-02, 1.655e-01,
     x 5.068e+02, 2.462e+01, 3.920e-01,
! band        6
     x 9.334e+00, 1.002e+00,-2.429e-02, 3.575e-02,
     x 2.751e-01,-1.189e-03, 9.593e-01,
! band        7
     x-1.222e+00, 5.423e-01,-2.327e-02, 5.197e-01,
     x 6.423e+02, 5.038e+01, 1.502e-01,
! band        7
     x 8.506e+00, 1.000e+00,-2.339e-02, 8.891e-03,
     x-6.805e-01,-1.639e-04, 9.917e-01,
! band        8
     x-3.638e+00, 8.534e-01,-1.344e-02, 6.816e-01,
     x 5.385e+02, 4.428e+01,-6.366e-03,
! band        8
     x 6.921e+00, 1.002e+00,-1.974e-02, 6.350e-02,
     x 6.838e-01,-1.121e-03, 9.237e-01,
! band        9
     x-2.329e+00, 7.893e-01,-2.588e-03, 1.017e+00,
     x 1.525e+02, 1.029e+01,-1.486e-01,
! band        9
     x 6.742e-01, 1.008e+00,-3.376e-03, 9.105e-01,
     x 1.074e+01,-3.307e-03, 5.741e-02,
! band       10
     x-1.677e+00, 9.173e-01,-5.780e-03, 1.504e+00,
     x 7.886e+02, 2.288e+01,-5.999e-01,
! band       10
     x 3.396e+00, 1.005e+00,-3.433e-03, 1.012e+00,
     x 7.635e+00, 3.010e-03,-2.418e-02,
! band       11
     x 7.943e-01, 9.260e-01,-5.050e-03, 1.141e+00,
     x 2.221e+02, 1.021e+01,-2.246e-01,
! band       11
     x 3.356e+00, 1.002e+00,-4.719e-03, 9.578e-01,
     x 6.164e+00, 1.186e-03, 2.264e-02,
! band       12
     x-5.874e+00, 7.060e-01,-1.532e-03, 1.141e+00,
     x 1.463e+02, 6.534e+00,-4.308e-01,
! band       12
     x 4.709e-01, 1.010e+00,-6.067e-03, 8.513e-01,
     x 1.161e+01,-6.629e-03, 8.885e-02
     &/

	
! These Regressions are more sensitive to pathlength
! So accomodations for very Thin or Thick layers are made.
	 dz1  =dz
	 factor=1.000
	   if ( dz .lt. 0.25  ) then
		factor = 0.25/dz
		 dz1   = 0.25
	   elseif (dz .gt. 1.50) then
	        factor = 1.50/dz
	        dz1    = 1.50
	   endif

	amnt1=amnt*factor

! Regression is now broken up into TWO parts one for small
! one for large	water vapor ammounts.

	ireg=1
	if ( log(amnt1) .gt. h2obnd(iband) ) ireg=2

       ph2o = amnt1 *(8.314d+07 *temp )/
     &              (dz1*1.0d+05*18.01534 *1.01325d+06)

	patmx = log(patm)

!	print'(8f8.3)',log(amnt1),temp,patmx,(ph2o),amnt1,log(ph2o),dz1
!	print*, ireg,iband,aa(1:7,ireg,iband)

 	tau_log	    = ck24_3(1,ireg,iband)	      +
     $		      ck24_3(2,ireg,iband)* log(amnt1)  +
     $		      ck24_3(3,ireg,iband)* temp       +
     $                ck24_3(4,ireg,iband)* patmx       +
     $		      ck24_3(5,ireg,iband)* (ph2o)     +
     $		      ck24_3(6,ireg,iband)* amnt1       +
     $		      ck24_3(7,ireg,iband)* log(ph2o) 
!	print*,tau_log
	parm_ckd24 = exp ( tau_log )
	parm_ckd24 = parm_ckd24/factor
	return
caky26102009	end function parm_ckd24
c not f77, but f90, thanks to Reto Stockli for pointing this out
	end 

!==============================================================


      subroutine gascon_ckd(ib)
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      integer i, m, ib
      real pp, pt, pa, ph, po, dz
      real tgm
      real log_average
      integer iwtas
      common /con/tgm(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /cont_tas/iwtas
      common /thick/dz(nvx)

      integer iflb(18)
      real dp, p, t, pres, sum, dvabs
      real war, wn2, wo2, wco2, wh2o
      integer ii, iv1, iv2, idv, nptabs
      integer ib1, ie1
      real fe, fr
      double precision v1abs, v2abs, vbnds(12,2)
      real absrb(1000), wnbr(1000)
      data vbnds/ 5, 280, 400, 540, 670, 800, 980, 1100, 1250, 1400, 
     1   1700, 1900, 280, 400, 540, 670, 800, 980, 1100, 1250, 1400, 
     2   1700, 1900, 2200/ 
      data dvabs/ 1.0/ 
      if (iwtas.lt.0 .or. iwtas.gt.4) stop 'iwtas'
      if (iwtas .eq. 0) then
 
         iflb(1) = 0
         iflb(2) = 0
         iflb(3) = 0
         iflb(4) = 0
         iflb(5) = 0
         iflb(6) = 0
         iflb(7) = 0
         iflb(8) = 0
         iflb(9) = 0
         iflb(10) = 9
         iflb(11) = 8
         iflb(12) = 7
         iflb(13) = 6
         iflb(14) = 5
         iflb(15) = 4
         iflb(16) = 3
         iflb(17) = 2
         iflb(18) = 0
      endif
      if (iwtas .gt. 0) then
         iflb(1) = 0
         iflb(2) = 0
         iflb(3) = 0
         iflb(4) = 0
         iflb(5) = 0
         iflb(6) = 0
         iflb(7) = 12
         iflb(8) = 11
         iflb(9) = 10
         iflb(10) = 9
         iflb(11) = 8
         iflb(12) = 7
         iflb(13) = 6
         iflb(14) = 5
         iflb(15) = 4
         iflb(16) = 3
         iflb(17) = 2
         iflb(18) = 1
      endif
 
      do m = 1, nv
         tgm(m) = 0.0
      end do
 
      if (iflb(ib) .eq. 0) return 
 
      v1abs = vbnds(iflb(ib),1)
      v2abs = vbnds(iflb(ib),2)
      iv1 = v1abs
      iv2 = v2abs
      idv = dvabs
 
      ii = 0
      do i = iv1, iv2, idv
         ii = ii + 1
         wnbr(ii) = i
      end do
      nptabs = ii
 
      do m = 1, nv
cbm         dp = pp(m+1) - pp(m)
cbm logarithmic interpolation instead of linear
cbm ??? this should be made consistent with uvspec ???
cbm      pres = (pp(m)+pp(m+1))*0.5
         pres = log_average(pp(m), pp(m+1))
         t = (pt(m)+pt(m+1))*0.5
cbm         dz = dz_quick(pp(m),pp(m+1),pt(m),pt(m+1),ph(m),ph(m+1))
c                                        ! g/g -> g/cm2 -->molecules/cm2
cbm      wh2o = 1.02*dp*(ph(m)+ph(m+1))*0.5*(6.023E+23/18.01534)
         wh2o = pa(m)*0.1*(ph(m)+ph(m+1))*0.5*(6.023E+23/18.01534)
 
         p = pres/1013.25
         war = p*1.01325E+06/(1.3805E-16*t)*0.00934*1.0E+05*dz(m)
         wn2 = p*1.01325E+06/(1.3805E-16*t)*0.78084*1.0E+05*dz(m)
         wo2 = p*1.01325E+06/(1.3805E-16*t)*0.20948*1.0E+05*dz(m)
         wco2 = p*1.01325E+06/(1.3805E-16*t)*3.50E-04*1.0E+05*dz(m)
 
         do ii = 1, 1000
            absrb(ii) = 0.0
         end do
 
         call clough (pres, t, war, wn2, wo2, wco2, wh2o, v1abs, v2abs, 
     1      dvabs, nptabs, absrb)
 
         ii = 0
         sum = 0
         do i = iv1, iv2, idv
            ib1 = i
            ie1 = i + idv
            if (iwtas.eq.3 .or. iwtas.eq.4) call subplank (t, iv1, iv2, 
     1         idv, ib1, ie1, fe, fr)
 
            ii = ii + 1
c                                                !! log noplank
            if (iwtas .eq. 0) sum = sum + log(absrb(ii))
 
c                                                !! lin noplank
            if (iwtas .eq. 1) sum = sum + absrb(ii)
c                                                !! log noplank
            if (iwtas .eq. 2) sum = sum + log(absrb(ii))
c                                                !! lin Plank wgt
            if (iwtas .eq. 3) sum = sum + absrb(ii)*fr
c                                                !! log Plank wgt
            if (iwtas .eq. 4) sum = sum + log(absrb(ii))*fr
         end do
c                                                !! log noplank
         if (iwtas .eq. 0) tgm(m) = exp(sum/float(ii))
 
         if (iwtas .eq. 1) tgm(m) = sum/float(ii)!! lin noplank
c                                                !! log noplank
         if (iwtas .eq. 2) tgm(m) = exp(sum/float(ii))
         if (iwtas .eq. 3) tgm(m) = sum          !! lin Plank wgt
         if (iwtas .eq. 4) tgm(m) = exp(sum)     !! log Plank wgt
      end do
 
      return 
      end 
 

c  ************************************************************
c  The Following routines are for the exact CKD_2.1 solution
c  From S.A.Clough with modifications for LWONLY by Dave Kratz.
c  ************************************************************
      subroutine clough(pres, temp, wa, wn2, wo2, wco2, wh2o, v1abs, 
     1   v2abs, dvabs, nptabs, absrb)
      implicit double precision(v)
      real absrb(*)
      real wk(60)
      data xlosmt/ 2.68675E+19/ 
      data wk/ 60*0/ 
 
c  *********************************************************************
c   this program calculates the continuum optical depth
c         for an homogeneous layer
c
c   the following quantities must be specified:
c
c          pressure                   pave (mb)
c          temperature                tave ( k)
c          column amount
c            nitrogen                 wn2    (molec/cm**2)
c            oxygen                   wk(7)  (molec/cm**2)
c            carbon dioxide           wk(2)  (molec/cm**2)
c            water vapor              wk(1)  (molec/cm**2)
c          number of molecules        nmol
c          beginning wavenumber       v1abs (cm-1)
c          ending wavenumber          v2abs (cm-1)
c          sampling interval          dvabs (cm-1)
c          number of values           nptabs
c
c   the results are in array absorb
c
c   note that for an atmospheric layer:
c
c     wtot   = xlosmt * (pave/1013) * (273/tave) * (path length)
c     wbroad = the column amount for all species not explicitly provided
c     wk(m)  = (volume mixing ratio) * (column of dry air)
c
c      note that continua are included for
c
c             m=1    water vapor    air & self         0 - 20,000 cm-1
c             wn2    nitrogen       air             2020 -   2800 cm-1
c
c  The following is an example for a 1 km path (see cntnout for results)
c  *********************************************************************
      pave = pres
      tave = temp
      wk(7) = wo2
      wk(2) = wco2
      wk(1) = wh2o
      wbroad = wn2 + wa
      nmol = 7
      jrad = 1
 
      call contnm (jrad, pave, tave, wk, wbroad, nmol, v1abs, v2abs, 
     1   dvabs, nptabs, absrb)
 
      return 
      end 
 
      subroutine xint(v1a,v2a,dva,a,afact,vft,dvr3,r3,n1r3,n2r3)
c  ****************************************************************
c     this subroutine interpolates the a array stored
c     from v1a to v2a in increments of dva using a multiplicative
c     factor afact, into the r3 array from location n1r3 to n2r3 in
c     increments of dvr3
c  ****************************************************************
      implicit double precision(v)
      real a(*), r3(*)
      data onepl/ 1.001/ 
      data onemi/ 0.999/ 
      data argmin/ 34./ 
 
      recdva = 1./dva
      ilo = (v1a + dva - vft)/dvr3 + 1. + onemi
      ilo = max(ilo,n1r3)
      ihi = (v2a - dva - vft)/dvr3 + onemi
      ihi = min(ihi,n2r3)
 
      do i = ilo, ihi
         vi = vft + dvr3*float(i - 1)
         j = (vi - v1a)*recdva + onepl
         vj = v1a + dva*float(j - 1)
         p = recdva*(vi - vj)
         c = (3. - 2.*p)*p*p
         b = 0.5*p*(1. - p)
         b1 = b*(1. - p)
         b2 = b*p
         conti=(-a(j-1)*b1)+a(j)*(1.-c+b2)+a(j+1)*(c+b1)-a(j+2)*b2
         r3(i) = r3(i) + conti*afact
      end do
 
      return 
      end 
 
      real function radfn (vi, xkt)
 
c
c     function radfn calculates the radiation term for the line shape
c  *********************************************************************
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
c  *********************************************************************
      implicit double precision(v)
      data onepl/ 1.001/ 
      data onemi/ 0.999/ 
      data argmin/ 34./ 
 
c  *******************************************
c  In the small xviokt region 0.5 is required.
c  *******************************************
      xvi = vi
      if (xkt .gt. 0.0) then
         xviokt = xvi/xkt
         if (xviokt .le. 0.01) then
            radfn = 0.5*xviokt*xvi
         else if (xviokt .le. 10.0) then
            expvkt = exp((-xviokt))
            radfn = xvi*(1. - expvkt)/(1. + expvkt)
         else
            radfn = xvi
         endif
      else
         radfn = xvi
      endif
 
      return 
      end 
 
      subroutine contnm(jrad, pave, tave, wk, wbroad, nmol, v1abs, v2abs
     1   , dvabs, nptabs, absrb)
c  *********************************************
c  Subroutine contnm contains the continuum data
c  which is interpolated into the array absrb
c  *********************************************
      implicit double precision(v)
      real absrb(*)
      real wk(60)
      real csh2o(2050), cfh2o(2050)
      real c(2050), c0(2050), c1(2050), c2(2050)
      real xfac(0:50)

      data pi/ 3.1415927/ 
      data planck/ 6.626176E-27/ 
      data boltz/ 1.380662E-16/ 
      data clight/ 2.99792458E10/ 
      data avog/ 6.022045E23/ 
      data p0/ 1013./ 
      data t0/ 296./ 
      data xlosmt/ 2.68675E+19/ 
      data (xfac(i),i=0,50)/ 1.00000, 1.01792, 1.03767, 1.05749, 1.07730
     1   , 1.09708, 1.10489, 1.11268, 1.12047, 1.12822, 1.13597, 1.14367
     2   , 1.15135, 1.15904, 1.16669, 1.17431, 1.18786, 1.20134, 1.21479
     3   , 1.22821, 1.24158, 1.26580, 1.28991, 1.28295, 1.27600, 1.26896
     4   , 1.25550, 1.24213, 1.22879, 1.21560, 1.20230, 1.18162, 1.16112
     5   , 1.14063, 1.12016, 1.10195, 1.09207, 1.08622, 1.08105, 1.07765
     6   , 1.07398, 1.06620, 1.05791, 1.04905, 1.03976, 1.02981, 1.00985
     7   , 1.00000, 1.00000, 1.00000, 1.00000/ 
 
      radcn1 = 2.*planck*clight*clight*1.E-07
      radcn2 = planck*clight/boltz
 
c  ************************************
c  Assign sccs version number to module
c      hvrcnt = '3.3'
c  ************************************
      rhoave = (pave/p0)*(t0/tave)
      xkt = tave/radcn2
      wtot = wbroad
      do m = 1, nmol
         wtot = wtot + wk(m)
      end do
      wtot = wtot*1.E-20
      patm = pave/p0
 
c  ******************************
c  ********    water    ********
c  ******************************
      call sl296 (v1c, v2c, dvc, nptc, c0, v1abs, v2abs, dvabs, nptabs)
      call sl260 (v1c, v2c, dvc, nptc, c1, v1abs, v2abs, dvabs, nptabs)
      call frn296 (v1c, v2c, dvc, nptc, c2, v1abs, v2abs, dvabs, nptabs)
 
      w1 = wk(1)*1.0E-20
c  *************************************************************
c  The factor of 1.e-20 is handled this way to avoid underflows.
c  *************************************************************
      ph2o = patm*(w1/wtot)
      rh2o = ph2o*(t0/tave)
      pfrgn = patm - ph2o
      rfrgn = pfrgn*(t0/tave)
      xkt = tave/radcn2
      tfac = (tave - t0)/(260. - t0)
 
c  ******************************
c  ********     self     ********
c  ******************************
      alpha2 = 200.**2
      alphs2 = 120.**2
      betas = 5.E-06
      v0s = 1310.
      factrs = 0.15
 
c  ******************************
c  ********    foreign   ********
c  ******************************
      hwsqf = 330.**2
      betaf = 8.E-11
      v0f = 1130.
      factrf = 0.97
 
      v0f2 = 1900.
      hwsqf2 = 150.**2
      beta2 = 3.E-06
 
 
      v1h = v1c
      dvh = dvc
      npth = nptc
 
      do j = 1, nptc
         vj = v1c + dvc*float(j - 1)
         vs2 = (vj - v0s)**2
         sh2o = 0.
         if (c0(j) .gt. 0.) then
            sh2o = c0(j)*(c1(j)/c0(j))**tfac
            sfac = 1.
            if (vj.ge.700. .and. vj.le.1200.) then
               jfac = (vj - 700.)/10. + 0.00001
               sfac = xfac(jfac)
            endif
 
c     ----------------------------------------------------------------
c     correction to self continuum (1 sept 85); factor of 0.78 at 1000
c                             and  .......
c     ----------------------------------------------------------------
            sh2o = sfac*sh2o*(1. - 0.2333*(alpha2/((vj - 1050.)**2 + 
     1         alpha2)))*(1. - factrs*(alphs2/(vs2 + betas*vs2**2 + 
     2         alphs2)))
         endif
 
c     -------------------------------
c     Correction to foreign continuum
c     -------------------------------
         vf2 = (vj - v0f)**2
         vf6 = vf2*vf2*vf2
         fscal = 1. - factrf*(hwsqf/(vf2 + betaf*vf6 + hwsqf))
         vf2 = (vj - v0f2)**2
         vf4 = vf2*vf2
         fscal = fscal*(1. - 0.6*(hwsqf2/(vf2 + beta2*vf4 + hwsqf2)))
 
         c2(j) = c2(j)*fscal
         c(j) = w1*(sh2o*rh2o + c2(j)*rfrgn)
 
         csh2o(j) = 1.E-20*sh2o
         cfh2o(j) = 1.E-20*c2(j)
 
c     ---------------
c     Radiation field
c     ---------------
         if (jrad .eq. 1) c(j) = c(j)*radfn(vj,xkt)
      end do
 
      call xint (v1c, v2c, dvc, c, 1.0, v1abs, dvabs, absrb, 1, nptabs)
 
      return 
      end 
 
 
      subroutine sl296(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
c  ******************************************
c               06/28/82
c      units of (cm**3/mol) * 1.e-20
c  ******************************************
      implicit double precision(v)
      real s(305)
      real c(*)
      data v1s, v2s, dvs, npts/ -20.0, 3020.0, 10.0, 305/ 

      data (s(i),i=1,2)/ 1.1109E-01, 1.0573E-01/ 
      data (s(i),i=3,52)/ 1.0162E-01, 1.0573E-01, 1.1109E-01, 1.2574E-01
     1   , 1.3499E-01, 1.4327E-01, 1.5065E-01, 1.5164E-01, 1.5022E-01, 
     2   1.3677E-01, 1.3115E-01, 1.2253E-01, 1.1271E-01, 1.0070E-01, 
     3   8.7495E-02, 8.0118E-02, 6.9940E-02, 6.2034E-02, 5.6051E-02, 
     4   4.7663E-02, 4.2450E-02, 3.6690E-02, 3.3441E-02, 3.0711E-02, 
     5   2.5205E-02, 2.2113E-02, 1.8880E-02, 1.6653E-02, 1.4626E-02, 
     6   1.2065E-02, 1.0709E-02, 9.1783E-03, 7.7274E-03, 6.7302E-03, 
     7   5.6164E-03, 4.9089E-03, 4.1497E-03, 3.5823E-03, 3.1124E-03, 
     8   2.6414E-03, 2.3167E-03, 2.0156E-03, 1.7829E-03, 1.5666E-03, 
     9   1.3928E-03, 1.2338E-03, 1.0932E-03, 9.7939E-04, 8.8241E-04, 
     .   7.9173E-04/ 
      data (s(i),i=53,102)/ 7.1296E-04, 6.4179E-04, 5.8031E-04, 
     1   5.2647E-04, 4.7762E-04, 4.3349E-04, 3.9355E-04, 3.5887E-04, 
     2   3.2723E-04, 2.9919E-04, 2.7363E-04, 2.5013E-04, 2.2876E-04, 
     3   2.0924E-04, 1.9193E-04, 1.7618E-04, 1.6188E-04, 1.4891E-04, 
     4   1.3717E-04, 1.2647E-04, 1.1671E-04, 1.0786E-04, 9.9785E-05, 
     5   9.2350E-05, 8.5539E-05, 7.9377E-05, 7.3781E-05, 6.8677E-05, 
     6   6.3993E-05, 5.9705E-05, 5.5788E-05, 5.2196E-05, 4.8899E-05, 
     7   4.5865E-05, 4.3079E-05, 4.0526E-05, 3.8182E-05, 3.6025E-05, 
     8   3.4038E-05, 3.2203E-05, 3.0511E-05, 2.8949E-05, 2.7505E-05, 
     9   2.6170E-05, 2.4933E-05, 2.3786E-05, 2.2722E-05, 2.1736E-05, 
     .   2.0819E-05, 1.9968E-05/ 
      data (s(i),i=103,152)/ 1.9178E-05, 1.8442E-05, 1.7760E-05, 
     1   1.7127E-05, 1.6541E-05, 1.5997E-05, 1.5495E-05, 1.5034E-05, 
     2   1.4614E-05, 1.4230E-05, 1.3883E-05, 1.3578E-05, 1.3304E-05, 
     3   1.3069E-05, 1.2876E-05, 1.2732E-05, 1.2626E-05, 1.2556E-05, 
     4   1.2544E-05, 1.2604E-05, 1.2719E-05, 1.2883E-05, 1.3164E-05, 
     5   1.3581E-05, 1.4187E-05, 1.4866E-05, 1.5669E-05, 1.6717E-05, 
     6   1.8148E-05, 2.0268E-05, 2.2456E-05, 2.5582E-05, 2.9183E-05, 
     7   3.3612E-05, 3.9996E-05, 4.6829E-05, 5.5055E-05, 6.5897E-05, 
     8   7.5360E-05, 8.7213E-05, 1.0046E-04, 1.1496E-04, 1.2943E-04, 
     9   1.5049E-04, 1.6973E-04, 1.8711E-04, 2.0286E-04, 2.2823E-04, 
     .   2.6780E-04, 2.8766E-04/ 
      data (s(i),i=153,202)/ 3.1164E-04, 3.3640E-04, 3.6884E-04, 
     1   3.9159E-04, 3.8712E-04, 3.7433E-04, 3.4503E-04, 3.1003E-04, 
     2   2.8027E-04, 2.5253E-04, 2.3408E-04, 2.2836E-04, 2.4442E-04, 
     3   2.7521E-04, 2.9048E-04, 3.0489E-04, 3.2646E-04, 3.3880E-04, 
     4   3.3492E-04, 3.0987E-04, 2.9482E-04, 2.8711E-04, 2.6068E-04, 
     5   2.2683E-04, 1.9996E-04, 1.7788E-04, 1.6101E-04, 1.3911E-04, 
     6   1.2013E-04, 1.0544E-04, 9.4224E-05, 8.1256E-05, 7.3667E-05, 
     7   6.2233E-05, 5.5906E-05, 5.1619E-05, 4.5140E-05, 4.0273E-05, 
     8   3.3268E-05, 3.0258E-05, 2.6440E-05, 2.3103E-05, 2.0749E-05, 
     9   1.8258E-05, 1.6459E-05, 1.4097E-05, 1.2052E-05, 1.0759E-05, 
     .   9.1400E-06, 8.1432E-06/ 
      data (s(i),i=203,252)/ 7.1460E-06, 6.4006E-06, 5.6995E-06, 
     1   4.9372E-06, 4.4455E-06, 3.9033E-06, 3.4740E-06, 3.1269E-06, 
     2   2.8059E-06, 2.5558E-06, 2.2919E-06, 2.0846E-06, 1.8983E-06, 
     3   1.7329E-06, 1.5929E-06, 1.4631E-06, 1.3513E-06, 1.2461E-06, 
     4   1.1519E-06, 1.0682E-06, 9.9256E-07, 9.2505E-07, 8.6367E-07, 
     5   8.0857E-07, 7.5674E-07, 7.0934E-07, 6.6580E-07, 6.2580E-07, 
     6   5.8853E-07, 5.5333E-07, 5.2143E-07, 4.9169E-07, 4.6431E-07, 
     7   4.3898E-07, 4.1564E-07, 3.9405E-07, 3.7403E-07, 3.5544E-07, 
     8   3.3819E-07, 3.2212E-07, 3.0714E-07, 2.9313E-07, 2.8003E-07, 
     9   2.6777E-07, 2.5628E-07, 2.4551E-07, 2.3540E-07, 2.2591E-07, 
     .   2.1701E-07, 2.0866E-07/ 
      data (s(i),i=253,302)/ 2.0082E-07, 1.9349E-07, 1.8665E-07, 
     1   1.8027E-07, 1.7439E-07, 1.6894E-07, 1.6400E-07, 1.5953E-07, 
     2   1.5557E-07, 1.5195E-07, 1.4888E-07, 1.4603E-07, 1.4337E-07, 
     3   1.4093E-07, 1.3828E-07, 1.3569E-07, 1.3270E-07, 1.2984E-07, 
     4   1.2714E-07, 1.2541E-07, 1.2399E-07, 1.2102E-07, 1.1878E-07, 
     5   1.1728E-07, 1.1644E-07, 1.1491E-07, 1.1305E-07, 1.1235E-07, 
     6   1.1228E-07, 1.1224E-07, 1.1191E-07, 1.1151E-07, 1.1098E-07, 
     7   1.1068E-07, 1.1109E-07, 1.1213E-07, 1.1431E-07, 1.1826E-07, 
     8   1.2322E-07, 1.3025E-07, 1.4066E-07, 1.5657E-07, 1.7214E-07, 
     9   1.9449E-07, 2.2662E-07, 2.6953E-07, 3.1723E-07, 3.7028E-07, 
     .   4.4482E-07, 5.3852E-07/ 
      data (s(i),i=303,305)/ 6.2639E-07, 7.2175E-07, 7.7626E-07/ 
 
      dvc = dvs
      v1c = v1abs - dvc
      v2c = v2abs + dvc
 
      i1 = (v1c - v1s)/dvs
      if (v1c .lt. v1s) i1 = i1 - 1
 
      v1c = v1s + dvs*float(i1)
      i2 = (v2c - v1s)/dvs
      nptc = i2 - i1 + 3
      v2c = v1c + dvs*float(nptc - 1)
      do j = 1, nptc
         i = i1 + j
         c(j) = 0.
         if (i.ge.1 .and. i.le.305) c(j) = s(i)
      end do
 
      return 
      end 
 
      subroutine sl260(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
c  **************************************
c               06/28/82
c      units of (cm**3/mol) * 1.e-20
c  **************************************
      implicit double precision(v)
      real s(305)
      real c(*)
      data v1s, v2s, dvs, npts/ -20.0, 3020.0, 10.0, 305/ 
      data (s(i),i=1,2)/ 1.7750E-01, 1.7045E-01/ 
      data (s(i),i=3,52)/ 1.6457E-01, 1.7045E-01, 1.7750E-01, 2.0036E-01
     1   , 2.1347E-01, 2.2454E-01, 2.3428E-01, 2.3399E-01, 2.3022E-01, 
     2   2.0724E-01, 1.9712E-01, 1.8317E-01, 1.6724E-01, 1.4780E-01, 
     3   1.2757E-01, 1.1626E-01, 1.0098E-01, 8.9033E-02, 7.9770E-02, 
     4   6.7416E-02, 5.9588E-02, 5.1117E-02, 4.6218E-02, 4.2179E-02, 
     5   3.4372E-02, 2.9863E-02, 2.5252E-02, 2.2075E-02, 1.9209E-02, 
     6   1.5816E-02, 1.3932E-02, 1.1943E-02, 1.0079E-02, 8.7667E-03, 
     7   7.4094E-03, 6.4967E-03, 5.5711E-03, 4.8444E-03, 4.2552E-03, 
     8   3.6953E-03, 3.2824E-03, 2.9124E-03, 2.6102E-03, 2.3370E-03, 
     9   2.1100E-03, 1.9008E-03, 1.7145E-03, 1.5573E-03, 1.4206E-03, 
     .   1.2931E-03/ 
      data (s(i),i=53,102)/ 1.1803E-03, 1.0774E-03, 9.8616E-04, 
     1   9.0496E-04, 8.3071E-04, 7.6319E-04, 7.0149E-04, 6.4637E-04, 
     2   5.9566E-04, 5.4987E-04, 5.0768E-04, 4.6880E-04, 4.3317E-04, 
     3   4.0037E-04, 3.7064E-04, 3.4325E-04, 3.1809E-04, 2.9501E-04, 
     4   2.7382E-04, 2.5430E-04, 2.3630E-04, 2.1977E-04, 2.0452E-04, 
     5   1.9042E-04, 1.7740E-04, 1.6544E-04, 1.5442E-04, 1.4425E-04, 
     6   1.3486E-04, 1.2618E-04, 1.1817E-04, 1.1076E-04, 1.0391E-04, 
     7   9.7563E-05, 9.1696E-05, 8.6272E-05, 8.1253E-05, 7.6607E-05, 
     8   7.2302E-05, 6.8311E-05, 6.4613E-05, 6.1183E-05, 5.8001E-05, 
     9   5.5048E-05, 5.2307E-05, 4.9761E-05, 4.7395E-05, 4.5197E-05, 
     .   4.3155E-05, 4.1256E-05/ 
      data (s(i),i=103,152)/ 3.9491E-05, 3.7849E-05, 3.6324E-05, 
     1   3.4908E-05, 3.3594E-05, 3.2374E-05, 3.1244E-05, 3.0201E-05, 
     2   2.9240E-05, 2.8356E-05, 2.7547E-05, 2.6814E-05, 2.6147E-05, 
     3   2.5551E-05, 2.5029E-05, 2.4582E-05, 2.4203E-05, 2.3891E-05, 
     4   2.3663E-05, 2.3531E-05, 2.3483E-05, 2.3516E-05, 2.3694E-05, 
     5   2.4032E-05, 2.4579E-05, 2.5234E-05, 2.6032E-05, 2.7119E-05, 
     6   2.8631E-05, 3.0848E-05, 3.3262E-05, 3.6635E-05, 4.0732E-05, 
     7   4.5923E-05, 5.3373E-05, 6.1875E-05, 7.2031E-05, 8.5980E-05, 
     8   9.8642E-05, 1.1469E-04, 1.3327E-04, 1.5390E-04, 1.7513E-04, 
     9   2.0665E-04, 2.3609E-04, 2.6220E-04, 2.8677E-04, 3.2590E-04, 
     .   3.8624E-04, 4.1570E-04/ 
      data (s(i),i=153,202)/ 4.5207E-04, 4.9336E-04, 5.4500E-04, 
     1   5.8258E-04, 5.8086E-04, 5.6977E-04, 5.3085E-04, 4.8020E-04, 
     2   4.3915E-04, 4.0343E-04, 3.7853E-04, 3.7025E-04, 3.9637E-04, 
     3   4.4675E-04, 4.7072E-04, 4.9022E-04, 5.2076E-04, 5.3676E-04, 
     4   5.2755E-04, 4.8244E-04, 4.5473E-04, 4.3952E-04, 3.9614E-04, 
     5   3.4086E-04, 2.9733E-04, 2.6367E-04, 2.3767E-04, 2.0427E-04, 
     6   1.7595E-04, 1.5493E-04, 1.3851E-04, 1.1874E-04, 1.0735E-04, 
     7   9.0490E-05, 8.1149E-05, 7.4788E-05, 6.5438E-05, 5.8248E-05, 
     8   4.8076E-05, 4.3488E-05, 3.7856E-05, 3.3034E-05, 2.9592E-05, 
     9   2.6088E-05, 2.3497E-05, 2.0279E-05, 1.7526E-05, 1.5714E-05, 
     .   1.3553E-05, 1.2145E-05/ 
      data (s(i),i=203,252)/ 1.0802E-05, 9.7681E-06, 8.8196E-06, 
     1   7.8291E-06, 7.1335E-06, 6.4234E-06, 5.8391E-06, 5.3532E-06, 
     2   4.9079E-06, 4.5378E-06, 4.1716E-06, 3.8649E-06, 3.5893E-06, 
     3   3.3406E-06, 3.1199E-06, 2.9172E-06, 2.7348E-06, 2.5644E-06, 
     4   2.4086E-06, 2.2664E-06, 2.1359E-06, 2.0159E-06, 1.9051E-06, 
     5   1.8031E-06, 1.7074E-06, 1.6185E-06, 1.5356E-06, 1.4584E-06, 
     6   1.3861E-06, 1.3179E-06, 1.2545E-06, 1.1951E-06, 1.1395E-06, 
     7   1.0873E-06, 1.0384E-06, 9.9250E-07, 9.4935E-07, 9.0873E-07, 
     8   8.7050E-07, 8.3446E-07, 8.0046E-07, 7.6834E-07, 7.3800E-07, 
     9   7.0931E-07, 6.8217E-07, 6.5648E-07, 6.3214E-07, 6.0909E-07, 
     .   5.8725E-07, 5.6655E-07/ 
      data (s(i),i=253,302)/ 5.4693E-07, 5.2835E-07, 5.1077E-07, 
     1   4.9416E-07, 4.7853E-07, 4.6381E-07, 4.5007E-07, 4.3728E-07, 
     2   4.2550E-07, 4.1450E-07, 4.0459E-07, 3.9532E-07, 3.8662E-07, 
     3   3.7855E-07, 3.7041E-07, 3.6254E-07, 3.5420E-07, 3.4617E-07, 
     4   3.3838E-07, 3.3212E-07, 3.2655E-07, 3.1865E-07, 3.1203E-07, 
     5   3.0670E-07, 3.0252E-07, 2.9749E-07, 2.9184E-07, 2.8795E-07, 
     6   2.8501E-07, 2.8202E-07, 2.7856E-07, 2.7509E-07, 2.7152E-07, 
     7   2.6844E-07, 2.6642E-07, 2.6548E-07, 2.6617E-07, 2.6916E-07, 
     8   2.7372E-07, 2.8094E-07, 2.9236E-07, 3.1035E-07, 3.2854E-07, 
     9   3.5481E-07, 3.9377E-07, 4.4692E-07, 5.0761E-07, 5.7715E-07, 
     .   6.7725E-07, 8.0668E-07/ 
      data (s(i),i=303,305)/ 9.3716E-07, 1.0797E-06, 1.1689E-06/ 
 
      dvc = dvs
      v1c = v1abs - dvc
      v2c = v2abs + dvc
 
      i1 = (v1c - v1s)/dvs
      if (v1c .lt. v1s) i1 = i1 - 1
 
      v1c = v1s + dvs*float(i1)
      i2 = (v2c - v1s)/dvs
      nptc = i2 - i1 + 3
      v2c = v1c + dvs*float(nptc - 1)
      do j = 1, nptc
         i = i1 + j
         c(j) = 0.
         if (i.ge.1 .and. i.le.305) c(j) = s(i)
      end do
 
      return 
      end 
 
      subroutine frn296(v1c,v2c,dvc,nptc,c,v1abs,v2abs,dvabs,nptabs)
c  **********************************
c               06/28/82
c      units of (cm**3/mol)*1.e-20
c  **********************************
      implicit double precision(v)
      real s(305)
      real c(*)
      data v1s, v2s, dvs, npts/ -20.0, 3020.0, 10.0, 305/ 
      data (s(i),i=1,2)/ 1.2859E-02, 1.1715E-02/ 
      data (s(i),i=3,52)/ 1.1038E-02, 1.1715E-02, 1.2859E-02, 1.5326E-02
     1   , 1.6999E-02, 1.8321E-02, 1.9402E-02, 1.9570E-02, 1.9432E-02, 
     2   1.7572E-02, 1.6760E-02, 1.5480E-02, 1.3984E-02, 1.2266E-02, 
     3   1.0467E-02, 9.4526E-03, 8.0485E-03, 6.9484E-03, 6.1416E-03, 
     4   5.0941E-03, 4.4836E-03, 3.8133E-03, 3.4608E-03, 3.1487E-03, 
     5   2.4555E-03, 2.0977E-03, 1.7266E-03, 1.4920E-03, 1.2709E-03, 
     6   9.8081E-04, 8.5063E-04, 6.8822E-04, 5.3809E-04, 4.4679E-04, 
     7   3.3774E-04, 2.7979E-04, 2.1047E-04, 1.6511E-04, 1.2993E-04, 
     8   9.3033E-05, 7.4360E-05, 5.6428E-05, 4.5442E-05, 3.4575E-05, 
     9   2.7903E-05, 2.1374E-05, 1.6075E-05, 1.3022E-05, 1.0962E-05, 
     .   8.5959E-06/ 
      data (s(i),i=53,102)/ 6.9125E-06, 5.3808E-06, 4.3586E-06, 
     1   3.6394E-06, 2.9552E-06, 2.3547E-06, 1.8463E-06, 1.6036E-06, 
     2   1.3483E-06, 1.1968E-06, 1.0333E-06, 8.4484E-07, 6.7195E-07, 
     3   5.0947E-07, 4.2343E-07, 3.4453E-07, 2.7830E-07, 2.3063E-07, 
     4   1.9951E-07, 1.7087E-07, 1.4393E-07, 1.2575E-07, 1.0750E-07, 
     5   8.2325E-08, 5.7524E-08, 4.4482E-08, 3.8106E-08, 3.4315E-08, 
     6   2.9422E-08, 2.5069E-08, 2.2402E-08, 1.9349E-08, 1.6152E-08, 
     7   1.2208E-08, 8.9660E-09, 7.1322E-09, 6.1028E-09, 5.2938E-09, 
     8   4.5350E-09, 3.4977E-09, 2.9511E-09, 2.4734E-09, 2.0508E-09, 
     9   1.8507E-09, 1.6373E-09, 1.5171E-09, 1.3071E-09, 1.2462E-09, 
     .   1.2148E-09, 1.2590E-09/ 
      data (s(i),i=103,152)/ 1.3153E-09, 1.3301E-09, 1.4483E-09, 
     1   1.6944E-09, 2.0559E-09, 2.2954E-09, 2.6221E-09, 3.2606E-09, 
     2   4.2392E-09, 5.2171E-09, 6.2553E-09, 8.2548E-09, 9.5842E-09, 
     3   1.1280E-08, 1.3628E-08, 1.7635E-08, 2.1576E-08, 2.4835E-08, 
     4   3.0014E-08, 3.8485E-08, 4.7440E-08, 5.5202E-08, 7.0897E-08, 
     5   9.6578E-08, 1.3976E-07, 1.8391E-07, 2.3207E-07, 2.9960E-07, 
     6   4.0408E-07, 5.9260E-07, 7.8487E-07, 1.0947E-06, 1.4676E-06, 
     7   1.9325E-06, 2.6587E-06, 3.4534E-06, 4.4376E-06, 5.8061E-06, 
     8   7.0141E-06, 8.4937E-06, 1.0186E-05, 1.2034E-05, 1.3837E-05, 
     9   1.6595E-05, 1.9259E-05, 2.1620E-05, 2.3681E-05, 2.7064E-05, 
     .   3.2510E-05, 3.5460E-05/ 
      data (s(i),i=153,202)/ 3.9109E-05, 4.2891E-05, 4.7757E-05, 
     1   5.0981E-05, 5.0527E-05, 4.8618E-05, 4.4001E-05, 3.7982E-05, 
     2   3.2667E-05, 2.7794E-05, 2.4910E-05, 2.4375E-05, 2.7316E-05, 
     3   3.2579E-05, 3.5499E-05, 3.8010E-05, 4.1353E-05, 4.3323E-05, 
     4   4.3004E-05, 3.9790E-05, 3.7718E-05, 3.6360E-05, 3.2386E-05, 
     5   2.7409E-05, 2.3626E-05, 2.0631E-05, 1.8371E-05, 1.5445E-05, 
     6   1.2989E-05, 1.1098E-05, 9.6552E-06, 8.0649E-06, 7.2365E-06, 
     7   5.9137E-06, 5.2759E-06, 4.8860E-06, 4.1321E-06, 3.5918E-06, 
     8   2.7640E-06, 2.4892E-06, 2.1018E-06, 1.7848E-06, 1.5855E-06, 
     9   1.3569E-06, 1.1986E-06, 9.4693E-07, 7.4097E-07, 6.3443E-07, 
     .   4.8131E-07, 4.0942E-07/ 
      data (s(i),i=203,252)/ 3.3316E-07, 2.8488E-07, 2.3461E-07, 
     1   1.7397E-07, 1.4684E-07, 1.0953E-07, 8.5396E-08, 6.9261E-08, 
     2   5.4001E-08, 4.5430E-08, 3.2791E-08, 2.5995E-08, 2.0225E-08, 
     3   1.5710E-08, 1.3027E-08, 1.0229E-08, 8.5277E-09, 6.5249E-09, 
     4   5.0117E-09, 3.9906E-09, 3.2332E-09, 2.7847E-09, 2.4570E-09, 
     5   2.3359E-09, 2.0599E-09, 1.8436E-09, 1.6559E-09, 1.4910E-09, 
     6   1.2794E-09, 9.8229E-10, 8.0054E-10, 6.0769E-10, 4.5646E-10, 
     7   3.3111E-10, 2.4428E-10, 1.8007E-10, 1.3291E-10, 9.7974E-11, 
     8   7.8271E-11, 6.3833E-11, 5.4425E-11, 4.6471E-11, 4.0209E-11, 
     9   3.5227E-11, 3.1212E-11, 2.8840E-11, 2.7762E-11, 2.7935E-11, 
     .   3.2012E-11, 3.9525E-11/ 
      data (s(i),i=253,302)/ 5.0303E-11, 6.8027E-11, 9.3954E-11, 
     1   1.2986E-10, 1.8478E-10, 2.5331E-10, 3.4827E-10, 4.6968E-10, 
     2   6.2380E-10, 7.9106E-10, 1.0026E-09, 1.2102E-09, 1.4146E-09, 
     3   1.6154E-09, 1.7510E-09, 1.8575E-09, 1.8742E-09, 1.8700E-09, 
     4   1.8582E-09, 1.9657E-09, 2.1204E-09, 2.0381E-09, 2.0122E-09, 
     5   2.0436E-09, 2.1213E-09, 2.0742E-09, 1.9870E-09, 2.0465E-09, 
     6   2.1556E-09, 2.2222E-09, 2.1977E-09, 2.1047E-09, 1.9334E-09, 
     7   1.7357E-09, 1.5754E-09, 1.4398E-09, 1.4018E-09, 1.5459E-09, 
     8   1.7576E-09, 2.1645E-09, 2.9480E-09, 4.4439E-09, 5.8341E-09, 
     9   8.0757E-09, 1.1658E-08, 1.6793E-08, 2.2694E-08, 2.9468E-08, 
     .   3.9278E-08, 5.2145E-08/ 
      data (s(i),i=303,305)/ 6.4378E-08, 7.7947E-08, 8.5321E-08/ 
 
      dvc = dvs
      v1c = v1abs - dvc
      v2c = v2abs + dvc
 
      i1 = (v1c - v1s)/dvs
      if (v1c .lt. v1s) i1 = i1 - 1
 
      v1c = v1s + dvs*float(i1)
      i2 = (v2c - v1s)/dvs
      nptc = i2 - i1 + 3
      v2c = v1c + dvs*float(nptc - 1)
 
      do j = 1, nptc
         i = i1 + j
         c(j) = 0.
         if (i.ge.1 .and. i.le.305) c(j) = s(i)
      end do
 
      return 
      end 
c---------- 4/2/97 (3) -- PREVIOUS 801 LINES
c==================================================
      subroutine subplank(temp, ib, ie, id, ib1, ie1, fe, fr)
      integer id_fine
      save id_fine
      integer ib0
      save ib0
      integer ie0
      save ie0
      integer id0
      save id0
      integer nnn
      save nnn
      real fluxb
      save fluxb
      real temp0
      save temp0
      data id_fine/ 1/ 
 
      if(ib0.ne.ib.or.ie0.ne.ie.or.id0.ne.id.or.temp.ne.temp0)then
         print *, 'initalize subplank'
         ib0 = ib
         ie0 = ie
         id0 = id
         temp0 = temp
         call plank (temp, ib, ie, id_fine, fluxb)
         nnn = (ie - ib)/float(id)
      endif
 
      call plank (temp, ib1, ie1, id_fine, flux)
      fr = flux/fluxb
      fe = fr*nnn
!     print*,nnn,flux,fluxb,fe
 
      return 
      end 
c==================================================
      subroutine plank(t, ib, ie, id, fluxo)
      double precision c1, c2, u, flux, el, h, c, rk
      real fluxo
      c = 2.99793E+10
      h = 6.62620E-27
      rk = 1.38062E-16
      c1 = 3.74E-16
      c2 = 1.44E-02
 
      c1 = 3.741951E-16
      c2 = 1.439E-02
 
      rw = 1.0E-06
 
      flux = 0
      do ic = ib, ie - id, id
         um1 = 1.0E+04/float(ic)
         um2 = 1.0E+04/float(ic + id)
         dum = um1 - um2
 
         rcm = ic + id/2.0                       !center
 
         um = 1.0E+04/rcm             !!float(ic) ! WAVENUMBER to Micron
 
         upw = um/rcm                            !!float(ic)
 
         u = um*1.0E-06
         el = c1/(u**5*(exp(c2/(u*t)) - 1.0))
 
c     rw=rwm*1.0e-06
         rw = dum*1.0E-6
 
         el2 = el*rw
 
         flux = flux + el2
         umo = um
      end do
 
      fluxo = flux
      return 
      end 
c Fu 07-08-98
      block data ice1new
      INCLUDE 'fl_radparams.inc'
      common /ic1new/ap(3,mbx), bps(4,mbsx), bpir(4,mbirx), cp(4,mbx), 
     1   dps(4,mbsx)
      data ap/ -2.9172062E-05, 2.5192544E+00, 0.0, -2.2948980E-05, 
     1   2.5212550E+00, 0.0, -2.9772840E-04, 2.5400320E+00, 0.0, 
     2   4.2668223E-04, 2.4933372E+00, 0.0, 4.3226531E-04, 2.4642946E+00
     3   , 0.0, 9.5918990E-05, 2.5232218E+00, 0.0, -2.308881E-03, 
     4   2.814002E+00, 1.072211E+00, -2.465236E-03, 2.833187E+00, 
     5   -4.227573E-01, -3.034573E-03, 2.900043E+00, -1.849911E+00, 
     6   -4.936610E-03, 3.087764E+00, -3.884262E+00, -8.178608E-03, 
     7   3.401245E+00, -8.812820E+00, -8.372696E-03, 3.455018E+00, 
     8   -1.516692E+01, -1.691632E-03, 2.765756E+00, -8.331033E+00, 
     9   -4.159424E-03, 3.047325E+00, -5.061568E+00, -9.524174E-03, 
     .   3.587742E+00, -1.068895E+01, -1.334860E-02, 4.043808E+00, 
     1   -2.171029E+01, 3.325756E-03, 2.601360E+00, -1.909602E+01, 
     2   4.919685E-03, 2.327741E+00,  - 1.390858E+01/ 
      data bps/ 1.3540265E-07, 9.9282217E-08, -7.3843168E-11, 
     1   3.3111862E-13, -2.1458450E-06, 2.1984010E-05, -4.4225520E-09, 
     2   1.0711940E-11, 1.4027890E-04, 1.3919010E-03, -5.1005610E-06, 
     3   1.4032930E-08, 5.7801650E-03, 2.4420420E-03, -1.1985030E-05, 
     4   3.3878720E-08, 2.7122737E-01, 1.9809794E-03, -1.5071269E-05, 
     5   5.0103900E-08, 1.6215025E-01, 6.3734393E-03, -5.7740959E-05, 
     6   1.9109300E-07/ 
      data bpir/ 4.346482E-01, 1.721457E-02, -1.623227E-04, 5.561523E-07
     1   , 7.428957E-01, 1.279601E-02, -1.391803E-04, 5.180104E-07, 
     2   8.862434E-01, 1.226538E-02, -1.523076E-04, 6.000892E-07, 
     3   7.152274E-01, 1.621734E-02, -1.868544E-04, 7.078738E-07, 
     4   5.874323E-01, 1.876628E-02, -2.045834E-04, 7.510080E-07, 
     5   5.409536E-01, 1.949649E-02, -2.050908E-04, 7.364680E-07, 
     6   1.195515E+00, 3.350616E-03, -5.266996E-05, 2.233377E-07, 
     7   1.466481E+00, -2.129226E-03, -1.361630E-05, 1.193649E-07, 
     8   9.551440E-01, 1.309792E-02, -1.793694E-04, 7.313392E-07, 
     9   3.003701E-01, 2.051529E-02, -1.931684E-04, 6.583031E-07, 
     .   2.005578E-01, 2.132614E-02, -1.751052E-04, 5.355885E-07, 
     1   8.869787E-01, 2.118409E-02, -2.781429E-04, 1.094562E-06/ 
      data cp/ 7.4812728E-01, 9.5684492E-04, -1.1151708E-06, 
     1   -8.1557303E-09, 7.5212480E-01, 1.1045100E-03, -2.9157100E-06, 
     2   -1.3429900E-09, 7.5320460E-01, 1.8845180E-03, -9.7571460E-06, 
     3   2.2428270E-08, 7.7381780E-01, 2.2260760E-03, -1.4052790E-05, 
     4   3.7896870E-08, 8.7020490E-01, 1.6645530E-03, -1.4886030E-05, 
     5   4.9867270E-08, 7.4212060E-01, 5.2621900E-03, -5.0877550E-05, 
     6   1.7307870E-07, 7.962716E-01, 3.003488E-03, -2.082376E-05, 
     7   5.366545E-08, 8.472918E-01, 2.559953E-03, -2.182660E-05, 
     8   6.879977E-08, 8.741665E-01, 2.455409E-03, -2.456935E-05, 
     9   8.641223E-08, 8.522816E-01, 2.523627E-03, -2.149196E-05, 
     .   6.685067E-08, 8.609604E-01, 2.200445E-03, -1.748105E-05, 
     1   5.176616E-08, 8.906280E-01, 1.903269E-03, -1.733552E-05, 
     2   5.855071E-08, 8.663385E-01, 2.797934E-03, -3.187011E-05, 
     3   1.217209E-07, 7.984021E-01, 3.977117E-03, -4.471984E-05, 
     4   1.694919E-07, 7.363466E-01, 4.798266E-03, -4.513292E-05, 
     5   1.525774E-07, 7.260484E-01, 2.664334E-03, -1.251136E-05, 
     6   2.243377E-08, 6.891414E-01, 6.192281E-03, -6.459514E-05, 
     7   2.436963E-07, 4.949276E-01, 1.186174E-02, -1.267629E-04, 
     8   4.603574E-07/ 
      data dps/ 1.1572963E-01, 2.5648064E-04, 1.9131293E-06, 
     1   -1.2460341E-08, 1.1360752E-01, 2.4156171E-04, 2.0185942E-06, 
     2   -1.2876106E-08, 1.1241170E-01, -1.7635186E-07, 2.1499248E-06, 
     3   -1.2949304E-08, 1.0855775E-01, -3.2496217E-04, 3.4207304E-06, 
     4   -1.6247759E-08, 5.7783360E-02, -4.1158260E-04, 4.2361240E-06, 
     5   -1.7204950E-08, 1.1367129E-01, -1.9711061E-03, 1.6078010E-05, 
     6    - 5.1736898E-08/ 
c **********************************************************************
      end 
 
      subroutine icenew(ib, unscaled)
      INCLUDE 'fl_radparams.inc'
      integer unscaled
      common /clouds/pre(nvx), plwc(nvx), pdge(nvx), piwc(nvx)
      common /ic1new/ap(3,mbx), bps(4,mbsx), bpir(4,mbirx), cp(4,mbx), 
     1   dps(4,mbsx)
      common /thick/dz(nvx)
      common /ic/ti(nvx), wi(nvx), wwi(nvx,4)
      do i = 1, nv
         if (piwc(i) .lt. 1.0E-5) then
            ti(i) = 0.0
            wi(i) = 0.0
            wwi(i,1) = 0.0
            wwi(i,2) = 0.0
            wwi(i,3) = 0.0
            wwi(i,4) = 0.0
         else
            fw1 = pdge(i)
            fw2 = fw1*pdge(i)
            fw3 = fw2*pdge(i)
            if (ib .le. mbs) then
               tau = dz(i)*1000.0*piwc(i)*(ap(1,ib)+ap(2,ib)/fw1)
               omega = 1.0 - (bps(1,ib)+bps(2,ib)*fw1+bps(3,ib)*fw2+bps(
     1            4,ib)*fw3)
               asy=cp(1,ib)+cp(2,ib)*fw1+cp(3,ib)*fw2+cp(4,ib)*fw3
               fd=dps(1,ib)+dps(2,ib)*fw1+dps(3,ib)*fw2+dps(4,ib)*fw3
               f = 0.5/omega + fd
               fw = f*omega
               if (unscaled .eq. 0) then
                  ti(i) = (1.0 - fw)*tau
                  wi(i) = (1.0 - f)*omega/(1.0 - fw)
                  gg = (asy - f)/(1.0 - f)
               else
                  ti(i) = tau
                  wi(i) = omega
                  gg    = asy
               endif
               x1 = gg
               x2 = x1*gg
               x3 = x2*gg
               x4 = x3*gg
               wwi(i,1) = 3.0*x1
               wwi(i,2) = 5.0*x2
               wwi(i,3) = 7.0*x3
               wwi(i,4) = 9.0*x4
            else
               ibr = ib - mbs
               betae = piwc(i)*(ap(1,ib)+ap(2,ib)/fw1+ap(3,ib)/fw2)
               betaa = piwc(i)/fw1*(bpir(1,ibr)+bpir(2,ibr)*fw1+bpir(3,
     1            ibr)*fw2+bpir(4,ibr)*fw3)
               asy=cp(1,ib)+cp(2,ib)*fw1+cp(3,ib)*fw2+cp(4,ib)*fw3
               ti(i) = dz(i)*1000.0*betae
               wi(i) = 1.0 - betaa/betae
               gg = asy
               x1 = gg
               x2 = x1*gg
               x3 = x2*gg
               x4 = x3*gg
               wwi(i,1) = 3.0*x1
               wwi(i,2) = 5.0*x2
               wwi(i,3) = 7.0*x3
               wwi(i,4) = 9.0*x4
            endif
         endif
      end do
      return 
      end 
c Fu 07-08-98
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      block data
      real a, u
      common /dis/a(4)
      common /point/u(4)
      data a/ 0.5, 0.5, 0.5, 0.5/ 
      data u/ -0.7886752, -0.2113247, 0.2113247, 0.7886752/ 
      end 
 
      BLOCK DATA ICE1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      REAL AP, BP, CPS, DPS, CPIR
      COMMON /IC1/AP(3,MBX), BP(4,MBX), CPS(4,4,MBSX), DPS(4,MBSX), CPIR
     1   (4,MBIRX)
      DATA AP/ -6.656E-3, 3.686, 0.00, -6.656E-3, 3.686, 0.00, -6.656E-3
     1   , 3.686, 0.00, -6.656E-3, 3.686, 0.00, -6.656E-3, 3.686, 0.00, 
     2   -6.656E-3, 3.686, 0.00, -7.770E-3, 3.734, 11.85, -8.088E-3, 
     3   3.717, 17.17, -8.441E-3, 3.715, 19.48, -9.061E-3, 3.741, 26.48
     4   , -9.609E-3, 3.768, 34.11, -1.153E-2, 4.109, 17.32, -8.294E-3, 
     5   3.925, 1.315, -1.026E-2, 4.105, 16.36, -1.151E-2, 4.182, 31.13
     6   , -1.704E-2, 4.830, 16.27, -1.741E-2, 5.541, -58.42, -7.752E-3
     7   , 4.624,  - 42.01/ 
      DATA BP/ .10998E-05, -.26101E-07, .10896E-08, -.47387E-11, 
     1   .20208E-04, .96483E-05, .83009E-07, -.32217E-09, .13590E-03, 
     2   .73453E-03, .28281E-05, -.18272E-07, -.16598E-02, .20933E-02, 
     3   -.13977E-05, -.18703E-07, .46180E+00, .24471E-03, -.27839E-05, 
     4   .10379E-07, .42362E-01, .86425E-02, -.75519E-04, .24056E-06, 
     5   .19960E+00, .37800E-02, -.14910E-04, .00000E+00, .30140E+00, 
     6   .26390E-02, -.11160E-04, .00000E+00, .39080E+00, .12720E-02, 
     7   -.55640E-05, .00000E+00, .31050E+00, .26030E-02, -.11390E-04, 
     8   .00000E+00, .20370E+00, .42470E-02, -.18100E-04, .00000E+00, 
     9   .23070E+00, .38300E-02, -.16160E-04, .00000E+00, .56310E+00, 
     .   -.14340E-02, .62980E-05, .00000E+00, .52070E+00, -.97780E-03, 
     1   .37250E-05, .00000E+00, .32540E+00, .34340E-02, -.30810E-04, 
     2   .91430E-07, .10280E+00, .50190E-02, -.20240E-04, .00000E+00, 
     3   .39640E+00, -.31550E-02, .64170E-04, -.29790E-06, .80790E+00, 
     4   -.70040E-02, .52090E-04,  - .14250E-06/ 
      DATA CPS/ .22110E+01, -.10398E-02, .65199E-04, -.34498E-06, 
     1   .32201E+01, .94227E-03, .80947E-04, -.47428E-06, .41610E+01, 
     2   .74396E-03, .82690E-04, -.45251E-06, .51379E+01, .51545E-02, 
     3   .11881E-04, -.15556E-06, .22151E+01, -.77982E-03, .63750E-04, 
     4   -.34466E-06, .31727E+01, .15597E-02, .82021E-04, -.49665E-06, 
     5   .40672E+01, .25800E-02, .71550E-04, -.43051E-06, .49882E+01, 
     6   .86489E-02, -.18318E-04, -.59275E-07, .22376E+01, .10293E-02, 
     7   .50842E-04, -.30135E-06, .31549E+01, .47115E-02, .70684E-04, 
     8   -.47622E-06, .39917E+01, .82830E-02, .53927E-04, -.41778E-06, 
     9   .48496E+01, .15998E-01, -.39320E-04, -.43862E-07, .23012E+01, 
     .   .33854E-02, .23528E-04, -.20068E-06, .31730E+01, .93439E-02, 
     1   .36367E-04, -.38390E-06, .39298E+01, .16424E-01, .10502E-04, 
     2   -.35086E-06, .47226E+01, .25872E-01, -.77542E-04, -.21999E-07, 
     3   .27975E+01, .29741E-02, -.32344E-04, .11636E-06, .43532E+01, 
     4   .11234E-01, -.12081E-03, .43435E-06, .56835E+01, .24681E-01, 
     5   -.26480E-03, .95314E-06, .68271E+01, .42788E-01, -.45615E-03, 
     6   .16368E-05, .19655E+01, .20094E-01, -.17067E-03, .50806E-06, 
     7   .28803E+01, .36091E-01, -.28365E-03, .79656E-06, .34613E+01, 
     8   .58525E-01, -.46455E-03, .13444E-05, .39568E+01, .81480E-01, 
     9   -.64777E-03, .19022E-05/ 
      DATA DPS/ .12495E+00, -.43582E-03, .14092E-04, -.69565E-07, 
     1   .12363E+00, -.44419E-03, .14038E-04, -.68851E-07, .12117E+00, 
     2   -.48474E-03, .12495E-04, -.62411E-07, .11581E+00, -.55031E-03, 
     3   .98776E-05, -.50193E-07, -.15968E-03, .10115E-04, -.12472E-06, 
     4   .48667E-09, .13830E+00, -.18921E-02, .12030E-04,  - .31698E-07
     5   / 
      DATA CPIR/ .79550, 2.524E-3, -1.022E-5, 0.000E+0, .86010, 1.599E-3
     1   , -6.465E-6, 0.000E+0, .89150, 1.060E-3, -4.171E-6, 0.000E+0, 
     2   .87650, 1.198E-3, -4.485E-6, 0.000E+0, .88150, 9.858E-4, 
     3   -3.116E-6, 0.000E+0, .91670, 5.499E-4, -1.507E-6, 0.000E+0, 
     4   .90920, 9.295E-4, -3.877E-6, 0.000E+0, .84540, 1.429E-3, 
     5   -5.859E-6, 0.000E+0, .76780, 2.571E-3, -1.041E-5, 0.000E+0, 
     6   .72900, 2.132E-3, -5.584E-6, 0.000E+0, .70240, 4.581E-3, 
     7   -3.054E-5, 6.684E-8, .22920, 1.724E-2, -1.573E-4, 4.995E-7/ 
      END 
 
      block data water1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      real re, fl, bz, wz, gz
 
      common/wat1/re(ncx),fl(ncx),bz(ncx,mbx),wz(ncx,mbx),gz(ncx,mbx)
      data re/ 4.18, 5.36, 5.89, 6.16, 9.27, 9.84, 12.10, 31.23/ 
      data fl/ 0.05, 0.14, 0.22, 0.28, 0.50, 0.47, 1.00, 2.50/ 
      data bz/ 15.11, 40.25, 59.81, 72.43, 83.69, 73.99, 128.17, 120.91
     1   , 15.74, 41.70, 61.52, 74.47, 85.78, 75.59, 130.46, 121.84, 
     2   16.38, 43.52, 64.84, 77.97, 87.31, 77.36, 134.30, 124.06, 17.57
     3   , 45.78, 66.44, 80.15, 90.49, 79.90, 137.56, 125.92, 18.19, 
     4   46.63, 69.39, 82.20, 91.46, 79.99, 138.21, 126.08, 21.30, 51.88
     5   , 77.77, 87.02, 94.91, 83.55, 143.46, 128.45, 22.44, 57.35, 
     6   84.41, 103.50, 103.49, 84.17, 152.77, 132.07, 18.32, 52.69, 
     7   76.67, 100.31, 105.46, 92.86, 157.82, 133.03, 17.27, 50.44, 
     8   74.18, 96.76, 105.32, 95.25, 158.07, 134.48, 13.73, 44.90, 
     9   67.70, 90.85, 109.16, 105.48, 163.11, 136.21, 10.30, 36.28, 
     .   57.23, 76.43, 106.45, 104.90, 161.73, 136.62, 7.16, 26.40, 
     1   43.51, 57.24, 92.55, 90.55, 149.10, 135.13, 6.39, 21.00, 33.81
     2   , 43.36, 66.90, 63.58, 113.83, 125.65, 10.33, 30.87, 47.63, 
     3   60.33, 79.54, 73.92, 127.46, 128.21, 11.86, 35.64, 54.81, 69.85
     4   , 90.39, 84.16, 142.49, 135.25, 10.27, 33.08, 51.81, 67.26, 
     5   93.24, 88.60, 148.71, 140.42, 6.72, 24.09, 39.42, 51.68, 83.34
     6   , 80.72, 140.14, 143.57, 3.92, 14.76, 25.32, 32.63, 60.85, 
     7   58.81, 112.30, 145.62/ 
      data wz/ .999999, .999999, .999999, .999999, .999998, .999999, 
     1   .999998, .999997, .999753, .999700, .999667, .999646, .999492, 
     2   .999470, .999344, .998667, .995914, .994967, .994379, .993842, 
     3   .991385, .990753, .988908, .974831, .983761, .978981, .976568, 
     4   .974700, .963466, .959934, .953865, .897690, .702949, .683241, 
     5   .679723, .669045, .642616, .632996, .629776, .588820, .947343, 
     6   .929619, .924806, .914557, .877169, .867047, .853661, .737426, 
     7   .919356, .896274, .885924, .881097, .812772, .781637, .775418, 
     8   .637341, .874717, .861122, .847850, .851677, .787171, .772952, 
     9   .753143, .618656, .764750, .752410, .736529, .743435, .671272, 
     .   .659392, .639492, .549941, .807536, .808700, .795994, .805489, 
     1   .750577, .755524, .709472, .571989, .753346, .772026, .767273, 
     2   .777079, .751264, .760973, .712536, .568286, .632722, .676332, 
     3   .684631, .693552, .707986, .717724, .682430, .552867, .288885, 
     4   .348489, .371653, .380367, .454540, .465769, .475409, .493881, 
     5   .261827, .306283, .321340, .333051, .392917, .406876, .417450, 
     6   .484593, .295804, .339929, .352494, .365502, .416229, .430369, 
     7   .435267, .491356, .301214, .354746, .369346, .381906, .433602, 
     8   .447397, .447406, .486968, .243714, .318761, .344642, .352770, 
     9   .427906, .438979, .445972, .477264, .109012, .187230, .226849, 
     .   .224976, .331382, .335917, .374882, .457067/ 
      data gz/ .838, .839, .844, .847, .849, .860, .853, .859, .809, 
     1   .810, .819, .823, .823, .849, .833, .843, .774, .787, .781, 
     2   .792, .812, .836, .815, .833, .801, .802, .793, .793, .814, 
     3   .829, .818, .832, .877, .873, .879, .880, .885, .899, .891, 
     4   .908, .783, .769, .777, .756, .764, .776, .770, .797, .818, 
     5   .805, .824, .830, .815, .801, .820, .845, .810, .802, .826, 
     6   .840, .829, .853, .840, .868, .774, .766, .799, .818, .815, 
     7   .869, .834, .869, .734, .728, .767, .797, .796, .871, .818, 
     8   .854, .693, .688, .736, .772, .780, .880, .808, .846, .643, 
     9   .646, .698, .741, .759, .882, .793, .839, .564, .582, .637, 
     .   .690, .719, .871, .764, .819, .466, .494, .546, .609, .651, 
     1   .823, .701, .766, .375, .410, .455, .525, .583, .773, .637, 
     2   .710, .262, .301, .334, .406, .485, .695, .545, .631, .144, 
     3   .181, .200, .256, .352, .562, .413, .517, .060, .077, .088, 
     4   .112, .181, .310, .222, .327/ 
      end 
 
      block data rayle1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      real ri, rix, riy
 
      common /ray1/ri(mbsx), rix(10), riy(10)
      data ri/ 0.9022E-5, 0.5282E-6, 0.5722E-7, 0.1433E-7, 0.4526E-8, 
     1   0.1529E-8/ 
      data rix/ 1.775E-04, 1.222E-04, 6.606E-05, 4.588E-05, 3.524E-05, 
     1   2.420E-05, 1.244E-05, 6.510E-06, 3.475E-06, 1.792E-06/ 
      data riy/ 1.761E-04, 1.214E-04, 6.512E-05, 4.583E-05, 3.510E-05, 
     1   2.401E-05, 1.209E-05, 6.510E-06, 3.475E-06, 1.792E-06/ 
 
      end 
 
      block data rain1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      real rwc, brn, wrnf, grn
      common /rai1/rwc, brn(mbx), wrnf(mbx), grn(mbx)
      data rwc/ 0.5/ 
      data brn/ 1.5377, 1.5377, 1.5379, 1.5385, 1.5396, 1.5417, 1.5454, 
     1   1.5478, 1.5512, 1.5559, 1.5600, 1.5642, 1.5647, 1.5741, 1.5862
     2   , 1.5993, 1.6149, 1.6765/ 
      data wrnf/ .999932, .97096, .74627, .56719, .53023, .53815, .53233
     1   , .52884, .53192, .52969, .52716, .52321, .51904, .53859, 
     2   .55169, .55488, .55334, .55218/ 
      data grn/ .88323, .89067, .92835, .96626, .97553, .96626, .97226, 
     1   .97663, .97216, .97467, .97745, .98156, .98584, .96374, .94218
     2   , .93266, .92990, .90729/ 
      end 
 
      block data graup1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      real gwc, bg, wgf, gg
      common /gra1/gwc, bg(mbx), wgf(mbx), gg(mbx)
      data gwc/ 0.5/ 
      data bg/ 0.83939, 0.83940, 0.83940, 0.83941, 0.83946, 0.83951, 
     1   0.83967, 0.83979, 0.83995, 0.84029, 0.84058, 0.84097, 0.84143, 
     2   0.84286, 0.84418, 0.84825, 0.85421, 0.87477/ 
      data wgf/ 0.999911, 0.97115, 0.56192, 0.53156, 0.52579, 0.53846, 
     1   0.53296, 0.53017, 0.53182, 0.53180, 0.52959, 0.52446, 0.52342, 
     2   0.54914, 0.55258, 0.54307, 0.53160, 0.55474/ 
      data gg/ 0.89218, 0.89940, 0.96820, 0.97816, 0.98141, 0.96373, 
     1   0.97173, 0.97559, 0.97330, 0.97327, 0.97626, 0.98274, 0.98396, 
     2   0.94673, 0.94213, 0.95539, 0.97097, 0.93183/ 
      end 
 
      block data legend
      real p0d, p1d, p2d, p3d
      common /legen/p0d(4), p1d(4), p2d(4), p3d(4)
      data p0d/ .100000E+01, .100000E+01, .100000E+01, .100000E+01/ 
      data p1d/ -.788675E+00, -.211325E+00, .211325E+00, .788675E+00/ 
      data p2d/ .433013E+00, -.433013E+00, -.433013E+00, .433013E+00/ 
      data p3d/ -.433940E-01, .293394E+00, -.293394E+00, .433940E-01/ 
      end 
 
      BLOCK DATA LEGENF
      REAL P11D, P22D, P33D
      COMMON /LEGEN1/P11D(4,4), P22D(4,4), P33D(4,4)
      DATA P11D/ .311004E+00, .833334E-01, -.833334E-01, -.311004E+00, 
     1   .833334E-01, .223291E-01, -.223291E-01, -.833334E-01, 
     2   -.833334E-01, -.223291E-01, .223291E-01, .833334E-01, 
     3   -.311004E+00, -.833334E-01, .833334E-01, .311004E+00/ 
      DATA P22D/ .937501E-01, -.937501E-01, -.937501E-01, .937501E-01, 
     1   -.937501E-01, .937501E-01, .937501E-01, -.937501E-01, 
     2   -.937501E-01, .937501E-01, .937501E-01, -.937501E-01, 
     3   .937501E-01, -.937501E-01, -.937501E-01, .937501E-01/ 
      DATA P33D/ .941520E-03, -.636577E-02, .636577E-02, -.941520E-03, 
     1   -.636577E-02, .430400E-01, -.430400E-01, .636577E-02, 
     2   .636577E-02, -.430400E-01, .430400E-01, -.636577E-02, 
     3   -.941520E-03, .636577E-02, -.636577E-02, .941520E-03/ 
      END 
 
c---------- 4/1/97 (6) -- Replaces old ckd1 block data.
      subroutine ckd1_init(isolar_spectrum, lray)
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 to
c one. fko3 is the corresponding ozone absorption coefficient in units
c of (cm-atm)**-1 (Fu, 1991). The spectral region is from 50000 cm**-1
c to 14500 cm**-1.
c *********************************************************************
c      common /band1/ hk(10), fko3(10)
c      data hk / .24, .16, .24, .28, .03,
c     1            .016, .01, .008, .008, .008 /
c      data fko3 / .2204e-08,.1207e-01,.4537e-01,.1032e+00,.1740e+00,
c     1           .1210e+01,.7367e+01,.2050e+02,.8100e+02,.2410e+03 /
c *********************************************************************
c The treatment of O3 absorption here is modified to consider  the
c dependence of aerosol single-scattering properties on wavelength
c in the first band.  We divide this band into 10 subintervals  as
c 175.439 - 224.719 - 243.902 - 285.714 - 298.507 - 322.500 -
c 357.500 - 437.500 - 497.500 - 595.000 - 692.500 following Grant
c (personal communication, 1996).
c *********************************************************************
c---------- 4/1/97 (6)
      common /band1/hk(10), fko3(10), sol_spect(0:7), fkh2o(10)
 
      real fko3_1(10), fko3_2(10), fko3_3(10), fko3_4(10)
      logical lray
      integer j1
      fko3_1(1) = .2751E+02
      fko3_1(2) = .1564E+03
      fko3_1(3) = .1746E+03
      fko3_1(4) = .2684E+02
      fko3_1(5) = .1962E+01
      fko3_1(6) = .7402E-01
      fko3_1(7) = .7295E-03
      fko3_1(8) = .1239E-01
      fko3_1(9) = .8011E-01
      fko3_1(10) = .8011E-01
 
! ! RECOMPUTED FROM 111 SUBBAND FU CODE.
! ! ENERGY WEIGHTED - Linear Average
! !    fko3_2 = ( /3.87E+01,1.66E+02,1.94E+02,2.93E+01,3.22E+00,
! !     1            8.38E-02,7.38E-04,1.25E-02,8.45E-02,7.88E-02/)
! ! ENERGY WEIGHTED - LOG Average
! !   fko3_3 =( /3.28E+01,1.56E+02,1.75E+02,2.70E+01,2.19E+00,
! !   &     2.98E-02,1.94E-04,1.02E-02,7.76E-02,7.19E-02/)
! ! ENERGY WEIGHTED - LOG(1-6) Lin(7-10) Average
! !   data fko3_4 /3.28E+01,1.56E+02,1.75E+02,2.70E+01,2.19E+00,
! !     &      2.98E-02,7.38E-04,1.25E-02,8.45E-02,7.88E-02/
 
! ! ENERGY WEIGHTED - LOG(1-5) Lin(6-10) Average
      fko3_4(1) = 3.28E+01
      fko3_4(2) = 1.56E+02
      fko3_4(3) = 1.75E+02
      fko3_4(4) = 2.70E+01
      fko3_4(5) = 2.19E+00
      fko3_4(6) = 8.38E-02
      fko3_4(7) = 7.38E-04
      fko3_4(8) = 1.25E-02
      fko3_4(9) = 8.45E-02
      fko3_4(10) = 7.88E-02
 
! CHOU
!       data fkh2o / 6*0.0,4*0.00075 /
! Q.Fu
      fkh2o(1) = 0.
      fkh2o(2) = 0.
      fkh2o(3) = 0.
      fkh2o(4) = 0.
      fkh2o(5) = 0.
      fkh2o(6) = 0.
      fkh2o(7) = 0.
      fkh2o(8) = 7.63E-05
      fkh2o(9) = 1.90E-03
      fkh2o(10) = 2.57E-03
 
! SOLAR SPECTRA Choice::
 
      if (lray) then
         do j1 = 1, 10
            fko3(j1) = fko3_4(j1)                !!! LOG & Lin wgt
         end do
      else
         do j1 = 1, 10
            fko3(j1) = fko3_1(j1)                !! OLD
         end do
      endif
!--------------------------------------------------------------------
      if (isolar_spectrum .eq. 1) then           !!!  Mod3.5_kurucz
         sol_spect(1-1) = 1373.013
         sol_spect(2-1) = 626.929
         sol_spect(3-1) = 493.654
         sol_spect(4-1) = 155.256
         sol_spect(5-1) = 51.183
         sol_spect(6-1) = 28.516
         sol_spect(7-1) = 5.550
         sol_spect(8-1) = 11.924
         hk(1) = 0.00117
         hk(2) = 0.00145
         hk(3) = 0.01060
         hk(4) = 0.01060
         hk(5) = 0.02618
         hk(6) = 0.05375
         hk(7) = 0.17448
         hk(8) = 0.18855
         hk(9) = 0.28872
         hk(10) = 0.24451
! --------------------------------------------------------
      else if (isolar_spectrum .eq. 2) then      !!!  Mod3.7_cebchkur
         sol_spect(1-1) = 1362.119
         sol_spect(2-1) = 626.734
         sol_spect(3-1) = 484.953
         sol_spect(4-1) = 153.697
         sol_spect(5-1) = 51.119
         sol_spect(6-1) = 28.420
         sol_spect(7-1) = 5.542
         sol_spect(8-1) = 11.910
         hk(1) = 0.00125
         hk(2) = 0.00154
         hk(3) = 0.01040
         hk(4) = 0.01009
         hk(5) = 0.02481
         hk(6) = 0.05374
         hk(7) = 0.17567
         hk(8) = 0.18909
         hk(9) = 0.28927
         hk(10) = 0.24413
! --------------------------------------------------------
      else if (isolar_spectrum .eq. 3) then      !!!  Mod3.7_chkur
         sol_spect(1-1) = 1359.751
         sol_spect(2-1) = 624.364
         sol_spect(3-1) = 484.953
         sol_spect(4-1) = 153.697
         sol_spect(5-1) = 51.119
         sol_spect(6-1) = 28.420
         sol_spect(7-1) = 5.542
         sol_spect(8-1) = 11.910
         hk(1) = 0.00128
         hk(2) = 0.00158
         hk(3) = 0.01108
         hk(4) = 0.01120
         hk(5) = 0.02336
         hk(6) = 0.05156
         hk(7) = 0.17472
         hk(8) = 0.18980
         hk(9) = 0.29037
         hk(10) = 0.24506
! --------------------------------------------------------
      else if (isolar_spectrum .eq. 4) then      !!!  Mod3.7_Newkur
         sol_spect(1-1) = 1367.997
         sol_spect(2-1) = 630.440
         sol_spect(3-1) = 487.123
         sol_spect(4-1) = 153.697
         sol_spect(5-1) = 51.119
         sol_spect(6-1) = 28.420
         sol_spect(7-1) = 5.542
         sol_spect(8-1) = 11.910
         hk(1) = 0.00127
         hk(2) = 0.00142
         hk(3) = 0.01142
         hk(4) = 0.01087
         hk(5) = 0.02440
         hk(6) = 0.05532
         hk(7) = 0.17325
         hk(8) = 0.19032
         hk(9) = 0.28889
         hk(10) = 0.24283
! --------------------------------------------------------
      else if (isolar_spectrum .eq. 5) then      !!!  Mod3.7_thkur
         sol_spect(1-1) = 1376.231
         sol_spect(2-1) = 634.524
         sol_spect(3-1) = 491.276
         sol_spect(4-1) = 153.697
         sol_spect(5-1) = 51.119
         sol_spect(6-1) = 28.420
         sol_spect(7-1) = 5.542
         sol_spect(8-1) = 11.910
         hk(1) = 0.00125
         hk(2) = 0.00152
         hk(3) = 0.01005
         hk(4) = 0.00970
         hk(5) = 0.02423
         hk(6) = 0.05515
         hk(7) = 0.17862
         hk(8) = 0.19247
         hk(9) = 0.28592
         hk(10) = 0.24111
! --------------------------------------------------------
      else if (isolar_spectrum .eq. 6) then      !!!
!SPECIAL FU_111 subband integration of vis-subband
         sol_spect(1-1) = 1351.0
         sol_spect(2-1) = 619.618
         sol_spect(3-1) = 484.295
         sol_spect(4-1) = 149.845
         sol_spect(5-1) = 48.7302
         sol_spect(6-1) = 31.6576
         sol_spect(7-1) = 5.79927
         sol_spect(8-1) = 11.0
         hk(1) = 0.001113
         hk(2) = 0.001452
         hk(3) = 0.010599
         hk(4) = 0.010586
         hk(5) = 0.026237
         hk(6) = 0.053743
         hk(7) = 0.174724
         hk(8) = 0.188548
         hk(9) = 0.295757
         hk(10) = 0.237242
 
 
 
      else
!ORIGINAL FU THEKEKERA
         sol_spect(1-1) = 1351.0
         sol_spect(2-1) = 619.618
         sol_spect(3-1) = 484.295
         sol_spect(4-1) = 149.845
         sol_spect(5-1) = 48.7302
         sol_spect(6-1) = 31.6576
         sol_spect(7-1) = 5.79927
         sol_spect(8-1) = 11.0
         hk(1) = .0012221
         hk(2) = .0014711
         hk(3) = .0100210
         hk(4) = .0097280
         hk(5) = .0240560
         hk(6) = .0528740
         hk(7) = .1766475
         hk(8) = .1897695
         hk(9) = .2870026
         hk(10) = .2472082
      endif
!     print*,'SOLSPEC',isolar_spectrum,sol_spect(0:7)
 
      end 
c==================================================
cbm calculate the average assuming exponential variation
      real function log_average (d1, d2)
      
      real d1, d2
      real smaller, larger, diff

      smaller = d1
      larger  = d2
      if (d2 .LE. d1) then
         smaller = d2
         larger  = d1
      endif
      
      diff = abs(d2-d1)

c linear interpolation to avoid instabilities      
      if (smaller .LE. 0 .OR. diff .LE. 0.001*larger) then
         log_average = 0.5 * (d1 + d2)
      else 
         log_average = (d2-d1) / log(d2/d1)
      endif

      return 
      end
c==================================================
      BLOCK DATA CKD2
      COMMON /BAND2/HK(8), COEH2O(3,11,8)
      DATA HK/ .71, .11, .06, .06, .04, .016, .0034, .0006/ 
      DATA (((COEH2O(K,J,I),I=1,8),J=1,11),K=1,3)/ -.1735E+02, 
     1   -.1407E+02, -.1268E+02, -.1131E+02, -.9261E+01, -.6666E+01, 
     2   -.3937E+01, -.5448E+00, -.1690E+02, -.1365E+02, -.1232E+02, 
     3   -.1101E+02, -.9058E+01, -.6574E+01, -.3914E+01, -.5529E+00, 
     4   -.1643E+02, -.1323E+02, -.1195E+02, -.1068E+02, -.8840E+01, 
     5   -.6475E+01, -.3889E+01, -.6143E+00, -.1598E+02, -.1282E+02, 
     6   -.1157E+02, -.1035E+02, -.8598E+01, -.6339E+01, -.3848E+01, 
     7   -.6636E+00, -.1551E+02, -.1241E+02, -.1119E+02, -.1001E+02, 
     8   -.8342E+01, -.6178E+01, -.3788E+01, -.8181E+00, -.1506E+02, 
     9   -.1201E+02, -.1082E+02, -.9692E+01, -.8073E+01, -.6017E+01, 
     .   -.3703E+01, -.9003E+00, -.1446E+02, -.1154E+02, -.1042E+02, 
     1   -.9332E+01, -.7810E+01, -.5846E+01, -.3576E+01, -.1083E+01, 
     2   -.1394E+02, -.1112E+02, -.1005E+02, -.8992E+01, -.7548E+01, 
     3   -.5674E+01, -.3477E+01, -.1266E+01, -.1351E+02, -.1076E+02, 
     4   -.9722E+01, -.8702E+01, -.7334E+01, -.5531E+01, -.3401E+01, 
     5   -.1524E+01, -.1311E+02, -.1044E+02, -.9422E+01, -.8423E+01, 
     6   -.7117E+01, -.5383E+01, -.3410E+01, -.1785E+01, -.1274E+02, 
     7   -.1015E+02, -.9162E+01, -.8190E+01, -.6949E+01, -.5236E+01, 
     8   -.3477E+01, -.2082E+01, .2407E-02, .2847E-02, .3768E-02, 
     9   .4626E-02, .5631E-02, .4542E-02, .3475E-02, -.3085E-02, 
     .   .2428E-02, .2805E-02, .3412E-02, .3893E-02, .4773E-02, 
     1   .3998E-02, .2742E-02, -.2556E-02, .2428E-02, .2721E-02, 
     2   .3077E-02, .3161E-02, .4019E-02, .3224E-02, .2512E-02, 
     3   -.1884E-02, .2449E-02, .2617E-02, .2763E-02, .2658E-02, 
     4   .3286E-02, .2617E-02, .1989E-02, -.1740E-02, .2512E-02, 
     5   .2470E-02, .2470E-02, .2282E-02, .2512E-02, .1926E-02, 
     6   .1465E-02, -.2612E-02, .2554E-02, .2303E-02, .2303E-02, 
     7   .1842E-02, .2030E-02, .1340E-02, .1068E-02, -.1413E-02, 
     8   .2449E-02, .2198E-02, .2030E-02, .1465E-02, .1528E-02, 
     9   .9838E-03, .1005E-02, -.1099E-02, .2868E-02, .2198E-02, 
     .   .1968E-02, .1382E-02, .1172E-02, .5652E-03, .6070E-03, 
     1   -.1662E-02, .3077E-02, .2219E-02, .1800E-02, .1277E-02, 
     2   .1005E-02, .3349E-03, .2512E-03, -.1195E-02, .3182E-02, 
     3   .2219E-02, .1758E-02, .1172E-02, .7326E-03, .4815E-03, 
     4   .6280E-04, -.1880E-02, .3265E-02, .2114E-02, .1696E-02, 
     5   .1298E-02, .4187E-03, .4187E-03, -.3768E-03, -.1467E-02, 
     6   -.1180E-04, -.1294E-04, -.1142E-04, -.7232E-05, -.8754E-05, 
     7   -.1484E-04, -.8373E-05, .1028E-04, -.1218E-04, -.1142E-04, 
     8   -.9515E-05, -.1522E-05, -.9134E-05, -.1484E-04, -.3425E-05, 
     9   .1142E-06, -.1294E-04, -.9895E-05, -.7231E-05, -.4187E-05, 
     .   -.7612E-05, -.3806E-05, .1522E-05, -.3882E-05, -.1256E-04, 
     1   -.8754E-05, -.7612E-05, -.6470E-05, -.4948E-05, -.3425E-05, 
     2   .4948E-05, -.1054E-04, -.1370E-04, -.6089E-05, -.8373E-05, 
     3   -.5709E-05, -.3045E-05, -.3806E-05, .5328E-05, .8678E-05, 
     4   -.1370E-04, -.6851E-05, -.8373E-05, -.1522E-05, -.3425E-05, 
     5   .0000E+00, .1256E-04, -.1572E-04, -.1484E-04, -.7231E-05, 
     6   -.7992E-05, -.4567E-05, -.2664E-05, -.3807E-06, -.1522E-05, 
     7   .2169E-05, -.1713E-04, -.9515E-05, -.6089E-05, -.6851E-05, 
     8   -.3045E-05, -.1142E-05, .1903E-05, .9363E-05, -.1560E-04, 
     9   -.9134E-05, -.5328E-05, -.4948E-05, .0000E+00, .7611E-06, 
     .   -.6851E-05, .1252E-04, -.1522E-04, -.8373E-05, -.6089E-05, 
     1   -.6089E-05, -.3805E-06, -.1142E-05, -.3807E-06, .2512E-05, 
     2   -.1599E-04, -.7231E-05, -.5709E-05, -.4567E-05, .1522E-05, 
     3   -.2284E-05, -.3941E-10, .5290E-05/ 
      END 
 
      BLOCK DATA CKD3
      COMMON /BAND3/HK(12), COEH2O(3,11,12)
      DATA HK/.34,.11,.1,.09,.12,.1,.06,.04,.026,.01,.0035,.0005/
      DATA (((COEH2O(K,J,I),I=1,12),J=1,11),K=1,3)/ -.1900E+02, 
     1   -.1515E+02, -.1344E+02, -.1224E+02, -.1081E+02, -.9337E+01, 
     2   -.7965E+01, -.6585E+01, -.4578E+01, -.2247E+01, .1747E+00, 
     3   .3083E+01, -.1854E+02, -.1471E+02, -.1300E+02, -.1181E+02, 
     4   -.1039E+02, -.8927E+01, -.7576E+01, -.6238E+01, -.4317E+01, 
     5   -.2119E+01, .1888E+00, .3033E+01, -.1808E+02, -.1426E+02, 
     6   -.1257E+02, -.1137E+02, -.9966E+01, -.8513E+01, -.7177E+01, 
     7   -.5885E+01, -.4053E+01, -.1977E+01, .2245E+00, .3005E+01, 
     8   -.1763E+02, -.1381E+02, -.1213E+02, -.1094E+02, -.9542E+01, 
     9   -.8094E+01, -.6779E+01, -.5524E+01, -.3788E+01, -.1796E+01, 
     .   .2961E+00, .2828E+01, -.1716E+02, -.1337E+02, -.1170E+02, 
     1   -.1051E+02, -.9116E+01, -.7677E+01, -.6381E+01, -.5153E+01, 
     2   -.3493E+01, -.1607E+01, .3850E+00, .2660E+01, -.1670E+02, 
     3   -.1295E+02, -.1127E+02, -.1008E+02, -.8690E+01, -.7265E+01, 
     4   -.5991E+01, -.4799E+01, -.3212E+01, -.1438E+01, .4582E+00, 
     5   .2588E+01, -.1596E+02, -.1231E+02, -.1067E+02, -.9501E+01, 
     6   -.8151E+01, -.6793E+01, -.5588E+01, -.4458E+01, -.2940E+01, 
     7   -.1257E+01, .4888E+00, .2260E+01, -.1530E+02, -.1184E+02, 
     8   -.1017E+02, -.8992E+01, -.7661E+01, -.6369E+01, -.5213E+01, 
     9   -.4145E+01, -.2701E+01, -.1108E+01, .4239E+00, .1974E+01, 
     .   -.1481E+02, -.1144E+02, -.9756E+01, -.8573E+01, -.7255E+01, 
     1   -.5994E+01, -.4868E+01, -.3829E+01, -.2485E+01, -.9738E+00, 
     2   .3343E+00, .1667E+01, -.1439E+02, -.1108E+02, -.9360E+01, 
     3   -.8183E+01, -.6885E+01, -.5646E+01, -.4559E+01, -.3555E+01, 
     4   -.2314E+01, -.8904E+00, .2169E+00, .1289E+01, -.1402E+02, 
     5   -.1073E+02, -.8987E+01, -.7817E+01, -.6551E+01, -.5335E+01, 
     6   -.4278E+01, -.3316E+01, -.2147E+01, -.8695E+00, .1587E-01, 
     7   .8658E+00, .1132E-01, .8855E-02, .6698E-02, .5296E-02, 
     8   .4396E-02, .3370E-02, .3245E-02, .4145E-02, .4731E-02, 
     9   .4756E-02, .3116E-02, -.2763E-02, .1135E-01, .8917E-02, 
     .   .6657E-02, .5170E-02, .4207E-02, .3056E-02, .2868E-02, 
     1   .3433E-02, .3726E-02, .4109E-02, .2836E-02, -.3119E-02, 
     2   .1135E-01, .8980E-02, .6615E-02, .5045E-02, .4061E-02, 
     3   .2847E-02, .2491E-02, .2847E-02, .2910E-02, .2671E-02, 
     4   .2396E-02, -.3245E-02, .1135E-01, .9043E-02, .6594E-02, 
     5   .4940E-02, .3914E-02, .2638E-02, .2156E-02, .2261E-02, 
     6   .2051E-02, .1978E-02, .1566E-02, -.3203E-02, .1139E-01, 
     7   .9085E-02, .6531E-02, .4835E-02, .3768E-02, .2428E-02, 
     8   .1842E-02, .1612E-02, .1591E-02, .1279E-02, .7201E-03, 
     9   -.2763E-02, .1143E-01, .9085E-02, .6447E-02, .4752E-02, 
     .   .3684E-02, .2261E-02, .1570E-02, .1235E-02, .1151E-02, 
     1   .7243E-03, .6489E-04, -.2240E-02, .1135E-01, .9001E-02, 
     2   .5694E-02, .4438E-02, .3412E-02, .1968E-02, .1235E-02, 
     3   .9420E-03, .8792E-03, .5045E-03, -.1821E-03, -.1936E-02, 
     4   .1174E-01, .9273E-02, .5882E-02, .4689E-02, .3454E-02, 
     5   .1947E-02, .1151E-02, .6070E-03, .6698E-03, .9420E-04, 
     6   -.6740E-03, -.2707E-02, .1218E-01, .9336E-02, .6050E-02, 
     7   .4731E-02, .3475E-02, .1863E-02, .1151E-02, .4605E-03, 
     8   .3768E-03, -.1214E-03, -.4396E-03, -.1903E-02, .1235E-01, 
     9   .9294E-02, .6029E-02, .4584E-02, .3370E-02, .1800E-02, 
     .   .1068E-02, .2303E-03, .1675E-03, -.4501E-03, -.7571E-03, 
     1   -.1149E-02, .1233E-01, .9315E-02, .6029E-02, .4438E-02, 
     2   .3203E-02, .1842E-02, .9629E-03, .0000E+00, -.2198E-03, 
     3   -.5338E-03, -.9721E-03, -.7661E-03, -.3692E-04, -.3844E-04, 
     4   -.2588E-04, -.1180E-04, -.1066E-04, -.3426E-05, -.2664E-05, 
     5   .7611E-06, .6089E-05, -.4568E-06, -.2077E-04, -.1142E-04, 
     6   -.3730E-04, -.3806E-04, -.2360E-04, -.1256E-04, -.1180E-04, 
     7   -.4567E-05, -.3425E-05, -.2284E-05, -.1522E-05, -.4225E-05, 
     8   -.9940E-05, -.4187E-05, -.3501E-04, -.3844E-04, -.2131E-04, 
     9   -.1256E-04, -.9896E-05, -.3806E-05, -.4186E-05, .7612E-06, 
     .   -.1903E-05, .4110E-05, .1789E-05, -.2169E-04, -.3425E-04, 
     1   -.3882E-04, -.1941E-04, -.1294E-04, -.9515E-05, -.4567E-05, 
     2   -.4186E-05, .1522E-05, -.4187E-10, .4605E-05, -.2588E-05, 
     3   .6470E-05, -.3501E-04, -.3730E-04, -.1751E-04, -.1332E-04, 
     4   -.1066E-04, -.3806E-05, -.4567E-05, -.1142E-05, -.3045E-05, 
     5   .1104E-05, -.1058E-04, .2816E-04, -.3578E-04, -.3501E-04, 
     6   -.1751E-04, -.1332E-04, -.1218E-04, -.3806E-05, -.3425E-05, 
     7   -.3806E-06, -.4187E-05, -.6090E-06, -.6965E-05, -.3463E-04, 
     8   -.3578E-04, -.3349E-04, -.1675E-04, -.9895E-05, -.9515E-05, 
     9   -.6090E-05, -.6470E-05, -.3807E-06, -.5328E-05, -.4186E-06, 
     .   -.3996E-05, .2074E-04, -.3540E-04, -.3083E-04, -.1789E-04, 
     1   -.9896E-05, -.1104E-04, -.6470E-05, -.5709E-05, .3425E-05, 
     2   -.4567E-05, .3463E-05, .5633E-05, -.3159E-05, -.3730E-04, 
     3   -.2740E-04, -.1484E-04, -.1066E-04, -.1142E-04, -.6470E-05, 
     4   -.6470E-05, .1522E-05, -.1522E-05, -.3045E-05, .3197E-05, 
     5   -.1039E-04, -.3425E-04, -.2284E-04, -.1370E-04, -.1028E-04, 
     6   -.1104E-04, -.8373E-05, -.4948E-05, .1903E-05, -.7612E-06, 
     7   -.1104E-05, .2455E-05, -.3805E-07, -.3235E-04, -.2093E-04, 
     8   -.1294E-04, -.1142E-04, -.1180E-04, -.6851E-05, -.3045E-05, 
     9   -.7611E-06, .1256E-05, -.7231E-06, .9924E-05, .3578E-05/ 
      END 
 
      BLOCK DATA CKD4
      COMMON /BAND4/HK(7), COEH2O(3,11,7)
      DATA HK/ .52, .21, .11, .1, .04, .015, .005/ 
      DATA (((COEH2O(K,J,I),I=1,7),J=1,11),K=1,3)/ -.1722E+02, 
     1   -.1402E+02, -.1202E+02, -.1001E+02, -.7702E+01, -.5273E+01, 
     2   -.6530E+00, -.1677E+02, -.1359E+02, -.1164E+02, -.9662E+01, 
     3   -.7419E+01, -.5001E+01, -.6040E+00, -.1630E+02, -.1316E+02, 
     4   -.1125E+02, -.9303E+01, -.7092E+01, -.4750E+01, -.5715E+00, 
     5   -.1584E+02, -.1274E+02, -.1086E+02, -.8939E+01, -.6751E+01, 
     6   -.4458E+01, -.4928E+00, -.1538E+02, -.1232E+02, -.1048E+02, 
     7   -.8579E+01, -.6399E+01, -.4191E+01, -.4683E+00, -.1493E+02, 
     8   -.1192E+02, -.1011E+02, -.8241E+01, -.6065E+01, -.3910E+01, 
     9   -.4310E+00, -.1440E+02, -.1145E+02, -.9643E+01, -.7873E+01, 
     .   -.5710E+01, -.3668E+01, -.3304E+00, -.1391E+02, -.1104E+02, 
     1   -.9238E+01, -.7479E+01, -.5367E+01, -.3387E+01, -.3604E+00, 
     2   -.1348E+02, -.1069E+02, -.8918E+01, -.7122E+01, -.5086E+01, 
     3   -.3152E+01, -.3030E+00, -.1310E+02, -.1037E+02, -.8626E+01, 
     4   -.6790E+01, -.4815E+01, -.2945E+01, -.4789E+00, -.1275E+02, 
     5   -.1011E+02, -.8347E+01, -.6484E+01, -.4584E+01, -.2788E+01, 
     6   -.5807E+00, .7934E-02, .9231E-02, .1005E-01, .9043E-02, 
     7   .8164E-02, .8980E-02, .6403E-02, .7954E-02, .9169E-02, 
     8   .9797E-02, .8687E-02, .7724E-02, .7954E-02, .6652E-02, 
     9   .7954E-02, .9043E-02, .9608E-02, .8499E-02, .7347E-02, 
     .   .7473E-02, .6382E-02, .7996E-02, .8980E-02, .9378E-02, 
     1   .8289E-02, .7264E-02, .6594E-02, .6674E-02, .8059E-02, 
     2   .8938E-02, .9294E-02, .8227E-02, .7201E-02, .6678E-02, 
     3   .7032E-02, .8122E-02, .8896E-02, .9189E-02, .8038E-02, 
     4   .7033E-02, .5987E-02, .5475E-02, .8268E-02, .9064E-02, 
     5   .8792E-02, .7975E-02, .6573E-02, .5087E-02, .4657E-02, 
     6   .8541E-02, .8980E-02, .9085E-02, .7996E-02, .6133E-02, 
     7   .4501E-02, .3860E-02, .8813E-02, .9043E-02, .9294E-02, 
     8   .8122E-02, .5861E-02, .4354E-02, .3964E-02, .8875E-02, 
     9   .8834E-02, .9797E-02, .8164E-02, .5463E-02, .4417E-02, 
     .   .3270E-02, .8938E-02, .8771E-02, .1005E-01, .8247E-02, 
     1   .5589E-02, .4835E-02, .3033E-02, -.1484E-04, -.2169E-04, 
     2   -.2436E-04, -.2588E-04, -.1142E-04, -.1142E-05, -.1519E-04, 
     3   -.1522E-04, -.2055E-04, -.2131E-04, -.2398E-04, -.4948E-05, 
     4   -.1675E-04, -.3593E-04, -.1522E-04, -.2055E-04, -.1865E-04, 
     5   -.2207E-04, -.4948E-05, -.1180E-04, -.1237E-04, -.1598E-04, 
     6   -.2017E-04, -.1903E-04, -.2284E-04, -.1028E-04, -.1865E-04, 
     7   -.2381E-04, -.1713E-04, -.2017E-04, -.1827E-04, -.2169E-04, 
     8   -.1218E-04, -.9515E-05, -.2415E-04, -.1827E-04, -.2093E-04, 
     9   -.1637E-04, -.1827E-04, -.9134E-05, -.8373E-05, -.1243E-04, 
     .   -.1560E-04, -.1865E-04, -.1599E-04, -.1256E-04, -.1066E-04, 
     1   -.1142E-05, -.2181E-04, -.1675E-04, -.1560E-04, -.1522E-04, 
     2   -.1675E-04, -.1865E-04, -.1865E-04, -.9522E-05, -.1332E-04, 
     3   -.1370E-04, -.1446E-04, -.2055E-04, -.1142E-04, -.2512E-04, 
     4   -.3343E-04, -.1294E-04, -.1294E-04, -.1751E-04, -.2512E-04, 
     5   -.1560E-04, -.2854E-04, -.7003E-05, -.8753E-05, -.1028E-04, 
     6   -.1751E-04, -.2512E-04, -.1713E-04, -.1713E-04,  - .1245E-04/ 
      END 
 
      BLOCK DATA CKD5
      COMMON /BAND5/HK(12), COEH2O(3,11,12)
      DATA HK/.13,.14,.13,.16,.18,.14,.07,.02,.016,.008,.004,.002/
      DATA (((COEH2O(K,J,I),I=1,12),J=1,11),K=1,3)/ -.1499E+02, 
     1   -.1267E+02, -.1118E+02, -.9696E+01, -.7992E+01, -.6323E+01, 
     2   -.4414E+01, -.2961E+01, -.1715E+01, -.1406E+00, .1612E+01, 
     3   .3689E+01, -.1454E+02, -.1223E+02, -.1075E+02, -.9277E+01, 
     4   -.7576E+01, -.5915E+01, -.4043E+01, -.2630E+01, -.1449E+01, 
     5   .2314E-01, .1708E+01, .3744E+01, -.1408E+02, -.1178E+02, 
     6   -.1031E+02, -.8851E+01, -.7154E+01, -.5503E+01, -.3666E+01, 
     7   -.2288E+01, -.1141E+01, .2772E+00, .1819E+01, .3788E+01, 
     8   -.1363E+02, -.1134E+02, -.9876E+01, -.8423E+01, -.6733E+01, 
     9   -.5091E+01, -.3286E+01, -.1938E+01, -.8649E+00, .5349E+00, 
     .   .1969E+01, .3795E+01, -.1318E+02, -.1091E+02, -.9452E+01, 
     1   -.8004E+01, -.6309E+01, -.4677E+01, -.2904E+01, -.1595E+01, 
     2   -.5641E+00, .7592E+00, .2109E+01, .3783E+01, -.1275E+02, 
     3   -.1048E+02, -.9028E+01, -.7585E+01, -.5892E+01, -.4267E+01, 
     4   -.2524E+01, -.1274E+01, -.2782E+00, .9376E+00, .2257E+01, 
     5   .3714E+01, -.1180E+02, -.9887E+01, -.8492E+01, -.7014E+01, 
     6   -.5390E+01, -.3834E+01, -.2156E+01, -.9775E+00, -.3129E-01, 
     7   .1151E+01, .2330E+01, .3592E+01, -.1114E+02, -.9367E+01, 
     8   -.8002E+01, -.6514E+01, -.4928E+01, -.3435E+01, -.1835E+01, 
     9   -.7064E+00, .2153E+00, .1309E+01, .2422E+01, .3488E+01, 
     .   -.1074E+02, -.8941E+01, -.7582E+01, -.6116E+01, -.4536E+01, 
     1   -.3072E+01, -.1521E+01, -.4651E+00, .4053E+00, .1465E+01, 
     2   .2374E+01, .3260E+01, -.1041E+02, -.8545E+01, -.7180E+01, 
     3   -.5745E+01, -.4177E+01, -.2735E+01, -.1245E+01, -.2356E+00, 
     4   .5786E+00, .1516E+01, .2263E+01, .3074E+01, -.1008E+02, 
     5   -.8149E+01, -.6804E+01, -.5409E+01, -.3855E+01, -.2427E+01, 
     6   -.9857E+00, -.4939E-01, .7060E+00, .1483E+01, .2159E+01, 
     7   .2745E+01, .9985E-02, .8373E-02, .7431E-02, .6866E-02, 
     8   .4584E-02, .2952E-02, .3098E-02, .3768E-02, .4013E-02, 
     9   .3960E-02, .3228E-02, .3203E-02, .1007E-01, .8436E-02, 
     .   .7368E-02, .6657E-02, .4375E-02, .2617E-02, .2742E-02, 
     1   .3286E-02, .3192E-02, .2992E-02, .2612E-02, .1968E-02, 
     2   .1019E-01, .8457E-02, .7264E-02, .6426E-02, .4187E-02, 
     3   .2365E-02, .2324E-02, .2614E-02, .2736E-02, .2068E-02, 
     4   .2085E-02, .1005E-02, .1028E-01, .8478E-02, .7138E-02, 
     5   .6259E-02, .3998E-02, .2156E-02, .1926E-02, .1953E-02, 
     6   .2250E-02, .1844E-02, .1869E-02, -.6489E-03, .1030E-01, 
     7   .8478E-02, .7033E-02, .6112E-02, .3852E-02, .1989E-02, 
     8   .1716E-02, .1763E-02, .1432E-02, .1193E-02, .1306E-02, 
     9   -.5861E-03, .1042E-01, .8499E-02, .6887E-02, .5987E-02, 
     .   .3768E-02, .1800E-02, .1549E-02, .1712E-02, .1287E-02, 
     1   .7389E-03, .7222E-03, -.1130E-02, .8227E-02, .7201E-02, 
     2   .6866E-02, .5903E-02, .3412E-02, .1591E-02, .1402E-02, 
     3   .1346E-02, .1041E-02, .8185E-03, .3349E-03, -.4815E-03, 
     4   .8268E-02, .6992E-02, .7159E-02, .6384E-02, .3286E-02, 
     5   .1591E-02, .1271E-02, .1202E-02, .9187E-03, .6531E-03, 
     6   -.4187E-03, -.7954E-03, .8478E-02, .7159E-02, .7117E-02, 
     7   .6447E-02, .3349E-02, .1528E-02, .9964E-03, .9210E-03, 
     8   .6112E-03, .6259E-03, -.3768E-03, -.1298E-02, .8520E-02, 
     9   .7075E-02, .7096E-02, .6405E-02, .3245E-02, .1528E-02, 
     .   .1011E-02, .7877E-03, .7536E-03, .9001E-04, -.6719E-03, 
     1   -.1026E-02, .8561E-02, .6950E-02, .7033E-02, .6280E-02, 
     2   .2993E-02, .1528E-02, .6698E-03, .5847E-03, .2847E-03, 
     3   -.6280E-04, -.9420E-03, -.1444E-02, -.1408E-04, -.2664E-04, 
     4   -.1180E-04, -.1903E-04, -.9515E-05, .3806E-06, -.6851E-05, 
     5   -.3806E-05, -.4834E-05, -.3239E-05, -.2284E-05, -.1028E-04, 
     6   -.1484E-04, -.2550E-04, -.1142E-04, -.1827E-04, -.9515E-05, 
     7   .3805E-06, -.4948E-05, .3806E-06, -.2664E-06, .1058E-04, 
     8   -.1012E-04, -.1142E-04, -.1560E-04, -.2512E-04, -.1256E-04, 
     9   -.1865E-04, -.9134E-05, .1142E-05, -.3425E-05, .2474E-05, 
     .   -.9781E-05, -.1519E-05, -.7916E-05, -.1294E-04, -.1560E-04, 
     1   -.2474E-04, -.1180E-04, -.2017E-04, -.7992E-05, .3805E-06, 
     2   -.2283E-05, -.4453E-05, -.1180E-05, -.5138E-05, -.4453E-05, 
     3   -.3425E-05, -.1522E-04, -.2550E-04, -.9896E-05, -.1903E-04, 
     4   -.9134E-05, -.1142E-05, -.7611E-06, -.5252E-05, -.4567E-06, 
     5   -.4643E-05, -.4567E-06, -.4567E-05, -.1294E-04, -.2512E-04, 
     6   -.1028E-04, -.2055E-04, -.9896E-05, -.4567E-05, -.2284E-05, 
     7   -.5100E-05, -.4339E-06, -.9515E-06, -.1252E-04, -.7612E-06, 
     8   -.2246E-04, -.1370E-04, -.1066E-04, -.1598E-04, -.8754E-05, 
     9   -.5328E-05, -.6622E-05, -.5138E-05, -.8754E-07, -.9515E-06, 
     .   .6090E-05, .4187E-05, -.3463E-04, -.1599E-04, -.1218E-04, 
     1   -.2093E-04, -.9515E-05, -.4567E-05, -.1104E-05, -.1903E-05, 
     2   -.1488E-05, -.3730E-05, -.4567E-05, .3045E-05, -.3463E-04, 
     3   -.1675E-04, -.1294E-04, -.1979E-04, -.1066E-04, -.4187E-05, 
     4   -.4034E-05, -.2893E-05, -.2588E-05, -.9401E-05, .2284E-05, 
     5   .3045E-05, -.2778E-04, -.1522E-04, -.1560E-04, -.1751E-04, 
     6   -.1256E-04, -.5709E-05, -.2474E-05, -.2577E-05, -.2284E-05, 
     7   -.4187E-06, .7650E-05, -.3425E-05, -.3083E-04, -.1827E-04, 
     8   -.1370E-04, -.1751E-04, -.1104E-04, -.9515E-05, -.6318E-05, 
     9   -.4358E-05, -.7613E-07, .4643E-05, .4415E-05, .1028E-04/ 
      END 
 
      BLOCK DATA CKD6
      COMMON /BAND6/HK(5), COEH2O(3,11,5)
      DATA HK/ .3, .2, .2, .2, .1/ 
      DATA (((COEH2O(K,J,I),I=1,5),J=1,11),K=1,3)/ -.1905E+02, 
     1   -.1602E+02, -.1472E+02, -.1307E+02, -.1024E+02, -.1823E+02, 
     2   -.1555E+02, -.1427E+02, -.1266E+02, -.9938E+01, -.1749E+02, 
     3   -.1508E+02, -.1381E+02, -.1225E+02, -.9641E+01, -.1684E+02, 
     4   -.1462E+02, -.1337E+02, -.1185E+02, -.9367E+01, -.1630E+02, 
     5   -.1417E+02, -.1294E+02, -.1145E+02, -.9123E+01, -.1578E+02, 
     6   -.1373E+02, -.1251E+02, -.1108E+02, -.8881E+01, -.1517E+02, 
     7   -.1327E+02, -.1209E+02, -.1072E+02, -.8653E+01, -.1463E+02, 
     8   -.1284E+02, -.1169E+02, -.1040E+02, -.8453E+01, -.1421E+02, 
     9   -.1244E+02, -.1133E+02, -.1014E+02, -.8312E+01, -.1382E+02, 
     .   -.1207E+02, -.1100E+02, -.9887E+01, -.8220E+01, -.1348E+02, 
     1   -.1173E+02, -.1071E+02, -.9685E+01, -.8220E+01, .1024E-01, 
     2   .1842E-02, .6908E-03, .1737E-02, .3517E-02, .8394E-02, 
     3   .2072E-02, .8164E-03, .1716E-02, .2805E-02, .8143E-02, 
     4   .2240E-02, .9001E-03, .1570E-02, .1800E-02, .8227E-02, 
     5   .2386E-02, .9420E-03, .1486E-02, .1068E-02, .8373E-02, 
     6   .2533E-02, .9210E-03, .1319E-02, .9420E-03, .8394E-02, 
     7   .2700E-02, .9629E-03, .1026E-02, .5233E-03, .8917E-02, 
     8   .2575E-02, .8792E-03, .7536E-03, .4187E-03, .9378E-02, 
     9   .2617E-02, .7955E-03, .6070E-03, .4815E-03, .9797E-02, 
     .   .2638E-02, .6908E-03, .5233E-03, .6280E-03, .1009E-01, 
     1   .2638E-02, .4815E-03, .2931E-03, .4815E-03, .1036E-01, 
     2   .2428E-02, .3140E-03, .3977E-03, .2093E-03, -.5366E-04, 
     3   -.1522E-04, -.5709E-05, -.2664E-05, .3806E-05, -.4301E-04, 
     4   -.1484E-04, -.4948E-05, -.7610E-06, .7610E-06, -.3920E-04, 
     5   -.1484E-04, -.4948E-05, .3804E-06, -.3806E-05, -.3920E-04, 
     6   -.1522E-04, -.4948E-05, .3425E-05, .1903E-05, -.3806E-04, 
     7   -.1484E-04, -.3045E-05, .2664E-05, .7993E-05, -.4148E-04, 
     8   -.1408E-04, -.3806E-05, .4187E-05, .7993E-05, -.5481E-04, 
     9   -.1180E-04, -.3045E-05, .3045E-05, .2284E-05, -.5709E-04, 
     .   -.1104E-04, -.2283E-05, -.2664E-05, -.1142E-05, -.6090E-04, 
     1   -.1218E-04, -.2664E-05, .3804E-06, .3045E-05, -.6698E-04, 
     2   -.1218E-04, -.2664E-05, .1523E-05, -.1142E-05, -.6508E-04, 
     3   -.1218E-04, -.3425E-05, .1903E-05, .7612E-06/ 
      END 
 
      BLOCK DATA CKD7
      COMMON /BAND7/HK(2), COEH2O(3,19,2)
      DATA HK/ 0.7, 0.3/ 
      DATA (((COEH2O(K,J,I),I=1,2),J=1,19),K=1,3)/ -.2008E+02, 
     1   -.1467E+02, -.2004E+02, -.1426E+02, -.2001E+02, -.1386E+02, 
     2   -.1998E+02, -.1345E+02, -.1995E+02, -.1304E+02, -.1992E+02, 
     3   -.1263E+02, -.1989E+02, -.1223E+02, -.1986E+02, -.1183E+02, 
     4   -.1984E+02, -.1143E+02, -.1758E+02, -.1038E+02, -.1602E+02, 
     5   -.9480E+01, -.1469E+02, -.8752E+01, -.1349E+02, -.8218E+01, 
     6   -.1255E+02, -.7677E+01, -.1174E+02, -.7184E+01, -.1110E+02, 
     7   -.6735E+01, -.1056E+02, -.6332E+01, -.1019E+02, -.5975E+01, 
     8   -.9874E+01, -.5644E+01, .2533E-02, .2269E-01, .2575E-02, 
     9   .2263E-01, .2554E-02, .2267E-01, .2491E-02, .2250E-01, 
     .   .2449E-02, .2244E-01, .2344E-02, .2234E-01, .2219E-02, 
     1   .2208E-01, .5694E-02, .2190E-01, .9650E-02, .2162E-01, 
     2   .3286E-01, .1848E-01, .2987E-01, .1578E-01, .2527E-01, 
     3   .1465E-01, .2175E-01, .1386E-01, .2056E-01, .1235E-01, 
     4   .1963E-01, .1116E-01, .1926E-01, .1040E-01, .2014E-01, 
     5   .1040E-01, .2024E-01, .1042E-01, .1972E-01, .1080E-01, 
     6   -.8754E-05, -.6698E-04, -.1104E-04, -.6432E-04, -.1142E-04, 
     7   -.6051E-04, -.1180E-04, -.6128E-04, -.1180E-04, -.6242E-04, 
     8   -.1218E-04, -.6280E-04, -.1218E-04, -.6204E-04, .5328E-04, 
     9   -.5709E-04, .1275E-03, -.5214E-04, -.1370E-03, -.4148E-04, 
     .   -.1100E-03, -.3045E-04, -.9248E-04, -.3197E-04, -.7346E-04, 
     1   -.2436E-04, -.5100E-04, -.2131E-04, -.5861E-04, -.2550E-04, 
     2   -.5328E-04, -.3311E-04, -.6090E-04, -.4225E-04, -.5443E-04, 
     3   -.4415E-04, -.4034E-04,  - .4339E-04/ 
      END 
 
      BLOCK DATA CKD8
      COMMON /BAND8/HK(3), COEH2O(3,19,3)
      DATA HK/ 0.2, 0.7, 0.1/ 
      DATA (((COEH2O(K,J,I),I=1,3),J=1,19),K=1,3)/ -.2283E+02, 
     1   -.1639E+02, -.6155E+01, -.2237E+02, -.1595E+02, -.5775E+01, 
     2   -.2191E+02, -.1551E+02, -.5381E+01, -.2145E+02, -.1507E+02, 
     3   -.5004E+01, -.2099E+02, -.1463E+02, -.4617E+01, -.2053E+02, 
     4   -.1419E+02, -.4218E+01, -.2025E+02, -.1375E+02, -.3806E+01, 
     5   -.2021E+02, -.1330E+02, -.3403E+01, -.2018E+02, -.1287E+02, 
     6   -.2993E+01, -.1998E+02, -.1091E+02, -.2586E+01, -.1744E+02, 
     7   -.9171E+01, -.2162E+01, -.1490E+02, -.7642E+01, -.1763E+01, 
     8   -.1303E+02, -.6526E+01, -.1373E+01, -.1113E+02, -.5846E+01, 
     9   -.9699E+00, -.9814E+01, -.5280E+01, -.5955E+00, -.8582E+01, 
     .   -.4787E+01, -.2510E+00, -.8020E+01, -.4350E+01, .2770E-01, 
     1   -.7571E+01, -.3942E+01, .2406E+00, -.7140E+01, -.3537E+01, 
     2   .3567E+00, .3722E-01, .1505E-01, .6615E-02, .3722E-01, 
     3   .1518E-01, .5840E-02, .3720E-01, .1526E-01, .5170E-02, 
     4   .3399E-01, .1530E-01, .4773E-02, .3012E-01, .1551E-01, 
     5   .4333E-02, .2625E-01, .1553E-01, .3956E-02, .2240E-01, 
     6   .1562E-01, .3454E-02, .1846E-01, .1574E-01, .3161E-02, 
     7   .1446E-01, .1572E-01, .3098E-02, .5924E-02, .8875E-02, 
     8   .2658E-02, .2204E-01, .7096E-02, .2504E-02, .1591E-01, 
     9   .5233E-02, .2292E-02, .8855E-02, .4249E-02, .2190E-02, 
     .   .5422E-02, .3496E-02, .2041E-02, .4919E-02, .3621E-02, 
     1   .2200E-02, .6657E-02, .3663E-02, .2248E-02, .8645E-02, 
     2   .3852E-02, .2118E-02, .8771E-02, .3873E-02, .2176E-02, 
     3   .9043E-02, .3747E-02, .2079E-02, -.1568E-03, -.4681E-04, 
     4   .4567E-05, -.1568E-03, -.4605E-04, -.3425E-05, -.1572E-03, 
     5   -.4605E-04, -.1104E-04, -.2154E-03, -.4453E-04, -.6851E-05, 
     6   -.2843E-03, -.4225E-04, -.7231E-05, -.3562E-03, -.4110E-04, 
     7   -.7231E-05, -.3692E-03, -.4110E-04, -.1028E-04, -.3007E-03, 
     8   -.4263E-04, -.6470E-05, -.2325E-03, -.3996E-04, -.8373E-05, 
     9   -.5290E-04, -.7612E-05, -.4948E-05, -.7422E-04, -.1256E-04, 
     .   -.8449E-05, -.3501E-04, -.1446E-04, -.4834E-05, .4529E-04, 
     1   -.2246E-04, -.2893E-05, .6470E-05, -.1789E-04, -.7498E-05, 
     2   -.4948E-05, -.1713E-04, -.8183E-05, -.5481E-04, -.1713E-04, 
     3   -.1447E-04, -.4986E-04, -.1903E-04, -.1353E-04, -.5138E-04, 
     4   -.1484E-04, -.1147E-04, -.5328E-04, -.1560E-04,  - .6588E-05/ 
      END 
 
      BLOCK DATA CKD9
      COMMON /BAND9/HK(4), COEH2O(3,19,4)
      DATA HK/ 0.22, 0.51, 0.22, 0.05/ 
      DATA (((COEH2O(K,J,I),I=1,4),J=1,19),K=1,3)/ -.2066E+02, 
     1   -.1464E+02, -.8301E+01, -.3548E+01, -.2025E+02, -.1419E+02, 
     2   -.7905E+01, -.3260E+01, -.2019E+02, -.1374E+02, -.7495E+01, 
     3   -.2927E+01, -.2013E+02, -.1329E+02, -.7078E+01, -.2584E+01, 
     4   -.2007E+02, -.1284E+02, -.6675E+01, -.2247E+01, -.2001E+02, 
     5   -.1239E+02, -.6268E+01, -.1890E+01, -.1996E+02, -.1194E+02, 
     6   -.5853E+01, -.1530E+01, -.1991E+02, -.1150E+02, -.5441E+01, 
     7   -.1133E+01, -.1987E+02, -.1105E+02, -.5022E+01, -.7447E+00, 
     8   -.1575E+02, -.9657E+01, -.4191E+01, -.3728E+00, -.1329E+02, 
     9   -.8133E+01, -.3638E+01, .1616E-01, -.1181E+02, -.6675E+01, 
     .   -.3178E+01, .4083E+00, -.1036E+02, -.5655E+01, -.2731E+01, 
     1   .7953E+00, -.8628E+01, -.4990E+01, -.2303E+01, .1153E+01, 
     2   -.7223E+01, -.4453E+01, -.1877E+01, .1454E+01, -.6567E+01, 
     3   -.3974E+01, -.1461E+01, .1663E+01, -.6077E+01, -.3551E+01, 
     4   -.1071E+01, .1800E+01, -.5651E+01, -.3136E+01, -.7005E+00, 
     5   .1809E+01, -.5241E+01, -.2726E+01, -.3859E+00, .1781E+01, 
     6   .1315E-01, .4542E-02, .3496E-02, .4877E-02, .9650E-02, 
     7   .4542E-02, .3098E-02, .3956E-02, .6154E-02, .4626E-02, 
     8   .2763E-02, .3077E-02, .2658E-02, .4626E-02, .2512E-02, 
     9   .2261E-02, .2658E-02, .4689E-02, .2219E-02, .1405E-02, 
     .   .2700E-02, .4752E-02, .1926E-02, .7473E-03, .2658E-02, 
     1   .4773E-02, .1737E-02, .5066E-03, .4668E-02, .4815E-02, 
     2   .1507E-02, .1842E-03, .8541E-02, .4794E-02, .1382E-02, 
     3   -.2156E-03, .1022E-01, .2198E-02, .3977E-03, -.2910E-03, 
     4   .5484E-02, .6698E-03, .0000E+00, -.2339E-03, .3349E-02, 
     5   .1068E-02, -.2512E-03, -.4228E-03, .1884E-02, .2093E-03, 
     6   -.3977E-03, -.6405E-03, -.8373E-04, -.5233E-03, -.4124E-03, 
     7   -.5945E-03, .7536E-03, -.6698E-03, -.4919E-03, -.4794E-03, 
     8   .3600E-02, -.4605E-03, -.4375E-03, -.3517E-03, .3873E-02, 
     9   -.5861E-03, -.3203E-03, -.4689E-03, .3935E-02, -.7326E-03, 
     .   -.2072E-03, -.4228E-03, .4124E-02, -.8582E-03, -.4187E-04, 
     1   -.5945E-03, -.8525E-04, .1865E-04, -.1142E-05, .2664E-05, 
     2   -.1313E-03, .1865E-04, .0000E+00, .1256E-04, -.6470E-04, 
     3   .1865E-04, -.3045E-05, .8754E-05, .3805E-06, .1789E-04, 
     4   -.6851E-05, .5328E-05, .1142E-05, .1827E-04, -.6090E-05, 
     5   .4148E-05, .1142E-05, .1865E-04, -.3806E-05, -.3768E-05, 
     6   -.1903E-05, .1751E-04, -.4948E-05, .3121E-05, .3159E-04, 
     7   .1979E-04, -.3045E-05, -.9896E-06, .1005E-03, .1789E-04, 
     8   -.6089E-05, -.1865E-05, -.2207E-04, .1941E-04, .1903E-05, 
     9   .2322E-05, -.1675E-04, .6090E-05, -.7611E-06, .4397E-05, 
     .   .3425E-04, .3806E-06, .1522E-05, .3806E-05, .4796E-04, 
     1   .1522E-05, -.3806E-06, .3654E-05, -.6851E-05, .2664E-05, 
     2   -.3920E-05, -.6850E-06, -.1370E-04, .5328E-05, -.6584E-05, 
     3   -.8716E-05, -.8374E-10, .1522E-05, -.6356E-05, .1294E-05, 
     4   -.9515E-05, .7612E-06, -.3235E-05, -.1066E-05, -.7612E-05, 
     5   .1142E-05, -.4529E-05, .3730E-05, -.2664E-05, -.3806E-06, 
     6   -.3501E-05,  - .5328E-06/ 
      END 
 
      BLOCK DATA CKD10
      COMMON /BAND10/HK(4), COEH2O(3,19,4), COECH4(3,19), COEN2O(3,19)
      DATA HK/ 0.28, 0.42, 0.25, 0.05/ 
      DATA (((COEH2O(K,J,I),I=1,4),J=1,19),K=1,3)/ -.2023E+02, 
     1   -.1641E+02, -.1171E+02, -.6090E+01, -.2016E+02, -.1595E+02, 
     2   -.1133E+02, -.5867E+01, -.2011E+02, -.1550E+02, -.1095E+02, 
     3   -.5660E+01, -.2005E+02, -.1504E+02, -.1055E+02, -.5407E+01, 
     4   -.2001E+02, -.1459E+02, -.1015E+02, -.5137E+01, -.1997E+02, 
     5   -.1413E+02, -.9749E+01, -.4852E+01, -.1993E+02, -.1367E+02, 
     6   -.9337E+01, -.4534E+01, -.1990E+02, -.1321E+02, -.8920E+01, 
     7   -.4211E+01, -.1987E+02, -.1276E+02, -.8506E+01, -.3889E+01, 
     8   -.1645E+02, -.1179E+02, -.7711E+01, -.3613E+01, -.1442E+02, 
     9   -.1081E+02, -.6942E+01, -.3316E+01, -.1308E+02, -.9950E+01, 
     .   -.6344E+01, -.2950E+01, -.1212E+02, -.9217E+01, -.5904E+01, 
     1   -.2577E+01, -.1131E+02, -.8559E+01, -.5519E+01, -.2256E+01, 
     2   -.1064E+02, -.7962E+01, -.5183E+01, -.1929E+01, -.1013E+02, 
     3   -.7447E+01, -.4833E+01, -.1643E+01, -.9712E+01, -.7071E+01, 
     4   -.4485E+01, -.1410E+01, -.9305E+01, -.6760E+01, -.4145E+01, 
     5   -.1249E+01, -.8966E+01, -.6477E+01, -.3820E+01, -.1114E+01, 
     6   .7913E-02, .8206E-02, .1509E-01, .1869E-01, .4228E-02, 
     7   .8247E-02, .1467E-01, .1783E-01, .2010E-02, .8227E-02, 
     8   .1442E-01, .1687E-01, .1947E-02, .8289E-02, .1394E-01, 
     9   .1568E-01, .1863E-02, .8289E-02, .1346E-01, .1484E-01, 
     .   .1842E-02, .8415E-02, .1310E-01, .1400E-01, .1800E-02, 
     1   .8457E-02, .1275E-01, .1377E-01, .1696E-02, .8478E-02, 
     2   .1220E-01, .1321E-01, .1842E-02, .8478E-02, .1189E-01, 
     3   .1250E-01, .1409E-01, .8624E-02, .1254E-01, .1214E-01, 
     4   .9043E-02, .1045E-01, .1225E-01, .1260E-01, .8561E-02, 
     5   .1202E-01, .1181E-01, .1296E-01, .1114E-01, .1235E-01, 
     6   .1191E-01, .1330E-01, .1199E-01, .1271E-01, .1195E-01, 
     7   .1371E-01, .1415E-01, .1315E-01, .1218E-01, .1361E-01, 
     8   .1478E-01, .1338E-01, .1296E-01, .1306E-01, .1518E-01, 
     9   .1375E-01, .1365E-01, .1334E-01, .1530E-01, .1411E-01, 
     .   .1392E-01, .1327E-01, .1547E-01, .1507E-01, .1390E-01, 
     1   .1264E-01, -.1089E-03, -.2740E-04, -.2017E-04, -.5519E-04, 
     2   -.4491E-04, -.2740E-04, -.1408E-04, -.5937E-04, -.6090E-05, 
     3   -.2702E-04, -.6470E-05, -.4719E-04, -.7232E-05, -.2740E-04, 
     4   -.6089E-05, -.4910E-04, -.7231E-05, -.2969E-04, -.4186E-05, 
     5   -.5366E-04, -.6090E-05, -.3045E-04, -.2284E-05, -.4986E-04, 
     6   -.4568E-05, -.3121E-04, -.4948E-05, -.5100E-04, -.3426E-05, 
     7   -.3007E-04, -.7993E-05, -.4910E-04, .1522E-05, -.2931E-04, 
     8   -.9896E-05, -.5366E-04, -.5823E-04, -.1599E-04, -.1713E-04, 
     9   -.4110E-04, -.3121E-04, -.1713E-04, -.3159E-04, -.3578E-04, 
     .   -.3996E-04, -.1598E-04, -.3958E-04, -.4605E-04, -.3349E-04, 
     1   -.1751E-04, -.3844E-04, -.5576E-04, -.2626E-04, -.2474E-04, 
     2   -.3920E-04, -.4464E-04, -.1979E-04, -.3045E-04, -.3958E-04, 
     3   -.5336E-04, -.2893E-04, -.3616E-04, -.3996E-04, -.4754E-04, 
     4   -.2398E-04, -.3083E-04, -.4415E-04, -.5119E-04, -.2702E-04, 
     5   -.2664E-04, -.4605E-04, -.4038E-04, -.2398E-04, -.2360E-04, 
     6   -.4948E-04,  - .5149E-04/ 
      DATA ((COECH4(K,J),J=1,19),K=1,3)/ -.8909E+01, -.8464E+01, 
     1   -.8018E+01, -.7573E+01, -.7133E+01, -.6687E+01, -.6240E+01, 
     2   -.5803E+01, -.5377E+01, -.4534E+01, -.3983E+01, -.3502E+01, 
     3   -.3062E+01, -.2648E+01, -.2265E+01, -.1896E+01, -.1568E+01, 
     4   -.1234E+01, -.9298E+00, .9629E-03, .9838E-03, .1088E-02, 
     5   .1172E-02, .1256E-02, .1402E-02, .1528E-02, .1633E-02, 
     6   .1716E-02, .4815E-03, -.3977E-03, -.5652E-03, -.5024E-03, 
     7   -.4605E-03, -.4563E-03, -.4438E-03, -.4521E-03, -.4312E-03, 
     8   -.3789E-03, -.1294E-04, -.1408E-04, -.1522E-04, -.1675E-04, 
     9   -.1751E-04, -.1941E-04, -.2246E-04, -.2207E-04, -.1827E-04, 
     .   -.1256E-04, -.9515E-05, -.6470E-05, -.3045E-05, -.3806E-05, 
     1   -.2055E-05, -.3730E-05, -.7612E-06, -.3806E-05, .1256E-05/ 
      DATA ((COEN2O(K,J),J=1,19),K=1,3)/ -.7863E+01, -.7412E+01, 
     1   -.6963E+01, -.6514E+01, -.6065E+01, -.5611E+01, -.5167E+01, 
     2   -.4720E+01, -.4283E+01, -.3454E+01, -.2858E+01, -.2404E+01, 
     3   -.1922E+01, -.1491E+01, -.1097E+01, -.7177E+00, -.3548E+00, 
     4   .1218E-01, .3088E+00, .4459E-02, .4542E-02, .4668E-02, 
     5   .4752E-02, .4815E-02, .4919E-02, .5087E-02, .5254E-02, 
     6   .5296E-02, .2324E-02, .2093E-02, .2294E-02, .2125E-02, 
     7   .2058E-02, .1920E-02, .1786E-02, .1689E-02, .1788E-02, 
     8   .2144E-02, -.7231E-05, -.7231E-05, -.7231E-05, -.6470E-05, 
     9   -.6851E-05, -.7231E-05, -.5709E-05, -.6470E-05, -.4186E-05, 
     .   .8754E-05, -.7612E-05, -.9134E-06, -.8640E-05, -.8487E-05, 
     1   -.8259E-05, -.9553E-05, -.8107E-05, -.1654E-04,  - .1858E-04/ 
      END 
 
      BLOCK DATA CKD11
      COMMON /BAND11/HK(3), COEH2O(3,19,3), COECH4(3,19), COEN2O(3,19)
      DATA HK/ 0.80, 0.15, 0.05/ 
      DATA (((COEH2O(K,J,I),I=1,3),J=1,19),K=1,3)/ -.2005E+02, 
     1   -.1548E+02, -.1021E+02, -.2001E+02, -.1504E+02, -.1001E+02, 
     2   -.1997E+02, -.1459E+02, -.9814E+01, -.1993E+02, -.1416E+02, 
     3   -.9595E+01, -.1989E+02, -.1373E+02, -.9349E+01, -.1985E+02, 
     4   -.1328E+02, -.9072E+01, -.1982E+02, -.1286E+02, -.8833E+01, 
     5   -.1957E+02, -.1243E+02, -.8566E+01, -.1911E+02, -.1200E+02, 
     6   -.8276E+01, -.1743E+02, -.1134E+02, -.7958E+01, -.1625E+02, 
     7   -.1078E+02, -.7629E+01, -.1524E+02, -.1036E+02, -.7334E+01, 
     8   -.1429E+02, -.9970E+01, -.7051E+01, -.1348E+02, -.9620E+01, 
     9   -.6749E+01, -.1282E+02, -.9270E+01, -.6505E+01, -.1229E+02, 
     .   -.8932E+01, -.6277E+01, -.1186E+02, -.8628E+01, -.6120E+01, 
     1   -.1148E+02, -.8345E+01, -.6049E+01, -.1112E+02, -.8066E+01, 
     2   -.5906E+01, .1842E-02, .2131E-01, .3033E-01, .1905E-02, 
     3   .2137E-01, .2841E-01, .1926E-02, .2135E-01, .2696E-01, 
     4   .1926E-02, .2133E-01, .2514E-01, .1884E-02, .2154E-01, 
     5   .2401E-01, .5589E-02, .2156E-01, .2321E-01, .9483E-02, 
     6   .2156E-01, .2210E-01, .1333E-01, .2150E-01, .2133E-01, 
     7   .1725E-01, .2154E-01, .2074E-01, .2254E-01, .1999E-01, 
     8   .2005E-01, .2118E-01, .1926E-01, .1978E-01, .1936E-01, 
     9   .1920E-01, .1963E-01, .1905E-01, .1911E-01, .1934E-01, 
     .   .1909E-01, .1903E-01, .1920E-01, .1922E-01, .1901E-01, 
     1   .1899E-01, .1934E-01, .1930E-01, .1974E-01, .1966E-01, 
     2   .1909E-01, .2014E-01, .1976E-01, .1905E-01, .1984E-01, 
     3   .1963E-01, .1940E-01, .1897E-01, -.1522E-05, -.6013E-04, 
     4   -.5062E-04, -.2665E-05, -.6204E-04, -.5519E-04, -.3806E-05, 
     5   -.6394E-04, -.5633E-04, -.4567E-05, -.6280E-04, -.5214E-04, 
     6   -.6090E-05, -.6128E-04, -.5290E-04, .6051E-04, -.6242E-04, 
     7   -.5823E-04, .1313E-03, -.6013E-04, -.5176E-04, .1336E-03, 
     8   -.5747E-04, -.4072E-04, .6318E-04, -.5671E-04, -.3996E-04, 
     9   -.5595E-04, -.3996E-04, -.4263E-04, -.3958E-04, -.4719E-04, 
     .   -.4453E-04, -.3387E-04, -.5138E-04, -.5100E-04, -.5252E-04, 
     1   -.4986E-04, -.4491E-04, -.5100E-04, -.4453E-04, -.4529E-04, 
     2   -.5176E-04, -.4795E-04, -.4453E-04, -.5557E-04, -.5176E-04, 
     3   -.5062E-04, -.5747E-04, -.4795E-04, -.5633E-04, -.5709E-04, 
     4   -.4643E-04, -.3806E-04, -.5481E-04, -.5671E-04,  - .4948E-04/ 
      DATA ((COECH4(K,J),J=1,19),K=1,3)/ -.1207E+02, -.1162E+02, 
     1   -.1116E+02, -.1070E+02, -.1024E+02, -.9777E+01, -.9319E+01, 
     2   -.8858E+01, -.8398E+01, -.7384E+01, -.6643E+01, -.6081E+01, 
     3   -.5602E+01, -.5188E+01, -.4822E+01, -.4479E+01, -.4184E+01, 
     4   -.3884E+01, -.3627E+01, .1036E-01, .1036E-01, .1040E-01, 
     5   .1040E-01, .1045E-01, .1047E-01, .1049E-01, .1055E-01, 
     6   .1059E-01, .1059E-01, .1026E-01, .1011E-01, .1024E-01, 
     7   .1049E-01, .1072E-01, .1089E-01, .1109E-01, .1153E-01, 
     8   .1191E-01, -.4910E-04, -.4834E-04, -.4910E-04, -.4910E-04, 
     9   -.4910E-04, -.4872E-04, -.4834E-04, -.4948E-04, -.5100E-04, 
     .   -.5633E-04, -.6166E-04, -.5595E-04, -.5366E-04, -.5366E-04, 
     1   -.5328E-04, -.5328E-04, -.4948E-04, -.5519E-04,  - .5595E-04/ 
      DATA ((COEN2O(K,J),J=1,19),K=1,3)/ -.9461E+01, -.9003E+01, 
     1   -.8543E+01, -.8084E+01, -.7629E+01, -.7166E+01, -.6707E+01, 
     2   -.6249E+01, -.5793E+01, -.5312E+01, -.4847E+01, -.4393E+01, 
     3   -.3974E+01, -.3587E+01, -.3231E+01, -.2885E+01, -.2602E+01, 
     4   -.2358E+01, -.2108E+01, .4710E-02, .4752E-02, .4773E-02, 
     5   .4773E-02, .4815E-02, .4877E-02, .4898E-02, .4982E-02, 
     6   .5066E-02, .5296E-02, .5149E-02, .5129E-02, .5024E-02, 
     7   .4752E-02, .4501E-02, .4270E-02, .4019E-02, .3646E-02, 
     8   .2759E-02, -.1484E-04, -.1408E-04, -.1446E-04, -.1446E-04, 
     9   -.1522E-04, -.1560E-04, -.1522E-04, -.1522E-04, -.1598E-04, 
     .   -.1484E-04, -.9895E-05, -.1028E-04, -.7612E-05, -.1903E-05, 
     1   .1903E-05, .0000E+00, .2283E-05, .6166E-05,  - .2740E-05/ 
      END 
 
      BLOCK DATA CKD12
      COMMON /BAND12/HK(5), COEO3(3,19,5), COEH2O(3,19)
      DATA HK/ 0.45, 0.30, 0.2, 0.04, 0.01/ 
      DATA (((COEO3(K,J,I),I=1,5),J=1,19),K=1,3)/ -.6590E+01, -.3912E+01
     1   , -.8513E+00, .2731E+01, .5515E+01, -.6157E+01, -.3583E+01, 
     2   -.7292E+00, .2740E+01, .5508E+01, -.5731E+01, -.3242E+01, 
     3   -.5800E+00, .2782E+01, .5485E+01, -.5301E+01, -.2901E+01, 
     4   -.4131E+00, .2805E+01, .5455E+01, -.4879E+01, -.2551E+01, 
     5   -.2288E+00, .2878E+01, .5416E+01, -.4449E+01, -.2201E+01, 
     6   -.2228E-01, .3000E+01, .5374E+01, -.4018E+01, -.1843E+01, 
     7   .2055E+00, .3143E+01, .5342E+01, -.3615E+01, -.1502E+01, 
     8   .4561E+00, .3288E+01, .5204E+01, -.3228E+01, -.1172E+01, 
     9   .7099E+00, .3396E+01, .5077E+01, -.2828E+01, -.8499E+00, 
     .   .9664E+00, .3463E+01, .4893E+01, -.2480E+01, -.5393E+00, 
     1   .1229E+01, .3493E+01, .4656E+01, -.2181E+01, -.2653E+00, 
     2   .1504E+01, .3456E+01, .4398E+01, -.1950E+01, -.1469E-01, 
     3   .1735E+01, .3387E+01, .4115E+01, -.1788E+01, .2517E+00, 
     4   .1919E+01, .3251E+01, .3832E+01, -.1677E+01, .5027E+00, 
     5   .2032E+01, .3088E+01, .3581E+01, -.1637E+01, .7373E+00, 
     6   .2100E+01, .2910E+01, .3364E+01, -.1650E+01, .9383E+00, 
     7   .2123E+01, .2793E+01, .3150E+01, -.1658E+01, .1091E+01, 
     8   .2112E+01, .2683E+01, .3021E+01, -.1654E+01, .1163E+01, 
     9   .2099E+01, .2602E+01, .2871E+01, .9498E-02, .8894E-02, 
     .   .1161E-01, .8828E-02, -.1669E-02, .9613E-02, .8347E-02, 
     1   .1053E-01, .8462E-02, -.1612E-02, .9700E-02, .7829E-02, 
     2   .9101E-02, .7915E-02, -.1439E-02, .9815E-02, .7167E-02, 
     3   .7981E-02, .7282E-02, -.1094E-02, .9671E-02, .6764E-02, 
     4   .6930E-02, .5613E-02, -.8347E-03, .9613E-02, .6312E-02, 
     5   .6225E-02, .4145E-02, -.1295E-02, .9728E-02, .6099E-02, 
     6   .5293E-02, .2965E-02, -.1756E-02, .9844E-02, .5915E-02, 
     7   .4496E-02, .1871E-02, -.2044E-02, .9930E-02, .5817E-02, 
     8   .3509E-02, .1324E-02, -.2044E-02, .9988E-02, .5535E-02, 
     9   .2711E-02, .6620E-03, -.1813E-02, .1034E-01, .5247E-02, 
     .   .1926E-02, -.2303E-03, -.1842E-02, .1058E-01, .4795E-02, 
     1   .1197E-02, -.9498E-03, -.2216E-02, .1084E-01, .4414E-02, 
     2   .6188E-03, -.1123E-02, -.2303E-02, .1079E-01, .3926E-02, 
     3   .1756E-03, -.1497E-02, -.2274E-02, .1039E-01, .3425E-02, 
     4   -.1900E-03, -.1353E-02, -.2389E-02, .9815E-02, .2769E-02, 
     5   -.6620E-03, -.1756E-02, -.1785E-02, .9818E-02, .2444E-02, 
     6   -.1016E-02, -.1410E-02, -.1698E-02, .1074E-01, .3218E-02, 
     7   -.1235E-02, -.1900E-02, -.2533E-02, .1145E-01, .3684E-02, 
     8   -.1364E-02, -.1353E-02, -.1957E-02, -.4030E-04, -.2375E-04, 
     9   -.3814E-05, -.4943E-04, -.3166E-04, -.3742E-04, -.1871E-04, 
     .   -.1137E-04, -.4317E-04, -.2878E-04, -.3526E-04, -.2015E-04, 
     1   -.1295E-04, -.4821E-04, -.2303E-04, -.3382E-04, -.2087E-04, 
     2   -.1519E-04, -.2231E-04, -.1871E-04, -.3454E-04, -.2087E-04, 
     3   -.8109E-05, -.6476E-05, -.1511E-04, -.3454E-04, -.1820E-04, 
     4   -.1269E-05, -.1439E-04, -.5037E-05, -.4173E-04, -.2598E-04, 
     5   .6645E-05, -.1943E-04, -.2087E-04, -.3454E-04, -.2267E-04, 
     6   .2159E-05, -.2231E-04, -.2159E-05, -.2950E-04, -.2080E-04, 
     7   .2159E-06, -.4317E-05, .1799E-04, -.3670E-04, -.1590E-04, 
     8   -.4461E-05, -.9354E-05, -.3598E-05, -.3216E-04, -.1475E-04, 
     9   -.2231E-05, -.1295E-04, -.2878E-05, -.3576E-04, -.7347E-05, 
     .   -.1022E-04, -.2159E-05, -.7915E-05, -.3015E-04, -.5230E-05, 
     1   -.5109E-05, -.6476E-05, -.7196E-05, -.2331E-04, -.1079E-04, 
     2   -.4102E-05, .1439E-05, -.1223E-04, -.2216E-04, -.1094E-04, 
     3   -.5325E-05, -.7196E-06, -.1655E-04, -.1036E-04, -.7627E-05, 
     4   -.2878E-05, .5037E-05, -.1295E-04, .1029E-04, -.1346E-04, 
     5   -.4821E-05, -.7915E-05, .7915E-05, .2835E-04, -.2893E-04, 
     6   -.1367E-05, -.7196E-05, -.1871E-04, .3965E-04, -.3310E-04, 
     7   -.3310E-05, -.7195E-06, .2303E-04/ 
      DATA ((COEH2O(K,J),J=1,19),K=1,3)/ -.1984E+02, -.1983E+02, 
     1   -.1982E+02, -.1981E+02, -.1963E+02, -.1917E+02, -.1871E+02, 
     2   -.1825E+02, -.1779E+02, -.1639E+02, -.1545E+02, -.1484E+02, 
     3   -.1433E+02, -.1387E+02, -.1345E+02, -.1305E+02, -.1268E+02, 
     4   -.1231E+02, -.1196E+02, .6071E-03, .2072E-02, .6196E-02, 
     5   .1030E-01, .1436E-01, .1846E-01, .2259E-01, .2667E-01, 
     6   .2993E-01, .2878E-01, .2803E-01, .2851E-01, .2864E-01, 
     7   .2874E-01, .2862E-01, .2859E-01, .2853E-01, .2868E-01, 
     8   .2887E-01, -.3808E-06, .2474E-04, .9895E-04, .1728E-03, 
     9   .1911E-03, .1165E-03, .4225E-04, -.3121E-04, -.8982E-04, 
     .   -.9553E-04, -.9705E-04, -.9591E-04, -.9287E-04, -.9172E-04, 
     1   -.9096E-04, -.9134E-04, -.9248E-04, -.1050E-03,  - .1031E-03/ 
      END 
 
      BLOCK DATA CKD13
      COMMON /BAND13/HK(2), COEH2O(3,19,2)
      DATA HK/ 0.95, 0.05/ 
      DATA (((COEH2O(K,J,I),I=1,2),J=1,19),K=1,3)/ -.1992E+02, 
     1   -.1446E+02, -.1992E+02, -.1405E+02, -.1991E+02, -.1363E+02, 
     2   -.1990E+02, -.1322E+02, -.1989E+02, -.1282E+02, -.1989E+02, 
     3   -.1242E+02, -.1988E+02, -.1201E+02, -.1987E+02, -.1159E+02, 
     4   -.1986E+02, -.1119E+02, -.1982E+02, -.1079E+02, -.1817E+02, 
     5   -.1039E+02, -.1659E+02, -.1000E+02, -.1537E+02, -.9623E+01, 
     6   -.1460E+02, -.9266E+01, -.1406E+02, -.8959E+01, -.1354E+02, 
     7   -.8676E+01, -.1309E+02, -.8411E+01, -.1267E+02, -.8232E+01, 
     8   -.1229E+02, -.8094E+01, .5024E-03, .3199E-01, .5652E-03, 
     9   .3199E-01, .6071E-03, .3211E-01, .6489E-03, .3199E-01, 
     .   .6699E-03, .3178E-01, .6908E-03, .3157E-01, .6908E-03, 
     1   .3109E-01, .6698E-03, .3075E-01, .6698E-03, .3054E-01, 
     2   .1474E-01, .3000E-01, .3085E-01, .2960E-01, .3659E-01, 
     3   .2935E-01, .3016E-01, .2920E-01, .2834E-01, .2895E-01, 
     4   .2780E-01, .2870E-01, .2753E-01, .2843E-01, .2755E-01, 
     5   .2820E-01, .2765E-01, .2732E-01, .2769E-01, .2705E-01, 
     6   .6299E-09, -.7993E-04, -.3802E-06, -.7992E-04, -.3802E-06, 
     7   -.8525E-04, -.3808E-06, -.8449E-04, -.7610E-06, -.7764E-04, 
     8   -.1142E-05, -.7231E-04, -.1142E-05, -.7345E-04, -.2284E-05, 
     9   -.8259E-04, -.2284E-05, -.8031E-04, .2436E-03, -.7878E-04, 
     .   .7612E-05, -.8525E-04, -.1248E-03, -.9439E-04, -.9477E-04, 
     1   -.9172E-04, -.8982E-04, -.8640E-04, -.7916E-04, -.6813E-04, 
     2   -.7574E-04, -.6090E-04, -.7612E-04, -.7117E-04, -.7498E-04, 
     3   -.7041E-04, -.7269E-04,  - .7992E-04/ 
      END 
 
      BLOCK DATA CKD14
      COMMON /BAND14/HK(10), COEHCA(3,19,10), COEHCB(3,19,10)
      DATA HK/ .3, .3, .2, .12, .06, .012, .004, .0025, .0011, .0004/ 
      DATA (((COEHCA(K,J,I),I=1,10),J=1,19),K=1,3)/-.1847E+02,-.1399E+02
     1,-.1106E+02,-.8539E+01,-.5852E+01,-.3295E+01,-.1208E+01,-.6272E-01
     2,.2055E+01,.6071E+01,-.1801E+02,-.1357E+02,-.1067E+02,-.8171E+01,
     3-.5562E+01,-.3071E+01,-.1073E+01,.1033E+00,.2055E+01,.6071E+01,
     4-.1755E+02,-.1314E+02,-.1027E+02,-.7798E+01,-.5224E+01,-.2823E+01,
     5-.9280E+00,.2723E+00,.2165E+01,.5969E+01,-.1709E+02,-.1272E+02,
     6-.9868E+01,-.7404E+01,-.4880E+01,-.2569E+01,-.6908E+00,.4453E+00,
     7.2241E+01,.5969E+01,-.1663E+02,-.1230E+02,-.9467E+01,-.7013E+01,
     8-.4535E+01,-.2297E+01,-.4408E+00,.6353E+00,.2359E+01,.5969E+01,
     9-.1617E+02,-.1188E+02,-.9050E+01,-.6619E+01,-.4160E+01,-.1967E+01,
     .-.1687E+00,.8213E+00,.2421E+01,.5969E+01,-.1571E+02,-.1147E+02,
     1-.8629E+01,-.6230E+01,-.3771E+01,-.1648E+01,.1573E+00,.1019E+01,
     2.2511E+01,.5884E+01,-.1525E+02,-.1106E+02,-.8215E+01,-.5841E+01,
     3-.3393E+01,-.1331E+01,.4013E+00,.1198E+01,.2654E+01,.5794E+01,
     4-.1480E+02,-.1066E+02,-.7800E+01,-.5454E+01,-.3032E+01,-.9870E+00,
     5.6323E+00,.1373E+01,.2905E+01,.5647E+01,-.1402E+02,-.9693E+01,
     6-.7206E+01,-.4846E+01,-.2656E+01,-.6540E+00,.8323E+00,.1530E+01,
     7.3211E+01,.5355E+01,-.1343E+02,-.9060E+01,-.6596E+01,-.4399E+01,
     8-.2294E+01,-.3519E+00,.9823E+00,.1673E+01,.3420E+01,.5083E+01,
     9-.1279E+02,-.8611E+01,-.5785E+01,-.4010E+01,-.1936E+01,-.1177E+00,
     ..1134E+01,.1974E+01,.3591E+01,.4770E+01,-.1230E+02,-.8174E+01,
     1-.5298E+01,-.3611E+01,-.1607E+01,.3636E-01,.1433E+01,.2260E+01,
     2.3539E+01,.4439E+01,-.1192E+02,-.7763E+01,-.4946E+01,-.3228E+01,
     3-.1321E+01,.1991E+00,.1720E+01,.2420E+01,.3383E+01,.4041E+01,
     4-.1154E+02,-.7377E+01,-.4576E+01,-.2851E+01,-.1093E+01,.4430E+00,
     5.1896E+01,.2462E+01,.3122E+01,.3620E+01,-.1118E+02,-.7003E+01,
     6-.4210E+01,-.2524E+01,-.8973E+00,.7490E+00,.1966E+01,.2363E+01,
     7.2818E+01,.3182E+01,-.1080E+02,-.6677E+01,-.3872E+01,-.2264E+01,
     8-.6846E+00,.9392E+00,.1867E+01,.2138E+01,.2505E+01,.2738E+01,
     9-.1031E+02,-.6353E+01,-.3596E+01,-.1938E+01,-.4537E+00,.1015E+01,
     ..1659E+01,.1830E+01,.2142E+01,.2287E+01,-.9695E+01,-.5977E+01,
     1-.3427E+01,-.1596E+01,-.1979E+00,.9458E+00,.1363E+01,.1545E+01,
     2.1743E+01,.1832E+01,.3628E-01,.2728E-01,.2213E-01,.1656E-01,
     3.1507E-01,.1564E-01,.1623E-01,.1419E-01,.1455E-01,.1089E-02,
     4.3632E-01,.2740E-01,.2164E-01,.1606E-01,.1369E-01,.1418E-01,
     5.1444E-01,.1275E-01,.1331E-01,.9210E-03,.3636E-01,.2746E-01,
     6.2114E-01,.1557E-01,.1239E-01,.1285E-01,.1237E-01,.1141E-01,
     7.1141E-01,.9210E-03,.3640E-01,.2748E-01,.2064E-01,.1516E-01,
     8.1141E-01,.1125E-01,.1092E-01,.1026E-01,.1011E-01,-.5652E-03,
     9.3646E-01,.2746E-01,.2024E-01,.1478E-01,.1036E-01,.9688E-02,
     ..9610E-02,.9305E-02,.9399E-02,-.6489E-03,.3651E-01,.2734E-01,
     1.1984E-01,.1438E-01,.9436E-02,.8486E-02,.8214E-02,.8995E-02,
     2.7892E-02,-.8582E-03,.3655E-01,.2723E-01,.1951E-01,.1402E-01,
     3.8716E-02,.7433E-02,.7169E-02,.8072E-02,.5443E-02,-.1172E-02,
     4.3659E-01,.2709E-01,.1911E-01,.1379E-01,.8107E-02,.6818E-02,
     5.6818E-02,.7033E-02,.3056E-02,-.1047E-02,.3670E-01,.2698E-01,
     6.1890E-01,.1363E-01,.7502E-02,.6371E-02,.6558E-02,.6489E-02,
     7-.5652E-03,-.1340E-02,.3592E-01,.2238E-01,.1804E-01,.1007E-01,
     8.6730E-02,.5512E-02,.6194E-02,.4375E-02,-.1109E-02,-.3559E-03,
     9.3609E-01,.2242E-01,.1526E-01,.8582E-02,.6284E-02,.5809E-02,
     ..4501E-02,.9420E-03,-.9001E-03,-.1005E-02,.3703E-01,.2196E-01,
     1.1281E-01,.7860E-02,.5861E-02,.5842E-02,.1800E-02,-.1591E-02,
     2-.1235E-02,-.9420E-03,.3728E-01,.2114E-01,.1347E-01,.6678E-02,
     3.5449E-02,.4837E-02,-.1084E-02,-.1361E-02,-.6699E-03,-.1256E-03,
     4.3683E-01,.2061E-01,.1350E-01,.6133E-02,.5449E-02,.2111E-02,
     5-.1386E-02,-.1235E-02,-.5652E-03,-.8373E-04,.3656E-01,.1988E-01,
     6.1348E-01,.5441E-02,.5149E-02,-.8813E-03,-.1116E-02,-.8373E-03,
     7-.3140E-03,-.6280E-04,.3669E-01,.1934E-01,.1363E-01,.5035E-02,
     8.3585E-02,-.1250E-02,-.9357E-03,-.8227E-03,-.3140E-03,-.4187E-04,
     9.3618E-01,.1856E-01,.1390E-01,.3836E-02,.1470E-02,-.1096E-02,
     .-.8080E-03,-.4480E-03,-.2093E-03,-.2093E-04,.3416E-01,.1741E-01,
     1.1431E-01,.1951E-02,-.2923E-04,-.9422E-03,-.4576E-03,-.2395E-03,
     2-.1565E-03,-.2799E-04,.3219E-01,.1674E-01,.1516E-01,.6652E-03,
     3-.5051E-03,-.7052E-03,-.2002E-03,-.2135E-03,-.7633E-04,-.7300E-04,
     4-.1290E-03,-.9934E-04,-.5595E-04,-.3996E-04,.1294E-04,-.9134E-05,
     5.1294E-05,-.3121E-05,-.4757E-04,-.1979E-04,-.1305E-03,-.9629E-04,
     6-.5481E-04,-.4301E-04,.1827E-04,-.9363E-05,.1777E-04,-.2185E-04,
     7-.1903E-04,-.1675E-04,-.1313E-03,-.9439E-04,-.5404E-04,-.4263E-04,
     8.9134E-05,-.1020E-04,.3524E-04,-.2599E-04,-.2093E-04,.1675E-04,
     9-.1313E-03,-.9172E-04,-.5252E-04,-.4567E-04,.4186E-05,-.3920E-05,
     ..2552E-04,-.2059E-04,-.2246E-04,-.1028E-04,-.1324E-03,-.9210E-04,
     1-.5138E-04,-.4491E-04,.6470E-05,-.2131E-05,.1496E-04,-.1572E-04,
     2-.3311E-04,-.8754E-05,-.1324E-03,-.9058E-04,-.5328E-04,-.4225E-04,
     3.1827E-05,-.8411E-06,.4719E-05,-.6813E-05,-.2474E-04,-.1256E-04,
     4-.1340E-03,-.8868E-04,-.5633E-04,-.4187E-04,-.4415E-05,.6055E-05,
     5-.1648E-04,-.1507E-04,.1979E-04,-.2131E-04,-.1340E-03,-.8373E-04,
     6-.5899E-04,-.3920E-04,-.4072E-05,.1491E-04,-.9781E-05,-.5328E-05,
     7.3578E-04,-.1979E-04,-.1321E-03,-.7954E-04,-.5899E-04,-.4072E-04,
     8.1066E-05,.5728E-05,-.5138E-05,-.8373E-05,.2626E-04,-.2436E-04,
     9-.1363E-03,-.6432E-04,-.5176E-04,-.3083E-04,.2169E-05,-.8944E-05,
     ..3159E-05,.6470E-05,-.4187E-05,.4948E-05,-.1302E-03,-.7802E-04,
     1-.3311E-04,-.1903E-04,.5328E-05,-.1884E-04,.1408E-04,.3311E-04,
     2.1142E-05,-.7613E-06,-.1473E-03,-.6737E-04,-.7536E-04,-.1085E-04,
     3-.1903E-05,-.1458E-04,.4034E-04,-.3941E-10,-.7992E-05,.2664E-05,
     4-.1361E-03,-.5709E-04,-.8550E-04,-.5709E-05,-.8640E-05,.6523E-05,
     5.1903E-05,-.8221E-05,-.3045E-05,-.9134E-05,-.1329E-03,-.5529E-04,
     6-.7107E-04,.2664E-05,-.9020E-05,.3320E-04,-.2131E-05,-.4187E-05,
     7-.7231E-05,-.3806E-05,-.1278E-03,-.5247E-04,-.6465E-04,.3806E-05,
     8-.6091E-05,.1245E-04,-.3844E-05,-.6090E-05,-.8754E-05,-.2664E-05,
     9-.1321E-03,-.5632E-04,-.5897E-04,.1012E-04,.1168E-04,-.4196E-06,
     .-.8411E-05,-.8868E-05,-.1484E-04,-.1522E-05,-.1252E-03,-.4907E-04,
     1-.5932E-04,.3245E-04,.1996E-04,-.3325E-05,-.5785E-05,-.6394E-05,
     2-.6851E-05,-.1142E-05,-.1093E-03,-.4731E-04,-.6761E-04,.1808E-04,
     3.1754E-04,-.5079E-05,-.5809E-05,-.5649E-05,-.3988E-05,-.5849E-06,
     4-.1151E-03,-.4965E-04,-.7163E-04,.7839E-05,.5505E-05,-.6084E-05,
     5-.3344E-05,-.3894E-05,-.1391E-05,-.1327E-05/
      DATA (((COEHCB(K,J,I),I=1,10),J=1,19),K=1,3)/-.9398E+01,-.5678E+01
     1,-.3606E+01,-.2192E+01,.2104E+01,.3044E+01,-.4587E+02,-.4587E+02,
     2-.4587E+02,-.4587E+02,-.9094E+01,-.5422E+01,-.3448E+01,-.1650E+01,
     3.2046E+01,.2749E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     4-.8760E+01,-.5270E+01,-.3329E+01,-.1147E+01,.2112E+01,.2709E+01,
     5-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.8537E+01,-.5152E+01,
     6-.3129E+01,-.9544E+00,.2254E+01,.2771E+01,-.4587E+02,-.4587E+02,
     7-.4587E+02,-.4587E+02,-.8176E+01,-.4936E+01,-.2680E+01,-.9259E+00,
     8.2247E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     9-.7836E+01,-.4676E+01,-.2378E+01,-.3550E+00,.1396E+01,.1976E+01,
     .-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.7419E+01,-.4122E+01,
     1-.2407E+01,-.1204E-01,.1744E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     2-.4587E+02,-.4587E+02,-.7124E+01,-.3727E+01,-.2160E+01,.6158E+00,
     3.1953E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     4-.6823E+01,-.3324E+01,-.1748E+01,-.9806E-01,.2319E+01,-.4587E+02,
     5-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.5957E+01,-.3017E+01,
     6-.1647E+01,.1398E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     7-.4587E+02,-.4587E+02,-.5115E+01,-.2290E+01,-.5273E+00,.5662E+00,
     8.1459E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     9-.4162E+01,-.1453E+01,.1116E+00,-.4587E+02,.9569E+00,-.4587E+02,
     .-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.3611E+01,-.9744E+00,
     1-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     2-.4587E+02,-.4587E+02,-.3075E+01,-.4176E+00,-.4587E+02,-.4587E+02,
     3-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     4-.3469E+01,-.9395E+00,.5092E+00,.6200E+00,-.4587E+02,-.4587E+02,
     5-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.3808E+01,-.1505E+01,
     6.3901E+00,.6264E+00,-.1155E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     7-.4587E+02,-.4587E+02,-.4058E+01,-.1818E+01,.2693E+00,.7087E+00,
     8.3820E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     9-.4262E+01,-.2097E+01,-.5711E-01,.5681E+00,.1310E+01,.7371E+00,
     .-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.3997E+01,-.1784E+01,
     1.4388E-01,.5167E+00,.6930E+00,-.6906E+00,-.4587E+02,-.4587E+02,
     2-.4587E+02,-.4587E+02,.2944E-01,.2723E-01,.1854E-01,.2023E-01,
     3.2254E-01,.3059E-02,.4788E+00,.3059E-02,.3059E-02,.3059E-02,
     4.3080E-01,.2549E-01,.1547E-01,.2225E-01,.2107E-01,.3059E-02,
     5.4737E+00,.3059E-02,.3059E-02,.3059E-02,.3269E-01,.2656E-01,
     6.2125E-01,.2179E-01,.2162E-01,.4589E+00,.4643E+00,.3059E-02,
     7.3059E-02,.3059E-02,.3322E-01,.2476E-01,.2075E-01,.2139E-01,
     8.1907E-01,.4501E+00,.4441E+00,.3059E-02,.3059E-02,.3059E-02,
     9.3387E-01,.2182E-01,.2665E-01,.1841E-01,.2506E-01,.3059E-02,
     ..3059E-02,.3059E-02,.3059E-02,.3059E-02,.3532E-01,.2091E-01,
     1.1995E-01,.2067E-01,.1949E-01,.4491E+00,.3059E-02,.3059E-02,
     2.3059E-02,.3059E-02,.3468E-01,.2075E-01,.2587E-01,.1401E-01,
     3.8646E-02,.3059E-02,.3059E-02,.3059E-02,.3059E-02,.3059E-02,
     4.3666E-01,.2430E-01,.1919E-01,.2007E-01,.3059E-02,.3059E-02,
     5.3059E-02,.3059E-02,.3059E-02,.3059E-02,.3613E-01,.2147E-01,
     6.1892E-01,.1361E-01,.3059E-02,.4506E+00,.3059E-02,.3059E-02,
     7.3059E-02,.3059E-02,.3129E-01,.1954E-01,.2442E-01,.1011E-01,
     8.4420E+00,.3059E-02,.3059E-02,.3059E-02,.3059E-02,.3059E-02,
     9.3177E-01,.2101E-01,.1526E-01,.4376E+00,.4379E+00,.3059E-02,
     ..3059E-02,.3059E-02,.3059E-02,.3059E-02,.2887E-01,.2044E-01,
     1.1285E-01,.3059E-02,-.4862E-03,.3059E-02,.3059E-02,.3059E-02,
     2.3059E-02,.3059E-02,.2759E-01,.2114E-01,.4303E+00,.3059E-02,
     3.3059E-02,.3059E-02,.3059E-02,.3059E-02,.3059E-02,.3059E-02,
     4.2880E-01,.1690E-01,-.4187E+00,.3059E-02,.3059E-02,.3059E-02,
     5.3059E-02,.3059E-02,.3059E-02,.3059E-02,.2852E-01,.2255E-01,
     6.2184E-01,.4334E+00,.4217E+00,.3059E-02,.3059E-02,.3059E-02,
     7.3059E-02,.3059E-02,.2840E-01,.2136E-01,.1644E-01,.2812E-01,
     8.4358E+00,.4288E+00,.3059E-02,.3059E-02,.3059E-02,.3059E-02,
     9.2809E-01,.2173E-01,.1708E-01,.3346E-01,.4225E-01,.4419E+00,
     ..3059E-02,.3059E-02,.3059E-02,.3059E-02,.2702E-01,.2260E-01,
     1.1607E-01,.2720E-01,.3982E-01,.4452E+00,.4365E+00,.4345E+00,
     2.4432E+00,.4623E+00,.2684E-01,.2328E-01,.2099E-01,.3040E-01,
     3.3867E-01,.4389E+00,.3132E-01,.3158E-01,.4083E-01,.4580E+00,
     4-.1581E-03,-.9707E-04,-.1250E-03,.2580E-03,.7378E-04,-.1617E-01,
     5.8646E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.1319E-03,-.9528E-04,
     6-.1710E-03,.7118E-04,.2076E-04,-.1608E-01,.8552E-02,-.4656E-05,
     7-.4656E-05,-.4656E-05,-.1721E-03,-.4680E-04,-.5522E-04,-.6242E-04,
     8.4517E-04,-.7777E-02,.8382E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     9-.1482E-03,-.4208E-04,-.5216E-04,-.6514E-04,-.8378E-04,-.7956E-02,
     ..8013E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.1501E-03,-.4002E-04,
     1-.1664E-03,.2272E-04,-.1888E-03,-.4656E-05,-.4656E-05,-.4656E-05,
     2-.4656E-05,-.4656E-05,-.1201E-03,-.4709E-04,-.5371E-04,-.1574E-03,
     3.1854E-03,-.7712E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     4-.1333E-03,-.1062E-03,.5785E-04,-.4150E-04,-.5717E-05,-.4656E-05,
     5-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1212E-03,-.8524E-04,
     6-.5895E-04,-.2884E-03,-.1581E-01,-.4656E-05,-.4656E-05,-.4656E-05,
     7-.4656E-05,-.4656E-05,-.8148E-04,-.9361E-04,-.2873E-03,.1883E-03,
     8-.1594E-01,.8133E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     9-.1221E-03,-.1430E-04,.6335E-04,-.2581E-03,.7977E-02,-.4656E-05,
     .-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.9257E-04,-.5008E-04,
     1.6389E-04,-.7455E-02,-.7745E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     2-.4656E-05,-.4656E-05,-.1186E-03,-.9037E-04,-.7461E-04,-.4656E-05,
     3.1168E-03,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     4-.8513E-04,-.5708E-04,.7763E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     5-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1124E-03,-.1228E-03,
     6.7663E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     7-.4656E-05,-.4656E-05,-.1015E-03,-.8369E-04,-.2167E-03,-.7548E-02,
     8.7608E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     9-.1049E-03,-.6414E-04,-.1384E-03,-.1644E-03,-.6919E-02,.7736E-02,
     .-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1008E-03,-.7047E-04,
     1-.1276E-03,-.2445E-03,-.1860E-03,.7975E-02,-.4656E-05,-.4656E-05,
     2-.4656E-05,-.4656E-05,-.9629E-04,-.1007E-03,-.1127E-03,-.1527E-03,
     3-.3238E-03,-.7373E-02,.7877E-02,.7840E-02,.7997E-02,.8345E-02,
     4-.8800E-04,-.1072E-03,-.1046E-03,-.1777E-03,-.2146E-03,-.7016E-02,
     5.1516E-01,.1532E-01,.1509E-01,.8268E-02/
      END 
 
      BLOCK DATA CKD15
      COMMON /BAND15/HK(12), COEHCA(3,19,12), COEHCB(3,19,12)
      DATA HK/.24,.36,.18,.1,.05,.02,.016,.012,.01,.006,.0039,.0021/
      DATA (((COEHCA(K,J,I),I=1,12),J=1,19),K=1,2)/ -.1921E+02, 
     1   -.1363E+02, -.1080E+02, -.8392E+01, -.6776E+01, -.5696E+01, 
     2   -.4572E+01, -.3752E+01, -.2382E+01, -.1110E+01, .6803E+00, 
     3   .3259E+01, -.1875E+02, -.1321E+02, -.1040E+02, -.8026E+01, 
     4   -.6449E+01, -.5401E+01, -.4316E+01, -.3498E+01, -.2141E+01, 
     5   -.9439E+00, .8103E+00, .3314E+01, -.1829E+02, -.1278E+02, 
     6   -.1000E+02, -.7646E+01, -.6089E+01, -.5085E+01, -.4047E+01, 
     7   -.3217E+01, -.1872E+01, -.7106E+00, .9573E+00, .3390E+01, 
     8   -.1783E+02, -.1236E+02, -.9596E+01, -.7264E+01, -.5735E+01, 
     9   -.4740E+01, -.3743E+01, -.2882E+01, -.1587E+01, -.4714E+00, 
     .   .1120E+01, .3425E+01, -.1737E+02, -.1195E+02, -.9193E+01, 
     1   -.6877E+01, -.5371E+01, -.4404E+01, -.3405E+01, -.2574E+01, 
     2   -.1298E+01, -.1747E+00, .1327E+01, .3547E+01, -.1691E+02, 
     3   -.1153E+02, -.8776E+01, -.6490E+01, -.4993E+01, -.4049E+01, 
     4   -.3039E+01, -.2256E+01, -.1012E+01, .1103E+00, .1530E+01, 
     5   .3651E+01, -.1644E+02, -.1112E+02, -.8360E+01, -.6105E+01, 
     6   -.4623E+01, -.3688E+01, -.2694E+01, -.1915E+01, -.6855E+00, 
     7   .3993E+00, .1714E+01, .3950E+01, -.1598E+02, -.1073E+02, 
     8   -.7943E+01, -.5723E+01, -.4236E+01, -.3314E+01, -.2338E+01, 
     9   -.1596E+01, -.3583E+00, .6963E+00, .1868E+01, .4127E+01, 
     .   -.1553E+02, -.1034E+02, -.7542E+01, -.5357E+01, -.3856E+01, 
     1   -.2942E+01, -.1986E+01, -.1299E+01, -.5472E-01, .9443E+00, 
     2   .2149E+01, .4261E+01, -.1485E+02, -.9661E+01, -.7008E+01, 
     3   -.4830E+01, -.3458E+01, -.2566E+01, -.1658E+01, -.9639E+00, 
     4   .2083E+00, .1182E+01, .2458E+01, .4452E+01, -.1427E+02, 
     5   -.9166E+01, -.6373E+01, -.4404E+01, -.3073E+01, -.2209E+01, 
     6   -.1349E+01, -.6648E+00, .4023E+00, .1452E+01, .2739E+01, 
     7   .4466E+01, -.1380E+02, -.8726E+01, -.5772E+01, -.3982E+01, 
     8   -.2732E+01, -.1874E+01, -.1052E+01, -.4403E+00, .5763E+00, 
     9   .1792E+01, .2999E+01, .4335E+01, -.1305E+02, -.8270E+01, 
     .   -.5304E+01, -.3586E+01, -.2392E+01, -.1568E+01, -.8299E+00, 
     1   -.2650E+00, .8584E+00, .2062E+01, .3141E+01, .4168E+01, 
     2   -.1269E+02, -.7900E+01, -.4956E+01, -.3205E+01, -.2065E+01, 
     3   -.1332E+01, -.6415E+00, -.7921E-01, .1170E+01, .2269E+01, 
     4   .3198E+01, .4066E+01, -.1227E+02, -.7536E+01, -.4576E+01, 
     5   -.2859E+01, -.1815E+01, -.1139E+01, -.4520E+00, .2272E+00, 
     6   .1371E+01, .2351E+01, .3150E+01, .3935E+01, -.1186E+02, 
     7   -.7159E+01, -.4223E+01, -.2538E+01, -.1619E+01, -.9324E+00, 
     8   -.1566E+00, .5151E+00, .1520E+01, .2339E+01, .3132E+01, 
     9   .3880E+01, -.1120E+02, -.6777E+01, -.3919E+01, -.2330E+01, 
     .   -.1387E+01, -.6737E+00, .1108E+00, .6991E+00, .1531E+01, 
     1   .2163E+01, .3150E+01, .3767E+01, -.9973E+01, -.6279E+01, 
     2   -.3638E+01, -.2048E+01, -.1098E+01, -.4407E+00, .3043E+00, 
     3   .7797E+00, .1424E+01, .2002E+01, .3122E+01, .3611E+01, 
     4   -.8483E+01, -.5607E+01, -.3357E+01, -.1744E+01, -.8884E+00, 
     5   -.2264E+00, .3800E+00, .7504E+00, .1245E+01, .2032E+01, 
     6   .3097E+01, .3546E+01, .3762E-01, .2372E-01, .1643E-01, 
     7   .1208E-01, .1170E-01, .1164E-01, .1214E-01, .1161E-01, 
     8   .1028E-01, .9185E-02, .7712E-02, .1001E-01, .3762E-01, 
     9   .2382E-01, .1593E-01, .1145E-01, .1059E-01, .1049E-01, 
     .   .1080E-01, .1057E-01, .8894E-02, .7807E-02, .7132E-02, 
     1   .1032E-01, .3764E-01, .2386E-01, .1555E-01, .1080E-01, 
     2   .9692E-02, .9231E-02, .9585E-02, .9644E-02, .7711E-02, 
     3   .6443E-02, .6223E-02, .9922E-02, .3764E-01, .2395E-01, 
     4   .1516E-01, .1028E-01, .8917E-02, .8415E-02, .8457E-02, 
     5   .8777E-02, .6436E-02, .5428E-02, .5499E-02, .8017E-02, 
     6   .3768E-01, .2399E-01, .1482E-01, .9692E-02, .8247E-02, 
     7   .7640E-02, .7582E-02, .7783E-02, .5432E-02, .4482E-02, 
     8   .4919E-02, .5903E-02, .3770E-01, .2401E-01, .1449E-01, 
     9   .9252E-02, .7620E-02, .6678E-02, .6845E-02, .6925E-02, 
     .   .4939E-02, .3471E-02, .4124E-02, .3873E-02, .3776E-01, 
     1   .2395E-01, .1419E-01, .8959E-02, .7096E-02, .6184E-02, 
     2   .6110E-02, .6075E-02, .4419E-02, .2891E-02, .3056E-02, 
     3   .1214E-02, .3780E-01, .2391E-01, .1392E-01, .8687E-02, 
     4   .6573E-02, .5733E-02, .5359E-02, .5009E-02, .4034E-02, 
     5   .2755E-02, .1968E-02, -.4187E-04, .3791E-01, .2382E-01, 
     6   .1373E-01, .8561E-02, .6060E-02, .5120E-02, .4618E-02, 
     7   .4713E-02, .3965E-02, .2481E-02, .8164E-03, -.1088E-02, 
     8   .3843E-01, .2148E-01, .1302E-01, .6384E-02, .5256E-02, 
     9   .4260E-02, .4077E-02, .4181E-02, .4132E-02, .2135E-02, 
     .   -.2931E-03, -.1151E-02, .3896E-01, .2081E-01, .1097E-01, 
     1   .5568E-02, .4475E-02, .3795E-02, .3828E-02, .3996E-02, 
     2   .3766E-02, .1193E-02, -.1089E-02, -.9420E-03, .3973E-01, 
     3   .2024E-01, .9943E-02, .4815E-02, .3820E-02, .3663E-02, 
     4   .3568E-02, .3881E-02, .2859E-02, .6698E-03, -.1549E-02, 
     5   -.6280E-03, .3635E-01, .1963E-01, .1061E-01, .3812E-02, 
     6   .3509E-02, .3429E-02, .3693E-02, .3316E-02, .1120E-02, 
     7   .6552E-03, -.1193E-02, -.1109E-02, .3631E-01, .1893E-01, 
     8   .1056E-01, .3172E-02, .3378E-02, .3164E-02, .2751E-02, 
     9   .1722E-02, .1112E-02, .4354E-03, -.7327E-03, -.1319E-02, 
     .   .3500E-01, .1828E-01, .1050E-01, .2831E-02, .2784E-02, 
     1   .2564E-02, .1469E-02, .7739E-03, .1209E-02, .7913E-03, 
     2   -.2512E-03, -.1758E-02, .3352E-01, .1763E-01, .1045E-01, 
     3   .2401E-02, .1928E-02, .1340E-02, .3753E-03, .5794E-03, 
     4   .9060E-03, .1042E-02, .1465E-03, -.2533E-02, .2880E-01, 
     5   .1729E-01, .1077E-01, .1347E-02, .1194E-02, -.1191E-03, 
     6   .2828E-03, .6606E-03, .9743E-03, .1002E-02, .0000E+00, 
     7   -.3140E-02, .2040E-01, .1585E-01, .1165E-01, .3871E-05, 
     8   .1509E-04, -.1046E-02, .2444E-03, .4359E-03, .1041E-02, 
     9   .2429E-02, -.1721E-03, -.2786E-02, .1737E-01, .1560E-01, 
     .   .1240E-01, -.2139E-03, -.1025E-02, -.1248E-02, -.6934E-04, 
     1   .1649E-03, .4062E-03, .1554E-02, -.4179E-03,  - .7795E-03/ 
      DATA (((COEHCA(K,J,I),I=1,12),J=1,19),K=3,3)/ -.1488E-03, 
     1   -.9248E-04, -.2322E-04, -.4187E-05, .1104E-04, .9895E-05, 
     2   -.2283E-05, .2512E-05, -.9058E-05, .8449E-05, .8297E-05, 
     3   -.3882E-04, -.1488E-03, -.9058E-04, -.2398E-04, -.5709E-05, 
     4   .1218E-04, .1180E-04, .1522E-05, .6927E-05, -.1161E-04, 
     5   .1714E-04, -.4948E-06, -.3540E-04, -.1500E-03, -.8830E-04, 
     6   -.2474E-04, -.8373E-05, .6470E-05, .7992E-05, .9096E-05, 
     7   .6737E-05, -.1485E-04, .1873E-04, -.4948E-06, -.4491E-04, 
     8   -.1500E-03, -.8601E-04, -.2664E-04, -.1028E-04, .6851E-05, 
     9   .6851E-05, .1294E-04, -.2550E-05, -.1520E-04, .2310E-04, 
     .   .4948E-06, -.2017E-04, -.1507E-03, -.8373E-04, -.2664E-04, 
     1   -.1256E-04, .4567E-05, .1028E-04, .9210E-05, -.2131E-05, 
     2   -.6995E-05, .7498E-05, -.1104E-04, -.2284E-05, -.1519E-03, 
     3   -.8183E-04, -.2816E-04, -.1142E-04, .7611E-06, .7231E-05, 
     4   .1751E-05, -.7612E-06, .8312E-05, .2436E-05, -.7231E-05, 
     5   .2398E-04, -.1530E-03, -.7992E-04, -.2893E-04, -.9896E-05, 
     6   .3806E-06, .8906E-05, .3159E-05, -.5328E-05, .3692E-05, 
     7   -.2093E-05, -.6851E-05, -.3045E-05, -.1538E-03, -.7536E-04, 
     8   -.3007E-04, -.8754E-05, -.3045E-05, .5138E-05, .9134E-06, 
     9   -.1979E-06, .1560E-05, -.1507E-04, .2284E-04, .9895E-05, 
     .   -.1541E-03, -.7688E-04, -.2969E-04, -.5709E-05, -.3996E-05, 
     1   .1142E-05, -.8373E-06, .1235E-04, -.7079E-05, -.6737E-05, 
     2   .1028E-04, .3578E-04, -.1560E-03, -.6851E-04, -.1903E-04, 
     3   -.4187E-05, -.4605E-05, -.1142E-06, .3878E-05, .3597E-05, 
     4   -.9591E-05, .5328E-05, .7612E-05, -.4948E-05, -.1587E-03, 
     5   -.6546E-04, -.2740E-04, -.7612E-06, -.3578E-05, .1713E-05, 
     6   .6064E-05, -.9781E-05, .1408E-05, .5709E-05, .8373E-05, 
     7   -.1256E-04, -.1484E-03, -.5823E-04, -.4301E-04, -.1522E-05, 
     8   .7498E-05, -.5328E-06, -.7855E-05, -.1599E-05, .1964E-04, 
     9   -.2284E-05, .7882E-10, .5328E-05, -.1238E-03, -.5700E-04, 
     .   -.5266E-04, .3286E-05, .4910E-05, -.8602E-05, .6090E-06, 
     1   .8454E-05, .1256E-05, -.4072E-05, -.1903E-05, .6470E-05, 
     2   -.1155E-03, -.5231E-04, -.4396E-04, .3626E-05, -.7051E-05, 
     3   -.1743E-05, .9667E-05, .2064E-04, -.2778E-05, -.6546E-05, 
     4   -.4948E-05, .1903E-05, -.1024E-03, -.5129E-04, -.4506E-04, 
     5   .7943E-06, .3074E-06, .3243E-05, .2754E-04, -.1479E-05, 
     6   .1661E-05, -.2969E-05, -.1066E-04, .7612E-06, -.8473E-04, 
     7   -.5418E-04, -.4674E-04, -.3418E-05, .9460E-05, .1151E-04, 
     8   .5714E-05, -.1069E-04, -.2022E-05, -.9061E-05, -.1104E-04, 
     9   -.3083E-04, -.4283E-04, -.5037E-04, -.4476E-04, .1951E-04, 
     .   .8922E-05, .1296E-04, -.4053E-05, -.4355E-05, -.2355E-05, 
     1   -.5004E-05, -.1218E-04, -.1522E-04, .6411E-05, -.5937E-04, 
     2   -.5331E-04, .1934E-04, .5284E-05, .1129E-04, -.2166E-05, 
     3   -.1484E-06, -.5407E-05, -.1364E-04, -.3115E-05, .3004E-04, 
     4   -.5074E-04, -.6256E-04, -.5097E-04, .2218E-04, .1228E-04, 
     5   -.1160E-05, -.1105E-05, .1618E-06, -.6089E-05, -.4216E-06, 
     6   -.5314E-05, .7903E-05/ 
      DATA (((COEHCB(K,J,I),I=1,12),J=1,19),K=1,2)/ -.9593E+01, 
     1   -.4078E+01, -.2812E+01, -.6506E+00, -.4123E+00, .2055E+01, 
     2   .4097E+01, .4671E+01, .4639E+01, -.4587E+02, -.4587E+02, 
     3   -.4587E+02, -.9276E+01, -.3757E+01, -.2467E+01, -.5784E+00, 
     4   .8833E-01, .2232E+01, .3826E+01, .4723E+01, .4942E+01, 
     5   .5135E+01, -.4587E+02, -.4587E+02, -.8968E+01, -.3508E+01, 
     6   -.2116E+01, -.1363E+00, .1662E+00, .2424E+01, .4220E+01, 
     7   .4513E+01, .1375E+01, .4601E+01, -.4587E+02, -.4587E+02, 
     8   -.8662E+01, -.3164E+01, -.1722E+01, .5178E-01, .7288E+00, 
     9   .2411E+01, .3805E+01, .4766E+01, .4342E+01, -.4587E+02, 
     .   -.4587E+02, -.4587E+02, -.8292E+01, -.2799E+01, -.1359E+01, 
     1   .3271E+00, .1650E+01, .2395E+01, .4192E+01, .4758E+01, 
     2   .2470E+01, -.4587E+02, -.4587E+02, -.4587E+02, -.7812E+01, 
     3   -.2404E+01, -.1085E+01, .7167E+00, .2202E+01, .2922E+01, 
     4   .4322E+01, .4591E+01, .4186E+01, -.4587E+02, -.4587E+02, 
     5   -.4587E+02, -.7441E+01, -.2066E+01, -.7142E+00, .1057E+01, 
     6   .2524E+01, .2946E+01, .4220E+01, .3607E+01, -.4587E+02, 
     7   -.4587E+02, -.4587E+02, -.4587E+02, -.7191E+01, -.1745E+01, 
     8   -.3487E+00, .1453E+01, .2739E+01, .3660E+01, .4114E+01, 
     9   .3245E+01, -.4587E+02, -.4587E+02, -.4587E+02, -.4587E+02, 
     .   -.6895E+01, -.1326E+01, -.3500E+00, .1647E+01, .2899E+01, 
     1   .4023E+01, .3361E+01, .3360E+01, -.4587E+02, -.4587E+02, 
     2   -.4587E+02, -.4587E+02, -.5876E+01, -.9573E+00, .2014E+00, 
     3   .2130E+01, .3493E+01, .4088E+01, -.4587E+02, -.4587E+02, 
     4   -.4587E+02, -.4587E+02, -.4587E+02, -.4587E+02, -.4429E+01, 
     5   -.3417E+00, .1204E+01, .2780E+01, .3843E+01, .3099E+01, 
     6   -.4587E+02, .3605E+01, -.4587E+02, -.4587E+02, -.4587E+02, 
     7   -.4587E+02, -.3122E+01, .2697E+00, .1866E+01, .3526E+01, 
     8   .3569E+01, .1025E+01, -.4587E+02, -.4587E+02, -.4587E+02, 
     9   -.4587E+02, -.4587E+02, -.4587E+02, -.2284E+01, .8186E+00, 
     .   .2754E+01, .3206E+01, .3704E+01, -.4587E+02, -.4587E+02, 
     1   .4625E+01, -.4587E+02, -.4587E+02, -.4587E+02, -.4587E+02, 
     2   -.1711E+01, .1220E+01, .3248E+01, -.4587E+02, .2565E+01, 
     3   .3297E+01, -.4587E+02, -.4587E+02, -.4587E+02, -.4587E+02, 
     4   -.4587E+02, -.4587E+02, -.1758E+01, .7970E+00, .2758E+01, 
     5   .2926E+01, .2613E+01, .1974E+01, -.4587E+02, .2310E+01, 
     6   -.4587E+02, -.4587E+02, -.4587E+02, -.4587E+02, -.1737E+01, 
     7   .3499E+00, .2246E+01, .2673E+01, .3308E+01, .3463E+01, 
     8   .3103E+01, .2611E+01, .2178E+01, -.4587E+02, -.4587E+02, 
     9   -.4587E+02, -.1559E+01, .2215E+00, .1875E+01, .2500E+01, 
     .   .3346E+01, .3585E+01, .3946E+01, .3533E+01, .3205E+01, 
     1   -.4587E+02, -.4587E+02, -.4587E+02, -.1601E+01, .5060E-01, 
     2   .1275E+01, .2176E+01, .3081E+01, .3649E+01, .3940E+01, 
     3   .4106E+01, .4112E+01, .4349E+01, .2292E+01, -.4587E+02, 
     4   -.1222E+01, .3199E+00, .1642E+01, .2380E+01, .3254E+01, 
     5   .3534E+01, .3687E+01, .3717E+01, .3402E+01, .3868E+01, 
     6   -.4587E+02, -.4587E+02, .2967E-01, .1697E-01, .1795E-01, 
     7   .1387E-01, .2032E-01, .1187E-01, .2560E-01, .1044E-01, 
     8   -.4560E+00, .3059E-02, .3059E-02, .3059E-02, .2998E-01, 
     9   .1586E-01, .1786E-01, .1521E-01, .1710E-01, .1061E-01, 
     .   .2030E-01, .1158E-01, .4452E+00, .3059E-02, .3059E-02, 
     1   .3059E-02, .2993E-01, .1551E-01, .1481E-01, .9846E-02, 
     2   .2443E-01, .1150E-01, .1865E-01, .1376E-01, .4617E+00, 
     3   .3059E-02, .3059E-02, .3059E-02, .3035E-01, .1417E-01, 
     4   .1438E-01, .1511E-01, .1901E-01, .8582E-02, .1746E-01, 
     5   .1450E-01, .4523E+00, .3059E-02, .3059E-02, .3059E-02, 
     6   .2970E-01, .1347E-01, .1322E-01, .1252E-01, .1665E-01, 
     7   .1037E-01, .1320E-01, .1199E-01, .4436E+00, .3059E-02, 
     8   .3059E-02, .3059E-02, .2949E-01, .1291E-01, .1671E-01, 
     9   .1111E-01, .1400E-01, .1318E-01, .1060E-01, .1046E-01, 
     .   .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3004E-01, 
     1   .1300E-01, .1413E-01, .9085E-02, .9764E-02, .2260E-01, 
     2   .9778E-02, .4671E+00, .3059E-02, .3059E-02, .3059E-02, 
     3   .3059E-02, .3086E-01, .1436E-01, .1205E-01, .1081E-01, 
     4   .4681E-02, .1479E-01, .1888E-01, .3494E-01, .3059E-02, 
     5   .3059E-02, .3059E-02, .3059E-02, .3094E-01, .1500E-01, 
     6   .1457E-01, .1060E-01, .8319E-02, .8983E-02, .3791E-01, 
     7   .2232E-01, .4631E+00, .3059E-02, .3059E-02, .3059E-02, 
     8   .3158E-01, .1585E-01, .1292E-01, .6531E-02, .1383E-01, 
     9   .4605E+00, .4662E+00, .3059E-02, .3059E-02, .3059E-02, 
     .   .3059E-02, .3059E-02, .3182E-01, .1586E-01, .8724E-02, 
     1   .5798E-02, .2454E-01, .4607E+00, .4560E+00, .4511E+00, 
     2   .3059E-02, .3059E-02, .3059E-02, .3059E-02, .2369E-01, 
     3   .1606E-01, .5477E-02, .1228E-01, .4579E+00, .4561E+00, 
     4   .4497E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, 
     5   .3059E-02, .2190E-01, .1779E-01, .6267E-02, .4535E+00, 
     6   .4533E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, 
     7   .3059E-02, .3059E-02, .3059E-02, .2100E-01, .1653E-01, 
     8   .7449E-02, .4543E+00, .4472E+00, .4439E+00, .3059E-02, 
     9   .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, 
     .   .1864E-01, .1771E-01, .7040E-02, .2877E-01, .3381E-01, 
     1   .2691E-01, .4466E+00, .3059E-02, .4613E+00, .3059E-02, 
     2   .3059E-02, .3059E-02, .1637E-01, .1641E-01, .8424E-02, 
     3   .1318E-01, .2060E-01, .3426E-01, .4122E-01, .4621E+00, 
     4   .4555E+00, .4525E+00, .3059E-02, .3059E-02, .1607E-01, 
     5   .1452E-01, .8013E-02, .1213E-01, .1482E-01, .2125E-01, 
     6   .3379E-01, .3562E-01, .4619E+00, .4569E+00, .3059E-02, 
     7   .3059E-02, .1698E-01, .1538E-01, .6616E-02, .1147E-01, 
     8   .1217E-01, .1696E-01, .1871E-01, .2273E-01, .4513E-01, 
     9   .4702E+00, .4617E+00, .4553E+00, .1700E-01, .1547E-01, 
     .   .6456E-02, .1324E-01, .1502E-01, .2095E-01, .2547E-01, 
     1   .2823E-01, .4107E-01, .4676E+00, .4583E+00, .4498E+00/ 
      DATA (((COEHCB(K,J,I),I=1,12),J=1,19),K=3,3)/ -.6747E-05, 
     1   -.2483E-04, .6575E-04, .1026E-03, .3888E-03, -.8519E-04, 
     2   -.1629E-03, -.1808E-04, -.8355E-02, -.4656E-05, -.4656E-05, 
     3   -.4656E-05, -.2270E-04, -.3427E-04, .5118E-04, .1218E-03, 
     4   .1245E-03, -.1245E-03, .3841E-05, -.4151E-04, -.8763E-02, 
     5   -.1687E-01, -.4656E-05, -.4656E-05, -.4557E-04, -.3023E-04, 
     6   .2286E-04, .5656E-04, .4113E-04, -.1407E-03, -.1301E-03, 
     7   .8503E-04, -.7284E-02, -.1669E-01, -.4656E-05, -.4656E-05, 
     8   -.5325E-04, -.5309E-04, -.1246E-04, .2244E-04, .5136E-04, 
     9   -.1272E-03, .4217E-04, -.1749E-04, -.8435E-02, -.4656E-05, 
     .   -.4656E-05, -.4656E-05, -.6857E-04, -.7217E-04, .1740E-05, 
     1   .3653E-04, -.1490E-03, -.4090E-04, -.2376E-04, .2047E-04, 
     2   -.7974E-02, -.4656E-05, -.4656E-05, -.4656E-05, -.1232E-03, 
     3   -.9826E-04, -.2849E-04, .1703E-04, -.1895E-03, -.3363E-03, 
     4   .7102E-04, -.1838E-05, -.1655E-01, -.4656E-05, -.4656E-05, 
     5   -.4656E-05, -.9896E-04, -.5127E-04, -.2704E-04, -.1218E-04, 
     6   -.1207E-03, -.5883E-04, .6893E-04, -.7924E-02, -.4656E-05, 
     7   -.4656E-05, -.4656E-05, -.4656E-05, -.7837E-04, -.4980E-04, 
     8   .6902E-05, -.1072E-03, -.4051E-04, -.1991E-05, -.1173E-03, 
     9   -.5195E-04, -.4656E-05, -.4656E-05, -.4656E-05, -.4656E-05, 
     .   -.8136E-04, -.8102E-04, .1254E-03, -.4658E-04, .3173E-04, 
     1   -.4461E-05, -.1558E-03, -.2036E-03, .8360E-02, -.4656E-05, 
     2   -.4656E-05, -.4656E-05, -.2232E-04, -.6411E-04, .9486E-04, 
     3   -.2322E-03, -.8282E-04, -.8202E-02, .8416E-02, -.4656E-05, 
     4   -.4656E-05, -.4656E-05, -.4656E-05, -.4656E-05, -.1398E-03, 
     5   -.7165E-04, -.4258E-04, -.3970E-04, -.2839E-03, -.7873E-02, 
     6   .8231E-02, -.8213E-02, -.4656E-05, -.4656E-05, -.4656E-05, 
     7   -.4656E-05, -.6754E-04, -.7469E-04, -.6898E-04, -.1702E-03, 
     8   -.8079E-02, -.7270E-02, .8116E-02, -.4656E-05, -.4656E-05, 
     9   -.4656E-05, -.4656E-05, -.4656E-05, -.2396E-04, -.2361E-04, 
     .   -.8664E-04, -.8038E-02, -.8207E-02, -.4656E-05, -.4656E-05, 
     1   -.1670E-01, -.4656E-05, -.4656E-05, -.4656E-05, -.4656E-05, 
     2   -.5479E-04, -.7593E-04, -.1005E-03, .8199E-02, -.7942E-02, 
     3   -.8244E-02, -.4656E-05, -.4656E-05, -.4656E-05, -.4656E-05, 
     4   -.4656E-05, -.4656E-05, -.3806E-04, -.5825E-04, -.1003E-03, 
     5   -.2925E-03, -.1506E-03, .3148E-04, .8060E-02, -.1593E-01, 
     6   .8327E-02, -.4656E-05, -.4656E-05, -.4656E-05, -.4706E-04, 
     7   -.3630E-04, -.7811E-04, -.6881E-04, -.1822E-03, -.3091E-03, 
     8   -.3033E-03, -.7684E-02, -.7663E-02, .8167E-02, -.4656E-05, 
     9   -.4656E-05, -.7669E-04, -.4610E-04, -.8063E-04, -.7250E-04, 
     .   -.1094E-03, -.1241E-03, -.2944E-03, -.1736E-03, -.7886E-02, 
     1   .8248E-02, -.4656E-05, -.4656E-05, -.7138E-04, -.4545E-04, 
     2   -.3653E-04, -.6075E-04, -.4528E-04, -.1077E-03, -.1119E-03, 
     3   -.1657E-03, -.4695E-03, -.8112E-02, -.7587E-02, .8217E-02, 
     4   -.6812E-04, -.4558E-04, -.6739E-04, -.8861E-04, -.9386E-04, 
     5   -.1334E-03, -.2007E-03, -.2179E-03, -.1650E-03, -.8001E-02, 
     6   .8273E-02, .8118E-02/ 
      END 
 
      BLOCK DATA CKD16
      COMMON /BAND16/HK(7), COEH2O(3,19,7)
      DATA HK/ .12, .24, .24, .20, .12, .06, .02/ 
      DATA (((COEH2O(K,J,I),I=1,7),J=1,19),K=1,3)/ -.2344E+02, 
     1   -.2016E+02, -.1986E+02, -.1655E+02, -.1243E+02, -.8437E+01, 
     2   -.4858E+01, -.2298E+02, -.2014E+02, -.1984E+02, -.1609E+02, 
     3   -.1198E+02, -.8020E+01, -.4548E+01, -.2252E+02, -.2012E+02, 
     4   -.1981E+02, -.1564E+02, -.1153E+02, -.7596E+01, -.4239E+01, 
     5   -.2206E+02, -.2009E+02, -.1957E+02, -.1517E+02, -.1111E+02, 
     6   -.7161E+01, -.3871E+01, -.2160E+02, -.2007E+02, -.1911E+02, 
     7   -.1472E+02, -.1065E+02, -.6721E+01, -.3479E+01, -.2113E+02, 
     8   -.2005E+02, -.1865E+02, -.1426E+02, -.1021E+02, -.6302E+01, 
     9   -.3081E+01, -.2067E+02, -.2003E+02, -.1819E+02, -.1379E+02, 
     .   -.9765E+01, -.5883E+01, -.2678E+01, -.2026E+02, -.2001E+02, 
     1   -.1773E+02, -.1333E+02, -.9332E+01, -.5443E+01, -.2253E+01, 
     2   -.2024E+02, -.1999E+02, -.1727E+02, -.1288E+02, -.8897E+01, 
     3   -.5029E+01, -.1858E+01, -.2026E+02, -.1959E+02, -.1481E+02, 
     4   -.1147E+02, -.7477E+01, -.4555E+01, -.1464E+01, -.2022E+02, 
     5   -.1632E+02, -.1305E+02, -.9885E+01, -.6689E+01, -.4108E+01, 
     6   -.1068E+01, -.1936E+02, -.1438E+02, -.1163E+02, -.8499E+01, 
     7   -.6146E+01, -.3673E+01, -.6816E+00, -.1675E+02, -.1281E+02, 
     8   -.1020E+02, -.7716E+01, -.5678E+01, -.3256E+01, -.3125E+00, 
     9   -.1510E+02, -.1124E+02, -.8821E+01, -.7140E+01, -.5243E+01, 
     .   -.2851E+01, -.2560E-01, -.1334E+02, -.9708E+01, -.8061E+01, 
     1   -.6611E+01, -.4842E+01, -.2459E+01, .1711E+00, -.1155E+02, 
     2   -.8798E+01, -.7440E+01, -.6123E+01, -.4439E+01, -.2089E+01, 
     3   .2480E+00, -.1020E+02, -.8154E+01, -.6945E+01, -.5681E+01, 
     4   -.4055E+01, -.1737E+01, .2390E+00, -.9464E+01, -.7677E+01, 
     5   -.6512E+01, -.5284E+01, -.3707E+01, -.1453E+01, .2015E+00, 
     6   -.9033E+01, -.7246E+01, -.6093E+01, -.4882E+01, -.3346E+01, 
     7   -.1264E+01, .1033E+00, .4658E-01, .5840E-02, .4626E-02, 
     8   .2688E-01, .2395E-01, .1804E-01, .2074E-01, .4660E-01, 
     9   .1884E-02, .8561E-02, .2690E-01, .2403E-01, .1788E-01, 
     .   .1934E-01, .4660E-01, .1800E-02, .1252E-01, .2694E-01, 
     1   .2393E-01, .1786E-01, .1825E-01, .4660E-01, .1779E-02, 
     2   .1649E-01, .2696E-01, .2397E-01, .1779E-01, .1765E-01, 
     3   .4348E-01, .1758E-02, .2043E-01, .2696E-01, .2393E-01, 
     4   .1748E-01, .1675E-01, .3944E-01, .1737E-02, .2445E-01, 
     5   .2698E-01, .2384E-01, .1752E-01, .1549E-01, .3538E-01, 
     6   .1654E-02, .2847E-01, .2702E-01, .2384E-01, .1714E-01, 
     7   .1565E-01, .3127E-01, .1570E-02, .3245E-01, .2705E-01, 
     8   .2374E-01, .1712E-01, .1514E-01, .2715E-01, .1444E-02, 
     9   .3540E-01, .2711E-01, .2363E-01, .1702E-01, .1446E-01, 
     .   .2960E-01, .1760E-01, .2977E-01, .2397E-01, .2087E-01, 
     1   .1618E-01, .1445E-01, .2466E-01, .3039E-01, .2428E-01, 
     2   .2217E-01, .1821E-01, .1593E-01, .1463E-01, .2640E-01, 
     3   .2545E-01, .2231E-01, .2060E-01, .1773E-01, .1555E-01, 
     4   .1473E-01, .3456E-01, .2135E-01, .2030E-01, .1844E-01, 
     5   .1740E-01, .1559E-01, .1428E-01, .3203E-01, .2047E-01, 
     6   .1809E-01, .1760E-01, .1725E-01, .1545E-01, .1541E-01, 
     7   .2137E-01, .1857E-01, .1616E-01, .1698E-01, .1700E-01, 
     8   .1537E-01, .1636E-01, .1338E-01, .1518E-01, .1580E-01, 
     9   .1658E-01, .1710E-01, .1518E-01, .1513E-01, .1570E-01, 
     .   .1614E-01, .1603E-01, .1673E-01, .1706E-01, .1497E-01, 
     1   .1439E-01, .1987E-01, .1731E-01, .1601E-01, .1675E-01, 
     2   .1681E-01, .1535E-01, .1425E-01, .2018E-01, .1723E-01, 
     3   .1597E-01, .1691E-01, .1666E-01, .1509E-01, .1446E-01, 
     4   -.2873E-03, -.8031E-04, .4225E-04, -.9287E-04, -.6013E-04, 
     5   -.4339E-04, -.2474E-04, -.2862E-03, -.8372E-05, .1146E-03, 
     6   -.9248E-04, -.6166E-04, -.3882E-04, -.1827E-04, -.2870E-03, 
     7   -.6851E-05, .1865E-03, -.9172E-04, -.6128E-04, -.3616E-04, 
     8   -.7612E-05, -.2877E-03, -.7231E-05, .1880E-03, -.9287E-04, 
     9   -.5671E-04, -.4110E-04, -.1104E-04, -.3429E-03, -.7612E-05, 
     .   .1149E-03, -.9287E-04, -.6356E-04, -.4529E-04, -.2436E-04, 
     1   -.4187E-03, -.7992E-05, .4339E-04, -.9325E-04, -.6280E-04, 
     2   -.4225E-04, -.3197E-04, -.4925E-03, -.8754E-05, -.2740E-04, 
     3   -.9477E-04, -.6432E-04, -.3768E-04, -.3361E-04, -.5511E-03, 
     4   -.8753E-05, -.9972E-04, -.9515E-04, -.6394E-04, -.3806E-04, 
     5   -.3787E-04, -.4792E-03, -.1028E-04, -.1534E-03, -.9477E-04, 
     6   -.6356E-04, -.3616E-04, -.2923E-04, -.5070E-03, .1922E-03, 
     7   -.1028E-03, -.5823E-04, -.7954E-04, -.2550E-04, -.3893E-04, 
     8   -.3776E-03, -.1043E-03, -.7993E-04, -.7422E-04, -.4948E-04, 
     9   -.3007E-04, -.3863E-04, .8335E-04, -.5709E-04, -.6090E-04, 
     .   -.7840E-04, -.3692E-04, -.3007E-04, -.4251E-04, -.6204E-04, 
     1   -.4872E-04, -.3806E-04, -.4681E-04, -.3463E-04, -.3007E-04, 
     2   -.4312E-04, -.1142E-04, -.5176E-04, -.5024E-04, -.3007E-04, 
     3   -.3730E-04, -.3037E-04, -.3888E-04, .2550E-04, -.6508E-04, 
     4   -.2512E-04, -.3083E-04, -.3197E-04, -.3041E-04, -.3750E-04, 
     5   .1484E-04, -.1941E-04, -.2626E-04, -.3349E-04, -.3463E-04, 
     6   -.2896E-04, -.1716E-04, -.7231E-04, -.3920E-04, -.2893E-04, 
     7   -.3540E-04, -.3311E-04, -.3734E-04, -.2550E-05, -.7650E-04, 
     8   -.3159E-04, -.2778E-04, -.3121E-04, -.2169E-04, -.4365E-04, 
     9   -.1546E-04, -.7916E-04, -.2931E-04, -.2854E-04, -.3654E-04, 
     .   -.1979E-04, -.4811E-04,  - .1435E-04/ 
      END 
 
      BLOCK DATA CKD17
      COMMON /BAND17/HK(7), COEH2O(3,19,7)
      DATA HK/ .12, .26, .22, .20, .10, .085, .015/ 
      DATA (((COEH2O(K,J,I),I=1,7),J=1,19),K=1,3)/ -.2255E+02, 
     1   -.2000E+02, -.1703E+02, -.1282E+02, -.9215E+01, -.5938E+01, 
     2   -.2009E+01, -.2209E+02, -.1997E+02, -.1657E+02, -.1236E+02, 
     3   -.8764E+01, -.5499E+01, -.1582E+01, -.2163E+02, -.1993E+02, 
     4   -.1611E+02, -.1191E+02, -.8324E+01, -.5061E+01, -.1170E+01, 
     5   -.2117E+02, -.1990E+02, -.1565E+02, -.1146E+02, -.7889E+01, 
     6   -.4631E+01, -.7737E+00, -.2071E+02, -.1987E+02, -.1519E+02, 
     7   -.1100E+02, -.7440E+01, -.4179E+01, -.3719E+00, -.2026E+02, 
     8   -.1985E+02, -.1473E+02, -.1054E+02, -.6995E+01, -.3721E+01, 
     9   .0000E+00, -.2024E+02, -.1982E+02, -.1426E+02, -.1009E+02, 
     .   -.6549E+01, -.3284E+01, .4053E+00, -.2022E+02, -.1980E+02, 
     1   -.1381E+02, -.9639E+01, -.6097E+01, -.2821E+01, .8375E+00, 
     2   -.2021E+02, -.1933E+02, -.1335E+02, -.9187E+01, -.5653E+01, 
     3   -.2379E+01, .1272E+01, -.2010E+02, -.1503E+02, -.1125E+02, 
     4   -.7665E+01, -.4492E+01, -.1893E+01, .1642E+01, -.1747E+02, 
     5   -.1278E+02, -.9547E+01, -.6120E+01, -.3756E+01, -.1443E+01, 
     6   .1995E+01, -.1529E+02, -.1095E+02, -.8107E+01, -.5036E+01, 
     7   -.3182E+01, -.1032E+01, .2429E+01, -.1370E+02, -.9303E+01, 
     8   -.6691E+01, -.4357E+01, -.2683E+01, -.6173E+00, .2805E+01, 
     9   -.1150E+02, -.7859E+01, -.5618E+01, -.3843E+01, -.2234E+01, 
     .   -.2171E+00, .2973E+01, -.9590E+01, -.6537E+01, -.4886E+01, 
     1   -.3355E+01, -.1805E+01, .1615E+00, .3157E+01, -.7530E+01, 
     2   -.5699E+01, -.4306E+01, -.2892E+01, -.1388E+01, .5448E+00, 
     3   .3155E+01, -.6758E+01, -.5112E+01, -.3809E+01, -.2464E+01, 
     4   -.9947E+00, .8713E+00, .3203E+01, -.6245E+01, -.4610E+01, 
     5   -.3376E+01, -.2058E+01, -.6166E+00, .1073E+01, .3109E+01, 
     6   -.5777E+01, -.4175E+01, -.2963E+01, -.1671E+01, -.2556E+00, 
     7   .1241E+01, .3014E+01, .4264E-01, .1968E-02, .1863E-01, 
     8   .1436E-01, .1101E-01, .1055E-01, .1281E-01, .4264E-01, 
     9   .1989E-02, .1861E-01, .1438E-01, .1095E-01, .1030E-01, 
     .   .1211E-01, .3996E-01, .1968E-02, .1861E-01, .1434E-01, 
     1   .1103E-01, .1019E-01, .1160E-01, .3600E-01, .1947E-02, 
     2   .1861E-01, .1442E-01, .1086E-01, .1003E-01, .1157E-01, 
     3   .3203E-01, .5756E-02, .1861E-01, .1444E-01, .1080E-01, 
     4   .9922E-02, .1151E-01, .2801E-01, .9713E-02, .1859E-01, 
     5   .1446E-01, .1070E-01, .9880E-02, .1066E-01, .2393E-01, 
     6   .1369E-01, .1859E-01, .1451E-01, .1057E-01, .9880E-02, 
     7   .1072E-01, .1987E-01, .1767E-01, .1863E-01, .1451E-01, 
     8   .1040E-01, .9880E-02, .1057E-01, .1572E-01, .2169E-01, 
     9   .1863E-01, .1442E-01, .1022E-01, .9742E-02, .1036E-01, 
     .   .3391E-02, .1884E-01, .1566E-01, .1105E-01, .1011E-01, 
     1   .1001E-01, .1017E-01, .1982E-01, .1444E-01, .1189E-01, 
     2   .1030E-01, .9859E-02, .9861E-02, .1038E-01, .1748E-01, 
     3   .1321E-01, .9922E-02, .1068E-01, .1013E-01, .9937E-02, 
     4   .9958E-02, .1346E-01, .9943E-02, .9566E-02, .1097E-01, 
     5   .9815E-02, .9964E-02, .1059E-01, .9817E-02, .7159E-02, 
     6   .8687E-02, .1114E-01, .1007E-01, .1014E-01, .1058E-01, 
     7   .3370E-02, .7264E-02, .9378E-02, .1112E-01, .9767E-02, 
     8   .1016E-01, .1101E-01, .2993E-02, .8017E-02, .9566E-02, 
     9   .1116E-01, .9738E-02, .1025E-01, .1086E-01, .8331E-02, 
     .   .8771E-02, .1001E-01, .1117E-01, .9847E-02, .1076E-01, 
     1   .1084E-01, .7850E-02, .9378E-02, .1001E-01, .1105E-01, 
     2   .9964E-02, .1113E-01, .1168E-01, .8038E-02, .9336E-02, 
     3   .9817E-02, .1096E-01, .1024E-01, .1175E-01, .1107E-01, 
     4   -.2188E-03, -.2283E-05, -.8069E-04, -.4415E-04, -.2284E-04, 
     5   -.4491E-04, -.4518E-04, -.2196E-03, -.2665E-05, -.8107E-04, 
     6   -.4301E-04, -.2398E-04, -.4795E-04, -.4693E-04, -.2683E-03, 
     7   -.3045E-05, -.8107E-04, -.4301E-04, -.2246E-04, -.4757E-04, 
     8   -.4152E-04, -.3403E-03, -.4187E-05, -.8031E-04, -.3996E-04, 
     9   -.1865E-04, -.4301E-04, -.4350E-04, -.4118E-03, .6584E-04, 
     .   -.8107E-04, -.4034E-04, -.1903E-04, -.4643E-04, -.4834E-04, 
     1   -.4803E-03, .1378E-03, -.8069E-04, -.4072E-04, -.1713E-04, 
     2   -.5176E-04, -.3460E-04, -.4099E-03, .2101E-03, -.8069E-04, 
     3   -.3920E-04, -.1713E-04, -.5024E-04, -.3524E-04, -.3391E-03, 
     4   .2809E-03, -.7992E-04, -.3616E-04, -.2017E-04, -.5633E-04, 
     5   -.4886E-04, -.2668E-03, .2078E-03, -.8069E-04, -.3768E-04, 
     6   -.2131E-04, -.5580E-04, -.5454E-04, -.2207E-04, -.8601E-04, 
     7   -.4643E-04, -.2436E-04, -.4148E-04, -.5458E-04, -.4579E-04, 
     8   -.5138E-04, -.2893E-04, -.3273E-04, -.3882E-04, -.3920E-04, 
     9   -.5035E-04, -.3170E-04, -.2169E-04, -.3007E-04, -.2740E-04, 
     .   -.5328E-04, -.4491E-04, -.4403E-04, -.6383E-04, .4834E-04, 
     1   -.2702E-04, -.4453E-04, -.4339E-04, -.4457E-04, -.4551E-04, 
     2   -.8133E-04, .3768E-04, -.7611E-06, -.2626E-04, -.4643E-04, 
     3   -.4305E-04, -.4840E-04, -.5149E-04, .7193E-04, -.2169E-04, 
     4   -.4491E-04, -.3996E-04, -.4483E-04, -.4487E-04, -.6698E-04, 
     5   -.4834E-04, -.3463E-04, -.4986E-04, -.4377E-04, -.4514E-04, 
     6   -.5377E-04, -.2626E-04, -.4187E-04, -.3692E-04, -.5100E-04, 
     7   -.4651E-04, -.4392E-04, -.5386E-04, -.4643E-04, -.4301E-04, 
     8   -.3578E-04, -.5176E-04, -.4594E-04, -.4551E-04, -.3920E-04, 
     9   -.3425E-04, -.4491E-04, -.3654E-04, -.5138E-04, -.4377E-04, 
     .   -.5614E-04, -.5758E-04,  - .3600E-04/ 
      END 
 
      BLOCK DATA CKD18
      COMMON /BAND18/HK(8), COEH2O(3,19,8)
      DATA HK/ .10, .1, .2, .25, .2, .1, .03, .02/ 
      DATA (((COEH2O(K,J,I),I=1,8),J=1,19),K=1,3)/ -.2121E+02, 
     1   -.2002E+02, -.1676E+02, -.1274E+02, -.8780E+01, -.5167E+01, 
     2   -.2692E+01, -.6275E+00, -.2075E+02, -.1996E+02, -.1630E+02, 
     3   -.1228E+02, -.8324E+01, -.4718E+01, -.2260E+01, -.2303E+00, 
     4   -.2029E+02, -.1990E+02, -.1584E+02, -.1182E+02, -.7868E+01, 
     5   -.4269E+01, -.1806E+01, .1645E+00, -.2022E+02, -.1985E+02, 
     6   -.1538E+02, -.1136E+02, -.7417E+01, -.3820E+01, -.1373E+01, 
     7   .5657E+00, -.2018E+02, -.1981E+02, -.1492E+02, -.1090E+02, 
     8   -.6965E+01, -.3369E+01, -.9319E+00, .9577E+00, -.2013E+02, 
     9   -.1937E+02, -.1446E+02, -.1044E+02, -.6512E+01, -.2917E+01, 
     .   -.4928E+00, .1376E+01, -.2009E+02, -.1891E+02, -.1400E+02, 
     1   -.9984E+01, -.6063E+01, -.2466E+01, -.6887E-01, .1768E+01, 
     2   -.2006E+02, -.1845E+02, -.1354E+02, -.9530E+01, -.5618E+01, 
     3   -.2024E+01, .3615E+00, .2196E+01, -.2003E+02, -.1800E+02, 
     4   -.1308E+02, -.9075E+01, -.5174E+01, -.1593E+01, .7820E+00, 
     5   .2600E+01, -.1827E+02, -.1464E+02, -.1097E+02, -.7525E+01, 
     6   -.3733E+01, -.1077E+01, .1204E+01, .3014E+01, -.1525E+02, 
     7   -.1210E+02, -.9275E+01, -.5876E+01, -.2768E+01, -.6286E+00, 
     8   .1622E+01, .3394E+01, -.1298E+02, -.1060E+02, -.7764E+01, 
     9   -.4462E+01, -.2154E+01, -.2001E+00, .2034E+01, .3756E+01, 
     .   -.1157E+02, -.8941E+01, -.5984E+01, -.3509E+01, -.1651E+01, 
     1   .2279E+00, .2422E+01, .4066E+01, -.9986E+01, -.7062E+01, 
     2   -.4794E+01, -.2818E+01, -.1196E+01, .6394E+00, .2791E+01, 
     3   .4283E+01, -.8064E+01, -.5512E+01, -.3933E+01, -.2274E+01, 
     4   -.7559E+00, .1036E+01, .3085E+01, .4444E+01, -.6440E+01, 
     5   -.4863E+01, -.3219E+01, -.1791E+01, -.3279E+00, .1427E+01, 
     6   .3304E+01, .4527E+01, -.5902E+01, -.4207E+01, -.2756E+01, 
     7   -.1350E+01, .7686E-01, .1776E+01, .3475E+01, .4550E+01, 
     8   -.5439E+01, -.3739E+01, -.2330E+01, -.9233E+00, .4612E+00, 
     9   .2066E+01, .3564E+01, .4502E+01, -.5006E+01, -.3316E+01, 
     .   -.1906E+01, -.5066E+00, .8352E+00, .2272E+01, .3587E+01, 
     1   .4419E+01, .2338E-01, .1968E-02, .9503E-02, .3412E-02, 
     2   .6280E-03, -.1109E-02, -.1089E-02, -.1026E-02, .1972E-01, 
     3   .2093E-02, .9503E-02, .3391E-02, .6489E-03, -.1172E-02, 
     4   -.1164E-02, -.1158E-02, .1603E-01, .3328E-02, .9524E-02, 
     5   .3391E-02, .6489E-03, -.1277E-02, -.1229E-02, -.1296E-02, 
     6   .1229E-01, .7138E-02, .9524E-02, .3370E-02, .6070E-03, 
     7   -.1319E-02, -.1264E-02, -.1610E-02, .8478E-02, .1095E-01, 
     8   .9566E-02, .3412E-02, .5652E-03, -.1382E-02, -.1266E-02, 
     9   -.1566E-02, .4563E-02, .1480E-01, .9566E-02, .3412E-02, 
     .   .5443E-03, -.1423E-02, -.1199E-02, -.1679E-02, .2261E-02, 
     1   .1865E-01, .9608E-02, .3454E-02, .4815E-03, -.1423E-02, 
     2   -.1296E-02, -.1555E-02, .2198E-02, .2250E-01, .9671E-02, 
     3   .3412E-02, .4187E-03, -.1426E-02, -.1472E-02, -.1800E-02, 
     4   .2072E-02, .2600E-01, .9734E-02, .3433E-02, .3977E-03, 
     5   -.1428E-02, -.1541E-02, -.1591E-02, .1987E-01, .8645E-02, 
     6   .6280E-02, .1298E-02, -.1151E-02, -.1509E-02, -.1662E-02, 
     7   -.1570E-02, .4668E-02, .8373E-02, .3956E-02, -.4187E-04, 
     8   -.1968E-02, -.1624E-02, -.1700E-02, -.1947E-02, .9231E-02, 
     9   .5694E-02, .1444E-02, -.2512E-03, -.1827E-02, -.1662E-02, 
     .   -.1576E-02, -.1633E-02, .8666E-02, .3077E-02, -.1737E-02, 
     1   -.1277E-02, -.1507E-02, -.1757E-02, -.1612E-02, -.1612E-02, 
     2   .8164E-03, -.4375E-02, -.1884E-02, -.1277E-02, -.1564E-02, 
     3   -.1853E-02, -.1591E-02, -.1486E-02, -.1486E-02, -.2596E-02, 
     4   -.1633E-02, -.1539E-02, -.1662E-02, -.1846E-02, -.1423E-02, 
     5   -.1277E-02, -.1423E-02, -.2617E-02, -.1005E-02, -.1379E-02, 
     6   -.1687E-02, -.1905E-02, -.1528E-02, -.1298E-02, -.1675E-03, 
     7   -.1947E-02, -.5024E-03, -.1325E-02, -.1696E-02, -.1698E-02, 
     8   -.1486E-02, -.1277E-02, .1047E-03, -.1109E-02, -.5861E-03, 
     9   -.1363E-02, -.1620E-02, -.1666E-02, -.1507E-02, -.9210E-03, 
     .   .1047E-03, -.1047E-02, -.8394E-03, -.1342E-02, -.1591E-02, 
     1   -.1323E-02, -.1340E-02, -.9420E-03, -.1085E-03, .2283E-05, 
     2   -.4719E-04, -.3807E-06, -.1522E-05, -.3425E-05, -.7612E-06, 
     3   .1751E-05, -.1766E-03, .1523E-05, -.4719E-04, -.7609E-06, 
     4   -.3807E-06, -.3045E-05, .1599E-05, .8723E-05, -.2443E-03, 
     5   .1941E-04, -.4757E-04, -.1522E-05, -.3806E-06, -.1903E-05, 
     6   -.2778E-05, .1294E-04, -.1838E-03, .8563E-04, -.4757E-04, 
     7   -.1903E-05, .1142E-05, -.2664E-05, -.6090E-06, .1321E-04, 
     8   -.1161E-03, .1526E-03, -.4757E-04, -.2664E-05, -.3805E-06, 
     9   -.3806E-05, -.2093E-05, .2253E-04, -.4795E-04, .9248E-04, 
     .   -.4757E-04, -.1903E-05, .0000E+00, -.3045E-05, -.7992E-06, 
     1   .1393E-04, -.9134E-05, .2246E-04, -.4834E-04, -.2664E-05, 
     2   .3804E-06, -.5328E-05, -.1510E-05, .1465E-04, -.1028E-04, 
     3   -.4757E-04, -.4948E-04, -.1142E-05, .7614E-06, -.4910E-05, 
     4   -.5709E-06, .1477E-04, -.1256E-04, -.1066E-03, -.4910E-04, 
     5   -.1523E-05, -.3805E-06, -.3121E-05, -.2512E-05, .1142E-04, 
     6   -.7878E-04, -.2664E-05, -.8373E-05, -.7612E-06, .1104E-04, 
     7   -.3311E-05, -.1979E-05, .5709E-05, -.2626E-04, -.4872E-04, 
     8   -.3808E-06, -.2283E-05, .2284E-05, -.3349E-05, -.4034E-05, 
     9   .7231E-05, -.4910E-04, .1599E-04, .1256E-04, -.7612E-05, 
     .   .1180E-05, -.1815E-05, -.7193E-05, .3045E-05, .1576E-09, 
     1   .6470E-05, -.1408E-04, -.1903E-05, .1522E-05, -.4746E-05, 
     2   -.4948E-05, .3806E-06, .9020E-04, .5214E-04, .6090E-05, 
     3   -.1104E-04, .1180E-05, -.2778E-05, -.6090E-05, -.2664E-05, 
     4   -.6737E-04, -.1218E-04, -.3806E-05, -.5214E-05, -.1066E-05, 
     5   -.1294E-05, -.3045E-05, -.2664E-05, -.4643E-04, .1713E-04, 
     6   -.1218E-04, -.6204E-05, -.2360E-05, -.1979E-05, -.1903E-05, 
     7   -.3806E-05, -.3045E-04, -.1256E-04, -.9134E-05, -.6508E-05, 
     8   -.1027E-05, -.7993E-06, -.1142E-05, -.7992E-05, -.3616E-04, 
     9   -.1028E-04, -.1066E-04, -.6051E-05, .1066E-05, -.1751E-05, 
     .   -.2284E-05, -.2284E-05, -.3920E-04, -.9895E-05, -.1321E-04, 
     1   -.3844E-05, -.2055E-05, -.2512E-05, -.3806E-05,  - .3425E-05/ 
      END 
 
c---------- 4/1/97 (7) -- NEXT 1142 LINES -- Replaces old
c                         aerosol1,aerosol2 block data.
      BLOCK DATA AEROSOL1
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      INTEGER I, J
c#      real a_ssax(mbx,nrh,naer),a_extx(mbx,nrh,naer)
c#      real a_asyx(mbx,nrh,naer)
      COMMON /AER_OPTX/A_SSAX, A_EXTX, A_ASYX
      REAL A_SSAX(mbx,nrh,naer), A_EXTX(mbx,nrh,naer)
      REAL A_ASYX(mbx,nrh,naer)

      DATA ((A_SSAX(I,J,1),I=1,mbx),J=1,nrh)/ .1000E+01, .9984E+00, 
     1   .9525E+00, .9053E+00, .7378E+00, .8873E+00, .8528E+00, 
     2   .8678E+00, .6329E+00, .7734E+00, .7571E+00, .7446E+00, 
     3   .5500E+00, .3973E+00, .4265E+00, .4511E+00, .4341E+00, 
     4   .3346E+00, .1000E+01, .9974E+00, .9586E+00, .9109E+00, 
     5   .7298E+00, .8807E+00, .8421E+00, .8447E+00, .6212E+00, 
     6   .7637E+00, .7352E+00, .7322E+00, .5276E+00, .3942E+00, 
     7   .4226E+00, .4474E+00, .4344E+00, .3404E+00, .1000E+01, 
     8   .9980E+00, .9691E+00, .9182E+00, .7075E+00, .8584E+00, 
     9   .8072E+00, .8201E+00, .5870E+00, .7255E+00, .6977E+00, 
     .   .6968E+00, .4866E+00, .3946E+00, .4212E+00, .4429E+00, 
     1   .4396E+00, .3688E+00, .1000E+01, .9988E+00, .9820E+00, 
     2   .9212E+00, .6840E+00, .8189E+00, .7384E+00, .7583E+00, 
     3   .5412E+00, .6484E+00, .6295E+00, .6340E+00, .4620E+00, 
     4   .4177E+00, .4341E+00, .4484E+00, .4522E+00, .4161E+00, 
     5   .1000E+01, .9989E+00, .9836E+00, .9178E+00, .6825E+00, 
     6   .8084E+00, .7180E+00, .7351E+00, .5334E+00, .6226E+00, 
     7   .6058E+00, .6108E+00, .4623E+00, .4255E+00, .4399E+00, 
     8   .4518E+00, .4559E+00, .4284E+00, .1000E+01, .9990E+00, 
     9   .9832E+00, .9107E+00, .6815E+00, .7994E+00, .7018E+00, 
     .   .7143E+00, .5313E+00, .6011E+00, .5836E+00, .5877E+00, 
     1   .4635E+00, .4341E+00, .4456E+00, .4551E+00, .4589E+00, 
     2   .4382E+00, .1000E+01, .9987E+00, .9813E+00, .8925E+00, 
     3   .6748E+00, .7865E+00, .6908E+00, .6951E+00, .5373E+00, 
     4   .5836E+00, .5624E+00, .5605E+00, .4652E+00, .4443E+00, 
     5   .4537E+00, .4598E+00, .4620E+00, .4474E+00, .1000E+01, 
     6   .9988E+00, .9800E+00, .8969E+00, .6654E+00, .7781E+00, 
     7   .6947E+00, .6954E+00, .5480E+00, .5842E+00, .5572E+00, 
     8   .5477E+00, .4642E+00, .4479E+00, .4572E+00, .4614E+00, 
     9   .4620E+00, .4495E+00/ 
      DATA ((A_EXTX(I,J,1),I=1,mbx),J=1,nrh)/ .2085E-03, .2085E-03, 
     1   .1753E-03, .1667E-03, .1655E-03, .1667E-03, .1721E-03, 
     2   .1735E-03, .1698E-03, .1700E-03, .1691E-03, .1647E-03, 
     3   .1267E-03, .1256E-03, .1477E-03, .1473E-03, .1320E-03, 
     4   .1206E-03, .2442E-03, .2391E-03, .1959E-03, .1850E-03, 
     5   .1841E-03, .1836E-03, .1895E-03, .1909E-03, .1867E-03, 
     6   .1895E-03, .1879E-03, .1794E-03, .1379E-03, .1395E-03, 
     7   .1642E-03, .1644E-03, .1482E-03, .1336E-03, .3488E-03, 
     8   .3479E-03, .3010E-03, .2796E-03, .2720E-03, .2663E-03, 
     9   .2693E-03, .2725E-03, .2678E-03, .2743E-03, .2717E-03, 
     .   .2589E-03, .2028E-03, .2152E-03, .2470E-03, .2496E-03, 
     1   .2322E-03, .2076E-03, .7848E-03, .7872E-03, .7928E-03, 
     2   .7466E-03, .7085E-03, .6744E-03, .6381E-03, .6362E-03, 
     3   .6401E-03, .6470E-03, .6477E-03, .6350E-03, .5307E-03, 
     4   .5726E-03, .6321E-03, .6438E-03, .6297E-03, .5842E-03, 
     5   .1112E-02, .1113E-02, .1148E-02, .1112E-02, .1057E-02, 
     6   .1004E-02, .9203E-03, .9076E-03, .9195E-03, .9147E-03, 
     7   .9172E-03, .9072E-03, .7833E-03, .8441E-03, .9175E-03, 
     8   .9317E-03, .9216E-03, .8724E-03, .1636E-02, .1619E-02, 
     9   .1667E-02, .1673E-02, .1619E-02, .1548E-02, .1385E-02, 
     .   .1345E-02, .1367E-02, .1335E-02, .1334E-02, .1324E-02, 
     1   .1184E-02, .1269E-02, .1366E-02, .1379E-02, .1373E-02, 
     2   .1323E-02, .2803E-02, .2748E-02, .2765E-02, .2829E-02, 
     3   .2862E-02, .2813E-02, .2508E-02, .2396E-02, .2421E-02, 
     4   .2312E-02, .2280E-02, .2252E-02, .2093E-02, .2240E-02, 
     5   .2390E-02, .2388E-02, .2373E-02, .2328E-02, .4213E-02, 
     6   .4113E-02, .4088E-02, .4098E-02, .4248E-02, .4287E-02, 
     7   .3951E-02, .3743E-02, .3733E-02, .3520E-02, .3416E-02, 
     8   .3331E-02, .3154E-02, .3390E-02, .3609E-02, .3580E-02, 
     9   .3527E-02, .3473E-02/ 
      DATA ((A_ASYX(I,J,1),I=1,mbx),J=1,nrh)/ .7972E+00, .8182E+00, 
     1   .8172E+00, .8200E+00, .8119E+00, .7766E+00, .8040E+00, 
     2   .8212E+00, .8646E+00, .8447E+00, .8440E+00, .8411E+00, 
     3   .8880E+00, .8602E+00, .7911E+00, .7291E+00, .6673E+00, 
     4   .5545E+00, .8017E+00, .8218E+00, .8187E+00, .8216E+00, 
     5   .8160E+00, .7809E+00, .8095E+00, .8488E+00, .8715E+00, 
     6   .8498E+00, .8488E+00, .8597E+00, .8958E+00, .8652E+00, 
     7   .7976E+00, .7375E+00, .6763E+00, .5685E+00, .7986E+00, 
     8   .8234E+00, .8312E+00, .8353E+00, .8296E+00, .7968E+00, 
     9   .8248E+00, .8507E+00, .8891E+00, .8648E+00, .8726E+00, 
     .   .8853E+00, .9177E+00, .8834E+00, .8248E+00, .7727E+00, 
     1   .7198E+00, .6379E+00, .7617E+00, .8120E+00, .8494E+00, 
     2   .8614E+00, .8610E+00, .8308E+00, .8540E+00, .8626E+00, 
     3   .9124E+00, .8874E+00, .9025E+00, .9183E+00, .9476E+00, 
     4   .9123E+00, .8683E+00, .8336E+00, .7948E+00, .7445E+00, 
     5   .7412E+00, .7992E+00, .8491E+00, .8673E+00, .8711E+00, 
     6   .8437E+00, .8652E+00, .8700E+00, .9176E+00, .8950E+00, 
     7   .9099E+00, .9256E+00, .9550E+00, .9187E+00, .8787E+00, 
     8   .8512E+00, .8183E+00, .7759E+00, .7144E+00, .7752E+00, 
     9   .8417E+00, .8684E+00, .8779E+00, .8554E+00, .8775E+00, 
     .   .8804E+00, .9226E+00, .9026E+00, .9169E+00, .9319E+00, 
     1   .9607E+00, .9236E+00, .8850E+00, .8645E+00, .8394E+00, 
     2   .8044E+00, .6858E+00, .7430E+00, .8251E+00, .8605E+00, 
     3   .8799E+00, .8649E+00, .8931E+00, .8955E+00, .9294E+00, 
     4   .9133E+00, .9253E+00, .9394E+00, .9660E+00, .9273E+00, 
     5   .8877E+00, .8751E+00, .8610E+00, .8366E+00, .6686E+00, 
     6   .7251E+00, .8155E+00, .8500E+00, .8752E+00, .8642E+00, 
     7   .9001E+00, .9040E+00, .9324E+00, .9183E+00, .9292E+00, 
     8   .9420E+00, .9677E+00, .9280E+00, .8855E+00, .8724E+00, 
     9   .8665E+00, .8517E+00/ 
      DATA ((A_SSAX(I,J,2),I=1,mbx),J=1,nrh)/ .9607E+00, .9253E+00, 
     1   .7650E+00, .3869E+00, .7830E+00, .8196E+00, .5468E+00, 
     2   .3954E+00, .2303E+00, .6683E-01, .8012E-01, .1274E+00, 
     3   .1627E+00, .9903E-01, .5161E-01, .4431E-01, .2697E-01, 
     4   .1631E-01, .9606E+00, .9252E+00, .7650E+00, .3872E+00, 
     5   .7821E+00, .8195E+00, .5486E+00, .3983E+00, .2330E+00, 
     6   .6891E-01, .8092E-01, .1285E+00, .1625E+00, .1015E+00, 
     7   .5113E-01, .4522E-01, .2781E-01, .1691E-01, .9632E+00, 
     8   .9301E+00, .7820E+00, .4110E+00, .7464E+00, .8202E+00, 
     9   .5511E+00, .4098E+00, .2105E+00, .7610E-01, .8126E-01, 
     .   .1259E+00, .1316E+00, .6796E-01, .4130E-01, .4058E-01, 
     1   .2661E-01, .1672E-01, .9730E+00, .9487E+00, .8461E+00, 
     2   .5175E+00, .7033E+00, .8338E+00, .5724E+00, .4600E+00, 
     3   .1834E+00, .1095E+00, .8760E-01, .1199E+00, .7362E-01, 
     4   .2678E-01, .2572E-01, .3075E-01, .2423E-01, .1656E-01, 
     5   .9820E+00, .9667E+00, .9056E+00, .6542E+00, .7047E+00, 
     6   .8543E+00, .6027E+00, .5125E+00, .1824E+00, .1479E+00, 
     7   .1006E+00, .1160E+00, .4699E-01, .1763E-01, .2012E-01, 
     8   .2466E-01, .2149E-01, .1529E-01, .9859E+00, .9745E+00, 
     9   .9303E+00, .7255E+00, .7137E+00, .8662E+00, .6230E+00, 
     .   .5426E+00, .1894E+00, .1718E+00, .1117E+00, .1168E+00, 
     1   .3984E-01, .1625E-01, .1928E-01, .2322E-01, .2079E-01, 
     2   .1493E-01, .9891E+00, .9808E+00, .9500E+00, .7911E+00, 
     3   .7245E+00, .8778E+00, .6552E+00, .5913E+00, .2128E+00, 
     4   .2271E+00, .1459E+00, .1391E+00, .4313E-01, .2095E-01, 
     5   .2574E-01, .3247E-01, .3267E-01, .2510E-01, .9914E+00, 
     6   .9853E+00, .9630E+00, .8391E+00, .7353E+00, .8871E+00, 
     7   .6780E+00, .6222E+00, .2295E+00, .2635E+00, .1717E+00, 
     8   .1553E+00, .4561E-01, .2404E-01, .2966E-01, .3769E-01, 
     9   .3967E-01, .3733E-01/ 
      DATA ((A_EXTX(I,J,2),I=1,mbx),J=1,nrh)/ .1067E-04, .5658E-05, 
     1   .1248E-05, .1317E-05, .2144E-06, .1635E-06, .1051E-06, 
     2   .1039E-06, .1074E-06, .1852E-06, .3665E-06, .2548E-06, 
     3   .8879E-07, .9337E-07, .1557E-06, .1269E-06, .1362E-06, 
     4   .1536E-06, .1067E-04, .5659E-05, .1250E-05, .1318E-05, 
     5   .2156E-06, .1645E-06, .1060E-06, .1047E-06, .1083E-06, 
     6   .1859E-06, .3671E-06, .2554E-06, .8921E-07, .9426E-07, 
     7   .1563E-06, .1274E-06, .1366E-06, .1539E-06, .1145E-04, 
     8   .6089E-05, .1366E-05, .1390E-05, .2705E-06, .1893E-06, 
     9   .1202E-06, .1160E-06, .1392E-06, .1982E-06, .3933E-06, 
     .   .2795E-06, .1146E-06, .1395E-06, .1987E-06, .1543E-06, 
     1   .1569E-06, .1719E-06, .1554E-04, .8394E-05, .2017E-05, 
     2   .1792E-05, .5780E-06, .3314E-06, .1978E-06, .1767E-06, 
     3   .3055E-06, .2637E-06, .5149E-06, .3873E-06, .2577E-06, 
     4   .4005E-06, .4371E-06, .3044E-06, .2654E-06, .2658E-06, 
     5   .2344E-04, .1308E-04, .3456E-05, .2666E-05, .1253E-05, 
     6   .6574E-06, .3619E-06, .3005E-06, .6400E-06, .3929E-06, 
     7   .7060E-06, .5523E-06, .5532E-06, .9393E-06, .9290E-06, 
     8   .6129E-06, .4778E-06, .4448E-06, .3004E-04, .1716E-04, 
     9   .4801E-05, .3491E-05, .1886E-05, .9781E-06, .5168E-06, 
     .   .4150E-06, .9341E-06, .5071E-06, .8499E-06, .6765E-06, 
     1   .8087E-06, .1406E-05, .1356E-05, .8808E-06, .6589E-06, 
     2   .5955E-06, .3935E-04, .2315E-04, .6913E-05, .4819E-05, 
     3   .2908E-05, .1535E-05, .8007E-06, .6334E-06, .1408E-05, 
     4   .7153E-06, .1081E-05, .8758E-06, .1200E-05, .2109E-05, 
     5   .2009E-05, .1301E-05, .9483E-06, .8356E-06, .5037E-04, 
     6   .3051E-04, .9659E-05, .6565E-05, .4238E-05, .2277E-05, 
     7   .1165E-05, .9083E-06, .1992E-05, .9700E-06, .1351E-05, 
     8   .1111E-05, .1676E-05, .2962E-05, .2801E-05, .1811E-05, 
     9   .1302E-05, .9664E-06/ 
      DATA ((A_ASYX(I,J,2),I=1,mbx),J=1,nrh)/ .6406E+00, .6057E+00, 
     1   .5447E+00, .4976E+00, .4323E+00, .4216E+00, .4084E+00, 
     2   .4038E+00, .3530E+00, .5334E+00, .4666E+00, .3619E+00, 
     3   .4654E+00, .5418E+00, .5190E+00, .4775E+00, .4633E+00, 
     4   .3869E+00, .6406E+00, .6057E+00, .5449E+00, .4982E+00, 
     5   .4338E+00, .4240E+00, .4135E+00, .4106E+00, .3639E+00, 
     6   .5480E+00, .4744E+00, .3681E+00, .4720E+00, .5471E+00, 
     7   .5244E+00, .4836E+00, .4694E+00, .3936E+00, .6514E+00, 
     8   .6161E+00, .5532E+00, .5076E+00, .4378E+00, .4297E+00, 
     9   .4202E+00, .4202E+00, .3811E+00, .5519E+00, .4816E+00, 
     .   .3768E+00, .4809E+00, .5492E+00, .5237E+00, .4896E+00, 
     1   .4800E+00, .4086E+00, .6892E+00, .6537E+00, .5854E+00, 
     2   .5436E+00, .4553E+00, .4509E+00, .4419E+00, .4494E+00, 
     3   .4211E+00, .5500E+00, .4860E+00, .4133E+00, .5070E+00, 
     4   .5397E+00, .5109E+00, .5082E+00, .5147E+00, .4602E+00, 
     5   .7238E+00, .6909E+00, .6212E+00, .5832E+00, .4777E+00, 
     6   .4729E+00, .4524E+00, .4581E+00, .4252E+00, .5112E+00, 
     7   .4636E+00, .4340E+00, .5076E+00, .4990E+00, .4697E+00, 
     8   .4977E+00, .5255E+00, .4928E+00, .7390E+00, .7084E+00, 
     9   .6399E+00, .6042E+00, .4924E+00, .4874E+00, .4601E+00, 
     .   .4628E+00, .4269E+00, .4922E+00, .4535E+00, .4423E+00, 
     1   .5024E+00, .4712E+00, .4434E+00, .4853E+00, .5258E+00, 
     2   .5079E+00, .7522E+00, .7245E+00, .6593E+00, .6270E+00, 
     3   .5142E+00, .5121E+00, .4943E+00, .5049E+00, .4707E+00, 
     4   .5452E+00, .5173E+00, .5260E+00, .5827E+00, .5439E+00, 
     5   .5142E+00, .5589E+00, .5958E+00, .5818E+00, .7620E+00, 
     6   .7371E+00, .6754E+00, .6456E+00, .5312E+00, .5297E+00, 
     7   .5112E+00, .5218E+00, .4872E+00, .5553E+00, .5396E+00, 
     8   .5633E+00, .6177E+00, .5629E+00, .5322E+00, .5821E+00, 
     9   .6221E+00, .6250E+00/ 
      DATA ((A_SSAX(I,J,3),I=1,mbx),J=1,nrh)/ .9371E+00, .8999E+00, 
     1   .7175E+00, .3628E+00, .6462E+00, .6564E+00, .4011E+00, 
     2   .2856E+00, .1754E+00, .3630E-01, .6500E-01, .8672E-01, 
     3   .8039E-01, .3570E-01, .1633E-01, .1202E-01, .4884E-02, 
     4   .2383E-02, .9365E+00, .8992E+00, .7160E+00, .3622E+00, 
     5   .6417E+00, .6511E+00, .3969E+00, .2828E+00, .1736E+00, 
     6   .3634E-01, .6489E-01, .8637E-01, .7938E-01, .3683E-01, 
     7   .1618E-01, .1207E-01, .4972E-02, .2454E-02, .9386E+00, 
     8   .9035E+00, .7316E+00, .3838E+00, .6180E+00, .6530E+00, 
     9   .3951E+00, .2856E+00, .1549E+00, .3982E-01, .6398E-01, 
     .   .8328E-01, .6367E-01, .2501E-01, .1367E-01, .1095E-01, 
     1   .4813E-02, .2445E-02, .9522E+00, .9265E+00, .8037E+00, 
     2   .4882E+00, .6214E+00, .7102E+00, .4367E+00, .3326E+00, 
     3   .1370E+00, .5985E-01, .6507E-01, .7825E-01, .3715E-01, 
     4   .1113E-01, .9892E-02, .9000E-02, .4744E-02, .2534E-02, 
     5   .9669E+00, .9502E+00, .8773E+00, .6289E+00, .6589E+00, 
     6   .7812E+00, .5069E+00, .4094E+00, .1465E+00, .9469E-01, 
     7   .7425E-01, .7846E-01, .2566E-01, .8722E-02, .9422E-02, 
     8   .8668E-02, .5022E-02, .2619E-02, .9733E+00, .9620E+00, 
     9   .9086E+00, .7041E+00, .6797E+00, .8128E+00, .5458E+00, 
     .   .4545E+00, .1573E+00, .1192E+00, .8364E-01, .8086E-01, 
     1   .2264E-01, .8812E-02, .9970E-02, .9079E-02, .5429E-02, 
     2   .2765E-02, .9790E+00, .9710E+00, .9339E+00, .7734E+00, 
     3   .6992E+00, .8399E+00, .5866E+00, .5038E+00, .1729E+00, 
     4   .1520E+00, .9926E-01, .8757E-01, .2187E-01, .1003E-01, 
     5   .1172E-01, .1103E-01, .7504E-02, .4245E-02, .9832E+00, 
     6   .9776E+00, .9511E+00, .8253E+00, .7164E+00, .8604E+00, 
     7   .6219E+00, .5461E+00, .1893E+00, .1831E+00, .1159E+00, 
     8   .9521E-01, .2208E-01, .1143E-01, .1350E-01, .1280E-01, 
     9   .9112E-02, .6206E-02/ 
      DATA ((A_EXTX(I,J,3),I=1,mbx),J=1,nrh)/ .6974E-05, .3689E-05, 
     1   .8308E-06, .8639E-06, .1517E-06, .1172E-06, .7603E-07, 
     2   .7487E-07, .7732E-07, .1226E-06, .2351E-06, .1602E-06, 
     3   .5585E-07, .5913E-07, .9781E-07, .7830E-07, .8414E-07, 
     4   .9499E-07, .6982E-05, .3693E-05, .8327E-06, .8656E-06, 
     5   .1529E-06, .1183E-06, .7697E-07, .7577E-07, .7826E-07, 
     6   .1234E-06, .2358E-06, .1609E-06, .5659E-07, .6006E-07, 
     7   .9840E-07, .7871E-07, .8445E-07, .9516E-07, .7505E-05, 
     8   .3978E-05, .9110E-06, .9159E-06, .1901E-06, .1358E-06, 
     9   .8773E-07, .8470E-07, .9992E-07, .1329E-06, .2542E-06, 
     .   .1781E-06, .7504E-07, .9126E-07, .1270E-06, .9669E-07, 
     1   .9796E-07, .1067E-06, .1017E-04, .5464E-05, .1328E-05, 
     2   .1172E-05, .3857E-06, .2239E-06, .1355E-06, .1219E-06, 
     3   .2064E-06, .1743E-06, .3321E-06, .2476E-06, .1704E-06, 
     4   .2630E-06, .2819E-06, .1932E-06, .1672E-06, .1659E-06, 
     5   .1528E-04, .8466E-05, .2240E-05, .1723E-05, .8128E-06, 
     6   .4251E-06, .2346E-06, .1956E-06, .4180E-06, .2532E-06, 
     7   .4515E-06, .3513E-06, .3619E-06, .6124E-06, .5990E-06, 
     8   .3904E-06, .3021E-06, .2789E-06, .1963E-04, .1111E-04, 
     9   .3099E-05, .2247E-05, .1219E-05, .6254E-06, .3299E-06, 
     .   .2652E-06, .6084E-06, .3245E-06, .5427E-06, .4307E-06, 
     1   .5313E-06, .9227E-06, .8816E-06, .5662E-06, .4201E-06, 
     2   .3763E-06, .2594E-04, .1505E-04, .4440E-05, .3077E-05, 
     3   .1861E-05, .9545E-06, .4851E-06, .3783E-06, .8998E-06, 
     4   .4362E-06, .6725E-06, .5438E-06, .7839E-06, .1386E-05, 
     5   .1306E-05, .8310E-06, .5971E-06, .5214E-06, .3343E-04, 
     6   .1989E-04, .6191E-05, .4176E-05, .2703E-05, .1404E-05, 
     7   .6921E-06, .5279E-06, .3548E-06, .5771E-06, .8264E-06, 
     8   .6774E-06, .1089E-05, .1946E-05, .1820E-05, .1152E-05, 
     9   .8146E-06, .5980E-06/ 
      DATA ((A_ASYX(I,J,3),I=1,mbx),J=1,nrh)/ .6381E+00, .6035E+00, 
     1   .5386E+00, .4849E+00, .3957E+00, .3761E+00, .3199E+00, 
     2   .3006E+00, .2684E+00, .2713E+00, .2763E+00, .2241E+00, 
     3   .2351E+00, .2725E+00, .2616E+00, .2607E+00, .3195E+00, 
     4   .3162E+00, .6381E+00, .6035E+00, .5385E+00, .4849E+00, 
     5   .3958E+00, .3764E+00, .3207E+00, .3018E+00, .2700E+00, 
     6   .2780E+00, .2836E+00, .2254E+00, .2366E+00, .2730E+00, 
     7   .2629E+00, .2658E+00, .3269E+00, .3239E+00, .6490E+00, 
     8   .6137E+00, .5468E+00, .4946E+00, .4020E+00, .3834E+00, 
     9   .3274E+00, .3089E+00, .2773E+00, .2849E+00, .2936E+00, 
     .   .2285E+00, .2374E+00, .2654E+00, .2526E+00, .2639E+00, 
     1   .3284E+00, .3327E+00, .6866E+00, .6512E+00, .5797E+00, 
     2   .5327E+00, .4280E+00, .4121E+00, .3550E+00, .3374E+00, 
     3   .3039E+00, .3043E+00, .2931E+00, .2444E+00, .2438E+00, 
     4   .2427E+00, .2265E+00, .2565E+00, .3291E+00, .3575E+00, 
     5   .7213E+00, .6884E+00, .6168E+00, .5756E+00, .4606E+00, 
     6   .4473E+00, .3883E+00, .3700E+00, .3315E+00, .3171E+00, 
     7   .2900E+00, .2650E+00, .2506E+00, .2220E+00, .2017E+00, 
     8   .2297E+00, .2957E+00, .3450E+00, .7362E+00, .7055E+00, 
     9   .6357E+00, .5977E+00, .4792E+00, .4674E+00, .4080E+00, 
     .   .3893E+00, .3483E+00, .3287E+00, .2983E+00, .2787E+00, 
     1   .2562E+00, .2159E+00, .1939E+00, .2169E+00, .2754E+00, 
     2   .3308E+00, .7488E+00, .7206E+00, .6539E+00, .6192E+00, 
     3   .4991E+00, .4892E+00, .4329E+00, .4163E+00, .3733E+00, 
     4   .3608E+00, .3302E+00, .3274E+00, .2936E+00, .2500E+00, 
     5   .2283E+00, .2661E+00, .3470E+00, .4161E+00, .7584E+00, 
     6   .7326E+00, .6694E+00, .6376E+00, .5167E+00, .5082E+00, 
     7   .4535E+00, .4373E+00, .3927E+00, .3798E+00, .3500E+00, 
     8   .3381E+00, .3171E+00, .2640E+00, .2404E+00, .2809E+00, 
     9   .3669E+00, .4590E+00/ 
      DATA ((A_SSAX(I,J,4),I=1,mbx),J=1,nrh)/ .9140E+00, .9726E+00, 
     1   .9759E+00, .9737E+00, .8492E+00, .8986E+00, .8344E+00, 
     2   .6125E+00, .2537E+00, .9996E-01, .3744E-01, .1756E+00, 
     3   .6959E-01, .3767E-01, .1425E-01, .1772E-01, .7060E-02, 
     4   .2826E-02, .9140E+00, .9726E+00, .9759E+00, .9737E+00, 
     5   .8492E+00, .8986E+00, .8344E+00, .6125E+00, .2537E+00, 
     6   .9996E-01, .3744E-01, .1756E+00, .6959E-01, .3767E-01, 
     7   .1425E-01, .1772E-01, .7060E-02, .2826E-02, .9140E+00, 
     8   .9726E+00, .9759E+00, .9737E+00, .8492E+00, .8986E+00, 
     9   .8344E+00, .6125E+00, .2537E+00, .9996E-01, .3744E-01, 
     .   .1756E+00, .6959E-01, .3767E-01, .1425E-01, .1772E-01, 
     1   .7060E-02, .2826E-02, .9140E+00, .9726E+00, .9759E+00, 
     2   .9737E+00, .8492E+00, .8986E+00, .8344E+00, .6125E+00, 
     3   .2537E+00, .9996E-01, .3744E-01, .1756E+00, .6959E-01, 
     4   .3767E-01, .1425E-01, .1772E-01, .7060E-02, .2826E-02, 
     5   .9140E+00, .9726E+00, .9759E+00, .9737E+00, .8492E+00, 
     6   .8986E+00, .8344E+00, .6125E+00, .2537E+00, .9996E-01, 
     7   .3744E-01, .1756E+00, .6959E-01, .3767E-01, .1425E-01, 
     8   .1772E-01, .7060E-02, .2826E-02, .9140E+00, .9726E+00, 
     9   .9759E+00, .9737E+00, .8492E+00, .8986E+00, .8344E+00, 
     .   .6125E+00, .2537E+00, .9996E-01, .3744E-01, .1756E+00, 
     1   .6959E-01, .3767E-01, .1425E-01, .1772E-01, .7060E-02, 
     2   .2826E-02, .9140E+00, .9726E+00, .9759E+00, .9737E+00, 
     3   .8492E+00, .8986E+00, .8344E+00, .6125E+00, .2537E+00, 
     4   .9996E-01, .3744E-01, .1756E+00, .6959E-01, .3767E-01, 
     5   .1425E-01, .1772E-01, .7060E-02, .2826E-02, .9140E+00, 
     6   .9726E+00, .9759E+00, .9737E+00, .8492E+00, .8986E+00, 
     7   .8344E+00, .6125E+00, .2537E+00, .9996E-01, .3744E-01, 
     8   .1756E+00, .6959E-01, .3767E-01, .1425E-01, .1772E-01, 
     9   .7060E-02, .2826E-02/ 
      DATA ((A_EXTX(I,J,4),I=1,mbx),J=1,nrh)/ .1013E+01, .1046E+01, 
     1   .7036E+00, .4361E+00, .1101E+00, .7263E-01, .3980E-01, 
     2   .3442E-01, .3402E-01, .3102E-01, .7158E-01, .1016E+00, 
     3   .5528E-01, .2937E-01, .3969E-01, .3820E-01, .2108E-01, 
     4   .1806E-01, .1013E+01, .1046E+01, .7036E+00, .4361E+00, 
     5   .1101E+00, .7263E-01, .3980E-01, .3442E-01, .3402E-01, 
     6   .3102E-01, .7158E-01, .1016E+00, .5528E-01, .2937E-01, 
     7   .3969E-01, .3820E-01, .2108E-01, .1806E-01, .1013E+01, 
     8   .1046E+01, .7036E+00, .4361E+00, .1101E+00, .7263E-01, 
     9   .3980E-01, .3442E-01, .3402E-01, .3102E-01, .7158E-01, 
     .   .1016E+00, .5528E-01, .2937E-01, .3969E-01, .3820E-01, 
     1   .2108E-01, .1806E-01, .1013E+01, .1046E+01, .7036E+00, 
     2   .4361E+00, .1101E+00, .7263E-01, .3980E-01, .3442E-01, 
     3   .3402E-01, .3102E-01, .7158E-01, .1016E+00, .5528E-01, 
     4   .2937E-01, .3969E-01, .3820E-01, .2108E-01, .1806E-01, 
     5   .1013E+01, .1046E+01, .7036E+00, .4361E+00, .1101E+00, 
     6   .7263E-01, .3980E-01, .3442E-01, .3402E-01, .3102E-01, 
     7   .7158E-01, .1016E+00, .5528E-01, .2937E-01, .3969E-01, 
     8   .3820E-01, .2108E-01, .1806E-01, .1013E+01, .1046E+01, 
     9   .7036E+00, .4361E+00, .1101E+00, .7263E-01, .3980E-01, 
     .   .3442E-01, .3402E-01, .3102E-01, .7158E-01, .1016E+00, 
     1   .5528E-01, .2937E-01, .3969E-01, .3820E-01, .2108E-01, 
     2   .1806E-01, .1013E+01, .1046E+01, .7036E+00, .4361E+00, 
     3   .1101E+00, .7263E-01, .3980E-01, .3442E-01, .3402E-01, 
     4   .3102E-01, .7158E-01, .1016E+00, .5528E-01, .2937E-01, 
     5   .3969E-01, .3820E-01, .2108E-01, .1806E-01, .1013E+01, 
     6   .1046E+01, .7036E+00, .4361E+00, .1101E+00, .7263E-01, 
     7   .3980E-01, .3442E-01, .3402E-01, .3102E-01, .7158E-01, 
     8   .1016E+00, .5528E-01, .2937E-01, .3969E-01, .3820E-01, 
     9   .2108E-01, .1806E-01/ 
      DATA ((A_ASYX(I,J,4),I=1,mbx),J=1,nrh)/ .6727E+00, .6788E+00, 
     1   .6599E+00, .6079E+00, .4306E+00, .3754E+00, .2599E+00, 
     2   .2139E+00, .1488E+00, .1066E+00, .8476E-01, .1280E+00, 
     3   .6212E-01, .4009E-01, .2821E-01, .2439E-01, .1238E-01, 
     4   .7042E-02, .6727E+00, .6788E+00, .6599E+00, .6079E+00, 
     5   .4306E+00, .3754E+00, .2599E+00, .2139E+00, .1488E+00, 
     6   .1066E+00, .8476E-01, .1280E+00, .6212E-01, .4009E-01, 
     7   .2821E-01, .2439E-01, .1238E-01, .7042E-02, .6727E+00, 
     8   .6788E+00, .6599E+00, .6079E+00, .4306E+00, .3754E+00, 
     9   .2599E+00, .2139E+00, .1488E+00, .1066E+00, .8476E-01, 
     .   .1280E+00, .6212E-01, .4009E-01, .2821E-01, .2439E-01, 
     1   .1238E-01, .7042E-02, .6727E+00, .6788E+00, .6599E+00, 
     2   .6079E+00, .4306E+00, .3754E+00, .2599E+00, .2139E+00, 
     3   .1488E+00, .1066E+00, .8476E-01, .1280E+00, .6212E-01, 
     4   .4009E-01, .2821E-01, .2439E-01, .1238E-01, .7042E-02, 
     5   .6727E+00, .6788E+00, .6599E+00, .6079E+00, .4306E+00, 
     6   .3754E+00, .2599E+00, .2139E+00, .1488E+00, .1066E+00, 
     7   .8476E-01, .1280E+00, .6212E-01, .4009E-01, .2821E-01, 
     8   .2439E-01, .1238E-01, .7042E-02, .6727E+00, .6788E+00, 
     9   .6599E+00, .6079E+00, .4306E+00, .3754E+00, .2599E+00, 
     .   .2139E+00, .1488E+00, .1066E+00, .8476E-01, .1280E+00, 
     1   .6212E-01, .4009E-01, .2821E-01, .2439E-01, .1238E-01, 
     2   .7042E-02, .6727E+00, .6788E+00, .6599E+00, .6079E+00, 
     3   .4306E+00, .3754E+00, .2599E+00, .2139E+00, .1488E+00, 
     4   .1066E+00, .8476E-01, .1280E+00, .6212E-01, .4009E-01, 
     5   .2821E-01, .2439E-01, .1238E-01, .7042E-02, .6727E+00, 
     6   .6788E+00, .6599E+00, .6079E+00, .4306E+00, .3754E+00, 
     7   .2599E+00, .2139E+00, .1488E+00, .1066E+00, .8476E-01, 
     8   .1280E+00, .6212E-01, .4009E-01, .2821E-01, .2439E-01, 
     9   .1238E-01, .7042E-02/ 
      DATA ((A_SSAX(I,J,5),I=1,mbx),J=1,nrh)/ .8498E+00, .9415E+00, 
     1   .9649E+00, .9728E+00, .9141E+00, .9502E+00, .9317E+00, 
     2   .8228E+00, .5514E+00, .3158E+00, .1352E+00, .3908E+00, 
     3   .2884E+00, .1955E+00, .8936E-01, .1136E+00, .5145E-01, 
     4   .2186E-01, .8498E+00, .9415E+00, .9649E+00, .9728E+00, 
     5   .9141E+00, .9502E+00, .9317E+00, .8228E+00, .5514E+00, 
     6   .3158E+00, .1352E+00, .3908E+00, .2884E+00, .1955E+00, 
     7   .8936E-01, .1136E+00, .5145E-01, .2186E-01, .8498E+00, 
     8   .9415E+00, .9649E+00, .9728E+00, .9141E+00, .9502E+00, 
     9   .9317E+00, .8228E+00, .5514E+00, .3158E+00, .1352E+00, 
     .   .3908E+00, .2884E+00, .1955E+00, .8936E-01, .1136E+00, 
     1   .5145E-01, .2186E-01, .8498E+00, .9415E+00, .9649E+00, 
     2   .9728E+00, .9141E+00, .9502E+00, .9317E+00, .8228E+00, 
     3   .5514E+00, .3158E+00, .1352E+00, .3908E+00, .2884E+00, 
     4   .1955E+00, .8936E-01, .1136E+00, .5145E-01, .2186E-01, 
     5   .8498E+00, .9415E+00, .9649E+00, .9728E+00, .9141E+00, 
     6   .9502E+00, .9317E+00, .8228E+00, .5514E+00, .3158E+00, 
     7   .1352E+00, .3908E+00, .2884E+00, .1955E+00, .8936E-01, 
     8   .1136E+00, .5145E-01, .2186E-01, .8498E+00, .9415E+00, 
     9   .9649E+00, .9728E+00, .9141E+00, .9502E+00, .9317E+00, 
     .   .8228E+00, .5514E+00, .3158E+00, .1352E+00, .3908E+00, 
     1   .2884E+00, .1955E+00, .8936E-01, .1136E+00, .5145E-01, 
     2   .2186E-01, .8498E+00, .9415E+00, .9649E+00, .9728E+00, 
     3   .9141E+00, .9502E+00, .9317E+00, .8228E+00, .5514E+00, 
     4   .3158E+00, .1352E+00, .3908E+00, .2884E+00, .1955E+00, 
     5   .8936E-01, .1136E+00, .5145E-01, .2186E-01, .8498E+00, 
     6   .9415E+00, .9649E+00, .9728E+00, .9141E+00, .9502E+00, 
     7   .9317E+00, .8228E+00, .5514E+00, .3158E+00, .1352E+00, 
     8   .3908E+00, .2884E+00, .1955E+00, .8936E-01, .1136E+00, 
     9   .5145E-01, .2186E-01/ 
      DATA ((A_EXTX(I,J,5),I=1,mbx),J=1,nrh)/ .1011E+01, .1126E+01, 
     1   .1274E+01, .1194E+01, .5876E+00, .4705E+00, .3210E+00, 
     2   .2489E+00, .1574E+00, .1099E+00, .2069E+00, .5297E+00, 
     3   .1960E+00, .9338E-01, .1105E+00, .1188E+00, .5688E-01, 
     4   .4516E-01, .1011E+01, .1126E+01, .1274E+01, .1194E+01, 
     5   .5876E+00, .4705E+00, .3210E+00, .2489E+00, .1574E+00, 
     6   .1099E+00, .2069E+00, .5297E+00, .1960E+00, .9338E-01, 
     7   .1105E+00, .1188E+00, .5688E-01, .4516E-01, .1011E+01, 
     8   .1126E+01, .1274E+01, .1194E+01, .5876E+00, .4705E+00, 
     9   .3210E+00, .2489E+00, .1574E+00, .1099E+00, .2069E+00, 
     .   .5297E+00, .1960E+00, .9338E-01, .1105E+00, .1188E+00, 
     1   .5688E-01, .4516E-01, .1011E+01, .1126E+01, .1274E+01, 
     2   .1194E+01, .5876E+00, .4705E+00, .3210E+00, .2489E+00, 
     3   .1574E+00, .1099E+00, .2069E+00, .5297E+00, .1960E+00, 
     4   .9338E-01, .1105E+00, .1188E+00, .5688E-01, .4516E-01, 
     5   .1011E+01, .1126E+01, .1274E+01, .1194E+01, .5876E+00, 
     6   .4705E+00, .3210E+00, .2489E+00, .1574E+00, .1099E+00, 
     7   .2069E+00, .5297E+00, .1960E+00, .9338E-01, .1105E+00, 
     8   .1188E+00, .5688E-01, .4516E-01, .1011E+01, .1126E+01, 
     9   .1274E+01, .1194E+01, .5876E+00, .4705E+00, .3210E+00, 
     .   .2489E+00, .1574E+00, .1099E+00, .2069E+00, .5297E+00, 
     1   .1960E+00, .9338E-01, .1105E+00, .1188E+00, .5688E-01, 
     2   .4516E-01, .1011E+01, .1126E+01, .1274E+01, .1194E+01, 
     3   .5876E+00, .4705E+00, .3210E+00, .2489E+00, .1574E+00, 
     4   .1099E+00, .2069E+00, .5297E+00, .1960E+00, .9338E-01, 
     5   .1105E+00, .1188E+00, .5688E-01, .4516E-01, .1011E+01, 
     6   .1126E+01, .1274E+01, .1194E+01, .5876E+00, .4705E+00, 
     7   .3210E+00, .2489E+00, .1574E+00, .1099E+00, .2069E+00, 
     8   .5297E+00, .1960E+00, .9338E-01, .1105E+00, .1188E+00, 
     9   .5688E-01, .4516E-01/ 
      DATA ((A_ASYX(I,J,5),I=1,mbx),J=1,nrh)/ .7338E+00, .6749E+00, 
     1   .6812E+00, .6876E+00, .6653E+00, .6352E+00, .5506E+00, 
     2   .5123E+00, .4335E+00, .3460E+00, .2780E+00, .2550E+00, 
     3   .2217E+00, .1555E+00, .1096E+00, .9265E-01, .5052E-01, 
     4   .2847E-01, .7338E+00, .6749E+00, .6812E+00, .6876E+00, 
     5   .6653E+00, .6352E+00, .5506E+00, .5123E+00, .4335E+00, 
     6   .3460E+00, .2780E+00, .2550E+00, .2217E+00, .1555E+00, 
     7   .1096E+00, .9265E-01, .5052E-01, .2847E-01, .7338E+00, 
     8   .6749E+00, .6812E+00, .6876E+00, .6653E+00, .6352E+00, 
     9   .5506E+00, .5123E+00, .4335E+00, .3460E+00, .2780E+00, 
     .   .2550E+00, .2217E+00, .1555E+00, .1096E+00, .9265E-01, 
     1   .5052E-01, .2847E-01, .7338E+00, .6749E+00, .6812E+00, 
     2   .6876E+00, .6653E+00, .6352E+00, .5506E+00, .5123E+00, 
     3   .4335E+00, .3460E+00, .2780E+00, .2550E+00, .2217E+00, 
     4   .1555E+00, .1096E+00, .9265E-01, .5052E-01, .2847E-01, 
     5   .7338E+00, .6749E+00, .6812E+00, .6876E+00, .6653E+00, 
     6   .6352E+00, .5506E+00, .5123E+00, .4335E+00, .3460E+00, 
     7   .2780E+00, .2550E+00, .2217E+00, .1555E+00, .1096E+00, 
     8   .9265E-01, .5052E-01, .2847E-01, .7338E+00, .6749E+00, 
     9   .6812E+00, .6876E+00, .6653E+00, .6352E+00, .5506E+00, 
     .   .5123E+00, .4335E+00, .3460E+00, .2780E+00, .2550E+00, 
     1   .2217E+00, .1555E+00, .1096E+00, .9265E-01, .5052E-01, 
     2   .2847E-01, .7338E+00, .6749E+00, .6812E+00, .6876E+00, 
     3   .6653E+00, .6352E+00, .5506E+00, .5123E+00, .4335E+00, 
     4   .3460E+00, .2780E+00, .2550E+00, .2217E+00, .1555E+00, 
     5   .1096E+00, .9265E-01, .5052E-01, .2847E-01, .7338E+00, 
     6   .6749E+00, .6812E+00, .6876E+00, .6653E+00, .6352E+00, 
     7   .5506E+00, .5123E+00, .4335E+00, .3460E+00, .2780E+00, 
     8   .2550E+00, .2217E+00, .1555E+00, .1096E+00, .9265E-01, 
     9   .5052E-01, .2847E-01/ 
      DATA ((A_SSAX(I,J,6),I=1,mbx),J=1,nrh)/ .7767E+00, .8913E+00, 
     1   .9229E+00, .9437E+00, .9070E+00, .9518E+00, .9450E+00, 
     2   .8785E+00, .7097E+00, .5202E+00, .2521E+00, .4713E+00, 
     3   .4974E+00, .4416E+00, .2753E+00, .3267E+00, .2329E+00, 
     4   .1353E+00, .7767E+00, .8913E+00, .9229E+00, .9437E+00, 
     5   .9070E+00, .9518E+00, .9450E+00, .8785E+00, .7097E+00, 
     6   .5202E+00, .2521E+00, .4713E+00, .4974E+00, .4416E+00, 
     7   .2753E+00, .3267E+00, .2329E+00, .1353E+00, .7767E+00, 
     8   .8913E+00, .9229E+00, .9437E+00, .9070E+00, .9518E+00, 
     9   .9450E+00, .8785E+00, .7097E+00, .5202E+00, .2521E+00, 
     .   .4713E+00, .4974E+00, .4416E+00, .2753E+00, .3267E+00, 
     1   .2329E+00, .1353E+00, .7767E+00, .8913E+00, .9229E+00, 
     2   .9437E+00, .9070E+00, .9518E+00, .9450E+00, .8785E+00, 
     3   .7097E+00, .5202E+00, .2521E+00, .4713E+00, .4974E+00, 
     4   .4416E+00, .2753E+00, .3267E+00, .2329E+00, .1353E+00, 
     5   .7767E+00, .8913E+00, .9229E+00, .9437E+00, .9070E+00, 
     6   .9518E+00, .9450E+00, .8785E+00, .7097E+00, .5202E+00, 
     7   .2521E+00, .4713E+00, .4974E+00, .4416E+00, .2753E+00, 
     8   .3267E+00, .2329E+00, .1353E+00, .7767E+00, .8913E+00, 
     9   .9229E+00, .9437E+00, .9070E+00, .9518E+00, .9450E+00, 
     .   .8785E+00, .7097E+00, .5202E+00, .2521E+00, .4713E+00, 
     1   .4974E+00, .4416E+00, .2753E+00, .3267E+00, .2329E+00, 
     2   .1353E+00, .7767E+00, .8913E+00, .9229E+00, .9437E+00, 
     3   .9070E+00, .9518E+00, .9450E+00, .8785E+00, .7097E+00, 
     4   .5202E+00, .2521E+00, .4713E+00, .4974E+00, .4416E+00, 
     5   .2753E+00, .3267E+00, .2329E+00, .1353E+00, .7767E+00, 
     6   .8913E+00, .9229E+00, .9437E+00, .9070E+00, .9518E+00, 
     7   .9450E+00, .8785E+00, .7097E+00, .5202E+00, .2521E+00, 
     8   .4713E+00, .4974E+00, .4416E+00, .2753E+00, .3267E+00, 
     9   .2329E+00, .1353E+00/ 
      DATA ((A_EXTX(I,J,6),I=1,mbx),J=1,nrh)/ .1004E+01, .1058E+01, 
     1   .1170E+01, .1268E+01, .1279E+01, .1229E+01, .1090E+01, 
     2   .9105E+00, .5986E+00, .3776E+00, .4888E+00, .1196E+01, 
     3   .6530E+00, .3654E+00, .3515E+00, .4897E+00, .2131E+00, 
     4   .1327E+00, .1004E+01, .1058E+01, .1170E+01, .1268E+01, 
     5   .1279E+01, .1229E+01, .1090E+01, .9105E+00, .5986E+00, 
     6   .3776E+00, .4888E+00, .1196E+01, .6530E+00, .3654E+00, 
     7   .3515E+00, .4897E+00, .2131E+00, .1327E+00, .1004E+01, 
     8   .1058E+01, .1170E+01, .1268E+01, .1279E+01, .1229E+01, 
     9   .1090E+01, .9105E+00, .5986E+00, .3776E+00, .4888E+00, 
     .   .1196E+01, .6530E+00, .3654E+00, .3515E+00, .4897E+00, 
     1   .2131E+00, .1327E+00, .1004E+01, .1058E+01, .1170E+01, 
     2   .1268E+01, .1279E+01, .1229E+01, .1090E+01, .9105E+00, 
     3   .5986E+00, .3776E+00, .4888E+00, .1196E+01, .6530E+00, 
     4   .3654E+00, .3515E+00, .4897E+00, .2131E+00, .1327E+00, 
     5   .1004E+01, .1058E+01, .1170E+01, .1268E+01, .1279E+01, 
     6   .1229E+01, .1090E+01, .9105E+00, .5986E+00, .3776E+00, 
     7   .4888E+00, .1196E+01, .6530E+00, .3654E+00, .3515E+00, 
     8   .4897E+00, .2131E+00, .1327E+00, .1004E+01, .1058E+01, 
     9   .1170E+01, .1268E+01, .1279E+01, .1229E+01, .1090E+01, 
     .   .9105E+00, .5986E+00, .3776E+00, .4888E+00, .1196E+01, 
     1   .6530E+00, .3654E+00, .3515E+00, .4897E+00, .2131E+00, 
     2   .1327E+00, .1004E+01, .1058E+01, .1170E+01, .1268E+01, 
     3   .1279E+01, .1229E+01, .1090E+01, .9105E+00, .5986E+00, 
     4   .3776E+00, .4888E+00, .1196E+01, .6530E+00, .3654E+00, 
     5   .3515E+00, .4897E+00, .2131E+00, .1327E+00, .1004E+01, 
     6   .1058E+01, .1170E+01, .1268E+01, .1279E+01, .1229E+01, 
     7   .1090E+01, .9105E+00, .5986E+00, .3776E+00, .4888E+00, 
     8   .1196E+01, .6530E+00, .3654E+00, .3515E+00, .4897E+00, 
     9   .2131E+00, .1327E+00/ 
      DATA ((A_ASYX(I,J,6),I=1,mbx),J=1,nrh)/ .8134E+00, .7428E+00, 
     1   .6826E+00, .6685E+00, .7403E+00, .7278E+00, .6859E+00, 
     2   .6902E+00, .6914E+00, .6552E+00, .5806E+00, .3985E+00, 
     3   .4769E+00, .4206E+00, .3269E+00, .2136E+00, .1687E+00, 
     4   .1122E+00, .8134E+00, .7428E+00, .6826E+00, .6685E+00, 
     5   .7403E+00, .7278E+00, .6859E+00, .6902E+00, .6914E+00, 
     6   .6552E+00, .5806E+00, .3985E+00, .4769E+00, .4206E+00, 
     7   .3269E+00, .2136E+00, .1687E+00, .1122E+00, .8134E+00, 
     8   .7428E+00, .6826E+00, .6685E+00, .7403E+00, .7278E+00, 
     9   .6859E+00, .6902E+00, .6914E+00, .6552E+00, .5806E+00, 
     .   .3985E+00, .4769E+00, .4206E+00, .3269E+00, .2136E+00, 
     1   .1687E+00, .1122E+00, .8134E+00, .7428E+00, .6826E+00, 
     2   .6685E+00, .7403E+00, .7278E+00, .6859E+00, .6902E+00, 
     3   .6914E+00, .6552E+00, .5806E+00, .3985E+00, .4769E+00, 
     4   .4206E+00, .3269E+00, .2136E+00, .1687E+00, .1122E+00, 
     5   .8134E+00, .7428E+00, .6826E+00, .6685E+00, .7403E+00, 
     6   .7278E+00, .6859E+00, .6902E+00, .6914E+00, .6552E+00, 
     7   .5806E+00, .3985E+00, .4769E+00, .4206E+00, .3269E+00, 
     8   .2136E+00, .1687E+00, .1122E+00, .8134E+00, .7428E+00, 
     9   .6826E+00, .6685E+00, .7403E+00, .7278E+00, .6859E+00, 
     .   .6902E+00, .6914E+00, .6552E+00, .5806E+00, .3985E+00, 
     1   .4769E+00, .4206E+00, .3269E+00, .2136E+00, .1687E+00, 
     2   .1122E+00, .8134E+00, .7428E+00, .6826E+00, .6685E+00, 
     3   .7403E+00, .7278E+00, .6859E+00, .6902E+00, .6914E+00, 
     4   .6552E+00, .5806E+00, .3985E+00, .4769E+00, .4206E+00, 
     5   .3269E+00, .2136E+00, .1687E+00, .1122E+00, .8134E+00, 
     6   .7428E+00, .6826E+00, .6685E+00, .7403E+00, .7278E+00, 
     7   .6859E+00, .6902E+00, .6914E+00, .6552E+00, .5806E+00, 
     8   .3985E+00, .4769E+00, .4206E+00, .3269E+00, .2136E+00, 
     9   .1687E+00, .1122E+00/ 
      DATA ((A_SSAX(I,J,7),I=1,mbx),J=1,nrh)/ .6979E+00, .8213E+00, 
     1   .8632E+00, .8896E+00, .8303E+00, .9082E+00, .9083E+00, 
     2   .8444E+00, .7304E+00, .6223E+00, .3484E+00, .4968E+00, 
     3   .5537E+00, .5626E+00, .4275E+00, .4462E+00, .4227E+00, 
     4   .3573E+00, .6979E+00, .8213E+00, .8632E+00, .8896E+00, 
     5   .8303E+00, .9082E+00, .9083E+00, .8444E+00, .7304E+00, 
     6   .6223E+00, .3484E+00, .4968E+00, .5537E+00, .5626E+00, 
     7   .4275E+00, .4462E+00, .4227E+00, .3573E+00, .6979E+00, 
     8   .8213E+00, .8632E+00, .8896E+00, .8303E+00, .9082E+00, 
     9   .9083E+00, .8444E+00, .7304E+00, .6223E+00, .3484E+00, 
     .   .4968E+00, .5537E+00, .5626E+00, .4275E+00, .4462E+00, 
     1   .4227E+00, .3573E+00, .6979E+00, .8213E+00, .8632E+00, 
     2   .8896E+00, .8303E+00, .9082E+00, .9083E+00, .8444E+00, 
     3   .7304E+00, .6223E+00, .3484E+00, .4968E+00, .5537E+00, 
     4   .5626E+00, .4275E+00, .4462E+00, .4227E+00, .3573E+00, 
     5   .6979E+00, .8213E+00, .8632E+00, .8896E+00, .8303E+00, 
     6   .9082E+00, .9083E+00, .8444E+00, .7304E+00, .6223E+00, 
     7   .3484E+00, .4968E+00, .5537E+00, .5626E+00, .4275E+00, 
     8   .4462E+00, .4227E+00, .3573E+00, .6979E+00, .8213E+00, 
     9   .8632E+00, .8896E+00, .8303E+00, .9082E+00, .9083E+00, 
     .   .8444E+00, .7304E+00, .6223E+00, .3484E+00, .4968E+00, 
     1   .5537E+00, .5626E+00, .4275E+00, .4462E+00, .4227E+00, 
     2   .3573E+00, .6979E+00, .8213E+00, .8632E+00, .8896E+00, 
     3   .8303E+00, .9082E+00, .9083E+00, .8444E+00, .7304E+00, 
     4   .6223E+00, .3484E+00, .4968E+00, .5537E+00, .5626E+00, 
     5   .4275E+00, .4462E+00, .4227E+00, .3573E+00, .6979E+00, 
     6   .8213E+00, .8632E+00, .8896E+00, .8303E+00, .9082E+00, 
     7   .9083E+00, .8444E+00, .7304E+00, .6223E+00, .3484E+00, 
     8   .4968E+00, .5537E+00, .5626E+00, .4275E+00, .4462E+00, 
     9   .4227E+00, .3573E+00/ 
      DATA ((A_EXTX(I,J,7),I=1,mbx),J=1,nrh)/ .1003E+01, .1034E+01, 
     1   .1088E+01, .1130E+01, .1278E+01, .1325E+01, .1396E+01, 
     2   .1369E+01, .1214E+01, .8716E+00, .7715E+00, .1332E+01, 
     3   .1231E+01, .9736E+00, .8431E+00, .1167E+01, .8008E+00, 
     4   .5356E+00, .1003E+01, .1034E+01, .1088E+01, .1130E+01, 
     5   .1278E+01, .1325E+01, .1396E+01, .1369E+01, .1214E+01, 
     6   .8716E+00, .7715E+00, .1332E+01, .1231E+01, .9736E+00, 
     7   .8431E+00, .1167E+01, .8008E+00, .5356E+00, .1003E+01, 
     8   .1034E+01, .1088E+01, .1130E+01, .1278E+01, .1325E+01, 
     9   .1396E+01, .1369E+01, .1214E+01, .8716E+00, .7715E+00, 
     .   .1332E+01, .1231E+01, .9736E+00, .8431E+00, .1167E+01, 
     1   .8008E+00, .5356E+00, .1003E+01, .1034E+01, .1088E+01, 
     2   .1130E+01, .1278E+01, .1325E+01, .1396E+01, .1369E+01, 
     3   .1214E+01, .8716E+00, .7715E+00, .1332E+01, .1231E+01, 
     4   .9736E+00, .8431E+00, .1167E+01, .8008E+00, .5356E+00, 
     5   .1003E+01, .1034E+01, .1088E+01, .1130E+01, .1278E+01, 
     6   .1325E+01, .1396E+01, .1369E+01, .1214E+01, .8716E+00, 
     7   .7715E+00, .1332E+01, .1231E+01, .9736E+00, .8431E+00, 
     8   .1167E+01, .8008E+00, .5356E+00, .1003E+01, .1034E+01, 
     9   .1088E+01, .1130E+01, .1278E+01, .1325E+01, .1396E+01, 
     .   .1369E+01, .1214E+01, .8716E+00, .7715E+00, .1332E+01, 
     1   .1231E+01, .9736E+00, .8431E+00, .1167E+01, .8008E+00, 
     2   .5356E+00, .1003E+01, .1034E+01, .1088E+01, .1130E+01, 
     3   .1278E+01, .1325E+01, .1396E+01, .1369E+01, .1214E+01, 
     4   .8716E+00, .7715E+00, .1332E+01, .1231E+01, .9736E+00, 
     5   .8431E+00, .1167E+01, .8008E+00, .5356E+00, .1003E+01, 
     6   .1034E+01, .1088E+01, .1130E+01, .1278E+01, .1325E+01, 
     7   .1396E+01, .1369E+01, .1214E+01, .8716E+00, .7715E+00, 
     8   .1332E+01, .1231E+01, .9736E+00, .8431E+00, .1167E+01, 
     9   .8008E+00, .5356E+00/ 
      DATA ((A_ASYX(I,J,7),I=1,mbx),J=1,nrh)/ .8694E+00, .8106E+00, 
     1   .7632E+00, .7289E+00, .7415E+00, .7160E+00, .6904E+00, 
     2   .7315E+00, .8001E+00, .8185E+00, .7903E+00, .5996E+00, 
     3   .6566E+00, .6460E+00, .5866E+00, .3647E+00, .3224E+00, 
     4   .2752E+00, .8694E+00, .8106E+00, .7632E+00, .7289E+00, 
     5   .7415E+00, .7160E+00, .6904E+00, .7315E+00, .8001E+00, 
     6   .8185E+00, .7903E+00, .5996E+00, .6566E+00, .6460E+00, 
     7   .5866E+00, .3647E+00, .3224E+00, .2752E+00, .8694E+00, 
     8   .8106E+00, .7632E+00, .7289E+00, .7415E+00, .7160E+00, 
     9   .6904E+00, .7315E+00, .8001E+00, .8185E+00, .7903E+00, 
     .   .5996E+00, .6566E+00, .6460E+00, .5866E+00, .3647E+00, 
     1   .3224E+00, .2752E+00, .8694E+00, .8106E+00, .7632E+00, 
     2   .7289E+00, .7415E+00, .7160E+00, .6904E+00, .7315E+00, 
     3   .8001E+00, .8185E+00, .7903E+00, .5996E+00, .6566E+00, 
     4   .6460E+00, .5866E+00, .3647E+00, .3224E+00, .2752E+00, 
     5   .8694E+00, .8106E+00, .7632E+00, .7289E+00, .7415E+00, 
     6   .7160E+00, .6904E+00, .7315E+00, .8001E+00, .8185E+00, 
     7   .7903E+00, .5996E+00, .6566E+00, .6460E+00, .5866E+00, 
     8   .3647E+00, .3224E+00, .2752E+00, .8694E+00, .8106E+00, 
     9   .7632E+00, .7289E+00, .7415E+00, .7160E+00, .6904E+00, 
     .   .7315E+00, .8001E+00, .8185E+00, .7903E+00, .5996E+00, 
     1   .6566E+00, .6460E+00, .5866E+00, .3647E+00, .3224E+00, 
     2   .2752E+00, .8694E+00, .8106E+00, .7632E+00, .7289E+00, 
     3   .7415E+00, .7160E+00, .6904E+00, .7315E+00, .8001E+00, 
     4   .8185E+00, .7903E+00, .5996E+00, .6566E+00, .6460E+00, 
     5   .5866E+00, .3647E+00, .3224E+00, .2752E+00, .8694E+00, 
     6   .8106E+00, .7632E+00, .7289E+00, .7415E+00, .7160E+00, 
     7   .6904E+00, .7315E+00, .8001E+00, .8185E+00, .7903E+00, 
     8   .5996E+00, .6566E+00, .6460E+00, .5866E+00, .3647E+00, 
     9   .3224E+00, .2752E+00/ 
      DATA ((A_SSAX(I,J,8),I=1,mbx),J=1,nrh)/ .6279E+00, .7298E+00, 
     1   .7835E+00, .8196E+00, .7267E+00, .8274E+00, .8228E+00, 
     2   .7316E+00, .6350E+00, .6091E+00, .4227E+00, .5355E+00, 
     3   .5210E+00, .5494E+00, .4857E+00, .4905E+00, .4786E+00, 
     4   .4606E+00, .6279E+00, .7298E+00, .7835E+00, .8196E+00, 
     5   .7267E+00, .8274E+00, .8228E+00, .7316E+00, .6350E+00, 
     6   .6091E+00, .4227E+00, .5355E+00, .5210E+00, .5494E+00, 
     7   .4857E+00, .4905E+00, .4786E+00, .4606E+00, .6279E+00, 
     8   .7298E+00, .7835E+00, .8196E+00, .7267E+00, .8274E+00, 
     9   .8228E+00, .7316E+00, .6350E+00, .6091E+00, .4227E+00, 
     .   .5355E+00, .5210E+00, .5494E+00, .4857E+00, .4905E+00, 
     1   .4786E+00, .4606E+00, .6279E+00, .7298E+00, .7835E+00, 
     2   .8196E+00, .7267E+00, .8274E+00, .8228E+00, .7316E+00, 
     3   .6350E+00, .6091E+00, .4227E+00, .5355E+00, .5210E+00, 
     4   .5494E+00, .4857E+00, .4905E+00, .4786E+00, .4606E+00, 
     5   .6279E+00, .7298E+00, .7835E+00, .8196E+00, .7267E+00, 
     6   .8274E+00, .8228E+00, .7316E+00, .6350E+00, .6091E+00, 
     7   .4227E+00, .5355E+00, .5210E+00, .5494E+00, .4857E+00, 
     8   .4905E+00, .4786E+00, .4606E+00, .6279E+00, .7298E+00, 
     9   .7835E+00, .8196E+00, .7267E+00, .8274E+00, .8228E+00, 
     .   .7316E+00, .6350E+00, .6091E+00, .4227E+00, .5355E+00, 
     1   .5210E+00, .5494E+00, .4857E+00, .4905E+00, .4786E+00, 
     2   .4606E+00, .6279E+00, .7298E+00, .7835E+00, .8196E+00, 
     3   .7267E+00, .8274E+00, .8228E+00, .7316E+00, .6350E+00, 
     4   .6091E+00, .4227E+00, .5355E+00, .5210E+00, .5494E+00, 
     5   .4857E+00, .4905E+00, .4786E+00, .4606E+00, .6279E+00, 
     6   .7298E+00, .7835E+00, .8196E+00, .7267E+00, .8274E+00, 
     7   .8228E+00, .7316E+00, .6350E+00, .6091E+00, .4227E+00, 
     8   .5355E+00, .5210E+00, .5494E+00, .4857E+00, .4905E+00, 
     9   .4786E+00, .4606E+00/ 
      DATA ((A_EXTX(I,J,8),I=1,mbx),J=1,nrh)/ .1002E+01, .1022E+01, 
     1   .1054E+01, .1076E+01, .1139E+01, .1159E+01, .1207E+01, 
     2   .1239E+01, .1280E+01, .1182E+01, .9419E+00, .1253E+01, 
     3   .1321E+01, .1320E+01, .1211E+01, .1372E+01, .1347E+01, 
     4   .1223E+01, .1002E+01, .1022E+01, .1054E+01, .1076E+01, 
     5   .1139E+01, .1159E+01, .1207E+01, .1239E+01, .1280E+01, 
     6   .1182E+01, .9419E+00, .1253E+01, .1321E+01, .1320E+01, 
     7   .1211E+01, .1372E+01, .1347E+01, .1223E+01, .1002E+01, 
     8   .1022E+01, .1054E+01, .1076E+01, .1139E+01, .1159E+01, 
     9   .1207E+01, .1239E+01, .1280E+01, .1182E+01, .9419E+00, 
     .   .1253E+01, .1321E+01, .1320E+01, .1211E+01, .1372E+01, 
     1   .1347E+01, .1223E+01, .1002E+01, .1022E+01, .1054E+01, 
     2   .1076E+01, .1139E+01, .1159E+01, .1207E+01, .1239E+01, 
     3   .1280E+01, .1182E+01, .9419E+00, .1253E+01, .1321E+01, 
     4   .1320E+01, .1211E+01, .1372E+01, .1347E+01, .1223E+01, 
     5   .1002E+01, .1022E+01, .1054E+01, .1076E+01, .1139E+01, 
     6   .1159E+01, .1207E+01, .1239E+01, .1280E+01, .1182E+01, 
     7   .9419E+00, .1253E+01, .1321E+01, .1320E+01, .1211E+01, 
     8   .1372E+01, .1347E+01, .1223E+01, .1002E+01, .1022E+01, 
     9   .1054E+01, .1076E+01, .1139E+01, .1159E+01, .1207E+01, 
     .   .1239E+01, .1280E+01, .1182E+01, .9419E+00, .1253E+01, 
     1   .1321E+01, .1320E+01, .1211E+01, .1372E+01, .1347E+01, 
     2   .1223E+01, .1002E+01, .1022E+01, .1054E+01, .1076E+01, 
     3   .1139E+01, .1159E+01, .1207E+01, .1239E+01, .1280E+01, 
     4   .1182E+01, .9419E+00, .1253E+01, .1321E+01, .1320E+01, 
     5   .1211E+01, .1372E+01, .1347E+01, .1223E+01, .1002E+01, 
     6   .1022E+01, .1054E+01, .1076E+01, .1139E+01, .1159E+01, 
     7   .1207E+01, .1239E+01, .1280E+01, .1182E+01, .9419E+00, 
     8   .1253E+01, .1321E+01, .1320E+01, .1211E+01, .1372E+01, 
     9   .1347E+01, .1223E+01/ 
      DATA ((A_ASYX(I,J,8),I=1,mbx),J=1,nrh)/ .9078E+00, .8641E+00, 
     1   .8296E+00, .8045E+00, .8136E+00, .7683E+00, .7318E+00, 
     2   .7706E+00, .8403E+00, .8819E+00, .8892E+00, .7333E+00, 
     3   .7707E+00, .7667E+00, .7564E+00, .5706E+00, .4955E+00, 
     4   .4467E+00, .9078E+00, .8641E+00, .8296E+00, .8045E+00, 
     5   .8136E+00, .7683E+00, .7318E+00, .7706E+00, .8403E+00, 
     6   .8819E+00, .8892E+00, .7333E+00, .7707E+00, .7667E+00, 
     7   .7564E+00, .5706E+00, .4955E+00, .4467E+00, .9078E+00, 
     8   .8641E+00, .8296E+00, .8045E+00, .8136E+00, .7683E+00, 
     9   .7318E+00, .7706E+00, .8403E+00, .8819E+00, .8892E+00, 
     .   .7333E+00, .7707E+00, .7667E+00, .7564E+00, .5706E+00, 
     1   .4955E+00, .4467E+00, .9078E+00, .8641E+00, .8296E+00, 
     2   .8045E+00, .8136E+00, .7683E+00, .7318E+00, .7706E+00, 
     3   .8403E+00, .8819E+00, .8892E+00, .7333E+00, .7707E+00, 
     4   .7667E+00, .7564E+00, .5706E+00, .4955E+00, .4467E+00, 
     5   .9078E+00, .8641E+00, .8296E+00, .8045E+00, .8136E+00, 
     6   .7683E+00, .7318E+00, .7706E+00, .8403E+00, .8819E+00, 
     7   .8892E+00, .7333E+00, .7707E+00, .7667E+00, .7564E+00, 
     8   .5706E+00, .4955E+00, .4467E+00, .9078E+00, .8641E+00, 
     9   .8296E+00, .8045E+00, .8136E+00, .7683E+00, .7318E+00, 
     .   .7706E+00, .8403E+00, .8819E+00, .8892E+00, .7333E+00, 
     1   .7707E+00, .7667E+00, .7564E+00, .5706E+00, .4955E+00, 
     2   .4467E+00, .9078E+00, .8641E+00, .8296E+00, .8045E+00, 
     3   .8136E+00, .7683E+00, .7318E+00, .7706E+00, .8403E+00, 
     4   .8819E+00, .8892E+00, .7333E+00, .7707E+00, .7667E+00, 
     5   .7564E+00, .5706E+00, .4955E+00, .4467E+00, .9078E+00, 
     6   .8641E+00, .8296E+00, .8045E+00, .8136E+00, .7683E+00, 
     7   .7318E+00, .7706E+00, .8403E+00, .8819E+00, .8892E+00, 
     8   .7333E+00, .7707E+00, .7667E+00, .7564E+00, .5706E+00, 
     9   .4955E+00, .4467E+00/ 
      DATA ((A_EXTX(I,J,9),I=1,mbx),J=1,1)/ 0.9992E+00, 0.1055E+01, 
     1   0.1097E+01, 0.9565E+00, 0.7209E+00, 0.8266E+00, 0.6757E+00, 
     2   0.4984E+00, 0.4294E+00, 0.4649E+00, 0.5541E+00, 0.8549E+00, 
     3   0.6774E+00, 0.5136E+00, 0.4909E+00, 0.4952E+00, 0.4213E+00, 
     4   0.3563E+00/ 
      DATA ((A_SSAX(I,J,9),I=1,mbx),J=1,1)/ 0.7289E+00, 0.7933E+00, 
     1   0.8553E+00, 0.8828E+00, 0.8465E+00, 0.8840E+00, 0.8537E+00, 
     2   0.7561E+00, 0.5914E+00, 0.6595E+00, 0.5205E+00, 0.5811E+00, 
     3   0.6361E+00, 0.6307E+00, 0.6348E+00, 0.5020E+00, 0.4057E+00, 
     4   0.3352E+00/ 
      DATA ((A_ASYX(I,J,9),I=1,mbx),J=1,1)/ 0.8317E+00, 0.7882E+00, 
     1   0.8003E+00, 0.8834E+00, 0.9145E+00, 0.8506E+00, 0.8563E+00, 
     2   0.8778E+00, 0.8615E+00, 0.8283E+00, 0.7892E+00, 0.6657E+00, 
     3   0.6808E+00, 0.6886E+00, 0.6387E+00, 0.5706E+00, 0.4973E+00, 
     4   0.3480E+00/ 
      DATA ((A_EXTX(I,J,10),I=1,mbx),J=1,8)/ 0.1015E+01, 0.4304E+00, 
     1   0.1407E+00, 0.4076E-01, 0.2793E-01, 0.9580E-02, 0.8208E-02, 
     2   0.8256E-02, 0.1281E-01, 0.1857E-01, 0.3590E-01, 0.2404E-01, 
     3   0.8152E-02, 0.8273E-02, 0.2205E-01, 0.1285E-01, 0.1426E-01, 
     4   0.2008E-01, 0.1015E+01, 0.4417E+00, 0.1461E+00, 0.4302E-01, 
     5   0.7731E-01, 0.1407E-01, 0.1064E-01, 0.1307E-01, 0.1784E-01, 
     6   0.1697E-01, 0.3190E-01, 0.2466E-01, 0.1843E-01, 0.3242E-01, 
     7   0.3835E-01, 0.2407E-01, 0.1947E-01, 0.2221E-01, 0.1014E+01, 
     8   0.4488E+00, 0.1501E+00, 0.4563E-01, 0.9123E-01, 0.1620E-01, 
     9   0.1167E-01, 0.1446E-01, 0.1949E-01, 0.1647E-01, 0.2965E-01, 
     .   0.2378E-01, 0.2184E-01, 0.4049E-01, 0.4408E-01, 0.2782E-01, 
     1   0.2101E-01, 0.2265E-01, 0.1014E+01, 0.4559E+00, 0.1546E+00, 
     2   0.4833E-01, 0.1014E+00, 0.1818E-01, 0.1258E-01, 0.1550E-01, 
     3   0.2072E-01, 0.1613E-01, 0.2769E-01, 0.2283E-01, 0.2429E-01, 
     4   0.4633E-01, 0.4824E-01, 0.3054E-01, 0.2207E-01, 0.2287E-01, 
     5   0.1013E+01, 0.4713E+00, 0.1650E+00, 0.5434E-01, 0.1170E+00, 
     6   0.2223E-01, 0.1438E-01, 0.1721E-01, 0.2264E-01, 0.1573E-01, 
     7   0.2436E-01, 0.2105E-01, 0.2776E-01, 0.5469E-01, 0.5422E-01, 
     8   0.3446E-01, 0.2350E-01, 0.2308E-01, 0.1012E+01, 0.4913E+00, 
     9   0.1795E+00, 0.6242E-01, 0.1310E+00, 0.2735E-01, 0.1663E-01, 
     .   0.1898E-01, 0.2444E-01, 0.1565E-01, 0.2154E-01, 0.1943E-01, 
     1   0.3035E-01, 0.6100E-01, 0.5878E-01, 0.3746E-01, 0.2454E-01, 
     2   0.2318E-01, 0.1011E+01, 0.5221E+00, 0.2038E+00, 0.7577E-01, 
     3   0.1474E+00, 0.3559E-01, 0.2032E-01, 0.2156E-01, 0.2687E-01, 
     4   0.1612E-01, 0.1925E-01, 0.1809E-01, 0.3271E-01, 0.6668E-01, 
     5   0.6308E-01, 0.4032E-01, 0.2557E-01, 0.2331E-01, 0.1010E+01, 
     6   0.5447E+00, 0.2230E+00, 0.8651E-01, 0.1583E+00, 0.4224E-01, 
     7   0.2339E-01, 0.2363E-01, 0.2874E-01, 0.1684E-01, 0.1855E-01, 
     8   0.1772E-01, 0.3399E-01, 0.6967E-01, 0.6553E-01, 0.4195E-01, 
     9   0.2624E-01, 0.2355E-01/ 
      DATA ((A_SSAX(I,J,10),I=1,mbx),J=1,8)/ 0.9633E+00, 0.8961E+00, 
     1   0.7687E+00, 0.7940E+00, 0.5192E+00, 0.7595E+00, 0.3996E+00, 
     2   0.2073E+00, 0.9201E-01, 0.1337E-01, 0.3585E-01, 0.4407E-01, 
     3   0.3943E-01, 0.1419E-01, 0.4977E-02, 0.3965E-02, 0.8396E-03, 
     4   0.1067E-03, 0.9776E+00, 0.9357E+00, 0.8539E+00, 0.8702E+00, 
     5   0.4324E+00, 0.7956E+00, 0.4373E+00, 0.2277E+00, 0.9151E-01, 
     6   0.3219E-01, 0.3617E-01, 0.4094E-01, 0.1815E-01, 0.4240E-02, 
     7   0.3668E-02, 0.3164E-02, 0.1059E-02, 0.1485E-03, 0.9820E+00, 
     8   0.9484E+00, 0.8829E+00, 0.8940E+00, 0.4503E+00, 0.8121E+00, 
     9   0.4596E+00, 0.2482E+00, 0.9977E-01, 0.4317E-01, 0.3868E-01, 
     .   0.4120E-01, 0.1566E-01, 0.3913E-02, 0.3780E-02, 0.3268E-02, 
     1   0.1219E-02, 0.1749E-03, 0.9850E+00, 0.9577E+00, 0.9041E+00, 
     2   0.9111E+00, 0.4686E+00, 0.8259E+00, 0.4808E+00, 0.2680E+00, 
     3   0.1088E+00, 0.5411E-01, 0.4245E-01, 0.4209E-01, 0.1441E-01, 
     4   0.3948E-02, 0.4031E-02, 0.3475E-02, 0.1394E-02, 0.2031E-03, 
     5   0.9895E+00, 0.9710E+00, 0.9349E+00, 0.9351E+00, 0.5035E+00, 
     6   0.8496E+00, 0.5230E+00, 0.3079E+00, 0.1297E+00, 0.7788E-01, 
     7   0.5266E-01, 0.4537E-01, 0.1346E-01, 0.4473E-02, 0.4808E-02, 
     8   0.4116E-02, 0.1823E-02, 0.2823E-03, 0.9928E+00, 0.9810E+00, 
     9   0.9578E+00, 0.9527E+00, 0.5393E+00, 0.8720E+00, 0.5707E+00, 
     .   0.3541E+00, 0.1579E+00, 0.1089E+00, 0.6886E-01, 0.5146E-01, 
     1   0.1377E-01, 0.5552E-02, 0.6102E-02, 0.5203E-02, 0.2488E-02, 
     2   0.4103E-03, 0.9956E+00, 0.9891E+00, 0.9761E+00, 0.9670E+00, 
     3   0.5805E+00, 0.8959E+00, 0.6318E+00, 0.4163E+00, 0.2022E+00, 
     4   0.1575E+00, 0.9848E-01, 0.6388E-01, 0.1577E-01, 0.7723E-02, 
     5   0.8602E-02, 0.7370E-02, 0.3795E-02, 0.6676E-03, 0.9968E+00, 
     6   0.9924E+00, 0.9834E+00, 0.9730E+00, 0.6038E+00, 0.9087E+00, 
     7   0.6692E+00, 0.4563E+00, 0.2352E+00, 0.1938E+00, 0.1231E+00, 
     8   0.7513E-01, 0.1802E-01, 0.9667E-02, 0.1084E-01, 0.9370E-02, 
     9   0.5020E-02, 0.9152E-03/ 
      DATA ((A_ASYX(I,J,10),I=1,mbx),J=1,8)/ 0.6143E+00, 0.5585E+00, 
     1   0.4813E+00, 0.4255E+00, 0.3592E+00, 0.3090E+00, 0.2571E+00, 
     2   0.2231E+00, 0.1912E+00, 0.1539E+00, 0.1519E+00, 0.1543E+00, 
     3   0.1251E+00, 0.9279E-01, 0.7617E-01, 0.6363E-01, 0.3473E-01, 
     4   0.1405E-01, 0.6722E+00, 0.6148E+00, 0.5341E+00, 0.4799E+00, 
     5   0.3881E+00, 0.3523E+00, 0.2969E+00, 0.2574E+00, 0.2247E+00, 
     6   0.1905E+00, 0.1774E+00, 0.1667E+00, 0.1328E+00, 0.9845E-01, 
     7   0.8187E-01, 0.6865E-01, 0.4277E-01, 0.1880E-01, 0.6904E+00, 
     8   0.6342E+00, 0.5549E+00, 0.5011E+00, 0.4022E+00, 0.3708E+00, 
     9   0.3143E+00, 0.2736E+00, 0.2399E+00, 0.2062E+00, 0.1893E+00, 
     .   0.1753E+00, 0.1392E+00, 0.1030E+00, 0.8623E-01, 0.7257E-01, 
     1   0.4688E-01, 0.2120E-01, 0.7042E+00, 0.6494E+00, 0.5709E+00, 
     2   0.5194E+00, 0.4152E+00, 0.3860E+00, 0.3294E+00, 0.2869E+00, 
     3   0.2526E+00, 0.2205E+00, 0.2002E+00, 0.1837E+00, 0.1454E+00, 
     4   0.1079E+00, 0.9064E-01, 0.7653E-01, 0.5081E-01, 0.2344E-01, 
     5   0.7254E+00, 0.6740E+00, 0.5991E+00, 0.5504E+00, 0.4394E+00, 
     6   0.4146E+00, 0.3578E+00, 0.3133E+00, 0.2782E+00, 0.2459E+00, 
     7   0.2224E+00, 0.2017E+00, 0.1597E+00, 0.1181E+00, 0.1001E+00, 
     8   0.8510E-01, 0.5879E-01, 0.2824E-01, 0.7433E+00, 0.6967E+00, 
     9   0.6270E+00, 0.5823E+00, 0.4667E+00, 0.4446E+00, 0.3888E+00, 
     .   0.3432E+00, 0.3068E+00, 0.2755E+00, 0.2484E+00, 0.2241E+00, 
     1   0.1777E+00, 0.1315E+00, 0.1119E+00, 0.9622E-01, 0.6871E-01, 
     2   0.3408E-01, 0.7612E+00, 0.7217E+00, 0.6594E+00, 0.6202E+00, 
     3   0.5032E+00, 0.4833E+00, 0.4295E+00, 0.3836E+00, 0.3454E+00, 
     4   0.3147E+00, 0.2852E+00, 0.2571E+00, 0.2045E+00, 0.1520E+00, 
     5   0.1300E+00, 0.1126E+00, 0.8357E-01, 0.4311E-01, 0.7701E+00, 
     6   0.7355E+00, 0.6784E+00, 0.6433E+00, 0.5268E+00, 0.5073E+00, 
     7   0.4559E+00, 0.4108E+00, 0.3712E+00, 0.3412E+00, 0.3106E+00, 
     8   0.2803E+00, 0.2242E+00, 0.1672E+00, 0.1436E+00, 0.1248E+00, 
     9   0.9454E-01, 0.5018E-01/ 
      DATA ((A_EXTX(I,J,11),I=1,mbx),J=1,1)/ 0.1017E+01, 0.5114E+00, 
     1   0.2718E+00, 0.1911E+00, 0.1446E+00, 0.1113E+00, 0.8555E-01, 
     2   0.7200E-01, 0.6089E-01, 0.5213E-01, 0.4566E-01, 0.4003E-01, 
     3   0.3392E-01, 0.2769E-01, 0.2254E-01, 0.1702E-01, 0.1213E-01, 
     4   0.7093E-02/ 
      DATA ((A_SSAX(I,J,11),I=1,mbx),J=1,1)/ 0.2102E+00, 0.1127E+00, 
     1   0.4250E-01, 0.2007E-01, 0.9655E-02, 0.5070E-02, 0.2738E-02, 
     2   0.1795E-02, 0.1192E-02, 0.8190E-03, 0.6116E-03, 0.4481E-03, 
     3   0.2918E-03, 0.1761E-03, 0.1049E-03, 0.5490E-04, 0.2308E-04, 
     4   0.5530E-05/ 
      DATA ((A_ASYX(I,J,11),I=1,mbx),J=1,1)/ 0.3375E+00, 0.2412E+00, 
     1   0.1541E+00, 0.1086E+00, 0.7644E-01, 0.5501E-01, 0.3917E-01, 
     2   0.3088E-01, 0.2415E-01, 0.1913E-01, 0.1581E-01, 0.1289E-01, 
     3   0.9727E-02, 0.6915E-02, 0.4858E-02, 0.3121E-02, 0.1759E-02, 
     4   0.1141E-02/ 
      DATA ((A_EXTX(I,J,12),I=1,mbx),J=1,8)/ 0.9977E+00, 0.9420E+00, 
     1   0.7044E+00, 0.4678E+00, 0.4148E+00, 0.2336E+00, 0.1493E+00, 
     2   0.7957E-01, 0.7247E-01, 0.4367E-01, 0.5551E-01, 0.4474E-01, 
     3   0.2241E-01, 0.1365E-01, 0.2654E-01, 0.2912E-01, 0.3158E-01, 
     4   0.1472E+00, 0.9989E+00, 0.1010E+01, 0.8520E+00, 0.6085E+00, 
     5   0.5675E+00, 0.4096E+00, 0.2546E+00, 0.1694E+00, 0.1644E+00, 
     6   0.1067E+00, 0.9235E-01, 0.7284E-01, 0.8140E-01, 0.1443E+00, 
     7   0.1495E+00, 0.1081E+00, 0.7115E-01, 0.1054E+00, 0.9984E+00, 
     8   0.1028E+01, 0.9029E+00, 0.6677E+00, 0.6213E+00, 0.4742E+00, 
     9   0.2993E+00, 0.2035E+00, 0.1981E+00, 0.1314E+00, 0.1092E+00, 
     .   0.8514E-01, 0.9789E-01, 0.1765E+00, 0.1827E+00, 0.1321E+00, 
     1   0.8473E-01, 0.1036E+00, 0.9991E+00, 0.1045E+01, 0.9466E+00, 
     2   0.7224E+00, 0.6693E+00, 0.5337E+00, 0.3425E+00, 0.2364E+00, 
     3   0.2303E+00, 0.1556E+00, 0.1262E+00, 0.9751E-01, 0.1125E+00, 
     4   0.2040E+00, 0.2118E+00, 0.1541E+00, 0.9772E-01, 0.1055E+00, 
     5   0.9984E+00, 0.1066E+01, 0.1022E+01, 0.8298E+00, 0.7614E+00, 
     6   0.6542E+00, 0.4369E+00, 0.3089E+00, 0.3009E+00, 0.2109E+00, 
     7   0.1665E+00, 0.1265E+00, 0.1427E+00, 0.2574E+00, 0.2702E+00, 
     8   0.2006E+00, 0.1267E+00, 0.1155E+00, 0.9990E+00, 0.1079E+01, 
     9   0.1088E+01, 0.9476E+00, 0.8625E+00, 0.7960E+00, 0.5610E+00, 
     .   0.4073E+00, 0.3966E+00, 0.2898E+00, 0.2262E+00, 0.1696E+00, 
     1   0.1816E+00, 0.3214E+00, 0.3426E+00, 0.2624E+00, 0.1681E+00, 
     2   0.1363E+00, 0.9994E+00, 0.1075E+01, 0.1141E+01, 0.1086E+01, 
     3   0.9872E+00, 0.9866E+00, 0.7578E+00, 0.5745E+00, 0.5599E+00, 
     4   0.4348E+00, 0.3419E+00, 0.2551E+00, 0.2490E+00, 0.4208E+00, 
     5   0.4585E+00, 0.3704E+00, 0.2469E+00, 0.1830E+00, 0.9998E+00, 
     6   0.1066E+01, 0.1149E+01, 0.1154E+01, 0.1059E+01, 0.1105E+01, 
     7   0.9090E+00, 0.7161E+00, 0.7000E+00, 0.5694E+00, 0.4561E+00, 
     8   0.3430E+00, 0.3119E+00, 0.5031E+00, 0.5565E+00, 0.4695E+00, 
     9   0.3264E+00, 0.2343E+00/ 
      DATA ((A_SSAX(I,J,12),I=1,mbx),J=1,8)/ 0.1000E+01, 0.9991E+00, 
     1   0.9957E+00, 0.9892E+00, 0.9560E+00, 0.9897E+00, 0.9825E+00, 
     2   0.9471E+00, 0.9222E+00, 0.8519E+00, 0.7555E+00, 0.8130E+00, 
     3   0.7537E+00, 0.5286E+00, 0.2979E+00, 0.1749E+00, 0.6170E-01, 
     4   0.1383E-01, 0.1000E+01, 0.9997E+00, 0.9976E+00, 0.9895E+00, 
     5   0.7743E+00, 0.9725E+00, 0.9208E+00, 0.8068E+00, 0.6949E+00, 
     6   0.7072E+00, 0.6425E+00, 0.5446E+00, 0.2100E+00, 0.9825E-01, 
     7   0.1218E+00, 0.1271E+00, 0.8758E-01, 0.1948E-01, 0.1000E+01, 
     8   0.9997E+00, 0.9978E+00, 0.9891E+00, 0.7705E+00, 0.9709E+00, 
     9   0.9168E+00, 0.8044E+00, 0.6898E+00, 0.7093E+00, 0.6428E+00, 
     .   0.5309E+00, 0.1996E+00, 0.1070E+00, 0.1312E+00, 0.1383E+00, 
     1   0.1012E+00, 0.2443E-01, 0.1000E+01, 0.9998E+00, 0.9979E+00, 
     2   0.9888E+00, 0.7694E+00, 0.9699E+00, 0.9150E+00, 0.8052E+00, 
     3   0.6901E+00, 0.7145E+00, 0.6474E+00, 0.5283E+00, 0.1984E+00, 
     4   0.1163E+00, 0.1411E+00, 0.1495E+00, 0.1139E+00, 0.2946E-01, 
     5   0.1000E+01, 0.9999E+00, 0.9979E+00, 0.9880E+00, 0.7688E+00, 
     6   0.9681E+00, 0.9134E+00, 0.8094E+00, 0.6964E+00, 0.7279E+00, 
     7   0.6622E+00, 0.5372E+00, 0.2072E+00, 0.1367E+00, 0.1631E+00, 
     8   0.1741E+00, 0.1408E+00, 0.4127E-01, 0.1000E+01, 0.9999E+00, 
     9   0.9978E+00, 0.9865E+00, 0.7684E+00, 0.9659E+00, 0.9123E+00, 
     .   0.8156E+00, 0.7064E+00, 0.7438E+00, 0.6826E+00, 0.5588E+00, 
     1   0.2270E+00, 0.1626E+00, 0.1908E+00, 0.2048E+00, 0.1747E+00, 
     2   0.5801E-01, 0.1000E+01, 0.9998E+00, 0.9975E+00, 0.9835E+00, 
     3   0.7657E+00, 0.9612E+00, 0.9093E+00, 0.8228E+00, 0.7186E+00, 
     4   0.7624E+00, 0.7101E+00, 0.5961E+00, 0.2647E+00, 0.2028E+00, 
     5   0.2331E+00, 0.2512E+00, 0.2267E+00, 0.9154E-01, 0.1000E+01, 
     6   0.9998E+00, 0.9970E+00, 0.9801E+00, 0.7611E+00, 0.9558E+00, 
     7   0.9043E+00, 0.8251E+00, 0.7236E+00, 0.7714E+00, 0.7266E+00, 
     8   0.6232E+00, 0.2972E+00, 0.2351E+00, 0.2662E+00, 0.2865E+00, 
     9   0.2670E+00, 0.1236E+00/ 
      DATA ((A_ASYX(I,J,12),I=1,mbx),J=1,8)/ 0.6925E+00, 0.7030E+00, 
     1   0.7037E+00, 0.7018E+00, 0.6290E+00, 0.6210E+00, 0.5823E+00, 
     2   0.5754E+00, 0.5304E+00, 0.5025E+00, 0.4631E+00, 0.4344E+00, 
     3   0.4025E+00, 0.3539E+00, 0.3069E+00, 0.2526E+00, 0.1773E+00, 
     4   0.5475E-01, 0.7710E+00, 0.7780E+00, 0.7844E+00, 0.7895E+00, 
     5   0.7592E+00, 0.7110E+00, 0.6965E+00, 0.6880E+00, 0.6458E+00, 
     6   0.6223E+00, 0.5924E+00, 0.5645E+00, 0.5042E+00, 0.4073E+00, 
     7   0.3515E+00, 0.3059E+00, 0.2449E+00, 0.1281E+00, 0.7783E+00, 
     8   0.7853E+00, 0.7928E+00, 0.8012E+00, 0.7771E+00, 0.7257E+00, 
     9   0.7169E+00, 0.7105E+00, 0.6705E+00, 0.6489E+00, 0.6225E+00, 
     .   0.5962E+00, 0.5341E+00, 0.4315E+00, 0.3729E+00, 0.3270E+00, 
     1   0.2664E+00, 0.1485E+00, 0.7840E+00, 0.7886E+00, 0.7979E+00, 
     2   0.8088E+00, 0.7893E+00, 0.7363E+00, 0.7324E+00, 0.7284E+00, 
     3   0.6893E+00, 0.6702E+00, 0.6456E+00, 0.6215E+00, 0.5591E+00, 
     4   0.4524E+00, 0.3921E+00, 0.3454E+00, 0.2847E+00, 0.1658E+00, 
     5   0.7933E+00, 0.7934E+00, 0.8035E+00, 0.8182E+00, 0.8105E+00, 
     6   0.7519E+00, 0.7564E+00, 0.7580E+00, 0.7218E+00, 0.7047E+00, 
     7   0.6859E+00, 0.6656E+00, 0.6055E+00, 0.4947E+00, 0.4303E+00, 
     8   0.3819E+00, 0.3207E+00, 0.2002E+00, 0.8009E+00, 0.7966E+00, 
     9   0.8065E+00, 0.8250E+00, 0.8289E+00, 0.7651E+00, 0.7777E+00, 
     .   0.7851E+00, 0.7526E+00, 0.7388E+00, 0.7254E+00, 0.7104E+00, 
     1   0.6555E+00, 0.5434E+00, 0.4756E+00, 0.4249E+00, 0.3629E+00, 
     2   0.2389E+00, 0.8136E+00, 0.8017E+00, 0.8068E+00, 0.8283E+00, 
     3   0.8495E+00, 0.7764E+00, 0.7991E+00, 0.8159E+00, 0.7883E+00, 
     4   0.7786E+00, 0.7723E+00, 0.7651E+00, 0.7221E+00, 0.6129E+00, 
     5   0.5424E+00, 0.4883E+00, 0.4249E+00, 0.2973E+00, 0.8244E+00, 
     6   0.8072E+00, 0.8066E+00, 0.8277E+00, 0.8617E+00, 0.7812E+00, 
     7   0.8112E+00, 0.8339E+00, 0.8106E+00, 0.8033E+00, 0.8011E+00, 
     8   0.7995E+00, 0.7670E+00, 0.6642E+00, 0.5932E+00, 0.5376E+00, 
     9   0.4736E+00, 0.3449E+00/ 
      DATA ((A_EXTX(I,J,13),I=1,mbx),J=1,8)/ 0.9980E+00, 0.1032E+01, 
     1   0.1084E+01, 0.1141E+01, 0.1180E+01, 0.1239E+01, 0.1267E+01, 
     2   0.1223E+01, 0.1223E+01, 0.1131E+01, 0.1187E+01, 0.1193E+01, 
     3   0.1003E+01, 0.7764E+00, 0.8981E+00, 0.8734E+00, 0.7043E+00, 
     4   0.1087E+01, 0.9993E+00, 0.1023E+01, 0.1062E+01, 0.1112E+01, 
     5   0.1116E+01, 0.1186E+01, 0.1231E+01, 0.1199E+01, 0.1211E+01, 
     6   0.1190E+01, 0.1166E+01, 0.1088E+01, 0.8417E+00, 0.8725E+00, 
     7   0.1002E+01, 0.1018E+01, 0.9030E+00, 0.8425E+00, 0.9997E+00, 
     8   0.1023E+01, 0.1056E+01, 0.1100E+01, 0.1103E+01, 0.1170E+01, 
     9   0.1217E+01, 0.1196E+01, 0.1210E+01, 0.1203E+01, 0.1176E+01, 
     .   0.1094E+01, 0.8508E+00, 0.9063E+00, 0.1030E+01, 0.1053E+01, 
     1   0.9529E+00, 0.8484E+00, 0.1000E+01, 0.1022E+01, 0.1054E+01, 
     2   0.1094E+01, 0.1097E+01, 0.1157E+01, 0.1208E+01, 0.1194E+01, 
     3   0.1209E+01, 0.1213E+01, 0.1186E+01, 0.1105E+01, 0.8661E+00, 
     4   0.9329E+00, 0.1052E+01, 0.1081E+01, 0.9940E+00, 0.8664E+00, 
     5   0.1000E+01, 0.1020E+01, 0.1045E+01, 0.1079E+01, 0.1084E+01, 
     6   0.1135E+01, 0.1183E+01, 0.1182E+01, 0.1200E+01, 0.1219E+01, 
     7   0.1201E+01, 0.1130E+01, 0.9010E+00, 0.9733E+00, 0.1084E+01, 
     8   0.1122E+01, 0.1062E+01, 0.9151E+00, 0.9997E+00, 0.1015E+01, 
     9   0.1039E+01, 0.1066E+01, 0.1070E+01, 0.1111E+01, 0.1154E+01, 
     .   0.1163E+01, 0.1181E+01, 0.1209E+01, 0.1205E+01, 0.1153E+01, 
     1   0.9422E+00, 0.1006E+01, 0.1107E+01, 0.1153E+01, 0.1123E+01, 
     2   0.9803E+00, 0.9997E+00, 0.1013E+01, 0.1034E+01, 0.1052E+01, 
     3   0.1058E+01, 0.1087E+01, 0.1119E+01, 0.1132E+01, 0.1150E+01, 
     4   0.1181E+01, 0.1193E+01, 0.1172E+01, 0.9974E+00, 0.1041E+01, 
     5   0.1125E+01, 0.1178E+01, 0.1183E+01, 0.1076E+01, 0.9997E+00, 
     6   0.1011E+01, 0.1027E+01, 0.1043E+01, 0.1048E+01, 0.1071E+01, 
     7   0.1097E+01, 0.1109E+01, 0.1124E+01, 0.1152E+01, 0.1169E+01, 
     8   0.1167E+01, 0.1027E+01, 0.1055E+01, 0.1127E+01, 0.1180E+01, 
     9   0.1203E+01, 0.1130E+01/ 
      DATA ((A_SSAX(I,J,13),I=1,mbx),J=1,8)/ 0.1000E+01, 0.9930E+00, 
     1   0.9727E+00, 0.9556E+00, 0.8560E+00, 0.9741E+00, 0.9710E+00, 
     2   0.9527E+00, 0.9273E+00, 0.9159E+00, 0.8371E+00, 0.8800E+00, 
     3   0.9069E+00, 0.8668E+00, 0.6916E+00, 0.5906E+00, 0.4595E+00, 
     4   0.2559E+00, 0.1000E+01, 0.9975E+00, 0.9823E+00, 0.9394E+00, 
     5   0.7044E+00, 0.9079E+00, 0.8498E+00, 0.8049E+00, 0.7115E+00, 
     6   0.7686E+00, 0.7483E+00, 0.7304E+00, 0.5310E+00, 0.3977E+00, 
     7   0.4240E+00, 0.4467E+00, 0.4340E+00, 0.2669E+00, 0.1000E+01, 
     8   0.9979E+00, 0.9830E+00, 0.9334E+00, 0.6994E+00, 0.8936E+00, 
     9   0.8287E+00, 0.7859E+00, 0.6880E+00, 0.7474E+00, 0.7322E+00, 
     .   0.7098E+00, 0.5087E+00, 0.4006E+00, 0.4249E+00, 0.4450E+00, 
     1   0.4379E+00, 0.2892E+00, 0.1000E+01, 0.9982E+00, 0.9833E+00, 
     2   0.9271E+00, 0.6955E+00, 0.8822E+00, 0.8119E+00, 0.7713E+00, 
     3   0.6711E+00, 0.7314E+00, 0.7197E+00, 0.6954E+00, 0.4980E+00, 
     4   0.4064E+00, 0.4286E+00, 0.4461E+00, 0.4420E+00, 0.3084E+00, 
     5   0.1000E+01, 0.9986E+00, 0.9832E+00, 0.9143E+00, 0.6881E+00, 
     6   0.8606E+00, 0.7812E+00, 0.7451E+00, 0.6433E+00, 0.7029E+00, 
     7   0.6968E+00, 0.6726E+00, 0.4893E+00, 0.4204E+00, 0.4386E+00, 
     8   0.4515E+00, 0.4504E+00, 0.3412E+00, 0.1000E+01, 0.9991E+00, 
     9   0.9817E+00, 0.8973E+00, 0.6799E+00, 0.8361E+00, 0.7480E+00, 
     .   0.7169E+00, 0.6158E+00, 0.6723E+00, 0.6710E+00, 0.6503E+00, 
     1   0.4885E+00, 0.4367E+00, 0.4511E+00, 0.4595E+00, 0.4593E+00, 
     2   0.3725E+00, 0.1000E+01, 0.9992E+00, 0.9791E+00, 0.8695E+00, 
     3   0.6676E+00, 0.8005E+00, 0.7024E+00, 0.6776E+00, 0.5821E+00, 
     4   0.6297E+00, 0.6327E+00, 0.6192E+00, 0.4926E+00, 0.4581E+00, 
     5   0.4687E+00, 0.4714E+00, 0.4702E+00, 0.4072E+00, 0.1000E+01, 
     6   0.9992E+00, 0.9761E+00, 0.8453E+00, 0.6583E+00, 0.7726E+00, 
     7   0.6703E+00, 0.6494E+00, 0.5614E+00, 0.6001E+00, 0.6039E+00, 
     8   0.5952E+00, 0.4959E+00, 0.4723E+00, 0.4807E+00, 0.4801E+00, 
     9   0.4769E+00, 0.4282E+00/ 
      DATA ((A_ASYX(I,J,13),I=1,mbx),J=1,8)/ 0.7964E+00, 0.7818E+00, 
     1   0.7631E+00, 0.7611E+00, 0.7325E+00, 0.7164E+00, 0.7131E+00, 
     2   0.7702E+00, 0.7411E+00, 0.7655E+00, 0.7186E+00, 0.6815E+00, 
     3   0.7176E+00, 0.7450E+00, 0.6774E+00, 0.6201E+00, 0.5680E+00, 
     4   0.2940E+00, 0.8469E+00, 0.8377E+00, 0.8242E+00, 0.8204E+00, 
     5   0.8874E+00, 0.7840E+00, 0.8116E+00, 0.8563E+00, 0.8519E+00, 
     6   0.8503E+00, 0.8474E+00, 0.8563E+00, 0.8891E+00, 0.8556E+00, 
     7   0.7968E+00, 0.7386E+00, 0.6789E+00, 0.5294E+00, 0.8506E+00, 
     8   0.8444E+00, 0.8319E+00, 0.8316E+00, 0.8958E+00, 0.7934E+00, 
     9   0.8231E+00, 0.8648E+00, 0.8639E+00, 0.8604E+00, 0.8618E+00, 
     .   0.8742E+00, 0.9029E+00, 0.8659E+00, 0.8120E+00, 0.7583E+00, 
     1   0.7014E+00, 0.5713E+00, 0.8570E+00, 0.8500E+00, 0.8384E+00, 
     2   0.8389E+00, 0.9015E+00, 0.8023E+00, 0.8313E+00, 0.8704E+00, 
     3   0.8722E+00, 0.8666E+00, 0.8712E+00, 0.8855E+00, 0.9122E+00, 
     4   0.8744E+00, 0.8241E+00, 0.7736E+00, 0.7190E+00, 0.6019E+00, 
     5   0.8604E+00, 0.8566E+00, 0.8500E+00, 0.8528E+00, 0.9109E+00, 
     6   0.8159E+00, 0.8449E+00, 0.8796E+00, 0.8850E+00, 0.8777E+00, 
     7   0.8852E+00, 0.9014E+00, 0.9259E+00, 0.8882E+00, 0.8437E+00, 
     8   0.7991E+00, 0.7488E+00, 0.6523E+00, 0.8629E+00, 0.8628E+00, 
     9   0.8581E+00, 0.8642E+00, 0.9207E+00, 0.8332E+00, 0.8599E+00, 
     .   0.8889E+00, 0.8985E+00, 0.8885E+00, 0.8968E+00, 0.9141E+00, 
     1   0.9375E+00, 0.9015E+00, 0.8625E+00, 0.8239E+00, 0.7785E+00, 
     2   0.6962E+00, 0.8680E+00, 0.8658E+00, 0.8674E+00, 0.8815E+00, 
     3   0.9323E+00, 0.8570E+00, 0.8826E+00, 0.9038E+00, 0.9169E+00, 
     4   0.9042E+00, 0.9111E+00, 0.9276E+00, 0.9504E+00, 0.9180E+00, 
     5   0.8843E+00, 0.8539E+00, 0.8157E+00, 0.7501E+00, 0.8671E+00, 
     6   0.8679E+00, 0.8732E+00, 0.8919E+00, 0.9398E+00, 0.8746E+00, 
     7   0.9008E+00, 0.9166E+00, 0.9314E+00, 0.9176E+00, 0.9221E+00, 
     8   0.9365E+00, 0.9578E+00, 0.9277E+00, 0.8974E+00, 0.8722E+00, 
     9   0.8395E+00, 0.7834E+00/ 
      DATA ((A_EXTX(I,J,14),I=1,mbx),J=1,1)/ 0.6970E+00, 0.3724E+00, 
     1   0.1420E+00, 0.6483E-01, 0.3820E-01, 0.1519E-01, 0.8261E-02, 
     2   0.1004E-01, 0.1296E-01, 0.1639E-01, 0.3130E-01, 0.2811E-01, 
     3   0.2521E-01, 0.1703E-01, 0.1988E-01, 0.1697E-01, 0.9376E-02, 
     4   0.9107E-02/ 
      DATA ((A_SSAX(I,J,14),I=1,mbx),J=1,1)/ 0.9647E+00, 0.9747E+00, 
     1   0.9551E+00, 0.9100E+00, 0.6865E+00, 0.7466E+00, 0.5914E+00, 
     2   0.2448E+00, 0.1105E+00, 0.2806E-01, 0.1280E-01, 0.9266E-01, 
     3   0.2898E-01, 0.1546E-01, 0.4009E-02, 0.7202E-02, 0.3913E-02, 
     4   0.4404E-03/ 
      DATA ((A_ASYX(I,J,14),I=1,mbx),J=1,1)/ 0.6649E+00, 0.6163E+00, 
     1   0.5404E+00, 0.4736E+00, 0.4018E+00, 0.3402E+00, 0.2763E+00, 
     2   0.2324E+00, 0.1920E+00, 0.1511E+00, 0.1181E+00, 0.1573E+00, 
     3   0.1089E+00, 0.8380E-01, 0.5339E-01, 0.5149E-01, 0.3630E-01, 
     4   0.1260E-01/ 
      DATA ((A_EXTX(I,J,15),I=1,mbx),J=1,1)/ 0.9984E+00, 0.1086E+01, 
     1   0.1096E+01, 0.9933E+00, 0.8202E+00, 0.6341E+00, 0.4556E+00, 
     2   0.3471E+00, 0.2878E+00, 0.1996E+00, 0.2565E+00, 0.6046E+00, 
     3   0.3391E+00, 0.2277E+00, 0.1790E+00, 0.2234E+00, 0.1218E+00, 
     4   0.6990E-01/ 
      DATA ((A_SSAX(I,J,15),I=1,mbx),J=1,1)/ 0.8711E+00, 0.9378E+00, 
     1   0.9463E+00, 0.9390E+00, 0.8556E+00, 0.9280E+00, 0.9132E+00, 
     2   0.7796E+00, 0.6446E+00, 0.3766E+00, 0.1883E+00, 0.4505E+00, 
     3   0.3751E+00, 0.3398E+00, 0.1766E+00, 0.2689E+00, 0.2345E+00, 
     4   0.6417E-01/ 
      DATA ((A_ASYX(I,J,15),I=1,mbx),J=1,1)/ 0.7372E+00, 0.6959E+00, 
     1   0.6875E+00, 0.6870E+00, 0.6976E+00, 0.6754E+00, 0.6587E+00, 
     2   0.6577E+00, 0.6356E+00, 0.6194E+00, 0.5500E+00, 0.3734E+00, 
     3   0.4415E+00, 0.4217E+00, 0.3678E+00, 0.2494E+00, 0.2310E+00, 
     4   0.1612E+00/ 
      DATA ((A_EXTX(I,J,16),I=1,mbx),J=1,1)/ 0.9996E+00, 0.1027E+01, 
     1   0.1068E+01, 0.1107E+01, 0.1148E+01, 0.1198E+01, 0.1233E+01, 
     2   0.1224E+01, 0.1191E+01, 0.1019E+01, 0.8557E+00, 0.1258E+01, 
     3   0.1215E+01, 0.1151E+01, 0.9892E+00, 0.1223E+01, 0.1120E+01, 
     4   0.8345E+00/ 
      DATA ((A_SSAX(I,J,16),I=1,mbx),J=1,1)/ 0.6601E+00, 0.7660E+00, 
     1   0.7855E+00, 0.7760E+00, 0.6695E+00, 0.8055E+00, 0.8213E+00, 
     2   0.7032E+00, 0.6400E+00, 0.5581E+00, 0.4115E+00, 0.5215E+00, 
     3   0.4952E+00, 0.5059E+00, 0.4368E+00, 0.4754E+00, 0.4695E+00, 
     4   0.3922E+00/ 
      DATA ((A_ASYX(I,J,16),I=1,mbx),J=1,1)/ 0.8973E+00, 0.8441E+00, 
     1   0.8113E+00, 0.7920E+00, 0.8221E+00, 0.7620E+00, 0.7560E+00, 
     2   0.8061E+00, 0.8300E+00, 0.8774E+00, 0.8754E+00, 0.6871E+00, 
     3   0.7447E+00, 0.7345E+00, 0.7462E+00, 0.5502E+00, 0.4931E+00, 
     4   0.4435E+00/ 
      DATA ((A_EXTX(I,J,17),I=1,mbx),J=1,1)/ 0.9986E+00, 0.1075E+01, 
     1   0.1146E+01, 0.1147E+01, 0.1071E+01, 0.9633E+00, 0.8081E+00, 
     2   0.6535E+00, 0.5407E+00, 0.3475E+00, 0.3921E+00, 0.9254E+00, 
     3   0.5998E+00, 0.4218E+00, 0.3044E+00, 0.4378E+00, 0.2417E+00, 
     4   0.1079E+00/ 
      DATA ((A_SSAX(I,J,17),I=1,mbx),J=1,1)/ 0.8289E+00, 0.9121E+00, 
     1   0.9248E+00, 0.9199E+00, 0.8342E+00, 0.9235E+00, 0.9212E+00, 
     2   0.8118E+00, 0.7005E+00, 0.4500E+00, 0.2373E+00, 0.4742E+00, 
     3   0.4401E+00, 0.4147E+00, 0.2391E+00, 0.3436E+00, 0.2894E+00, 
     4   0.6579E-01/ 
      DATA ((A_ASYX(I,J,17),I=1,mbx),J=1,1)/ 0.7784E+00, 0.7216E+00, 
     1   0.6970E+00, 0.6933E+00, 0.7183E+00, 0.6974E+00, 0.7035E+00, 
     2   0.7204E+00, 0.7080E+00, 0.6903E+00, 0.6241E+00, 0.4321E+00, 
     3   0.5088E+00, 0.4685E+00, 0.3810E+00, 0.2530E+00, 0.1987E+00, 
     4   0.7565E-01/ 
      DATA ((A_EXTX(I,J,18),I=1,mbx),J=1,8)/ 0.1009E+01, 0.5315E+00, 
     1   0.1916E+00, 0.7670E-01, 0.7874E-01, 0.9590E-01, 0.6253E-01, 
     2   0.7484E-01, 0.4925E-01, 0.9857E-01, 0.1587E+00, 0.9879E-01, 
     3   0.4700E-01, 0.2393E-01, 0.3030E-01, 0.5892E-02, 0.9433E-02, 
     4   0.7848E-02, 0.1006E+01, 0.5952E+00, 0.2531E+00, 0.1037E+00, 
     5   0.1446E+00, 0.7541E-01, 0.4527E-01, 0.5074E-01, 0.4011E-01, 
     6   0.5088E-01, 0.8283E-01, 0.6455E-01, 0.4619E-01, 0.5258E-01, 
     7   0.5271E-01, 0.2665E-01, 0.1876E-01, 0.1624E-01, 0.1006E+01, 
     8   0.6212E+00, 0.2790E+00, 0.1177E+00, 0.1633E+00, 0.7724E-01, 
     9   0.4479E-01, 0.4768E-01, 0.4012E-01, 0.4406E-01, 0.6866E-01, 
     .   0.5540E-01, 0.4558E-01, 0.5987E-01, 0.5860E-01, 0.3197E-01, 
     1   0.2137E-01, 0.1845E-01, 0.1005E+01, 0.6426E+00, 0.3008E+00, 
     2   0.1302E+00, 0.1777E+00, 0.8106E-01, 0.4587E-01, 0.4669E-01, 
     3   0.4100E-01, 0.4064E-01, 0.6060E-01, 0.4993E-01, 0.4548E-01, 
     4   0.6509E-01, 0.6299E-01, 0.3575E-01, 0.2328E-01, 0.2004E-01, 
     5   0.1004E+01, 0.6814E+00, 0.3426E+00, 0.1554E+00, 0.2042E+00, 
     6   0.9223E-01, 0.5051E-01, 0.4744E-01, 0.4419E-01, 0.3753E-01, 
     7   0.5094E-01, 0.4315E-01, 0.4619E-01, 0.7387E-01, 0.7072E-01, 
     8   0.4208E-01, 0.2658E-01, 0.2275E-01, 0.1002E+01, 0.7316E+00, 
     9   0.4003E+00, 0.1929E+00, 0.2399E+00, 0.1134E+00, 0.6077E-01, 
     .   0.5232E-01, 0.5093E-01, 0.3767E-01, 0.4499E-01, 0.3886E-01, 
     1   0.4868E-01, 0.8496E-01, 0.8106E-01, 0.5002E-01, 0.3089E-01, 
     2   0.2623E-01, 0.1001E+01, 0.8035E+00, 0.4928E+00, 0.2593E+00, 
     3   0.2985E+00, 0.1572E+00, 0.8413E-01, 0.6619E-01, 0.6603E-01, 
     4   0.4387E-01, 0.4423E-01, 0.3828E-01, 0.5504E-01, 0.1028E+00, 
     5   0.9851E-01, 0.6287E-01, 0.3814E-01, 0.3191E-01, 0.1000E+01, 
     6   0.8626E+00, 0.5794E+00, 0.3284E+00, 0.3569E+00, 0.2076E+00, 
     7   0.1130E+00, 0.8469E-01, 0.8470E-01, 0.5441E-01, 0.4914E-01, 
     8   0.4181E-01, 0.6305E-01, 0.1209E+00, 0.1170E+00, 0.7631E-01, 
     9   0.4595E-01, 0.3786E-01/ 
      DATA ((A_SSAX(I,J,18),I=1,mbx),J=1,8)/ 0.1000E+01, 0.1000E+01, 
     1   0.9976E+00, 0.9708E+00, 0.4906E+00, 0.1774E+00, 0.1232E+00, 
     2   0.6446E-01, 0.5671E-01, 0.1527E-01, 0.2826E-01, 0.4081E-01, 
     3   0.4101E-01, 0.2760E-01, 0.3588E-01, 0.7646E-01, 0.4123E-02, 
     4   0.9668E-03, 0.1000E+01, 0.1000E+01, 0.9983E+00, 0.9788E+00, 
     5   0.5812E+00, 0.4741E+00, 0.3407E+00, 0.1657E+00, 0.1505E+00, 
     6   0.5469E-01, 0.4088E-01, 0.4982E-01, 0.3349E-01, 0.1680E-01, 
     7   0.1471E-01, 0.1438E-01, 0.6970E-02, 0.1126E-02, 0.1000E+01, 
     8   0.1000E+01, 0.9985E+00, 0.9808E+00, 0.6104E+00, 0.5805E+00, 
     9   0.4312E+00, 0.2226E+00, 0.1945E+00, 0.8345E-01, 0.5229E-01, 
     .   0.5785E-01, 0.3432E-01, 0.1778E-01, 0.1680E-01, 0.1577E-01, 
     1   0.8189E-02, 0.1382E-02, 0.1000E+01, 0.1000E+01, 0.9986E+00, 
     2   0.9821E+00, 0.6299E+00, 0.6530E+00, 0.4985E+00, 0.2720E+00, 
     3   0.2305E+00, 0.1115E+00, 0.6457E-01, 0.6601E-01, 0.3573E-01, 
     4   0.1918E-01, 0.1896E-01, 0.1750E-01, 0.9418E-02, 0.1639E-02, 
     5   0.1000E+01, 0.1000E+01, 0.9987E+00, 0.9839E+00, 0.6588E+00, 
     6   0.7549E+00, 0.6029E+00, 0.3635E+00, 0.2943E+00, 0.1715E+00, 
     7   0.9478E-01, 0.8494E-01, 0.3964E-01, 0.2273E-01, 0.2375E-01, 
     8   0.2170E-01, 0.1224E-01, 0.2239E-02, 0.1000E+01, 0.1000E+01, 
     9   0.9988E+00, 0.9856E+00, 0.6855E+00, 0.8382E+00, 0.7003E+00, 
     .   0.4687E+00, 0.3675E+00, 0.2569E+00, 0.1474E+00, 0.1163E+00, 
     1   0.4669E-01, 0.2869E-01, 0.3127E-01, 0.2874E-01, 0.1698E-01, 
     2   0.3285E-02, 0.1000E+01, 0.1000E+01, 0.9988E+00, 0.9871E+00, 
     3   0.7126E+00, 0.9027E+00, 0.7880E+00, 0.5850E+00, 0.4549E+00, 
     4   0.3788E+00, 0.2437E+00, 0.1737E+00, 0.6041E-01, 0.3976E-01, 
     5   0.4490E-01, 0.4222E-01, 0.2651E-01, 0.5546E-02, 0.1000E+01, 
     6   0.1000E+01, 0.9988E+00, 0.9880E+00, 0.7286E+00, 0.9306E+00, 
     7   0.8321E+00, 0.6533E+00, 0.5139E+00, 0.4686E+00, 0.3312E+00, 
     8   0.2294E+00, 0.7508E-01, 0.5140E-01, 0.5907E-01, 0.5688E-01, 
     9   0.3758E-01, 0.8392E-02/ 
      DATA ((A_ASYX(I,J,18),I=1,mbx),J=1,8)/ 0.7172E+00, 0.6760E+00, 
     1   0.6086E+00, 0.5473E+00, 0.4571E+00, 0.3765E+00, 0.3163E+00, 
     2   0.2661E+00, 0.2370E+00, 0.1704E+00, 0.1353E+00, 0.1478E+00, 
     3   0.1519E+00, 0.1258E+00, 0.9522E-01, 0.8572E-01, 0.5257E-01, 
     4   0.2308E-01, 0.7690E+00, 0.7404E+00, 0.6846E+00, 0.6391E+00, 
     5   0.5324E+00, 0.4853E+00, 0.4211E+00, 0.3673E+00, 0.3292E+00, 
     6   0.2748E+00, 0.2325E+00, 0.2197E+00, 0.1904E+00, 0.1480E+00, 
     7   0.1204E+00, 0.1042E+00, 0.7121E-01, 0.3384E-01, 0.7779E+00, 
     8   0.7541E+00, 0.7040E+00, 0.6637E+00, 0.5565E+00, 0.5153E+00, 
     9   0.4523E+00, 0.3990E+00, 0.3583E+00, 0.3077E+00, 0.2645E+00, 
     .   0.2459E+00, 0.2086E+00, 0.1601E+00, 0.1321E+00, 0.1139E+00, 
     1   0.7974E-01, 0.3877E-01, 0.7837E+00, 0.7632E+00, 0.7173E+00, 
     2   0.6810E+00, 0.5745E+00, 0.5355E+00, 0.4752E+00, 0.4226E+00, 
     3   0.3808E+00, 0.3328E+00, 0.2889E+00, 0.2660E+00, 0.2231E+00, 
     4   0.1705E+00, 0.1414E+00, 0.1221E+00, 0.8694E-01, 0.4311E-01, 
     5   0.7900E+00, 0.7755E+00, 0.7368E+00, 0.7068E+00, 0.6048E+00, 
     6   0.5693E+00, 0.5128E+00, 0.4620E+00, 0.4186E+00, 0.3739E+00, 
     7   0.3307E+00, 0.3022E+00, 0.2507E+00, 0.1908E+00, 0.1600E+00, 
     8   0.1377E+00, 0.1008E+00, 0.5131E-01, 0.7948E+00, 0.7864E+00, 
     9   0.7564E+00, 0.7336E+00, 0.6385E+00, 0.6040E+00, 0.5541E+00, 
     .   0.5070E+00, 0.4625E+00, 0.4216E+00, 0.3796E+00, 0.3465E+00, 
     1   0.2864E+00, 0.2176E+00, 0.1836E+00, 0.1587E+00, 0.1193E+00, 
     2   0.6287E-01, 0.7975E+00, 0.7959E+00, 0.7766E+00, 0.7635E+00, 
     3   0.6804E+00, 0.6449E+00, 0.6051E+00, 0.5647E+00, 0.5193E+00, 
     4   0.4834E+00, 0.4441E+00, 0.4081E+00, 0.3392E+00, 0.2584E+00, 
     5   0.2196E+00, 0.1909E+00, 0.1477E+00, 0.8069E-01, 0.7967E+00, 
     6   0.8007E+00, 0.7892E+00, 0.7825E+00, 0.7117E+00, 0.6735E+00, 
     7   0.6423E+00, 0.6084E+00, 0.5632E+00, 0.5303E+00, 0.4951E+00, 
     8   0.4587E+00, 0.3850E+00, 0.2951E+00, 0.2518E+00, 0.2202E+00, 
     9   0.1738E+00, 0.9837E-01/ 
      END 
 
      BLOCK DATA AEROSOL2
      INCLUDE 'fl_radparams.inc'
c##      include 'rad_0698.h'
      INTEGER I, J
      real a_ssay(mby,nrh,naer),a_exty(mby,nrh,naer)
      real a_asyy(mby,nrh,naer)
      COMMON /AER_OPTY/A_SSAY, A_EXTY, A_ASYY
      DATA ((A_SSAY(I,J,1),I=1,mby),j=1,nrh)/ .1000E+01, .1000E+01, 
     1   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     2   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     3   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     4   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     5   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     6   .1000E+01, .1000E+01, .1000E+01, .9999E+00, .9999E+00, 
     7   .9999E+00, .9999E+00, .9999E+00, .1000E+01, .1000E+01, 
     8   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     9   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     .   .1000E+01, .1000E+01, .1000E+01, .9993E+00, .9993E+00, 
     1   .9993E+00, .9993E+00, .9993E+00, .1000E+01, .1000E+01, 
     2   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     3   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     4   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     5   .1000E+01, .1000E+01, .1000E+01, .1000E+01, .1000E+01, 
     6   .1000E+01, .1000E+01, .1000E+01/ 
      DATA ((A_EXTY(I,J,1),I=1,mby),j=1,nrh)/ .2071E-03, .2071E-03, 
     1   .2071E-03, .2071E-03, .2071E-03, .2084E-03, .2081E-03, 
     2   .2065E-03, .2071E-03, .2101E-03, .2448E-03, .2448E-03, 
     3   .2448E-03, .2448E-03, .2448E-03, .2459E-03, .2452E-03, 
     4   .2437E-03, .2427E-03, .2447E-03, .3519E-03, .3519E-03, 
     5   .3519E-03, .3519E-03, .3519E-03, .3499E-03, .3503E-03, 
     6   .3510E-03, .3486E-03, .3468E-03, .7975E-03, .7975E-03, 
     7   .7975E-03, .7975E-03, .7975E-03, .7928E-03, .7874E-03, 
     8   .7863E-03, .7813E-03, .7843E-03, .1135E-02, .1135E-02, 
     9   .1135E-02, .1135E-02, .1135E-02, .1122E-02, .1120E-02, 
     .   .1113E-02, .1113E-02, .1106E-02, .1685E-02, .1685E-02, 
     1   .1685E-02, .1685E-02, .1685E-02, .1671E-02, .1656E-02, 
     2   .1644E-02, .1632E-02, .1626E-02, .2879E-02, .2879E-02, 
     3   .2879E-02, .2879E-02, .2879E-02, .2872E-02, .2855E-02, 
     4   .2832E-02, .2806E-02, .2770E-02, .4241E-02, .4241E-02, 
     5   .4241E-02, .4241E-02, .4241E-02, .4274E-02, .4256E-02, 
     6   .4255E-02, .4223E-02, .4171E-02/ 
      DATA ((A_ASYY(I,J,1),I=1,mby),j=1,nrh)/ .7513E+00, .7513E+00, 
     1   .7513E+00, .7513E+00, .7513E+00, .7721E+00, .7842E+00, 
     2   .7893E+00, .7963E+00, .8072E+00, .7568E+00, .7568E+00, 
     3   .7568E+00, .7568E+00, .7568E+00, .7792E+00, .7907E+00, 
     4   .7940E+00, .8002E+00, .8113E+00, .7412E+00, .7412E+00, 
     5   .7412E+00, .7412E+00, .7412E+00, .7662E+00, .7783E+00, 
     6   .7912E+00, .8007E+00, .8096E+00, .6857E+00, .6857E+00, 
     7   .6857E+00, .6857E+00, .6857E+00, .7078E+00, .7249E+00, 
     8   .7462E+00, .7600E+00, .7868E+00, .6639E+00, .6639E+00, 
     9   .6639E+00, .6639E+00, .6639E+00, .6845E+00, .7070E+00, 
     .   .7252E+00, .7393E+00, .7655E+00, .6515E+00, .6515E+00, 
     1   .6515E+00, .6515E+00, .6515E+00, .6620E+00, .6810E+00, 
     2   .6925E+00, .7165E+00, .7380E+00, .6220E+00, .6220E+00, 
     3   .6220E+00, .6220E+00, .6220E+00, .6424E+00, .6525E+00, 
     4   .6656E+00, .6848E+00, .7081E+00, .6129E+00, .6129E+00, 
     5   .6129E+00, .6129E+00, .6129E+00, .6290E+00, .6397E+00, 
     6   .6509E+00, .6676E+00, .6865E+00/ 
      DATA ((A_SSAY(I,J,2),I=1,mby),j=1,nrh)/ .9419E+00, .9419E+00, 
     1   .9419E+00, .9419E+00, .9419E+00, .9634E+00, .9640E+00, 
     2   .9652E+00, .9628E+00, .9566E+00, .9418E+00, .9418E+00, 
     3   .9418E+00, .9418E+00, .9418E+00, .9633E+00, .9640E+00, 
     4   .9652E+00, .9627E+00, .9565E+00, .9460E+00, .9460E+00, 
     5   .9460E+00, .9460E+00, .9460E+00, .9650E+00, .9667E+00, 
     6   .9673E+00, .9650E+00, .9595E+00, .9596E+00, .9596E+00, 
     7   .9596E+00, .9596E+00, .9596E+00, .9744E+00, .9754E+00, 
     8   .9760E+00, .9742E+00, .9703E+00, .9722E+00, .9722E+00, 
     9   .9722E+00, .9722E+00, .9722E+00, .9827E+00, .9833E+00, 
     .   .9838E+00, .9828E+00, .9805E+00, .9776E+00, .9776E+00, 
     1   .9776E+00, .9776E+00, .9776E+00, .9861E+00, .9869E+00, 
     2   .9872E+00, .9865E+00, .9846E+00, .9823E+00, .9823E+00, 
     3   .9823E+00, .9823E+00, .9823E+00, .9892E+00, .9895E+00, 
     4   .9900E+00, .9896E+00, .9883E+00, .9857E+00, .9857E+00, 
     5   .9857E+00, .9857E+00, .9857E+00, .9912E+00, .9917E+00, 
     6   .9921E+00, .9919E+00, .9907E+00/ 
      DATA ((A_EXTY(I,J,2),I=1,mby),j=1,nrh)/ .1763E-04, .1763E-04, 
     1   .1763E-04, .1763E-04, .1763E-04, .1574E-04, .1402E-04, 
     2   .1248E-04, .1055E-04, .8482E-05, .1763E-04, .1763E-04, 
     3   .1763E-04, .1763E-04, .1763E-04, .1574E-04, .1402E-04, 
     4   .1249E-04, .1055E-04, .8483E-05, .1890E-04, .1890E-04, 
     5   .1890E-04, .1890E-04, .1890E-04, .1689E-04, .1504E-04, 
     6   .1339E-04, .1132E-04, .9110E-05, .2535E-04, .2535E-04, 
     7   .2535E-04, .2535E-04, .2535E-04, .2270E-04, .2027E-04, 
     8   .1811E-04, .1538E-04, .1244E-04, .3707E-04, .3707E-04, 
     9   .3707E-04, .3707E-04, .3707E-04, .3347E-04, .3014E-04, 
     .   .2714E-04, .2326E-04, .1903E-04, .4636E-04, .4636E-04, 
     1   .4636E-04, .4636E-04, .4636E-04, .4215E-04, .3817E-04, 
     2   .3459E-04, .2986E-04, .2465E-04, .5890E-04, .5890E-04, 
     3   .5890E-04, .5890E-04, .5890E-04, .5402E-04, .4933E-04, 
     4   .4501E-04, .3919E-04, .3269E-04, .7312E-04, .7312E-04, 
     5   .7312E-04, .7312E-04, .7312E-04, .6769E-04, .6224E-04, 
     6   .5721E-04, .5027E-04, .4240E-04/ 
      DATA ((A_ASYY(I,J,2),I=1,mby),j=1,nrh)/ .6740E+00, .6740E+00, 
     1   .6740E+00, .6740E+00, .6740E+00, .6635E+00, .6570E+00, 
     2   .6507E+00, .6414E+00, .6293E+00, .6740E+00, .6740E+00, 
     3   .6740E+00, .6740E+00, .6740E+00, .6635E+00, .6570E+00, 
     4   .6507E+00, .6414E+00, .6293E+00, .6809E+00, .6809E+00, 
     5   .6809E+00, .6809E+00, .6809E+00, .6740E+00, .6678E+00, 
     6   .6616E+00, .6523E+00, .6403E+00, .7167E+00, .7167E+00, 
     7   .7167E+00, .7167E+00, .7167E+00, .7097E+00, .7046E+00, 
     8   .6988E+00, .6904E+00, .6785E+00, .7447E+00, .7447E+00, 
     9   .7447E+00, .7447E+00, .7447E+00, .7407E+00, .7371E+00, 
     .   .7325E+00, .7251E+00, .7146E+00, .7561E+00, .7561E+00, 
     1   .7561E+00, .7561E+00, .7561E+00, .7534E+00, .7508E+00, 
     2   .7468E+00, .7404E+00, .7308E+00, .7656E+00, .7656E+00, 
     3   .7656E+00, .7656E+00, .7656E+00, .7643E+00, .7622E+00, 
     4   .7589E+00, .7536E+00, .7451E+00, .7723E+00, .7723E+00, 
     5   .7723E+00, .7723E+00, .7723E+00, .7715E+00, .7706E+00, 
     6   .7678E+00, .7635E+00, .7559E+00/ 
      DATA ((A_SSAY(I,J,3),I=1,mby),j=1,nrh)/ .9180E+00, .9180E+00, 
     1   .9180E+00, .9180E+00, .9180E+00, .9394E+00, .9404E+00, 
     2   .9417E+00, .9391E+00, .9333E+00, .9174E+00, .9174E+00, 
     3   .9174E+00, .9174E+00, .9174E+00, .9388E+00, .9397E+00, 
     4   .9411E+00, .9384E+00, .9327E+00, .9210E+00, .9210E+00, 
     5   .9210E+00, .9210E+00, .9210E+00, .9400E+00, .9421E+00, 
     6   .9428E+00, .9403E+00, .9353E+00, .9377E+00, .9377E+00, 
     7   .9377E+00, .9377E+00, .9377E+00, .9527E+00, .9543E+00, 
     8   .9551E+00, .9533E+00, .9500E+00, .9553E+00, .9553E+00, 
     9   .9553E+00, .9553E+00, .9553E+00, .9663E+00, .9675E+00, 
     .   .9685E+00, .9676E+00, .9659E+00, .9630E+00, .9630E+00, 
     1   .9630E+00, .9630E+00, .9630E+00, .9722E+00, .9736E+00, 
     2   .9743E+00, .9739E+00, .9728E+00, .9702E+00, .9702E+00, 
     3   .9702E+00, .9702E+00, .9702E+00, .9776E+00, .9786E+00, 
     4   .9795E+00, .9795E+00, .9788E+00, .9756E+00, .9756E+00, 
     5   .9756E+00, .9756E+00, .9756E+00, .9816E+00, .9827E+00, 
     6   .9836E+00, .9837E+00, .9832E+00/ 
      DATA ((A_EXTY(I,J,3),I=1,mby),j=1,nrh)/ .1160E-04, .1160E-04, 
     1   .1160E-04, .1160E-04, .1160E-04, .1033E-04, .9185E-05, 
     2   .8166E-05, .6890E-05, .5530E-05, .1161E-04, .1161E-04, 
     3   .1161E-04, .1161E-04, .1161E-04, .1034E-04, .9196E-05, 
     4   .8175E-05, .6897E-05, .5536E-05, .1248E-04, .1248E-04, 
     5   .1248E-04, .1248E-04, .1248E-04, .1112E-04, .9879E-05, 
     6   .8785E-05, .7413E-05, .5952E-05, .1675E-04, .1675E-04, 
     7   .1675E-04, .1675E-04, .1675E-04, .1494E-04, .1331E-04, 
     8   .1187E-04, .1005E-04, .8106E-05, .2446E-04, .2446E-04, 
     9   .2446E-04, .2446E-04, .2446E-04, .2199E-04, .1972E-04, 
     .   .1772E-04, .1514E-04, .1235E-04, .3079E-04, .3079E-04, 
     1   .3079E-04, .3079E-04, .3079E-04, .2783E-04, .2509E-04, 
     2   .2265E-04, .1948E-04, .1601E-04, .3977E-04, .3977E-04, 
     3   .3977E-04, .3977E-04, .3977E-04, .3616E-04, .3281E-04, 
     4   .2978E-04, .2579E-04, .2139E-04, .4994E-04, .4994E-04, 
     5   .4994E-04, .4994E-04, .4994E-04, .4577E-04, .4176E-04, 
     6   .3815E-04, .3329E-04, .2788E-04/ 
      DATA ((A_ASYY(I,J,3),I=1,mby),j=1,nrh)/ .6710E+00, .6710E+00, 
     1   .6710E+00, .6710E+00, .6710E+00, .6606E+00, .6543E+00, 
     2   .6481E+00, .6390E+00, .6271E+00, .6711E+00, .6711E+00, 
     3   .6711E+00, .6711E+00, .6711E+00, .6607E+00, .6543E+00, 
     4   .6481E+00, .6389E+00, .6270E+00, .6811E+00, .6811E+00, 
     5   .6811E+00, .6811E+00, .6811E+00, .6713E+00, .6652E+00, 
     6   .6590E+00, .6498E+00, .6379E+00, .7143E+00, .7143E+00, 
     7   .7143E+00, .7143E+00, .7143E+00, .7072E+00, .7020E+00, 
     8   .6962E+00, .6878E+00, .6760E+00, .7425E+00, .7425E+00, 
     9   .7425E+00, .7425E+00, .7425E+00, .7383E+00, .7346E+00, 
     .   .7299E+00, .7225E+00, .7121E+00, .7541E+00, .7541E+00, 
     1   .7541E+00, .7541E+00, .7541E+00, .7510E+00, .7482E+00, 
     2   .7440E+00, .7375E+00, .7279E+00, .7637E+00, .7637E+00, 
     3   .7637E+00, .7637E+00, .7637E+00, .7618E+00, .7593E+00, 
     4   .7557E+00, .7501E+00, .7414E+00, .7707E+00, .7707E+00, 
     5   .7707E+00, .7707E+00, .7707E+00, .7691E+00, .7677E+00, 
     6   .7645E+00, .7598E+00, .7519E+00/ 
      DATA ((A_SSAY(I,J,4),I=1,mby),j=1,nrh)/ .7035E+00, .7035E+00, 
     1   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     2   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     3   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     4   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     5   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     6   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     7   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     8   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     9   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     .   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     1   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     2   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     3   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     4   .8779E+00, .9276E+00, .9653E+00, .7035E+00, .7035E+00, 
     5   .7035E+00, .7035E+00, .7035E+00, .7798E+00, .8284E+00, 
     6   .8779E+00, .9276E+00, .9653E+00/ 
      DATA ((A_EXTY(I,J,4),I=1,mby),j=1,nrh)/ .8783E+00, .8783E+00, 
     1   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     2   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     3   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     4   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     5   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     6   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     7   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     8   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     9   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     .   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     1   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     2   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     3   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     4   .9674E+00, .1015E+01, .1067E+01, .8783E+00, .8783E+00, 
     5   .8783E+00, .8783E+00, .8783E+00, .9056E+00, .9356E+00, 
     6   .9674E+00, .1015E+01, .1067E+01/ 
      DATA ((A_ASYY(I,J,4),I=1,mby),j=1,nrh)/ .7678E+00, .7678E+00, 
     1   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     2   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     3   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     4   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     5   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     6   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     7   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     8   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     9   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     .   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     1   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     2   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     3   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     4   .6754E+00, .6626E+00, .6622E+00, .7678E+00, .7678E+00, 
     5   .7678E+00, .7678E+00, .7678E+00, .7230E+00, .6963E+00, 
     6   .6754E+00, .6626E+00, .6622E+00/ 
      DATA ((A_SSAY(I,J,5),I=1,mby),j=1,nrh)/ .6142E+00, .6142E+00, 
     1   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     2   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     3   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     4   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     5   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     6   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     7   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     8   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     9   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     .   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     1   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     2   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     3   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     4   .7920E+00, .8629E+00, .9255E+00, .6142E+00, .6142E+00, 
     5   .6142E+00, .6142E+00, .6142E+00, .6812E+00, .7317E+00, 
     6   .7920E+00, .8629E+00, .9255E+00/ 
      DATA ((A_EXTY(I,J,5),I=1,mby),j=1,nrh)/ .9410E+00, .9410E+00, 
     1   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     2   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     3   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     4   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     5   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     6   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     7   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     8   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     9   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     .   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     1   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     2   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     3   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     4   .9848E+00, .1008E+01, .1040E+01, .9410E+00, .9410E+00, 
     5   .9410E+00, .9410E+00, .9410E+00, .9556E+00, .9700E+00, 
     6   .9848E+00, .1008E+01, .1040E+01/ 
      DATA ((A_ASYY(I,J,5),I=1,mby),j=1,nrh)/ .8661E+00, .8661E+00, 
     1   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     2   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     3   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     4   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     5   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     6   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     7   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     8   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     9   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     .   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     1   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     2   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     3   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     4   .7654E+00, .7285E+00, .6931E+00, .8661E+00, .8661E+00, 
     5   .8661E+00, .8661E+00, .8661E+00, .8265E+00, .7970E+00, 
     6   .7654E+00, .7285E+00, .6931E+00/ 
      DATA ((A_SSAY(I,J,6),I=1,mby),j=1,nrh)/ .5631E+00, .5631E+00, 
     1   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     2   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     3   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     4   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     5   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     6   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     7   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     8   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     9   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     .   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     1   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     2   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     3   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     4   .6988E+00, .7839E+00, .8715E+00, .5631E+00, .5631E+00, 
     5   .5631E+00, .5631E+00, .5631E+00, .6011E+00, .6403E+00, 
     6   .6988E+00, .7839E+00, .8715E+00/ 
      DATA ((A_EXTY(I,J,6),I=1,mby),j=1,nrh)/ .9650E+00, .9650E+00, 
     1   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     2   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     3   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     4   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     5   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     6   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     7   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     8   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     9   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     .   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     1   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     2   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     3   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     4   .9916E+00, .1004E+01, .1019E+01, .9650E+00, .9650E+00, 
     5   .9650E+00, .9650E+00, .9650E+00, .9749E+00, .9831E+00, 
     6   .9916E+00, .1004E+01, .1019E+01/ 
      DATA ((A_ASYY(I,J,6),I=1,mby),j=1,nrh)/ .9183E+00, .9183E+00, 
     1   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     2   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     3   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     4   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     5   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     6   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     7   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     8   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     9   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     .   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     1   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     2   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     3   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     4   .8466E+00, .8097E+00, .7725E+00, .9183E+00, .9183E+00, 
     5   .9183E+00, .9183E+00, .9183E+00, .8957E+00, .8745E+00, 
     6   .8466E+00, .8097E+00, .7725E+00/ 
      DATA ((A_SSAY(I,J,7),I=1,mby),j=1,nrh)/ .5495E+00, .5495E+00, 
     1   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     2   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     3   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     4   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     5   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     6   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     7   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     8   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     9   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     .   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     1   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     2   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     3   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     4   .6141E+00, .6914E+00, .7949E+00, .5495E+00, .5495E+00, 
     5   .5495E+00, .5495E+00, .5495E+00, .5603E+00, .5775E+00, 
     6   .6141E+00, .6914E+00, .7949E+00/ 
      DATA ((A_EXTY(I,J,7),I=1,mby),j=1,nrh)/ .9779E+00, .9779E+00, 
     1   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     2   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     3   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     4   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     5   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     6   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     7   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     8   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     9   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     .   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     1   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     2   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     3   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     4   .9948E+00, .1002E+01, .1012E+01, .9779E+00, .9779E+00, 
     5   .9779E+00, .9779E+00, .9779E+00, .9839E+00, .9894E+00, 
     6   .9948E+00, .1002E+01, .1012E+01/ 
      DATA ((A_ASYY(I,J,7),I=1,mby),j=1,nrh)/ .9364E+00, .9364E+00, 
     1   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     2   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     3   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     4   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     5   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     6   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     7   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     8   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     9   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     .   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     1   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     2   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     3   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     4   .9026E+00, .8702E+00, .8309E+00, .9364E+00, .9364E+00, 
     5   .9364E+00, .9364E+00, .9364E+00, .9298E+00, .9204E+00, 
     6   .9026E+00, .8702E+00, .8309E+00/ 
      DATA ((A_SSAY(I,J,8),I=1,mby),j=1,nrh)/ .5507E+00, .5507E+00, 
     1   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     2   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     3   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     4   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     5   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     6   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     7   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     8   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     9   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     .   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     1   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     2   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     3   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     4   .5663E+00, .6106E+00, .6996E+00, .5507E+00, .5507E+00, 
     5   .5507E+00, .5507E+00, .5507E+00, .5512E+00, .5542E+00, 
     6   .5663E+00, .6106E+00, .6996E+00/ 
      DATA ((A_EXTY(I,J,8),I=1,mby),j=1,nrh)/ .9859E+00, .9859E+00, 
     1   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     2   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     3   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     4   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     5   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     6   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     7   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     8   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     9   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     .   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     1   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     2   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     3   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     4   .9967E+00, .1002E+01, .1007E+01, .9859E+00, .9859E+00, 
     5   .9859E+00, .9859E+00, .9859E+00, .9896E+00, .9932E+00, 
     6   .9967E+00, .1002E+01, .1007E+01/ 
      DATA ((A_ASYY(I,J,8),I=1,mby),j=1,nrh)/ .9401E+00, .9401E+00, 
     1   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     2   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     3   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     4   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     5   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     6   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     7   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     8   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     9   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     .   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     1   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     2   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     3   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     4   .9327E+00, .9136E+00, .8795E+00, .9401E+00, .9401E+00, 
     5   .9401E+00, .9401E+00, .9401E+00, .9399E+00, .9385E+00, 
     6   .9327E+00, .9136E+00, .8795E+00/ 
      DATA ((A_EXTY(I,J,9),I=1,mby),J=1,1)/ 0.9429E+00, 0.9453E+00, 
     1   0.9513E+00, 0.9558E+00, 0.9595E+00, 0.9650E+00, 0.9753E+00, 
     2   0.9866E+00, 0.9992E+00, 0.1015E+01/ 
      DATA ((A_SSAY(I,J,9),I=1,mby),J=1,1)/ 0.4624E+00, 0.5098E+00, 
     1   0.6053E+00, 0.6511E+00, 0.6674E+00, 0.6750E+00, 0.6918E+00, 
     2   0.7100E+00, 0.7289E+00, 0.7486E+00/ 
      DATA ((A_ASYY(I,J,9),I=1,mby),J=1,1)/ 0.9788E+00, 0.9585E+00, 
     1   0.9122E+00, 0.8872E+00, 0.8773E+00, 0.8720E+00, 0.8597E+00, 
     2   0.8458E+00, 0.8317E+00, 0.8163E+00/ 
      DATA ((A_EXTY(I,J,10),I=1,mby),J=1,8)/ 0.2636E+01, 0.2534E+01, 
     1   0.2300E+01, 0.2142E+01, 0.2022E+01, 0.1846E+01, 0.1535E+01, 
     2   0.1261E+01, 0.1015E+01, 0.7892E+00, 0.2564E+01, 0.2466E+01, 
     3   0.2240E+01, 0.2089E+01, 0.1975E+01, 0.1808E+01, 0.1513E+01, 
     4   0.1253E+01, 0.1015E+01, 0.7959E+00, 0.2500E+01, 0.2407E+01, 
     5   0.2193E+01, 0.2049E+01, 0.1940E+01, 0.1781E+01, 0.1498E+01, 
     6   0.1246E+01, 0.1014E+01, 0.8001E+00, 0.2434E+01, 0.2346E+01, 
     7   0.2144E+01, 0.2008E+01, 0.1904E+01, 0.1754E+01, 0.1482E+01, 
     8   0.1240E+01, 0.1014E+01, 0.8042E+00, 0.2291E+01, 0.2215E+01, 
     9   0.2039E+01, 0.1920E+01, 0.1829E+01, 0.1694E+01, 0.1448E+01, 
     .   0.1225E+01, 0.1013E+01, 0.8130E+00, 0.2121E+01, 0.2059E+01, 
     1   0.1914E+01, 0.1815E+01, 0.1738E+01, 0.1623E+01, 0.1407E+01, 
     2   0.1208E+01, 0.1012E+01, 0.8240E+00, 0.1896E+01, 0.1852E+01, 
     3   0.1746E+01, 0.1673E+01, 0.1614E+01, 0.1525E+01, 0.1350E+01, 
     4   0.1182E+01, 0.1011E+01, 0.8404E+00, 0.1752E+01, 0.1718E+01, 
     5   0.1638E+01, 0.1581E+01, 0.1534E+01, 0.1460E+01, 0.1311E+01, 
     6   0.1165E+01, 0.1010E+01, 0.8518E+00/ 
      DATA ((A_SSAY(I,J,10),I=1,mby),J=1,8)/ 0.6646E+00, 0.7653E+00, 
     1   0.8981E+00, 0.9419E+00, 0.9569E+00, 0.9666E+00, 0.9685E+00, 
     2   0.9688E+00, 0.9633E+00, 0.9558E+00, 0.7687E+00, 0.8417E+00, 
     3   0.9345E+00, 0.9639E+00, 0.9735E+00, 0.9796E+00, 0.9808E+00, 
     4   0.9810E+00, 0.9776E+00, 0.9730E+00, 0.8036E+00, 0.8665E+00, 
     5   0.9456E+00, 0.9704E+00, 0.9784E+00, 0.9835E+00, 0.9844E+00, 
     6   0.9847E+00, 0.9820E+00, 0.9782E+00, 0.8288E+00, 0.8845E+00, 
     7   0.9538E+00, 0.9751E+00, 0.9819E+00, 0.9861E+00, 0.9871E+00, 
     8   0.9873E+00, 0.9850E+00, 0.9820E+00, 0.8696E+00, 0.9126E+00, 
     9   0.9655E+00, 0.9816E+00, 0.9868E+00, 0.9900E+00, 0.9907E+00, 
     .   0.9909E+00, 0.9895E+00, 0.9874E+00, 0.9009E+00, 0.9345E+00, 
     1   0.9749E+00, 0.9869E+00, 0.9906E+00, 0.9929E+00, 0.9935E+00, 
     2   0.9938E+00, 0.9928E+00, 0.9915E+00, 0.9319E+00, 0.9553E+00, 
     3   0.9831E+00, 0.9913E+00, 0.9939E+00, 0.9954E+00, 0.9959E+00, 
     4   0.9961E+00, 0.9956E+00, 0.9949E+00, 0.9465E+00, 0.9651E+00, 
     5   0.9870E+00, 0.9934E+00, 0.9954E+00, 0.9966E+00, 0.9970E+00, 
     6   0.9972E+00, 0.9968E+00, 0.9964E+00/ 
      DATA ((A_ASYY(I,J,10),I=1,mby),J=1,8)/ 0.7099E+00, 0.6998E+00, 
     1   0.6762E+00, 0.6623E+00, 0.6551E+00, 0.6486E+00, 0.6386E+00, 
     2   0.6267E+00, 0.6143E+00, 0.5985E+00, 0.7400E+00, 0.7344E+00, 
     3   0.7212E+00, 0.7131E+00, 0.7084E+00, 0.7033E+00, 0.6946E+00, 
     4   0.6841E+00, 0.6722E+00, 0.6570E+00, 0.7487E+00, 0.7443E+00, 
     5   0.7339E+00, 0.7276E+00, 0.7241E+00, 0.7201E+00, 0.7118E+00, 
     6   0.7021E+00, 0.6904E+00, 0.6755E+00, 0.7554E+00, 0.7517E+00, 
     7   0.7427E+00, 0.7373E+00, 0.7343E+00, 0.7310E+00, 0.7239E+00, 
     8   0.7151E+00, 0.7042E+00, 0.6902E+00, 0.7629E+00, 0.7609E+00, 
     9   0.7561E+00, 0.7529E+00, 0.7507E+00, 0.7478E+00, 0.7426E+00, 
     .   0.7349E+00, 0.7254E+00, 0.7123E+00, 0.7685E+00, 0.7678E+00, 
     1   0.7659E+00, 0.7644E+00, 0.7633E+00, 0.7615E+00, 0.7576E+00, 
     2   0.7512E+00, 0.7433E+00, 0.7320E+00, 0.7721E+00, 0.7725E+00, 
     3   0.7736E+00, 0.7740E+00, 0.7738E+00, 0.7732E+00, 0.7715E+00, 
     4   0.7669E+00, 0.7612E+00, 0.7518E+00, 0.7749E+00, 0.7755E+00, 
     5   0.7768E+00, 0.7777E+00, 0.7784E+00, 0.7789E+00, 0.7779E+00, 
     6   0.7750E+00, 0.7701E+00, 0.7619E+00/ 
      DATA ((A_EXTY(I,J,11),I=1,mby),J=1,1)/ 0.2564E+01, 0.2504E+01, 
     1   0.2357E+01, 0.2233E+01, 0.2108E+01, 0.1900E+01, 0.1552E+01, 
     2   0.1267E+01, 0.1017E+01, 0.8078E+00/ 
      DATA ((A_SSAY(I,J,11),I=1,mby),J=1,1)/ 0.3009E+00, 0.3045E+00, 
     1   0.3124E+00, 0.3138E+00, 0.3091E+00, 0.2954E+00, 0.2666E+00, 
     2   0.2384E+00, 0.2102E+00, 0.1789E+00/ 
      DATA ((A_ASYY(I,J,11),I=1,mby),J=1,1)/ 0.5324E+00, 0.5169E+00, 
     1   0.4811E+00, 0.4589E+00, 0.4449E+00, 0.4272E+00, 0.3957E+00, 
     2   0.3664E+00, 0.3375E+00, 0.3079E+00/ 
      DATA ((A_EXTY(I,J,12),I=1,mby),J=1,8)/ 0.8629E+00, 0.8715E+00, 
     1   0.8930E+00, 0.9073E+00, 0.9173E+00, 0.9311E+00, 0.9576E+00, 
     2   0.9787E+00, 0.9977E+00, 0.1004E+01, 0.8793E+00, 0.8838E+00, 
     3   0.8954E+00, 0.9047E+00, 0.9137E+00, 0.9279E+00, 0.9519E+00, 
     4   0.9771E+00, 0.9989E+00, 0.1019E+01, 0.8777E+00, 0.8831E+00, 
     5   0.8965E+00, 0.9062E+00, 0.9141E+00, 0.9262E+00, 0.9497E+00, 
     6   0.9737E+00, 0.9984E+00, 0.1021E+01, 0.8839E+00, 0.8886E+00, 
     7   0.9002E+00, 0.9086E+00, 0.9155E+00, 0.9268E+00, 0.9508E+00, 
     8   0.9735E+00, 0.9991E+00, 0.1024E+01, 0.8990E+00, 0.9021E+00, 
     9   0.9098E+00, 0.9160E+00, 0.9219E+00, 0.9322E+00, 0.9531E+00, 
     .   0.9738E+00, 0.9984E+00, 0.1025E+01, 0.9067E+00, 0.9113E+00, 
     1   0.9224E+00, 0.9293E+00, 0.9335E+00, 0.9397E+00, 0.9585E+00, 
     2   0.9766E+00, 0.9990E+00, 0.1023E+01, 0.9296E+00, 0.9329E+00, 
     3   0.9409E+00, 0.9461E+00, 0.9497E+00, 0.9550E+00, 0.9668E+00, 
     4   0.9813E+00, 0.9994E+00, 0.1018E+01, 0.9437E+00, 0.9462E+00, 
     5   0.9524E+00, 0.9566E+00, 0.9596E+00, 0.9641E+00, 0.9729E+00, 
     6   0.9851E+00, 0.9998E+00, 0.1013E+01/ 
      DATA ((A_SSAY(I,J,12),I=1,mby),J=1,8)/ 0.9998E+00, 0.9998E+00, 
     1   0.9998E+00, 0.9999E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9995E+00, 0.9998E+00, 
     3   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9995E+00, 0.9998E+00, 
     5   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     7   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     8   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     9   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     .   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     1   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     3   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     5   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01/ 
      DATA ((A_ASYY(I,J,12),I=1,mby),J=1,8)/ 0.7300E+00, 0.7254E+00, 
     1   0.7146E+00, 0.7078E+00, 0.7037E+00, 0.6998E+00, 0.6978E+00, 
     2   0.6951E+00, 0.6925E+00, 0.6965E+00, 0.7850E+00, 0.7850E+00, 
     3   0.7846E+00, 0.7830E+00, 0.7798E+00, 0.7750E+00, 0.7740E+00, 
     4   0.7708E+00, 0.7710E+00, 0.7724E+00, 0.8029E+00, 0.8004E+00, 
     5   0.7945E+00, 0.7909E+00, 0.7889E+00, 0.7866E+00, 0.7831E+00, 
     6   0.7799E+00, 0.7783E+00, 0.7800E+00, 0.8050E+00, 0.8040E+00, 
     7   0.8015E+00, 0.7996E+00, 0.7980E+00, 0.7950E+00, 0.7895E+00, 
     8   0.7869E+00, 0.7840E+00, 0.7844E+00, 0.8178E+00, 0.8169E+00, 
     9   0.8146E+00, 0.8127E+00, 0.8108E+00, 0.8072E+00, 0.8002E+00, 
     .   0.7963E+00, 0.7933E+00, 0.7912E+00, 0.8304E+00, 0.8287E+00, 
     1   0.8246E+00, 0.8224E+00, 0.8216E+00, 0.8198E+00, 0.8116E+00, 
     2   0.8061E+00, 0.8009E+00, 0.7985E+00, 0.8380E+00, 0.8380E+00, 
     3   0.8378E+00, 0.8367E+00, 0.8345E+00, 0.8309E+00, 0.8271E+00, 
     4   0.8189E+00, 0.8136E+00, 0.8086E+00, 0.8448E+00, 0.8449E+00, 
     5   0.8450E+00, 0.8444E+00, 0.8431E+00, 0.8406E+00, 0.8373E+00, 
     6   0.8299E+00, 0.8244E+00, 0.8185E+00/ 
      DATA ((A_EXTY(I,J,13),I=1,mby),J=1,8)/ 0.9630E+00, 0.9648E+00, 
     1   0.9695E+00, 0.9730E+00, 0.9758E+00, 0.9790E+00, 0.9821E+00, 
     2   0.9899E+00, 0.9980E+00, 0.1007E+01, 0.9727E+00, 0.9743E+00, 
     3   0.9783E+00, 0.9809E+00, 0.9829E+00, 0.9856E+00, 0.9907E+00, 
     4   0.9950E+00, 0.9993E+00, 0.1006E+01, 0.9761E+00, 0.9773E+00, 
     5   0.9805E+00, 0.9827E+00, 0.9844E+00, 0.9866E+00, 0.9897E+00, 
     6   0.9943E+00, 0.9997E+00, 0.1006E+01, 0.9762E+00, 0.9780E+00, 
     7   0.9824E+00, 0.9849E+00, 0.9862E+00, 0.9877E+00, 0.9924E+00, 
     8   0.9962E+00, 0.1000E+01, 0.1007E+01, 0.9811E+00, 0.9821E+00, 
     9   0.9844E+00, 0.9861E+00, 0.9874E+00, 0.9893E+00, 0.9931E+00, 
     .   0.9972E+00, 0.1000E+01, 0.1005E+01, 0.9831E+00, 0.9837E+00, 
     1   0.9852E+00, 0.9861E+00, 0.9868E+00, 0.9880E+00, 0.9925E+00, 
     2   0.9962E+00, 0.9997E+00, 0.1004E+01, 0.9858E+00, 0.9865E+00, 
     3   0.9882E+00, 0.9894E+00, 0.9904E+00, 0.9918E+00, 0.9948E+00, 
     4   0.9969E+00, 0.9997E+00, 0.1003E+01, 0.9901E+00, 0.9891E+00, 
     5   0.9872E+00, 0.9870E+00, 0.9887E+00, 0.9911E+00, 0.9895E+00, 
     6   0.9966E+00, 0.9997E+00, 0.1002E+01/ 
      DATA ((A_SSAY(I,J,13),I=1,mby),J=1,8)/ 0.9974E+00, 0.9980E+00, 
     1   0.9990E+00, 0.9994E+00, 0.9996E+00, 0.9998E+00, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9995E+00, 0.9995E+00, 
     3   0.9995E+00, 0.9997E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9994E+00, 0.9995E+00, 
     5   0.9997E+00, 0.9999E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9994E+00, 0.9995E+00, 
     7   0.9997E+00, 0.9999E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     8   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9996E+00, 0.9997E+00, 
     9   0.9998E+00, 0.9999E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     .   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9998E+00, 0.9998E+00, 
     1   0.9998E+00, 0.9999E+00, 0.9999E+00, 0.1000E+01, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9995E+00, 0.9998E+00, 
     3   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.9995E+00, 0.9998E+00, 
     5   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01/ 
      DATA ((A_ASYY(I,J,13),I=1,mby),J=1,8)/ 0.8085E+00, 0.8092E+00, 
     1   0.8108E+00, 0.8106E+00, 0.8085E+00, 0.8049E+00, 0.8045E+00, 
     2   0.8017E+00, 0.7964E+00, 0.7932E+00, 0.8452E+00, 0.8456E+00, 
     3   0.8467E+00, 0.8476E+00, 0.8486E+00, 0.8498E+00, 0.8505E+00, 
     4   0.8504E+00, 0.8469E+00, 0.8464E+00, 0.8483E+00, 0.8481E+00, 
     5   0.8480E+00, 0.8486E+00, 0.8500E+00, 0.8530E+00, 0.8581E+00, 
     6   0.8557E+00, 0.8506E+00, 0.8505E+00, 0.8373E+00, 0.8422E+00, 
     7   0.8535E+00, 0.8586E+00, 0.8582E+00, 0.8542E+00, 0.8541E+00, 
     8   0.8559E+00, 0.8570E+00, 0.8549E+00, 0.8441E+00, 0.8465E+00, 
     9   0.8523E+00, 0.8548E+00, 0.8546E+00, 0.8539E+00, 0.8614E+00, 
     .   0.8626E+00, 0.8604E+00, 0.8586E+00, 0.8390E+00, 0.8405E+00, 
     1   0.8443E+00, 0.8471E+00, 0.8495E+00, 0.8536E+00, 0.8641E+00, 
     2   0.8666E+00, 0.8629E+00, 0.8652E+00, 0.8385E+00, 0.8383E+00, 
     3   0.8382E+00, 0.8405E+00, 0.8455E+00, 0.8545E+00, 0.8586E+00, 
     4   0.8622E+00, 0.8680E+00, 0.8681E+00, 0.8070E+00, 0.8121E+00, 
     5   0.8244E+00, 0.8326E+00, 0.8384E+00, 0.8465E+00, 0.8575E+00, 
     6   0.8605E+00, 0.8671E+00, 0.8673E+00/ 
      DATA ((A_EXTY(I,J,14),I=1,mby),J=1,1)/ 0.1013E+01, 0.1007E+01, 
     1   0.9897E+00, 0.9760E+00, 0.9624E+00, 0.9367E+00, 0.8702E+00, 
     2   0.7911E+00, 0.6970E+00, 0.5919E+00/ 
      DATA ((A_SSAY(I,J,14),I=1,mby),J=1,1)/ 0.7841E+00, 0.7925E+00, 
     1   0.8134E+00, 0.8331E+00, 0.8543E+00, 0.8839E+00, 0.9211E+00, 
     2   0.9491E+00, 0.9647E+00, 0.9740E+00/ 
      DATA ((A_ASYY(I,J,14),I=1,mby),J=1,1)/ 0.7398E+00, 0.7359E+00, 
     1   0.7260E+00, 0.7185E+00, 0.7118E+00, 0.7018E+00, 0.6876E+00, 
     2   0.6759E+00, 0.6649E+00, 0.6521E+00/ 
      DATA ((A_EXTY(I,J,15),I=1,mby),J=1,1)/ 0.8958E+00, 0.8998E+00, 
     1   0.9097E+00, 0.9170E+00, 0.9230E+00, 0.9325E+00, 0.9520E+00, 
     2   0.9736E+00, 0.9984E+00, 0.1028E+01/ 
      DATA ((A_SSAY(I,J,15),I=1,mby),J=1,1)/ 0.5676E+00, 0.5719E+00, 
     1   0.5844E+00, 0.6014E+00, 0.6254E+00, 0.6655E+00, 0.7403E+00, 
     2   0.8148E+00, 0.8711E+00, 0.9093E+00/ 
      DATA ((A_ASYY(I,J,15),I=1,mby),J=1,1)/ 0.9144E+00, 0.9086E+00, 
     1   0.8938E+00, 0.8802E+00, 0.8654E+00, 0.8411E+00, 0.8007E+00, 
     2   0.7640E+00, 0.7372E+00, 0.7162E+00/ 
      DATA ((A_EXTY(I,J,16),I=1,mby),J=1,1)/ 0.9717E+00, 0.9730E+00, 
     1   0.9761E+00, 0.9783E+00, 0.9801E+00, 0.9828E+00, 0.9880E+00, 
     2   0.9935E+00, 0.9996E+00, 0.1007E+01/ 
      DATA ((A_SSAY(I,J,16),I=1,mby),J=1,1)/ 0.5462E+00, 0.5457E+00, 
     1   0.5448E+00, 0.5454E+00, 0.5477E+00, 0.5523E+00, 0.5714E+00, 
     2   0.6059E+00, 0.6601E+00, 0.7118E+00/ 
      DATA ((A_ASYY(I,J,16),I=1,mby),J=1,1)/ 0.9433E+00, 0.9447E+00, 
     1   0.9478E+00, 0.9490E+00, 0.9485E+00, 0.9466E+00, 0.9380E+00, 
     2   0.9212E+00, 0.8973E+00, 0.8741E+00/ 
      DATA ((A_EXTY(I,J,17),I=1,mby),J=1,1)/ 0.9243E+00, 0.9273E+00, 
     1   0.9348E+00, 0.9403E+00, 0.9447E+00, 0.9516E+00, 0.9657E+00, 
     2   0.9810E+00, 0.9986E+00, 0.1019E+01/ 
      DATA ((A_SSAY(I,J,17),I=1,mby),J=1,1)/ 0.5535E+00, 0.5553E+00, 
     1   0.5615E+00, 0.5724E+00, 0.5898E+00, 0.6205E+00, 0.6873E+00, 
     2   0.7635E+00, 0.8289E+00, 0.8763E+00/ 
      DATA ((A_ASYY(I,J,17),I=1,mby),J=1,1)/ 0.9371E+00, 0.9340E+00, 
     1   0.9257E+00, 0.9168E+00, 0.9057E+00, 0.8860E+00, 0.8473E+00, 
     2   0.8087E+00, 0.7784E+00, 0.7532E+00/ 
      DATA ((A_EXTY(I,J,18),I=1,mby),J=1,8)/ 0.1618E+01, 0.1603E+01, 
     1   0.1566E+01, 0.1534E+01, 0.1500E+01, 0.1438E+01, 0.1299E+01, 
     2   0.1155E+01, 0.1009E+01, 0.8521E+00, 0.1357E+01, 0.1352E+01, 
     3   0.1341E+01, 0.1329E+01, 0.1315E+01, 0.1285E+01, 0.1207E+01, 
     4   0.1115E+01, 0.1006E+01, 0.8812E+00, 0.1272E+01, 0.1271E+01, 
     5   0.1268E+01, 0.1263E+01, 0.1255E+01, 0.1235E+01, 0.1175E+01, 
     6   0.1100E+01, 0.1006E+01, 0.8928E+00, 0.1210E+01, 0.1211E+01, 
     7   0.1215E+01, 0.1215E+01, 0.1211E+01, 0.1198E+01, 0.1152E+01, 
     8   0.1089E+01, 0.1005E+01, 0.9023E+00, 0.1120E+01, 0.1125E+01, 
     9   0.1137E+01, 0.1143E+01, 0.1144E+01, 0.1139E+01, 0.1112E+01, 
     .   0.1069E+01, 0.1004E+01, 0.9180E+00, 0.1033E+01, 0.1040E+01, 
     1   0.1057E+01, 0.1067E+01, 0.1073E+01, 0.1077E+01, 0.1070E+01, 
     2   0.1047E+01, 0.1002E+01, 0.9381E+00, 0.9487E+00, 0.9569E+00, 
     3   0.9770E+00, 0.9906E+00, 0.1000E+01, 0.1012E+01, 0.1022E+01, 
     4   0.1020E+01, 0.1001E+01, 0.9637E+00, 0.9078E+00, 0.9155E+00, 
     5   0.9348E+00, 0.9483E+00, 0.9585E+00, 0.9724E+00, 0.9917E+00, 
     6   0.1002E+01, 0.1000E+01, 0.9817E+00/ 
      DATA ((A_SSAY(I,J,18),I=1,mby),J=1,8)/ 0.1000E+01, 0.1000E+01, 
     1   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     3   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     5   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     7   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     8   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     9   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     .   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     1   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     2   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     3   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     4   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     5   0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 
     6   0.1000E+01, 0.1000E+01, 0.1000E+01/ 
      DATA ((A_ASYY(I,J,18),I=1,mby),J=1,8)/ 0.6955E+00, 0.6978E+00, 
     1   0.7034E+00, 0.7076E+00, 0.7112E+00, 0.7165E+00, 0.7227E+00, 
     2   0.7231E+00, 0.7172E+00, 0.7080E+00, 0.7408E+00, 0.7444E+00, 
     3   0.7532E+00, 0.7586E+00, 0.7617E+00, 0.7650E+00, 0.7706E+00, 
     4   0.7717E+00, 0.7690E+00, 0.7634E+00, 0.7498E+00, 0.7529E+00, 
     5   0.7606E+00, 0.7656E+00, 0.7690E+00, 0.7729E+00, 0.7779E+00, 
     6   0.7797E+00, 0.7779E+00, 0.7743E+00, 0.7559E+00, 0.7590E+00, 
     7   0.7664E+00, 0.7709E+00, 0.7734E+00, 0.7762E+00, 0.7818E+00, 
     8   0.7840E+00, 0.7837E+00, 0.7803E+00, 0.7601E+00, 0.7631E+00, 
     9   0.7703E+00, 0.7748E+00, 0.7776E+00, 0.7810E+00, 0.7867E+00, 
     .   0.7894E+00, 0.7900E+00, 0.7891E+00, 0.7657E+00, 0.7684E+00, 
     1   0.7749E+00, 0.7789E+00, 0.7814E+00, 0.7842E+00, 0.7890E+00, 
     2   0.7928E+00, 0.7948E+00, 0.7949E+00, 0.7749E+00, 0.7765E+00, 
     3   0.7802E+00, 0.7825E+00, 0.7836E+00, 0.7853E+00, 0.7906E+00, 
     4   0.7940E+00, 0.7975E+00, 0.7987E+00, 0.7795E+00, 0.7808E+00, 
     5   0.7838E+00, 0.7856E+00, 0.7865E+00, 0.7876E+00, 0.7909E+00, 
     6   0.7937E+00, 0.7967E+00, 0.7995E+00/ 
!=====================================================================
      END 
 
c===================================================================
      block data cont_coef
      parameter (ncoef = 7, nband = 12)
      real aa0, aa1, aa2, aa3, aa4
 
      common /pcont_aa0/aa0(ncoef,nband)
      common /pcont_aa1/aa1(ncoef,nband)
      common /pcont_aa2/aa2(ncoef,nband)
      common /pcont_aa3/aa3(ncoef,nband)
      common /pcont_aa4/aa4(ncoef,nband)
      data aa0/ -1.00E+01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     1   0.000E+00, 0.000E+00, 3.457E+00, 9.977E-01, -6.202E-03, 
     2   1.822E+00, 8.200E+00, 8.343E-03, 1.413E-02, 2.327E+00, 
     3   1.011E+00, -7.050E-03, 1.368E+00, 1.508E+01, -1.747E-03, 
     4   1.271E-01, 6.139E+00, 1.015E+00, -1.505E-02, 5.499E-01, 
     5   6.284E+00, -7.911E-03, 5.933E-01, 9.234E+00, 1.003E+00, 
     6   -2.254E-02, 1.204E-01, 3.237E-01, -1.697E-03, 9.105E-01, 
     7   9.426E+00, 9.999E-01, -2.438E-02, 2.410E-02, -7.441E-01, 
     8   4.155E-05, 9.880E-01, 8.408E+00, 1.000E+00, -2.319E-02, 
     9   2.111E-02, -8.137E-01, 6.080E-06, 9.933E-01, 7.285E+00, 
     .   1.001E+00, -2.049E-02, 5.285E-02, -5.057E-01, -3.719E-04, 
     1   9.647E-01, 1.096E+00, 1.008E+00, -6.852E-03, 1.021E+00, 
     2   8.523E+00, -5.127E-03, 2.513E-01, -1.00E+01, 0.000E+00, 
     3   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     4   -1.00E+01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     5   0.000E+00, 0.000E+00, -1.00E+01, 0.000E+00, 0.000E+00, 
     6   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/ 
      data aa1/ 6.762E+00, 9.996E-01, -7.259E-03, 9.977E-01, 1.603E+01, 
     1   1.709E-01, -1.394E-03, 5.261E+00, 9.971E-01, -6.219E-03, 
     2   9.868E-01, 1.561E+01, 1.275E-01, 2.789E-03, 3.235E+00, 
     3   9.872E-01, -6.404E-03, 8.690E-01, 2.107E+01, 1.222E-01, 
     4   6.908E-02, 5.942E+00, 9.982E-01, -1.368E-02, 3.916E-01, 
     5   9.076E+00, 2.730E-02, 5.256E-01, 9.217E+00, 1.003E+00, 
     6   -2.221E-02, 8.940E-02, 6.127E-01, -6.341E-04, 8.978E-01, 
     7   9.478E+00, 1.000E+00, -2.441E-02, 1.837E-02, -7.174E-01, 
     8   -8.151E-06, 9.865E-01, 8.451E+00, 1.000E+00, -2.323E-02, 
     9   1.559E-02, -8.067E-01, 5.546E-06, 9.932E-01, 7.281E+00, 
     .   1.001E+00, -2.035E-02, 4.216E-02, -4.622E-01, -4.440E-04, 
     1   9.615E-01, 1.081E+00, 9.983E-01, -4.259E-03, 8.882E-01, 
     2   9.883E+00, 1.679E-02, 6.666E-02, 3.738E+00, 9.973E-01, 
     3   -3.824E-03, 9.842E-01, 8.171E+00, 5.865E-02, 4.520E-03, 
     4   2.879E+00, 9.953E-01, -3.841E-03, 9.721E-01, 7.980E+00, 
     5   5.344E-02, 1.106E-02, 7.045E-01, 1.008E+00, -5.968E-03, 
     6   7.708E-01, 7.788E+00, -5.452E-03, 1.683E-01/ 
      data aa2/ 6.405E+00, 9.986E-01, -7.277E-03, 9.963E-01, 1.304E+01, 
     1   1.588E-01, -3.569E-04, 5.039E+00, 9.967E-01, -6.107E-03, 
     2   9.820E-01, 1.819E+01, 1.188E-01, 4.073E-03, 3.333E+00, 
     3   9.832E-01, -6.811E-03, 8.141E-01, 2.107E+01, 1.416E-01, 
     4   1.102E-01, 6.498E+00, 1.000E+00, -1.486E-02, 3.513E-01, 
     5   7.720E+00, 2.083E-02, 5.794E-01, 9.328E+00, 1.003E+00, 
     6   -2.251E-02, 8.079E-02, 4.519E-01, -5.421E-04, 9.084E-01, 
     7   9.448E+00, 9.999E-01, -2.438E-02, 1.711E-02, -7.349E-01, 
     8   5.084E-05, 9.879E-01, 8.424E+00, 1.000E+00, -2.318E-02, 
     9   1.598E-02, -8.060E-01, -2.452E-05, 9.931E-01, 7.343E+00, 
     .   1.001E+00, -2.051E-02, 3.936E-02, -5.032E-01, -3.829E-04, 
     1   9.655E-01, 1.962E+00, 1.004E+00, -6.678E-03, 7.064E-01, 
     2   9.220E+00, 5.062E-03, 2.399E-01, 3.616E+00, 9.972E-01, 
     3   -3.827E-03, 9.825E-01, 8.823E+00, 5.525E-02, 5.027E-03, 
     4   2.481E+00, 9.947E-01, -3.944E-03, 9.568E-01, 9.687E+00, 
     5   5.039E-02, 1.942E-02, 2.141E+00, 1.009E+00, -1.016E-02, 
     6   5.728E-01, 5.100E+00, -5.114E-03, 3.688E-01/ 
      data aa3/ 6.899E+00, 9.989E-01, -7.230E-03, 9.979E-01, 1.680E+01, 
     1   2.287E-01, -9.470E-04, 5.300E+00, 9.971E-01, -6.528E-03, 
     2   9.864E-01, 1.615E+01, 1.283E-01, 2.595E-03, 3.351E+00, 
     3   9.867E-01, -6.864E-03, 8.681E-01, 2.138E+01, 1.249E-01, 
     4   6.981E-02, 5.981E+00, 9.983E-01, -1.384E-02, 3.944E-01, 
     5   9.316E+00, 2.791E-02, 5.220E-01, 9.253E+00, 1.003E+00, 
     6   -2.230E-02, 9.128E-02, 6.481E-01, -6.383E-04, 8.963E-01, 
     7   9.555E+00, 1.000E+00, -2.460E-02, 1.921E-02, -6.885E-01, 
     8   -3.569E-05, 9.855E-01, 8.495E+00, 1.000E+00, -2.333E-02, 
     9   1.484E-02, -7.828E-01, -1.622E-06, 9.927E-01, 7.328E+00, 
     .   1.001E+00, -2.046E-02, 3.898E-02, -5.016E-01, -3.945E-04, 
     1   9.649E-01, 8.158E-01, 9.974E-01, -3.582E-03, 8.775E-01, 
     2   9.817E+00, 1.744E-02, 7.603E-02, 3.482E+00, 9.971E-01, 
     3   -3.205E-03, 9.827E-01, 8.142E+00, 5.638E-02, 5.275E-03, 
     4   3.262E+00, 9.969E-01, -4.671E-03, 9.760E-01, 8.702E+00, 
     5   4.119E-02, 8.187E-03, 1.006E+00, 1.008E+00, -6.560E-03, 
     6   8.004E-01, 8.480E+00, -5.247E-03, 1.390E-01/ 
      data aa4/ 6.883E+00, 1.007E+00, -7.231E-03, 1.006E+00, 1.711E+01, 
     1   2.401E-01, -8.842E-04, 5.128E+00, 1.006E+00, -6.531E-03, 
     2   9.905E-01, 1.907E+01, 1.219E-01, 4.032E-03, 3.529E+00, 
     3   9.903E-01, -7.494E-03, 8.184E-01, 2.150E+01, 1.434E-01, 
     4   1.137E-01, 6.623E+00, 1.008E+00, -1.525E-02, 3.557E-01, 
     5   7.867E+00, 2.121E-02, 5.826E-01, 9.438E+00, 1.010E+00, 
     6   -2.281E-02, 8.279E-02, 4.923E-01, -4.263E-04, 9.138E-01, 
     7   9.574E+00, 1.005E+00, -2.472E-02, 1.761E-02, -7.078E-01, 
     8   1.032E-04, 9.917E-01, 8.530E+00, 1.007E+00, -2.347E-02, 
     9   1.612E-02, -7.979E-01, 1.613E-05, 1.000E+00, 7.429E+00, 
     .   1.006E+00, -2.076E-02, 3.668E-02, -5.353E-01, -2.511E-04, 
     1   9.738E-01, 1.838E+00, 1.013E+00, -6.037E-03, 6.632E-01, 
     2   8.190E+00, 1.126E-04, 2.812E-01, 3.284E+00, 9.994E-01, 
     3   -3.042E-03, 9.830E-01, 8.944E+00, 4.588E-02, 5.890E-03, 
     4   2.977E+00, 9.987E-01, -5.050E-03, 9.663E-01, 9.910E+00, 
     5   4.875E-02, 1.525E-02, 2.177E+00, 1.010E+00, -1.035E-02, 
     6   6.566E-01, 6.613E+00, -5.513E-03, 2.837E-01/ 
      end 
 
 
c===================================================================
