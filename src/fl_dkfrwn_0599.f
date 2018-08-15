      subroutine ck_dkfr_win_hyb(ib, iga, hka, tg)
      parameter (maxig = 5, ngas = 6, maxjp = 10)
      INCLUDE 'fl_radparams.inc'
C#        include 'rad_0698.h'
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
 
      real fkg(101)
      common /band11/hk11(3), c11h2o(3,19,3), c11ch4(3,19), c11n2o(3,19)
      common /band12/hk12(5), c12o3(3,19,5), c12h2o(3,19)
      common /band13/hk13(2), c13h2o(3,19,2)
 
      real tg(nvx)
 
      real tgall(nvx)
 
      integer ngv(11:13,6), ng(6)
      save ngv, ng
 
      integer igs1(maxjp,11:13), igs2(maxjp,11:13), igs3(maxjp,11:13),
     1    igs4(maxjp,11:13), igs5(maxjp,11:13), igs6(maxjp,11:13)
      save igs1, igs2, igs3, igs4, igs5, igs6
      real tgx(nvx,maxig,ngas), hkx(maxig,ngas)
      save tgx, hkx
      character glab(6)*5
      integer j1, j2, j3
      real r1, r2
      data glab/ 'H2O', 'O3', 'CO2', 'N2O', 'CH4', 'CFC'/ 
      ngv(1+10,1) = 5
      ngv(2+10,1) = 1
      ngv(3+10,1) = 2
      ngv(1+10,2) = 1
      ngv(2+10,2) = 5
      ngv(3+10,2) = 1
      ngv(1+10,3) = 0
      ngv(2+10,3) = 1
      ngv(3+10,3) = 1
      ngv(1+10,4) = 1
      ngv(2+10,4) = 0
      ngv(3+10,4) = 0
      ngv(1+10,5) = 1
      ngv(2+10,5) = 0
      ngv(3+10,5) = 0
      ngv(1+10,6) = 1
      ngv(2+10,6) = 1
      ngv(3+10,6) = 1
! 5,5,2=12
!-------------------------------------------------------------------
!
      if (iga .eq. 1) then
         do j2 = 1, 6
            do j1 = 1, 5
               hkx(j1,j2) = 1.000
            end do
         end do
         do j3 = 1, 6
            do j2 = 1, 5
               do j1 = 1, 100
                  tgx(j1,j2,j3) = 0.000
               end do
            end do
         end do
 
         do kgas = 1, 6
            if (ngv(ib,kgas) .ne. 0) then
!      print*, 'HYB',glab(kgas)
               do ig = 1, ngv(ib,kgas)
 
                  if (kgas .eq. 1) then
                     if (ib .eq. 12) then
                        call qki (c12h2o, fkg)
                        call qoph2o (fkg, tgx(1,ig,1))
                        hkx(ig,1) = 1.0
                     else if (ib .eq. 13) then
                        call qki (c13h2o(1,1,ig), fkg)
                        call qoph2o (fkg, tgx(1,ig,1))
                        hkx(ig,1) = hk13(ig)
                     else if (ib .eq. 11) then
                        call ck_dkfr_win_h2o (ib, ig, hkx(ig,1), tgx(1,
     1                     ig,1))
                     endif
 
                  else if (kgas .eq. 2) then
                     if (ib .eq. 11) then
                        call ck_dkfr_win_o3_1k (ib, ig, hkx(ig,2), tgx(1
     1                     ,ig,2))
 
                     else if (ib .eq. 12) then
                        call qkio3 (c12o3(1,1,ig), fkg)
                        call qopo3i (fkg, tgx(1,ig,2))
                        hkx(ig,2) = hk12(ig)
 
                     else if (ib .eq. 13) then
                        call ck_dkfr_win_o3 (ib, ig, hkx(ig,2), tgx(1,ig
     1                     ,2))
                     endif
 
                  else if (kgas .eq. 3) then
                     call ck_dkfr_win_co2_thin (ib, ig, hkx(ig,3), tgx(1
     1                  ,ig,3))
 
                  else if (kgas .eq. 4) then
 
                     call qki (c11n2o, fkg)
                     call qopn2o (fkg, tgx(1,ig,4))
                     do j1 = 1, nv
                        tgx(j1,ig,4) = tgx(j1,ig,4)/0.28*umn2o
                     end do
                     hkx(ig,4) = 1.0
 
!     call ck_dkfr_win_n2o(ib,ig, hkx(ig,4 ),tgx(1,ig,4) )
 
                  else if (kgas .eq. 5) then
 
                     call qki (c11ch4, fkg)
                     call qopch4 (fkg, tgx(1,ig,5))
                     do j1 = 1, nv
                        tgx(j1,ig,5) = tgx(j1,ig,5)/1.6*umch4
                     end do
                     hkx(ig,5) = 1.0
 
!        call ck_dkfr_win_ch4(ib,ig, hkx(ig,5 ),tgx(1,ig,5) )
 
                  else if (kgas .eq. 6) then
                     call ck_dkfr_win_cfc(ib,ig,hkx(ig,6),tgx(1,ig,6))
                  endif
!     print*, tgx(1:nv,ig,kgas)
               end do
!     print*, hkx(1:ngv(ib,kgas),kgas )
 
            endif
         end do
 
         kk = 0
         sumhka1 = 0
         do kgas = 1, 6
            ng(kgas) = ngv(ib,kgas)
            if (ngv(ib,kgas) .eq. 0) ng(kgas) = 1
         end do
         tgsumw = 0
         do ig1 = 1, ng(1)
            do ig2 = 1, ng(2)
               do ig3 = 1, ng(3)
                  do ig4 = 1, ng(4)
                     do ig5 = 1, ng(5)
                        do ig6 = 1, ng(6)
                           kk = kk + 1
                           igs1(kk,ib) = ig1
                           igs2(kk,ib) = ig2
                           igs3(kk,ib) = ig3
                           igs4(kk,ib) = ig4
                           igs5(kk,ib) = ig5
                           igs6(kk,ib) = ig6
 
                           hka1 = hkx(ig1,1)*hkx(ig2,2)*hkx(ig3,3)*hkx(
     1                        ig4,4)*hkx(ig5,5)*hkx(ig6,6)
                           sumhka1 = sumhka1 + hka1
                           do j1 = 1, 100
                              tgall(j1) = 0.0
                           end do
 
                           do j1 = 1, nv
                              tgall(j1) = tgx(j1,ig1,1) + tgx(j1,ig2,2)
     1                            + tgx(j1,ig3,3) + tgx(j1,ig4,4) + tgx(
     2                           j1,ig5,5) + tgx(j1,ig6,6)
                           end do
                           r1 = 0
                           do j1 = 1, nv
                              r1 = r1 + tgall(j1)
                           end do
                           tgsum = r1
                           r2 = 0
                           do j1 = 1, nv
                              r2 = r2 + tgall(j1)*hka1
                           end do
                           tgsumw = tgsumw + r2
 
!      print'(I5,6I2,6f8.5,2f10.5,2x,4f8.4)'
!     &       ,kk,ig1,ig2,ig3,ig4,ig5,ig6,
!     &         hkx(ig1,1),hkx(ig2,2),hkx(ig3,3),
!     &         hkx(ig4,4),hkx(ig5,5),hkx(ig6,6),
!     &  hka1,sumhka1,tgsum,tgsumw
 
                        end do
                     end do
                  end do
               end do
            end do
         end do
      endif
      hka = hkx(igs1(iga,ib),1)*hkx(igs2(iga,ib),2)*hkx(igs3(iga,ib),3)*
     1   hkx(igs4(iga,ib),4)*hkx(igs5(iga,ib),5)*hkx(igs6(iga,ib),6)
      do j1 = 1, nv
         tg(j1) = tgx(j1,igs1(iga,ib),1) + tgx(j1,igs2(iga,ib),2) + tgx(
     1      j1,igs3(iga,ib),3) + tgx(j1,igs4(iga,ib),4) + tgx(j1,igs5(
     2      iga,ib),5) + tgx(j1,igs6(iga,ib),6)
      end do
      return 
      end 
c *********************************************************************
      subroutine ck_dkfr_win_co2_thin(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x), dz
      common /thick/ dz(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
      integer j1
 
      do m = 1, nv
cbm         dz = gpt(pp(m+1),pp(m),pt(m+1),pt(m),ph(m+1),ph(m))
! ! CO2 UNITS
!       u(m)=pp(m)*1.01325e+06/(1.3805e-16*pt(m))*3.50e-04*
!     *          1.0e+05/6.023e+23*44.00995
 
cbm         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umco2*1.0E-06*dz*
cbm     1      1.0E+05/6.023E+23*44.00995
         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umco2*1.0E-06*dz(m)*
     1      1.0E+05/6.023E+23*44.00995
      end do
 
      if (ig .gt. 1) stop ' co2thin ig>1'
      do j1 = 1, nv
         tg(j1) = 0.0
      end do
 
      if (ib .eq. 12) then
         hk = 1
         rk1 = 0.03479137
         do m = 1, nv
            tg(m) = 0.837764*rk1*(1.0000 + 3.4268E-02*(pt(m)-250.0)+
     1         3.7401E-04*(pt(m)-250.0)**2)*u(m)
         end do
 
      else if (ib .eq. 13) then
         hk = 1.0
         rk1 = 0.01407591
         do m = 1, nv
            tg(m) = 0.865828*rk1*(1.0000 + 3.7154E-02*(pt(m)-250.0)+
     1         4.3205E-04*(pt(m)-250.0)**2)*u(m)
         end do
      endif
 
 
 
      return 
      end 
c *********************************************************************
      subroutine ck_dkfr_win_o3_1k(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      integer j1
 
c     o3 1100-1250 cm^{-1} BAND 11
 
      hk = 1.000000
      rk1 = 21.305812
      do j1 = 1, nv
         tg(j1) = 0.0
      end do
      if (ig .gt. 1) stop ' ig >1 o3 '
      if (ib .eq. 11) then
         do m = 1, nv
! ! O3 UNITS
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
cbm         u(m) = 1.02*(po(m)+po(m+1))*0.5*(pp(m+1)-pp(m))
            u(m) = (po(m)+po(m+1))*0.5*( pa(i)*0.1 )
 
            tg(m) = 0.778043*rk1*(1.0000 + 1.4500E-03*(pt(m)-250.0)+
     1         6.3265E-06*(pt(m)-250.0)**2)*u(m)
         end do
      endif
      return 
      end 
c *********************************************************************
!================================================================
!================================================================
      subroutine ck_dkfr_win_all(ib, iga, hka, tg)
      parameter (maxig = 5, ngas = 6, maxjp = 80)
      INCLUDE 'fl_radparams.inc'
c#        include 'rad_0698.h'
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
      real tg(nvx)
 
      real tgall(nvx)
 
      integer ngvkr(11:13,6), ng(6)
      save ngvkr, ng
      integer igs1(maxjp,11:13), igs2(maxjp,11:13), igs3(maxjp,11:13), 
     1   igs4(maxjp,11:13), igs5(maxjp,11:13), igs6(maxjp,11:13)
      save igs1, igs2, igs3, igs4, igs5, igs6
      real tgx(nvx,maxig,ngas),hkx(maxig,ngas)
      save tgx, hkx
      character glab(6)*5
      integer j1, j2, j3
      real r1, r2
      data glab/ 'H2O', 'O3', 'CO2', 'N2O', 'CH4', 'CFC'/ 
      ngvkr(1+10,1) = 5
      ngvkr(2+10,1) = 4
      ngvkr(3+10,1) = 5
      ngvkr(1+10,2) = 2
      ngvkr(2+10,2) = 5
      ngvkr(3+10,2) = 1
      ngvkr(1+10,3) = 0
      ngvkr(2+10,3) = 2
      ngvkr(3+10,3) = 2
      ngvkr(1+10,4) = 2
      ngvkr(2+10,4) = 0
      ngvkr(3+10,4) = 0
      ngvkr(1+10,5) = 4
      ngvkr(2+10,5) = 0
      ngvkr(3+10,5) = 0
      ngvkr(1+10,6) = 1
      ngvkr(2+10,6) = 1
      ngvkr(3+10,6) = 1
!-------------------------------------------------------------------
!
      if (iga .eq. 1) then
         do j2 = 1, 6
            do j1 = 1, 5
               hkx(j1,j2) = 1.000
            end do
         end do
         do j3 = 1, 6
            do j2 = 1, 5
               do j1 = 1, 100
                  tgx(j1,j2,j3) = 0.000
               end do
            end do
         end do
 
         do kgas = 1, 6
            if (ngvkr(ib,kgas) .ne. 0) then
!      print*, 'Kratz',glab(kgas)
               do ig = 1, ngvkr(ib,kgas)
 
                  if (kgas .eq. 1) then
                     call ck_dkfr_win_h2o(ib,ig,hkx(ig,1),tgx(1,ig,1))
                  else if (kgas .eq. 2) then
                     call ck_dkfr_win_o3(ib,ig,hkx(ig,2),tgx(1,ig,2))
                  else if (kgas .eq. 3) then
                     call ck_dkfr_win_co2(ib,ig,hkx(ig,3),tgx(1,ig,3))
                  else if (kgas .eq. 4) then
                     call ck_dkfr_win_n2o(ib,ig,hkx(ig,4),tgx(1,ig,4))
                  else if (kgas .eq. 5) then
                     call ck_dkfr_win_ch4(ib,ig,hkx(ig,5),tgx(1,ig,5))
                  else if (kgas .eq. 6) then
                     call ck_dkfr_win_cfc(ib,ig,hkx(ig,6),tgx(1,ig,6))
                  endif
!     print*, tgx(1:nv,ig,kgas)
               end do
!     print*, hkx(1:ngvkr(ib,kgas),kgas )
 
            endif
         end do
 
         kk = 0
         sumhka1 = 0
         do kgas = 1, 6
            ng(kgas) = ngvkr(ib,kgas)
            if (ngvkr(ib,kgas) .eq. 0) ng(kgas) = 1
         end do
         tgsumw = 0
         do ig1 = 1, ng(1)
            do ig2 = 1, ng(2)
               do ig3 = 1, ng(3)
                  do ig4 = 1, ng(4)
                     do ig5 = 1, ng(5)
                        do ig6 = 1, ng(6)
                           kk = kk + 1
                           igs1(kk,ib) = ig1
                           igs2(kk,ib) = ig2
                           igs3(kk,ib) = ig3
                           igs4(kk,ib) = ig4
                           igs5(kk,ib) = ig5
                           igs6(kk,ib) = ig6
 
                           hka1 = hkx(ig1,1)*hkx(ig2,2)*hkx(ig3,3)*hkx(
     1                        ig4,4)*hkx(ig5,5)*hkx(ig6,6)
                           sumhka1 = sumhka1 + hka1
                           do j1 = 1, 100
                              tgall(j1) = 0.0
                           end do
 
                           do j1 = 1, nv
                              tgall(j1) = tgx(j1,ig1,1) + tgx(j1,ig2,2)
     1                            + tgx(j1,ig3,3) + tgx(j1,ig4,4) + tgx(
     2                           j1,ig5,5) + tgx(j1,ig6,6)
                           end do
                           r1 = 0
                           do j1 = 1, nv
                              r1 = r1 + tgall(j1)
                           end do
                           tgsum = r1
                           r2 = 0
                           do j1 = 1, nv
                              r2 = r2 + tgall(j1)*hka1
                           end do
                           tgsumw = tgsumw + r2
 
!      print'(I5,6I2,6f8.5,2f10.5,2x,4f8.4)'
!     &       ,kk,ig1,ig2,ig3,ig4,ig5,ig6,
!     &         hkx(ig1,1),hkx(ig2,2),hkx(ig3,3),
!     &         hkx(ig4,4),hkx(ig5,5),hkx(ig6,6),
!     &  hka1,sumhka1,tgsum,tgsumw
 
                        end do
                     end do
                  end do
               end do
            end do
         end do
      endif
      hka = hkx(igs1(iga,ib),1)*hkx(igs2(iga,ib),2)*hkx(igs3(iga,ib),3)*
     1   hkx(igs4(iga,ib),4)*hkx(igs5(iga,ib),5)*hkx(igs6(iga,ib),6)
      do j1 = 1, nv
         tg(j1) = tgx(j1,igs1(iga,ib),1) + tgx(j1,igs2(iga,ib),2) + tgx(
     1      j1,igs3(iga,ib),3) + tgx(j1,igs4(iga,ib),4) + tgx(j1,igs5(
     2      iga,ib),5) + tgx(j1,igs6(iga,ib),6)
      end do
      return 
      end 
c *********************************************************************
      subroutine ck_dkfr_win_h2o(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp,pt,ph,po,tg(nvx), u(nv1x),coefk(3,19)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
      common /band13dkfr/hk13h2o(5), c13h2o(5,3,19)
      common /band12dkfr/hk12h2o(4), c12h2o(4,3,19)
      common /band11dkfr/hk11h2o(5), c11h2o(5,3,19)
      integer j1, j2
 
      if (ib .eq. 13) then
         ni = 5
         rk1 = 0.00207622908
         rk = 6.0
         hk = hk13h2o(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c13h2o(ig,j1,j2)
            end do
         end do
      else if (ib .eq. 12) then
         ni = 4
         rk1 = 0.00238376099
         rk = 10.0
         hk = hk12h2o(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c12h2o(ig,j1,j2)
            end do
         end do
      else if (ib .eq. 11) then
         ni = 5
         rk1 = 0.004186672
         rk = 9.0
         hk = hk11h2o(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c11h2o(ig,j1,j2)
            end do
         end do
      else
         stop ' bad band ib not 11,12,13'
      endif
      if (ig .gt. ni) stop ' BAD ig'
 
      do m = 1, nv
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
cbm      u(m) = 1.02*(ph(m)+ph(m+1))*0.5*(pp(m+1)-pp(m))
         u(m) = (ph(m)+ph(m+1))*0.5*( pa(i)*0.1 )
      end do
 
      call stp_correction (ig, rk1, rk, coefk, u, tg)
 
      return 
      end 
 
c===============================================================
 
 
c *********************************************************************
      subroutine ck_dkfr_win_o3(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x),coefk(3,19)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
! OZONE
      common /band13dkfr_o3/hk13o3(1), c13o3(1,3,19)
      common /band12dkfr_o3/hk12o3(5), c12o3(5,3,19)
      common /band11dkfr_o3/hk11o3(2), c11o3(2,3,19)
      integer j1, j2
 
      if (ib .eq. 13) then
         ni = 1
         rk1 = 1.0
         rk = 1.0
         hk = hk13o3(1)
         do jp = 1, 19
            do j1 = 1, 3
               coefk(j1,jp) = c13o3(1,j1,1)
            end do
         end do
      else if (ib .eq. 12) then
         ni = 5
         rk1 = 9.172445
         rk = 9.0
         hk = hk12o3(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c12o3(ig,j1,j2)
            end do
         end do
      else if (ib .eq. 11) then
         ni = 2
         rk1 = 11.0838052
         rk = 36.0
         hk = hk11o3(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c11o3(ig,j1,j2)
            end do
         end do
      else
         stop ' bad band ib not 11,12,13'
      endif
 
      if (ig .gt. ni) stop ' BAD ig'
 
      do m = 1, nv
! ! O3 UNITS
cbm modified original Fu code where the mass of the air in the layer 
cbm was calculated via the hydrostatic equation; replaced by the 
cbm average density which is consistently calculated in uvspec         
cbm      u(m) = 1.02*(po(m)+po(m+1))*0.5*(pp(m+1)-pp(m))
         u(m) = (po(m)+po(m+1))*0.5*( pa(i)*0.1 )
      end do
 
      call stp_correction (ig, rk1, rk, coefk, u, tg)
 
      return 
      end 
 
c *********************************************************************
      subroutine ck_dkfr_win_ch4(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x), coefk(3,19), dz
      common /thick/ dz(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
!CH4
      common /band11dkfr_ch4/hk11ch4(4), c11ch4(4,3,19)
      integer j1, j2
 
      if (ib .eq. 11) then
 
         ni = 4
         rk1 = 4.289433
         rk = 8.0
         hk = hk11ch4(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c11ch4(ig,j1,j2)
            end do
         end do
      else if (ib.eq.12 .or. ib.eq.13) then
         do j1 = 1, 100
            tg(j1) = 0
         end do
         return 
 
      else
         stop ' bad band ib not 11,12,13'
      endif
 
      if (ig .gt. ni) stop ' BAD ig'
 
      do m = 1, nv
cbm         dz = gpt(pp(m+1),pp(m),pt(m+1),pt(m),ph(m+1),ph(m))
 
! ! CH4 UNITS
!         uch4(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*1.75e-06*fac(m)*
!      *          1.0e+05/6.023e+23*16.04303
 
cbm         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umch4*1.0E-06*dz*
cbm     1      1.0E+05/6.023E+23*16.04303
         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umch4*1.0E-06*dz(m)*
     1      1.0E+05/6.023E+23*16.04303
      end do
 
      call stp_correction (ig, rk1, rk, coefk, u, tg)
 
      return 
      end 
c *********************************************************************
      subroutine ck_dkfr_win_n2o(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x),coefk(3,19), dz
      common /thick/ dz(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
 
!N2O
      common /band11dkfr_n2o/hk11n2o(2), c11n2o(2,3,19)
      integer j1, j2
 
      if (ib .eq. 11) then
 
         ni = 2
         rk1 = 15.071363
         rk = 23.0
         hk = hk11n2o(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c11n2o(ig,j1,j2)
            end do
         end do
 
      else if (ib.eq.12 .or. ib.eq.13) then
         do j1 = 1, 100
            tg(j1) = 0
         end do
         return 
      else
         stop ' bad band ib not 11,12,13'
      endif
 
      if (ig .gt. ni) stop ' BAD ig'
 
      do m = 1, nv
cbm         dz = gpt(pp(m+1),pp(m),pt(m+1),pt(m),ph(m+1),ph(m))
! ! N2O UNITS
! !      un2o(m)=p(m)*1.01325e+06/(1.3805e-16*t(m))*3.10e-07*fac(m)*
! !     *          1.0e+05/6.023e+23*44.0128
 
cbm         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umn2o*1.0E-06*dz*
cbm     1      1.0E+05/6.023E+23*44.0128
         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umn2o*1.0E-06*dz(m)*
     1        1.0E+05/6.023E+23*44.0128
 
      end do
 
      call stp_correction (ig, rk1, rk, coefk, u, tg)
 
      return 
      end 
 
c *********************************************************************
      subroutine ck_dkfr_win_co2(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), u(nv1x), coefk(3,19), dz
      common /thick/ dz(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /umcon/umco2, umch4, umn2o
 
!CO2
      common /band13dkfr_co2/hk13co2(2), c13co2(2,3,19)
      common /band12dkfr_co2/hk12co2(2), c12co2(2,3,19)
      integer j1, j2
 
      if (ib .eq. 13) then
 
         ni = 2
         rk1 = 0.0059926849
         rk = 60.
         hk = hk13co2(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c13co2(ig,j1,j2)
            end do
         end do
 
      else if (ib .eq. 12) then
         ni = 2
         rk1 = 0.0130594055
         rk = 44.0
         hk = hk12co2(ig)
         do j2 = 1, 19
            do j1 = 1, 3
               coefk(j1,j2) = c12co2(ig,j1,j2)
            end do
         end do
      else if (ib .eq. 11) then
         do j1 = 1, 100
            tg(j1) = 0
         end do
         return 
      else
         stop ' bad band ib not 11,12,13'
      endif
 
      if (ig .gt. ni) stop ' BAD ig'
 
      do m = 1, nv
cbm         dz = gpt(pp(m+1),pp(m),pt(m+1),pt(m),ph(m+1),ph(m))
! ! CO2 UNITS
!       u(m)=pp(m)*1.01325e+06/(1.3805e-16*pt(m))*3.50e-04*
!     *          1.0e+05/6.023e+23*44.00995
 
cbm         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umco2*1.0E-06*dz*
cbm     1      1.0E+05/6.023E+23*44.00995
         u(m) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*umco2*1.0E-06*dz(m)*
     1      1.0E+05/6.023E+23*44.00995
      end do
 
      call stp_correction (ig, rk1, rk, coefk, u, tg)
 
      return 
      end 
c---------------------------------------------------------------
c===============================================================
      subroutine stp_correction(ig, rk1, rk, coefk, ug, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), ug(nv1x),coefk(3,19)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
 
      real k(nvx),fkg(nv1x) 
! in the original code is written k(100), but 
! k(nvx) seems more logical, as tg and k run from
! m=1, nv
 
      real stp(19)
      data stp/ 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 3.98, 6.31, 10.0
     1   , 15.8, 25.1, 39.8, 63.1, 100.0, 158.0, 251.0, 398.0, 631.0, 
     2   1000.0/ 
 
c - - - - - - - - - - - - - - - - - - -
      k(1) = rk1
      do iii = 2, ig
         k(iii) = rk*k(iii-1)
!     print*,iii,k1,k(iii)
      end do
 
      do m = 1, nv
!     ug(m) = 1.02 *(ph(m)+ph(m+1))*0.5 * (pp(m+1) - pp(m))
         ml = 1
 
         if (pp(m) .lt. stp(1)) then
            x1 = coefk(1,1) + coefk(2,1)*(pt(m)-250.0) + coefk(3,1)*(pt(
     1         m)-250.0)**2
            fkg(m) = x1*pp(m)/stp(1)
         else if (pp(m) .gt. stp(19)) then
            x1 = coefk(1,18) + coefk(2,18)*(pt(m)-250.0) + coefk(3,18)*(
     1         pt(m)-250.0)**2
            x2 = coefk(1,19) + coefk(2,19)*(pt(m)-250.0) + coefk(3,19)*(
     1         pt(m)-250.0)**2
            fkg(m) = x1 + (x2 - x1)/(stp(19)-stp(18))*(pp(m)-stp(18))
         else
    2       continue
            if (pp(m) .ge. stp(ml)) then
               ml = ml + 1
               go to 2
            endif
            x1 = coefk(1,ml-1) + coefk(2,ml-1)*(pt(m)-250.0) + coefk(3,
     1         ml-1)*(pt(m)-250.0)**2
            x2 = coefk(1,ml) + coefk(2,ml)*(pt(m)-250.0) + coefk(3,ml)*(
     1         pt(m)-250.0)**2
            fkg(m)=x1+(x2-x1)/(stp(ml)-stp(ml-1))*(pp(m)-stp(ml-1))
         endif
 
         tg(m) = k(ig)*ug(m)*fkg(m)
 
      end do
 
!        do m=1,nv
!        tg(m)=k(ig)*ug(m)*fkg(m)
!      print'(3i3,3f10.7)',ib,ig,m,k(ig),u(m),fkg(m)
!        end do
      return 
      end 
c---------------------------------------------------------------
cbm      real function gpt (p1, p2, t1, t2, q1, q2)
cbm 
cbm      qbar = (q1 + q2)*0.5
cbm      tbar = (t1 + t2)*0.5
cbm      avt = tbar*(1.0 + 0.61*qbar)
cbm      h = .001*29.3*avt
cbm      gpt = h*log(p1/p2)
cbm      return 
cbm      end 
cbm 
c===============================================================
c *********************************************************************
      subroutine ck_dkfr_win_cfc(ib, ig, hk, tg)
      INCLUDE 'fl_radparams.inc'
      real pp, pt, pa, ph, po, tg(nvx), dz
      common /thick/ dz(nvx)
      common /atmos/pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)
      common /cfcs/cfc_conc(3)
      real cfcc(3,3,11:13), corc(3,11:13), rat(3,11:13), cfcmw(3)
      real uf(3), fk(3)
      data cfcmw/ 137.3685, 120.9139, 86.4689/ 
      data cfcc/ 37.940, 9.9010E-03, 4.2008E-05, 1822.228, 2.1404E-04, 
     1   4.4274E-06, 2698.829, -1.8124E-04, 6.7193E-07, 891.794, 
     2   1.9494E-04, 7.0829E-06, 800.771, 1.3996E-03, 1.7760E-06, 
     3   417.009, 2.8908E-03, 7.1388E-06, 5611.777, 1.2661E-03, 
     4   3.5594E-06, 2250.491, 8.7737E-04, 5.8844E-06, 1295.441, 
     5   9.6513E-04, 1.3128E-05/ 
      data corc/ 1.110020, 1.091854, 1.111199, 0.891109, 0.869315, 
     1   0.902223, 0.277778, 0.666667, 0.555556/ 
      data rat/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.277778, 0.666667, 
     1   0.555556/ 
 
 
      if (ib.lt.11 .or. ib.gt.13) stop 
 
      hk = 1
      if (ig .gt. 1) stop ' BAD ig'
 
      do m = 1, nv
cbm         dz = gpt(pp(m+1),pp(m),pt(m+1),pt(m),ph(m+1),ph(m))
 
         do ifc = 1, 3
cbm            uf(ifc) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*cfc_conc(ifc)*dz*
cbm     1         1.0E+05/6.023E+23*cfcmw(ifc)
            uf(ifc) = pp(m)*1.0E+03/(1.3805E-16*pt(m))*cfc_conc(ifc)
     1         *dz(m)*1.0E+05/6.023E+23*cfcmw(ifc)
 
            fk(ifc) = corc(ifc,ib)*rat(ifc,ib)*cfcc(1,ifc,ib)*(1.0000 + 
     1         cfcc(2,ifc,ib)*(pt(m)-250.0)+cfcc(3,ifc,ib)*(pt(m)-250.0)
     2         **2)
 
         end do
 
         tg(m) = fk(1)*uf(1) + fk(2)*uf(2) + fk(3)*uf(3)
 
      end do
 
      return 
      end 
 
 
 
c---------------------------------------------------------------
      BLOCK DATA DKFR13
      COMMON /BAND13DKFR/HK13H2O(5), C13H2O(5,3,19)
! OZONE
      COMMON /BAND13DKFR_O3/HK13O3(1), C13O3(1,3,19)
!CO2
      COMMON /BAND13DKFR_CO2/HK13CO2(2), C13CO2(2,3,19)
      DATA HK13H2O/ 0.874946, 0.028572, 0.067003, 0.021822, 0.007657/ 
      DATA ((C13H2O(1,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00000, 0.00000, 0.00000, 0.00013, 0.00105, 0.00509
     2   , 0.01786, 0.04702, 0.09765, 0.18006, 0.31177, 0.51197, 0.81269
     3   , 1.25906, 1.90230, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 0.000E+00, 0.000E+00, 0.000E+00, 5.375E-06, 4.287E-05, 
     5   2.076E-04, 6.810E-04, 1.616E-03, 3.120E-03, 5.532E-03, 
     6   9.357E-03, 1.515E-02, 2.361E-02, 3.602E-02, 5.343E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     8   0.000E+00, 0.000E+00, 6.812E-08, 5.406E-07, 2.611E-06, 
     9   8.206E-06, 1.806E-05, 3.347E-05, 5.794E-05, 9.710E-05, 
     .   1.560E-04, 2.385E-04, 3.588E-04, 5.217E-04/ 
      DATA ((C13H2O(2,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00000, 0.00000, 0.00011, 0.00129, 0.00602, 0.02093
     2   , 0.05523, 0.09989, 0.15814, 0.24711, 0.37581, 0.56387, 0.83740
     3   , 1.20265, 1.67662, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 0.000E+00, 0.000E+00, 4.513E-06, 5.184E-05, 2.485E-04, 
     5   7.961E-04, 1.801E-03, 3.121E-03, 4.914E-03, 7.687E-03, 
     6   1.156E-02, 1.708E-02, 2.491E-02, 3.588E-02, 4.880E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     8   0.000E+00, 5.406E-08, 6.447E-07, 3.169E-06, 9.547E-06, 
     9   1.889E-05, 3.294E-05, 5.286E-05, 8.252E-05, 1.228E-04, 
     .   1.793E-04, 2.563E-04, 3.768E-04, 5.087E-04/ 
      DATA ((C13H2O(3,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00014, 0.00094, 0.00394, 0.01208, 0.02774, 0.05169
     2   , 0.08407, 0.12949, 0.19531, 0.29073, 0.42682, 0.61650, 0.86541
     3   , 1.17320, 1.52430, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 5.388E-06, 3.800E-05, 1.634E-04, 4.536E-04, 9.603E-04, 
     5   1.688E-03, 2.688E-03, 4.092E-03, 6.089E-03, 8.962E-03, 
     6   1.313E-02, 1.890E-02, 2.647E-02, 3.556E-02, 4.414E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.156E-08, 
     8   4.769E-07, 2.089E-06, 5.409E-06, 1.083E-05, 1.848E-05, 
     9   2.926E-05, 4.432E-05, 6.505E-05, 9.451E-05, 1.386E-04, 
     .   1.993E-04, 2.806E-04, 3.767E-04, 4.511E-04/ 
      DATA ((C13H2O(4,JT,JP),JP=1,19),JT=1,3)/ 0.00120, 0.00128, 0.00242
     1   , 0.00553, 0.01117, 0.01998, 0.03244, 0.04899, 0.07138, 0.10317
     2   , 0.15001, 0.21618, 0.30731, 0.42906, 0.57909, 0.75300, 0.92892
     3   , 1.08421, 1.18233, 1.316E-04, 1.392E-04, 1.634E-04, 2.889E-04
     4   , 5.016E-04, 7.824E-04, 1.160E-03, 1.660E-03, 2.346E-03, 
     5   3.319E-03, 4.714E-03, 6.728E-03, 9.492E-03, 1.318E-02, 
     6   1.765E-02, 2.217E-02, 2.616E-02, 2.976E-02, 3.187E-02, 
     7   2.578E-06, 2.721E-06, 2.733E-06, 4.265E-06, 6.791E-06, 
     8   9.797E-06, 1.375E-05, 1.903E-05, 2.647E-05, 3.687E-05, 
     9   5.102E-05, 7.200E-05, 1.007E-04, 1.398E-04, 1.872E-04, 
     .   2.286E-04, 2.576E-04, 2.814E-04, 2.908E-04/ 
      DATA ((C13H2O(5,JT,JP),JP=1,19),JT=1,3)/ 1.77726, 1.77579, 1.77489
     1   , 1.77406, 1.77231, 1.76948, 1.76735, 1.76208, 1.75493, 1.74598
     2   , 1.72806, 1.69572, 1.63730, 1.54934, 1.42998, 1.28142, 1.10596
     3   , 0.91165, 0.71913, 5.213E-02, 5.209E-02, 5.203E-02, 5.197E-02
     4   , 5.187E-02, 5.172E-02, 5.157E-02, 5.126E-02, 5.082E-02, 
     5   5.023E-02, 4.936E-02, 4.807E-02, 4.614E-02, 4.337E-02, 
     6   3.979E-02, 3.564E-02, 3.083E-02, 2.559E-02, 2.038E-02, 
     7   5.478E-04, 5.475E-04, 5.464E-04, 5.454E-04, 5.440E-04, 
     8   5.420E-04, 5.393E-04, 5.352E-04, 5.289E-04, 5.206E-04, 
     9   5.093E-04, 4.933E-04, 4.725E-04, 4.419E-04, 4.043E-04, 
     .   3.635E-04, 3.168E-04, 2.693E-04, 2.206E-04/ 
 
! OZONE
!      COMMON /BAND13DKFR_O3/HK13O3(1), C13O3(1,3,19)
      DATA HK13O3/ 1.0/ 
      DATA (C13O3(1,JT,1),JT=1,3)/ 1.246705, 2.6904E-02, 2.7135E-04/ 
 
!CO2
!      COMMON /BAND13DKFR_CO2/HK13CO2(2), C13CO2(2,3,19)
      DATA HK13CO2/ 0.972025, 0.027975/ 
      DATA ((C13CO2(1,JT,JP),JP=1,19),JT=1,3)/ 0.00045, 0.00099, 0.00210
     1   , 0.00438, 0.00860, 0.01503, 0.02430, 0.03752, 0.05664, 0.08348
     2   , 0.12152, 0.17345, 0.24203, 0.33027, 0.45103, 0.61785, 0.85883
     3   , 1.15553, 1.47206, 3.141E-05, 6.461E-05, 1.278E-04, 2.514E-04
     4   , 4.454E-04, 7.285E-04, 1.132E-03, 1.704E-03, 2.511E-03, 
     5   3.616E-03, 5.135E-03, 7.211E-03, 9.914E-03, 1.338E-02, 
     6   1.774E-02, 2.311E-02, 3.138E-02, 4.135E-02, 5.177E-02, 
     7   5.347E-07, 1.064E-06, 2.035E-06, 3.884E-06, 6.516E-06, 
     8   1.028E-05, 1.561E-05, 2.321E-05, 3.374E-05, 4.788E-05, 
     9   6.663E-05, 9.232E-05, 1.255E-04, 1.686E-04, 2.194E-04, 
     .   2.739E-04, 3.639E-04, 4.687E-04, 5.743E-04/ 
      DATA ((C13CO2(2,JT,JP),JP=1,19),JT=1,3)/ 1.63918, 1.63805, 1.63609
     1   , 1.63448, 1.62868, 1.62197, 1.61184, 1.59748, 1.57952, 1.55626
     2   , 1.52820, 1.49593, 1.45372, 1.40012, 1.32752, 1.21144, 1.08034
     3   , 0.91050, 0.72074, 5.840E-02, 5.836E-02, 5.828E-02, 5.820E-02
     4   , 5.797E-02, 5.771E-02, 5.731E-02, 5.681E-02, 5.608E-02, 
     5   5.522E-02, 5.410E-02, 5.275E-02, 5.111E-02, 4.912E-02, 
     6   4.651E-02, 4.222E-02, 3.826E-02, 3.243E-02, 2.560E-02, 
     7   6.567E-04, 6.563E-04, 6.552E-04, 6.540E-04, 6.515E-04, 
     8   6.483E-04, 6.436E-04, 6.383E-04, 6.294E-04, 6.197E-04, 
     9   6.057E-04, 5.882E-04, 5.685E-04, 5.452E-04, 5.148E-04, 
     .   4.644E-04, 4.276E-04, 3.641E-04, 2.865E-04/ 
 
      END 
 
c===================================================================
      BLOCK DATA DKFR12
      COMMON /BAND12DKFR/HK12H2O(4), C12H2O(4,3,19)
! OZONE
      COMMON /BAND12DKFR_O3/HK12O3(5), C12O3(5,3,19)
!CO2
      COMMON /BAND12DKFR_CO2/HK12CO2(2), C12CO2(2,3,19)
      DATA HK12H2O/ 0.914825, 0.045491, 0.036386, 0.003298/ 
      DATA ((C12H2O(1,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00000, 0.00000, 0.00009, 0.00078, 0.00452, 0.01483
     2   , 0.03426, 0.06639, 0.11794, 0.19817, 0.32318, 0.51457, 0.80714
     3   , 1.23416, 1.80137, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 0.000E+00, 0.000E+00, 3.087E-06, 2.678E-05, 1.455E-04, 
     5   4.626E-04, 1.070E-03, 2.083E-03, 3.708E-03, 6.252E-03, 
     6   1.015E-02, 1.610E-02, 2.492E-02, 3.792E-02, 5.561E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     8   0.000E+00, 3.219E-08, 2.975E-07, 1.519E-06, 4.734E-06, 
     9   1.090E-05, 2.135E-05, 3.804E-05, 6.464E-05, 1.046E-04, 
     .   1.658E-04, 2.524E-04, 3.817E-04, 5.685E-04/ 
      DATA ((C12H2O(2,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00005, 0.00057, 0.00374, 0.01246, 0.02747, 0.04804
     2   , 0.07634, 0.11725, 0.17851, 0.26999, 0.40956, 0.60833, 0.85850
     3   , 1.14913, 1.39638, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 1.700E-06, 1.995E-05, 1.177E-04, 3.789E-04, 8.418E-04, 
     5   1.517E-03, 2.429E-03, 3.740E-03, 5.666E-03, 8.523E-03, 
     6   1.261E-02, 1.852E-02, 2.637E-02, 3.494E-02, 4.316E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.688E-08, 
     8   2.206E-07, 1.194E-06, 3.778E-06, 8.425E-06, 1.556E-05, 
     9   2.510E-05, 3.895E-05, 5.897E-05, 8.877E-05, 1.274E-04, 
     .   1.836E-04, 2.658E-04, 3.456E-04, 4.354E-04/ 
      DATA ((C12H2O(3,JT,JP),JP=1,19),JT=1,3)/ 0.02572, 0.02678, 0.02958
     1   , 0.03599, 0.04565, 0.06011, 0.08145, 0.11230, 0.15518, 0.21712
     2   , 0.30136, 0.40518, 0.52890, 0.66382, 0.79336, 0.90048, 0.97465
     3   , 1.01230, 0.98289, 1.362E-03, 1.393E-03, 1.450E-03, 1.638E-03
     4   , 1.924E-03, 2.335E-03, 2.941E-03, 3.848E-03, 5.123E-03, 
     5   6.871E-03, 9.322E-03, 1.241E-02, 1.614E-02, 2.032E-02, 
     6   2.429E-02, 2.793E-02, 3.071E-02, 3.224E-02, 3.241E-02, 
     7   2.040E-05, 2.072E-05, 2.110E-05, 2.297E-05, 2.579E-05, 
     8   2.946E-05, 3.465E-05, 4.291E-05, 5.505E-05, 7.050E-05, 
     9   9.342E-05, 1.237E-04, 1.601E-04, 2.030E-04, 2.433E-04, 
     .   2.842E-04, 3.160E-04, 3.324E-04, 3.490E-04/ 
      DATA ((C12H2O(4,JT,JP),JP=1,19),JT=1,3)/ 2.71173, 2.70906, 2.69995
     1   , 2.68591, 2.66765, 2.64146, 2.59903, 2.53778, 2.45614, 2.34626
     2   , 2.21089, 2.05945, 1.88700, 1.69715, 1.50251, 1.30286, 1.10293
     3   , 0.90151, 0.71955, 8.568E-02, 8.552E-02, 8.527E-02, 8.498E-02
     4   , 8.437E-02, 8.368E-02, 8.244E-02, 8.084E-02, 7.844E-02, 
     5   7.541E-02, 7.146E-02, 6.686E-02, 6.154E-02, 5.576E-02, 
     6   4.970E-02, 4.319E-02, 3.651E-02, 3.006E-02, 2.363E-02, 
     7   8.699E-04, 8.665E-04, 8.643E-04, 8.640E-04, 8.569E-04, 
     8   8.513E-04, 8.400E-04, 8.293E-04, 8.071E-04, 7.826E-04, 
     9   7.463E-04, 7.006E-04, 6.477E-04, 5.909E-04, 5.275E-04, 
     .   4.576E-04, 3.863E-04, 3.228E-04, 2.508E-04/ 
 
! OZONE
!      COMMON /BAND12DKFR_O3/HK12O3(5), C12O3(5,3,19)
      DATA HK12O3/ 0.303817, 0.420361, 0.200442, 0.059613, 0.015767/ 
      DATA ((C12O3(1,JT,JP),JP=1,19),JT=1,3)/ 0.00979, 0.02336, 0.04790
     1   , 0.09069, 0.16051, 0.27301, 0.43742, 0.67745, 1.04000, 1.45151
     2   , 2.05584, 2.83048, 3.73868, 4.69546, 5.56884, 6.26397, 6.73334
     3   , 7.07776, 7.46791, 2.488E-04, 4.436E-04, 7.564E-04, 1.272E-03
     4   , 2.068E-03, 3.314E-03, 5.033E-03, 7.398E-03, 1.078E-02, 
     5   1.542E-02, 2.158E-02, 2.948E-02, 3.868E-02, 4.834E-02, 
     6   5.633E-02, 6.131E-02, 6.424E-02, 6.715E-02, 6.959E-02, 
     7   2.181E-06, 2.691E-06, 3.966E-06, 6.009E-06, 8.725E-06, 
     8   1.186E-05, 1.672E-05, 2.224E-05, 8.728E-06, 4.355E-05, 
     9   5.922E-05, 7.806E-05, 9.862E-05, 1.250E-04, 1.489E-04, 
     .   1.626E-04, 1.971E-04, 2.556E-04, 2.880E-04/ 
      DATA ((C12O3(2,JT,JP),JP=1,19),JT=1,3)/ 0.06620, 0.08860, 0.12406
     1   , 0.17583, 0.25003, 0.35240, 0.50090, 0.71273, 1.04000, 1.39262
     2   , 1.89454, 2.52711, 3.33143, 4.34682, 5.58255, 6.97494, 8.39649
     3   , 9.69907, 10.73610, 1.341E-03, 1.543E-03, 1.816E-03, 2.133E-03
     4   , 2.548E-03, 3.105E-03, 3.869E-03, 4.897E-03, 6.356E-03, 
     5   8.401E-03, 1.094E-02, 1.395E-02, 1.701E-02, 2.004E-02, 
     6   2.315E-02, 2.693E-02, 3.068E-02, 3.452E-02, 3.773E-02, 
     7   9.797E-06, 9.625E-06, 9.444E-06, 9.709E-06, 9.672E-06, 
     8   9.862E-06, 9.825E-06, 9.816E-06, -1.221E-05, 1.395E-05, 
     9   1.855E-05, 2.357E-05, 2.492E-05, 2.418E-05, 1.727E-05, 
     .   7.997E-06, -6.288E-06, -1.890E-05,  - 4.795E-05/ 
      DATA ((C12O3(3,JT,JP),JP=1,19),JT=1,3)/ 0.25965, 0.27977, 0.30863
     1   , 0.35051, 0.40986, 0.49103, 0.60846, 0.77464, 1.04000, 1.32456
     2   , 1.75906, 2.30693, 2.92641, 3.54370, 4.06031, 4.43662, 4.65689
     3   , 4.74587, 4.77463, 3.488E-03, 3.439E-03, 3.397E-03, 3.410E-03
     4   , 3.492E-03, 3.645E-03, 3.881E-03, 4.122E-03, 4.333E-03, 
     5   4.246E-03, 3.837E-03, 3.322E-03, 2.773E-03, 2.215E-03, 
     6   1.372E-03, 2.910E-04, -8.093E-04, -1.853E-03, -2.843E-03, 
     7   1.886E-05, 1.845E-05, 1.774E-05, 1.646E-05, 1.524E-05, 
     8   1.421E-05, 1.372E-05, 1.383E-05, -8.053E-06, 1.462E-05, 
     9   1.117E-05, 6.956E-06, 3.922E-06, -3.063E-07, 1.575E-06, 
     .   1.906E-06, 4.144E-06, 5.428E-06, 1.528E-05/ 
      DATA ((C12O3(4,JT,JP),JP=1,19),JT=1,3)/ 0.56031, 0.56833, 0.58032
     1   , 0.59986, 0.63098, 0.67687, 0.75268, 0.86332, 1.04000, 1.16430
     2   , 1.29948, 1.37607, 1.37588, 1.30907, 1.19139, 1.05733, 0.93554
     3   , 0.85314, 0.79321, 4.727E-03, 4.675E-03, 4.568E-03, 4.404E-03
     4   , 4.150E-03, 3.796E-03, 3.303E-03, 2.660E-03, 1.944E-03, 
     5   1.260E-03, 7.381E-04, 2.411E-04, -1.656E-04, -4.431E-04, 
     6   -5.450E-04, -5.321E-04, -5.521E-04, -5.740E-04, -5.761E-04, 
     7   1.656E-06, 1.663E-06, 2.081E-06, 2.925E-06, 3.434E-06, 
     8   3.553E-06, 3.359E-06, 3.138E-06, -2.039E-05, -3.812E-07, 
     9   6.906E-07, 6.906E-07, 1.609E-06, 5.031E-07, 1.387E-06, 
     .   1.516E-06, 1.859E-06, 1.206E-06, 4.594E-07/ 
      DATA ((C12O3(5,JT,JP),JP=1,19),JT=1,3)/ 1.37137, 1.36422, 1.35091
     1   , 1.33042, 1.29804, 1.24585, 1.18554, 1.10850, 1.04000, 0.87896
     2   , 0.73774, 0.59572, 0.46651, 0.35956, 0.27922, 0.22430, 0.18639
     3   , 0.15853, 0.13974, -2.193E-03, -2.136E-03, -2.075E-03, 
     4   -1.948E-03, -1.765E-03, -1.479E-03, -1.287E-03, -1.134E-03, 
     5   -9.195E-04, -7.087E-04, -5.296E-04, -3.786E-04, -2.871E-04, 
     6   -2.111E-04, -1.944E-04, -1.777E-04, -1.585E-04, -1.308E-04, 
     7   -1.138E-04, 5.094E-07, -9.377E-08, -2.844E-07, -6.031E-07, 
     8   -1.406E-06, -2.656E-06, -2.278E-06, -1.478E-06, -2.396E-05, 
     9   -1.088E-06, -2.844E-07, 1.656E-07, 1.719E-07, 4.156E-07, 
     .   4.844E-07, 2.375E-07, 3.562E-07, 4.250E-07, 4.250E-07/ 
 
!CO2
!      COMMON /BAND12DKFR_CO2/HK12CO2(2), C12CO2(2,3,19)
      DATA HK12CO2/ 0.961324, 0.038676/ 
      DATA ((C12CO2(1,JT,JP),JP=1,19),JT=1,3)/ 0.00012, 0.00043, 0.00111
     1   , 0.00256, 0.00577, 0.01151, 0.01911, 0.03057, 0.04735, 0.07152
     2   , 0.10639, 0.15660, 0.22750, 0.32466, 0.45192, 0.62695, 0.86151
     3   , 1.15048, 1.45610, 8.013E-06, 2.554E-05, 6.049E-05, 1.357E-04
     4   , 2.728E-04, 4.795E-04, 7.884E-04, 1.237E-03, 1.897E-03, 
     5   2.829E-03, 4.110E-03, 5.950E-03, 8.587E-03, 1.225E-02, 
     6   1.695E-02, 2.269E-02, 3.052E-02, 3.991E-02, 5.045E-02, 
     7   1.334E-07, 4.066E-07, 9.128E-07, 2.027E-06, 3.781E-06, 
     8   6.025E-06, 9.991E-06, 1.554E-05, 2.375E-05, 3.518E-05, 
     9   5.013E-05, 7.142E-05, 1.023E-04, 1.458E-04, 2.020E-04, 
     .   2.619E-04, 3.453E-04, 4.411E-04, 5.584E-04/ 
      DATA ((C12CO2(2,JT,JP),JP=1,19),JT=1,3)/ 1.46378, 1.46528, 1.46749
     1   , 1.47035, 1.47393, 1.48012, 1.48753, 1.49545, 1.50304, 1.50490
     2   , 1.49706, 1.47418, 1.43550, 1.38056, 1.31013, 1.20658, 1.07377
     3   , 0.91729, 0.74149, 5.042E-02, 5.045E-02, 5.049E-02, 5.053E-02
     4   , 5.064E-02, 5.071E-02, 5.086E-02, 5.107E-02, 5.123E-02, 
     5   5.125E-02, 5.093E-02, 5.017E-02, 4.878E-02, 4.666E-02, 
     6   4.406E-02, 4.049E-02, 3.623E-02, 3.330E-02, 2.516E-02, 
     7   5.531E-04, 5.531E-04, 5.533E-04, 5.533E-04, 5.547E-04, 
     8   5.538E-04, 5.542E-04, 5.563E-04, 5.563E-04, 5.561E-04, 
     9   5.517E-04, 5.439E-04, 5.280E-04, 5.023E-04, 4.711E-04, 
     .   4.318E-04, 3.889E-04, 3.905E-04, 2.717E-04/ 
      END 
 
c===================================================================
      BLOCK DATA DKFR11
      COMMON /BAND11DKFR/HK11H2O(5), C11H2O(5,3,19)
! OZONE
      COMMON /BAND11DKFR_O3/HK11O3(2), C11O3(2,3,19)
! CH4
      COMMON /BAND11DKFR_CH4/HK11CH4(4), C11CH4(4,3,19)
!N2O
      COMMON /BAND11DKFR_N2O/HK11N2O(2), C11N2O(2,3,19)
      DATA HK11H2O/ 0.444198, 0.385220, 0.126262, 0.037311, 0.007009/ 
      DATA ((C11H2O(1,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00008, 0.00086
     2   , 0.00363, 0.01052, 0.03001, 0.07842, 0.18661, 0.41396, 0.76228
     3   , 1.30036, 2.06106, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.562E-06, 
     5   1.941E-05, 8.462E-05, 2.456E-04, 6.580E-04, 1.644E-03, 
     6   3.748E-03, 7.797E-03, 1.456E-02, 2.503E-02, 3.915E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     8   0.000E+00, 0.000E+00, 0.000E+00, 3.437E-09, 1.378E-07, 
     9   6.531E-07, 2.066E-06, 5.325E-06, 1.293E-05, 2.824E-05, 
     .   4.874E-05, 9.226E-05, 1.576E-04, 2.403E-04/ 
      DATA ((C11H2O(2,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00000
     1   , 0.00000, 0.00000, 0.00000, 0.00041, 0.00157, 0.00488, 0.01337
     2   , 0.03167, 0.06774, 0.12891, 0.22321, 0.36117, 0.55392, 0.82533
     3   , 1.21923, 1.76842, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00
     4   , 0.000E+00, 0.000E+00, 9.488E-06, 3.591E-05, 1.075E-04, 
     5   2.814E-04, 6.704E-04, 1.390E-03, 2.508E-03, 4.255E-03, 
     6   6.841E-03, 1.047E-02, 1.544E-02, 2.242E-02, 3.177E-02, 
     7   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 
     8   0.000E+00, 7.219E-08, 2.891E-07, 8.313E-07, 2.047E-06, 
     9   5.128E-06, 1.039E-05, 1.720E-05, 2.786E-05, 4.424E-05, 
     .   6.838E-05, 9.892E-05, 1.392E-04, 1.899E-04/ 
      DATA ((C11H2O(3,JT,JP),JP=1,19),JT=1,3)/ 0.00000, 0.00000, 0.00009
     1   , 0.00042, 0.00127, 0.00332, 0.00786, 0.01660, 0.03295, 0.05885
     2   , 0.09495, 0.14408, 0.21230, 0.30777, 0.43984, 0.62226, 0.86117
     3   , 1.16200, 1.51320, 0.000E+00, 0.000E+00, 2.700E-06, 1.196E-05
     4   , 3.451E-05, 8.425E-05, 1.880E-04, 3.823E-04, 7.209E-04, 
     5   1.236E-03, 1.947E-03, 2.900E-03, 4.223E-03, 6.060E-03, 
     6   8.547E-03, 1.193E-02, 1.638E-02, 2.238E-02, 2.939E-02, 
     7   0.000E+00, 0.000E+00, 2.813E-08, 1.153E-07, 3.172E-07, 
     8   7.562E-07, 1.613E-06, 3.231E-06, 5.766E-06, 9.300E-06, 
     9   1.415E-05, 2.055E-05, 2.922E-05, 4.091E-05, 5.612E-05, 
     .   7.652E-05, 1.039E-04, 1.446E-04, 1.910E-04/ 
      DATA ((C11H2O(4,JT,JP),JP=1,19),JT=1,3)/ 0.01451, 0.01494, 0.01790
     1   , 0.02372, 0.03162, 0.04273, 0.05752, 0.07786, 0.10566, 0.14347
     2   , 0.19622, 0.27002, 0.36739, 0.48830, 0.63343, 0.78459, 0.92817
     3   , 1.07507, 1.17461, 5.337E-04, 5.444E-04, 5.690E-04, 6.828E-04
     4   , 8.466E-04, 1.056E-03, 1.346E-03, 1.729E-03, 2.245E-03, 
     5   2.966E-03, 3.963E-03, 5.316E-03, 7.084E-03, 9.453E-03, 
     6   1.236E-02, 1.552E-02, 1.820E-02, 2.039E-02, 2.194E-02, 
     7   6.587E-06, 6.703E-06, 6.375E-06, 6.969E-06, 8.166E-06, 
     8   9.406E-06, 1.145E-05, 1.382E-05, 1.680E-05, 2.156E-05, 
     9   2.793E-05, 3.610E-05, 4.620E-05, 6.149E-05, 8.040E-05, 
     .   1.038E-04, 1.208E-04, 1.250E-04, 1.247E-04/ 
      DATA ((C11H2O(5,JT,JP),JP=1,19),JT=1,3)/ 1.82788, 1.82672, 1.82619
     1   , 1.82564, 1.82672, 1.82556, 1.82564, 1.82332, 1.81516, 1.80014
     2   , 1.76685, 1.70987, 1.63467, 1.53211, 1.40742, 1.26146, 1.09063
     3   , 0.90847, 0.71754, 3.468E-02, 3.466E-02, 3.461E-02, 3.460E-02
     4   , 3.453E-02, 3.442E-02, 3.431E-02, 3.414E-02, 3.385E-02, 
     5   3.342E-02, 3.274E-02, 3.170E-02, 3.024E-02, 2.829E-02, 
     6   2.584E-02, 2.293E-02, 1.979E-02, 1.677E-02, 1.326E-02, 
     7   2.034E-04, 2.034E-04, 2.029E-04, 2.033E-04, 2.011E-04, 
     8   2.002E-04, 1.985E-04, 1.964E-04, 1.938E-04, 1.885E-04, 
     9   1.834E-04, 1.780E-04, 1.678E-04, 1.567E-04, 1.400E-04, 
     .   1.198E-04, 1.020E-04, 8.899E-05, 7.262E-05/ 
 
 
! OZONE
!      COMMON /BAND11DKFR_O3/HK11O3(2), C11O3(2,3,19)
      DATA HK11O3/ 0.973683, 0.026317/ 
      DATA ((C11O3(1,JT,JP),JP=1,19),JT=1,3)/ 0.11160, 0.11694, 0.12485
     1   , 0.13708, 0.15490, 0.18035, 0.21504, 0.26357, 0.33159, 0.42643
     2   , 0.54808, 0.69060, 0.84504, 1.00000, 1.14124, 1.26234, 1.36787
     3   , 1.45446, 1.52682, 1.565E-03, 1.554E-03, 1.549E-03, 1.550E-03
     4   , 1.558E-03, 1.566E-03, 1.591E-03, 1.643E-03, 1.749E-03, 
     5   1.881E-03, 1.968E-03, 2.034E-03, 2.118E-03, 2.236E-03, 
     6   2.412E-03, 2.523E-03, 2.606E-03, 2.693E-03, 2.761E-03, 
     7   8.206E-06, 7.972E-06, 8.069E-06, 7.944E-06, 7.925E-06, 
     8   7.778E-06, 7.797E-06, 7.200E-06, 6.806E-06, 7.622E-06, 
     9   7.987E-06, 8.516E-06, 8.937E-06, 8.841E-06, 9.491E-06, 
     .   1.026E-05, 1.036E-05, 1.059E-05, 1.026E-05/ 
      DATA ((C11O3(2,JT,JP),JP=1,19),JT=1,3)/ 1.96365, 1.95785, 1.94512
     1   , 1.92795, 1.89909, 1.85230, 1.77961, 1.69148, 1.60713, 1.53430
     2   , 1.43732, 1.30709, 1.15729, 1.00000, 0.85407, 0.72868, 0.61960
     3   , 0.52682, 0.45076, 8.716E-04, 8.880E-04, 9.111E-04, 9.630E-04
     4   , 1.046E-03, 1.102E-03, 1.269E-03, 1.354E-03, 1.358E-03, 
     5   1.239E-03, 1.052E-03, 9.302E-04, 8.200E-04, 6.619E-04, 
     6   4.876E-04, 3.455E-04, 2.603E-04, 1.853E-04, 9.575E-05, 
     7   8.147E-06, 6.481E-06, 7.809E-06, 7.037E-06, 6.797E-06, 
     8   5.803E-06, 6.141E-06, 5.719E-06, 5.050E-06, 5.925E-06, 
     9   4.594E-06, 5.062E-06, 4.119E-06, 3.866E-06, 3.291E-06, 
     .   2.494E-06, 1.806E-06, 1.712E-06, 1.619E-06/ 
 
! CH4
!      COMMON /BAND11DKFR_CH4/HK11CH4(4), C11CH4(4,3,19)
      DATA HK11CH4/ 0.865020, 0.060018, 0.061265, 0.013697/ 
      DATA ((C11CH4(1,JT,JP),JP=1,19),JT=1,3)/ 0.00175, 0.00186, 0.00236
     1   , 0.00368, 0.00585, 0.00947, 0.01530, 0.02524, 0.04176, 0.06855
     2   , 0.10814, 0.16470, 0.24366, 0.35005, 0.48602, 0.66037, 0.88063
     3   , 1.15237, 1.48382, 5.966E-05, 6.169E-05, 5.980E-05, 7.925E-05
     4   , 1.093E-04, 1.556E-04, 2.256E-04, 3.337E-04, 4.966E-04, 
     5   7.509E-04, 1.138E-03, 1.722E-03, 2.543E-03, 3.719E-03, 
     6   5.353E-03, 7.563E-03, 1.047E-02, 1.427E-02, 1.895E-02, 
     7   6.503E-07, 6.559E-07, 5.644E-07, 6.750E-07, 8.688E-07, 
     8   1.097E-06, 1.459E-06, 1.906E-06, 2.516E-06, 3.159E-06, 
     9   4.234E-06, 6.203E-06, 9.141E-06, 1.401E-05, 2.062E-05, 
     .   2.791E-05, 3.566E-05, 4.864E-05, 6.190E-05/ 
      DATA ((C11CH4(2,JT,JP),JP=1,19),JT=1,3)/ 0.01050, 0.01070, 0.01170
     1   , 0.01435, 0.01842, 0.02516, 0.03542, 0.05029, 0.06951, 0.09404
     2   , 0.12764, 0.17138, 0.22980, 0.31581, 0.44617, 0.62437, 0.86205
     3   , 1.18590, 1.61479, 2.913E-04, 2.931E-04, 2.859E-04, 3.178E-04
     4   , 3.660E-04, 4.407E-04, 5.510E-04, 6.879E-04, 8.481E-04, 
     5   1.047E-03, 1.328E-03, 1.753E-03, 2.388E-03, 3.420E-03, 
     6   4.829E-03, 6.539E-03, 8.871E-03, 1.205E-02, 1.613E-02, 
     7   2.619E-06, 2.622E-06, 2.428E-06, 2.456E-06, 2.669E-06, 
     8   2.844E-06, 3.250E-06, 3.597E-06, 3.722E-06, 4.097E-06, 
     9   4.534E-06, 6.256E-06, 9.581E-06, 1.288E-05, 1.313E-05, 
     .   1.214E-05, 1.304E-05, 1.223E-05, 1.058E-05/ 
      DATA ((C11CH4(3,JT,JP),JP=1,19),JT=1,3)/ 0.04628, 0.04671, 0.04856
     1   , 0.05280, 0.05843, 0.06590, 0.07557, 0.08820, 0.10576, 0.13018
     2   , 0.16622, 0.21913, 0.29579, 0.40519, 0.55081, 0.72897, 0.91824
     3   , 1.09785, 1.24313, 7.931E-04, 7.956E-04, 7.733E-04, 8.039E-04
     4   , 8.394E-04, 8.905E-04, 9.610E-04, 1.052E-03, 1.195E-03, 
     5   1.413E-03, 1.751E-03, 2.267E-03, 3.048E-03, 4.213E-03, 
     6   5.740E-03, 7.539E-03, 9.501E-03, 1.118E-02, 1.237E-02, 
     7   5.384E-06, 5.391E-06, 5.194E-06, 5.266E-06, 5.403E-06, 
     8   5.488E-06, 5.744E-06, 5.806E-06, 5.884E-06, 6.200E-06, 
     9   6.344E-06, 6.313E-06, 6.097E-06, 6.766E-06, 8.466E-06, 
     .   8.053E-06, 1.122E-05, 1.293E-05, 1.237E-05/ 
      DATA ((C11CH4(4,JT,JP),JP=1,19),JT=1,3)/ 1.78344, 1.77799, 1.76799
     1   , 1.75348, 1.73023, 1.69947, 1.65701, 1.60388, 1.54595, 1.49093
     2   , 1.45025, 1.42143, 1.39841, 1.36918, 1.31127, 1.20896, 1.07866
     3   , 0.93319, 0.78806, 1.787E-02, 1.777E-02, 1.764E-02, 1.745E-02
     4   , 1.714E-02, 1.671E-02, 1.616E-02, 1.548E-02, 1.481E-02, 
     5   1.426E-02, 1.392E-02, 1.389E-02, 1.394E-02, 1.386E-02, 
     6   1.330E-02, 1.220E-02, 1.076E-02, 9.265E-03, 7.711E-03, 
     7   1.669E-05, 1.618E-05, 1.572E-05, 1.508E-05, 1.539E-05, 
     8   1.365E-05, 1.323E-05, 1.234E-05, 1.162E-05, 1.214E-05, 
     9   1.118E-05, 1.287E-05, 1.518E-05, 1.561E-05, 1.248E-05, 
     .   1.177E-05, 8.472E-06, 6.391E-06, 2.587E-06/ 

!N2O
!      COMMON /BAND11DKFR_N2O/HK11N2O(2), C11N2O(2,3,19)
      DATA HK11N2O/ 0.930764, 0.069236/ 
      DATA ((C11N2O(1,JT,JP),JP=1,19),JT=1,3)/ 0.00050, 0.00114, 0.00233
     1   , 0.00476, 0.00908, 0.01587, 0.02595, 0.04097, 0.06286, 0.09246
     2   , 0.13543, 0.19313, 0.27407, 0.38342, 0.51938, 0.68028, 0.88705
     3   , 1.14064, 1.41256, 1.893E-05, 2.733E-05, 4.675E-05, 7.575E-05
     4   , 1.148E-04, 1.660E-04, 2.398E-04, 3.465E-04, 5.000E-04, 
     5   7.157E-04, 9.908E-04, 1.342E-03, 1.798E-03, 2.432E-03, 
     6   3.199E-03, 3.929E-03, 4.523E-03, 5.232E-03, 5.988E-03, 
     7   2.587E-07, 2.581E-07, 4.000E-07, 4.562E-07, 5.750E-07, 
     8   7.063E-07, 9.250E-07, 1.194E-06, 1.513E-06, 2.750E-06, 
     9   3.450E-06, 4.894E-06, 6.013E-06, 7.222E-06, 8.538E-06, 
     .   1.260E-05, 1.383E-05, 1.305E-05, 1.275E-05/ 
      DATA ((C11N2O(2,JT,JP),JP=1,19),JT=1,3)/ 1.65497, 1.65363, 1.65134
     1   , 1.64740, 1.64078, 1.63262, 1.61814, 1.58872, 1.57649, 1.54641
     2   , 1.51493, 1.47966, 1.43172, 1.36766, 1.29030, 1.19705, 1.07703
     3   , 0.93204, 0.77538, 6.922E-03, 6.915E-03, 6.918E-03, 6.918E-03
     4   , 6.950E-03, 6.960E-03, 6.978E-03, 7.006E-03, 6.963E-03, 
     5   6.893E-03, 6.772E-03, 6.544E-03, 6.250E-03, 5.907E-03, 
     6   5.453E-03, 5.034E-03, 4.693E-03, 4.293E-03, 3.837E-03, 
     7   3.273E-05, 3.249E-05, 3.229E-05, 3.205E-05, 3.246E-05, 
     8   3.092E-05, 3.154E-05, 3.730E-05, 2.915E-05, 2.854E-05, 
     9   2.840E-05, 2.623E-05, 2.578E-05, 2.549E-05, 2.447E-05, 
     .   2.185E-05, 2.198E-05, 2.152E-05, 2.233E-05/ 
      END 
