!============================================================
      SUBROUTINE TRMM_WINDOW_TRMM_WNFLT(COSVZA)
      INCLUDE 'fl_trmm_window.inc'
      REAL TRMM_WINDOW_TRMM_FULIOU_FILTER
!fu_sr      = FU-liou Radiance in each window band
!fu_sf      = FU-liou Flux in each window band
 
! sat_flt_sr = emmulation of TRMM FILTERED RADIANCE in each window band
! sat_unf_sr = emmulation of TRMM UNFILTERED RADIANCE in each window band
 
! sat_flt_r = emmulation of TRMM FILTERED RADIANCE
! sat_unf_r = emmulation of TRMM UNFILTERED RADIANCE
! sat_f     = emmulation of TRMM Window Flux
 
      SAT_FLT_R = 0.0
      SAT_UNF_R = 0.0
      SAT_F = 0.0
C *** Do construct BAND ***
      DO IB = 11, 13
 
         ADM(IB) = FU_SF(IB)/(3.14159*FU_SR(IB))
         FF = TRMM_WINDOW_TRMM_FULIOU_FILTER(IB,0,IWNCLD(IB),FU_SR(IB),
     1      COSVZA)
         UFF = TRMM_WINDOW_TRMM_FULIOU_FILTER(IB,1,IWNCLD(IB),FU_SR(IB),
     1      COSVZA)
 
         SAT_FLT_SR(IB) = FU_SR(IB)*FF
         SAT_UNF_SR(IB) = SAT_FLT_SR(IB)*UFF
         SAT_FLT_R = SAT_FLT_R + SAT_FLT_SR(IB)
         SAT_UNF_R = SAT_UNF_R + SAT_UNF_SR(IB)
! Flux method 1 ( SPECTRAL ADM  )
         SAT_F = SAT_F + SAT_UNF_SR(IB)*(3.14159*ADM(IB))
 
 
!  FLUX METHOD 2 ( DIFFUSE ANGLE Aprox)
!     ff  = TRMM_FULIOU_FILTER(ib,0,iwncld(ib),fu_sr(ib),0.6)
!     uff = TRMM_FULIOU_FILTER(ib,1,iwncld(ib),fu_sr(ib),0.6)
!     sat_flt_sf(ib) =  fu_sf(ib) * ff
!     sat_unf_sf(ib) = sat_flt_sf(ib) * uff
!     sat_f     = sat_f     + sat_unf_sf(ib)
! !
 
      END DO
 
!     print* , 'SCCWN',sum( sat_unf_sr ) / sum( sat_flt_sr )
 
      RETURN 
      END 
!============================================================
      REAL FUNCTION TRMM_WINDOW_TRMM_FULIOU_FILTER (IB, IDIR, ICLD, 
     1   UNFILTERED_RADIANCE, COSVZA)
      PARAMETER (IBW = 1)
      INCLUDE 'fl_trmm_window.inc'
      REAL COEFS(4,3,0:1,0:1,2)
      INTEGER IBT(11:13)
      DATA IBT/ 3, 2, 1/ 
      DATA COEFS/ 5.137E-01, 5.002E-04, -4.887E-03, 3.724E-05, 7.156E-01
     1   , -2.305E-04, 9.184E-04, 1.272E-05, 6.067E-01, 6.400E-04, 
     2   -9.066E-03, 3.633E-05, 1.382E+00, -8.756E-05, 7.637E-04, 
     3   -3.077E-06, 1.398E+00, 4.567E-04, -1.826E-03, -2.524E-05, 
     4   1.463E+00, -3.787E-04, 1.142E-03, -4.245E-05, 4.946E-01, 
     5   3.123E-03, -2.532E-04, -6.721E-05, 7.187E-01, -1.150E-03, 
     6   2.028E-04, 9.825E-05, 6.101E-01, -2.080E-03, -2.730E-03, 
     7   2.223E-04, 1.384E+00, -3.153E-04, 1.790E-04, 6.385E-06, 
     8   1.391E+00, 2.233E-03, -3.908E-04, -1.896E-04, 1.460E+00, 
     9   8.689E-04, -2.329E-04, -1.286E-04, 5.137E-01, 5.002E-04, 
     .   -4.887E-03, 3.724E-05, 7.156E-01, -2.305E-04, 9.184E-04, 
     1   1.272E-05, 6.067E-01, 6.400E-04, -9.066E-03, 3.633E-05, 
     2   1.346E+00, 3.661E-05, 5.024E-04, -1.379E-06, 1.398E+00, 
     3   4.567E-04, -1.826E-03, -2.524E-05, 1.390E+00, 4.193E-04, 
     4   -8.215E-03, 5.476E-05, 4.946E-01, 3.123E-03, -2.532E-04, 
     5   -6.721E-05, 7.187E-01, -1.150E-03, 2.028E-04, 9.825E-05, 
     6   6.101E-01, -2.080E-03, -2.730E-03, 2.223E-04, 1.344E+00, 
     7   2.939E-04, 1.753E-04, -1.015E-05, 1.391E+00, 2.233E-03, 
     8   -3.908E-04, -1.896E-04, 1.395E+00, -2.624E-03, -2.174E-03, 
     9   2.697E-04/ 
 
      IF (IB.LT.11 .OR. IB.GT.13) THEN
         TRMM_WINDOW_TRMM_FULIOU_FILTER = 1.0000
         RETURN 
      ENDIF
 
      IF((IDIR.EQ.0.OR.IDIR.EQ.1).AND.(ICLD.EQ.0.OR.ICLD.EQ.1))THEN
      ELSE
         STOP ' ERROR idir or icld NE to 0 or 1'
      ENDIF
 
      ICLD_USE = ICLD
      IF (ICLD_USE .EQ. 1) THEN
         IF (IB.EQ.11 .AND. UNFILTERED_RADIANCE.GT.8.0 .OR. IB.EQ.12
     1       .AND. UNFILTERED_RADIANCE.GT.8.0 .OR. IB.EQ.13 .AND. 
     2      UNFILTERED_RADIANCE.GT.20.0) ICLD_USE = 0
      ENDIF
 
      TRMM_WINDOW_TRMM_FULIOU_FILTER = COEFS(1,IBT(IB),IDIR,ICLD_USE,1)
     1    + UNFILTERED_RADIANCE*COEFS(2,IBT(IB),IDIR,ICLD_USE,1) + 
     2   COSVZA*COEFS(3,IBT(IB),IDIR,ICLD_USE,1) + UNFILTERED_RADIANCE**
     3   2*COEFS(4,IBT(IB),IDIR,ICLD_USE,1)
 
      RETURN 
 
      END 
