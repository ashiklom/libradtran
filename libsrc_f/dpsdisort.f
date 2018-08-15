      SUBROUTINE  sdisort( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $     wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, umu, 
     $     nphi, phi, fbeam, beta, nil, umu0, phi0, newgeo, zd, spher, 
     $     radius, fisot, albedo, btemp, ttemp, temis, deltam, planck, 
     $     onlyfl, accur, quiet, ierror, prnt, header,
     $     maxcly, maxulv, maxumu, maxcmu, maxphi, 
     $     rfldir, rfldn, flup, dfdt, uavg, uu, u0u, nsca,
     $     uavgdn, uavgso, uavgup,
     $     nrefrac, ichap, vn,
     $     ndenssza, denssza, denssig, denstab, dtauc_mb )
C
C+---------------------------------------------------------------------+
C------------------    I/O VARIABLE SPECIFICATIONS     -----------------
C+---------------------------------------------------------------------+
C
      INCLUDE 'DISORT.MXD'

      CHARACTER  HEADER*127
      LOGICAL  DELTAM, newgeo, quiet, planck, onlyfl, PRNT(7), 
     $     SPHER, usrang, USRTAU
      INTEGER  ierror(33), MAXCLY, MAXULV, MAXCMU, maxphi, maxumu, 
     $     NLYR, nphi, NSTR, NTAU, numu, ndenssza
      REAL*4     accur, ALBEDO, beta(0:maxcly), BTEMP,
     $     DTAUC( MAXCLY ), FBEAM, FISOT, phi0, 
     $     PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ), RADIUS,
     $     SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     $     WVNMLO, WVNMHI, UMU0, UTAU( MAXULV ), 
     $     U0U(MAXUMU, MAXULV), UU( MAXUMU, MAXULV, MAXPHI ),
     $     ZD( 0:MAXCLY ), VN(0:MAXCLY), DTAUC_mb( MAXCLY )
C     
      REAL*4     RFLDIR( MAXULV ), RFLDN( MAXULV ),
     $     FLUP( MAXULV ),
     $     UAVG( MAXULV ), DFDT( MAXULV ), umu( maxumu ),
     $     UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )
*
      LOGICAL noplnk
      REAL*8     dpaccur, dpALBEDO, dpbeta(0:mxcly), dpBTEMP,
     $     dpDTAUC( MXCLY ), dpFBEAM, dpFISOT, dpphi0, 
     $     dpPHI( MXPHI ), dpPMOM( 0:MXCMU, MXCLY ), dpRADIUS,
     $     dpSSALB( MXCLY ), dpTEMPER( 0:MXCLY ), dpTEMIS, dpTTEMP,
     $     dpWVNMLO, dpWVNMHI, dpUMU0, dpUTAU( MXULV ), 
     $     dpU0U(MXUMU, MXULV), dpUU( MXUMU, MXULV, MXPHI ),
     $     dpZD( 0:MXCLY ), dpVN(0:MXCLY), dpDTAUC_mb( MXCLY )
C     
      REAL*8     dpRFLDIR( MXULV ), dpRFLDN( MXULV ),
     $     dpFLUP( MXULV ),
     $     dpUAVG( MXULV ), dpDFDT( MXULV ), dpumu( mxumu ),
     $     dpUAVGDN( MXULV ), dpUAVGUP( MXULV ), dpUAVGSO( MXULV )
*
      INTEGER nrefrac, nsca, ichap, nil
*
      INTEGER lc, lu, iq, j, iu
*--------------------------------------------------------
*     Variables that will be eventually input to sdisort
*
      REAL*4 denssza(*), denssig(*), denstab(0:maxcly,*)
      REAL*8 dpdenssza(MXSZA), dpdenssig(0:mxcly), 
     $     dpdenstab(mxsza,0:mxcly)
*
      INTEGER brosza
*
*     First check that the arys are large enough
*
      IF( MXSZA.LT.ndenssza ) THEN
         write (0,*) '* Error, too many solar zenith angles in denssza.' 
         write (0,*) '* Increase MXSZA in libsrc_f/DISORT.MXD to at'
         write (0,*) '* least',ndenssza,' and recompile libRadtran!'
         RETURN
      ENDIF
*
* For input convert between REAL*4 (uvspec working environment) 
* and REAL*8 for double precision sdisort
*
      noplnk     = .NOT. planck
      dpaccur    = dble(accur)
      dpalbedo   = dble(albedo)
      dpbtemp    = dble(btemp)
      dpfbeam    = dble(fbeam)
      dpfisot    = dble(fisot)
      dpphi0     = dble(phi0)
      dpradius   = dble(radius)
      dptemis    = dble(temis)
      dpttemp    = dble(ttemp)
      dpwvnmlo   = dble(wvnmlo)
      dpwvnmhi   = dble(wvnmhi)
      dpumu0     = dble(umu0)

* brosza = 0: density is a function of altitude, that is a profile
* brosza = 1: density is a function of altitude and sza, that is a matrix

      brosza = 0 ! This is the default

      IF ( ndenssza .GT. 0 ) THEN
         brosza = 1
         DO j = 1, ndenssza
            dpdenssza(j) = dble(denssza(j))
         ENDDO
         DO lc = 1, maxcly
            dpdenssig(lc) = dble(denssig(lc))
         ENDDO
         DO j = 1, ndenssza
            DO lc = 0, maxcly
               dpdenstab(j,lc) = dble(denstab(lc,j))
            ENDDO
         ENDDO
      ENDIF
      DO j = 1, maxphi
         dpphi(j) = dble(phi(j))
      ENDDO
      DO lc = 1, maxcly
         dpdtauc(lc)       = dble(dtauc(lc))
         dpdtauc_mb(lc)    = dble(dtauc_mb(lc))
         dpssalb(lc)       = dble(ssalb(lc))        
         DO iq = 0, maxcmu
            dppmom(iq,lc) = dble(pmom(iq,lc))
         ENDDO
      ENDDO
      DO lc = 0, maxcly
         dpbeta(lc)     = dble(beta(lc))
         dptemper(lc)   = dble(temper(lc))
         dpzd(lc)       = dble(zd(lc))    
         dpVN(lc)       = dble(vn(lc))
      ENDDO
      DO lu = 1, maxulv
         dputau(lu)     = dble(utau(lu))
      ENDDO
      DO iu = 1, maxumu
         dpumu(iu)      = dble(umu(iu))
      ENDDO

      CALL  dpsdisort
     $     ( ICHAP, nlyr, dpdtauc, dpssalb, dppmom, dptemper, dpwvnmlo,
     $     dpwvnmhi, usrtau, ntau, dputau, nstr, usrang, numu, dpumu, 
     $     nphi, dpphi, dpfbeam, dpbeta, nil, dpumu0, dpphi0, newgeo, 
     $     dpzd, spher, dpradius, dpfisot, 
     $     dpalbedo, dpbtemp, dpttemp, dptemis, deltam, noplnk, 
     $     onlyfl, dpaccur, quiet, ierror, prnt, header,
     $     mxcly, mxulv, mxumu, mxcmu, mxphi, 
     $     dprfldir, dprfldn, dpflup, dpdfdt, dpuavg, dpuu, dpu0u, nsca,
     $     brosza, dpdtauc_mb, dpdenssig, dpdenstab, dpdenssza,
     $     ndenssza, dpVN, nrefrac,
     $     dpuavgdn, dpuavgso, dpuavgup)

*
* For output convert between REAL*4 (uvspec working environment) 
* and REAL*8 for double precision sdisort
*

      DO lu = 1, maxulv
         rfldir(lu)   = real(dprfldir(lu))
         rfldn(lu)    = real(dprfldn(lu))
         flup(lu)     = real(dpflup(lu))
         uavg(lu)     = real(dpuavg(lu))
         uavgso(lu)   = real(dpuavgso(lu))
         uavgup(lu)   = real(dpuavgup(lu))
         uavgdn(lu)   = real(dpuavgdn(lu))
        DO iu = 1, maxumu
            u0u(iu,lu) = real(dpu0u(iu,lu))
            DO j = 1, maxphi
               uu(iu,lu, j) = real(dpuu(iu,lu,j))
caky               write(*,*)"sdisort", iu, lu, j, uu(iu,lu, j)
            ENDDO
         ENDDO         
      ENDDO

      RETURN
      END


      SUBROUTINE  dpsdisort
     $     ( ICHAP, nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $     wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, umu, 
     $     nphi, phi, fbeam, beta, nil, umu0, phi0, newgeo, zd, spher, 
     $     radius, fisot, albedo, btemp, ttemp, temis, deltam, noplnk, 
     $     onlyfl, accur, quiet, ierror, prnt, header,
     $     maxcly, maxulv, maxumu, maxcmu, maxphi, 
     $     rfldir, rfldn, flup, dfdt, uavg, uu, u0u, nsca,
     $     brosza, dtauc_mb, sigbro, DENSBRO, SZA_BRO,
     $     ndenssza, VN, nrefrac,
     $     uavgdn, uavgso, uavgup )
     
*+-----------------------------------------------------------------------+
* Copyright (C) 1992, 94, 95 Arve Kylling
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 1, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* To obtain a copy of the GNU General Public License write to the
* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
* USA.
*
*+-----------------------------------------------------------------------+
***********************************************************************
*
*     Strongly modified version of disort, this version solves the
*     radiative transfer equation in pseudo-spherical or spherical 
*     geometry. The main changes are:
*
*     1) The medium may be taken to be pseudo-spherical (spher=.true.
*        and nil=0) or spherical (spher=.true. and nil.gt.0). 
*
*     2) Only lambertian surface is allowed
*
*     3) The present implementation of sdisort.f does not work for 
*        thermal radiation. This is not a problem for UVspec, but 
*        may be of concern for those who want to use sdisort.f some 
*        other place.
*
*     If run in spherical geometry the code should be compiled and
*     linked in double precision. I have not tried the spherical
*     version in single precision and do not know if it does what it
*     should do. If you wonder how to build a double precision version
*     of this routine, see the file Makefile_double that comes with the
*     uvspec distribution.
*
*    ( See disort.doc and twostr.doc for further documentation )
*
*     Bug reports and other feedback to arve@kaja.gi.alaska.edu
*
************************************************************************
C
C+---------------------------------------------------------------------+
C------------------    I/O VARIABLE SPECIFICATIONS     -----------------
C+---------------------------------------------------------------------+
C
      implicit double precision (A-H,O-Z)
      INTEGER MXCLY, MXULV, MXCMU, MXUMU, MXPHI,
     $        MI, MI9M2, NNLYRI 
      INCLUDE 'DISORT.MXD'
      INTEGER brosza
      CHARACTER  HEADER*127
      LOGICAL  DELTAM, newgeo, quiet, NOPLNK, onlyfl, PRNT(7), 
     $     SPHER, usrang, USRTAU
      INTEGER  ierror(33), MAXCLY, MAXULV, MAXCMU, maxphi, maxumu, 
     $     NLYR, nphi, NSTR, NTAU, numu
      real*8    accur, ALBEDO, beta(0:maxcly), BTEMP,
     $     DTAUC( MAXCLY ), FBEAM, FISOT, phi0, 
     $     PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ), RADIUS,
     $     SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     $     WVNMLO, WVNMHI, UMU0, UTAU( MAXULV ), 
     $     U0U(MAXUMU, MAXULV), UU( MAXUMU, MAXULV, MAXPHI ),
     $     ZD( 0:MAXCLY ), DTAUC_MB(MAXCLY), SIGBRO(0:*), 
     $     DENSBRO(mxsza,0:*), VN(0:MAXCLY)
C     
      real*8     RFLDIR( MAXULV ), 
     $           RFLDN( MAXULV ),  
     $           FLUP( MAXULV ),
     $     UAVG( MAXULV ), DFDT( MAXULV ), umu( maxumu ),
     $     UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )
      integer nil
C
C+---------------------------------------------------------------------+
C      ROUTINES CALLED (IN ORDER):  SLFTST, ZEROAL, schekin, ssetdis,
C                                   sprtinp, slepoly, ssoleig, SETBCO,
C                                   supbeam, SETPCO, supisot, sterpev,
C                                   sterpso, ssetmtx, ssolve0, sfluxes,
*                                   susrint
C+---------------------------------------------------------------------+
C
C  INDEX CONVENTIONS (FOR ALL DO-LOOPS AND ALL VARIABLE DESCRIPTIONS):
*
*     il     :  For iteration over the perturbation source.
C
C  IQ,JQ,KQ  :  FOR COMPUTATIONAL POLAR ANGLES ('QUADRATURE ANGLES')
C
C   IQ/2     :  FOR HALF THE COMPUTATIONAL POLAR ANGLES (JUST THE ONES
C               IN EITHER 0-90 DEGREES, OR 90-180 DEGREES)
C
C     K,L    :  FOR LEGENDRE EXPANSION COEFFICIENTS OR, ALTERNATIVELY,
C               SUBSCRIPTS OF ASSOCIATED LEGENDRE POLYNOMIALS
C
C     LU     :  FOR USER LEVELS
C
C     LC     :  FOR COMPUTATIONAL LAYERS (EACH HAVING A DIFFERENT
C               SINGLE-SCATTER ALBEDO AND/OR PHASE FUNCTION)
C
C    LEV     :  FOR COMPUTATIONAL LEVELS
C
C+---------------------------------------------------------------------+
C               I N T E R N A L    V A R I A B L E S
C
C   AMB(IQ/2,IQ/2)    FIRST MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'ssoleig')
C
C   APB(IQ/2,IQ/2)    SECOND MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'ssoleig')
C
C   ARRAY(IQ,IQ)      SCRATCH MATRIX FOR 'ssoleig', 'supbeam' AND 'supisot'
C                     (SEE EACH SUBROUTINE FOR DEFINITION)
C
C   B()               RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
C                     *ssolve0*;  RETURNS AS SOLUTION VECTOR
C                     VECTOR CAPITAL-L, THE CONSTANTS OF INTEGsratioN
C
C   BDR(IQ/2,0:IQ/2)  BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A
C                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS
C                     REFERS TO A COMPUTATIONAL ANGLE.  SECOND INDEX:
C                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE -UMU0-;
C                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.
C
C   BPLANK            INTENSITY EMITTED FROM BOTTOM BOUNDARY
C
C   CBAND()           MATRIX OF LEFT-HAND SIDE OF THE LINEAR SYSTEM
C                     EQ. SC(5), SCALED BY EQ. SC(12);  IN BANDED
C                     FORM REQUIRED BY LINPACK SOLUTION ROUTINES
C
C   CC(IQ,IQ)         CAPITAL-C-SUB-IJ IN EQ. SS(5)
C
C   CH(LC)            THE CHAPMAN-FACTOR TO CORRECT FOR PSEUDO-
C                     SPHERICAL GEOMETRY IN THE DIRECT BEAM.
C
C   CMU(IQ)           COMPUTATIONAL POLAR ANGLES (GAUSSIAN)
C
C   CWT(IQ)           QUADRATURE WEIGHTS CORRESP. TO -CMU-
C
C   DELM0             KRONECKER DELTA, DELTA-SUB-M0, WHERE 'M' = MAZ
C                     IS THE NUMBER OF THE FOURIER COMPONENT IN THE
C                     AZIMUTH COSINE EXPANSION
C
C   EVAL(IQ)          TEMPORARY STORAGE FOR EIGENVALUES OF EQ. SS(12)
C
C   EVECC(IQ,IQ)      COMPLETE EIGENVECTORS OF SS(7) ON RETURN FROM
C                     *ssoleig* ; STORED PERMANENTLY IN -GC-
C
C   EXPBEA(LC)        TRANSMISSION OF DIRECT BEAM IN DELTA-M OPTICAL
C                     DEPTH COORDINATES
C
C   FLDN(LU)          DIFFUSE DOWN FLUX (DELTA-M SCALED)
C
C   FLDIR(LU)         DIRECT BEAM FLUX (DELTA-M SCALED)
C
C   FLYR(LC)          TRUNCATED FRACTION IN DELTA-M METHOD
C
C   GL(K,LC)          PHASE FUNCTION LEGENDRE POLY. EXPANSION
C                     COEFFICIENTS, CALCULATED FROM 'PMOM' BY
C                     INCLUDING SINGLE-SCATTERING ALBEDO, FACTOR
C                     2K+1, AND (IF DELTAM=TRUE) THE DELTA-M
C                     SCALING
C
C   GC(IQ,IQ,LC)      EIGENVECTORS AT POLAR QUADRATURE ANGLES,
C                     LITTLD-G  IN EQ. SC(1)
C
*   gl(iu,iq,lc)      Eigenvectors interpolated to user polar angles
*                     ( littlD-g in Eqs. SC(3) and S1(8-9), i.e.
*                     capital-G without the capital-L factor )
*
C   IPVT(LC*IQ)       INTEGER VECTOR OF PIVOT INDICES FOR LINPACK
C                     ROUTINES
C
C   KK(IQ,LC)         EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C
C   LAYRU(LU)         COMPUTATIONAL LAYER IN WHICH USER OUTPUT LEVEL
C                     -UTAU(LU)- IS LOCATED
C
C   LL(IQ,LC)         CONSTANTS OF INTEGsratioN CAPITAL-L IN EQ. SC(1),
C                     OBTAINED BY SOLVING SCALED VERSION OF EQ. SC(5)
C
C   LYRCUT            TRUE, RADIATION IS ASSUMED ZERO BELOW LAYER
C                     -NCUT- BECAUSE OF ALMOST COMPLETE ABSORPTION
C
C   NCUT              COMPUTATIONAL LAYER NUMBER IN WHICH ABSORPTION
C                     OPTICAL DEPTH FIRST EXCEEDS -ABSCUT-
C
C   OPRIM(LC)         SINGLE SCATTERING ALBEDO AFTER DELTA-M SCALING
C
C   PASS1             TRUE ON FIRST ENTRY, FALSE THEREAFTER
C
C   PKAG(0:LC)        INTEGRATED PLANCK FUNCTION FOR INTERNAL EMISSION
C                     AT LAYER BOUNDARIES
C
C   PKAGC(LC)         INTEGRATED PLANCK FUNCTION FOR INTERNAL EMISSION
C                     AT LAYER CENTER
C
*   psi0(iq),         Sum just after square bracket in Eq. SD(9), for
*    psi1(iq)         z0 and z1 (thermal source) respectively  (zbs0
*                     and zbs1 for beam source)
*
C   TAUC(0:LC)        CUMULATIVE OPTICAL DEPTH (UN-DELTA-M-SCALED)
C
C   TAUCPR(0:LC)      CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED IF
C                     DELTAM = TRUE, OTHERWISE EQUAL TO -TAUC-)
C
C   TPLANK            INTENSITY EMITTED FROM TOP BOUNDARY
C
C   U0C(IQ,LU)        AZIMUTHALLY-AVERAGED INTENSITY
C
C   UTAUPR(LU)        OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                     COORDINATES;  EQUAL TO  -UTAU(LU)- IF NO DELTA-M
C
C   WK()              SCRATCH ARRAY
C
C   xb0(IQ,LC)        X-SUB-ZERO COEFFICIENT IN EQ. KS(4) FOR PSEUDO-
C                     SPHERICAL BEAM SOURCE KS(15).
C
C   xb1(IQ,LC)        X-SUB-ONE COEFFICIENT IN EQ. KS(4) FOR PSEUDO-
C                     SPHERICAL BEAM SOURCE KS(15).
C
C   xba(LC)           ALFA COEFFICIENT IN EQ. KS(4) FOR PSEUDO-
C                     SPHERICAL BEAM SOURCE KS(15).
C
C   XR0(LC)           X-SUB-ZERO IN EXPANSION OF THERMAL SOURCE FUNC-
C                     TION; SEE EQS. KS(4) (HAS NO (MU) DEPENDENCE)
C
C   XR1(LC)           X-SUB-ONE IN EXPANSION OF THERMAL SOURCE FUNC-
C                     TION; SEE EQ. KS(4) (HAS NO (MU) DEPENDENCE)
C
C   XRA(LC)           ALFA IN EXPONENT IN EXPANSION OF THERMAL SOURCE 
C                     FUNCTION; SEE EQ. KS(4) (HAS NO (MU) DEPENDENCE)
C
C   YLM0(L)           NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                     OF SUBSCRIPT 'L' AT THE BEAM ANGLE (NOT SAVED
C                     AS FUNCTION OF SUPERSCIPT 'M')
C
C   YLMC(L,IQ)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                     OF SUBSCRIPT 'L' AT THE COMPUTATIONAL ANGLES
C                     (NOT SAVED AS FUNCTION OF SUPERSCIPT 'M')
C
*   ylmu(l,iq)        Normalized associated legendre polynomial
*                     of subscript 'l' at the user angles
*                     (NOT saved as function of superscipt 'm')
*
C   Z()               SCRATCH ARRAY USED IN * ssolve0* TO SOLVE A
C                     LINEAR SYSTEM FOR THE CONSTANTS OF INTEGsratioN
C
C   Z0(IQ)            SOLUTION VECTORS Z-SUB-ZERO OF EQ. KS(7-8),
C                     THERMAL RADIATION
C
C   Z1(IQ)            SOLUTION VECTORS Z-SUB-ONE  OF EQ. KS(7-8),
C                     THERMAL RADIATION
C
C   ZA                ALFA COEFFICIENT IN EQ. KS(4), THERMAL SOURCE
*
*   zbau(iu,lc)       zbsa interpolated to user angles 
*
*   zb0u(iu,lc)       zbs0 interpolated to user angles 
*
*   zb1u(iu,lc)       zbs1 interpolated to user angles 
*
*   zpau(iu,lc)       zpsa interpolated to user angles 
*
*   zp0u(iu,lc)       zps0 interpolated to user angles 
*
*   zp1u(iu,lc)       zps1 interpolated to user angles 
C
C   ZBEAMA(LC)        PERMANENT STORAGE FOR THE ALFA COEFFICIENT IN 
C                     EQ. KS(6) FOR PSEUDO-SPHERICAL BEAM SOURCE.
C
C   ZBEAM0(IQ,LC)     PERMANENT STORAGE FOR THE Y-SUB-ZERO COEFFICIENT 
C                     IN EQ. KS(6) FOR PSEUDO-SPHERICAL BEAM SOURCE 
C                     KS(15). OBTAINED BY SOLVING EQS. KS(7-8)
C
C   ZBEAM1(IQ,LC)     PERMANENT STORAGE FOR THE Y-SUB-ONE COEFFICIENT 
C                     IN EQ. KS(6) FOR PSEUDO-SPHERICAL BEAM SOURCE 
C                     KS(15). FOUND BY SOLVING EQS. KS(7-8)
C
C   ZBS0(IQ)          SOLUTION VECTORS Y-SUB-ZERO OF KS(7-8)
C
C   ZBS1(IQ)          SOLUTION VECTORS Y-SUB-ONE OF KS(7-8)
*
*   zbsa(iq)          Alfa coefficient corresponding to zbs0 and zbs1
C
C   ZJ(IQ)            RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN
C                     EQ. STWJ(6D)
C
C   ZPLK0(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE
C                     VECTORS  -Z0-  OBTAINED BY SOLVING  EQ. KS(7-8)
C
C   ZPLK1(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE
C                     VECTORS  -Z1-  OBTAINED BY SOLVING  EQ. KS(7-8)
C
C   ZPLKA(LC)         PERMANENT STORAGE FOR THE THERMAL SOURCE
C                     COEFFICIENT -ZA-; ALFA IN EQ. KS(7-8)
C
C+---------------------------------------------------------------------+
C   LOCAL SYMBOLIC DIMENSIONS:
C
C       MXCLY  = Max no. of computational layers
C       MXULV  = Max no. of output levels
C       MXCMU  = Max no. of computation polar angles
C       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
C+---------------------------------------------------------------------+
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY )

*
      DOUBLE PRECISION d1mach
      real*8 dpsratio
      LOGICAL LYRCUT
      INTEGER layru( mxulv )
      real*8 acmu(mxcmu), azerr, azterm, bdr( mi,0:mi ),
     $     bplank, ch( mxcly ), chtau( 0:2*mxcly+1 ),
     $     cosphi, cmu( mxcmu ), cwt( mxcmu ), delm0, cwt0(mxcmu),
     $     deltamu, epsil,expbea( 0:mxcly ),
     $     flyr( mxcly ), fldn( mxulv ), 
     $     fldir( mxulv ), gl( 0:mxcmu,mxcly ), oprim( mxcly ), 
     $     PHIRAD( MXPHI ), pi, pkag( 0:mxcly ), pkagc(mxcly), rpd, 
     $     tauc( 0:mxcly ), taucpr( 0:mxcly ), tplank,
     $     u0uls(mxumu, mxulv), uum( mxumu, mxulv, 0:mxcmu), 
     $     utaupr( mxulv ),sza_bro(mxsza)
*
*     Variables for the internal source at computational angles
*
      INTEGER ntau_c, numu_c, layru_c( mxulv )
      real*8  utau_c(mxcly+1),
     $     umu_c( 2*mxcmu ),  utaupr_c( mxulv )
*
*     Variables for the internal source at user angles
*
      INTEGER numu_u
      real*8  umu_u( 2*mxcmu )
*
      DATA  nerr / 33 /
*
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     bro density is a function of SZA => no delta-m method,
C                                         so  that dtaucpr=dtauc
      if (brosza.eq.1) then
        deltam=.false.
      endif
      
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      pi = 2.D0* DASIN(1.0D0)
    
      epsil = 10.D0*d1mach(4)

      
      rpd = pi / 180.D0
*
      IF ( PRNT(1) )  WRITE( *,1010 )  HEADER
*
*     Zero some arrays (not strictly necessary, but otherwise unused
*     parts of arrays collect garbage)
*
      DO i = 1, nerr
         ierror(i) = 0
      ENDDO
      CALL  dpszeroit( CH,   MXCLY )
      CALL  dpszeroit( CMU,  MXCMU )
      CALL  dpszeroit( CWT,  MXCMU )
      CALL  dpszeroit( CWT0, MXCMU )
*     
*     Calculate cumulative optical depth and dither SINGLE-scatter
*     albedo to improve numerical behavior of eigenvalue/vector
*     computation
*

    
      TAUC( 0 ) = 0.D0
      CALL  dpszeroit( TAUC(0), MXCLY+1 )
      DO 20  LC = 1, NLYR
*
*aky 27.11.2008 Added factor of 10 because otherwise I got nan on my MacOX
*
         IF( SSALB(LC) .GT. 1.0-epsil )  SSALB(LC) = 1.D0 - 10.*epsil
         TAUC(LC) = TAUC(LC-1) + DTAUC(LC)
 20   CONTINUE
*
*     Check input dimensions and variables
*
     	 
      CALL  dpschekin( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $     wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, 
     $     umu, fbeam, umu0, spher, fisot, albedo, btemp, ttemp, 
     $     temis, noplnk, onlyfl, ierror, quiet,
     $     maxcly, maxulv, maxumu, maxcmu, 
     $     mxcly, mxulv,  mxcmu, mxumu, tauc, zd )
*
      iret = 0
      DO ierr = 1, nerr
         IF ( ierror(ierr) .NE. 0 ) THEN
            iret = 1
            WRITE(*,'(/,A,I4,/)')  "SDISORT REPORTS FATAL ERROR: ",
     $           ierr
         ENDIF
      ENDDO
       	
      IF ( iret .EQ. 1 ) RETURN
*
*     Perform various setup operations
*
      
      CALL  dpssetdis
     $  (ICHAP, acmu, albedo, bdr, beta, bplank,btemp,ch,chtau, 
     $     cmu, cwt, deltam, deltamu, dtauc,EXPbea, fbeam, 
     $     flyr, gl, layru, layru_c, lyrcut, maxumu, maxcmu,
     $     ncut, newgeo, 
     $     nlyr, ntau, ntau_c, nn, nstr, noplnk, numu, onlyfl,
     $     oprim, pkag, pkagc, pmom, radius, spher, ssalb, tauc,
     $     taucpr,temis, temper, tplank, ttemp, utau, utau_c,
     $     utaupr, utaupr_c, umu0, umu, usrtau, usrang,
     $     wvnmlo, wvnmhi, 
     $     zd, numu_c, umu_c, numu_u, umu_u, BROSZA,SIGBRO,DENSBRO,
     $     dtauc_mb, sza_bro, ndenssza, VN, nrefrac)
    
      
      
*
*     Print input information
*     
      IF ( prnt(1) ) THEN
         write(*,'(A)') 'Layer  Chapman'
         write(*,'(A)') 'number function'
         DO lc = 1, nlyr
            write(*,'(I4,4x,F10.5)') lc, ch(lc)
         ENDDO
      ENDIF

      IF ( prnt(1) ) THEN
            IF (nsca .EQ. 1) THEN
               WRITE(*,*) "Running sdisort in single scattering mode"
            ELSE IF (nsca .EQ. 2) THEN
               WRITE(*,*) "Running sdisort in multiple scattering mode"
            ENDIF
            CALL dpsprtinp( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $           wvnmhi, ntau, utau, nstr, numu, umu, fbeam, umu0,
     $           spher,nil, fisot, albedo, btemp, ttemp, temis, deltam,
     $           noplnk, onlyfl, flyr, lyrcut, oprim, tauc, taucpr, 
     $           maxcmu, prnt(7) )
         ENDIF
C     
C ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
C ========  (EQ STWJ 5)
C
      KCONV = 0
      NAZ = NSTR-1
*
*                                            ** AZIMUTH-INDEPENDENT CASE
*
      IF (FBEAM.EQ.0.0 .OR. (1.D0-UMU0).LT.1.D-5 .OR. ONLYFL.OR.
     $     (NUMU.EQ.1.AND.(1.-UMU(1)).LT.1.D-5 ) ) then
          NAZ = 0
      endif
*
      CALL  dpszeroit( UU, MAXUMU*MAXULV*MAXPHI )
    
caky      write(*,*)"naz", naz
      DO MAZ = 0, NAZ

         IF ( MAZ.EQ.0 )  DELM0 = 1.0D0
         IF ( MAZ.GT.0 )  DELM0 = 0.0D0
*
*   
*        SINGLE SCATTERING
*
         do I=1,mxcmu
           cwt0(I)=cwt(I)
	 enddo
*	  
         IF (NSCA.EQ.1) THEN 
            CALL dpszeroit( CWT, MXCMU )
         ENDIF
*      
       
         CALL dpfoucom(cwt0,
     $        accur, acmu, albedo, bdr, beta, bplank, btemp,
     $        ch, chtau, cmu, cwt, delm0, deltam, deltamu, dtauc, 
     $        EXPbea, fbeam, fisot, flyr, gl, layru, 
     $        layru_c, lyrcut, maxcly, maxcmu, maxulv, maxumu, maz, 
     $        ncut, nil, nlyr, noplnk, nn, nstr, ntau, ntau_c, 
     $        numu, numu_c, numu_u, onlyfl, oprim, pi, pkag, pkagc, 
     $        pmom, prnt, radius, spher, ssalb, tauc,  
     $        taucpr, temis, temper, tplank, ttemp, umu0, umu, 
     $        umu_c, umu_u, usrang, utau, utau_c, utaupr, 
     $        utaupr_c, utaupr_u, wvnmhi, wvnmlo, zd,
     $        dfdt, flup, fldn, fldir, rfldir,rfldn, uavg, 
     $        u0u, uum, u0uls, uavgdn, uavgso, uavgup )
        
          
 
      IF( MAZ.EQ.0 ) THEN
C
         DO  140  J = 1, NPHI
            PHIRAD( J ) = RPD * ( PHI(J) - PHI0 )
 140     CONTINUE
C                               ** SAVE AZIMUTHALLY AVERAGED INTENSITIES
         DO 160  LU = 1, NTAU
            DO 160  IU = 1, NUMU
               U0U( IU,LU ) = UUM( IU,LU,0 )
 160     CONTINUE
C
      END IF
C                                ** INCREMENT INTENSITY BY CURRENT
C                                ** AZIMUTHAL COMPONENT (FOURIER
C                                ** COSINE SERIES);  EQ SD(2)
      AZERR = 0.D0
      DO  J = 1, NPHI
         COSPHI = DCOS( MAZ * PHIRAD(J) )
         DO  LU = 1, NTAU
            DO IU = 1, NUMU
               AZTERM = UUM( IU,LU,MAZ ) * COSPHI
               UU( IU,LU,J ) = UU( IU,LU,J ) + AZTERM
caky           write(*,*)"sdisort gaba", maz,j, lu, iu, uu(iu,lu, j), azterm
caky     $             ,UUM( IU,LU,MAZ )
               AZERR = dmax1(dpsratio( DABS(AZTERM), DABS(UU(IU,LU,J))),
     $                        AZERR )
            ENDDO
         ENDDO
      ENDDO
      IF ( AZERR.LE.ACCUR )  KCONV = KCONV + 1
      IF ( KCONV.GE.2 )   then
       GOTO 210
      endif 
C
      ENDDO
*
210   CONTINUE
*     
      
      RETURN
C
1010  FORMAT ( ////, 1X, 120('*'), /, 25X,
     $  'Spherical version of disort',
     $  /, 1X, A, /, 1X, 120('*') )
      END
*
      SUBROUTINE dpfoucom(cwt0, accur, acmu, albedo, bdr, beta, bplank,
     $     btemp, ch, chtau, cmu, cwt, delm0, deltam, deltamu, dtauc, 
     $     EXPbea, fbeam, fisot, flyr, gl, layru, layru_c, 
     $     lyrcut, maxcly, maxcmu, maxulv, maxumu, maz, ncut,
     $     nil, nlyr, noplnk, nn, nstr, ntau, ntau_c, numu, numu_c,
     $     numu_u, onlyfl, oprim, pi, pkag, pkagc, pmom, prnt, radius, 
     $     spher, ssalb, tauc, taucpr,  temis, temper, 
     $     tplank, ttemp, umu0, umu, umu_c, umu_u, usrang, utau, 
     $     utau_c, utaupr, utaupr_c, utaupr_u, wvnmhi, wvnmlo, 
     $     zd, dfdt, flup, fldn, fldir, rfldir, rfldn,
     $     uavg, u0u, uum, u0uls, uavgdn, uavgso, uavgup )
*
*     Solve Eq. 6a (STWJ) for each FOUrier COMponent
*     
      implicit double precision (A-H, O-Z)

C+---------------------------------------------------------------------+
C   LOCAL SYMBOLIC DIMENSIONS:
C
C       MXCLY  = Max no. of computational layers
C       MXULV  = Max no. of output levels
C       MXCMU  = Max no. of computation polar angles
C       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
C+---------------------------------------------------------------------+
      INTEGER MXCLY, MXULV, MXCMU, MXUMU, MXPHI,
     $        MI, MI9M2, NNLYRI 
      INCLUDE 'DISORT.MXD'
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY )
*
      INTEGER ipvt(nnlyri), layru(mxulv)
      LOGICAL lamber, LYRCUT
      LOGICAL  DELTAM, NOPLNK, onlyfl, PRNT(7), SPHER, usrang
      real*8     accur, ALBEDO, beta(0:maxcly), 
     $     bplank, BTEMP, delm0, deltamu, DTAUC( MAXCLY ), FBEAM,
     $     FISOT, pi, PMOM( 0:MAXCMU, MAXCLY ), RADIUS,
     $     SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, tplank, TTEMP,
     $     WVNMLO, WVNMHI, UMU0, UTAU( MAXULV ), 
     $     U0U(MAXUMU, MAXULV), ZD( 0:MAXCLY )
      real*8     RFLDIR( MAXULV ), RFLDN( MAXULV ),
     $     FLUP( MAXULV ), 
     $     UAVG( MAXULV ), DFDT( MAXULV ), umu( maxumu ),
     $     UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )
*
      real*8   ANGCOS(1), acmu(mxcmu), albedos,
     $     amb( mi,mi ), apb( mi,mi ),
     $     array( mxcmu,mxcmu ), b( nnlyri ), bdr( mi,0:mi ), 
     $     cband( mi9m2,nnlyri ), cc( mxcmu,mxcmu ), cerr, cdiff,
     $     ch( mxcly ), chtau(0:*), cmu( mxcmu ), cwt( mxcmu ),
     $     cwt0(mxcmu),  
     $     eval(mi), evecc( mxcmu, mxcmu ), 
     $     EXPbea( 0:mxcly ), dfdts(mxulv), fbeams, flups(mxulv), 
     $     flyr( mxcly ), fldn( mxulv ), fldir( mxulv ), 
     $     gl( 0:mxcmu,mxcly ), gc(mxcmu,mxcmu,mxcly ), 
     $     gu( mxumu, mxcmu, mxcly ), kk( mxcmu, mxcly ), 
     $     ll( mxcmu,mxcly ), oprim( mxcly ), pkag( 0:mxcly ),
     $     pkagc(mxcly), prev_cerr, psi0( mxcmu ), psi1( mxcmu ),
     $     rfldns(mxulv), rfldirs(mxulv), sgn, tauc( 0:mxcly ),
     $     taucpr( 0:mxcly ), tynt, 
     $     u0c( mxcmu,mxulv ), uavgs( mxulv ),
     $     u0uls(mxumu, mxulv), u0utmp(mxumu, mxulv),
     $     uum( mxumu, mxulv,0:mxcmu), utaupr( mxulv ), wk( mxcmu ),
     $     xr0(mxcly ), xr1( mxcly ), xra(mxcly), ylm0( 0:mxcmu ),
     $     ylmc(0:mxcmu,mxcmu ),ylmu(0:mxcmu, mxumu), z( nnlyri ),
     $     xba(mxcly ), xb0( mxcmu,mxcly ), xb1( mxcmu,mxcly ),
     $     zbeama(mxcly ), zbeam0( mxcmu,mxcly ), 
     $     zbeam1( mxcmu,mxcly ), zbs0(mxcmu ), zbsa, 
     $     zbs1( mxcmu ), z0( 2*mxcmu ),
     $     z1(2*mxcmu ), za,zj0(mxcmu ),  zj2( mxcmu ),
     $     zj( mxcmu ),  zju( mxumu), zplk0( mxcmu,mxcly ),
     $     zplk1( mxcmu,mxcly ), zplka(mxcly), zb0u( mxumu, mxulv ),
     $     zb1u( mxumu, mxulv ), zbau( mxumu,mxulv),
     $     zp0u( mxumu, mxulv), zp1u( mxumu, mxulv ),
     $     zpau(mxumu,mxulv)
*
*     Variables for the internal source at computational angles
*
      INTEGER ntau_c, numu_c, layru_c( mxulv ), nil
      real*8 gu_c( 2*mxumu, mxcmu, mxcly ),
     $     parmu_c(mxcmu, mxulv), utau_c(mxcly+1),
     $     umu_c( 2*mxcmu ),  utaupr_c( mxulv ),
     $     ylmu_c(0:mxcmu, 2*mxumu),
     $     zb0u_c( 2*mxumu, mxulv ), zb1u_c( 2*mxumu, mxulv ), 
     $     zbau_c( 2*mxumu,mxulv), uum_c( 2*mxumu, mxulv, 0:mxcmu)
*
*     Variables for the internal source at user angles
*
      INTEGER numu_u
      real*8 gu_u( 2*mxumu, mxcmu, mxcly ), 
     $     parmu_u(mxcmu, mxulv), 
     $     umu_u( 2*mxcmu ), 
     $     xba_u(mxcly ), xb0_u( mxcmu,mxcly ), xb1_u( mxcmu,mxcly ),
     $     ylmu_u(0:mxcmu, 2*mxumu),
     $     zb0u_u( 2*mxumu, mxulv ), zb1u_u( 2*mxumu, mxulv ), 
     $     zbau_u( 2*mxumu,mxulv),
     $     uum_u( 2*mxumu, mxulv, 0:mxcmu),
     $     zj0_u(mxcmu ), zj2_u( mxcmu ), zj_u( mxcmu )
*
*     This are just dummies, never used for anything, but
*     to avoid memory havocking when calling *sterpso*
*
      real*8   zp0u_c( 2*mxumu, mxulv), 
     $    zp1u_c( 2*mxumu, mxulv ), zpau_c( 2*mxumu,mxulv)
*
      DOUBLE PRECISION   AAD( MI,MI ), EVALD( MI ) , EVECCD( MI,MI ),
     $     WKD( MXCMU )
      DOUBLE PRECISION d1mach
*
*     No bidirectional reflectance here, so set lamber .true.
*
      lamber = .true.
      
   
*
*     Zero some stuff just to make sure....
*
      
      DO 10 I = 1, NNLYRI
         IPVT(I) = 0
10    CONTINUE
      CALL dpszeroit( dfdts, mxulv )
      CALL dpszeroit( flups, mxulv )
      CALL dpszeroit( rfldns, mxulv )
      CALL dpszeroit( rfldirs, mxulv )
      CALL dpszeroit( uavgs, mxulv )
      CALL dpszeroit( u0uls, mxumu*mxulv )
      CALL dpszeroit( u0utmp, mxumu*mxulv )
      CALL  dpszeroit( AMB, MI**2 )
      CALL  dpszeroit( APB, MI**2 )
      CALL  dpszeroit( ARRAY, MXCMU**2 )
      CALL  dpszeroit( CC,    MXCMU**2 )
      CALL  dpszeroit( EVAL  ,MI )
      CALL  dpszeroit( EVECC, MXCMU**2 )
      CALL  dpszeroit( gc, mxcmu*MXCMU*MXCLY )
      CALL  dpszeroit( GU, mxumu*MXCMU*MXCLY )
      CALL  dpszeroit( KK,     MXCMU*MXCLY )
      CALL  dpszeroit( LL,     MXCMU*MXCLY )
      CALL  dpszeroit( WK    ,  MXCMU )
      CALL  dpszeroit( xba, MXCLY )
      CALL  dpszeroit( xb0, MXCMU*MXCLY )
      CALL  dpszeroit( xb1, MXCMU*MXCLY )
      CALL  dpszeroit( XR0, MXCLY )
      CALL  dpszeroit( XR1, MXCLY )
      CALL  dpszeroit( XRA, MXCLY )
      CALL  dpszeroit( Z, NNLYRI )
      CALL  dpszeroit( Z0    ,  MXCMU )
      CALL  dpszeroit( Z1    ,  MXCMU )
      ZA = 0
      CALL  dpszeroit( ZBEAMA, MXCLY )
      CALL  dpszeroit( ZBEAM0, MXCMU*MXCLY )
      CALL  dpszeroit( ZBEAM1, MXCMU*MXCLY )
      CALL  dpszeroit( ZBS0, MXCMU )
      CALL  dpszeroit( ZBS1, MXCMU )
      CALL  dpszeroit( ZJ,  MXCMU )
      CALL  dpszeroit( ZPLK0,  MXCMU*MXCLY )
      CALL  dpszeroit( ZPLK1,  MXCMU*MXCLY )
      CALL  dpszeroit( ZPLKA,  MXCLY )
      CALL  dpszeroit( zbau, mxumu*mxulv )
      CALL  dpszeroit( zb0u, mxumu*mxulv )
      CALL  dpszeroit( zb1u, mxumu*mxulv )
*
*     Only zero these the first time, otherwise EXPect a downer
*
      IF ( maz .EQ. 0 ) THEN
         CALL  dpszeroit( YLM0, MXCMU+1 )
         CALL  dpszeroit( YLMC, (MXCMU+1)*MXCMU )
         CALL  dpszeroit( ylmu, (mxcmu+1)*mxumu )
      ENDIF
*
*
*     Get normalized associated legendre polynomials for incident beam
*     angle cosine
*
      IF ( FBEAM .GT. 0.0D0 ) THEN
         ANGCOS(1) = -UMU0
         CALL  dpslepoly( 1, MAZ, MXCMU, NSTR-1, ANGCOS, YLM0 )
      ENDIF
*
*     Get normalized associated legendre polynomials for computational
*     polar angle cosines
*
      IF ( .NOT. onlyfl .AND. usrang ) THEN
         CALL dpslepoly( numu, maz, mxcmu, nstr-1, umu, ylmu )
         CALL dpslepoly( numu_u, maz, mxcmu, nstr-1, umu_u, ylmu_u )
      ENDIF
      CALL dpslepoly( numu_c, maz, mxcmu, nstr-1, umu_c, ylmu_c )
      CALL  dpslepoly( NN,   MAZ, MXCMU, NSTR-1, CMU, YLMC )
      
*
*     Evaluate normalized associated legendre polynomials with negative
*     -cmu- from those with positive -cmu-; Dave/Armstrong eq. (15)
*
      SGN  = - 1.0D0
      DO L = MAZ, NSTR-1
         SGN = - SGN
         DO  IQ = NN+1, NSTR
            YLMC( L,IQ ) = SGN * YLMC( L,IQ-NN )
         ENDDO
      ENDDO
*     
*     
*     
      albedos  = albedo
      fbeams   = fbeam
      prev_cerr = d1mach(2)
     
*
* ===== Begin loop over each term in Eq. 10, Dahlback and Stamnes (1991)
*
      DO 99 il = 0, nil
*
*     For all but the pseudo-spherical approximation (il=0) set
*     fbeam=0.0 and albedo=0.0.
*
         IF ( il .GT. 0 ) THEN
            albedo = 0.0D0
            fbeam  = 0.0D0
         ENDIF
*     
C
C ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============
C
         DO 100  LC = 1, NCUT
*
*     Solve eigenfunction problem in eq. stwj(8b); return eigenvalues
*     and eigenvectors
*
            CALL dpssoleig( amb, apb, array, cmu, cwt, gl(0,lc),mi, maz,
     $           mxcmu, nn, nstr, wk, ylmc, cc, evecc, eval,
     $           kk(1,lc), gc(1,1,lc), aad, wkd, eveccd, evald, umu0 )
*
*     Set coefficients in ks(7)
*
*     For the beam source or for the perturbation source
*
            IF ( FBEAM.GT.0.D0 .OR. il.GT.0 ) THEN
*
               
               IF ( il .EQ. 0 ) THEN
	        
*
*     Set coefficients for the internal source for il=0
*               
              
                  CALL dpsetbco( ch, chtau, cmu, delm0, fbeam,
     $                 gl, lc, maz, mxcmu, nstr, 
     $                 taucpr, xba, xb0, xb1, zj, 
     $                 ylmc, ylm0)

                  
*
*     Get coefficients at cmu+-deltamu
*                  
                  IF ( nil .GT. 0 ) THEN
		    
                     CALL dpintbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, 2*mxumu, nstr, numu_c,
     $                    taucpr, zb0u_c, zb1u_c, xba,
     $                    zju, ylm0, ylmu_c)
                 
                  ENDIF
*     
                  IF ( usrang ) THEN
*     
*     Get coefficients at umu
*     
                     CALL dpintbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, mxumu, nstr, numu,
     $                    taucpr, zb0u, zb1u, xba,
     $                    zju, ylm0, ylmu)
*     
*     Get coefficients at umu+-deltamu
*     
                     CALL dpintbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, 2*mxumu, nstr, numu_u,
     $                    taucpr, zb0u_u, zb1u_u, xba,
     $                    zju, ylm0, ylmu_u)
                     
		    
                  ENDIF
*
               ELSE
	          
*
*     Set coefficients for the internal source for il > 0
*
                 
                  CALL dpsetpar( acmu, beta, lc, mxcmu, nstr, parmu_c,
     $                 radius,taucpr, xba, xb0, xb1, zd, zj, zj0, zj2 )
               
                  
*
*     Do some index gymnastic to satisfy the fantastic representation
*     of angles in cmu.
*
*     ooooh, bananas.....
*
                  DO iq = 1, nn
                     z0(2*iq-1)      = xb0(nstr+1-iq,lc)
                     z0(2*iq)        = xb0(nstr+1-iq,lc)
                     z0(2*nn+2*iq-1) = xb0(iq,lc)
                     z0(2*nn+2*iq)   = xb0(iq,lc)
                     z1(2*iq-1)      = xb1(nstr+1-iq,lc)
                     z1(2*iq)        = xb1(nstr+1-iq,lc)
                     z1(2*nn+2*iq-1) = xb1(iq,lc)
                     z1(2*nn+2*iq)   = xb1(iq,lc)
                  ENDDO
*     
                  DO iu = 1, numu_c
                     zb0u_c(iu,lc) = z0(iu)
                     zb1u_c(iu,lc) = z1(iu)
                  ENDDO
*
                 
                  IF ( usrang) THEN
                     CALL dpsetpar( umu_u, beta, lc, mxcmu, numu_u,
     $                    parmu_u, radius, taucpr, xba_u, xb0_u,
     $                    xb1_u, zd, zj_u, zj0_u, zj2_u )
*     
                     DO iu = 1, numu
                        zb0u(iu,lc) = xb0_u(iu,lc)
                        zb1u(iu,lc) = xb1_u(iu,lc)
                     ENDDO
                  ENDIF
               ENDIF
*
*     Calculate particular solutions of Eqs. ks(10-11) for incident
*     beam source in pseudo-spherical geometry
*

                              
               CALL dpsupbeam( array, cc, cmu, ipvt, mxcmu, nn, nstr, 
     $              wk, xb0(1,lc), xb1(1,lc), xba(lc), 
     $              zbs0, zbs1, zbsa, zbeam0(1,lc), zbeam1(1,lc), 
     $              zbeama(lc) )
                 
*
*     and now a nip-up....
*     I got this from Arne Dahlback, numerical problems may happen
*     (even in double precision) if the perturbation source is not
*     ignored for optically thin layers. The value for tynt is from
*     Arne. (and tynt is norwegian and means thin).
*
               tynt = 0.001D0
               IF ( tauc(lc).LT.tynt .AND. il .GT. 0 ) THEN
                  CALL dpszeroit(zbs0,mxcmu)
                  CALL dpszeroit(zbs1,mxcmu)
                  DO iq = 1, nstr
                     zbeam0(iq,lc) = 0.D0
                     zbeam1(iq,lc) = 0.D0
                  ENDDO
               ENDIF
            ENDIF
C     
          
            IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 ) THEN
	      
*
*     Set coefficients in ks(7) for thermal source
*
               CALL dpsetpco( lc, pkag, pkagc, taucpr, xra, xr0, xr1)
*
*     Calculate particular solutions of Eq. ks(11) for thermal emission
*     source
*
               CALL dpsupisot( array, cc, cmu, ipvt, mxcmu, nn, nstr,
     $              oprim(lc), wk, xr0(lc), xr1(lc), xra(lc), z0,
     $              z1, za, zplk0(1,lc), zplk1(1,lc), zplka(lc))
            END IF
*
*     Interpolate eigenvectors and source terms to angles needed to
*     get the internal source, cf. Eq. 11c, DS.
*
            IF ( nil .GT. 0 ) THEN
               CALL dpsterpev( cwt, evecc, gl(0,lc), gu_c(1,1,lc), maz,
     $              mxcmu, 2*mxumu, nn, nstr, numu_c, wk, ylmc, ylmu_c)
            ENDIF
*
            IF (il .GT. 0 )  THEN
               fbeam = 1.0 D0                ! Just to get beam-source 
            ENDIF                          ! interpolated for il .gt. 0
*                                          
*     
*     Standard disort way of interpolating source terms
*
*
*     First, interpolate to cmu+-deltamu to be able to calculate
*            derivative term
*     
            IF ( nil .GT. 0 ) THEN
	     
               CALL  dpsterpso( cwt, fbeam, gl(0,lc), maz,
     $              mxcmu, noplnk, numu_c, nstr, oprim(lc),
     $              ylmc, ylmu_c, psi0, psi1, 
     $              xr0(lc), xr1(lc), xra(lc), z0, z1, 
     $              zbs0, zbs1, zbsa,
     $              zb0u_c(1,lc), zb1u_c(1,lc), zbau_c(1,lc), 
     $              zp0u(1,lc), zp1u(1,lc), zpau(1,lc))   
            ENDIF
*
            IF ( usrang ) THEN
*
*     Second, interpolate to umu+-deltamu to be able to calculate
*            derivative term for user angles
*     
               CALL  dpsterpso( cwt, fbeam, gl(0,lc), maz,
     $              mxcmu, noplnk, numu_u, nstr, oprim(lc),
     $              ylmc, ylmu_u, psi0, psi1, 
     $              xr0(lc), xr1(lc), xra(lc), z0, z1,
     $              zbs0, zbs1, zbsa,
     $              zb0u_u(1,lc), zb1u_u(1,lc), zbau_u(1,lc), 
     $              zp0u(1,lc), zp1u(1,lc), zpau(1,lc))   
            ENDIF
*
            IF ( il .GT. 0 ) THEN
               fbeam = 0.0D0      ! Back to zero as it should be here
            ENDIF
*     
            IF ( .NOT.onlyfl .AND. usrang ) THEN
*
*     Interpolate eigenvectors to user angles
*
               CALL  dpsterpev( CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZ,
     $              MXCMU,MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU)
               CALL dpsterpev( cwt, evecc, gl(0,lc), gu_u(1,1,lc), maz,
     $              mxcmu, 2*mxumu, nn, nstr, numu_u, wk, ylmc, ylmu_u)

*
*     Interpolate source terms to user angles
*
               IF (il .gt. 0 )  THEN
                  fbeam = 1.0D0   ! Just to get beam-source interpolated
               ENDIF
*
*     Standard disort way of interpolating source terms
*     
               CALL  dpsterpso( cwt, fbeam, gl(0,lc), maz,
     $              mxcmu, noplnk, numu, nstr, oprim(lc),
     $              ylmc, ylmu, psi0, psi1, 
     $              xr0(lc), xr1(lc), xra(lc), z0, z1,
     $              zbs0, zbs1, zbsa,
     $              zb0u(1,lc), zb1u(1,lc), zbau(1,lc), 
     $              zp0u(1,lc), zp1u(1,lc), zpau(1,lc))   
*     
               IF (il .gt. 0 ) THEN
                  fbeam = 0.D0   ! Back to zero as it should be here
               ENDIF
*     
            END IF
C
C
 100     CONTINUE
*
* ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============
*
*     Set coefficient matrix of equations combining boundary and layer
*     interface conditions
*
         IF ( il .GT. 0 ) CALL dpszeroit(bdr,mi*(mi+1))
         CALL  dpssetmtx( bdr, cband, cmu, cwt, delm0, gc, kk, 
     $        lyrcut, mi, mi9m2, mxcmu, ncol, ncut, 
     $        nnlyri, nn, nstr, taucpr, wk )
*     
*     Solve for constants of integration in homo-geneous solution 
*     (general boundary conditions)
*
         CALL  dpssolve0( albedo, b, bplank, cband, cmu, cwt, 
     $         EXPbea, fbeam, fisot, il, ipvt, ll, lyrcut,
     $         maz, mi, mi9m2, mxcmu, ncol, ncut, nn, noplnk, 
     $         nstr,
     $         nnlyri, pi, tplank, taucpr, umu0, z, zbeam0,
     $         zbeam1, zbeama, zplk0, zplk1, zplka )
*
*     Compute upward and downward sfluxes
*
         IF ( maz .eq. 0 ) THEN
            CALL dpsfluxes( ch, cmu, cwt0, fbeam, gc, il, kk, 
     $           layru, ll, lyrcut, mxcmu, mxulv, ncut, nn, 
     $           nstr, ntau, pi, prnt, ssalb, taucpr, umu0,
     $           utau, utaupr, xr0, xr1, xra, zplk0,
     $           zplk1, zplka, zbeama, zbeam0, zbeam1, dfdt, 
     $           flup, fldn, fldir, rfldir, rfldn, 
     $           uavg, u0c, maxulv, uavgdn, uavgso, uavgup)
*     
         ENDIF
*
*     Compute intensity components at
*     angles needed to calculate the derivative term in Eq. 11c
*     at computational angles.
*
     
         IF ( nil .GT. 0 ) THEN

            CALL  dpsusrint( albedo, bplank, cmu, cwt, delm0, EXPbea,
     $           fbeam, fisot, gc, gu_c, il, kk, lamber, layru_c, ll,
     $           lyrcut, maz, mxcmu, mxulv, 2*mxumu, ncut,
     $           nlyr, nn, nstr, noplnk, numu_c, ntau_c, pi, 
     $           taucpr, tplank, umu_c, umu0, utaupr_c, wk,
     $           zb0u_c, zb1u_c, zbau_c, zbeam0, zbeam1, zbeama,
     $           zp0u_c, zp1u_c, zpau_c, zplk0, zplk1, zplka, uum_c )
         ENDIF
*     
         IF ( usrang ) THEN
*
*     Compute azimuthal averaged intensity components at
*     user angles.
*          
            CALL  dpsusrint( albedo, bplank, cmu, cwt, delm0, EXPbea,
     $           fbeam, fisot, gc, gu, il, kk, lamber, layru, ll,
     $           lyrcut, maz, mxcmu, mxulv, mxumu, ncut,
     $           nlyr, nn, nstr, noplnk, numu, ntau, pi, 
     $           taucpr, tplank, umu, umu0, utaupr, wk,
     $           zb0u, zb1u, zbau, zbeam0, zbeam1, zbeama,
     $           zp0u, zp1u, zpau, zplk0, zplk1, zplka, uum )
*
           
            DO lu = 1, ntau
               DO iu = 1, numu
                  u0utmp( iu,lu ) = u0utmp( iu,lu ) + uum( iu,lu, maz )
               ENDDO
            ENDDO
         ELSE 
*                                     ** Compute azimuthal averaged 
*                                     ** intensity components at
*                                     ** quadrature angles 
            DO lu = 1, ntau
               DO 130 iu = 1, nstr
                  u0utmp( iu,lu ) = u0utmp( iu,lu ) + uum( iu,lu, maz )
 130           CONTINUE
            ENDDO
         END IF
*
*     Check for convergence in Eq. 10, Dahlback and Stamnes (1991).
*      
         cerr  = 0.D0
         DO lu = 1, ntau
            IF ( il .GT. 0 ) THEN
               DO iu = 1, numu
                  IF (DABS(uum(iu,lu,maz)) .GT. 1.0D-09 ) THEN
                     cdiff = DABS( uum(iu,lu,maz) / u0utmp(iu,lu))
                  ELSE
                     cdiff = 0.D0
                  ENDIF
                  cerr = dmax1( cdiff, cerr )
               ENDDO
            ENDIF
         ENDDO
*
*     If cerr is not increasing save results from iteration
*
         IF ( cerr .LT. prev_cerr ) THEN
            IF ( maz .eq. 0 ) THEN
               DO lu = 1, ntau
                  dfdts(lu)    =  dfdts(lu)   + dfdt(lu)
                  flups(lu)    =  flups(lu)   + flup(lu)
                  rfldns(lu)   =  rfldns(lu)  + rfldn(lu)
                  rfldirs(lu)  =  rfldirs(lu) + rfldir(lu)
                  uavgs(lu)    =  uavgs(lu)   + uavg(lu)
               ENDDO
            ENDIF
            IF ( usrang ) THEN
               DO lu = 1, ntau
                  DO iu = 1, numu
                     u0uls( iu,lu ) = u0uls( iu,lu ) + uum( iu,lu, maz )
                  ENDDO
               ENDDO
            ELSE 
               DO lu = 1, ntau
                  DO  iu = 1, nstr
                     u0uls( iu,lu ) = u0uls( iu,lu ) + uum( iu,lu, maz )
                  ENDDO
               ENDDO
            END IF            
         ENDIF
         IF ( ( ( cerr .LT. accur ) .OR.
     $        ( cerr .GT. prev_cerr ) )
     $        .AND. il .GT. 0 ) GOTO 98
         IF ( il .GT. 0 ) prev_cerr = cerr
*     
*     Calculate the partial derivative of the intensity with respect to
*     the cosine of the polar angle for angles needed to calculate the
*     derivative term in Eq. 11c at computational angles.
*     
         IF ( nil .GT. 0 ) THEN
            CALL dpPARDER( deltamu, maz, mxcmu, mxulv, 2*mxumu, ntau_c,
     $           nstr, parmu_c, uum_c )
         ENDIF
*
         IF ( usrang ) THEN
            CALL dpPARDER( deltamu, maz, mxcmu, mxulv, 2*mxumu, ntau_c,
     $           numu_u/2, parmu_u, uum_u )
         ENDIF
*
 99   CONTINUE
*
 98   CONTINUE
*      
      DO lu = 1, ntau
       DO iu = 1, numu
        uum(iu, lu, maz ) = u0uls(iu, lu)
       ENDDO
      ENDDO
*
      IF ( maz .EQ. 0 ) THEN
         DO lu = 1, ntau
            dfdt(lu)  = dfdts(lu)
            uavg(lu)  = uavgs(lu)
            flup(lu)  = flups(lu)
            rfldn(lu) = rfldns(lu)
            rfldir(lu)= rfldirs(lu)
         ENDDO
      ENDIF
      albedo = albedos
      fbeam  = fbeams
*
      RETURN
      END
*
      SUBROUTINE  dpsasymtx( A, EVEC, EVAL, M, IA, IEVEC, IER, WK,
     $                    AAD, EVECD, EVALD, WKD, UMU0 )
C
C    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======
C
C       SOLVES EIGENFUNCTION PROBLEM FOR real*8 ASYMMETRIC MATRIX
C       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE real*8.
C
C       THIS IS AN ADAPTATION OF A SUBROUTINE 'EIGRF' IN THE IMSL
C       LIBRARY TO USE real*8 INSTEAD OF COMPLEX ARITHMETIC, ACCOUNTING
C       FOR THE KNOWN FACT THAT THE EIGENVALUES AND EIGENVECTORS IN
C       THE DISCRETE ORDINATE SOLUTION ARE real*8.  OTHER CHANGES INCLUDE
C       PUTTING ALL THE CALLED SUBROUTINES IN-LINE, IN DELETING
C       THE PERFORMANCE INDEX CALCULATION, IN UPDATING MANY DO-LOOPS
C       TO FORTRAN77, AND IN CALCULATING THE MACHINE PRECISION
C       'TOL' INSTEAD OF SPECIFYING IT IN A DATA STATEMENT.
C
C       'EIGRF' IS BASED PRIMARILY ON 'EISPACK' ROUTINES.  THE MATRIX
C       IS FIRST BALANCED USING THE PARLETT-REINSCH ALGORITHM.  THEN
C       THE MARTIN-WILKINSON ALGORITHM IS APPLIED.
C
C       REFERENCES:
C          DONGARRA, J. AND C. MOLER, EISPACK -- A PACKAGE FOR SOLVING
C             MATRIX EIGENVALUE PROBLEMS, IN COWELL, ED., 1984:
C             SOURCES AND DEVELOPMENT OF MATHEMATICAL SOFTWARE,
C             PRENTICD-HALL, ENGLEWOOD CLIFFS, NJ
C         PARLETT AND REINSCH, 1969: BALANCING A MATRIX FOR CALCULATION
C             OF EIGENVALUES AND EIGENVECTORS, NUM. MATH. 13, 293-304
C         WILKINSON, J., 1965: THE ALGEBRAIC EIGENVALUE PROBLEM,
C             CLARENDON PRESS, OXFORD
C
C   I N P U T    V A R I A B L E S:
C
C        A    :  INPUT ASYMMETRIC MATRIX, DESTROYED AFTER SOLVED
C        M    :  ORDER OF -A-
C       IA    :  FIRST DIMENSION OF -A-
C    IEVEC    :  FIRST DIMENSION OF -EVEC-
C
C   O U T P U T    V A R I A B L E S:
C
C       EVEC  :  (UNNORMALIZED) EIGENVECTORS OF -A-
C                   ( COLUMN J CORRESPONDS TO EVAL(J) )
C
C       EVAL  :  (UNORDERED) EIGENVALUES OF -A- ( DIMENSION AT LEAST M )
C
C       IER   :  IF .NE. 0, SIGNALS THAT EVAL(IER) FAILED TO CONVERGE;
C                   IN THAT CASE EIGENVALUES IER+1,IER+2,...,M  ARE
C                   CORRECT BUT EIGENVALUES 1,...,IER ARE SET TO ZERO.
C
C   S C R A T C H   V A R I A B L E S:
C
C       WK    :  WORK AREA ( DIMENSION AT LEAST 2*M )
C       AAD    :  DOUBLE PRECISION STAND-IN FOR -A-
C       EVECD :  DOUBLE PRECISION STAND-IN FOR -EVEC-
C       EVALD :  DOUBLE PRECISION STAND-IN FOR -EVAL-
C       WKD   :  DOUBLE PRECISION STAND-IN FOR -WK-
C+---------------------------------------------------------------------+
C
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      real*8               A( IA,* ),  WK(*),
     $     EVAL(*), EVEC( IEVEC,* )
      DOUBLE PRECISION  AAD( IA,* ), WKD(*), EVALD(*), EVECD( IA,* )
      real*8  C1, DISCRI, Q, P, R, RNORM, S, SGN, TOL, W, X, Z
      real*8 zero, umu0
      DOUBLE PRECISION  D1MACH
      LOGICAL           NOCONV, NOTLAS
      DATA     C1 / 0.4375D0 /, C2/ 0.5D0 /, C3/ 0.75D0 /, C4/ 0.95D0 /,
     $         C5/ 16.D0 /, C6/ 256.D0 /, ZERO / 0.D0 /, ONE / 1.D0 /
C
C
      TOL = D1MACH(3)   

      IER = 0
         
      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )
     $     CALL ERRMSG( 'dpsasymtx--BAD INPUT VARIABLE(S)', .TRUE. )
C
C                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVAL(1) = A(1,1)
         EVEC(1,1) = 1.0D0
         RETURN
C
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( A(1,1) - A(2,2) )**2 + 4.D0 * A(1,2) * A(2,1)
         IF ( DISCRI.LT.0.0 )
     $      CALL ERRMSG( 'dpsasymtx--COMPLEX EVALS IN 2X2 CASE', .TRUE.)
         SGN = 1.D0
         IF ( A(1,1).LT.A(2,2) )  SGN = - 1.0D0
         EVAL(1) = 0.5D0 * ( A(1,1) + A(2,2) + SGN*SQRT(DISCRI) )
         EVAL(2) = 0.5D0 * ( A(1,1) + A(2,2) - SGN*SQRT(DISCRI) )
	
         EVEC(1,1) = 1.0D0
         EVEC(2,2) = 1.0D0
         IF ( A(1,1).EQ.A(2,2) .AND.
     $        (A(2,1).EQ.0.0.OR.A(1,2).EQ.0.0) )
     $        THEN
            RNORM = DABS(A(1,1))+ DABS(A(1,2))
     $           + DABS(A(2,1)) + DABS(A(2,2))
            W = TOL * RNORM
            EVEC(2,1) = A(2,1) / W
            EVEC(1,2) = - A(1,2) / W
         ELSE
            EVEC(2,1) = A(2,1) / ( EVAL(1) - A(2,2) )
            EVEC(1,2) = A(1,2) / ( EVAL(2) - A(1,1) )
         ENDIF
         RETURN
      END IF
C                               ** PUT S.P. MATRIX INTO D.P. MATRIX
      DO 1  J = 1, M
         DO 1  K = 1, M
            AAD( J,K ) = DBLE( A(J,K) )
    1 CONTINUE
C                                        ** INITIALIZE OUTPUT VARIABLES
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE
C                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
C                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
C                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M
C
30    KKK = K
         DO 70  J = KKK, 1, -1
            ROW = ZERO
            DO 40 I = 1, K
               IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )
40          CONTINUE
            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               GO TO 30
            END IF
70       CONTINUE
C                                     ** SEARCH FOR COLUMNS ISOLATING AN
C                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL = L
         DO 120 J = LLL, K
            COL = ZERO
            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )
90          CONTINUE
            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               GO TO 80
            END IF
120      CONTINUE
C                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE
C
140   NOCONV = .FALSE.
         DO 200 I = L, K
            COL = ZERO
            ROW = ZERO
            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + DABS( AAD(J,I) )
                  ROW = ROW + DABS( AAD(I,J) )
               END IF
150         CONTINUE
            F = ONE
            G = ROW / C5
            H = COL + ROW
160         IF ( COL.LT.G ) THEN
               F   = F * C5
               COL = COL * C6
               GO TO 160
            END IF
            G = ROW * C5
170         IF ( COL.GE.G ) THEN
               F   = F / C5
               COL = COL / C6
               GO TO 170
            END IF
C                                                         ** NOW BALANCE
            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE
C
      IF ( NOCONV ) GO TO 140
C                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1 .LT. L+1 ) GO TO 350
C                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H        = ZERO
         WKD(N+M) = ZERO
         SCALE    = ZERO
C                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + DABS(AAD(I,N-1))
210      CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WKD(I+M) = AAD(I,N-1) / SCALE
               H = H + WKD(I+M)**2
220         CONTINUE
            G = - SIGN( SQRT(H), WKD(N+M) )
            H = H - WKD(N+M) * G
            WKD(N+M) = WKD(N+M) - G
C                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WKD(I+M) * AAD(I,J)
230            CONTINUE
               DO 240 I = N, K
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240            CONTINUE
250         CONTINUE
C                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WKD(J+M) * AAD(I,J)
260            CONTINUE
               DO 270 J = N, K
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WKD(N+M)  = SCALE * WKD(N+M)
            AAD(N,N-1) = SCALE * G
         END IF
290   CONTINUE
C
      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AAD(N1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WKD(N1+M)
            DO 300 I = N2, K
               WKD(I+M) = AAD(I,N)
300         CONTINUE
            IF ( N1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N1, K
                     G = G + WKD(I+M) * EVECD(I,J)
310               CONTINUE
                  G = G / F
                  DO 320 I = N1, K
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE
C
350   CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + DABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) then
	   EVALD(I) = AAD(I,I)
	 endif  
370   CONTINUE
      N = K
      T = ZERO
C                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
C                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410
400   CONTINUE
C
410   X = AAD(N,N)
      IF ( LB.EQ.N ) THEN
C                                        ** ONE EIGENVALUE FOUND
         AAD(N,N)  = X + T
         EVALD(N) = AAD(N,N)
         N = N1
         GO TO 380
      END IF
C
      Y = AAD(N1,N1)
      W = AAD(N,N1) * AAD(N1,N)
      IF ( LB.EQ.N1 ) THEN
C                                        ** TWO EIGENVALUES FOUND
         P = (Y-X) * C2
         Q = P**2 + W
         Z = DSQRT( DABS(Q) )
         AAD(N,N) = X + T
         X = AAD(N,N)
         AAD(N1,N1) = Y + T
C                                        ** real*8 PAIR
         Z = P + SIGN(Z,P)
         EVALD(N1) = X + Z
         EVALD(N)  = EVALD(N1)
         IF ( Z.NE.ZERO ) EVALD(N) = X - (W / Z)

         X = AAD(N,N1)
C                                  ** EMPLOY SCALE FACTOR IN CASE
C                                  ** X AND Z ARE VERY SMALL
         R = DSQRT( X*X + Z*Z )
         P = X / R
         Q = Z / R
C                                             ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AAD(N1,J)
            AAD(N1,J) = Q * Z + P * AAD(N,J)
            AAD(N,J)  = Q * AAD(N,J) - P * Z
420      CONTINUE
C                                             ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AAD(I,N1)
            AAD(I,N1) = Q * Z + P * AAD(I,N)
            AAD(I,N)  = Q * AAD(I,N) - P * Z
430      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVECD(I,N1)
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440      CONTINUE
C
         N = N2
         GO TO 380
      END IF
C
      IF ( IN.EQ.30 ) THEN
C                    ** NO CONVERGENCE AFTER 30 iterationS; SET ERROR
C                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE
         IER = N
         GO TO 670
      END IF
C                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AAD(I,I) = AAD(I,I) - X
450      CONTINUE
         S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
         X = C3 * S
         Y = X
         W = - C1 * S**2
      END IF
C
      IN = IN + 1
C                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
C
      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AAD(I,I)
         R = X - Z
         S = Y - Z
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
         Q = AAD(I+1,I+1) - Z - R - S
         R = AAD(I+2,I+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
         VV = DABS(P)*( DABS(AAD(I-1,I-1)) + DABS(Z) +
     $        DABS(AAD(I+1,I+1)) )
         IF ( UU .LE. TOL*VV ) GO TO 470
460   CONTINUE
C
470   CONTINUE
      AAD(I+2,I) = ZERO
      DO 480 J = I+3, N
         AAD(J,J-2) = ZERO
         AAD(J,J-3) = ZERO
480   CONTINUE
C
C             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N
C
      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
         ELSE
            P = AAD(KA,KA-1)
            Q = AAD(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)
            X = DABS(P) + DABS(Q) + DABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            AAD(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
C                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AAD(KA,J) + Q * AAD(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AAD(KA+2,J)
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z
            END IF
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y
            AAD(KA,J)   = AAD(KA,J)   - P * X
490      CONTINUE
C                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AAD(II,KA+2)
               AAD(II,KA+2) = AAD(II,KA+2) - P * R
            END IF
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q
            AAD(II,KA)   = AAD(II,KA) - P
500      CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVECD(II,KA+2)
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
            END IF
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
            EVECD(II,KA)   = EVECD(II,KA) - P
510      CONTINUE
C
520   CONTINUE
      GO TO 390
C                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE real*8 VECTOR
530   CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE
C                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS
C
         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE
C                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 600  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600         CONTINUE
         END IF
C
      END IF
C
      DO 620 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620   CONTINUE
C                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = WKD(I)
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE
C
      DO 660 I = K+1, M
         J = WKD(I)
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE
C                         ** PUT RESULTS INTO OUTPUT ARRAYS
  670 CONTINUE
      DO 680 J = 1, M
         EVAL( J ) = EVALD(J)
      DO 680 K = 1, M
            EVEC( J,K ) = EVECD(J,K)
680   CONTINUE
     
C
      RETURN
      END


     
      
      SUBROUTINE  dpschekin( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     $     WVNMHI, USRTAU, NTAU, UTAU, NSTR, 
     $     usrang, numu, umu, FBEAM, 
     $     UMU0, SPHER, FISOT, ALBEDO, BTEMP,
     $     TTEMP, TEMIS, NOPLNK, onlyfl, ierror, quiet,
     $     MAXCLY, MAXULV, maxumu, MAXCMU, MXCLY, MXULV,  MXCMU, 
     $     mxumu,  TAUC, ZD )
C
C           CHECKS THE INPUT DIMENSIONS AND VARIABLES
C
      implicit double precision (A-H, O-Z)
      LOGICAL  WRTBADQ, WRTDIMQ
      LOGICAL  NOPLNK, onlyfl, SPHER, usrang, USRTAU, INPERR,
     $     quiet
      INTEGER  MAXCLY, MAXULV, MAXCMU, NLYR, NSTR, NTAU, numu,
     $         MXCMU, MXCLY, MXULV, ierror(*), maxumu
      real*8     ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     $         PMOM( 0:MAXCMU, MAXCLY ), SSALB( MAXCLY ), 
     $         TEMPER( 0:MAXCLY ), TEMIS, TTEMP, WVNMLO, 
     $         WVNMHI, UMU0, UTAU( MAXULV ), TAUC( 0:* ), ZD( 0:* ),
     $         umu(maxumu), umumin
*
      inperr = .FALSE.
      IF ( nlyr.LT.1 ) THEN
         inperr = wrtbadq(quiet, 'nlyr' )
         ierror(1) = 1
      ENDIF
      IF ( nlyr.GT.maxcly ) THEN
         inperr = wrtbadq(quiet, 'maxcly' )
         ierror(2) = 1
      ENDIF
*
      DO 10  LC = 1, NLYR
         IF ( dtauc(lc).LT.0.0D0 )  THEN
            inperr = wrtbadq( quiet, 'dtauc' )
            ierror(3) = ierror(3) + 1
         ENDIF
         IF ( ssalb(lc).LT.0.0 .OR. ssalb(lc).GT.1.0 ) THEN
            inperr = wrtbadq( quiet, 'ssalb' )
            ierror(4) = ierror(4) + 1
         ENDIF
         IF ( .NOT.  noplnk )  THEN
            IF( lc.EQ.1 .AND. temper(0).LT.0.0 ) THEN
               inperr = wrtbadq( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
            IF( temper(lc).LT.0.0 ) THEN 
               inperr = wrtbadq( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
         ENDIF
         DO k = 0, nstr
            IF( pmom(k,lc).LT.-1.0 .OR. pmom(k,lc).GT.1.0 ) THEN
               inperr = wrtbadq( quiet, 'pmom' )
               ierror(6) = ierror(6) + 1
            ENDIF
         ENDDO
10    CONTINUE
*
      IF ( spher ) THEN
         DO 11 lc = 1, nlyr
            IF ( zd(lc) .GT. zd(lc-1) ) THEN
               inperr = wrtbadq( quiet, 'zd' )
               ierror(7) = ierror(7) + 1
            ENDIF
 11      CONTINUE     
      ENDIF
*
      IF ( USRTAU )  THEN
         IF ( ntau.LT.1 ) THEN
            inperr = wrtbadq( quiet, 'ntau' )
            ierror(8) = 1
         ENDIF
         IF ( maxulv.LT.ntau ) THEN
            inperr = wrtbadq( quiet, 'maxulv' )
            ierror(9) = 1
         ENDIF
         DO 20  lu = 1, ntau
            IF( DABS(utau(lu)-tauc(nlyr)).LE.1.D-4)
     $           utau(lu) = tauc(nlyr)
cbm   Added a relative check to the absolute 
cbm   (required for large optical thickness)
            IF ( tauc(nlyr) .GT. 0.D0 ) THEN
               IF( DABS((utau(lu)-tauc(nlyr))/tauc(nlyr)).LE.1.D-6)
     $              utau(lu) = tauc(nlyr)
            ENDIF
            IF( utau(lu).LT.0.0 .OR. utau(lu).GT.tauc(nlyr) ) THEN
               inperr = wrtbadq( quiet, 'utau' )
               ierror(10) = ierror(10) + 1
            ENDIF
20       CONTINUE
      ELSE
         IF ( maxulv.LT.nlyr+1 ) THEN
            inperr = wrtbadq( quiet, 'maxulv' )
            ierror(11) = 1
         ENDIF
      END IF
      IF ( NSTR.LT.2 .OR. MOD(NSTR,2).NE.0 ) THEN
         inperr = wrtbadq( quiet, 'nstr' )
         ierror(23) = 1
      ENDIF
      IF ( NSTR.GT.MAXCMU ) THEN
         inperr = wrtbadq( quiet, 'maxcmu' )
         ierror(24) = 1
      ENDIF
*
      IF ( usrang ) THEN
         IF ( numu .LT. 0 ) THEN
           inperr = wrtbadq( quiet, 'numu' )
           ierror(25) = 1
        ENDIF
        IF ( .NOT. onlyfl .AND. numu .EQ. 0 ) THEN
           inperr = wrtbadq( quiet, 'numu' )
           ierror(26) = 1
        ENDIF
        IF ( numu .GT. maxumu ) THEN
           inperr = wrtbadq( quiet, 'maxumu' )
           ierror(27) = 1
        ENDIF
        DO 30 iu = 1, numu
           IF( umu(iu).LT.-1.0 .OR. umu(iu).GT.1.0 .OR.
     $          umu(iu).EQ. 0.D0) THEN
           
              inperr = wrtbadq( quiet, 'umu' )
              ierror(28) = ierror(28) + 1
           ENDIF
           IF( iu.GT.1 .AND. umu(iu).LT.umu(iu-1) ) THEN
              inperr = wrtbadq( quiet, 'umu' )
              ierror(29) = ierror(29) + 1
           ENDIF
 30     CONTINUE
      ELSE
         IF( maxumu .LT. nstr ) THEN
           inperr = wrtbadq( quiet, 'maxumu' )
           ierror(30) = 1
        ENDIF
      ENDIF
*
      IF ( fbeam.LT.0.0 ) THEN
         inperr = wrtbadq( quiet, 'fbeam' )
         ierror(12) = 1
      ENDIF
      umumin = 0.D0
      IF ( spher ) umumin = -1.D0
      IF ( fbeam.GT.0.0 .AND. ( umu0.LE.umumin .OR. umu0.GT.1.0 ) )
     $     THEN
         inperr = wrtbadq( quiet, 'umu0' )
         ierror(13) = 1
      ENDIF
C
cbm   Be a bit more conservative; SZAs up to 0.1 degrees may cause troubles
      IF ( FBEAM.GT. 0.0 .AND. umu0 .GT. 0.9999985D0) then 
          CALL ERRMSG( 'Does not work for umu0=1.0', .TRUE. )
      ENDIF	   
cbm   Original code:
cbm      IF ( FBEAM.GT. 0.0 .AND. umu0.EQ. 1.D0) then 
cbm          CALL ERRMSG( 'Does not work for umu0=1.0', .TRUE. )
cbm      ENDIF	   
C
      IF ( fisot.LT.0.D0 ) THEN
         inperr = wrtbadq( quiet, 'fisot' )
         ierror(14) = 1
      ENDIF
      IF ( albedo.LT.0.0 .OR. albedo.GT.1.0 ) THEN
         inperr = wrtbadq( quiet, 'albedo' )
         ierror(15) = 1
      ENDIF
*
      IF ( .NOT.NOPLNK )  THEN
         IF ( wvnmlo.LT.0.0 .OR. wvnmhi.LT.wvnmlo ) THEN
            inperr = wrtbadq( quiet, 'wvnmlo,hi' )
            ierror(16) = 1
         ENDIF
         IF ( temis.LT.0.0 .OR. temis.GT.1.0 ) THEN
            inperr = wrtbadq( quiet, 'temis' )
            ierror(17) = 1
         ENDIF
         IF ( btemp.LT.0.0 ) THEN
            inperr = wrtbadq( quiet, 'btemp' )
            ierror(18) = 1
         ENDIF
         IF ( ttemp.LT.0.0 ) THEN
            inperr = wrtbadq( quiet, 'ttemp' )
            ierror(19) = 1
         ENDIF
      END IF
*
      IF ( mxcly.LT.nlyr ) THEN
         inperr = wrtdimq( quiet, 'mxcly', nlyr )
         ierror(20) = 1
      ENDIF
      IF ( usrtau .AND. mxulv.LT.ntau ) THEN
         inperr = wrtdimq( quiet, 'mxulv', ntau )
         ierror(21) = 1
      ENDIF
      IF ( .NOT.usrtau .AND. mxulv.LT.nlyr+1 ) THEN
         inperr = wrtdimq( quiet, 'mxulv', nlyr+1 )
         ierror(22) = 1
      ENDIF
      IF ( MXCMU.LT.NSTR ) THEN
         INPERR = WRTDIMQ( quiet, 'MXCMU', NSTR )
         ierror(31) = 1
      ENDIF
      IF ( usrang .AND. mxumu.LT.numu ) THEN
         inperr = wrtdimq( quiet, 'mxumu', numu )
         ierror(32) = 1
      ENDIF
      IF ( .NOT. usrang .AND. mxumu.LT.nstr ) THEN
         inperr = wrtdimq( quiet, 'mxumu', nstr )
         ierror(33) = 1
      ENDIF
*
      IF ( inperr )
     $   CALL ERRMSG( 'SDISORT--INPUT AND/OR DIMENSION ERRORS', .TRUE.)
*
      DO 100  LC = 1, NLYR
         IF ((.NOT.noplnk .AND.DABS(temper(lc)-temper(lc-1)).GT.50.0) 
     $          .AND. .NOT. quiet )
     $          CALL errmsg( 'schekin--vertical temperature step may'
     $          // ' be too large for good accuracy', .FALSE. )
100   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE  dpsfluxes(CH, CMU, CWT, FBEAM, GC, il, KK, LAYRU, LL,
     $                    LYRCUT, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,
     $                    PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                    XR0, XR1, XRA, ZPLK0, ZPLK1, ZPLKA, ZBEAMA, 
     $                    ZBEAM0, ZBEAM1, DFDT, FLUP, FLDN, FLDIR, 
     $                    RFLDIR, RFLDN, UAVG, U0C, MAXULV,
     $                    uavgdn, uavgso, uavgup)
*
*     Modified from original disort sfluxes routine to include il in the
*     parameter list and to use the EXPonential-linear approximation
*     for the beam source.
*
C
C       CALCULATES THE RADIATIVE FLUXES, MEAN INTENSITY, AND FLUX
C       DERIVATIVE WITH RESPECT TO OPTICAL DEPTH FROM THE M=0 INTENSITY
C       COMPONENTS (THE AZIMUTHALLY-AVERAGED INTENSITY)
C
C    I N P U T     V A R I A B L E S:
C
C      
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL       :  CONSTANTS OF INTEGsratioN IN EQ. SC(1), OBTAINED
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       NN       :  ORDER OF DOUBLD-GAUSS QUADRATURE (NSTR/2)
C       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
C                     OPTICAL DEPTH EXCEEDS -ABSCUT-
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
C       XR0      :  EXPANSION OF THERMAL SOURCE FUNCTION IN EQ. KS(4)
C       XR1      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. KS(4)
C       XRA      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. KS(4)
C       ZBEAMA   :  BEAM SOURCE COEFFICIENT ALFA IN EQ.KS(6)
C       ZBEAM0   :  BEAM SOURCE VECTORS Y-SUB-ZERO IN EQ.KS(6)
C       ZBEAM1   :  BEAM SOURCE VECTORS Y-SUB-ONE IN EQ.KS(6)
C       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. KS(7-8)
C       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. KS(7-8)
C       ZPLKA    :  THERMAL SOURCE COEFFICIENT -ZA-, 
C                   ALFA IN EQ. KS(6)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       U0C      :  AZIMUTHALLY AVERAGED INTENSITIES
C                   ( AT POLAR QUADRATURE ANGLES )
C       (RFLDIR, RFLDN, FLUP, DFDT, UAVG ARE 'DISORT' OUTPUT VARIABLES)
C
C   I N T E R N A L       V A R I A B L E S:
C
C       DIRINT   :  DIRECT INTENSITY ATTENUATED
C       FDNTOT   :  TOTAL DOWNWARD FLUX (DIRECT + DIFFUSE)
C       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)
C       FLDN     :  DIFFUSE DOWN-FLUX (DELTA-M SCALED)
C       FNET     :  NET FLUX (TOTAL-DOWN - DIFFUSD-UP)
C       FACT     :  EXP( - UTAUPR / UMU0 )
C       PLSORC   :  PLANCK SOURCE FUNCTION (THERMAL)
C       ZINT     :  INTENSITY OF m = 0 CASE, IN EQ. SC(1)
C+---------------------------------------------------------------------+
C
      implicit double precision (A-H, O-Z)
      LOGICAL LYRCUT, PRNT(*)
      real*8    DFDT(*), FLUP(*), FLDIR(*), FLDN(*),
     $     RFLDIR(*), RFLDN(* ), 
     $     U0C( MXCMU,MXULV ), UAVG(*), uavgdn(*), uavgso(*),
     $     uavgup(*), umu0
      INTEGER LAYRU(*)
      real*8    CH(*), CMU(*), CWT(*),
     $     FBEAM, GC( MXCMU,MXCMU,* ), KK( MXCMU,* ),
     $     LL( MXCMU,* ), pi, SSALB(*), TAUCPR( 0:* ),
     $     UTAU(*), UTAUPR(*), XR0(*), XR1(*), XRA(*), 
     $     ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), ZPLKA( * ) ,
     $     ZBEAMA(*), ZBEAM0(MXCMU,*), ZBEAM1(MXCMU,*) 
      real*8 ang1, ang2, dirint, fact, fdntot,
     $     fnet, plsorc, zint
C
C
      
      
      IF ( PRNT(2) )  WRITE( *,1010 )
C                                          ** ZERO DISORT OUTPUT ARRAYS
      CALL  dpszeroit( U0C, MXULV*MXCMU )
      CALL  dpszeroit( RFLDIR, MAXULV )
      CALL  dpszeroit( FLDIR,  MXULV )
      CALL  dpszeroit( RFLDN,  MAXULV )
      CALL  dpszeroit( FLDN,   MXULV )
      CALL  dpszeroit( FLUP,   MAXULV )
      CALL  dpszeroit( UAVG,   MAXULV )
      CALL  dpsZEROIT( UAVGdn, MAXULV )
      CALL  dpsZEROIT( UAVGso, MAXULV )
      CALL  dpsZEROIT( UAVGup, MAXULV )
      CALL  dpszeroit( DFDT,   MAXULV )
C                                        ** LOOP OVER USER LEVELS
      DO 100  LU = 1, NTAU
C
         LYU = LAYRU(LU)
C
         IF ( LYRCUT .AND. LYU.GT.NCUT ) THEN
C                                                ** NO RADIATION REACHES
C                                                ** THIS LEVEL
            FDNTOT = 0.0D0
            FNET   = 0.0D0
            PLSORC = 0.0D0
            GO TO 90
         END IF
*
         IF ( FBEAM.GT.0.D0 )  THEN
            FACT  = DEXP( - UTAUPR(LU) / CH(LYU) )
            DIRINT = FBEAM * FACT
            FLDIR(  LU ) = DABS(UMU0) * ( FBEAM * FACT )
            RFLDIR( LU ) = DABS(UMU0)*FBEAM*DEXP( - UTAU( LU ) /
     $                                           CH(LYU) )         
*           
                               	    
         ELSE
            DIRINT = 0.D0
            FLDIR(  LU ) = 0.D0
            RFLDIR( LU ) = 0.D0
	   
         END IF
C
         DO 20  IQ = 1, NN
C
            ZINT = 0.D0
            DO 10  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $               DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
10          CONTINUE
            DO 11  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                DEXP( - KK(JQ,LYU) *
     $              (UTAUPR(LU) - TAUCPR(LYU-1)) )
11          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.D0 .OR. il.GT. 0 )  U0C( IQ,LU ) = ZINT + 
     $                           DEXP(-ZBEAMA(LYU)*UTAUPR(LU)) *
     $                           ( ZBEAM0(IQ,LYU) +
     $                             ZBEAM1(IQ,LYU)*UTAUPR(LU) )
            U0C( IQ,LU ) = U0C( IQ,LU ) +
     $                             DEXP(-ZPLKA(LYU)*UTAUPR(LU))*
     $                         (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
            UAVG(LU) = UAVG(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )
            UAVGdn(LU) = UAVGdn(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )
            FLDN(LU) =FLDN(LU) + CWT(NN+1-IQ)*CMU(NN+1-IQ) * U0C(IQ,LU)
          
20       CONTINUE
C
         DO 40  IQ = NN+1, NSTR
C
            ZINT = 0.D0
            DO 30  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                DEXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
30          CONTINUE
            DO 31  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                DEXP( - KK(JQ,LYU) *
     $              (UTAUPR(LU) - TAUCPR(LYU-1)) )
31          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.0D0 .OR. il.GT.0 ) U0C( IQ,LU )=ZINT + 
     $                           DEXP(-ZBEAMA(LYU)*UTAUPR(LU))*
     $                           ( ZBEAM0(IQ,LYU) +
     $                             ZBEAM1(IQ,LYU)*UTAUPR(LU) )
            U0C( IQ,LU ) = U0C( IQ,LU ) + 
     $                             DEXP(-ZPLKA(LYU)*UTAUPR(LU))*
     $                       (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
            UAVG(LU) = UAVG(LU) + CWT(IQ-NN) * U0C( IQ,LU )
            UAVGup(LU) = UAVGup(LU) + CWT(IQ-NN) * U0C( IQ,LU )    
            FLUP(LU)=FLUP(LU) + CWT(IQ-NN)*CMU(IQ-NN)*U0C( IQ,LU )
         
40       CONTINUE
C
         FLUP( LU )  = 2.D0 * PI * FLUP( LU )
         FLDN( LU )  = 2.D0 * PI * FLDN( LU )
         FDNTOT = FLDN( LU ) + FLDIR( LU )
         FNET   = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU ) = (2.D0 * PI * UAVG(LU) + DIRINT)/( 4.D0*PI )
         UAVGdn( LU ) = ( 2.0 * PI * UAVGdn(LU) ) / ( 4.*PI )
         UAVGso( LU ) =  DIRINT  / ( 4.*PI )
         UAVGup( LU ) = ( 2.0 * PI * UAVGup(LU) ) / ( 4.*PI )
         PLSORC = DEXP(-XRA(LYU)*UTAUPR(LU))*
     $                          (XR0(LYU) + XR1(LYU)* UTAUPR(LU))
         DFDT( LU ) = ( 1.D0-SSALB(LYU) ) * 4.D0*PI*
     $        ( UAVG(LU) - PLSORC )
 90      IF( PRNT(2) )  WRITE( *,1020 ) UTAU(LU), LYU, RFLDIR(LU),
     $                                 RFLDN(LU), FDNTOT, FLUP(LU),
     $                                 FNET, UAVG(LU), PLSORC, DFDT(LU)
100   CONTINUE
C
      IF ( PRNT(3) )  THEN
         WRITE ( *,1100 )
         DO 200  LU = 1, NTAU
            WRITE( *,1110 )  UTAU( LU )
            DO  200  IQ = 1, NN
               ANG1 = 180.D0/PI * DACOS( CMU(2*NN-IQ+1) )
               ANG2 = 180.D0/PI * DACOS( CMU(IQ) )
               WRITE( *,1120 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                         ANG2, CMU(IQ),        U0C(IQ+NN,LU)
200      CONTINUE
      END IF
C
1010  FORMAT( //, 21X,
     $ '<----------------------- sfluxes ----------------------->', /,
     $ '   OPTICAL  COMPU    DOWNWARD    DOWNWARD    DOWNWARD     ',
     $ ' UPWARD                    MEAN      PLANCK   D(NET FLUX)', /,
     $ '     DEPTH  LAYER      DIRECT     DIFFUSE       TOTAL     ',
     $ 'DIFFUSE         NET   INTENSITY      SOURCE   / D(OP DEP)', / )
1020  FORMAT( F10.4, I7, 1P,7(E12.3), E14.3 )
1100  FORMAT( //, ' ******** AZIMUTHALLY AVERAGED INTENSITIES',
     $      ' ( AT POLAR QUADRATURE ANGLES ) *******' )
1110  FORMAT( /, ' OPTICAL DEPTH =', F10.4, //,
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY',
     $  '     ANGLE (DEG)   COS(ANGLE)     INTENSITY' )
1120  FORMAT( 2( 0P,F16.4, F13.5, 1P,E14.3 ) )
C
      RETURN
      END
      SUBROUTINE  dpslepoly( NMU, M, MAXMU, TWONM1, MU, YLM )
C
C       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
C       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C       PLM = P-SUB-L-SUPER-M AS
C
C             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU)
C
C       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1.
C       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE
C       FROM A PRIOR CALL TO THE ROUTINE.
C
C       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
C                  High-Order Associated Legendre Polynomials,
C                  J. Quant. Spectrosc. Radiat. Transfer 10,
C                  557-562, 1970.  (hereafter D/A)
C
C       METHOD: Varying degree recurrence relationship.
C
C       NOTE  The D/A formulas are transformed by
C               setting  M = n-1; L = k-1.
C       NOTE 2: Assumes that routine is called first with  M = 0,
C               then with  M = 1, etc. up to  M = TWONM1.
C       NOTE 3: Loops are written in such a way as to vectorize.
C
C  I N P U T     V A R I A B L E S:
C
C       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
C       M      :  ORDER OF -YLM-
C       MAXMU  :  FIRST DIMENSION OF -YLM-
C       TWONM1 :  MAX DEGREE OF -YLM-
C       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-
C       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED
C
C  O U T P U T     V A R I A B L E:
C
C       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
C                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
C+---------------------------------------------------------------------+
      implicit double precision (A-H, O-Z)
      real*8     MU(*), YLM( 0:MAXMU,* )
      INTEGER  M, NMU, TWONM1
      PARAMETER  ( MAXSQT = 1000 )
      real*8     SQT( MAXSQT )
      real*8 tmp1, tmp2
      LOGICAL  PASS1
      SAVE  SQT, PASS1
      DATA  PASS1 / .TRUE. /
C
C
      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         DO 1  NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT(NS) )
    1    CONTINUE
      ENDIF
C
      IF ( 2*TWONM1 .GT. MAXSQT )
     $   CALL ERRMSG( 'slepoly--NEED TO INCREASE PARAM MAXSQT', .TRUE. )
C
      IF ( M .EQ. 0 )  THEN
C                             ** UPWARD RECURRENCE FOR ORDINARY
C                             ** LEGENDRE POLYNOMIALS
         DO  10  I = 1, NMU
            YLM( 0,I ) = 1.D0
            YLM( 1,I ) = MU( I )
10       CONTINUE
         DO  20  L = 2, TWONM1
            DO  20  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $                      - ( L-1 ) * YLM( L-2,I ) ) / L
20       CONTINUE
C
      ELSE
C
         DO  30  I = 1, NMU
C                               ** Y-SUB-M-SUPER-M; DERIVED FROM
C                               ** D/A EQS. (11,12)
C
            YLM( M,I) = - SQT( 2*M-1 ) / SQT( 2*M )
     $                  * SQRT( 1. - MU(I)**2 ) * YLM( M-1,I )
C
C                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM
C                              ** D/A EQS. (13,14) USING EQS. (11,12)
C
            YLM( M+1,I ) = SQT( 2*M+1 ) * MU(I) * YLM( M,I )
30       CONTINUE
C                                   ** UPWARD RECURRENCE; D/A EQ. (10)
         DO  40  L = M+2, TWONM1
            TMP1 = SQT( L-M ) * SQT( L+M )
            TMP2 = SQT( L-M-1 ) * SQT( L+M-1 )
            DO  40  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $                        - TMP2 * YLM( L-2,I ) ) / TMP1
40       CONTINUE
C
      END IF
C
      RETURN
      END
      SUBROUTINE dpsprtinp( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     $                    WVNMHI, NTAU, UTAU, NSTR, numu, umu,
     $                    FBEAM, UMU0, spher, nil,
     $                    FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, 
     $                    NOPLNK, onlyfl, FLYR, LYRCUT, OPRIM, 
     $                    TAUC, TAUCPR, MAXCMU, PRTMOM )
C
C        PRINT VALUES OF INPUT VARIABLES
C
      implicit double precision (A-H, O-Z)
      LOGICAL  DELTAM, LYRCUT, NOPLNK, onlyfl, PRTMOM, spher
      real*8 albedo, btemp, DTAUC(*), fbeam,
     $     fisot, FLYR(*), OPRIM(*), 
     $     PMOM( 0:MAXCMU,* ), SSALB(*), UTAU(*), TAUC( 0:* ),
     $     TAUCPR( 0:* ), temis, TEMPER( 0:* ), ttemp, 
     $     umu(*), umu0, wvnmhi, wvnmlo
      real*8 yessct
      integer nil
C
C
      WRITE( *,1010 )  NSTR, NLYR
      WRITE( *,1030 )  NTAU, (UTAU(LU), LU = 1, NTAU)
      IF ( .NOT. onlyfl )
     $      WRITE( *,1040 ) numu, ( umu(iu), iu = 1, numu )
      IF ( NOPLNK  )  WRITE( *,1100 )
      WRITE( *,1060 ) FBEAM, UMU0,  FISOT
      WRITE( *,1080 ) ALBEDO
      IF ( .NOT.NOPLNK )  WRITE( *,1110 ) WVNMLO, WVNMHI, BTEMP,
     $                                    TTEMP, TEMIS
      IF ( spher .AND. nil.EQ.0 ) WRITE(*,1111)
      IF ( spher .AND. nil.GT.0 ) WRITE(*,1112)
      IF ( .NOT. spher )          WRITE(*,1113)
      IF ( DELTAM )      WRITE( *,1120 )
      IF ( .NOT.DELTAM ) WRITE( *,1130 )
      IF ( LYRCUT )  WRITE( *,1170 )
      IF( .NOT.NOPLNK )  WRITE ( *,1190 )
      IF(      NOPLNK )  WRITE ( *,1191 )
      YESSCT = 0.D0
      DO 10 LC = 1, NLYR
         YESSCT = YESSCT + SSALB(LC)
         IF( .NOT.NOPLNK )
     $       WRITE( *,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),
     $                    OPRIM(LC), PMOM(1,LC), TEMPER(LC-1)
         IF( NOPLNK )
     $       WRITE( *,1200 )  LC, DTAUC(LC), TAUC(LC), SSALB(LC),
     $                    FLYR(LC), TAUCPR(LC)-TAUCPR(LC-1), TAUCPR(LC),
     $                    OPRIM(LC), PMOM(1,LC)
 10   CONTINUE
      IF( .NOT.NOPLNK )  WRITE( *,1210 ) TEMPER(NLYR)
C
      IF( PRTMOM .AND. YESSCT.GT.0.D0 )  THEN
         WRITE( *, '(/,A)' )  ' LAYER   PHASE FUNCTION MOMENTS'
         DO 20 LC = 1, NLYR
            IF( SSALB(LC).GT.0.D0 )
     $          WRITE( *,1300 )  LC, ( PMOM(K,LC), K = 0, NSTR )
 20      CONTINUE
      ENDIF
C
      RETURN
C
 1010 FORMAT ( /, ' NO. STREAMS =', I4,
     $     '     NO. COMPUTATIONAL LAYERS =', I4 )
 1030 FORMAT( I4,' USER OPTICAL DEPTHS :',10F10.4, /, (26X,10F10.4) )
 1040 FORMAT( I4,' User polar angles cosines:',10F9.5,/,(31X,10F9.5) )
 1060 FORMAT( '    INCIDENT BEAM WITH INTENSITY =', 1P,E11.3, ' AND',
     $     ' POLAR ANGLE COSINE = ', 0P,F8.5,
     $     /,'    PLUS ISOTROPIC INCIDENT INTENSITY =', 1P,E11.3 )
 1080 FORMAT( '    BOTTOM ALBEDO (LAMBERTIAN) =', 0P,F8.4 )
 1100 FORMAT( ' NO THERMAL EMISSION' )
 1110 FORMAT( '    THERMAL EMISSION IN WAVENUMBER INTERVAL :',
     $     2F14.4,/
     $     ,'    BOTTOM TEMPERATURE =',
     $     F10.2, '     TOP TEMPERATURE ='
     $     ,F10.2,'    TOP EMISSIVITY =', F8.4 )
 1111 FORMAT( ' Uses pseudo-spherical geometry' )
 1112 FORMAT( ' Uses spherical geometry' )
 1113 FORMAT( ' Uses planD-parallel geometry' )
 1120 FORMAT( ' USES DELTA-M METHOD' )
 1130 FORMAT( ' DOES NOT USE DELTA-M METHOD' )
 1150 FORMAT( ' CALCULATE sfluxes AND INTENSITIES' )
 1170 FORMAT( ' SETS RADIATION = 0 BELOW ABSORPTION OPTICAL DEPTH 10' )
 1190 FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,
     $ '                   TOTAL    SINGLE                           ',
     $     'TOTAL    SINGLE', /,
     $     '       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',
     $     'OPTICAL   OPTICAL   SCATTER    ASYMM', /,
     $     '         DEPTH     DEPTH    ALBEDO    FRACTION     ',
     $     'DEPTH     DEPTH    ALBEDO   FACTOR   TEMPERATURE' )
 1191 FORMAT( /, 37X, '<------------- DELTA-M --------------->', /,
     $'                   TOTAL    SINGLE                           ',
     $     'TOTAL    SINGLE', /,
     $     '       OPTICAL   OPTICAL   SCATTER   TRUNCATED   ',
     $     'OPTICAL   OPTICAL   SCATTER    ASYMM', /,
     $     '         DEPTH     DEPTH    ALBEDO    FRACTION     ',
     $     'DEPTH     DEPTH    ALBEDO   FACTOR' )
 1200 FORMAT( I4, 2F10.4, F10.5, F12.5, 2F10.4, F10.5, F9.4,F14.3 )
 1210 FORMAT( 85X, F14.3 )
 1300 FORMAT( I6, 10F11.6, /, (6X,10F11.6) )
C     
      END
      SUBROUTINE  dpsqgausn( M, GMU, GWT )
C
C       COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE
C       (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)
C
C       REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
C                   Integration, Academic Press, New York, pp. 87, 1975.
C
C          METHOD:  Compute the abscissae as roots of the Legendre
C                   Polynomial P-SUB-N using a cubically convergent
C                   refinement of Newton's method.  Compute the
C                   weights from EQ. 2.7.3.8 of Davis/Rabinowitz.
C
C        ACCURACY:  at least 13 significant digits
C
C
C  I N P U T :    M       ORDER OF QUADRATURE RULE
C
C  O U T P U T :  GMU(I)  I = 1 TO M,    ARRAY OF ABSCISSAE
C                 GWT(I)  I = 1 TO M,    ARRAY OF WEIGHTS
C
C  I N T E R N A L    V A R I A B L E S:
C
C    PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS
C    PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL
C    P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL
C    TOL       : CONVERGENCE CRITERION FOR LEGENDRE POLY ROOT ITERATION
C    X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-
C                CONVERGENT VERSION OF NEWTON'S METHOD
C                ( SEEKING ROOTS OF LEGENDRE POLYNOMIAL )
C+---------------------------------------------------------------------+
      implicit double precision (A-H, O-Z)
      real*8     CONA, GMU(*), GWT(*), PI, T
      INTEGER  LIM, M, NP1
      DOUBLE   PRECISION  D1MACH
      DOUBLE   PRECISION  EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD,
     $                    TMP, TOL, X, XI
      SAVE tol, pi
      DATA     PI / 0.0D0 /
C
C
      IF ( PI.EQ.0.0D0 )  THEN
         PI = 2.D0 * DASIN(1.D0)
         TOL = 10.D0 * D1MACH(3)  
      END IF
C
      IF ( M.LE.1 )  THEN
         M = 1
         GMU( 1 ) = 0.5D0
         GWT( 1 ) = 1.D0
         RETURN
      END IF
C
      EN   = M
      NP1  = M + 1
      NNP1 = M * NP1
      CONA = FLOAT( M-1 ) / ( 8 * M**3 )
C                                        ** INITIAL GUESS FOR K-TH ROOT
C                                        ** OF LEGENDRE POLYNOMIAL, FROM
C                                        ** DAVIS/RABINOWITZ (2.7.3.3A)
      LIM  = M / 2
      DO 30  K = 1, LIM
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
         X = DCOS ( T + CONA / TAN( T ) )
C                                        ** RECURSION RELATION FOR
C                                        ** LEGENDRE POLYNOMIALS
10       PM2 = 1.D0
         PM1 = X
         DO 20 NN = 2, M
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
            PM2 = PM1
            PM1 = P
20       CONTINUE
C
         TMP   = 1.D0 / ( 1.D0 - X**2 )
         PPR   = EN * ( PM2 - X * P ) * TMP
         P2PRI = ( 2.D0 * X * PPR - NNP1 * P ) * TMP
         XI    = X - ( P / PPR ) * ( 1.D0 +
     $               ( P / PPR ) * P2PRI / ( 2.D0 * PPR ) )
C
C                                              ** CHECK FOR CONVERGENCE
         IF ( DABS(XI-X) .GT. TOL ) THEN
            X = XI
            GO TO 10
         END IF
C                          ** iteration FINISHED--CALC. WEIGHTS,
C                          ** ABSCISSAE FOR (-1,1)
         GMU( K ) = - X
         GWT( K ) = 2.D0 / ( TMP * ( EN * PM2 )**2 )
         GMU( NP1 - K ) = - GMU( K )
         GWT( NP1 - K ) =   GWT( K )
30    CONTINUE
C                                    ** SET MIDDLE ABSCISSA AND WEIGHT
C                                    ** FOR RULES OF ODD ORDER
      IF ( MOD( M,2 ) .NE. 0 )  THEN
         GMU( LIM + 1 ) = 0.D0
         PROD = 1.D0
         DO 40 K = 3, M, 2
            PROD = PROD * K / ( K-1 )
40       CONTINUE
         GWT( LIM + 1 ) = 2.D0 / PROD**2
      END IF
C                                        ** CONVERT FROM (-1,1) TO (0,1)
      DO 50  K = 1, M
         GMU( K ) = 0.5D0 * GMU( K ) + 0.5D0
         GWT( K ) = 0.5D0 * GWT( K )
50    CONTINUE
C
      RETURN
      END
      
      real*8 FUNCTION  dpsratio( A, B )
C
C        CALCULATE dpsratio  A/B  WITH OVER- AND UNDER-FLOW PROTECTION
C
      implicit double precision (A-H, O-Z)
      real*8 a, b
*
         IF ( DABS(A).LT.1.0D-8 .AND. DABS(B).LT.1.0D-8 )  THEN
            dpsratio = 1.0D0
         ELSE IF ( B.EQ.0.0 )  THEN
            dpsratio = 1.D+20
         ELSE
            dpsratio = A / B
         END IF
C
      RETURN
      END
*
      SUBROUTINE dpsetbco( ch, chtau, cmu, delm0, fbeam,
     $     gl, lc, maz, mxcmu, nstr, 
     $     taucpr,  xba, xb0, xb1,  
     $     zj, ylmc, ylm0)
C
C       SET COEFFICIENTS IN KS(7) FOR BEAM SOURCE
C
C       INPUT :  CMU         COMPUTATIONAL POLAR ANGLES
C                DELMO       KRONECKER DELTA, DELTA-SUB-M0
C                GL          PHASE FUNCTION LEGENDRE COEFFICIENTS MULTI-
C                              PLIED BY (2L+1) AND SINGLE-SCATTER ALBEDO
C                MAZ         ORDER OF AZIMUTHAL COMPONENT
C                TAUCPR      DELTA-M-SCALED OPTICAL DEPTH
C                YLMC        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                            AT THE QUADRATURE ANGLES -CMU-
C                YLM0        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                            AT THE BEAM ANGLE
C
C       OUTPUT:  xba         ALFA IN EQ. KS(7) 
C                xb0         X-SUB-ZERO IN KS(7)
C                xb1         X-SUB-ONE IN KS(7)
C
C
      implicit double precision (A-H, O-Z)
      DOUBLE PRECISION d1mach
      real*8 big, biggest, ch(*), chtau(0:*), cmu(*),
     $     delm0, fbeam, gl(0:mxcmu,*), q0, q2, 
     $     taucpr(0:*),  xba(*), xb0(mxcmu,*),
     $     xb1(mxcmu,*), zj(*), ylmc(0:mxcmu,*), ylm0(0:*)
      real*8 deltat, pi, q0a, 
     $     q2a, sum


      pi = 2.D0 * DASIN(1.0D0)
      biggest  = d1mach(2)  
      big   = dsqrt(biggest) / 1.D+10
*
*     Calculate x-sub-zero in STWJ(6d) 
*
      iqsav = 1
      DO 100 IQ = 1, NSTR
         SUM = 0.D0
         DO 110 K = MAZ, NSTR-1
            SUM = SUM + GL(K,LC)*YLMC(K,IQ)*YLM0(K)
 110     CONTINUE
         ZJ(IQ) = (2.D0-DELM0)*FBEAM*SUM / (4.D0*PI)
         IF ( zj(iq) .GT. zj(iqsav) ) iqsav = iq
 100  CONTINUE
*     
      q0a = DEXP( -chtau(lc-1) )
      q2a = DEXP( -chtau(lc) )
*     
*     Calculate alfa coefficient for one fixed angle
*     
      q0 = q0a*zj(iqsav)
      q2 = q2a*zj(iqsav)
  
      deltat = taucpr(lc) - taucpr(lc-1)
     
      xba(lc) = 1.D0/ch(lc)
        
      IF(DABS(xba(lc)).GT.big.AND.taucpr(lc).GT.1.D0)THEN
         xba(lc) = 0.0D0
      ENDIF
caky 20110915 Factor 0.25 added after problem with optically very thin cloud at 
caky 20110915 9-10 km that made xba too large and thus xb0 and xb1 too large
caky 20110915 and the solution completely rubbish
caky     IF( DABS(xba(lc)*taucpr(lc)) .GT. DLOG(big))THEN
      IF( DABS(xba(lc)*taucpr(lc)) .GT. 0.25*DLOG(big))THEN
	  xba(lc) = 0.0D0
      ENDIF	  
*
*     Dither alfa if it is close to one of the quadrature angles
*
      IF (  DABS(xba(lc)) .gt. 0.00001d0 ) THEN
         DO IQ = 1, NSTR/2
            IF (DABS((DABS(xba(LC))-1.D0/CMU(IQ))/xba(lc) ) .LT. 
     $           0.05D0 )then
C                  0.05                  001             
                xba(LC) = xba(LC) * 1.001D0

	    ENDIF	
         ENDDO
      ENDIF
      DO iq = 1, nstr
         q0 = q0a * zj(iq)
         q2 = q2a * zj(iq)
*     
*     x-sub-zero and x-sub-one in Eqs. KS(48-49)
*     
       
         xb1(iq,lc) = (1.D0/deltat)*(q2*DEXP(xba(lc)*taucpr(lc))
     $                     -q0*DEXP(xba(lc)*taucpr(lc-1)))
    
	                                       
         xb0(iq,lc) = q0 * DEXP(xba(lc)*taucpr(lc-1))-
     $        xb1(iq,lc)*taucpr(lc-1)
        
      ENDDO
     
*     
      RETURN
      END
*
      SUBROUTINE dpintbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $     mxcmu, mxumu, nstr, numu, taucpr, zb0u, zb1u, xba, 
     $     zju, ylm0, ylmu)
*
      implicit double precision (A-H, O-Z)
      real*8  chtau(0:*), delm0, fbeam,
     $     gl(0:mxcmu,*), q0, q2, taucpr(0:*), xba(*), 
     $     zb0u(mxumu,*), zb1u(mxumu,*),
     $     zju(*), ylm0(0:*), ylmu(0:mxcmu,*)
      real*8 deltat, pi, q0a, q2a, sum
    
*
*     Find coefficients at user angle, necessary for later use in
*     *sterpso*
*
      pi = 2.D0 * DASIN(1.0D0)
*
*     Calculate x-sub-zero in STWJ(6d) 
*
      deltat = taucpr(lc) - taucpr(lc-1)
*
      q0a = DEXP(-chtau(lc-1) )
      q2a = DEXP(-chtau(lc) )
*     
*
      DO 1000 iu = 1, numu 
         sum = 0.0d0
         DO 1100 k = maz, nstr-1
            sum = sum + gl(k,lc)*ylmu(k,iu)*ylm0(k)
 1100    CONTINUE
         zju(iu) = (2.D0-delm0)*fbeam*sum / (4.D0*pi)
 1000 CONTINUE
*     
      DO 5000 iu = 1, numu
         q0 = q0a * zju(iu)
         q2 = q2a * zju(iu)
*     
*     x-sub-zero and x-sub-one in Eqs. KS(48-49)
*

         zb1u(iu,lc)=(1.d0/deltat)*(q2*DEXP(xba(lc)*taucpr(lc))
     $        -q0*DEXP(xba(lc)*taucpr(lc-1)))
         zb0u(iu,lc) = q0 * DEXP(xba(lc)*taucpr(lc-1))-
     $        zb1u(iu,lc)*taucpr(lc-1)
         
      
 5000 CONTINUE
*     
      RETURN
      END
*
      SUBROUTINE  dpssetdis
     $     ( ICHAP,acmu, albedo, bdr, beta, bplank, btemp, ch,
     $     chtau, cmu, cwt, deltam, deltamu, dtauc,EXPbea,
     $     fbeam, flyr, gl, layru, layru_c, lyrcut, maxumu, maxcmu,
     $     ncut, 
     $     newgeo, nlyr, ntau, ntau_c, nn, nstr, noplnk, numu,
     $     onlyfl, oprim, pkag, pkagc, pmom, radius, spher, ssalb,
     $     tauc, taucpr, temis, temper, tplank, ttemp, utau,
     $     utau_c, utaupr, utaupr_c, umu0, umu, usrtau,
     $     usrang,  wvnmlo, 
     $     wvnmhi, zd, numu_c, umu_c, numu_u, umu_u,
     $     BROSZA,SIGBRO,DENSBRO, dtauc_mb,
     $     sza_bro, ndenssza, VN, nrefrac)
*
*     Modified from original disort setdis routine to do various setup
*     operations required by the spherical version of disort.
*
C
C          PERFORM MISCELLANEOUS SETTING-UP operations
C
C       ROUTINES CALLED:  ERRMSG, sqgausn, dpszeroit
C
C       INPUT :  ALL ARE 'DISORT' INPUT VARIABLES (SEE DOC FILE)
C                BROSZA: if =1, BrO density is a function of SZA 
C
C       OUTPUT:  NTAU,UTAU   IF USRTAU = FALSE
C                CMU,CWT     COMPUTATIONAL POLAR ANGLES AND
C                               CORRESPONDING QUADRATURE WEIGHTS
C                EXPBEA      TRANSMISSION OF DIRECT BEAM
C                FLYR        TRUNCATED FRACTION IN DELTA-M METHOD
C                GL          PHASE FUNCTION LEGENDRE COEFFICIENTS MULTI-
C                              PLIED BY (2L+1) AND SINGLE-SCATTER ALBEDO
C                LAYRU       COMPUTATIONAL LAYER IN WHICH -UTAU- FALLS
C                LYRCUT      FLAG AS TO WHETHER RADIATION WILL BE ZEROED
C                              BELOW LAYER -NCUT-
C                NCUT        COMPUTATIONAL LAYER WHERE ABSORPTION
C                              OPTICAL DEPTH FIRST EXCEEDS -ABSCUT-
C                NN          NSTR / 2
C                OPRIM       DELTA-M-SCALED SINGLE-SCATTER ALBEDO
C                TAUCPR      DELTA-M-SCALED OPTICAL DEPTH
C                UMU         OUTPUT ANGLES, EQUAL TO THE 
C                            COMPUTATIONAL POLAR ANGLES, BUT SORTED
C                            IN ASCENDING ORDER
C                UTAUPR      DELTA-M-SCALED VERSION OF -UTAU-
C
C      INTERNAL VARIABLES
C                TEMPC       TEMPERATURE AT CENTER OF LAYER, ASSUMED
C                            TO BE AVERAGE OF LAYER BOUNDARY TEMPERATURES
C
      INTEGER MXCLY, MXULV, MXCMU, MXUMU, MXPHI,
     $        MI, MI9M2, NNLYRI 
      INCLUDE "DISORT.MXD"
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY )

      external chpman2
      LOGICAL  DELTAM, LYRCUT, newgeo, NOPLNK, onlyfl, SPHER, usrang,
     $     USRTAU
      INTEGER  layru(*),  layru_c(*), ntau_c, BROSZA
      real*8 acmu(*), albedo, bdr( mi,0:* ), beta(0:*),
     $     bplank, btemp, ch(*), chtau(0:*), cmu(*), cwt(*), dtauc(*),
     $     EXPbea(0:*), fbeam, flyr(*), gl(0:mxcmu,*),
     $     oprim(*), pkag(0:*), pkagc(*), pmom(0:maxcmu,*), radius,
     $     ssalb(*), tauc(0:*), taucpr(0:*), taup, tempc, 
     $     temis, temper(0:*), tplank, ttemp, umu(*), umu0,
     $     utau(*), utau_c(*), utaupr(*),
     $     utaupr_c(*), wvnmhi, wvnmlo,SIGBRO(0:*),
     $     DENSBRO(mxsza,0:*),dtauc_mb(*),sza_bro(mxsza)

      INTEGER numu_c, numu_u
      real*8 umu_c(*), umu_u(*)
      real*8 abstau, deltamu, f, pi
      real*8 dpsplkavg
      real*8 dtaucpr(mxcly)
*      
      real*8 CHP1(mxcly),CHP2(mxcly),vdbr(mxcly)
      real*8 
     $       DTAUBRO_1(mxcly,mxcly),DTAUBRO_2(mxcly,mxcly),
     $       DBRO(mxcly,0:mxcly),
     $       DBRO_2(mxcly,0:mxcly)
*      real*8 
*     $       DTAUBRO_1(nlyr,nlyr),DTAUBRO_2(nlyr,nlyr),
*     $       DBRO(nlyr,0:nlyr),
*     $       DBRO_2(nlyr,0:nlyr)
     
 
*     inputs and outputs of geofac     
 
      integer NFAC_11(mxcly), NFAC_12(mxcly)
      real*8  fac_11(mxcly,mxcly), fac_12(mxcly,mxcly)
      integer NFAC_21(mxcly), NFAC_22(mxcly)
      real*8  fac_21(mxcly,mxcly), fac_22(mxcly,mxcly)
      
      real*8  SZALOC_1(mxcly,0:mxcly), 
     $        SZALOC_2(mxcly,0:mxcly)
      
*
      real*8  z_lay, zd(0:mxcly), zenang, rearth,
     $        VN(0:mxcly)
     
* 
*     save geometric factors (no dependent on wavelength)
      SAVE nfac_11,nfac_12,fac_11,fac_12
      SAVE nfac_21,nfac_22,fac_21,fac_22
*     
      DATA  ABSCUT / 400.D0 /
      rearth=radius
           
      PI = 2.D0 * DASIN( 1.D0 )

      CALL  dpszeroit( BDR, MI*(MI+1) )

    
      NN = NSTR / 2
      DO 20 IQ = 1, NN
         DO 20 JQ = 0, NN
            BDR(IQ,JQ) = ALBEDO
20    CONTINUE
C
*
*     Set usrtau to include every level when doing the spherical
*     correction.
*
      ntau_c = nlyr + 1
      DO  lc = 0, ntau_c-1
         utau_c(lc+1) = tauc(lc)
      ENDDO
*
      IF ( .NOT.USRTAU ) THEN
C                              ** SET OUTPUT LEVELS AT COMPUTATIONAL
C                              ** LAYER BOUNDARIES
         NTAU = NLYR + 1
         DO 30  LC = 0, NTAU-1
            UTAU(LC+1) = TAUC(LC)
30       CONTINUE
      END IF
C
C                        ** APPLY DELTA-M SCALING AND MOVE DESCRIPTION
C                        ** OF COMPUTATIONAL LAYERS TO LOCAL VARIABLES

      CALL  dpszeroit( TAUCPR(0), MXCLY+1 )
      CALL  dpszeroit( EXPBEA(0), MXCLY+1 )
      CALL  dpszeroit( FLYR, MXCLY )
      CALL  dpszeroit( GL, (MXCMU+1)*MXCLY )
      CALL  dpszeroit( OPRIM, MXCLY )
      ABSTAU = 0.0D0
      chtau(0) = 0.0D0
*
      DO  60  LC = 1, NLYR
      
         PMOM(0,LC) = 1.D0
         IF ( ABSTAU.LT.ABSCUT )  NCUT = LC
         ABSTAU = ABSTAU + ( 1.D0 - SSALB(LC) ) * DTAUC(LC)
C     
         IF ( .NOT.DELTAM )  THEN
            OPRIM(LC) = SSALB(LC)
            TAUCPR(LC) = TAUC(LC)
            DO 40  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) * PMOM(K,LC)
 40         CONTINUE
            dtaucpr(lc)= dtauc(lc)
            F = 0.0D0
         ELSE
C                             ** DO DELTA-M TRANSFORMATION
            F = PMOM( NSTR,LC )
            beta(lc)  = beta(lc) * ( 1.D0 - f * ssalb(lc) )
            OPRIM(LC) = SSALB(LC) * ( 1.D0 - F ) /
     $           ( 1.D0 - F * SSALB(LC) )
            TAUCPR(LC) = TAUCPR(LC-1) +
     $           ( 1. - F*SSALB(LC) ) * DTAUC(LC)
            DO 50  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) *
     $              (PMOM(K,LC)-F) / (1.D0-F)
 50         CONTINUE
            dtaucpr(lc)= taucpr(lc) - taucpr(lc-1)
	   
         ENDIF
         FLYR(LC) = F
 60   CONTINUE
 
   
*
      IF (BROSZA.EQ.0) THEN      
*           
         IF ( fbeam .GT. 0.D0 ) THEN
*
            zenang      = DACOS(umu0) * 180.D0/pi
            IF ( spher ) THEN
*	   
               IF(NEWGEO)THEN
*     GEOMETRIC FACTORS	     
                  z_lay = 0.0D0
              
	
                  if(nrefrac.lt.2)then
*     Refraction with fast method (or no refraction at all)	      
                     CALL geofast
     $                    (brosza,sza_bro, fac_11,nfac_11,szaloc_1,
     $                    fac_12,nfac_12,szaloc_2,
     $                    z_lay,nlyr,zd,zenang,rearth,VN)
                  else
*     Refraction with slow (more precise) method	      
                     CALL geofacpr
     $                    (brosza,sza_bro, fac_11,nfac_11,
     $                    szaloc_1, fac_12,nfac_12,szaloc_2,
     $                    z_lay,nlyr,zd,zenang,rearth,VN)
                  endif
                  rearth=rearth*1.D-5
               ENDIF
*     
*     CHAPMAN FACTORS 
               call chpman2(nlyr,zenang,dtaucpr,nfac_11,fac_11,
     $              nfac_12,fac_12,chp1, mxcly)
     
               IF(NEWGEO)THEN
*     GEOMETRIC FACTORS	 	     

                  z_lay = 0.5D0

                  if(nrefrac.lt.2)then
                     CALL geofast
     $                    (brosza,sza_bro, fac_21,nfac_21,szaloc_1,
     $                    fac_22,nfac_22,szaloc_2,
     $                    z_lay,nlyr,zd,zenang,rearth,VN)
                  else
*     Refraction with slow (more precise) method	      
                     CALL geofacpr
     $                    (brosza,sza_bro, fac_21,nfac_21,
     $                    szaloc_1, fac_22,nfac_22,szaloc_2,
     $                    z_lay,nlyr,zd,zenang,rearth,VN)
                  endif
               ENDIF
*     CHAPMAN FACTORS	      
               call chpman2(nlyr,zenang,dtauc,nfac_21,fac_21,
     $              nfac_22,fac_22,chp2,mxcly)

               EXPbea( 0 ) = 1.D0
               IF( umu0 .LT. 0.D0 ) THEN	     
                  EXPbea(0) =DEXP(-chp1(1))
               ENDIF
*
               DO lc = 1, ncut
                  taup       = tauc(lc-1) + dtauc(lc)/2.D0
                  chtau(lc)=chp1(lc)  
                  ch(lc)     = taup/chp2(lc)
                  EXPbea(lc) = DEXP(-chp1(lc))
               ENDDO
            ENDIF 
            
            IF( ICHAP .EQ. 0 ) THEN
               
*             ********************************************************	   
*             PSEUDO SPHERICAL model - CHAPMAN = PLANE PARALLEL model!
*             ********************************************************

               DO lc = 1, ncut            
                  ch(lc) = umu0
                  chtau(lc)= taucpr(lc) / umu0 
                  EXPbea(lc) = DEXP( - taucpr(lc) / umu0 )
               ENDDO
            ENDIF
*     
         ENDIF
*     
      ELSE
*        
*     BrO density is function of altitude and local SZA

         zenang      = DACOS(umu0) * 180.D0/pi
	
*     
         IF ( spher ) THEN
*	
            IF(newgeo) then	
               z_lay = 0.0D0
               
               if(nrefrac.eq.2)then
                  CALL geofacpr(brosza,sza_bro, fac_11,nfac_11,
     $                 szaloc_1, fac_12,nfac_12,szaloc_2,
     $                 z_lay,nlyr,zd,zenang,rearth,VN)
               else
                  CALL geofast(brosza,sza_bro, fac_11,nfac_11,
     $                 szaloc_1,fac_12,nfac_12,szaloc_2,
     $                 z_lay,nlyr,zd,zenang,rearth,VN)
               endif
               rearth=rearth*1.D-5
            ENDIF

            CALL dtbro(zd,zenang,sigbro,densbro,sza_bro,
     $           NFAC_11, NFAC_12,
     $           SZALOC_1, SZALOC_2, 
     $           nlyr, ndenssza,
     $           DTAUBRO_1, DTAUBRO_2,DBRO,DBRO_2,
     $           vdbr)
            
            CALL chapman4(NLYR,zenang,
     $           NFAC_11,FAC_11,DTAUC_MB,DTAUBRO_1,
     $           NFAC_12,FAC_12,DTAUBRO_2, mxcly,
     $           CHP1)
            
*	     
            IF(newgeo) then	     
*	     
               z_lay = 0.5D0
	     
               if(nrefrac.eq.2)then
                  CALL geofacpr
     $                 (brosza,sza_bro, fac_21,nfac_21,szaloc_1,
     $                 fac_22,nfac_22,szaloc_2,
     $                 z_lay,nlyr,zd,zenang,rearth,VN)
               else
                  CALL geofast
     $                 (brosza,sza_bro,fac_21,nfac_21,szaloc_1,
     $                 fac_22,nfac_22,szaloc_2,
     $                 z_lay,nlyr,zd,zenang,rearth,VN)
               endif
               
            ENDIF    
               
            CALL dtbro(zd,zenang,sigbro,densbro,sza_bro,
     $           nfac_21, nfac_22,
     $           szaloc_1, szaloc_2, 
     $           nlyr, ndenssza,
     $           dtaubro_1,dtaubro_2,DBRO,dbro_2,
     $           vdbr)
            
            
            CALL chapman4(NLYR,zenang,
     $           nfac_21,fac_21,DTAUC_MB,dtaubro_1,
     $           nfac_22,fac_22,dtaubro_2, mxcly,
     $           CHP2)
            
            
            EXPbea( 0 ) = 1.D0
            IF( umu0 .LT. 0.D0 ) THEN	     
               EXPbea(0) = DEXP(-CHP1(1))
            ENDIF
*     
            
            DO lc = 1, nlyr
               
               taup = taucpr(lc-1) + dtaucpr(lc)*0.5d0                      
               chtau(lc)  = CHP1(lc) 
               ch(lc)     = taup/CHP2(lc)
               EXPbea(lc) = DEXP(-CHP1(lc))
               
            ENDDO 
            
            IF( ICHAP .EQ. 0 ) THEN
               
*     ********************************************************	   
*     PSEUDO SPHERICAL model - CHAPMAN = PLANE PARALLEL model!
*     ********************************************************
               
               DO lc = 1, ncut            
                  ch(lc) = umu0
                  chtau(lc)= taucpr(lc) / umu0 
                  EXPbea(lc) = DEXP( - taucpr(lc) / umu0 )
               ENDDO
            ENDIF 
         ENDIF 
         
      ENDIF      
         
     
*

C                      ** IF NO THERMAL EMISSION, CUT OFF MEDIUM BELOW
C                      ** ABSORPTION OPTICAL DEPTH = ABSCUT ( NOTE THAT
C                      ** DELTA-M TRANSFORMATION LEAVES ABSORPTION
C                      ** OPTICAL DEPTH INVARIANT ).  NOT WORTH THE
C                      ** TROUBLE FOR ONE-LAYER PROBLEMS, THOUGH.
      LYRCUT = .false.
      IF ( ABSTAU.GE.ABSCUT .AND. NOPLNK
     $     .AND. NLYR.GT.1 )  LYRCUT =.TRUE.
      IF ( .NOT.LYRCUT )  NCUT = NLYR
C
C                             ** SET ARRAYS DEFINING LOCATION OF USER
C                             ** OUTPUT LEVELS WITHIN DELTA-M-SCALED
C                             ** COMPUTATIONAL MESH
      DO 90  LU = 1, NTAU
         DO 70 LC = 1, NLYR
            IF ( UTAU(LU).GE.TAUC(LC-1) .AND. UTAU(LU).LE.TAUC(LC) )
     $           GO TO 80
70       CONTINUE
         LC = NLYR

80       UTAUPR(LU) = UTAU(LU)
         IF(DELTAM) UTAUPR(LU) = TAUCPR(LC-1) +
     $        (1.-SSALB(LC)*FLYR(LC)) * (UTAU(LU) - TAUC(LC-1))
         LAYRU(LU) = LC
90    CONTINUE
*
      DO lu = 1, ntau_c
         DO lc = 1, nlyr
            IF ( utau_c(lu).GE.tauc(lc-1) .AND. utau_c(lu).LE.tauc(lc) )
     $           GO TO 88
         ENDDO
         lc = nlyr
*
 88      utaupr_c(lu) = utau_c(lu)
         IF ( deltam ) utaupr_c(lu) = taucpr(lc-1) +
     $        (1.-ssalb(lc)*flyr(lc)) * (utau_c(lu) - tauc(lc-1))
         layru_c(lu) = lc
      ENDDO
*
C                      ** CALCULATE COMPUTATIONAL POLAR ANGLE COSINES
C                      ** AND ASSOCIATED QUADRATURE WEIGHTS FOR GAUSSIAN
C                      ** QUADRATURE ON THE INTERVAL (0,1) (UPWARD)
      CALL  dpsqgausn( NN, CMU, CWT )
C                                  ** DOWNWARD (NEG) ANGLES AND WEIGHTS
      DO 100  IQ = 1, NN
         CMU(IQ+NN) = - CMU(IQ)
         CWT(IQ+NN) =   CWT(IQ)
100   CONTINUE
      
   
*
* Set computational angles at where derivatives are needed.
*
      DO iq = 1, nn
         acmu(iq) = - cmu(nn+1-iq)
      ENDDO
      DO iq = nn+1, nstr
         acmu(iq) = cmu(iq-nn)
      ENDDO
*
      numu_c   =  2*nstr
      deltamu  =  0.001D0
      jq       =  1
*
      DO iq = 1, 2*nstr, 2   ! Next all derivative angles
         umu_c(iq)   = acmu(jq) - deltamu
         umu_c(iq+1) = acmu(jq) + deltamu
         jq         = jq + 1
      ENDDO
*
* Set user angles at where derivatives are needed.
*
      IF ( usrang ) THEN
         numu_u   =  2*numu
         deltamu  =  0.001D0
         ju       =  1
*
         DO iu = 1, numu_u, 2   ! Next all derivative angles
            umu_u(iu)   = umu(ju) - deltamu
            umu_u(iu+1) = umu(ju) + deltamu
            IF ( umu_u(iu) .LT. -1.D0 ) THEN
               umu_u(iu) = -1.D0
               umu_u(iu+1) = -1.D0 + 2.0D0*deltamu
            ENDIF
            IF ( umu_u(iu+1) .GT. 1.D0 ) THEN
               umu_u(iu+1) = 1.D0
               umu_u(iu) = 1.D0 - 2.D0*deltamu
            ENDIF
            ju         = ju + 1
         ENDDO
      ENDIF
*
      IF ( .NOT.usrang .OR. (onlyfl .AND. maxumu.GE.nstr )) THEN
*
*                                   ** SET OUTPUT POLAR ANGLES TO
*                                   ** COMPUTATIONAL POLAR ANGLES
         numu = nstr 
         DO 120  IU = 1, NN
            UMU(IU) = - CMU(NN+1-IU)
 120     CONTINUE
         DO 121  IU = NN+1, NSTR
            UMU(IU) = CMU(IU-NN)
 121     CONTINUE
      ENDIF 
      
C
C
C                                   ** CALCULATE PLANCK FUNCTIONS
      IF ( NOPLNK )  THEN
         BPLANK = 0.0D0
         TPLANK = 0.0D0
         CALL  dpszeroit( PKAG, MXCLY+1 )
         CALL  dpszeroit( PKAGC, MXCLY )
      ELSE
         TPLANK = TEMIS * dpsplkavg( WVNMLO, WVNMHI, TTEMP )
         BPLANK =         dpsplkavg( WVNMLO, WVNMHI, BTEMP )
         DO 180  LEV = 0, NLYR
            PKAG( LEV ) = dpsplkavg( WVNMLO, WVNMHI, TEMPER(LEV) )
180      CONTINUE
         DO 190 LC = 1, NLYR
            TEMPC = 0.5D0 * ( TEMPER(LC-1) + TEMPER(LC) )
            PKAGC( LC ) = dpsplkavg( WVNMLO, WVNMHI, TEMPC )
190      CONTINUE
      END IF
C
      RETURN
      END
*
      SUBROUTINE  dpssetmtx( BDR, CBAND, CMU, CWT, DELM0, GC, KK, 
     $                    LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,
     $                    NN, NSTR, TAUCPR, WK )
C
C        CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS
C        OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-
C        OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE
C        SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES
C
C     ROUTINES CALLED:  dpszeroit
C
C     I N P U T      V A R I A B L E S:
C
C       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0    :  KRONECKER DELTA, DELTA-SUB-M0
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       NN       :  NUMBER OF STREAMS IN A HEMISPHERE (NSTR/2)
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
C                   BY LINPACK SOLUTION ROUTINES
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C
C   I N T E R N A L    V A R I A B L E S:
C
C       IROW     :  POINTS TO ROW IN  -CBAND-
C       JCOL     :  POINTS TO POSITION IN LAYER BLOCK
C       LDA      :  ROW DIMENSION OF -CBAND-
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C       NSHIFT   :  FOR POSITIONING NUMBER OF ROWS IN BAND STORAGE
C       WK       :  TEMPORARY STORAGE FOR 'EXP' EVALUATIONS
C ---------------------------------------------------------------------+
      implicit double precision (A-H, O-Z)
      LOGICAL LYRCUT
      real*8    BDR( MI,0:* ), CBAND( MI9M2,NNLYRI ),
     $     CMU(*), CWT(*), delm0, GC( MXCMU,MXCMU,* ),
     $     KK( MXCMU,* ), TAUCPR( 0:* ), WK(*)
      real*8 EXPa, sum
C
C
      CALL  dpszeroit( CBAND, MI9M2*NNLYRI )
      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
C                         ** USE CONTINUITY CONDITIONS OF EQ. STWJ(17)
C                         ** TO FORM COEFFICIENT MATRIX IN STWJ(20);
C                         ** EMPLOY SCALING TRANSFORMATION STWJ(22)
      DO 30  LC = 1, NCUT
C
         DO 4  IQ = 1, NN
            WK(IQ) = DEXP( KK(IQ,LC) * (TAUCPR(LC) - TAUCPR(LC-1)) )
 4       CONTINUE
C
         JCOL = 0
         DO 10  IQ = 1, NN
            NCOL = NCOL + 1
            IROW = NSHIFT - JCOL
            DO 5  JQ = 1, NSTR
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC)
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC) * WK(IQ)
               IROW = IROW + 1
 5          CONTINUE
            JCOL = JCOL + 1
10       CONTINUE
C
         DO 20  IQ = NN+1, NSTR
            NCOL = NCOL + 1
            IROW = NSHIFT - JCOL
            DO 15  JQ = 1, NSTR
               CBAND(IROW+NSTR,NCOL) =   GC(JQ,IQ,LC) * WK(NSTR+1-IQ)
               CBAND(IROW,     NCOL) = - GC(JQ,IQ,LC)
               IROW = IROW + 1
15          CONTINUE
            JCOL = JCOL + 1
20       CONTINUE
C
30    CONTINUE
C                  ** USE TOP BOUNDARY CONDITION OF STWJ(20A) FOR
C                  ** FIRST LAYER
      JCOL = 0
      DO 40  IQ = 1, NN
         EXPA = DEXP( KK(IQ,1) * TAUCPR(1) )
         IROW = NSHIFT - JCOL + NN
         DO 35  JQ = NN, 1, -1
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1) *EXPA
            IROW = IROW+1
35       CONTINUE
         JCOL = JCOL+1
40    CONTINUE
C
      DO 50  IQ = NN+1, NSTR
         IROW = NSHIFT - JCOL + NN
         DO 45  JQ = NN, 1, -1
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1)
            IROW = IROW+1
45       CONTINUE
         JCOL = JCOL+1
50    CONTINUE
C                           ** USE BOTTOM BOUNDARY CONDITION OF
C                           ** STWJ(20C) FOR LAST LAYER
      NNCOL = NCOL - NSTR
      JCOL  = 0
      DO 70  IQ = 1, NN
         NNCOL = NNCOL + 1
         IROW  = NSHIFT - JCOL + NSTR
C
         DO 60  JQ = NN+1, NSTR
            IF ( LYRCUT .OR. DELM0.EQ.0 ) THEN
C
C                          ** NO AZIMUTHAL-DEPENDENT INTENSITY IF LAM-
C                          ** BERT SURFACE; NO INTENSITY COMPONENT IF
C                          ** TRUNCATED BOTTOM LAYER
C
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT)
            ELSE
               SUM = 0.D0
               DO 55  K = 1, NN
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)
     $                        * GC(NN+1-K,IQ,NCUT)
55             CONTINUE
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) - (1.+DELM0) * SUM
            END IF
C
            IROW = IROW + 1
60       CONTINUE
         JCOL = JCOL + 1
70    CONTINUE
C
      DO 90  IQ = NN+1, NSTR
         NNCOL = NNCOL + 1
         IROW  = NSHIFT - JCOL + NSTR
         EXPA = WK(NSTR+1-IQ)
C
         DO 80  JQ = NN+1, NSTR
C
            IF ( LYRCUT .OR. DELM0.EQ.0) THEN
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) *EXPA
            ELSE
               SUM = 0.0D0
               DO 75  K = 1, NN
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)
     $                        * GC(NN+1-K,IQ,NCUT)
75             CONTINUE
               CBAND(IROW,NNCOL) = ( GC(JQ,IQ,NCUT)
     $                               - (1.D0+DELM0) * SUM ) *EXPA
            END IF
C
            IROW = IROW + 1
80       CONTINUE
         JCOL = JCOL + 1
90    CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  dpSETPCO( LC, PKAG, PKAGC, TAUCPR, XRA, XR0, XR1)
C
C       SET COEFFICIENTS IN KS(4) FOR THERMAL SOURCE
C
C       INPUT :    PKAG        INTEGRATED PLANCK FUNCTION 
C                              AT LAYER BOUNDARIES
C                  PKAGC       INTEGRATED PLANCK FUNCTION 
C                              AT LAYER CENTER
C                  TAUCPR      OPTICAL DEPTH (DELTA-M-SCALED)
C
C       OUTPUT:    XR0         X-SUB-ZERO IN EQ. KS(4) 
C                  XR1         X-SUB-ONE IN EQ. KS(4)
C                  XRA         ALFA IN EQ. KS(4)
C
      implicit double precision (A-H, O-Z)
      DOUBLE PRECISION d1mach
      real*8 arg, deltat, depsil
      real*8 PKAG(0:*), PKAGC(*), TAUCPR(0:*), XRA(*),
     $     XR0(*), XR1(*)
C
      DELTAT = TAUCPR(LC) - TAUCPR(LC-1)
      XR1( LC ) = 0.0D0
      XRA( LC ) = 0.0D0
      DEPSIL = 1000.D0 * d1mach(4) 
      
      
C                                  **   CHECK IF DELTAT .GT. THAN
C                                  **   SMAALEST NUMBER ALLOWED ON 
C                                  **   MACHINE, IF DEPSIL = 0.0 IS
C                                  **   USED, ARITHMETIC OVERFLOW
C                                  **   MAY BE THE RESULT LATER, SEE
C                                  **   EQ. FOR XRA(LC).
      IF ( DELTAT.GT.DEPSIL ) THEN
         
         IF ( PKAG(LC).GT.0.0D0 ) THEN
C                                                  EQ. KS(51)
C
            ARG = (PKAGC(LC)/PKAG(LC))**2.D0-PKAG(LC-1)/PKAG(LC)
            IF ( ARG .LT. 0.0D0 ) ARG = 0.0D0
C
C                                   **  USE PLUS (NEGATIVE) SOLUTION IF      
C                                   **  DE(IN)CREASING TEMPERATURE
C
            IF ( PKAG(LC) .GT. PKAG(LC-1) ) THEN
              XRA(LC) = (2.D0/DELTAT) *DLOG((PKAGC(LC)
     $              /PKAG(LC))+DSQRT(ARG))
            ELSE
              XRA(LC) = (2.D0/DELTAT) *DLOG((PKAGC(LC)
     $              /PKAG(LC))-DSQRT(ARG))
            ENDIF
         ELSEIF ( PKAG(LC).EQ.0.0D0 ) THEN
            IF (PKAGC(LC) .GT. 0.0D0 )
     $            XRA(LC) = (2.D0/DELTAT) *
     $           DLOG(PKAG(LC-1)/(2.D0*PKAGC(LC)))
         ENDIF
C                                                  EQ. KS(52)
C
         IF (DABS(XRA(LC)*TAUCPR(LC)).GT.DLOG(d1mach(2)/100.D0)) 
     $                                   XRA(LC) = 0.0D0
         XR1(LC)=(1.D0/DELTAT)*(PKAG(LC)*DEXP(XRA(LC)*TAUCPR(LC))-
     $        PKAG(LC-1)*DEXP(XRA(LC)*TAUCPR(LC-1)))
      ENDIF
C                                                  EQ. KS(53)
C
      XR0( LC ) = PKAG(LC-1)*DEXP(XRA(LC)*TAUCPR(LC-1))
     $                                   -XR1(LC)*TAUCPR(LC-1)
      RETURN
      END
*
      SUBROUTINE  dpssoleig( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZ,
     $                    MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,
     $                    KK, GC, AAD, WKD, EVECCD, EVALD, UMU0 )
C
C         SOLVES EIGENVALUE/VECTOR PROBLEM NECESSARY TO CONSTRUCT
C         HOMOGENEOUS PART OF DISCRETE ORDINATE SOLUTION; STWJ(8B)
C         ** NOTE ** EIGENVALUE PROBLEM IS DEGENERATE WHEN SINGLE
C                    SCATTERING ALBEDO = 1;  PRESENT WAY OF DOING IT
C                    SEEMS NUMERICALLY MORE STABLE THAN ALTERNATIVE
C                    METHODS THAT WE TRIED
C
C     ROUTINES CALLED:  sasymtx
C
C   I N P U T     V A R I A B L E S:
C
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND single-SCATTER ALBEDO)
C       CMU    :  COMPUTATIONAL POLAR ANGLE COSINES
C       CWT    :  WEIGHTS FOR QUADRATURE OVER POLAR ANGLE COSINE
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NN     :  HALF THE TOTAL NUMBER OF STREAMS
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES -CMU-
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T    V A R I A B L E S:
C
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5); NEEDED IN SS(15&18)
C       EVAL   :  -NN- EIGENVALUES OF EQ. SS(12) ON RETURN FROM 'sasymtx'
C                    BUT THEN SQUARE ROOTS TAKEN
C       EVECC  :  -NN- EIGENVECTORS  (G+) - (G-)  ON RETURN
C                    FROM 'sasymtx' ( COLUMN J CORRESPONDS TO -EVAL(J)- )
C                    BUT THEN  (G+) + (G-)  IS CALCULATED FROM SS(10),
C                    G+  AND  G-  ARE SEPARATED, AND  G+  IS STACKED ON
C                    TOP OF  G-  TO FORM -NSTR- EIGENVECTORS OF SS(7)
C       GC     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVECTORS, BUT
C                    IN AN ORDER CORRESPONDING TO -KK-
C       KK     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVALUES OF SS(7),
C                    BUT RD-ORDERED WITH NEGATIVE VALUES FIRST ( SQUARE
C                    ROOTS OF -EVAL- TAKEN AND NEGATIVES ADDED )
C
C   I N T E R N A L   V A R I A B L E S:
C
C       AMB,APB :  MATRICES (ALPHA-BETA), (ALPHA+BETA) IN REDUCED
C                    EIGENVALUE PROBLEM
C       ARRAY   :  COMPLETE COEFFICIENT MATRIX OF REDUCED EIGENVALUE
C                    PROBLEM: (ALFA+BETA)*(ALFA-BETA)
C       GPPLGM  :  (G+) + (G-) (CF. EQS. SS(10-11))
C       GPMIGM  :  (G+) - (G-) (CF. EQS. SS(10-11))
C       WK      :  SCRATCH ARRAY REQUIRED BY 'sasymtx'
C+---------------------------------------------------------------------+
      implicit double precision (A-H, O-Z)
      real*8    AMB( MI,* ), APB( MI,* ),
     $     ARRAY( MI,* ), CC( MXCMU,* ),
     $        CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ), GC( MXCMU,* ),
     $        GL(0:*), KK(*), WK(*), YLMC( 0:MXCMU,* )
      real*8 alpha, beta, gpmigm, gpplgm, sum, umu0
      DOUBLE PRECISION   EVECCD( MI,* ), EVALD(*), WKD(*), AAD( MI,* )
C
C
C                             ** CALCULATE QUANTITIES IN EQS. SS(5-6)
      DO 40 IQ  = 1, NN
C
         DO 20  JQ = 1, NSTR
            SUM = 0.0D0
            DO 10  L = MAZ, NSTR-1
               SUM = SUM + GL(L) * YLMC(L,IQ) * YLMC(L,JQ)
10          CONTINUE
            CC(IQ,JQ) = 0.5D0 * SUM * CWT(JQ)
20       CONTINUE
C
         DO 30  JQ = 1, NN
C                             ** FILL REMAINDER OF ARRAY USING SYMMETRY
C                             ** RELATIONS  C(-MUI,MUJ) = C(MUI,-MUJ)
C                             ** AND        C(-MUI,-MUJ) = C(MUI,MUJ)
C
            CC(IQ+NN,JQ) = CC(IQ,JQ+NN)
            CC(IQ+NN,JQ+NN) = CC(IQ,JQ)
C                                      ** GET FACTORS OF COEFF. MATRIX
C                                      ** OF REDUCED EIGENVALUE PROBLEM
            ALPHA =   CC(IQ,JQ) / CMU(IQ)
            BETA = CC(IQ,JQ+NN) / CMU(IQ)
            AMB(IQ,JQ) = ALPHA - BETA
            APB(IQ,JQ) = ALPHA + BETA
30       CONTINUE
         AMB(IQ,IQ) = AMB(IQ,IQ) - 1.0D0 / CMU(IQ)
         APB(IQ,IQ) = APB(IQ,IQ) - 1.0D0 / CMU(IQ)
C
40    CONTINUE
C                      ** FINISH CALCULATION OF COEFFICIENT MATRIX OF
C                      ** REDUCED EIGENVALUE PROBLEM:  GET MATRIX
C                      ** PRODUCT (ALFA+BETA)*(ALFA-BETA); SS(12)
      DO 70  IQ = 1, NN
         DO 70  JQ = 1, NN
            SUM = 0.D0
            DO 60  KQ = 1, NN
               SUM = SUM + APB(IQ,KQ) * AMB(KQ,JQ)
60          CONTINUE
            ARRAY(IQ,JQ) = SUM
70    CONTINUE
C                      ** FIND (real*8) EIGENVALUES AND EIGENVECTORS
C
      CALL  dpsasymtx( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WK,
     $              AAD, EVECCD, EVALD, WKD, umu0 )
      
      
C
      IF ( IER.GT.0 )  THEN
         WRITE( *, '(//,A,I4,A)' )  ' dpsasymtx--EIGENVALUE NO. ', IER,
     $     '  DIDNT CONVERGE.  LOWER-NUMBERED EIGENVALUES WRONG.'
         CALL  ERRMSG( 'dpsasymtx--CONVERGENCE PROBLEMS', .TRUE. )
      END IF
C
      DO 75  IQ = 1, NN
         EVAL(IQ) = DSQRT( DABS( EVAL(IQ) ) )
	 
         KK( IQ+NN ) = EVAL(IQ)
C                                             ** ADD NEGATIVE EIGENVALUE
         KK( NN+1-IQ ) = - EVAL(IQ)
75    CONTINUE
C                          ** FIND EIGENVECTORS (G+) + (G-) FROM SS(10)
C                          ** AND STORE TEMPORARILY IN -APB- ARRAY
      DO 90  JQ = 1, NN
         DO 90  IQ = 1, NN
            SUM = 0.D0
            DO 80  KQ = 1,NN
               SUM = SUM + AMB(IQ,KQ) * EVECC(KQ,JQ)
80          CONTINUE
            APB(IQ,JQ) = SUM / EVAL(JQ)	    
90    CONTINUE
C
      DO 100  JQ = 1, NN
         DO 100  IQ = 1, NN
            GPPLGM = APB(IQ,JQ)
            GPMIGM = EVECC(IQ,JQ)
C                                ** RECOVER EIGENVECTORS G+,G- FROM
C                                ** THEIR SUM AND DIFFERENCE; STACK THEM
C                                ** TO GET EIGENVECTORS OF FULL SYSTEM
C                                ** SS(7) (JQ = EIGENVECTOR NUMBER)
C
            EVECC(IQ,      JQ) = 0.5D0 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,   JQ) = 0.5D0 * ( GPPLGM - GPMIGM )
C
C                                ** EIGENVECTORS CORRESPONDING TO
C                                ** NEGATIVE EIGENVALUES (CORRESP. TO
C                                ** REVERSING SIGN OF 'K' IN SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5D0 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5D0 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )
	    
	    
100   CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  dpssolve0( ALBEDO, B, BPLANK, CBAND, CMU, CWT, 
     $     EXPBEA, FBEAM, FISOT, il, IPVT, LL, LYRCUT,
     $     MAZ, MI, MI9M2, MXCMU, NCOL, NCUT, NN, noplnk, NSTR,
     $     NNLYRI, PI, TPLANK, TAUCPR, UMU0, Z, ZBEAM0,
     $     ZBEAM1, ZBEAMA, ZPLK0, ZPLK1, ZPLKA )
C
C        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR GENERAL BOUNDARY
C        CONDITIONS STWJ(17) AND SOLVE SYSTEM OF EQUATIONS OBTAINED
C        FROM THE BOUNDARY CONDITIONS AND THE
C        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS.
C        THERMAL EMISSION CONTRIBUTES ONLY IN AZIMUTHAL INDEPENDENCE.
C
C     ROUTINES CALLED:  SGBCO, SGBSL, dpszeroit
C
C     I N P U T      V A R I A B L E S:
C
C       BPLANK   :  BOTTOM BOUNDARY THERMAL EMISSION
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
C                   BY LINPACK SOLUTION ROUTINES
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       EXPBEA   :  TRANSMISSION OF INCIDENT BEAM,DEXP(-TAUCPR/UMU0)
*       il       :  iteration index
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ      :  ORDER OF AZIMUTHAL COMPONENT
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C       NN       :  ORDER OF DOUBLD-GAUSS QUADRATURE (NSTR/2)
C       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       TPLANK   :  TOP BOUNDARY THERMAL EMISSION
C       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       ZBEAMA   :  BEAM SOURCE COEFFICIENT ALFA IN EQ.KS(6)
C       ZBEAM0   :  BEAM SOURCE VECTORS Y-SUB-ZERO IN EQ.KS(6)
C       ZBEAM1   :  BEAM SOURCE VECTORS Y-SUB-ONE IN EQ.KS(6)
C       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. FS(5)
C       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. FS(5)
C       ZPLKA    :  THERMAL SOURCE COEFFICIENT -ZA-, 
C                   ALFA IN EQ. FS(5)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T     V A R I A B L E S:
C
C       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
C                   *SGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.
C                   SC(12), CONSTANTS OF INTEGsratioN WITHOUT
C                   EXPONENTIAL TERM
C      LL        :  PERMANENT STORAGE FOR -B-, BUT RD-ORDERED
C
C   I N T E R N A L    V A R I A B L E S:
C
C       IPVT     :  INTEGER VECTOR OF PIVOT INDICES
C       IT       :  POINTER FOR POSITION IN  -B-
C       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
C       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-
C       Z        :  SCRATCH ARRAY REQUIRED BY *SGBCO*
C+---------------------------------------------------------------------+
C
      implicit double precision (A-H, O-Z)
      LOGICAL  LYRCUT, noplnk
      INTEGER  IPVT(*)
      real*8 albedo, B(*), bplank,
     $     CBAND( MI9M2,NNLYRI ),
     $     CMU(*), CWT(*),EXPBEA(0:*), fbeam, fisot, LL( MXCMU,* ),
     $     pi, TAUCPR( 0:* ), tplank, umu0, Z(*), ZPLK0( MXCMU,* ),
     $     ZPLK1( MXCMU,* ), ZPLKA( * ), ZBEAM0( MXCMU,* ),
     $     ZBEAM1( MXCMU,* ), ZBEAMA( * ) 
      real*8 rcond, sum
C
C
      CALL  dpszeroit( B, NNLYRI )
      IF ( pi .eq. 0.0D0 ) pi = 2.D0 * DASIN(1.D0)
*
C                             ** CONSTRUCT -B-,  STWJ(20A,C) FOR
C                             ** PARALLEL BEAM + BOTTOM REFLECTION +
C                             ** THERMAL EMISSION AT TOP AND/OR BOTTOM
C
      IF ( MAZ.GT.0 .AND. FBEAM.GT.0.0 )  THEN
C
C                                         ** AZIMUTH-DEPENDENT CASE
C                                         ** (NEVER CALLED IF FBEAM = 0)
C               ** NO AZIMUTHAL-DEPENDENT INTENSITY FOR LAMBERT SURFACE;
C               ** NO INTENSITY COMPONENT FOR TRUNCATED BOTTOM LAYER
C
            DO 10  IQ = 1, NN
C                                                     ** TOP BOUNDARY
               B(IQ) = - ZBEAM0(NN+1-IQ,1)
C                                                     ** BOTTOM BOUNDARY
               B(NCOL-NN+IQ) = -DEXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $              (ZBEAM0(IQ+NN,NCUT)+
     $              ZBEAM1(IQ+NN,NCUT)*TAUCPR(NCUT))
             
10          CONTINUE
C
C                             ** CONTINUITY CONDITION FOR LAYER
C                             ** INTERFACES OF EQ. STWJ(20B)
         IT = NN
         DO 40  LC = 1, NCUT-1
            DO 30  IQ = 1, NSTR
               IT    = IT + 1
               B(IT) =
     $              + DEXP(-ZBEAMA(LC+1)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC+1)+ZBEAM1(IQ,LC+1)*TAUCPR(LC))
     $              -  DEXP(-ZBEAMA(LC)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC)+ZBEAM1(IQ,LC)*TAUCPR(LC))
               
 30         CONTINUE
 40      CONTINUE
C
C                                   ** AZIMUTH-INDEPENDENT CASE
*
      ELSE IF ( FBEAM.EQ.0.0 .AND. il .EQ. 0)   THEN
C
         DO 50 IQ = 1, NN
C                                      ** TOP BOUNDARY
C
            B(IQ) = - ZPLK0(NN+1-IQ,1) + FISOT + TPLANK
 50      CONTINUE
C
         IF ( LYRCUT ) THEN
C                               ** NO INTENSITY COMPONENT FOR TRUNCATED
C                               ** BOTTOM LAYER
            DO 60 IQ = 1, NN
C                                      ** BOTTOM BOUNDARY
C
               B(NCOL-NN+IQ) = 
     $              -  DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $              (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
         
 60            continue
C
         ELSE
C
            DO 80 IQ = 1, NN
C
               SUM = 0.D0
               DO 70 JQ = 1, NN
                  SUM = SUM + CWT(JQ) * CMU(JQ) * ALBEDO
     $                 *(DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $                 (ZPLK0(NN+1-JQ,NCUT)+ZPLK1(NN+1-JQ,NCUT)
     $                 *TAUCPR(NCUT)))
 70            CONTINUE
               B(NCOL-NN+IQ) = 2.*SUM + (1.-ALBEDO) * BPLANK
     $              - DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $              (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
            
 80         CONTINUE
         END IF
C                             ** CONTINUITY CONDITION FOR LAYER
C                             ** INTERFACES, STWJ(20B)
         IT = NN
         DO 100  LC = 1, NCUT-1
            DO 90  IQ = 1, NSTR
               IT    = IT + 1
               B(IT) =
     $              +  DEXP(-ZPLKA(LC+1)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC+1)+ZPLK1(IQ,LC+1)*TAUCPR(LC))
     $              -  DEXP(-ZPLKA(LC)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC)+ZPLK1(IQ,LC)*TAUCPR(LC))
           
 90         CONTINUE
 100     CONTINUE
C     
      ELSE IF ( FBEAM .GT. 0.0D0 .OR. il.GT.0 ) THEN
C
         DO 150 IQ = 1, NN
C                                      ** TOP BOUNDARY
C
            B(IQ) = - ZBEAM0(NN+1-IQ,1)
     $           - ZPLK0(NN+1-IQ,1) +FISOT +TPLANK
         
 150     CONTINUE
C
         IF ( LYRCUT ) THEN
C
            DO 160 IQ = 1, NN
C                                      ** BOTTOM BOUNDARY
               B(NCOL-NN+IQ) = - DEXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $              (ZBEAM0(IQ+NN,NCUT)+ZBEAM1(IQ+NN,NCUT)*TAUCPR(NCUT)
     $              )- DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*(ZPLK0(IQ+NN,NCUT
     $              )+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
            
 160        CONTINUE
C
         ELSE
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            DO 180 IQ = 1, NN
C
               SUM = 0.D0
               DO 170 JQ = 1, NN
	          
                  SUM = SUM + CWT(JQ) * CMU(JQ) * ALBEDO
     $                 * ( 
     $                 DEXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $                 (ZBEAM0(NN+1-JQ,NCUT)+ZBEAM1(NN+1-JQ,NCUT)
     $                 *TAUCPR(NCUT))+  DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))
     $                 *(ZPLK0(NN+1-JQ,NCUT)+ZPLK1(NN+1-JQ,NCUT)
     $                 *TAUCPR(NCUT)))
                 
 170           CONTINUE
              
              
               B(NCOL-NN+IQ) = 2.D0*SUM + ( ALBEDO * UMU0*FBEAM/PI
     $              ) *expbea(NCUT)
     $              + (1.D0-ALBEDO) * BPLANK
     $              -  DEXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $              (ZBEAM0(IQ+NN,NCUT)+ZBEAM1(IQ+NN,NCUT)*TAUCPR(NCUT)
     $              )-  DEXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*(ZPLK0(IQ+NN
     $              ,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
             
	       
  
 180        CONTINUE
         END IF
C
         IT = NN
C                             ** CONTINUITY CONDITION FOR LAYER
C                             ** INTERFACES, STWJ(20B)
         DO 200  LC = 1, NCUT-1
            DO 190  IQ = 1, NSTR
               IT    = IT + 1
               B(IT) =
     $              +  DEXP(-ZBEAMA(LC+1)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC+1)+ZBEAM1(IQ,LC+1)*TAUCPR(LC))
     $              -  DEXP(-ZBEAMA(LC)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC)+ZBEAM1(IQ,LC)*TAUCPR(LC))
     $              +  DEXP(-ZPLKA(LC+1)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC+1)+ZPLK1(IQ,LC+1)*TAUCPR(LC))
     $              -  DEXP(-ZPLKA(LC)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC)+ZPLK1(IQ,LC)*TAUCPR(LC))
         	 
 190        CONTINUE
 200     CONTINUE
C     
      END IF
C
C                     ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
C                     ** OF BAND MATRIX -CBAND- AND TEST IF IT IS NEARLY
C                     ** SINGULAR (NOTE: -CBAND- IS DESTROYED)
C                     ** (-CBAND- IS IN LINPACK PACKED FORMAT)
      RCOND = 0.0D0
      NCD = 3*NN - 1
      CALL  dpSGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )
      IF ( 1.0D0+RCOND .EQ. 1.0D0 )  CALL  ERRMSG
     $   ( 'ssolve0--SGBCO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
C                   ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -CBAND-
C                   ** AND R.H. SIDE(S) -B- AFTER -CBAND- HAS BEEN L-U
C                   ** DECOMPOSED.  SOLUTION IS RETURNED IN -B-.
C
      CALL  dpSGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )
C
C                   ** ZERO -CBAND- (IT MAY CONTAIN 'FOREIGN'
C                   ** ELEMENTS UPON RETURNING FROM LINPACK);
C                   ** NECESSARY TO PREVENT ERRORS
C
      CALL  dpszeroit( CBAND, MI9M2*NNLYRI )
C
      DO 220  LC = 1, NCUT
         IPNT = LC*NSTR - NN
         DO 220  IQ = 1, NN
            LL(NN+1-IQ,LC) = B(IPNT+1-IQ)
            LL(IQ+NN,  LC) = B(IQ+IPNT)
220   CONTINUE

      RETURN
      END
*
      SUBROUTINE  dpsterpev( CWT, EVECC, GL, GU, MAZ, MXCMU, MXUMU,
     $                    NN, NSTR, NUMU, WK, YLMC, YLMU )
C
C         INTERPOLATE EIGENVECTORS TO USER ANGLES; EQ SD(8)
C
      implicit double precision (A-H, O-Z)
      integer l, iq, numu, nstr, nn, mxumu, mxcmu, maz, iu, jq
      real*8  CWT(*), EVECC( MXCMU,* ), GL(0:*),
     $     GU(  MXUMU,* ), WK(*),
     $      YLMC(  0:MXCMU,* ), YLMU(  0:MXCMU,* )
      real*8 sum
C
C
      DO 50  IQ = 1, NSTR
C
         DO 20  L = MAZ, NSTR-1
C                                       ** INNER SUM IN SD(8) TIMES ALL
C                                   ** FACTORS IN OUTER SUM BUT PLM(MU)
            SUM = 0.D0
            DO 10  JQ = 1, NSTR
               SUM = SUM + CWT(JQ) * YLMC(L,JQ) * EVECC(JQ,IQ)      
10          CONTINUE
            WK(L+1) = 0.5D0 * GL(L) * SUM
20       CONTINUE
C                                    ** FINISH OUTER SUM IN SD(8)
C                                    ** AND STORE EIGENVECTORS
         DO 40  IU = 1, NUMU
            SUM = 0.D0
            DO 30  L = MAZ, NSTR-1
               SUM = SUM + WK(L+1) * YLMU(L,IU)
30          CONTINUE
            IF ( IQ.LE.NN )  GU( IU, IQ+NN     ) = SUM
            IF ( IQ.GT.NN )  GU( IU, NSTR+1-IQ ) = SUM
40       CONTINUE
C
50    CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  dpsterpso( cwt, fbeam, gl, maz, mxcmu, noplnk,  
     $     numu, nstr, oprim, ylmc, ylmu, psi0, psi1, 
     $     xr0, xr1, xra, z0, z1, zbs0, zbs1, zbsa,
     $     zb0u, zb1u, zbau, zp0u, zp1u, zpau )   
C
C         INTERPOLATES SOURCE FUNCTIONS TO USER ANGLES
C
C    I N P U T      V A R I A B L E S:
C
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       OPRIM  :  SINGLE SCATTERING ALBEDO
C       XR0    :  EXPANSION OF THERMAL SOURCE FUNCTION
C       XR1    :  EXPANSION OF THERMAL SOURCE FUNCTION EQS.SS(14-16)
C       XR2    :  EXPANSION OF THERMAL SOURCE FUNCTION EQS.SS(14-16)
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES
C       YLMU   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE USER ANGLES
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)
C       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)
C       ZJ     :  SOLUTION VECTOR CAPITAL -Z-SUB-ZERO AFTER SOLVING
C                 EQ. SS(19)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
*
*    O U T P U T     V A R I A B L E S:
*
*       zb0u, zb1u:  Incident-beam source function at user angles,
*        zbau        EXPonential-linear-in-optical-depth approximation.
*       zp0u, zp1u:  Components of a EXPonential-linear-in-optical-
*        zpau        depth-dependent source (approximating the 
*                    Planck emission source), at users angles.
C
C   I N T E R N A L       V A R I A B L E S:
C
C       PSI0,  :  SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)
C        PSI1     WITH Z0 AND Z1 RESPECTIVELY
C+---------------------------------------------------------------------+
      implicit double precision (A-H, O-Z)
      LOGICAL  NOPLNK
      real*8     CWT(*), GL(0:*), fbeam, oprim,
     $     PSI0(*), PSI1(*), xr0, xr1, xra,  
     $     YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,*), Z0(*), Z1(*), 
     $     zbs0(*), zbs1(*), zbsa,
     $     zb0u(*), zb1u(*), zbau(*), zp0u(*), zp1u(*), zpau(*)
      real*8 psum0, psum1, sum0, sum1
C     
C
      IF ( FBEAM.GT.0.0 )  THEN
C                                  ** BEAM SOURCE TERMS; EQ. SD(9)
         DO 20  IQ = MAZ, NSTR-1
            psum0 = 0.D0
            psum1 = 0.D0
            DO 10  JQ = 1, NSTR
               psum0 = psum0 + cwt(jq) * ylmc(iq,jq) * zbs0(jq)
               psum1 = psum1 + cwt(jq) * ylmc(iq,jq) * zbs1(jq)
10          CONTINUE
            psi0(iq+1) = 0.5D0 * gl(iq) * psum0
            psi1(iq+1) = 0.5D0 * gl(iq) * psum1
20       CONTINUE
C
         DO 40  IU = 1, numu
            sum0 = 0.D0
            sum1 = 0.D0
            DO 30 IQ = MAZ, NSTR-1
               sum0 = sum0 + ylmu(iq,iu) * psi0(iq+1) 
               sum1 = sum1 + ylmu(iq,iu) * psi1(iq+1)
30          CONTINUE
            zb0u(iu) = sum0 + zb0u(iu)
            zb1u(iu) = sum1 + zb1u(iu)
            zbau(iu) = zbsa 
40       CONTINUE
      END IF
C
      IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 )  THEN
C
C                                ** THERMAL SOURCE TERMS, STWJ(27C)
         DO 80  IQ = MAZ, NSTR-1
            PSUM0 = 0.D0
            PSUM1 = 0.D0
            DO 70  JQ = 1, NSTR
               PSUM0 = PSUM0 + CWT(JQ) * YLMC(IQ,JQ) * Z0(JQ)
               PSUM1 = PSUM1 + CWT(JQ) * YLMC(IQ,JQ) * Z1(JQ)
 70         CONTINUE
            PSI0(IQ+1) = 0.5D0 * GL(IQ) * PSUM0
            PSI1(IQ+1) = 0.5D0 * GL(IQ) * PSUM1
 80       CONTINUE
C
          DO 100  IU = 1, numu
            SUM0 = 0.0D0
            SUM1 = 0.0D0
            DO 90   IQ = MAZ, NSTR-1
               SUM0 = SUM0 + YLMU(IQ,IU) * PSI0(IQ+1)
               SUM1 = SUM1 + YLMU(IQ,IU) * PSI1(IQ+1)
90          CONTINUE
            zp0u(iu) = sum0 + (1.D0-oprim) * xr0
            zp1u(iu) = sum1 + (1.D0-oprim) * xr1
            zpau(iu) = xra
100      CONTINUE
C
      END IF
C
      RETURN
      END
      SUBROUTINE  dpsupbeam( array, cc, cmu, ipvt, mxcmu, nn, nstr, 
     $                    wk, xb0, xb1, xba, zbs0, zbs1, zbsa, zbeam0, 
     $                    zbeam1, zbeama )
*
*       Finds the particular solution of beam source KS(10-11)
*
*     Routines called:  sgeco, sgesl
*
*   I N P U T     V A R I A B L E S:
*
*       cc     :  capital-c-sub-ij in Eq. SS(5)
*       cmu    :  abscissae for gauss quadrature over angle cosine
*       xb0    :  EXPansion of beam source function Eq. KS(7)
*       xb1    :  EXPansion of beam source function Eq. KS(7)
*       xba    :  EXPansion of beam source function Eq. KS(7)
*       (remainder are 'disort' input variables)
*
*    O U T P U T    V A R I A B L E S:
*
*       zbs0     :  solution vectors z-sub-zero of Eq. KS(10-11)
*       zbs1     :  solution vectors z-sub-one  of Eq. KS(10-11)
*       zbsa     :  alfa coefficient in Eq. KS(7)
*       zbeam0, :  permanent storage for -zbs0,zbs1,zbsa-, but rD-ordered
*        zbeam1,
*        zbeama
* 
*   I N T E R N A L    V A R I A B L E S:
*
*       array  :  coefficient matrix in left-hand side of Eq. KS(10)
*       ipvt   :  integer vector of pivot indices required by *linpack*
*       wk     :  scratch array required by *linpack*
*+---------------------------------------------------------------------+
*
      implicit double precision (A-H, O-Z)
      INTEGER ipvt(*)
      real*8    ARRAY( MXCMU,* ), CC( MXCMU,* ),
     $     CMU(*), WK(*), XB0(*), XB1(*), XBA, ZBS0(*), ZBS1(*), ZBSA, 
     $     ZBEAM0(*), ZBEAM1(*), ZBEAMA
      real*8 RCOND, RMIN
*
*
      CALL dpszeroit( ARRAY, MXCMU*MXCMU )
*
      DO IQ = 1, NSTR
         DO JQ = 1, NSTR
            ARRAY(IQ,JQ) = - CC(IQ,JQ)
         ENDDO
      
         ARRAY(IQ,IQ) = 1.D0 + XBA*CMU(IQ) + ARRAY(IQ,IQ)
	
         ZBSA     = XBA
         ZBS1(IQ) = XB1(IQ)
      ENDDO

     
*     SOLVE LINEAR EQUATIONS: 
*
      RCOND = 0.0D0
     
      
      CALL  dpSGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )
     
      RMIN = 1.0D-4

      IF ( RCOND .LT. RMIN ) THEN
*
*     Dither alfa if rcond to small
*               
         if(xba.eq.0.D0)then
* This is not working, do not how to fix. AK 21.01.2009
            xba=0.000000005D0  
         endif
         xba = xba * (1.00000005D0)

	 DO iq = 1, nstr   
	    DO jq = 1, nstr
	       array(iq,jq) = - cc(iq,jq)
	    ENDDO
	    array(iq,iq) = 1.D0 + xba*cmu(iq) + array(iq,iq)
	    zbsa     = xba
	    zbs1(iq) = xb1(iq)
         ENDDO
	
*
*     Solve linear equations KS(10-11)
*
         rcond = 0.0D0

         CALL  dpsgeco( array, mxcmu, nstr, ipvt, rcond, wk )
	  
      ENDIF

*
      IF ( 1.0D0+rcond .EQ. 1.0D0 )  CALL  errmsg
     $   ( ' supbeam--sgeco says matrix near singular',.FALSE.)
*
            
      CALL  dpsgesl( array, mxcmu, nstr, ipvt, zbs1, 0 )
          
      DO iq = 1, nstr
         zbs0(iq) = xb0(iq) + cmu(iq) * zbs1(iq)
      ENDDO
      CALL  dpsgesl( array, mxcmu, nstr, ipvt, zbs0, 0 )
*
*     ... and now some index gymnastic for the inventive ones...
*
      zbeama            = zbsa
                  
      DO iq = 1, nn
         zbeam0( iq+nn )   = zbs0( iq )
         zbeam1( iq+nn )   = zbs1( iq )
         zbeam0( nn+1-iq ) = zbs0( iq+nn )
         zbeam1( nn+1-iq ) = zbs1( iq+nn )
      ENDDO
*
      RETURN
      END
      SUBROUTINE dpsupisot( ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR, 
     $     OPRIM,
     $                    WK, XR0, XR1, XRA, Z0, Z1, ZA, ZPLK0, ZPLK1, 
     $                    ZPLKA )
C
C       FINDS THE PARTICULAR SOLUTION OF THERMAL RADIATION OF SS(15)
C
C     ROUTINES CALLED:  SGECO, SGESL
C
C   I N P U T     V A R I A B L E S:
C
C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       OPRIM  :  DELTA-M SCALED SINGLE SCATTERING ALBEDO
C       XR0    :  EXPANSION OF THERMAL SOURCE FUNCTION EQ. FS(1)
C       XR1    :  EXPANSION OF THERMAL SOURCE FUNCTION EQ. FS(1)
C       XRA    :  EXPANSION OF THERMAL SOURCE FUNCTION EQ. FS(1)
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C    O U T P U T    V A R I A B L E S:
C
C       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. FS(5)
C       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. FS(5)
C       ZA     :  ALFA COEFFICIENT IN FS(5)
C       ZPLK0, :  PERMANENT STORAGE FOR -Z0,Z1,ZA-, BUT RD-ORDERED
C        ZPLK1,
C        ZPLKA
C
C   I N T E R N A L    V A R I A B L E S:
C
C       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. FS(5)
C       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
C       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
C+---------------------------------------------------------------------+
C
      implicit double precision (A-H, O-Z)
      INTEGER IPVT(*)
      real*8    ARRAY( MXCMU,* ), CC( MXCMU,* ),
     $     CMU(*), oprim, WK(*), xr0, xr1, xra, 
     $     Z0(*), Z1(*), ZA, ZPLK0(*), ZPLK1(*), ZPLKA
      real*8 rcond
C
C
      DO 20 IQ = 1, NSTR
C
         DO 10 JQ = 1, NSTR
            ARRAY(IQ,JQ) = - CC(IQ,JQ)
10       CONTINUE
         ARRAY(IQ,IQ) = 1.0D0 + XRA*CMU(IQ) + ARRAY(IQ,IQ)
C
         ZA     = XRA
         Z1(IQ) = ( 1.D0 - OPRIM ) * XR1
20    CONTINUE
C                       ** SOLVE LINEAR EQUATIONS: SAME AS IN *supbeam*,
C                       ** EXCEPT -ZJ- REPLACED BY -Z0-
      RCOND = 0.0D0
      CALL  dpSGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )
      IF ( 1.0D0+RCOND .EQ. 1.D0 )  CALL  ERRMSG
     $   ( 'supisot--SGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
      CALL  dpSGESL( ARRAY, MXCMU, NSTR, IPVT, Z1, 0 )
C
      DO 30 IQ = 1, NSTR
         Z0(IQ) = (1.D0-OPRIM) * XR0 + CMU(IQ) * Z1(IQ)
30    CONTINUE
      CALL  dpSGESL( ARRAY, MXCMU, NSTR, IPVT, Z0, 0 )
C
      DO 40  IQ = 1, NN
         ZPLK0( IQ+NN )   = Z0( IQ )
         ZPLK1( IQ+NN )   = Z1( IQ )
         ZPLKA            = ZA
         ZPLK0( NN+1-IQ ) = Z0( IQ+NN )
         ZPLK1( NN+1-IQ ) = Z1( IQ+NN )
         ZPLKA            = ZA
40    CONTINUE
C
      RETURN
      END
*
      SUBROUTINE dpsusrint( albedo, BPLANK, CMU, CWT, DELM0,EXPBEA,
     $     FBEAM, FISOT, GC, GU, il, KK, LAMBER, LAYRU, LL,
     $     LYRCUT, MAZ, MXCMU, MXULV, MXUMU, NCUT,
     $     NLYR, NN, NSTR, NOPLNK, NUMU, NTAU, PI, 
     $     TAUCPR, TPLANK, UMU, UMU0, UTAUPR, WK,
     $     zb0u, zb1u, zbau, zbeam0, zbeam1, zbeama,
     $     zp0u, zp1u, zpau, ZPLK0, ZPLK1, zplka, uum )
C
C       COMPUTES INTENSITY COMPONENTS AT USER OUTPUT ANGLES
C       FOR AZIMUTHAL EXPANSION TERMS IN EQ. SD(2)
C
C   I N P U T    V A R I A B L E S:
C
C       BPLANK :  INTEGRATED PLANCK FUNCTION FOR EMISSION FROM
C                 BOTTOM BOUNDARY
C       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT    :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
C       emu    :  SURFACE DIRECTIONAL EMISSIVITY
C      expbea :  TRANSMISSION OF INCIDENT BEAM,DEXP(-TAUCPR/UMU0)
C       GC     :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
C                 (i.e., g IN EQ. SC(1) )
C       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU  :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL     :  CONSTANTS OF INTEGsratioN IN EQ. SC(1), OBTAINED
C                 BY SOLVING SCALED VERSION OF EQ. SC(5);
C                 exponentIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NCUT   :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       NN     :  ORDER OF DOUBLD-GAUSS QUADRATURE (NSTR/2)
C       rmu    :  SURFACE BIDIRECTIONAL REFLECTIVITY
C       TAUCPR :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
C       TPLANK :  INTEGRATED PLANCK FUNCTION FOR EMISSION FROM
C                 TOP BOUNDARY
C       UTAUPR :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
C                    COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
C       zp0u   :  Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C       zp1u   :  Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C       zp2u   :  Z-SUB-TWO IN EQ. SS(16) INTERPOLATED TO USER
C                 ANGLES FROM AN EQUATION DERIVED FROM SS(16)
C       ZPLK0  :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
C       ZPLK1  :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
C       zplka  :  THERMAL SOURCE VECTORS -Z2 , BY SOLVING EQ. SS(16)
C       ZBEAM  :  INCIDENT-BEAM SOURCE VECTORS
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)
C
C   O U T P U T    V A R I A B L E S:
C
C       UUM  :  AZIMUTHAL COMPONENTS OF THE INTENSITY IN EQ. STWJ(5)
C
C   I N T E R N A L    V A R I A B L E S:
C
C       BNDDIR :  DIRECT INTENSITY DOWN AT THE BOTTOM BOUNDARY
C       BNDDFU :  DIFFUSE INTENSITY DOWN AT THE BOTTOM BOUNDARY
C       BNDINT :  INTENSITY ATTENUATED AT BOTH BOUNDARIES, STWJ(25-6)
C       DTAU   :  OPTICAL DEPTH OF A COMPUTATIONAL LAYER
C       LYREND :  END LAYER OF INTEGRATION
C       LYRSTR :  START LAYER OF INTEGRATION
C       PALINT :  INTENSITY COMPONENT FROM PARALLEL BEAM
C       PLKINT :  INTENSITY COMPONENT FROM PLANCK SOURCE
C       WK     :  SCRATCH VECTOR FOR SAVING 'EXP' EVALUATIONS
C       ALL THE exponentIAL FACTORS ( exp1, expn,... etc.)
C       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGsratioN IN
C       EQ. SC(12) INTO EQS. S1(8-9).  THEY ALL HAVE NEGATIVE
C       ARGUMENTS SO THERE SHOULD NEVER BE OVERFLOW PROBLEMS.
C+---------------------------------------------------------------------+
C
      implicit double precision (A-H, O-Z)
      LOGICAL  LAMBER, LYRCUT, NOPLNK, NEGUMU
      INTEGER  LAYRU(*)
      real*8     albedo, bplank, CMU(*), CWT(*), delm0, denom, 
     $     dtau, dtau1, dtau2,expbea(0:*), fbeam, fisot,
     $     GC( MXCMU,MXCMU,* ),
     $     GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),
     $     pi, TAUCPR( 0:* ), tplank, uum( mxumu,mxulv, 0:* ),
     $     UMU(*), umu0, UTAUPR(*), WK(*),
     $     zb0u( mxumu,* ), zb1u( mxumu,* ), zbau( mxumu,* ),
     $     zp0u( mxumu,* ), zp1u( mxumu,* ), zpau( mxumu,* ),
     $     zbeam0( mxcmu,* ), zbeam1( mxcmu,* ), zbeama( * ),
     $     ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), zplka( * ) 
      real*8 alfa, bnddfu, bnddir, bndint, dfuint,
     $     emu, exp1, exp2, expn, palint, plkint, rmu, sgn
C
C
      CALL  dpszeroit( UUM, MXUMU*MXULV*(mxcmu+1) )
C>>>>>>>>>>>>>>
       
    
      emu = 1.d0 - albedo
      rmu = albedo 
C                          ** INCORPORATE CONSTANTS OF INTEGRATION INTO
C                          ** INTERPOLATED EIGENVECTORS
      DO 10  LC = 1, NCUT
         DO  10  IQ = 1, NSTR
            DO 10  IU = 1, NUMU 
               GU(IU,IQ,LC) = GU(IU,IQ,LC)* LL(IQ,LC)
              
10    CONTINUE
      
C                           ** LOOP OVER LEVELS AT WHICH INTENSITIES
C                           ** ARE DESIRED ('USER OUTPUT LEVELS')
      DO 200  LU = 1, NTAU
C
         LYU = LAYRU(LU)
C                              ** LOOP OVER POLAR ANGLES AT WHICH
C                              ** INTENSITIES ARE DESIRED
         DO 100  IU = 1, NUMU
            IF ( LYRCUT .AND. LYU.GT.NCUT )  GO TO 100
            NEGUMU = UMU(IU).LT.0.0D0
            IF( NEGUMU )  THEN
               LYRSTR = 1
               LYREND = LYU - 1
               SGN = - 1.0D0
            ELSE
               LYRSTR = LYU + 1
               LYREND = NCUT
               SGN = 1.0D0
            END IF
C                          ** FOR DOWNWARD INTENSITY, INTEGRATE FROM TOP
C                          ** TO 'LYU-1' IN EQ. S1(8); FOR UPWARD,
C                          ** INTEGRATE FROM BOTTOM TO 'LYU+1' IN S1(9)
            PALINT = 0.0D0
            PLKINT = 0.0D0

            DO 30  LC = LYRSTR, LYREND
C             
               DTAU = TAUCPR(LC) - TAUCPR(LC-1)
C	         
               exp1 = DEXP( (UTAUPR(LU) - TAUCPR(LC-1)) / UMU(IU) )
               exp2 = DEXP( (UTAUPR(LU) - TAUCPR( LC )) / UMU(IU) )
C
                 
               IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 ) THEN
                 denom  =  SGN*1.D0/(zpau(iu,lc)*UMU(IU)+1.D0)
                 PLKINT = PLKINT +
     $                    (zp0u(iu,lc)*denom*
     $                        (DEXP(-zpau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                        -DEXP(-zpau(iu,lc)*TAUCPR(LC)) *EXP2 )
     $                    +zp1u(iu,lc)*denom*
     $                        ( (TAUCPR(LC-1)+SGN*denom*UMU(IU) )
     $                        *DEXP(-zpau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                         -(TAUCPR(LC)+SGN*denom*UMU(IU) )
     $                      *DEXP(-zpau(iu,lc)*TAUCPR(LC))* exp2) )
               ENDIF
C
         
	       IF ( FBEAM.GT.0.0 .OR. il.GT.0 )  THEN
		 denom  =  SGN*1.D0/(zbau(iu,lc)*UMU(IU)+1.D0)
		
		 palint = palint +
     $  		 ( zb0u(iu,lc)*denom*
     $  		     (DEXP(-zbau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $  		     -DEXP(-zbau(iu,lc)*TAUCPR(LC)) *EXP2 ) 
     $  		+zb1u(iu,lc)*denom*
     $  		    ( (TAUCPR(LC-1)+SGN*denom*UMU(IU) )
     $  		   *DEXP(-zbau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $  		     -(TAUCPR(LC)+SGN*denom*UMU(IU) )
     $  		   *DEXP(-zbau(iu,lc)*TAUCPR(LC))* exp2) )
		 
               ENDIF
               
C                                                   ** -KK- IS NEGATIVE
 1715          DO 20  IQ = 1, NN
                  WK(IQ) = DEXP( KK(IQ,LC) * DTAU )
		  
                  DENOM = 1.d0 + UMU(IU) * KK(IQ,LC)

                  IF ( DABS(DENOM).LT.0.0001D0 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     expn = DTAU / UMU(IU) * exp2
		    
                  ELSE
                     expn = SGN * ( exp1 * WK(IQ) - exp2 ) / DENOM

                  END IF

                  PALINT = PALINT + GU(IU,IQ,LC) * expn
                  
20             CONTINUE
               
C  ** -KK- IS POSITIVE
                 
               DO 21  IQ = NN+1, NSTR
                  DENOM = 1.0D0 + UMU(IU) * KK(IQ,LC)
		 
                  IF (DABS(DENOM).LT.0.0001D0) THEN
C                                                   ** L'HOSPITAL LIMIT
                     expn = - DTAU / UMU(IU) * exp1
		    
		   
                  ELSE
                     expn = SGN *( exp1 - exp2 * WK(NSTR+1-IQ) ) / DENOM
                  END IF
		 
                  PALINT = PALINT + GU(IU,IQ,LC) * expn
		 
21             CONTINUE
C
30          CONTINUE
C                           ** CALCULATE CONTRIBUTION FROM USER
C                           ** OUTPUT LEVEL TO NEXT COMPUTATIONAL LEVEL
C
            DTAU1 = UTAUPR(LU) - TAUCPR(LYU-1)
            DTAU2 = UTAUPR(LU) - TAUCPR(LYU)
            IF(DABS(DTAU1).LT.1.D-6 .AND. NEGUMU ) then
	       GO TO 50
	    endif
            IF(DABS(DTAU2).LT.1.D-6 .AND. (.NOT.NEGUMU) )then
	       GO TO 50
	    endif
            IF( NEGUMU ) then
	       exp1 = DEXP( DTAU1 / UMU(IU) )
	     
	    endif   
            IF( .NOT.NEGUMU) then	   
	       exp2 = DEXP( DTAU2 / UMU(IU) )
	    endif   
	   
C
            IF ( FBEAM.GT.0.0 .OR. il.GT.0 )  THEN
              IF ( NEGUMU ) THEN
                 expn = exp1
                 ALFA = zbau(iu,lyu)
                 denom = (-1.D0/(ALFA*UMU(IU)+1.D0))
	        
                 palint = palint + 
     $                    zb0u(iu,lyu)*denom*
     $  		       (-DEXP(-ALFA*UTAUPR(LU))
     $  			+ expn*DEXP(-alfa*taucpr(lyu-1)) )
     $  	       +  zb1u(iu,lyu)*denom*
     $  		      ( -(UTAUPR(LU)-UMU(IU)*denom)*
     $  			DEXP(-ALFA*UTAUPR(LU))
     $  			+(TAUCPR(LYU-1)-UMU(IU)*denom) 
     $  			*EXPN*DEXP(-alfa*taucpr(lyu-1)) )
	        
	      ELSE
		 expn = exp2
		 ALFA = zbau(iu,lyu)
		 denom = (1.D0/(ALFA*UMU(IU)+1.D0))
		 palint = palint + 
     $  		  zb0u(iu,lyu)*denom*
     $  		      (DEXP(-ALFA*UTAUPR(LU))
     $  		       -DEXP(-ALFA*TAUCPR(LYU))*EXPN)
     $  	       +  zb1u(iu,lyu)*denom*
     $  		      ( (UTAUPR(LU) +UMU(IU)*denom)
     $  			 *DEXP(-ALFA*UTAUPR(LU))
     $  		       -(TAUCPR(LYU)+UMU(IU)*denom)
     $  			 *DEXP(-ALFA*TAUCPR(LYU)) *EXPN )
	          
	      ENDIF
            ENDIF
	    
C                                                   ** -KK- IS NEGATIVE
            DTAU = TAUCPR(LYU) - TAUCPR(LYU-1)
            DO 40  IQ = 1, NN
               DENOM = 1.D0 + UMU(IU) * KK(IQ,LYU)
               IF (DABS(DENOM).LT.0.0001D0 ) THEN
                  expn = - DTAU2 / UMU(IU) * exp2
		 
               ELSE IF ( NEGUMU ) THEN
                  expn = ( DEXP( - KK(IQ,LYU) * DTAU2 ) -
     $                   DEXP( KK(IQ,LYU) * DTAU ) * exp1 ) / DENOM
               ELSE
                  expn = (DEXP( - KK(IQ,LYU) * DTAU2 ) - exp2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * expn
	       
40          CONTINUE
C                                                   ** -KK- IS POSITIVE
	    
            DO 41  IQ = NN+1, NSTR
               DENOM = 1.D0 + UMU(IU) * KK(IQ,LYU)
               IF (DABS(DENOM).LT.0.0001D0 ) THEN
                  expn = - DTAU1 / UMU(IU) * exp1
		 
               ELSE IF ( NEGUMU ) THEN
                  expn = (DEXP(- KK(IQ,LYU) * DTAU1 ) - exp1 ) / DENOM
               ELSE
                  expn = (DEXP( - KK(IQ,LYU) * DTAU1 ) -
     $                    DEXP( - KK(IQ,LYU) * DTAU ) * exp2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * expn
	  
41          CONTINUE
C
    
            IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 )  THEN
              IF ( NEGUMU ) THEN
                 expn = exp1
                 ALFA = zpau(iu,lyu)
                 denom = (-1.D0/(ALFA*UMU(IU)+1.D0))
                 PLKINT = PLKINT + 
     $                    zp0u(iu,lyu)*denom*
     $                        (-DEXP(-ALFA*UTAUPR(LU))
     $                        + expn*DEXP(-alfa*taucpr(lyu-1)) )
     $                 +  zp1u(iu,lyu)*denom*
     $                        ( -(UTAUPR(LU)-UMU(IU)*denom)*
     $                            DEXP(-ALFA*UTAUPR(LU))
     $                          +(TAUCPR(LYU-1)-UMU(IU)*denom) 
     $                           *EXPN*DEXP(-alfa*taucpr(lyu-1)) )
              ELSE
                 expn = exp2
                 ALFA = zpau(iu,lyu)
                 denom = (1.D0/(ALFA*UMU(IU)+1.D0))
                 PLKINT = PLKINT + 
     $                    zp0u(iu,lyu)*denom*
     $                        (DEXP(-ALFA*UTAUPR(LU))
     $                          -DEXP(-ALFA*TAUCPR(LYU))*EXPN)
     $                 +  zp1u(iu,lyu)*denom*
     $                        ( (UTAUPR(LU) +UMU(IU)*denom)
     $                           *DEXP(-ALFA*UTAUPR(LU))
     $                         -(TAUCPR(LYU)+UMU(IU)*denom)
     $                          *DEXP(-ALFA*TAUCPR(LYU)) *EXPN )
              END IF
            END IF
C                            ** CALCULATE INTENSITY COMPONENTS
C                            ** ATTENUATED AT BOTH BOUNDARIES.
C                            ** NOTE:: NO AZIMUTHAL INTENSITY
C                            ** COMPONENT FOR ISOTROPIC SURFACE
50          BNDINT = 0.0D0
            IF ( il .EQ. 0 ) THEN
               IF ( NEGUMU .AND. MAZ.EQ.0 ) THEN
                  BNDINT = ( FISOT + TPLANK ) * DEXP( UTAUPR(LU) /
     $                 UMU(IU) )
               ELSE IF ( .NOT.NEGUMU ) THEN
                  IF ( LYRCUT .OR. (LAMBER .AND. MAZ.GT.0) )  GO TO 90
                  DO 60  JQ = NN+1, NSTR
                     WK(JQ) = DEXP(-KK(JQ,NLYR)*
     $                    (TAUCPR(NLYR)-TAUCPR(NLYR-1)))
 60               CONTINUE
                  BNDDFU = 0.D0
                  DO 80  IQ = NN, 1, -1
                     DFUINT = 0.D0
                     DO 70  JQ = 1, NN
                        DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
 70                  CONTINUE
                     DO 71  JQ = NN+1, NSTR
                        DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
     $                       * WK(JQ)
 71                  CONTINUE
                     IF ( FBEAM.GT.0.D0 )
     $                    DFUINT = DFUINT + 
     $                    DEXP(-zbeama(nlyr)*taucpr(nlyr)) *
     $                    (zbeam0(iq,nlyr)+
     $                    zbeam1(iq,nlyr)*taucpr(nlyr))
                     DFUINT = DFUINT + DELM0 * ( 
     $                    DEXP(-zplka(nlyr)*TAUCPR(NLYR))*
     $                    (ZPLK0(IQ,NLYR)+
     $                    ZPLK1(IQ,NLYR)*TAUCPR(NLYR)))
                     BNDDFU = BNDDFU + ( 1.D0 + DELM0 ) * rmu
     $                    * CMU(NN+1-IQ) * CWT(NN+1-IQ) * DFUINT
 80               CONTINUE
C     
                  BNDDIR = 0.0D0
                  IF (FBEAM.GT.0.0D0  .OR. umu0.GT.0 ) BNDDIR = UMU0
     $                 * FBEAM / PI * rmu *expbea(NLYR)
                  BNDINT = ( BNDDFU + BNDDIR + DELM0 * emu * BPLANK )
     $                 *DEXP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )
                 
               END IF
            END IF
C
          
90          UUM( IU, LU, maz ) =  PALINT+PLKINT + BNDINT
         
         
                       
	  
C
100     CONTINUE
200   CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  dpszeroit( A, LENGTH )
C
C         ZEROS A real*8 ARRAY -A- HAVING -LENGTH- ELEMENTS
C
      implicit double precision (A-H, O-Z)
      real*8  A(*)
C
      DO 10  L = 1, LENGTH
         A( L ) = 0.D0
10    CONTINUE
C
      RETURN
      END
      real*8 FUNCTION  dpsplkavg ( WNUMLO, WNUMHI, T )
C
C        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS,
*        except if wnunlo .EQ. wnmuhi, then the Planck function at 
*        wnumlo is returned

C
C  NOTE ** CHANGE 'd1mach' TO 'D1MACH' TO RUN IN DOUBLE PRECISION
C
C  I N P U T :  WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL
C                           INTERVAL
C               WNUMHI : UPPER WAVENUMBER
C               T      : TEMPERATURE (K)
C
C  O U T P U T :  dpsplkavg : INTEGRATED PLANCK FUNCTION ( WATTS/SQ M )
C                           = INTEGRAL (WNUMLO TO WNUMHI) OF
C                              2H C**2  NU**3 / (DEXP(HC NU/KT) - 1)
C                              (WHERE H=PLANCKS CONSTANT, C=SPEED OF
C                              LIGHT, NU=WAVENUMBER, T=TEMPERATURE,
C                              AND K = BOLTZMANN CONSTANT)
C
C  REFERENCE : SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
C                 OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.,
C                 JAN. 1974
C
C  METHOD :  FOR  -WNUMLO-  CLOSE TO  -WNUMHI-, A SIMPSON-RULE
C            QUADRATURE IS DONE TO AVOID ILL-CONDITIONING; OTHERWISE
C
C            (1)  FOR WAVENUMBER (WNUMLO OR WNUMHI) SMALL,
C                 INTEGRAL(0 TO WNUM) IS CALCULATED BY EXPANDING
C                 THE INTEGRAND IN A POWER SERIES AND INTEGRATING
C                 TERM BY TERM;
C
C            (2)  OTHERWISE, INTEGRAL(WNUMLO/HI TO INFINITY) IS
C                 CALCULATED BY expandING THE DENOMINATOR OF THE
C                 INTEGRAND IN POWERS OF THE EXPONENTIAL AND
C                 INTEGRATING TERM BY TERM.
C
C  ACCURACY :  AT LEAST 6 SIGNIFICANT DIGITS, ASSUMING THE
C              PHYSICAL CONSTANTS ARE INFINITELY ACCURATE
C
C  ERRORS WHICH ARE NOT TRAPPED:
C
C      * POWER OR EXPONENTIAL SERIES MAY UNDERFLOW, GIVING NO
C        SIGNIFICANT DIGITS.  THIS MAY OR MAY NOT BE OF CONCERN,
C        DEPENDING ON THE APPLICATION.
C
C      * SIMPSON-RULE SPECIAL CASE IS SKIPPED WHEN DENOMINATOR OF
C        INTEGRAND WILL CAUSE OVERFLOW.  IN THAT CASE THE NORMAL
C        PROCEDURE IS USED, WHICH MAY BE INACCURATE IF THE
C        WAVENUMBER LIMITS (WNUMLO, WNUMHI) ARE CLOSE TOGETHER.
C ----------------------------------------------------------------------
C        
                               
      implicit double precision (A-H, O-Z)
      real*8     T, WNUMLO, WNUMHI
C                                   *** LOCAL VARIABLES
C
C        A1,2,... :  POWER SERIES COEFFICIENTS
C        c1       :  First radiation constant ( 2 * h * c**2 )
C        C2       :  H * C / K, IN UNITS CM*K (H = PLANCKS CONSTANT,
C                      C = SPEED OF LIGHT, K = BOLTZMANN CONSTANT)
C        D(I)     :  exponentIAL SERIESDEXPANSION OF INTEGRAL OF
C                       PLANCK FUNCTION FROM WNUMLO (I=1) OR WNUMHI
C                       (I=2) TO INFINITY
C        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 ON
C                       COMPUTER
C        EX       : DEXP( - V(I) )
C        EXM      :  EX**M
C        MMAX     :  NO. OF TERMS TO TAKE IN exponentIAL SERIES
C        MV       :  MULTIPLES OF 'V(I)'
C        P(I)     :  POWER SERIESDEXPANSION OF INTEGRAL OF
C                       PLANCK FUNCTION FROM ZERO TO WNUMLO (I=1) OR
C                       WNUMHI (I=2)
C        PI       :  3.14159...
C        SIGMA    :  STEFAN-BOLTZMANN CONSTANT (W/M**2/K**4)
C        SIGDPI   :  SIGMA / PI
C        SMALLV   :  NUMBER OF TIMES THE POWER SERIES IS USED (0,1,2)
C        V(I)     :  C2 * (WNUMLO(I=1) OR WNUMHI(I=2)) / TEMPERATURE
C        VCUT     :  POWER-SERIES CUTOFF POINT
C        VCP      :  exponentIAL SERIES CUTOFF POINTS
C        VMAX     :  LARGEST ALLOWABLE ARGUMENT OF 'EXP' FUNCTION
C
      real*8 a1, a2, a3, a4, a5, a6
      PARAMETER (A1 = 1.D0/3.D0,A2 = -1.D0/8.D0,A3 = 1.D0/60.D0,
     $     A4 = -1.D0/5040.D0,
     $     A5 = 1.D0/272160.D0, A6 = -1.D0/13305600.D0)
      INTEGER  SMALLV
      real*8 arg, c1, C2, CONC, D(2), del, EPSIL,
     $      EX, exm, hh, MV, oldval, P(2), pi, SIGMA, SIGDPI,
     $     val, val0, V(2), VCUT, VCP(7), vmax, VSQ, wvn, x
      DOUBLE PRECISION   D1MACH
      real*8 F
      SAVE     CONC, VMAX, EPSIL, SIGDPI
      DATA     c1 /1.1911D-8/
      DATA     C2 / 1.438786D0 /,  SIGMA / 5.67032D-8 /,
     $     VCUT / 1.5D0 /,
     $     VCP / 10.25D0, 5.7D0, 3.9D0, 2.9D0, 2.3D0, 1.9D0, 0.0D0/
      DATA     PI / 0.D0 /
      F(X) = X**3 / ( DEXP(X) - 1 )
C
C
      IF ( PI.EQ.0.0D0 )  THEN
         PI = 2.D0 * DASIN( 1.D0 )
         VMAX = DLOG( d1mach(2) )
         EPSIL = d1mach(4)
	
         SIGDPI = SIGMA / PI
         CONC = 15.D0 / PI**4
      END IF
C
      IF( T.LT.0.0D0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. )
     $  CALL ERRMSG( 'dpsplkavg--TEMPERATURE OR WAVENUMS. WRONG',.TRUE.)
C
      IF ( T.LT.1.D-4 )  THEN
         dpsplkavg = 0.0D0
         RETURN
      ENDIF
*
      IF ( wnumhi .eq. wnumlo ) THEN
         wvn  =  wnumhi
         arg  = DEXP( - C2 * wvn / T )
         dpsplkavg = c1 * (wvn**3.D0) * arg / ( 1.D0 - arg )
         RETURN
      ENDIF
*     
C
      V(1) = C2 * WNUMLO / T
      V(2) = C2 * WNUMHI / T
      IF ( V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND.
     $     (WNUMHI-WNUMLO)/WNUMHI .LT. 1.D-2 )  THEN
C
C                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
C                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
         HH = V(2) - V(1)
         OLDVAL = 0.0D0
         VAL0 = F( V(1) ) + F( V(2) )
C
         DO  2  N = 1, 10
            DEL = HH / (2*N)
            VAL = VAL0
            DO  1  K = 1, 2*N-1
               VAL = VAL + 2*(1+MOD(K,2)) * F( V(1) + K*DEL )
    1       CONTINUE
            VAL = DEL/3.D0 * VAL
            IF ( DABS( (VAL-OLDVAL)/VAL ) .LE. 1.D-6 )  GO TO 3
            OLDVAL = VAL
    2    CONTINUE
         CALL ERRMSG( 'dpsplkavg--SIMPSON RULE DIDNT CONVERGE', .FALSE.)
C
    3    dpsplkavg = SIGDPI * T**4 * CONC * VAL
         RETURN
      END IF
C
      SMALLV = 0
      DO  50  I = 1, 2
C
         IF( V(I).LT.VCUT )  THEN
C                                   ** USE POWER SERIES
            SMALLV = SMALLV + 1
            VSQ = V(I)**2
            P(I) =  CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) *
     $                ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) )
         ELSE
C                    ** USE exponentIAL SERIES
            MMAX = 0
C                                ** FIND UPPER LIMIT OF SERIES
   20       MMAX = MMAX + 1
               IF ( V(I).LT.VCP( MMAX ) )  GO TO 20
C
            EX = DEXP( - V(I) )
            EXM = 1.0D0
            D(I) = 0.0D0
C
            DO  30  M = 1, MMAX
               MV = M * V(I)
               EXM = EX * EXM
               D(I) = D(I) +
     $              EXM * ( 6.D0 +
     $              MV*( 6.D0 + MV*( 3.D0 + MV ) ) ) / M**4
   30       CONTINUE
C
            D(I) = CONC * D(I)
         END IF
C
   50 CONTINUE
C
      IF ( SMALLV .EQ. 2 ) THEN
C                                    ** WNUMLO AND WNUMHI BOTH SMALL
         dpsplkavg = P(2) - P(1)
C
      ELSE IF ( SMALLV .EQ. 1 ) THEN
C                                    ** WNUMLO SMALL, WNUMHI LARGE
         dpsplkavg = 1.D0 - P(1) - D(2)
C
      ELSE
C                                    ** WNUMLO AND WNUMHI BOTH LARGE
         dpsplkavg = D(1) - D(2)
C
      END IF
C
      dpsplkavg = SIGDPI * T**4 * dpsplkavg
      IF( dpsplkavg.EQ.0.0 )
     $    CALL ERRMSG( 'dpsplkavg--RETURNS ZERO; POSSIBLE UNDERFLOW',
     $                 .FALSE. )
C
      RETURN
      END
*
      SUBROUTINE dpPARDER( deltamu, maz, mxcmu, mxulv, mxumu, ntau,
     $      nstr, parmu, uu )
*
*     Calculate the partial derivative of the intensity with respect to
*     the cosine of the polar angle, the parital derivative in the
*     first term on the right-hand side of Eq. 13 Dahlback and Stamnes
*     (1991).
*
*     References:
*     Dahlback, A. and K. Stamnes, 1991, 'A new spherical model for
*     computing the radiation field available for photolysis and
*     heating at twilight', Planet. Space Sci, vol. 39, pp. 671-683.
*     
*
      implicit double precision (A-H, O-Z)
      real*8 deltamu, parmu(mxcmu, *),
     $     uu(mxumu, mxulv, 0:*)
*
      jq = 1
      DO iq = 1, 2*nstr, 2
         DO lu = 1, ntau
            parmu(jq, lu)=(uu(iq+1,lu, maz)-uu(iq,lu, maz))/deltamu
         ENDDO
         jq = jq + 1
      ENDDO
*
      RETURN
      END
*
      SUBROUTINE dpsetpar( acmu, beta, lc, mxcmu, nstr, parmu, radius,
     $     taucpr, xba, xb0, xb1, zd, zj, zj0, zj2 )
*
*     Set coefficients in Eq. 7, Kylling and Stamnes (1992) for
*     perturbation source. Alfa is set equal to zero so a linear
*     approximation is used.
*
*     References:
*
*     Dahlback, A. and K. Stamnes, 1991, 'A new spherical model for
*     computing the radiation field available for photolysis and
*     heating at twilight', Planet. Space Sci, vol. 39, pp. 671-683.
*
*     Kylling, A. and K. Stamnes, 1992, 'Efficient yet accurate
*     solution of the linear transport equation in the presence of
*     internal sources: the exponential-linear-in-depth approximation',
*     Journal of Computational Physics, vol. 102, pp. 265-276.
*
      implicit double precision (A-H, O-Z)
      real*8 acmu(*), beta(0:*), parmu(mxcmu, *), q0, q2, 
     $     radius, taucpr(0:*),xba(*), xb0(mxcmu,*), xb1(mxcmu,*), 
     $     zd(0:*), zj(*), zj0(*), zj2(*)
      real*8 betafact0, betafact1, deltat, fact
*
      fact  = 1.D+05         ! Convert from km to cm
      nn    = nstr/2
*
*     Calculate the internal source, second term Eq. 11c,
*     Dahlback and Stamnes (1991)
*
      betafact0 = beta(lc-1)*(zd(lc-1)+radius)*fact
      betafact1 = beta(lc)*(zd(lc)+radius)*fact
      DO iq = 1, nstr
*
*     For the very top level beta might be zero, avoid problems
*
         IF (beta(lc-1) .EQ. 0.D0 ) THEN
            zj0(iq) = 0.0D0 
         ELSE
            zj0(iq)  =  -((1.D0-acmu(iq)*acmu(iq)) / betafact0)
     $           * parmu( iq, lc )
         ENDIF
         zj2(iq)  = -((1.D0-acmu(iq)*acmu(iq)) / betafact1)
     $        * parmu( iq, lc+1 )
*
      ENDDO      
*
*     Rearrange everything so that internal source agree with Disort
*     weird cmu angles, yac.....
*
      DO iq = 1, nn
         zj(iq)         = zj0(nn+iq)
         zj(nn+iq)      = zj0(nn+1-iq)
      ENDDO
      DO iq = 1, nstr
         zj0(iq)        = zj(iq)
      ENDDO
      DO iq = 1, nn
         zj(iq)         = zj2(nn+iq)
         zj(nn+iq)      = zj2(nn+1-iq)
      ENDDO
      DO iq = 1, nstr
         zj2(iq)        = zj(iq)
      ENDDO
*
*     Calculate coefficients in Eq. 7, Kylling and Stamnes (1992),
*     actually same as coefficents in Eq. 12, Dahlback and Stamnes
*     (1991) since the linear approximation is used.
*
      deltat = taucpr(lc) - taucpr(lc-1)
      xba(lc) = 0.D0
      DO iq = 1, nstr
         q0 = zj0(iq)
         q2 = zj2(iq)
*     
*     X-sub-zero and X-sub-one in Eqs. 48-49, Kylling and Stamnes
*     (1992) for alfa=0, i.e. linear approximation.
*     
         xb1(IQ,LC) = 0.0D0
         IF ( deltat .GT. 1.0D-09 )then
	   xb1(iq,lc) =(1.D0/deltat)*(q2-q0)
	 ENDIF
         xb0(iq,lc) = q0 - xb1(iq,lc)*taucpr(lc-1)
      ENDDO
*
      RETURN
*
      END
