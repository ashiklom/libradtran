      SUBROUTINE  spsdisort( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $     wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, umu, 
     $     nphi, phi, fbeam, beta, nil, umu0, phi0, newgeo, zd, spher, 
     $     radius, fisot, albedo, btemp, ttemp, temis, deltam, planck, 
     $     onlyfl, accur, quiet, ierror, prnt, header,
     $     maxcly, maxulv, maxumu, maxcmu, maxphi, 
     $     rfldir, rfldn, flup, dfdt, uavg, uu, u0u,
     $     uavgdn, uavgso, uavgup)
*+-----------------------------------------------------------------------+
* 04.01.2001: Added old asymtx, lepoly and soleig to this file. 
*             Disort 2.0 subs not compatible with old ones. AKY.
*
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
*     If run in spherical geometry the code should be compiled and
*     linked in double precision. I have not tried the spherical
*     version in single precision and do not know if it does what it
*     should do. If you wonder how to build a double precision version
*     of this routine, see the file Makefile_double that comes with the
*     uvspec distribution.
*
*    ( See disort.doc and twostr.doc for further documentation )
*
*     Bug reports and other feedback to arve.kylling@nilu.no
*
************************************************************************
C
C+---------------------------------------------------------------------+
C------------------    I/O VARIABLE SPECIFICATIONS     -----------------
C+---------------------------------------------------------------------+
C
      CHARACTER  HEADER*127
      LOGICAL  DELTAM, newgeo, quiet, planck, onlyfl, PRNT(7), 
     $     SPHER, usrang, USRTAU
      INTEGER  ierror(33), MAXCLY, MAXULV, MAXCMU, maxphi, maxumu, 
     $     NLYR, nphi, NSTR, NTAU, numu
      REAL     accur, ALBEDO, beta(0:maxcly), BTEMP,
     $     DTAUC( MAXCLY ), FBEAM, FISOT, phi0, 
     $     PHI( MAXPHI ), PMOM( 0:MAXCMU, MAXCLY ), RADIUS,
     $     SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     $     WVNMLO, WVNMHI, UMU0, UTAU( MAXULV ), 
     $     U0U(MAXUMU, MAXULV), UU( MAXUMU, MAXULV, MAXPHI ),
     $     ZD( 0:MAXCLY )
C     
      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ),
     $     FLUP( MAXULV ),
     $     UAVG( MAXULV ), DFDT( MAXULV ), umu( maxumu ),
     $     UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )
C
C+---------------------------------------------------------------------+
C      ROUTINES CALLED (IN ORDER):  SLFTST, ZEROAL, CHEKIN, SETDIS,
C                                   PRTINP, LEPOLY, SOLEIG, SETBCO,
C                                   UPBEAM, SETPCO, UPISOT, terpev,
C                                   terpso, SETMTX, SOLVE0, FLUXES,
*                                   USRINT
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
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'SOLEIG')
C
C   APB(IQ/2,IQ/2)    SECOND MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
C                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN 'SOLEIG')
C
C   ARRAY(IQ,IQ)      SCRATCH MATRIX FOR 'SOLEIG', 'UPBEAM' AND 'UPISOT'
C                     (SEE EACH SUBROUTINE FOR DEFINITION)
C
C   B()               RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
C                     *SOLVE0*;  RETURNS AS SOLUTION VECTOR
C                     VECTOR CAPITAL-L, THE CONSTANTS OF INTEGRATION
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
C                     *SOLEIG* ; STORED PERMANENTLY IN -GC-
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
C                     LITTLE-G  IN EQ. SC(1)
C
*   gl(iu,iq,lc)      Eigenvectors interpolated to user polar angles
*                     ( little-g in Eqs. SC(3) and S1(8-9), i.e.
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
C   LL(IQ,LC)         CONSTANTS OF INTEGRATION CAPITAL-L IN EQ. SC(1),
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
C   Z()               SCRATCH ARRAY USED IN * SOLVE0* TO SOLVE A
C                     LINEAR SYSTEM FOR THE CONSTANTS OF INTEGRATION
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
c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
C+---------------------------------------------------------------------+
      
      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
      INCLUDE 'DISORT.MXD'
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY, MXIL = 20 )
*
      DOUBLE PRECISION d1mach
      REAL ratio
      LOGICAL LYRCUT, noplnk
      INTEGER layru( mxulv )
      REAL acmu(mxcmu), azerr, azterm, bdr( mi,0:mi ),
     $     bplank, ch( mxcly ), chtau( 0:2*mxcly+1 ),
     $     cosphi, cmu( mxcmu ), cwt( mxcmu ), delm0,
     $     deltamu, epsil, expbea( 0:mxcly ),
     $     flyr( mxcly ), fldn( mxulv ), 
     $     fldir( mxulv ), gl( 0:mxcmu,mxcly ), oprim( mxcly ), 
     $     PHIRAD( MXPHI ), pi, pkag( 0:mxcly ), pkagc(mxcly), rpd, 
     $     tauc( 0:mxcly ), taucpr( 0:mxcly ), 
     $     tplank,
     $     u0uls(mxumu, mxulv), uum( mxumu, mxulv, 0:mxcmu), 
     $     utaupr( mxulv )
*
*     Variables for the internal source at computational angles
*
      INTEGER ntau_c, numu_c, layru_c( mxulv )
      REAL  utau_c(mxcly+1),
     $     umu_c( 2*mxcmu ),  utaupr_c( mxulv )
*
*     Variables for the internal source at user angles
*
      INTEGER numu_u
      REAL  umu_u( 2*mxcmu )
*
      DATA  nerr / 33 /
*
*
      pi = 2. * ASIN(1.0)
      epsil = 10.*r1mach(4)
      rpd = pi / 180.0
      noplnk = .NOT. planck
*
*     Zero some arrays (not strictly necessary, but otherwise unused
*     parts of arrays collect garbage)
*
      DO i = 1, nerr
         ierror(i) = 0
      ENDDO
      CALL  ZEROIT( CH    , MXCLY )
      CALL  ZEROIT( CMU   , MXCMU )
      CALL  ZEROIT( CWT   , MXCMU )
*     
*     Calculate cumulative optical depth and dither single-scatter
*     albedo to improve numerical behavior of eigenvalue/vector
*     computation
*
      TAUC( 0 ) = 0.
      CALL  ZEROIT( TAUC(0), MXCLY+1 )
      DO 20  LC = 1, NLYR
         epsil = 0.001
         IF( SSALB(LC) .GT. 1.0-epsil )  SSALB(LC) = 1.0 - epsil
         TAUC(LC) = TAUC(LC-1) + DTAUC(LC)
 20   CONTINUE
*
*     Check input dimensions and variables
*
      CALL  schekin( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
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
      CALL  ssetdis( acmu, albedo, bdr, beta, bplank, btemp, ch, 
     $     chtau, cmu, cwt, deltam, deltamu, dtauc, expbea, fbeam, 
     $     flyr, gl, layru, layru_c, lyrcut, maxumu, maxcmu,
     $     ncut, newgeo, 
     $     nlyr, ntau, ntau_c, nn, nstr, noplnk, numu, onlyfl,
     $     oprim, pkag, pkagc, pmom, radius, spher, ssalb, tauc,
     $     taucpr, temis, temper, tplank, ttemp, utau, utau_c,
     $     utaupr, utaupr_c, umu0, umu, usrtau, usrang,
     $     wvnmlo, wvnmhi, 
     $     zd, numu_c, umu_c, numu_u, umu_u)
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
      IF ( prnt(1) )
     $     CALL sprtinp( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $     wvnmhi, ntau, utau, nstr, numu, umu, fbeam, umu0,
     $     spher, nil, fisot, albedo, btemp, ttemp, temis, deltam,
     $     noplnk, onlyfl, flyr, lyrcut, oprim, tauc, taucpr, 
     $     maxcmu, prnt(7) )
C
C ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
C ========  (EQ STWJ 5)
C
      KCONV = 0
      NAZ = NSTR-1
*
*                                            ** AZIMUTH-INDEPENDENT CASE
*
      IF ( FBEAM.EQ.0.0 .OR. (1.-UMU0).LT.1.E-5 .OR. ONLYFL .OR.
     $     (NUMU.EQ.1.AND.(1.-UMU(1)).LT.1.E-5 ) )
     $     NAZ = 0
*
      CALL  ZEROIT( UU, MAXUMU*MAXULV*MAXPHI )
      DO MAZ = 0, NAZ
*      DO MAZ = 0, 0
*
         IF ( MAZ.EQ.0 )  DELM0 = 1.0
         IF ( MAZ.GT.0 )  DELM0 = 0.0
*
*     
*
         CALL foucom( accur, acmu, albedo, bdr, beta, bplank, btemp,
     $        ch, chtau, cmu, cwt, delm0, deltam, deltamu, dtauc, 
     $        expbea, fbeam, fisot, flyr, gl, layru, 
     $        layru_c, lyrcut, maxcly, maxcmu, maxulv, maxumu, maz, 
     $        ncut, newgeo, nil, nlyr, noplnk, nn, nstr, ntau, ntau_c, 
     $        numu, numu_c, numu_u, onlyfl, oprim, pi, pkag, pkagc, 
     $        pmom, prnt, radius, spher, ssalb, tauc,  
     $        taucpr, temis, temper, tplank, ttemp, umu0, umu, 
     $        umu_c, umu_u, usrang, utau, utau_c, utaupr, 
     $        utaupr_c, utaupr_u, wvnmhi, wvnmlo, zd,
     $        dfdt, flup, fldn, fldir, rfldir, rfldn, uavg, u0u,
     $        uum, u0uls, uavgdn, uavgso, uavgup )
C
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
      AZERR = 0.0
      DO  J = 1, NPHI
         COSPHI = COS( MAZ * PHIRAD(J) )
         DO  LU = 1, NTAU
            DO IU = 1, NUMU
               AZTERM = UUM( IU,LU,MAZ ) * COSPHI
               UU( IU,LU,J ) = UU( IU,LU,J ) + AZTERM
               AZERR = AMAX1( RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ),
     $                        AZERR )
            ENDDO
         ENDDO
      ENDDO
      IF ( AZERR.LE.ACCUR )  KCONV = KCONV + 1
      IF ( KCONV.GE.2 )      GOTO 210
C
*
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
      SUBROUTINE foucom( accur, acmu, albedo, bdr, beta, bplank,
     $     btemp, ch, chtau, cmu, cwt, delm0, deltam, deltamu, dtauc, 
     $     expbea, fbeam, fisot, flyr, gl, layru, layru_c, 
     $     lyrcut, maxcly, maxcmu, maxulv, maxumu, maz, ncut, newgeo,
     $     nil, nlyr, noplnk, nn, nstr, ntau, ntau_c, numu, numu_c,
     $     numu_u, onlyfl, oprim, pi, pkag, pkagc, pmom, prnt, radius, 
     $     spher, ssalb, tauc, taucpr, temis, temper, 
     $     tplank, ttemp, umu0, umu, umu_c, umu_u, usrang, utau, 
     $     utau_c, utaupr, utaupr_c, utaupr_u, wvnmhi, wvnmlo, 
     $     zd, dfdt, flup, fldn, fldir, rfldir, rfldn,
     $     uavg, u0u, uum, u0uls, uavgdn, uavgso, uavgup)
*
*     Solve Eq. 6a (STWJ) for each FOUrier COMponent
*     
C+---------------------------------------------------------------------+
C   LOCAL SYMBOLIC DIMENSIONS:
C
c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
C+---------------------------------------------------------------------+
      
      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
      INCLUDE 'DISORT.MXD'
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY, MXIL = 20 )
*
      INTEGER ipvt(nnlyri), layru(mxulv)
      LOGICAL lamber, LYRCUT
      LOGICAL  DELTAM, newgeo, NOPLNK, onlyfl, PRNT(7), SPHER, usrang
      REAL     accur, ALBEDO, beta(0:maxcly), 
     $     bplank, BTEMP, delm0, deltamu, DTAUC( MAXCLY ), FBEAM,
     $     FISOT, pi, PMOM( 0:MAXCMU, MAXCLY ), RADIUS,
     $     SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, tplank, TTEMP,
     $     WVNMLO, WVNMHI, UMU0, UTAU( MAXULV ), 
     $     U0U(MAXUMU, MAXULV), ZD( 0:MAXCLY )
      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ),
     $     FLUP( MAXULV ),
     $     UAVG( MAXULV ), DFDT( MAXULV ), umu( maxumu ),
     $     UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )
*
      REAL    ANGCOS(1), acmu(mxcmu), albedos,
     $     amb( mi,mi ), apb( mi,mi ),
     $     array( mxcmu,mxcmu ), b( nnlyri ), bdr( mi,0:mi ), 
     $     cband( mi9m2,nnlyri ), cc( mxcmu,mxcmu ), cerr, cdiff,
     $     ch( mxcly ), chtau(0:*), cmu( mxcmu ), cwt( mxcmu ), 
     $     eval(mi), evecc( mxcmu, mxcmu ), 
     $     expbea( 0:mxcly ), dfdts(mxcly), fbeams, flups(mxcly), 
     $     flyr( mxcly ), fldn( mxulv ), fldir( mxulv ), 
     $     gl( 0:mxcmu,mxcly ), gc(mxcmu,mxcmu,mxcly ), 
     $     gu( mxumu, mxcmu, mxcly ), kk( mxcmu, mxcly ), 
     $     ll( mxcmu,mxcly ), oprim( mxcly ), pkag( 0:mxcly ),
     $     pkagc(mxcly), prev_cerr, psi0( mxcmu ), psi1( mxcmu ),
     $     rfldns(mxcly), rfldirs(mxcly), sgn, tauc( 0:mxcly ),
     $     taucpr( 0:mxcly ), tynt, 
     $     u0c( mxcmu,mxulv ), uavgs( mxulv ),
     $     u0uls(mxumu, mxulv), u0utmp(mxumu, mxulv),
     $     uum( mxumu, mxulv,0:mxcmu), utaupr( mxulv ), wk( mxcmu ),
     $     xr0(mxcly ), xr1( mxcly ), xra(mxcly), ylm0( 0:mxcmu ),
     $     ylmc(0:mxcmu,mxcmu ),ylmu(0:mxcmu, mxumu), z( nnlyri )
      REAL
     $     xba(mxcly ), xb0( mxcmu,mxcly ), xb1( mxcmu,mxcly ),
     $     zbeama(mxcly ), zbeam0( mxcmu,mxcly ), 
     $     zbeam1( mxcmu,mxcly ), zbs0(mxcmu ), zbsa, 
     $     zbs1( mxcmu ), z0( 2*mxcmu ),
     $     z1(2*mxcmu ), za, zb0( 2*mxcmu ),
     $     zb1(2*mxcmu ), zba, zj0(mxcmu ),  zj2( mxcmu ),
     $     zj( mxcmu ),  zju( mxumu), zplk0( mxcmu,mxcly ),
     $     zplk1( mxcmu,mxcly ), zb0u( mxumu, mxulv ),
     $     zb1u( mxumu, mxulv ), zbau( mxumu,mxulv),
     $     zp0u( mxumu, mxulv), zp1u( mxumu, mxulv ),
     $     zpau(mxumu,mxulv)
*
*     Variables for the internal source at computational angles
*
      INTEGER ntau_c, numu_c, layru_c( mxulv )
      REAL gu_c( 2*mxumu, mxcmu, mxcly ),
     $     parmu_c(mxcmu, mxulv), utau_c(mxcly+1),
     $     umu_c( 2*mxcmu ),  utaupr_c( mxulv ),
     $     ylmu_c(0:mxcmu, 2*mxumu),
     $     zb0u_c( 2*mxumu, mxulv ), zb1u_c( 2*mxumu, mxulv ), 
     $     zbau_c( 2*mxumu,mxulv), uum_c( 2*mxumu, mxulv, 0:mxcmu)
*
*     Variables for the internal source at user angles
*
      INTEGER numu_u
      REAL gu_u( 2*mxumu, mxcmu, mxcly ),
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
*     to avoid memory havocking when calling *terpso*
*
      REAL   zp0u_c( 2*mxumu, mxulv), 
     $    zp1u_c( 2*mxumu, mxulv ), zpau_c( 2*mxumu,mxulv)
*
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
      CALL zeroit( dfdts, mxulv )
      CALL zeroit( flups, mxulv )
      CALL zeroit( rfldns, mxulv )
      CALL zeroit( rfldirs, mxulv )
      CALL zeroit( uavgs, mxulv )
      CALL zeroit( u0uls, mxumu*mxulv )
      CALL zeroit( u0utmp, mxumu*mxulv )
      CALL  ZEROIT( AMB, MI**2 )
      CALL  ZEROIT( APB, MI**2 )
      CALL  ZEROIT( ARRAY, MXCMU**2 )
      CALL  ZEROIT( CC,    MXCMU**2 )
      CALL  ZEROIT( EVAL  ,MXCMU )
      CALL  ZEROIT( EVECC, MXCMU**2 )
      CALL  ZEROIT( gc, mxcmu*MXCMU*MXCLY )
      CALL  ZEROIT( GU, mxumu*MXCMU*MXCLY )
      CALL  ZEROIT( KK,     MXCMU*MXCLY )
      CALL  ZEROIT( LL,     MXCMU*MXCLY )
      CALL  ZEROIT( WK    ,  MXCMU )
      CALL  ZEROIT( xba, MXCLY )
      CALL  ZEROIT( xb0, MXCMU*MXCLY )
      CALL  ZEROIT( xb1, MXCMU*MXCLY )
      CALL  ZEROIT( XR0, MXCLY )
      CALL  ZEROIT( XR1, MXCLY )
      CALL  ZEROIT( XRA, MXCLY )
      CALL  ZEROIT( Z, NNLYRI )
      CALL  ZEROIT( Z0    ,  2*MXCMU )
      CALL  ZEROIT( Z1    ,  2*MXCMU )
      CALL  ZEROIT( Zb0    ,  2*MXCMU )
      CALL  ZEROIT( Zb1    ,  2*MXCMU )
      ZA = 0
      ZbA = 0
      CALL  ZEROIT( ZBEAMA, MXCLY )
      CALL  ZEROIT( ZBEAM0, MXCMU*MXCLY )
      CALL  ZEROIT( ZBEAM1, MXCMU*MXCLY )
      CALL  ZEROIT( ZBS0, MXCMU )
      CALL  ZEROIT( ZBS1, MXCMU )
      CALL  ZEROIT( ZJ,  MXCMU )
      CALL  ZEROIT( ZPLK0,  MXCMU*MXCLY )
      CALL  ZEROIT( ZPLK1,  MXCMU*MXCLY )
*
*     Only zero these the first time, otherwise expect a downer
*
      IF ( maz .EQ. 0 ) THEN
         CALL  ZEROIT( YLM0, MXCMU+1 )
         CALL  ZEROIT( YLMC, (MXCMU+1)*MXCMU )
         CALL  ZEROIT( ylmu, (mxcmu+1)*mxumu )
      ENDIF
*
*
*     Get normalized associated legendre polynomials for incident beam
*     angle cosine
*
      IF ( FBEAM .GT. 0.0 ) THEN
         ANGCOS(1) = -UMU0
         CALL  LEPOLYs( 1, MAZ, MXCMU, NSTR-1, ANGCOS, YLM0 )
      ENDIF
*
*     Get normalized associated legendre polynomials for computational
*     polar angle cosines
*
      IF ( .NOT. onlyfl .AND. usrang ) THEN
         CALL lepolys( numu, maz, mxcmu, nstr-1, umu, ylmu )
         CALL lepolys( numu_u, maz, mxcmu, nstr-1, umu_u, ylmu_u )
      ENDIF
      CALL lepolys( numu_c, maz, mxcmu, nstr-1, umu_c, ylmu_c )
      CALL  LEPOLYs( NN,   MAZ, MXCMU, NSTR-1, CMU, YLMC )
*
*     Evaluate normalized associated legendre polynomials with negative
*     -cmu- from those with positive -cmu-; Dave/Armstrong eq. (15)
*
      SGN  = - 1.0
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
      prev_cerr = r1mach(2)
*
* ===== Begin loop over each term in Eq. 10, Dahlback and Stamnes (1991)
*
      DO 99 il = 0, nil
*
*     For all but the pseudo-spherical approximation (il=0) set
*     fbeam=0.0 and albedo=0.0.
*
         IF ( il .GT. 0 ) THEN
            albedo = 0.0
            fbeam  = 0.0
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
            CALL  soleigs( amb, apb, array, cmu, cwt, gl(0,lc), mi, maz,
     $           mxcmu, nn, nstr, wk, ylmc, cc, evecc, eval,
     $           kk(1,lc), gc(1,1,lc) )
*
*     Set coefficients in ks(7)
*
*     For the beam source or for the perturbation source
*
            IF ( FBEAM.GT.0.0 .OR. il.GT.0 ) THEN
*
               IF ( il .EQ. 0 ) THEN
*
*     Set coefficients for the internal source for il=0
*
                  CALL setbco( ch, chtau, cmu, delm0, fbeam,
     $                 gl, lc, maz, mxcmu, nstr, 
     $                 taucpr, xba, xb0, xb1, zj, 
     $                 ylmc, ylm0 )
*
*     Get coefficients at cmu+-deltamu
*                  
                  IF ( nil .GT. 0 ) THEN
                     CALL intbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, 2*mxumu, nstr, numu_c,
     $                    taucpr, zb0u_c, zb1u_c, xba,
     $                    zju, ylm0, ylmu_c )
                  ENDIF
*     
                  IF ( usrang ) THEN
*     
*     Get coefficients at umu
*     
                     CALL intbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, mxumu, nstr, numu,
     $                    taucpr, zb0u, zb1u, xba,
     $                    zju, ylm0, ylmu )
*     
*     Get coefficients at umu+-deltamu
*     
                     CALL intbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $                    mxcmu, 2*mxumu, nstr, numu_u,
     $                    taucpr, zb0u_u, zb1u_u, xba,
     $                    zju, ylm0, ylmu_u )
                  ENDIF
*
               ELSE
*
*     Set coefficients for the internal source for il > 0
*
                  CALL setpar( acmu, beta, lc, mxcmu, nstr, parmu_c,
     $                 radius,taucpr, xba, xb0, xb1, zd, zj, zj0, zj2 )
*
*     Do some index gymnastic to satisfy the fantastic representation
*     of angles in cmu.
*
*     ooooh, bananas.....
*
                  DO iq = 1, nn
                     zb0(2*iq-1)      = xb0(nstr+1-iq,lc)
                     zb0(2*iq)        = xb0(nstr+1-iq,lc)
                     zb0(2*nn+2*iq-1) = xb0(iq,lc)
                     zb0(2*nn+2*iq)   = xb0(iq,lc)
                     zb1(2*iq-1)      = xb1(nstr+1-iq,lc)
                     zb1(2*iq)        = xb1(nstr+1-iq,lc)
                     zb1(2*nn+2*iq-1) = xb1(iq,lc)
                     zb1(2*nn+2*iq)   = xb1(iq,lc)
                  ENDDO
*     
                  DO iu = 1, numu_c
                     zb0u_c(iu,lc) = zb0(iu)
                     zb1u_c(iu,lc) = zb1(iu)
                  ENDDO
*
                  IF ( usrang ) THEN
                     CALL setpar( umu_u, beta, lc, mxcmu, numu_u,
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
               CALL  supbeam( array, cc, cmu, ipvt, mxcmu, nn, nstr, 
     $              wk, xb0(1,lc), xb1(1,lc), xba(lc), 
     $              zbs0, zbs1, zbsa, zbeam0(1,lc), zbeam1(1,lc), 
     $              zbeama(lc) )
cakyp
*               write(*,*)lc, zbeam0(1,lc), zbeam0(1,lc),
*     $              zbeam1(2,lc), zbeam1(2,lc)
cakyp
*
*     and now a nip-up....
*     I got this from Arne Dahlback, numerical problems may happen
*     (even in double precision) if the perturbation source is not
*     ignored for optically thin layers. The value for tynt is from
*     Arne. (and tynt is norwegian and means thin).
*
               tynt = 0.001
               IF ( tauc(lc).LT.tynt .AND. il .GT. 0 ) THEN
                  CALL zeroit(zbs0,mxcmu)
                  CALL zeroit(zbs1,mxcmu)
                  DO iq = 1, nstr
                     zbeam0(iq,lc) = 0.0
                     zbeam1(iq,lc) = 0.0
                  ENDDO
               ENDIF
            ENDIF
C     
            IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 ) THEN
*
*     Set coefficients in ks(7) for thermal source
*
	      CALL SETCOE( KK(1, LC), LC, NCUT, NN, PKAG(LC-1), 
     $             PKAGC(LC), PKAG(LC), TAUCPR, XR0(LC), XR1(LC), 
     $             XRA(LC) )
*
*     Calculate particular solutions of Eq. ks(11) for thermal emission
*     source
*
               CALL upisot( array, cc, cmu, ipvt, mxcmu, nn, nstr,
     $              oprim(lc), wk, xr0(lc), xr1(lc), xra(lc), z0,
     $              z1, zplk0(1,lc), zplk1(1,lc) )
            END IF
*
*     Interpolate eigenvectors and source terms to angles needed to
*     get the internal source, cf. Eq. 11c, DS.
*
            IF ( nil .GT. 0 ) THEN
               CALL terpev( cwt, evecc, gl(0,lc), gu_c(1,1,lc), maz,
     $              mxcmu, 2*mxumu, nn, nstr, numu_c, wk, ylmc, ylmu_c)
            ENDIF
*
            IF (il .GT. 0 )  THEN
               fbeam = 1.0                 ! Just to get beam-source 
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
               CALL  sterpso( cwt, fbeam, gl(0,lc), maz,
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
               CALL  sterpso( cwt, fbeam, gl(0,lc), maz,
     $              mxcmu, noplnk, numu_u, nstr, oprim(lc),
     $              ylmc, ylmu_u, psi0, psi1, 
     $              xr0(lc), xr1(lc), xra(lc), z0, z1,
     $              zbs0, zbs1, zbsa,
     $              zb0u_u(1,lc), zb1u_u(1,lc), zbau_u(1,lc), 
     $              zp0u(1,lc), zp1u(1,lc), zpau(1,lc))   
            ENDIF
*
            IF ( il .GT. 0 ) THEN
               fbeam = 0.0      ! Back to zero as it should be here
            ENDIF
*     
            IF ( .NOT.onlyfl .AND. usrang ) THEN
*
*     Interpolate eigenvectors to user angles
*
               CALL  TERPEV( CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZ,
     $              MXCMU,MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU)
               CALL terpev( cwt, evecc, gl(0,lc), gu_u(1,1,lc), maz,
     $              mxcmu, 2*mxumu, nn, nstr, numu_u, wk, ylmc, ylmu_u)
*
*     Interpolate source terms to user angles
*
               IF (il .gt. 0 )  THEN
                  fbeam = 1.0   ! Just to get beam-source interpolated
               ENDIF
*
*     Standard disort way of interpolating source terms
*     
               CALL  sterpso( cwt, fbeam, gl(0,lc), maz,
     $              mxcmu, noplnk, numu, nstr, oprim(lc),
     $              ylmc, ylmu, psi0, psi1, 
     $              xr0(lc), xr1(lc), xra(lc), z0, z1,
     $              zbs0, zbs1, zbsa,
     $              zb0u(1,lc), zb1u(1,lc), zbau(1,lc), 
     $              zp0u(1,lc), zp1u(1,lc), zpau(1,lc))   
*     
               IF (il .gt. 0 ) THEN
                  fbeam = 0.0   ! Back to zero as it should be here
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
         IF ( il .GT. 0 ) CALL zeroit(bdr,mi*(mi+1))
         CALL  ssetmtx( bdr, cband, cmu, cwt, delm0, gc, kk, 
     $        lyrcut, mi, mi9m2, mxcmu, ncol, ncut, 
     $        nnlyri, nn, nstr, taucpr, wk )
*     
*     Solve for constants of integration in homo-geneous solution 
*     (general boundary conditions)
*
         CALL  ssolve0( albedo, b, bplank, cband, cmu, cwt, 
     $        expbea, fbeam, fisot, il, ipvt, ll, lyrcut,
     $        maz, mi, mi9m2, mxcmu, ncol, ncut, nn, noplnk, nstr,
     $        nnlyri, pi, tplank, taucpr, umu0, z, zbeam0,
     $        zbeam1, zbeama, zplk0, zplk1, xra )
*
*     Compute upward and downward fluxes
*
         IF ( maz .eq. 0 ) THEN
            CALL sfluxes( ch, cmu, cwt, fbeam, gc, il, kk, 
     $           layru, ll, lyrcut, mxcmu, mxulv, ncut, nn, 
     $           nstr, ntau, pi, prnt, ssalb, taucpr, umu0,
     $           utau, utaupr, xr0, xr1, xra, zplk0,
     $           zplk1, xra, zbeama, zbeam0, zbeam1, dfdt, 
     $           flup, fldn, fldir, rfldir, rfldn, uavg,
     $           u0c, maxulv, uavgdn, uavgso, uavgup)
*     
         ENDIF
*
*     Compute intensity components at
*     angles needed to calculate the derivative term in Eq. 11c
*     at computational angles.
*
         IF ( nil .GT. 0 ) THEN
            CALL  susrint( albedo, bplank, cmu, cwt, delm0, expbea,
     $           fbeam, fisot, gc, gu_c, il, kk, lamber, layru_c, ll,
     $           lyrcut, maz, mxcmu, mxulv, 2*mxumu, ncut,
     $           nlyr, nn, nstr, noplnk, numu_c, ntau_c, pi, 
     $           taucpr, tplank, umu_c, umu0, utaupr_c, wk,
     $           zb0u_c, zb1u_c, zbau_c, zbeam0, zbeam1, zbeama,
     $           zp0u_c, zp1u_c, zpau_c, zplk0, zplk1, xra, uum_c )
         ENDIF
*     
         IF ( usrang ) THEN
*
*     Compute azimuthal averaged intensity components at
*     angles for the derivative computations.
*
            CALL  susrint( albedo, bplank, cmu, cwt, delm0, expbea,
     $           fbeam, fisot, gc, gu_u, il, kk, lamber, layru_c, ll,
     $           lyrcut, maz, mxcmu, mxulv, 2*mxumu, ncut,
     $           nlyr, nn, nstr, noplnk, numu_u, ntau_c, pi, 
     $           taucpr, tplank, umu_u, umu0, utaupr_c, wk,
     $           zb0u_u, zb1u_u, zbau_u, zbeam0, zbeam1, zbeama,
     $           zp0u_c, zp1u_c, zpau_c, zplk0, zplk1, xra, uum_u )
*
*     Compute azimuthal averaged intensity components at
*     user angles.
*
            CALL  susrint( albedo, bplank, cmu, cwt, delm0, expbea,
     $           fbeam, fisot, gc, gu, il, kk, lamber, layru, ll,
     $           lyrcut, maz, mxcmu, mxulv, mxumu, ncut,
     $           nlyr, nn, nstr, noplnk, numu, ntau, pi, 
     $           taucpr, tplank, umu, umu0, utaupr, wk,
     $           zb0u, zb1u, zbau, zbeam0, zbeam1, zbeama,
     $           zp0u, zp1u, zpau, zplk0, zplk1, xra, uum )
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
         cerr  = 0.0
         DO lu = 1, ntau
            IF ( il .GT. 0 ) THEN
               DO iu = 1, numu
                  IF ( ABS(uum(iu,lu,maz)) .GT. 1.0E-09 ) THEN
                     cdiff = ABS( uum(iu,lu,maz) / u0utmp(iu,lu) )
                  ELSE
                     cdiff = 0.0
                  ENDIF
                  cerr = AMAX1( cdiff, cerr )
               ENDDO
            ENDIF
         ENDDO
*
*     If cerr is not increasing save results from iteration
*
         IF ( cerr .LT. prev_cerr ) THEN
            IF ( maz .eq. 0 ) THEN
               DO lu = 1, ntau
                  dfdts(lu)    =  dfdts(lu) + dfdt(lu)
                  flups(lu)    =  flups(lu) + flup(lu)
                  rfldns(lu)   =  rfldns(lu) + rfldn(lu)
                  rfldirs(lu)  =  rfldirs(lu) + rfldir(lu)
                  uavgs(lu)    =  uavgs(lu) + uavg(lu)
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
         IF ( nil .GT. 0 .AND. fbeam .GT. 0.0 ) THEN
            CALL PARDER( deltamu, maz, mxcmu, mxulv, 2*mxumu, ntau_c,
     $           nstr, parmu_c, uum_c )
         ENDIF
*
         IF ( usrang .AND. fbeam .GT. 0.0 ) THEN
            CALL PARDER( deltamu, maz, mxcmu, mxulv, 2*mxumu, ntau_c,
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
      SUBROUTINE intbco( chtau, delm0,  fbeam, gl, lc, maz, 
     $     mxcmu, mxumu, nstr, numu, taucpr, zb0u, zb1u, xba, 
     $     zju, ylm0, ylmu )
*
      REAL  chtau(0:*), delm0, fbeam,
     $     gl(0:mxcmu,*), q0, q2, taucpr(0:*), xba(*), 
     $     zb0u(mxumu,*), zb1u(mxumu,*),
     $     zju(*), ylm0(0:*), ylmu(0:mxcmu,*)
      REAL deltat, pi, q0a, q2a, sum
*
*     Find coefficients at user angle, necessary for later use in
*     *terpso*
*
      pi = 2. * ASIN(1.0)
*
*     Calculate x-sub-zero in STWJ(6d) 
*
      deltat = taucpr(lc) - taucpr(lc-1)
*
      q0a = EXP(-chtau(lc-1) )
      q2a = EXP(-chtau(lc) )
*     
*
      DO 1000 iu = 1, numu 
         sum = 0.0
         DO 1100 k = maz, nstr-1
            sum = sum + gl(k,lc)*ylmu(k,iu)*ylm0(k)
 1100    CONTINUE
         zju(iu) = (2.-delm0)*fbeam*sum / (4.0*pi)
 1000 CONTINUE
*     
      DO 5000 iu = 1, numu
         q0 = q0a * zju(iu)
         q2 = q2a * zju(iu)
*     
*     x-sub-zero and x-sub-one in Eqs. KS(48-49)
*
         zb1u(iu,lc) = (1./deltat)*(q2*EXP(xba(lc)*taucpr(lc))
     $        -q0*EXP(xba(lc)*taucpr(lc-1)))
         zb0u(iu,lc) = q0 * EXP(xba(lc)*taucpr(lc-1))-
     $        zb1u(iu,lc)*taucpr(lc-1)
 5000 CONTINUE
*     
      RETURN
      END
*
      SUBROUTINE PARDER( deltamu, maz, mxcmu, mxulv, mxumu, ntau,
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
      REAL deltamu, parmu(mxcmu, *),
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
      SUBROUTINE  sCHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     $     WVNMHI, USRTAU, NTAU, UTAU, NSTR, 
     $     usrang, numu, umu, FBEAM, 
     $     UMU0, SPHER, FISOT, ALBEDO, BTEMP,
     $     TTEMP, TEMIS, NOPLNK, onlyfl, ierror, quiet,
     $     MAXCLY, MAXULV, maxumu, MAXCMU, MXCLY, MXULV,  MXCMU, 
     $     mxumu,  TAUC, ZD )
C
C           CHECKS THE INPUT DIMENSIONS AND VARIABLES
C
      LOGICAL  WRTBADQ, WRTDIMQ
      LOGICAL  NOPLNK, onlyfl, SPHER, usrang, USRTAU, INPERR,
     $     quiet
      INTEGER  MAXCLY, MAXULV, MAXCMU, NLYR, NSTR, NTAU, numu,
     $         MXCMU, MXCLY, MXULV, ierror(*)
      REAL     ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
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
         IF ( dtauc(lc).LT.0.0 )  THEN
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
            IF( ABS(utau(lu)-tauc(nlyr)).LE.1.E-4 .OR.
     $             (ABS(UTAU(LU)-TAUC(NLYR)).LE.1.E-6*TAUC(NLYR)))
     $           utau(lu) = tauc(nlyr)
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
     $          umu(iu).EQ. 0.0) THEN
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
      umumin = 0.0
      IF ( spher ) umumin = - 1.0
      IF ( fbeam.GT.0.0 .AND. ( umu0.LE.umumin .OR. umu0.GT.1.0 ) )
     $     THEN
         inperr = wrtbadq( quiet, 'umu0' )
         ierror(13) = 1
      ENDIF
Cgg
      IF ( FBEAM.GT. 0.0 .AND. umu0.EQ. 1.0) 
     $     CALL ERRMSG( 'Does not work for umu0=1.0', .TRUE. )
Cgg
      IF ( fisot.LT.0.0 ) THEN
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
         IF ((.NOT.noplnk .AND. ABS(temper(lc)-temper(lc-1)).GT.50.0) 
     $          .AND. .NOT. quiet )
     $          CALL errmsg( 'chekin--vertical temperature step may'
     $          // ' be too large for good accuracy', .FALSE. )
100   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE setbco( ch, chtau, cmu, delm0, fbeam,
     $     gl, lc, maz, mxcmu, nstr, 
     $     taucpr, xba, xb0, xb1,  
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
      DOUBLE PRECISION d1mach
      REAL big, biggest, ch(*), chtau(0:*), cmu(*),
     $     delm0, fbeam, gl(0:mxcmu,*), q0, q2, 
     $     taucpr(0:*), xba(*), xb0(mxcmu,*),
     $     xb1(mxcmu,*), zj(*), ylmc(0:mxcmu,*), ylm0(0:*)
      REAL deltat, pi, q0a, 
     $     q2a, sum
*
      pi = 2. * ASIN(1.0)
      biggest  = r1mach(2)
      big   = sqrt(biggest) / 1.E+10
*
*     Calculate x-sub-zero in STWJ(6d) 
*
      iqsav = 1
      DO 100 IQ = 1, NSTR
         SUM = 0.0
         DO 110 K = MAZ, NSTR-1
            SUM = SUM + GL(K,LC)*YLMC(K,IQ)*YLM0(K)
 110     CONTINUE
         ZJ(IQ) = (2.-DELM0)*FBEAM*SUM / (4.0*PI)
         IF ( zj(iq) .GT. zj(iqsav) ) iqsav = iq
 100  CONTINUE
*     
      q0a = EXP( -chtau(lc-1) )
      q2a = EXP( -chtau(lc) )
*     
*     Calculate alfa coefficient for one fixed angle
*     
      q0 = q0a*zj(iqsav)
      q2 = q2a*zj(iqsav)
      deltat = taucpr(lc) - taucpr(lc-1)
      xba(lc) = 1./ch(lc)
      IF( ABS(xba(lc)).GT.big .AND. taucpr(lc).GT. 1.0 )
     $     xba(lc) = 0.0
      IF( ABS(xba(lc)*taucpr(lc)) .GT. LOG(big))
     $     xba(lc) = 0.0
*
*     Dither alfa if it is close to one of the quadrature angles
*
      IF (  ABS(xba(lc)) .gt. 0.00001 ) THEN
         DO IQ = 1, NSTR/2
            IF ( ABS(( ABS(xba(LC)) - 1./CMU(IQ))/xba(lc) ) .LT. 0
     $           .05 )
     $           xba(LC) = xba(LC) * 1.001
         ENDDO
      ENDIF
      DO iq = 1, nstr
         q0 = q0a * zj(iq)
         q2 = q2a * zj(iq)
*     
*     x-sub-zero and x-sub-one in Eqs. KS(48-49)
*     
         xb1(iq,lc) = (1./deltat)*(q2*EXP(xba(lc)*taucpr(lc))
     $        -q0*EXP(xba(lc)*taucpr(lc-1)))
         xb0(iq,lc) = q0 * EXP(xba(lc)*taucpr(lc-1))-
     $        xb1(iq,lc)*taucpr(lc-1)
      ENDDO
*     
      RETURN
      END
*
      SUBROUTINE setpar( acmu, beta, lc, mxcmu, nstr, parmu, radius,
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
      REAL acmu(*), beta(0:*), parmu(mxcmu, *), q0, q2, 
     $     radius, taucpr(0:*),xba(*), xb0(mxcmu,*), xb1(mxcmu,*), 
     $     zd(0:*), zj(*), zj0(*), zj2(*)
      REAL betafact0, betafact1, deltat, fact
*
      fact  = 1.E+05         ! Convert from km to cm
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
         IF (beta(lc-1) .EQ. 0.0 ) THEN
            zj0(iq) = 0.0 
         ELSE
            zj0(iq)  =  -((1.0-acmu(iq)*acmu(iq)) / betafact0)
     $           * parmu( iq, lc )
         ENDIF
         zj2(iq)  = -((1.0-acmu(iq)*acmu(iq)) / betafact1)
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
      xba(lc) = 0.0
      DO iq = 1, nstr
         q0 = zj0(iq)
         q2 = zj2(iq)
*     
*     X-sub-zero and X-sub-one in Eqs. 48-49, Kylling and Stamnes
*     (1992) for alfa=0, i.e. linear approximation.
*     
         xb1(IQ,LC) = 0.0
         IF ( deltat .GT. 1.0E-09 )  xb1(iq,lc) =
     $        (1./deltat)*(q2-q0)
         xb0(iq,lc) = q0 - xb1(iq,lc)*taucpr(lc-1)
      ENDDO
*
      RETURN
*
      END
      SUBROUTINE  sFLUXES( CH, CMU, CWT, FBEAM, GC, il, KK, LAYRU, LL, 
     $                    LYRCUT, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,
     $                    PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                    XR0, XR1, XRA, ZPLK0, ZPLK1, ZPLKA, ZBEAMA, 
     $                    ZBEAM0, ZBEAM1, DFDT, FLUP, FLDN, FLDIR, 
     $                    RFLDIR, RFLDN, UAVG, U0C, MAXULV,
     $                    uavgdn, uavgso, uavgup)
*
*     Modified from original disort fluxes routine to include il in the
*     parameter list and to use the exponential-linear approximation
*     for the beam source.
*
C
C       CALCULATES THE RADIATIVE FLUXES, MEAN INTENSITY, AND FLUX
C       DERIVATIVE WITH RESPECT TO OPTICAL DEPTH FROM THE M=0 INTENSITY
C       COMPONENTS (THE AZIMUTHALLY-AVERAGED INTENSITY)
C
C    I N P U T     V A R I A B L E S:
C
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
C                   BY SOLVING SCALED VERSION OF EQ. SC(5);
C                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
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
C       FNET     :  NET FLUX (TOTAL-DOWN - DIFFUSE-UP)
C       FACT     :  EXP( - UTAUPR / UMU0 )
C       PLSORC   :  PLANCK SOURCE FUNCTION (THERMAL)
C       ZINT     :  INTENSITY OF m = 0 CASE, IN EQ. SC(1)
C+---------------------------------------------------------------------+
C
      LOGICAL LYRCUT, PRNT(*)
      REAL    DFDT(*), FLUP(*), FLDIR(*), FLDN(*),
     $     RFLDIR(*), RFLDN(* ),
     $     U0C( MXCMU,MXULV ), UAVG(*), uavgdn(*), uavgso(*),
     $     uavgup(*), umu0
      INTEGER LAYRU(*)
      REAL    CH(*), CMU(*), CWT(*),
     $     FBEAM, GC( MXCMU,MXCMU,* ), KK( MXCMU,* ),
     $     LL( MXCMU,* ), pi, SSALB(*), TAUCPR( 0:* ),
     $     UTAU(*), UTAUPR(*), XR0(*), XR1(*), XRA(*), 
     $     ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), ZPLKA( * ) ,
     $     ZBEAMA(*), ZBEAM0(MXCMU,*), ZBEAM1(MXCMU,*) 
      REAL ang1, ang2, dirint, fact, fdntot,
     $     fnet, plsorc, zint
C
C
      IF ( PRNT(2) )  WRITE( *,1010 )
C                                          ** ZERO DISORT OUTPUT ARRAYS
      CALL  ZEROIT( U0C, MXULV*MXCMU )
      CALL  ZEROIT( RFLDIR, MAXULV )
      CALL  ZEROIT( FLDIR,  MXULV )
      CALL  ZEROIT( RFLDN,  MAXULV )
      CALL  ZEROIT( FLDN,   MXULV )
      CALL  ZEROIT( FLUP,   MAXULV )
      CALL  ZEROIT( UAVG,   MAXULV )
      CALL  ZEROIT( UAVGdn, MAXULV )
      CALL  ZEROIT( UAVGso, MAXULV )
      CALL  ZEROIT( UAVGup, MAXULV )
      CALL  ZEROIT( DFDT,   MAXULV )
C                                        ** LOOP OVER USER LEVELS
      DO 100  LU = 1, NTAU
C
         LYU = LAYRU(LU)
C
         IF ( LYRCUT .AND. LYU.GT.NCUT ) THEN
C                                                ** NO RADIATION REACHES
C                                                ** THIS LEVEL
            FDNTOT = 0.0
            FNET   = 0.0
            PLSORC = 0.0
            GO TO 90
         END IF
*
         IF ( FBEAM.GT.0.0 )  THEN
            FACT  = EXP( - UTAUPR(LU) / CH(LYU) )
            DIRINT = FBEAM * FACT
            FLDIR(  LU ) = ABS(UMU0) * ( FBEAM * FACT )
            RFLDIR( LU ) = ABS(UMU0) * FBEAM * EXP( - UTAU( LU ) /
     $           CH(LYU) )
         ELSE
            DIRINT = 0.0
            FLDIR(  LU ) = 0.0
            RFLDIR( LU ) = 0.0
         END IF
C
         DO 20  IQ = 1, NN
C
            ZINT = 0.0
            DO 10  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
10          CONTINUE
            DO 11  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXP( - KK(JQ,LYU) *
     $              (UTAUPR(LU) - TAUCPR(LYU-1)) )
11          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.0 .OR. il.GT. 0 )  U0C( IQ,LU ) = ZINT + 
     $                           EXP(-ZBEAMA(LYU)*UTAUPR(LU)) *
     $                           ( ZBEAM0(IQ,LYU) +
     $                             ZBEAM1(IQ,LYU)*UTAUPR(LU) )
            IF ( il .EQ. 0 ) U0C( IQ,LU ) = U0C( IQ,LU ) +
     $           EXP(-ZPLKA(LYU)*UTAUPR(LU))*
     $           (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
            UAVG(LU) = UAVG(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )
            UAVGdn(LU) = UAVGdn(LU) + CWT(NN+1-IQ) * U0C( IQ,LU )
            FLDN(LU) = FLDN(LU) + CWT(NN+1-IQ)*CMU(NN+1-IQ) * U0C(IQ,LU)
20       CONTINUE
C
         DO 40  IQ = NN+1, NSTR
C
            ZINT = 0.0
            DO 30  JQ = 1, NN
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXP( - KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)) )
30          CONTINUE
            DO 31  JQ = NN+1, NSTR
               ZINT = ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *
     $                EXP( - KK(JQ,LYU) *
     $              (UTAUPR(LU) - TAUCPR(LYU-1)) )
31          CONTINUE
C
            U0C( IQ,LU ) = ZINT
            IF ( FBEAM.GT.0.0 .OR. il.GT.0 )  U0C( IQ,LU ) = ZINT + 
     $                           EXP(-ZBEAMA(LYU)*UTAUPR(LU))*
     $                           ( ZBEAM0(IQ,LYU) +
     $                             ZBEAM1(IQ,LYU)*UTAUPR(LU) )
            IF ( il .EQ. 0 ) U0C( IQ,LU ) = U0C( IQ,LU ) + 
     $           EXP(-ZPLKA(LYU)*UTAUPR(LU))*
     $           (ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
            UAVG(LU) = UAVG(LU) + CWT(IQ-NN) * U0C( IQ,LU )
            UAVGup(LU) = UAVGup(LU) + CWT(IQ-NN) * U0C( IQ,LU )
            FLUP(LU) = FLUP(LU) + CWT(IQ-NN) * CMU(IQ-NN) * U0C( IQ,LU )
40       CONTINUE
C
         FLUP( LU )  = 2.0 * PI * FLUP( LU )
         FLDN( LU )  = 2.0 * PI * FLDN( LU )
         FDNTOT = FLDN( LU ) + FLDIR( LU )
         FNET   = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU ) = ( 2.0 * PI * UAVG(LU) + DIRINT ) / ( 4.*PI )
         UAVGdn( LU ) = ( 2.0 * PI * UAVGdn(LU) ) / ( 4.*PI )
         UAVGso( LU ) =  DIRINT  / ( 4.*PI )
         UAVGup( LU ) = ( 2.0 * PI * UAVGup(LU) ) / ( 4.*PI )
         IF ( il .EQ. 0 )  THEN
            PLSORC = EXP(-XRA(LYU)*UTAUPR(LU))*
     $           (XR0(LYU) + XR1(LYU)* UTAUPR(LU))
         ELSE
            plsorc = 0.0
         ENDIF
         DFDT( LU ) = ( 1.0-SSALB(LYU) ) * 4.*PI*
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
               ANG1 = 180./PI * ACOS( CMU(2*NN-IQ+1) )
               ANG2 = 180./PI * ACOS( CMU(IQ) )
               WRITE( *,1120 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                         ANG2, CMU(IQ),        U0C(IQ+NN,LU)
200      CONTINUE
      END IF
C
1010  FORMAT( //, 21X,
     $ '<----------------------- FLUXES ----------------------->', /,
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
      SUBROUTINE  sPRTINP( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     $                    WVNMHI, NTAU, UTAU, NSTR, numu, umu,
     $                    FBEAM, UMU0, spher, nil,
     $                    FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, DELTAM, 
     $                    NOPLNK, onlyfl, FLYR, LYRCUT, OPRIM, 
     $                    TAUC, TAUCPR, MAXCMU, PRTMOM )
C
C        PRINT VALUES OF INPUT VARIABLES
C
      LOGICAL  DELTAM, LYRCUT, NOPLNK, onlyfl, PRTMOM, spher
      REAL albedo, btemp, DTAUC(*), fbeam,
     $     fisot, FLYR(*), OPRIM(*), 
     $     PMOM( 0:MAXCMU,* ), SSALB(*), UTAU(*), TAUC( 0:* ),
     $     TAUCPR( 0:* ), temis, TEMPER( 0:* ), ttemp, 
     $     umu(*), umu0, wvnmhi, wvnmlo
      REAL yessct
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
      YESSCT = 0.0
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
      IF( PRTMOM .AND. YESSCT.GT.0.0 )  THEN
         WRITE( *, '(/,A)' )  ' LAYER   PHASE FUNCTION MOMENTS'
         DO 20 LC = 1, NLYR
            IF( SSALB(LC).GT.0.0 )
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
 1113 FORMAT( ' Uses plane-parallel geometry' )
 1120 FORMAT( ' USES DELTA-M METHOD' )
 1130 FORMAT( ' DOES NOT USE DELTA-M METHOD' )
 1150 FORMAT( ' CALCULATE FLUXES AND INTENSITIES' )
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
      SUBROUTINE  ssetdis( acmu, albedo, bdr, beta, bplank, btemp, ch,
     $     chtau, cmu, cwt, deltam, deltamu, dtauc, expbea,
     $     fbeam, flyr, gl, layru, layru_c, lyrcut, maxumu, maxcmu,
     $     ncut, 
     $     newgeo, nlyr, ntau, ntau_c, nn, nstr, noplnk, numu,
     $     onlyfl, oprim, pkag, pkagc, pmom, radius, spher, ssalb,
     $     tauc, taucpr, temis, temper, tplank, ttemp, utau,
     $     utau_c, utaupr, utaupr_c, umu0, umu, usrtau,
     $     usrang,  wvnmlo, 
     $     wvnmhi, zd, numu_c, umu_c, numu_u, umu_u)
*
*     Modified from original disort setdis routine to do various setup
*     operations required by the spherical version of disort.
*
C
C          PERFORM MISCELLANEOUS SETTING-UP OPERATIONS
C
C       ROUTINES CALLED:  ERRMSG, QGAUSN, ZEROIT
C
C       INPUT :  ALL ARE 'DISORT' INPUT VARIABLES (SEE DOC FILE)
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
C+---------------------------------------------------------------------+
C   LOCAL SYMBOLIC DIMENSIONS:
C
c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
C+---------------------------------------------------------------------+
      
      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
      INCLUDE 'DISORT.MXD'
      PARAMETER ( MI = MXCMU/2, MI9M2 = 9*MI-2,
     $     NNLYRI = MXCMU*MXCLY, MXIL = 20 )
*


      LOGICAL  DELTAM, LYRCUT, newgeo, NOPLNK, onlyfl, SPHER, usrang,
     $     USRTAU
      INTEGER  layru(*),  layru_c(*), ntau_c
      REAL acmu(*), albedo, bdr( mi,0:* ), beta(0:*),
     $     bplank, btemp, ch(*), chtau(0:*), cmu(*), cwt(*), dtauc(*),
     $     expbea(0:*), fbeam, flyr(*), gl(0:mxcmu,*),
     $     oprim(*), pkag(0:*), pkagc(*), pmom(0:maxcmu,*), radius,
     $     ssalb(*), tauc(0:*), taucpr(0:*), taup, tempc, 
     $     temis, temper(0:*), tplank, ttemp, umu(*), umu0,
     $     utau(*), utau_c(*), utaupr(*),
     $     utaupr_c(*), wvnmhi, wvnmlo, zd(0:*)
*
      INTEGER numu_c, numu_u
      REAL umu_c(*), umu_u(*)
      REAL abstau, deltamu, f, pi, zenang
      REAL plkavg
      INTEGER nfac1(mxcly), nfac2(mxcly), nfac3(mxcly)
      REAL dtaucpr(mxcly),
     $     fac1(mxcly,2*mxcly), fac2(mxcly,2*mxcly),
     $     fac3(mxcly,2*mxcly), z_lay
      REAL chapman
*
      SAVE nfac1, nfac2, nfac3, fac1, fac2, fac3
*     
      DATA  ABSCUT / 400. /
C
C
      PI = 2. * ASIN( 1.0 )
C
      CALL  ZEROIT( BDR, MI*(MI+1) )
C
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
      CALL  ZEROIT( TAUCPR(0), MXCLY+1 )
      CALL  ZEROIT( EXPBEA(0), MXCLY+1 )
      CALL  ZEROIT( FLYR, MXCLY )
      CALL  ZEROIT( GL, (MXCMU+1)*MXCLY )
      CALL  ZEROIT( OPRIM, MXCLY )
      ABSTAU = 0.0
      chtau(0) = 0.0
*
      DO  60  LC = 1, NLYR
         PMOM(0,LC) = 1.0
         IF ( ABSTAU.LT.ABSCUT )  NCUT = LC
         ABSTAU = ABSTAU + ( 1. - SSALB(LC) ) * DTAUC(LC)
C     
         IF ( .NOT.DELTAM )  THEN
            OPRIM(LC) = SSALB(LC)
            TAUCPR(LC) = TAUC(LC)
            DO 40  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) * PMOM(K,LC)
 40         CONTINUE
            dtaucpr(lc)= dtauc(lc)
            F = 0.0
         ELSE
C                                    ** DO DELTA-M TRANSFORMATION
            F = PMOM( NSTR,LC )
            beta(lc)  = beta(lc) * ( 1. - f * ssalb(lc) )
            OPRIM(LC) = SSALB(LC) * ( 1. - F ) /
     $           ( 1. - F * SSALB(LC) )
            TAUCPR(LC) = TAUCPR(LC-1) +
     $           ( 1. - F*SSALB(LC) ) * DTAUC(LC)
            DO 50  K = 0, NSTR-1
               GL(K,LC) = (2*K+1) * OPRIM(LC) *
     $              (PMOM(K,LC)-F) / (1.-F)
 50         CONTINUE
            dtaucpr(lc)= taucpr(lc) - taucpr(lc-1)
         ENDIF
         FLYR(LC) = F
 60   CONTINUE
*
      IF ( fbeam .GT. 0.0 ) THEN
*
         zenang      = ACOS(umu0) * 180. / pi
         IF ( spher .and. newgeo ) THEN
            z_lay = 0.0
            CALL geofac( fac1, nfac1, z_lay, nlyr, zd, zenang, radius )
            z_lay = 0.5
            CALL geofac( fac2, nfac2, z_lay, nlyr, zd, zenang, radius )
            z_lay = 1.0
            CALL geofac( fac3, nfac3, z_lay, nlyr, zd, zenang, radius )
         ENDIF
         IF ( spher ) THEN
            expbea( 0 ) = 1.0
            IF( umu0 .LT. 0.0 ) THEN 
               expbea(0) =
     $              EXP(-chapman( 1, nfac1, fac1, dtaucpr ))
cakym               chtau(0) =  chapman( 1, nfac3, fac3, dtaucpr )
cakym               expbea(0) = EXP( -chtau(0) )
cakyp
*               write(11,*) chapman( 1, nfac1, fac1, dtaucpr ),
*     $              chtau(0), expbea(0)
cakyp
            ENDIF
*
            DO lc = 1, ncut
               taup       = tauc(lc-1) + dtauc(lc)/2.0
               chtau(lc)  = chapman( lc, nfac1, fac1, dtaucpr )
               ch(lc)     = taup/chapman( lc, nfac2, fac2, dtauc )
cakyp
*               write(*,*)lc,taup,tauc(lc-1),tauc(lc), chtau(lc-1),
*     $              chtau(lc),
*     $              chapman( lc, nfac2, fac2, dtauc ), ch(lc)
cakyp
               expbea(lc) = EXP(-chapman( lc, nfac1, fac1, dtaucpr ) )
cakyp
*               write(*,'(i2, 5(2x,F10.5))') lc, taup, chtau(lc), ch(lc), 
*     $              chapman( lc, nfac2, fac2, dtauc ), expbea(lc)
cakyp
            ENDDO
         ELSE IF ( .NOT. spher ) THEN
            DO lc = 1, ncut            
               ch(lc) = umu0
               chtau(lc)= taucpr(lc) / umu0 
               expbea(lc) = EXP( - taucpr(lc) / umu0 )
            ENDDO
         ENDIF
      ENDIF
*
C                      ** IF NO THERMAL EMISSION, CUT OFF MEDIUM BELOW
C                      ** ABSORPTION OPTICAL DEPTH = ABSCUT ( NOTE THAT
C                      ** DELTA-M TRANSFORMATION LEAVES ABSORPTION
C                      ** OPTICAL DEPTH INVARIANT ).  NOT WORTH THE
C                      ** TROUBLE FOR ONE-LAYER PROBLEMS, THOUGH.
      LYRCUT = .FALSE.
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
C
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
      CALL  QGAUSN( NN, CMU, CWT )
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
      deltamu  =  0.001
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
         deltamu  =  0.001
         ju       =  1
*
         DO iu = 1, numu_u, 2   ! Next all derivative angles
            umu_u(iu)   = umu(ju) - deltamu
            umu_u(iu+1) = umu(ju) + deltamu
            IF ( umu_u(iu) .LT. -1.0 ) THEN
               umu_u(iu) = -1.0
               umu_u(iu+1) = -1.0 + 2*deltamu
            ENDIF
            IF ( umu_u(iu+1) .GT. 1.0 ) THEN
               umu_u(iu+1) = 1.0
               umu_u(iu) = 1.0 - 2*deltamu
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
         BPLANK = 0.0
         TPLANK = 0.0
         CALL  ZEROIT( PKAG, MXCLY+1 )
         CALL  ZEROIT( PKAGC, MXCLY )
      ELSE
         TPLANK = TEMIS * PLKAVG( WVNMLO, WVNMHI, TTEMP )
         BPLANK =         PLKAVG( WVNMLO, WVNMHI, BTEMP )
         DO 180  LEV = 0, NLYR
            PKAG( LEV ) = PLKAVG( WVNMLO, WVNMHI, TEMPER(LEV) )
180      CONTINUE
         DO 190 LC = 1, NLYR
            TEMPC = 0.5 * ( TEMPER(LC-1) + TEMPER(LC) )
            PKAGC( LC ) = PLKAVG( WVNMLO, WVNMHI, TEMPC )
190      CONTINUE
      END IF
C
      RETURN
      END
*
      SUBROUTINE  sSETMTX( BDR, CBAND, CMU, CWT, DELM0, GC, KK, 
     $                    LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,
     $                    NN, NSTR, TAUCPR, WK )
C
C        CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS
C        OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-
C        OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE
C        SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES
C
C     ROUTINES CALLED:  ZEROIT
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
      LOGICAL LYRCUT
      REAL    BDR( MI,0:* ), CBAND( MI9M2,NNLYRI ),
     $     CMU(*), CWT(*), delm0, GC( MXCMU,MXCMU,* ),
     $     KK( MXCMU,* ), TAUCPR( 0:* ), WK(*)
      REAL expa, sum
C
C
      CALL  ZEROIT( CBAND, MI9M2*NNLYRI )
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
            WK(IQ) = EXP( KK(IQ,LC) * (TAUCPR(LC) - TAUCPR(LC-1)) )
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
         EXPA = EXP( KK(IQ,1) * TAUCPR(1) )
         IROW = NSHIFT - JCOL + NN
         DO 35  JQ = NN, 1, -1
            CBAND(IROW,JCOL+1) = GC(JQ,IQ,1) * EXPA
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
               SUM = 0.0
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
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) * EXPA
            ELSE
               SUM = 0.0
               DO 75  K = 1, NN
                  SUM = SUM + CWT(K) * CMU(K) * BDR(JQ-NN,K)
     $                        * GC(NN+1-K,IQ,NCUT)
75             CONTINUE
               CBAND(IROW,NNCOL) = ( GC(JQ,IQ,NCUT)
     $                               - (1.+DELM0) * SUM ) * EXPA
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
      SUBROUTINE  sSOLVE0( ALBEDO, B, BPLANK, CBAND, CMU, CWT, 
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
C     ROUTINES CALLED:  SGBCO, SGBSL, ZEROIT
C
C     I N P U T      V A R I A B L E S:
C
C       BPLANK   :  BOTTOM BOUNDARY THERMAL EMISSION
C       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
C                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
C                   BY LINPACK SOLUTION ROUTINES
C       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
C       EXPBEA   :  TRANSMISSION OF INCIDENT BEAM, EXP(-TAUCPR/UMU0)
*       il       :  iteration index
C       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ      :  ORDER OF AZIMUTHAL COMPONENT
C       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
C       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
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
C                   SC(12), CONSTANTS OF INTEGRATION WITHOUT
C                   EXPONENTIAL TERM
C      LL        :  PERMANENT STORAGE FOR -B-, BUT RE-ORDERED
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
      LOGICAL  LYRCUT, noplnk
      INTEGER  IPVT(*)
      REAL albedo, B(*), bplank,
     $     CBAND( MI9M2,NNLYRI ),
     $     CMU(*), CWT(*), EXPBEA(0:*), fbeam, fisot, LL( MXCMU,* ),
     $     pi, TAUCPR( 0:* ), tplank, umu0, Z(*), ZPLK0( MXCMU,* ),
     $     ZPLK1( MXCMU,* ), ZPLKA( * ), ZBEAM0( MXCMU,* ),
     $     ZBEAM1( MXCMU,* ), ZBEAMA( * ) 
      REAL rcond, sum
C
C
      CALL  ZEROIT( B, NNLYRI )
      IF ( pi .eq. 0.0 ) pi = 2.0 * ASIN(1.0)
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
               B(NCOL-NN+IQ) = - EXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
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
     $              +  EXP(-ZBEAMA(LC+1)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC+1)+ZBEAM1(IQ,LC+1)*TAUCPR(LC))
     $              -  EXP(-ZBEAMA(LC)*TAUCPR(LC))*
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
     $              -  EXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $              (ZPLK0(IQ+NN,NCUT)+ZPLK1(IQ+NN,NCUT)*TAUCPR(NCUT))
 60         CONTINUE
C
         ELSE
C
            DO 80 IQ = 1, NN
C
               SUM = 0.
               DO 70 JQ = 1, NN
                  SUM = SUM + CWT(JQ) * CMU(JQ) * ALBEDO
     $                 *(    EXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
     $                 (ZPLK0(NN+1-JQ,NCUT)+ZPLK1(NN+1-JQ,NCUT)
     $                 *TAUCPR(NCUT)))
 70            CONTINUE
               B(NCOL-NN+IQ) = 2.*SUM + (1.-ALBEDO) * BPLANK
     $              - EXP(-ZPLKA(NCUT)*TAUCPR(NCUT))*
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
     $              +  EXP(-ZPLKA(LC+1)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC+1)+ZPLK1(IQ,LC+1)*TAUCPR(LC))
     $              -  EXP(-ZPLKA(LC)*TAUCPR(LC))*
     $              (ZPLK0(IQ,LC)+ZPLK1(IQ,LC)*TAUCPR(LC))
 90         CONTINUE
 100     CONTINUE
C     
      ELSE IF ( FBEAM .GT. 0.0 .OR. il.GT.0 ) THEN
C
         DO 150 IQ = 1, NN
C                                      ** TOP BOUNDARY
C
            B(IQ) = - ZBEAM0(NN+1-IQ,1)
 150     CONTINUE
C
         IF ( LYRCUT ) THEN
C
            DO 160 IQ = 1, NN
C                                      ** BOTTOM BOUNDARY
               B(NCOL-NN+IQ) = - EXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $              (ZBEAM0(IQ+NN,NCUT)+ZBEAM1(IQ+NN,NCUT)*TAUCPR(NCUT)
     $              )
 160        CONTINUE
C
         ELSE
C
            DO 180 IQ = 1, NN
C
               SUM = 0.
               DO 170 JQ = 1, NN
                  SUM = SUM + CWT(JQ) * CMU(JQ) * ALBEDO
     $                 * ( 
     $                 EXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $                 (ZBEAM0(NN+1-JQ,NCUT)+ZBEAM1(NN+1-JQ,NCUT)
     $                 *TAUCPR(NCUT))
     $                 )
 170           CONTINUE
               B(NCOL-NN+IQ) = 2.*SUM + ( ALBEDO * UMU0*FBEAM/PI
     $              ) * EXPBEA(NCUT)
     $              -  EXP(-ZBEAMA(NCUT)*TAUCPR(NCUT))*
     $              (ZBEAM0(IQ+NN,NCUT)+ZBEAM1(IQ+NN,NCUT)*TAUCPR(NCUT)
     $              )
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
     $              +  EXP(-ZBEAMA(LC+1)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC+1)+ZBEAM1(IQ,LC+1)*TAUCPR(LC))
     $              -  EXP(-ZBEAMA(LC)*TAUCPR(LC))*
     $              (ZBEAM0(IQ,LC)+ZBEAM1(IQ,LC)*TAUCPR(LC))
 190        CONTINUE
 200     CONTINUE
C     
      END IF
C
C                     ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
C                     ** OF BAND MATRIX -CBAND- AND TEST IF IT IS NEARLY
C                     ** SINGULAR (NOTE: -CBAND- IS DESTROYED)
C                     ** (-CBAND- IS IN LINPACK PACKED FORMAT)
      RCOND = 0.0
      NCD = 3*NN - 1
      CALL  SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )
      IF ( 1.0+RCOND .EQ. 1.0 )  CALL  ERRMSG
     $   ( 'SOLVE0--SGBCO SAYS MATRIX NEAR SINGULAR',.FALSE.)
C
C                   ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -CBAND-
C                   ** AND R.H. SIDE(S) -B- AFTER -CBAND- HAS BEEN L-U
C                   ** DECOMPOSED.  SOLUTION IS RETURNED IN -B-.
C
      CALL  SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )
C
C                   ** ZERO -CBAND- (IT MAY CONTAIN 'FOREIGN'
C                   ** ELEMENTS UPON RETURNING FROM LINPACK);
C                   ** NECESSARY TO PREVENT ERRORS
C
      CALL  ZEROIT( CBAND, MI9M2*NNLYRI )
C
      DO 220  LC = 1, NCUT
         IPNT = LC*NSTR - NN
         DO 220  IQ = 1, NN
            LL(NN+1-IQ,LC) = B(IPNT+1-IQ)
            LL(IQ+NN,  LC) = B(IQ+IPNT)
220   CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  sterpso( cwt, fbeam, gl, maz, mxcmu, noplnk,  
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
*        zbau        exponential-linear-in-optical-depth approximation.
*       zp0u, zp1u:  Components of a exponential-linear-in-optical-
*        zpau        depth-dependent source (approximating the 
*                    Planck emission source), at users angles.
C
C   I N T E R N A L       V A R I A B L E S:
C
C       PSI0,  :  SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)
C        PSI1     WITH Z0 AND Z1 RESPECTIVELY
C+---------------------------------------------------------------------+
      LOGICAL  NOPLNK
      REAL     CWT(*), GL(0:*), fbeam, oprim,
     $     PSI0(*), PSI1(*), xr0, xr1, xra,  
     $     YLMC( 0:MXCMU,* ), YLMU( 0:MXCMU,*), Z0(*), Z1(*), 
     $     zbs0(*), zbs1(*), zbsa,
     $     zb0u(*), zb1u(*), zbau(*), zp0u(*), zp1u(*), zpau(*)
      REAL psum0, psum1, sum0, sum1
C     
C
      IF ( FBEAM.GT.0.0 )  THEN
C                                  ** BEAM SOURCE TERMS; EQ. SD(9)
         DO 20  IQ = MAZ, NSTR-1
            psum0 = 0.
            psum1 = 0.
            DO 10  JQ = 1, NSTR
               psum0 = psum0 + cwt(jq) * ylmc(iq,jq) * zbs0(jq)
               psum1 = psum1 + cwt(jq) * ylmc(iq,jq) * zbs1(jq)
10          CONTINUE
            psi0(iq+1) = 0.5 * gl(iq) * psum0
            psi1(iq+1) = 0.5 * gl(iq) * psum1
20       CONTINUE
C
         DO 40  IU = 1, numu
            sum0 = 0.
            sum1 = 0.
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
            PSUM0 = 0.0
            PSUM1 = 0.0
            DO 70  JQ = 1, NSTR
               PSUM0 = PSUM0 + CWT(JQ) * YLMC(IQ,JQ) * Z0(JQ)
               PSUM1 = PSUM1 + CWT(JQ) * YLMC(IQ,JQ) * Z1(JQ)
 70         CONTINUE
            PSI0(IQ+1) = 0.5 * GL(IQ) * PSUM0
            PSI1(IQ+1) = 0.5 * GL(IQ) * PSUM1
 80       CONTINUE
C
          DO 100  IU = 1, numu
            SUM0 = 0.0
            SUM1 = 0.0
            DO 90   IQ = MAZ, NSTR-1
               SUM0 = SUM0 + YLMU(IQ,IU) * PSI0(IQ+1)
               SUM1 = SUM1 + YLMU(IQ,IU) * PSI1(IQ+1)
90          CONTINUE
            zp0u(iu) = sum0 + (1.-oprim) * xr0
            zp1u(iu) = sum1 + (1.-oprim) * xr1
            zpau(iu) = xra
100      CONTINUE
C
      END IF
C
      RETURN
      END
      SUBROUTINE  supbeam( array, cc, cmu, ipvt, mxcmu, nn, nstr, 
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
*       xb0    :  expansion of beam source function Eq. KS(7)
*       xb1    :  expansion of beam source function Eq. KS(7)
*       xba    :  expansion of beam source function Eq. KS(7)
*       (remainder are 'disort' input variables)
*
*    O U T P U T    V A R I A B L E S:
*
*       zbs0     :  solution vectors z-sub-zero of Eq. KS(10-11)
*       zbs1     :  solution vectors z-sub-one  of Eq. KS(10-11)
*       zbsa     :  alfa coefficient in Eq. KS(7)
*       zbeam0, :  permanent storage for -zbs0,zbs1,zbsa-, but re-ordered
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
      INTEGER ipvt(*)
      REAL    ARRAY( MXCMU,* ), CC( MXCMU,* ),
     $     CMU(*), WK(*), XB0(*), XB1(*), XBA, ZBS0(*), ZBS1(*), ZBSA, 
     $     ZBEAM0(*), ZBEAM1(*), ZBEAMA
      REAL RCOND, RMIN
*
*
      CALL ZEROIT( ARRAY, MXCMU*MXCMU )
*
      DO IQ = 1, NSTR
         DO JQ = 1, NSTR
            ARRAY(IQ,JQ) = - CC(IQ,JQ)
         ENDDO
         ARRAY(IQ,IQ) = 1.0 + XBA*CMU(IQ) + ARRAY(IQ,IQ)
         ZBSA     = XBA
         ZBS1(IQ) = XB1(IQ)
      ENDDO
*
*     SOLVE LINEAR EQUATIONS: 
*
      RCOND = 0.0
      CALL  SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )
      RMIN = 1.0E-04
      IF ( RCOND .LT. RMIN ) THEN
*
*     Dither alfa if rcond to small
*
         xba = xba * (1.005) 
         DO iq = 1, nstr   
            DO jq = 1, nstr
               array(iq,jq) = - cc(iq,jq)
            ENDDO
            array(iq,iq) = 1.0 + xba*cmu(iq) + array(iq,iq)
            zbsa     = xba
            zbs1(iq) = xb1(iq)
         ENDDO
*
*     Solve linear equations KS(10-11)
*
         rcond = 0.0
         CALL  sgeco( array, mxcmu, nstr, ipvt, rcond, wk )
      ENDIF
*
      IF ( 1.0+rcond .EQ. 1.0 )  CALL  errmsg
     $   ( ' Upbeam--sgeco says matrix near singular',.FALSE.)
*
      CALL  sgesl( array, mxcmu, nstr, ipvt, zbs1, 0 )
*
      DO iq = 1, nstr
         zbs0(iq) = xb0(iq) + cmu(iq) * zbs1(iq)
      ENDDO
      CALL  sgesl( array, mxcmu, nstr, ipvt, zbs0, 0 )
*
*     ... and now some index gymnastic for the inventive ones...
*
      DO iq = 1, nn
         zbeam0( iq+nn )   = zbs0( iq )
         zbeam1( iq+nn )   = zbs1( iq )
         zbeama            = zbsa
         zbeam0( nn+1-iq ) = zbs0( iq+nn )
         zbeam1( nn+1-iq ) = zbs1( iq+nn )
      ENDDO
*
      RETURN
      END
      SUBROUTINE sUSRINT( albedo, BPLANK, CMU, CWT, DELM0, EXPBEA,
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
C       EXPBEA :  TRANSMISSION OF INCIDENT BEAM, EXP(-TAUCPR/UMU0)
C       GC     :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
C       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
C                 (i.e., g IN EQ. SC(1) )
C       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
C       LAYRU  :  LAYER NUMBER OF USER LEVEL -UTAU-
C       LL     :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
C                 BY SOLVING SCALED VERSION OF EQ. SC(5);
C                 EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
C       LYRCUT :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
C       MAZ    :  ORDER OF AZIMUTHAL COMPONENT
C       NCUT   :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
C       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
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
C       ALL THE EXPONENTIAL FACTORS ( EXP1, EXPN,... etc.)
C       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGRATION IN
C       EQ. SC(12) INTO EQS. S1(8-9).  THEY ALL HAVE NEGATIVE
C       ARGUMENTS SO THERE SHOULD NEVER BE OVERFLOW PROBLEMS.
C+---------------------------------------------------------------------+
C
      LOGICAL  LAMBER, LYRCUT, NOPLNK, NEGUMU
      INTEGER  LAYRU(*)
      REAL     albedo, bplank, CMU(*), CWT(*), delm0, denom, 
     $     dtau, dtau1, dtau2, EXPBEA(0:*), fbeam, fisot,
     $     GC( MXCMU,MXCMU,* ),
     $     GU( MXUMU,MXCMU,* ), KK( MXCMU,* ), LL( MXCMU,* ),
     $     pi, TAUCPR( 0:* ), tplank, uum( mxumu,mxulv, 0:* ),
     $     UMU(*), umu0, UTAUPR(*), WK(*),
     $     zb0u( mxumu,* ), zb1u( mxumu,* ), zbau( mxumu,* ),
     $     zp0u( mxumu,* ), zp1u( mxumu,* ), zpau( mxumu,* ),
     $     zbeam0( mxcmu,* ), zbeam1( mxcmu,* ), zbeama( * ),
     $     ZPLK0( MXCMU,* ), ZPLK1( MXCMU,* ), zplka( * ) 
      REAL alfa, bnddfu, bnddir, bndint, dfuint,
     $     emu, exp1, exp2, expn, palint, plkint, rmu, sgn
C
C
      CALL  ZEROIT( UUM, MXUMU*MXULV*(mxcmu+1) )
C
      emu = 1. - albedo
      rmu = albedo 
C                          ** INCORPORATE CONSTANTS OF INTEGRATION INTO
C                          ** INTERPOLATED EIGENVECTORS
      DO 10  LC = 1, NCUT
         DO  10  IQ = 1, NSTR
            DO 10  IU = 1, NUMU
               GU(IU,IQ,LC) = GU(IU,IQ,LC) * LL(IQ,LC)
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
            NEGUMU = UMU(IU).LT.0.0
            IF( NEGUMU )  THEN
               LYRSTR = 1
               LYREND = LYU - 1
               SGN = - 1.0
            ELSE
               LYRSTR = LYU + 1
               LYREND = NCUT
               SGN = 1.0
            END IF
C                          ** FOR DOWNWARD INTENSITY, INTEGRATE FROM TOP
C                          ** TO 'LYU-1' IN EQ. S1(8); FOR UPWARD,
C                          ** INTEGRATE FROM BOTTOM TO 'LYU+1' IN S1(9)
            PALINT = 0.0
            PLKINT = 0.0
            DO 30  LC = LYRSTR, LYREND
C
               DTAU = TAUCPR(LC) - TAUCPR(LC-1)
               EXP1 =  EXP( (UTAUPR(LU) - TAUCPR(LC-1)) / UMU(IU) )
               EXP2 =  EXP( (UTAUPR(LU) - TAUCPR( LC )) / UMU(IU) )
C
               IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 ) THEN
                 denom  =  SGN*1./(zpau(iu,lc)*UMU(IU)+1.)
                 PLKINT = PLKINT +
     $                    (zp0u(iu,lc)*denom*
     $                         (EXP(-zpau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                         -EXP(-zpau(iu,lc)*TAUCPR(LC)) *EXP2 )
     $                    +zp1u(iu,lc)*denom*
     $                        ( (TAUCPR(LC-1)+SGN*denom*UMU(IU) )
     $                          *EXP(-zpau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                         -(TAUCPR(LC)+SGN*denom*UMU(IU) )
     $                          *EXP(-zpau(iu,lc)*TAUCPR(LC))* EXP2) )
               ENDIF
C
               IF ( FBEAM.GT.0.0 .OR. il.GT.0 )  THEN
                 denom  =  SGN*1./(zbau(iu,lc)*UMU(IU)+1.)
                 palint = palint +
     $                    (zb0u(iu,lc)*denom*
     $                         (EXP(-zbau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                         -EXP(-zbau(iu,lc)*TAUCPR(LC)) *EXP2 )
     $                    +zb1u(iu,lc)*denom*
     $                        ( (TAUCPR(LC-1)+SGN*denom*UMU(IU) )
     $                          *EXP(-zbau(iu,lc)*TAUCPR(LC-1)) *EXP1
     $                         -(TAUCPR(LC)+SGN*denom*UMU(IU) )
     $                          *EXP(-zbau(iu,lc)*TAUCPR(LC))* EXP2) )
               ENDIF
C                                                   ** -KK- IS NEGATIVE
               DO 20  IQ = 1, NN
                  WK(IQ) = EXP( KK(IQ,LC) * DTAU )
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)
                  IF ( ABS(DENOM).LT.0.0001 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     EXPN = DTAU / UMU(IU) * EXP2
                  ELSE
                     EXPN = SGN * ( EXP1 * WK(IQ) - EXP2 ) / DENOM
                  END IF
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN
20             CONTINUE
C                                                   ** -KK- IS POSITIVE
               DO 21  IQ = NN+1, NSTR
                  DENOM = 1.0 + UMU(IU) * KK(IQ,LC)
                  IF ( ABS(DENOM).LT.0.0001 ) THEN
C                                                   ** L'HOSPITAL LIMIT
                     EXPN = - DTAU / UMU(IU) * EXP1
                  ELSE
                     EXPN = SGN *( EXP1 - EXP2 * WK(NSTR+1-IQ) ) / DENOM
                  END IF
                  PALINT = PALINT + GU(IU,IQ,LC) * EXPN
21             CONTINUE
C
30          CONTINUE
C                           ** CALCULATE CONTRIBUTION FROM USER
C                           ** OUTPUT LEVEL TO NEXT COMPUTATIONAL LEVEL
C
            DTAU1 = UTAUPR(LU) - TAUCPR(LYU-1)
            DTAU2 = UTAUPR(LU) - TAUCPR(LYU)
            IF( ABS(DTAU1).LT.1.E-6 .AND. NEGUMU )  GO TO 50
            IF( ABS(DTAU2).LT.1.E-6 .AND. (.NOT.NEGUMU) )  GO TO 50
            IF( NEGUMU ) EXP1 = EXP( DTAU1 / UMU(IU) )
            IF( .NOT.NEGUMU ) EXP2 = EXP( DTAU2 / UMU(IU) )
C
            IF ( FBEAM.GT.0.0 .OR. il.GT.0 )  THEN
              IF ( NEGUMU ) THEN
                 EXPN = EXP1
                 ALFA = zbau(iu,lyu)
                 denom = (-1./(ALFA*UMU(IU)+1.))
                 palint = palint + 
     $                    zb0u(iu,lyu)*denom*
     $                          (-EXP(-ALFA*UTAUPR(LU))
     $                          + EXPN*EXP(-alfa*taucpr(lyu-1)) )
     $                 +  zb1u(iu,lyu)*denom*
     $                        ( -(UTAUPR(LU)-UMU(IU)*denom)*
     $                            EXP(-ALFA*UTAUPR(LU))
     $                          +(TAUCPR(LYU-1)-UMU(IU)*denom) 
     $                           *EXPN*EXP(-alfa*taucpr(lyu-1)) )
              ELSE
                 EXPN = EXP2
                 ALFA = zbau(iu,lyu)
                 denom = (1./(ALFA*UMU(IU)+1.))
                 palint = palint + 
     $                    zb0u(iu,lyu)*denom*
     $                        ( EXP(-ALFA*UTAUPR(LU))
     $                           -EXP(-ALFA*TAUCPR(LYU))*EXPN)
     $                 +  zb1u(iu,lyu)*denom*
     $                        ( (UTAUPR(LU) +UMU(IU)*denom)
     $                             *EXP(-ALFA*UTAUPR(LU))
     $                         -(TAUCPR(LYU)+UMU(IU)*denom)
     $                             *EXP(-ALFA*TAUCPR(LYU)) *EXPN )
              END IF
            ENDIF
C                                                   ** -KK- IS NEGATIVE
            DTAU = TAUCPR(LYU) - TAUCPR(LYU-1)
            DO 40  IQ = 1, NN
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)
               IF ( ABS(DENOM).LT.0.0001 ) THEN
                  EXPN = - DTAU2 / UMU(IU) * EXP2
               ELSE IF ( NEGUMU ) THEN
                  EXPN = ( EXP( - KK(IQ,LYU) * DTAU2 ) -
     $                     EXP( KK(IQ,LYU) * DTAU ) * EXP1 ) / DENOM
               ELSE
                  EXPN = ( EXP( - KK(IQ,LYU) * DTAU2 ) - EXP2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN
40          CONTINUE
C                                                   ** -KK- IS POSITIVE
            DO 41  IQ = NN+1, NSTR
               DENOM = 1.0 + UMU(IU) * KK(IQ,LYU)
               IF ( ABS(DENOM).LT.0.0001 ) THEN
                  EXPN = - DTAU1 / UMU(IU) * EXP1
               ELSE IF ( NEGUMU ) THEN
                  EXPN = ( EXP(- KK(IQ,LYU) * DTAU1 ) - EXP1 ) / DENOM
               ELSE
                  EXPN = ( EXP( - KK(IQ,LYU) * DTAU1 ) -
     $                     EXP( - KK(IQ,LYU) * DTAU ) * EXP2 ) / DENOM
               END IF
               PALINT = PALINT + GU(IU,IQ,LYU) * EXPN
41          CONTINUE
C
            IF ( .NOT.NOPLNK .AND. MAZ.EQ.0 )  THEN
              IF ( NEGUMU ) THEN
                 EXPN = EXP1
                 ALFA = zpau(iu,lyu)
                 denom = (-1./(ALFA*UMU(IU)+1.))
                 PLKINT = PLKINT + 
     $                    zp0u(iu,lyu)*denom*
     $                        (-EXP(-ALFA*UTAUPR(LU))
     $                        + EXPN*EXP(-alfa*taucpr(lyu-1)) )
     $                 +  zp1u(iu,lyu)*denom*
     $                        ( -(UTAUPR(LU)-UMU(IU)*denom)*
     $                            EXP(-ALFA*UTAUPR(LU))
     $                          +(TAUCPR(LYU-1)-UMU(IU)*denom) 
     $                           *EXPN*EXP(-alfa*taucpr(lyu-1)) )
              ELSE
                 EXPN = EXP2
                 ALFA = zpau(iu,lyu)
                 denom = (1./(ALFA*UMU(IU)+1.))
                 PLKINT = PLKINT + 
     $                    zp0u(iu,lyu)*denom*
     $                        ( EXP(-ALFA*UTAUPR(LU))
     $                           -EXP(-ALFA*TAUCPR(LYU))*EXPN)
     $                 +  zp1u(iu,lyu)*denom*
     $                        ( (UTAUPR(LU) +UMU(IU)*denom)
     $                             *EXP(-ALFA*UTAUPR(LU))
     $                         -(TAUCPR(LYU)+UMU(IU)*denom)
     $                             *EXP(-ALFA*TAUCPR(LYU)) *EXPN )
              END IF
            END IF
C                            ** CALCULATE INTENSITY COMPONENTS
C                            ** ATTENUATED AT BOTH BOUNDARIES.
C                            ** NOTE:: NO AZIMUTHAL INTENSITY
C                            ** COMPONENT FOR ISOTROPIC SURFACE
50          BNDINT = 0.0
            IF ( il .EQ. 0 ) THEN
               IF ( NEGUMU .AND. MAZ.EQ.0 ) THEN
                  BNDINT = ( FISOT + TPLANK ) * EXP( UTAUPR(LU) /
     $                 UMU(IU) )
               ELSE IF ( .NOT.NEGUMU ) THEN
                  IF ( LYRCUT .OR. (LAMBER .AND. MAZ.GT.0) )  GO TO 90
                  DO 60  JQ = NN+1, NSTR
                     WK(JQ) = EXP(-KK(JQ,NLYR)*
     $                    (TAUCPR(NLYR)-TAUCPR(NLYR-1)))
 60               CONTINUE
                  BNDDFU = 0.0
                  DO 80  IQ = NN, 1, -1
                     DFUINT = 0.0
                     DO 70  JQ = 1, NN
                        DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
 70                  CONTINUE
                     DO 71  JQ = NN+1, NSTR
                        DFUINT = DFUINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
     $                       * WK(JQ)
 71                  CONTINUE
                     IF ( FBEAM.GT.0.0 )
     $                    DFUINT = DFUINT + 
     $                    EXP(-zbeama(nlyr)*taucpr(nlyr)) *
     $                    (zbeam0(iq,nlyr)+
     $                    zbeam1(iq,nlyr)*taucpr(nlyr))
                     DFUINT = DFUINT + DELM0 * ( 
     $                    EXP(-zplka(nlyr)*TAUCPR(NLYR))*
     $                    (ZPLK0(IQ,NLYR)+
     $                    ZPLK1(IQ,NLYR)*TAUCPR(NLYR)))
                     BNDDFU = BNDDFU + ( 1. + DELM0 ) * rmu
     $                    * CMU(NN+1-IQ) * CWT(NN+1-IQ) * DFUINT
 80               CONTINUE
C     
                  BNDDIR = 0.0
                  IF (FBEAM.GT.0.0  .OR. umu0.GT.0 ) BNDDIR = UMU0
     $                 * FBEAM / PI * rmu * EXPBEA(NLYR)
                  BNDINT = ( BNDDFU + BNDDIR + DELM0 * emu * BPLANK )
     $                 * EXP( (UTAUPR(LU)-TAUCPR(NLYR)) / UMU(IU) )
               END IF
            END IF
C
90          UUM( IU, LU, maz ) = PALINT + PLKINT + BNDINT
C
100      CONTINUE
200   CONTINUE
C
      RETURN
      END
*
      SUBROUTINE  SETCOE( KK, LC, NCUT, NN, Q0, Q1, Q2, TAUCPR,
     $     XR0, XR1, XRA ) 

C       Set coefficients for expansion of thermal source function, 
C       eq. KS(15) 

C       Routines called:  None

C       INPUT :  
C                Q0,Q1,Q2    Source at top, center and bottom of layer
C                TAUCPR      delta-M-scaled optical depth

C       OUTPUT:  
C                XR0         expansion of thermal source function
C                XR1         expansion of thermal source function
C                XRA         expansion of thermal source function eq. KS(15)

      INTEGER LC, NCUT
      LOGICAL FIRST
      SAVE big, biggest, first, small, smallest
      REAL big, biggest, small, smallest
      REAL   KK(*), TAUCPR( 0:* ), XR0, XR1, XRA 
C
      DATA FIRST /.TRUE./
C
C The calculation of the particular solutions require some care, small
C and big have been set so that no problems should occurr on 32-bits
C machine running single precision
C
      IF ( FIRST ) THEN
         SMALLEST = R1MACH(1)
         BIGGEST  = R1MACH(2)	 
	 BIG   = SQRT(BIGGEST)
	 SMALL = SMALLEST
	 FIRST = .FALSE.
      ENDIF
C
      XR1 = 0.0
      XRA = 0.0
C
C Set coefficients in KS(15) for thermal source
C
C Calculate alpha coefficient 
C
      DELTAT = TAUCPR(LC) - TAUCPR(LC-1)
C
C Case 1: source small at bottom layer
C
      IF ( (Q2.LT.(Q0*1.E-02) .OR. Q2.LE.SMALL )
     $     .AND. Q1.GT.SMALL .AND. Q0.GT.SMALL ) THEN
C
C alpha Eq. KS(50)
C
         XRA = (2./DELTAT) * LOG( Q0/Q1 )
	 IF ( XRA .GT. BIG )    XRA = BIG
	 IF ( XRA*TAUCPR(LC-1) .GE. ALOG(BIG) ) THEN
	    XR0 =  BIG
	 ELSE
	    XR0 = Q0
	 ENDIF
	 XR1 = 0.0
C
C Case 2: Source small at center and bottom of layer
C
      ELSE IF ( (Q2.LE.(Q1*1.E-02) .OR. Q2.LE.SMALL ) .AND.
     $     ((Q1.LE.(Q0*1.E-02)) .OR. (Q1.LE.SMALL))
     $     .AND. (Q0.GT.SMALL) ) THEN
C       
	 XRA  =   BIG / TAUCPR(NCUT)
	 XR0 = Q0
	 XR1 = 0.0
C
C Case 3:All sources zero
C
	ELSE IF ( Q2.LE.SMALL .AND. Q1.LE.SMALL
     $                           .AND. Q0.LE.SMALL) THEN
	   XRA = 0.0
	   XR0 = 0.0
	   XR1 = 0.0
C
C Case 4: Sources same at center, bottom and top of layer
C         or layer optically very thin
C
	ELSE IF ( (ABS((Q2-Q0)/Q2).LT.1.E-04) .AND.
     $          (ABS((Q2-Q1)/Q2).LT.1.E-04)
     $             .OR. DELTAT.LT. 1.E-04           ) THEN

	   XRA = 0.0
	   XR0 = Q0
	   XR1 = 0.0

C                         **  Case 5: Normal case
	ELSE

	 ARG = (Q1/Q2)**2. - Q0/Q2
	 IF ( ARG .LT. 0.0 ) ARG = 0.0

C alpha Eq. KS(44). For source that has its maximum the top of the
C layer, use negative solution

	 SGN = 1.0
	 IF ( Q0 .GT. Q2 ) SGN = -1.
	 FACT3 = LOG(Q1/Q2 + SGN*SQRT(ARG) )
	 IF ( ABS(FACT3) .LE. 0.005 ) THEN ! BE CAREFUL WITH LOG OF
	    Q1 = 0.99 * Q1	! NUMBERS CLOSE TO ONE
	    FACT3 = LOG(Q1/Q2 + SGN*SQRT(ARG) )
	 ENDIF 
	 XRA = (2./DELTAT) * FACT3
	 IF(ABS(XRA*TAUCPR(LC)) .GT. 
     $             (LOG(BIGGEST)-LOG(Q0*100.) ) )    XRA = 0.0 

C Dither alpha if it is close to an eigenvalue
	 
	 DO IQ = 1, NN
	    IF ( ABS((XRA-KK(IQ))/XRA) .LT. 1.E-4 ) THEN
	       XRA = 1.01*XRA
	    ENDIF
	 ENDDO

C
C Set constants in Eqs. KS(45-46)
C
	 XR1 = (1./DELTAT)*(Q2*EXP(XRA*TAUCPR(LC)) 
     $                                 -Q0*EXP(XRA*TAUCPR(LC-1)))
	 XR0 = Q0 * EXP(XRA*TAUCPR(LC-1)) -
     $                                   XR1*TAUCPR(LC-1)
      ENDIF

      RETURN
      END
      SUBROUTINE  LEPOLYs( NMU, M, MAXMU, TWONM1, MU, YLM )

C       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
C       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C       PLM = P-SUB-L-SUPER-M AS

C             YLM(MU) = SQRT( (L-M)!/(L+M)! ) * PLM(MU)

C       FOR FIXED ORDER -M- AND ALL DEGREES FROM L = M TO TWONM1.
C       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE
C       FROM A PRIOR CALL TO THE ROUTINE.

C       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
C                  High-Order Associated Legendre Polynomials,
C                  J. Quant. Spectrosc. Radiat. Transfer 10,
C                  557-562, 1970.  (hereafter D/A)

C       METHOD: Varying degree recurrence relationship.

C       NOTE 1: The D/A formulas are transformed by
C               setting  M = n-1; L = k-1.
C       NOTE 2: Assumes that routine is called first with  M = 0,
C               then with  M = 1, etc. up to  M = TWONM1.
C       NOTE 3: Loops are written in such a way as to vectorize.

C  I N P U T     V A R I A B L E S:

C       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
C       M      :  ORDER OF -YLM-
C       MAXMU  :  FIRST DIMENSION OF -YLM-
C       TWONM1 :  MAX DEGREE OF -YLM-
C       MU(I)  :  I = 1 TO NMU, ARGUMENTS OF -YLM-
C       IF M.GT.0, YLM(M-1,I) FOR I = 1 TO NMU IS REQUIRED

C  O U T P U T     V A R I A B L E:

C       YLM(L,I) :  L = M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
C                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
C+---------------------------------------------------------------------+
      REAL     MU(*), YLM( 0:MAXMU,* )
      INTEGER  M, NMU, TWONM1
      PARAMETER  ( MAXSQT = 1000 )
      REAL     SQT( MAXSQT )
      LOGICAL  PASS1
      SAVE  SQT, PASS1
      DATA  PASS1 / .TRUE. /

      IF ( PASS1 )  THEN
         PASS1 = .FALSE.
         DO 1  NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT(NS) )
    1    CONTINUE
      ENDIF

      IF ( 2*TWONM1 .GT. MAXSQT )
     $     CALL ERRMSG( 'LEPOLYs--NEED TO INCREASE PARAM MAXSQT',.TRUE.)

      IF ( M .EQ. 0 )  THEN
C                             ** UPWARD RECURRENCE FOR ORDINARY
C                             ** LEGENDRE POLYNOMIALS
         DO  10  I = 1, NMU
            YLM( 0,I ) = 1.
            YLM( 1,I ) = MU( I )
 10      CONTINUE

         DO  20  L = 2, TWONM1
            DO  20  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $              - ( L-1 ) * YLM( L-2,I ) ) / L
 20      CONTINUE

      ELSE

         DO  30  I = 1, NMU
C                               ** Y-SUB-M-SUPER-M; DERIVED FROM
C                               ** D/A EQS. (11,12)

            YLM( M,I) = - SQT( 2*M-1 ) / SQT( 2*M )
     $           * SQRT( 1. - MU(I)**2 ) * YLM( M-1,I )

C                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM
C                              ** D/A EQS. (13,14) USING EQS. (11,12)

            YLM( M+1,I ) = SQT( 2*M+1 ) * MU(I) * YLM( M,I )
 30      CONTINUE
C                                   ** UPWARD RECURRENCE; D/A EQ. (10)
         DO  40  L = M+2, TWONM1
            TMP1 = SQT( L-M ) * SQT( L+M )
            TMP2 = SQT( L-M-1 ) * SQT( L+M-1 )
            DO  40  I = 1, NMU
               YLM( L,I ) = ( ( 2*L-1 ) * MU(I) * YLM( L-1,I )
     $              - TMP2 * YLM( L-2,I ) ) / TMP1
40       CONTINUE

      END IF
      
      RETURN
      END
      SUBROUTINE  SOLEIGs( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM,
     $     MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,
     $     KK, GC )

C         SOLVES EIGENVALUE/VECTOR PROBLEM NECESSARY TO CONSTRUCT
C         HOMOGENEOUS PART OF DISCRETE ORDINATE SOLUTION; STWJ(8B)
C         ** NOTE ** EIGENVALUE PROBLEM IS DEGENERATE WHEN SINGLE
C                    SCATTERING ALBEDO = 1;  PRESENT WAY OF DOING IT
C                    SEEMS NUMERICALLY MORE STABLE THAN ALTERNATIVE
C                    METHODS THAT WE TRIED

C     ROUTINES CALLED:  ASYMTX

C   I N P U T     V A R I A B L E S:

C       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
C                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
C       CMU    :  COMPUTATIONAL POLAR ANGLE COSINES
C       CWT    :  WEIGHTS FOR QUADRATURE OVER POLAR ANGLE COSINE
C       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT
C       NN     :  HALF THE TOTAL NUMBER OF STREAMS
C       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
C                 AT THE QUADRATURE ANGLES -CMU-
C       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

C   O U T P U T    V A R I A B L E S:

C       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5); NEEDED IN SS(15&18)
C       EVAL   :  -NN- EIGENVALUES OF EQ. SS(12) ON RETURN FROM 'ASYMTX'
C                    BUT THEN SQUARE ROOTS TAKEN
C       EVECC  :  -NN- EIGENVECTORS  (G+) - (G-)  ON RETURN
C                    FROM 'ASYMTX' ( COLUMN J CORRESPONDS TO -EVAL(J)- )
C                    BUT THEN  (G+) + (G-)  IS CALCULATED FROM SS(10),
C                    G+  AND  G-  ARE SEPARATED, AND  G+  IS STACKED ON
C                    TOP OF  G-  TO FORM -NSTR- EIGENVECTORS OF SS(7)
C       GC     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVECTORS, BUT
C                    IN AN ORDER CORRESPONDING TO -KK-
C       KK     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVALUES OF SS(7),
C                    BUT RE-ORDERED WITH NEGATIVE VALUES FIRST ( SQUARE
C                    ROOTS OF -EVAL- TAKEN AND NEGATIVES ADDED )

C   I N T E R N A L   V A R I A B L E S:

C       AMB,APB :  MATRICES (ALPHA-BETA), (ALPHA+BETA) IN REDUCED
C                    EIGENVALUE PROBLEM
C       ARRAY   :  COMPLETE COEFFICIENT MATRIX OF REDUCED EIGENVALUE
C                    PROBLEM: (ALFA+BETA)*(ALFA-BETA)
C       GPPLGM  :  (G+) + (G-) (CF. EQS. SS(10-11))
C       GPMIGM  :  (G+) - (G-) (CF. EQS. SS(10-11))
C       WK      :  SCRATCH ARRAY REQUIRED BY 'ASYMTX'
C+---------------------------------------------------------------------+
      REAL    AMB( MI,* ), APB( MI,* ), ARRAY( MI,* ), CC( MXCMU,* ),
     $     CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ), GC( MXCMU,* ),
     $     GL(0:*), KK(*), WK(*), YLMC( 0:MXCMU,* )

C                             ** CALCULATE QUANTITIES IN EQS. SS(5-6)
      DO 40 IQ  = 1, NN
         
         DO 20  JQ = 1, NSTR
            SUM = 0.0
            DO 10  L = MAZIM, NSTR-1
               SUM = SUM + GL(L) * YLMC(L,IQ) * YLMC(L,JQ)
caky               write(*,*)iq,jq,l,gl(l),ylmc(l,iq),ylmc(l,jq)
 10         CONTINUE
            CC(IQ,JQ) = 0.5 * SUM * CWT(JQ)
caky            write(*,*)iq,jq,cc(iq,jq),sum,cwt(jq)
 20      CONTINUE

         DO 30  JQ = 1, NN
C                             ** FILL REMAINDER OF ARRAY USING SYMMETRY
C                             ** RELATIONS  C(-MUI,MUJ) = C(MUI,-MUJ)
C                             ** AND        C(-MUI,-MUJ) = C(MUI,MUJ)

            CC(IQ+NN,JQ) = CC(IQ,JQ+NN)
            CC(IQ+NN,JQ+NN) = CC(IQ,JQ)
C                                      ** GET FACTORS OF COEFF. MATRIX
C                                      ** OF REDUCED EIGENVALUE PROBLEM
            ALPHA =   CC(IQ,JQ) / CMU(IQ)
            BETA = CC(IQ,JQ+NN) / CMU(IQ)
            AMB(IQ,JQ) = ALPHA - BETA
            APB(IQ,JQ) = ALPHA + BETA
 30      CONTINUE
         AMB(IQ,IQ) = AMB(IQ,IQ) - 1.0 / CMU(IQ)
         APB(IQ,IQ) = APB(IQ,IQ) - 1.0 / CMU(IQ)

 40   CONTINUE
C                      ** FINISH CALCULATION OF COEFFICIENT MATRIX OF
C                      ** REDUCED EIGENVALUE PROBLEM:  GET MATRIX
C                      ** PRODUCT (ALFA+BETA)*(ALFA-BETA); SS(12)
      DO 70  IQ = 1, NN
         DO 70  JQ = 1, NN
            SUM = 0.
            DO 60  KQ = 1, NN
               SUM = SUM + APB(IQ,KQ) * AMB(KQ,JQ)
 60         CONTINUE
            ARRAY(IQ,JQ) = SUM
70    CONTINUE
C                      ** FIND (REAL) EIGENVALUES AND EIGENVECTORS


      CALL  ASYMTXs( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WK )
      
      IF ( IER.GT.0 )  THEN
         WRITE( *, '(//,A,I4,A)' )  ' ASYMTX--EIGENVALUE NO. ', IER,
     $        '  DIDNT CONVERGE.  LOWER-NUMBERED EIGENVALUES WRONG.'
         CALL  ERRMSG( 'ASYMTX--CONVERGENCE PROBLEMS', .TRUE. )
      END IF

CDIR$ IVDEP
      DO 75  IQ = 1, NN
         EVAL(IQ) = SQRT( ABS( EVAL(IQ) ) )
         KK( IQ+NN ) = EVAL(IQ)
C                                             ** ADD NEGATIVE EIGENVALUE
         KK( NN+1-IQ ) = - EVAL(IQ)
75    CONTINUE
C                          ** FIND EIGENVECTORS (G+) + (G-) FROM SS(10)
C                          ** AND STORE TEMPORARILY IN -APB- ARRAY
      DO 90  JQ = 1, NN
         DO 90  IQ = 1, NN
            SUM = 0.
            DO 80  KQ = 1,NN
               SUM = SUM + AMB(IQ,KQ) * EVECC(KQ,JQ)
 80         CONTINUE
            APB(IQ,JQ) = SUM / EVAL(JQ)
 90   CONTINUE

      DO 100  JQ = 1, NN
CDIR$ IVDEP
         DO 100  IQ = 1, NN
            GPPLGM = APB(IQ,JQ)
            GPMIGM = EVECC(IQ,JQ)
C                                ** RECOVER EIGENVECTORS G+,G- FROM
C                                   THEIR SUM AND DIFFERENCE; STACK THEM
C                                   TO GET EIGENVECTORS OF FULL SYSTEM
C                                   SS(7) (JQ = EIGENVECTOR NUMBER)

            EVECC(IQ,      JQ) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,   JQ) = 0.5 * ( GPPLGM - GPMIGM )

C                                ** EIGENVECTORS CORRESPONDING TO
C                                ** NEGATIVE EIGENVALUES (CORRESP. TO
C                                ** REVERSING SIGN OF 'K' IN SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )
100   CONTINUE

      RETURN
      END
      SUBROUTINE  ASYMTXs( AA, EVEC, EVAL, M, IA, IEVEC, IER, WK )

C    =======  S I N G L E    P R E C I S I O N    V E R S I O N  ======

C       SOLVES EIGENFUNCTION PROBLEM FOR REAL ASYMMETRIC MATRIX
C       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE REAL.

C           ( C O M M E N T S    O M I T T E D )
C+---------------------------------------------------------------------+

      REAL              AA( IA,* ),  WK(*),  EVAL(*), EVEC( IEVEC,* )
      LOGICAL           NOCONV, NOTLAS
      DATA     C1/ 0.4375 /, C2/ 0.5 /, C3/ 0.75 /, C4/ 0.95 /,
     $     C5/ 16.0 /, C6/ 256.0 /, ZERO / 0.0 /, ONE / 1.0 /


      IER = 0
      TOL = R1MACH(3)
      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )
     $     CALL ERRMSG( 'ASYMTX--BAD INPUT VARIABLE(S)', .TRUE. )

C                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVAL(1) = AA(1,1)
         EVEC(1,1) = ONE
         RETURN

      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AA(1,1) - AA(2,2) )**2 + 4. * AA(1,2) * AA(2,1)
         IF ( DISCRI.LT.ZERO )
     $        CALL ERRMSG( 'ASYMTX--COMPLEX EVALS IN 2X2 CASE', .TRUE. )
         SGN = ONE
         IF ( AA(1,1).LT.AA(2,2) )  SGN = - ONE
         EVAL(1) = 0.5 * ( AA(1,1) + AA(2,2) + SGN*SQRT(DISCRI) )
         EVAL(2) = 0.5 * ( AA(1,1) + AA(2,2) - SGN*SQRT(DISCRI) )
         EVEC(1,1) = ONE
         EVEC(2,2) = ONE
         IF ( AA(1,1).EQ.AA(2,2) .AND. 
     $        (AA(2,1).EQ.ZERO.OR.AA(1,2).EQ.ZERO) )  THEN
            RNORM =   ABS(AA(1,1)) + ABS(AA(1,2)) + ABS(AA(2,1))
     $           + ABS(AA(2,2))
            W = TOL * RNORM
            EVEC(2,1) =   AA(2,1) / W
            EVEC(1,2) = - AA(1,2) / W
         ELSE
            EVEC(2,1) = AA(2,1) / ( EVAL(1) - AA(2,2) )
            EVEC(1,2) = AA(1,2) / ( EVAL(2) - AA(1,1) )
         ENDIF
         
         RETURN

	END IF
C                                        ** INITIALIZE OUTPUT VARIABLES
	IER = 0
	DO 20 I = 1, M
	   EVAL(I) = ZERO
	   DO 10 J = 1, M
	      EVEC(I,J) = ZERO
10       CONTINUE
         EVEC(I,I) = ONE
20    CONTINUE
C                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
C                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
C                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M

30    KKK = K
      DO 70  J = KKK, 1, -1
         ROW = ZERO
         DO 40 I = 1, K
            IF ( I.NE.J ) ROW = ROW + ABS( AA(J,I) )
 40      CONTINUE
         IF ( ROW.EQ.ZERO ) THEN
            WK(K) = J
            IF ( J.NE.K ) THEN
               DO 50 I = 1, K
                  REPL   = AA(I,J)
                  AA(I,J) = AA(I,K)
                  AA(I,K) = REPL
50             CONTINUE
               DO 60 I = L, M
                  REPL   = AA(J,I)
                  AA(J,I) = AA(K,I)
                  AA(K,I) = REPL
60             CONTINUE
            END IF
            K = K - 1
            GO TO 30
         END IF
70    CONTINUE
C                                     ** SEARCH FOR COLUMNS ISOLATING AN
C                                       ** EIGENVALUE AND PUSH THEM LEFT
 80   LLL = L
      DO 120 J = LLL, K
         COL = ZERO
         DO 90 I = L, K
            IF ( I.NE.J ) COL = COL + ABS( AA(I,J) )
90       CONTINUE
         IF ( COL.EQ.ZERO ) THEN
            WK(L) = J
            IF ( J.NE.L ) THEN
               DO 100 I = 1, K
                  REPL   = AA(I,J)
                  AA(I,J) = AA(I,L)
                  AA(I,L) = REPL
100            CONTINUE
               DO 110 I = L, M
                  REPL   = AA(J,I)
                  AA(J,I) = AA(L,I)
                  AA(L,I) = REPL
110            CONTINUE
            END IF
            L = L + 1
            GO TO 80
         END IF
 120  CONTINUE
C                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WK(I) = ONE
130   CONTINUE

140   NOCONV = .FALSE.
      DO 200 I = L, K
         COL = ZERO
         ROW = ZERO
         DO 150 J = L, K
            IF ( J.NE.I ) THEN
               COL = COL + ABS( AA(J,I) )
               ROW = ROW + ABS( AA(I,J) )
            END IF
 150     CONTINUE
         F = ONE
         G = ROW / C5
         H = COL + ROW
 160     IF ( COL.LT.G ) THEN
            F   = F * C5
            COL = COL * C6
            GO TO 160
         END IF
         G = ROW * C5
 170     IF ( COL.GE.G ) THEN
            F   = F / C5
            COL = COL / C6
            GO TO 170
         END IF
C                                                         ** NOW BALANCE
         IF ( (COL+ROW) / F .LT. C4 * H ) THEN
            WK(I)  = WK(I) * F
            NOCONV = .TRUE.
            DO 180 J = L, M
               AA(I,J) = AA(I,J) / F
 180        CONTINUE
            DO 190 J = 1, K
               AA(J,I) = AA(J,I) * F
 190        CONTINUE
         END IF
 200  CONTINUE

      IF ( NOCONV ) GO TO 140
C                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1.LT.L+1 ) GO TO 350
C                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H       = ZERO
         WK(N+M) = ZERO
         SCALE   = ZERO
C                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + ABS(AA(I,N-1))
 210     CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WK(I+M) = AA(I,N-1) / SCALE
               H = H + WK(I+M) * WK(I+M)
 220        CONTINUE
            G = - SIGN( SQRT(H),WK(N+M) )
            H = H - WK(N+M) * G
            WK(N+M) = WK(N+M) - G
C                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WK(I+M) * AA(I,J)
 230           CONTINUE
               DO 240 I = N, K
                  AA(I,J) = AA(I,J) - WK(I+M) * F / H
 240           CONTINUE
 250        CONTINUE
C                                  ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WK(J+M) * AA(I,J)
 260           CONTINUE
               DO 270 J = N, K
                  AA(I,J) = AA(I,J) - WK(J+M) * F / H
 270           CONTINUE
 280        CONTINUE
            WK(N+M)  = SCALE * WK(N+M)
            AA(N,N-1) = SCALE * G
         END IF
 290  CONTINUE

      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AA(N+1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WK(N+1+M)
            DO 300 I = N+2, K
               WK(I+M) = AA(I,N)
 300        CONTINUE
            IF ( N+1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N+1, K
                     G = G + WK(I+M) * EVEC(I,J)
 310              CONTINUE
                  G = G / F
                  DO 320 I = N+1, K
                     EVEC(I,J) = EVEC(I,J) + G * WK(I+M)
 320              CONTINUE
 330           CONTINUE
            END IF
         END IF
 340  CONTINUE

 350  CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + ABS(AA(I,J))
 360     CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVAL(I) = AA(I,I)
 370  CONTINUE
      N = K
      T = ZERO
C                            ** SEARCH FOR NEXT EIGENVALUES
 380  IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
C                        ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
 390  CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = ABS( AA(LB-1,LB-1)) + ABS(AA(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( ABS(AA(LB,LB-1)) .LE. TOL * S ) GO TO 410
 400  CONTINUE
C                                                      ** ONE EVAL FOUND
 410  X = AA(N,N)
      IF ( LB.EQ.N  ) THEN
         AA(N,N)  = X + T
         EVAL(N) = AA(N,N)
         N = N1
         GO TO 380
      END IF
C                                                     ** TWO EVALS FOUND
      Y = AA(N1,N1)
      W = AA(N,N1) * AA(N1,N)
      IF ( LB.EQ.N1 ) THEN
         P = (Y-X) * C2
         Q = P * P + W
         Z = SQRT(ABS(Q))
         AA(N,N) = X + T
         X = AA(N,N)
         AA(N1,N1) = Y + T
C                                                           ** REAL PAIR
         Z = P + SIGN(Z,P)
         EVAL(N1) = X + Z
         EVAL(N)  = EVAL(N1)
         IF ( Z.NE.ZERO ) EVAL(N) = X - W / Z
         X = AA(N,N1)
C                                  ** EMPLOY SCALE FACTOR IN CASE
C                                  ** X AND Z ARE VERY SMALL
         R = SQRT(X * X + Z * Z)
         P = X / R
         Q = Z / R
C                                                    ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AA(N1,J)
            AA(N1,J) = Q * Z + P * AA(N,J)
            AA(N,J)  = Q * AA(N,J) - P * Z
 420     CONTINUE
C                                                 ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AA(I,N1)
            AA(I,N1) = Q * Z + P * AA(I,N)
            AA(I,N)  = Q * AA(I,N) - P * Z
 430     CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVEC(I,N1)
            EVEC(I,N1) = Q * Z + P * EVEC(I,N)
            EVEC(I,N)  = Q * EVEC(I,N) - P * Z
 440     CONTINUE
         N = N2
         GO TO 380
      END IF
C                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
C                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE

      IF ( IN.EQ.30 ) THEN
         IER = 128 + N
         RETURN
      END IF
C                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AA(I,I) = AA(I,I) - X
 450     CONTINUE
         S = ABS(AA(N,N1)) + ABS(AA(N1,N2))
         X = C3 * S
         Y = X
         W = -C1 * S * S
      END IF

      IN = IN + 1
C                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AA(I,I)
         R = X - Z
         S = Y - Z
         P = (R * S-W) / AA(I+1,I) + AA(I,I+1)
         Q = AA(I+1,I+1) - Z - R - S
         R = AA(I+2,I+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = ABS( AA(I,I-1) ) * ( ABS(Q) + ABS(R) )
         VV = ABS(P) * ( ABS(AA(I-1,I-1)) + ABS(Z) + ABS(AA(I+1,I+1)) )
         IF ( UU .LE. TOL*VV ) GO TO 470
 460  CONTINUE

 470  CONTINUE
      AA(I+2,I) = ZERO
      DO 480 J = I+3, N
         AA(J,J-2) = ZERO
         AA(J,J-3) = ZERO
 480  CONTINUE

C             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N
      
      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( SQRT( P*P + Q*Q + R*R ),P )
            IF ( LB.NE.I ) AA(KA,KA-1) = - AA(KA,KA-1)
         ELSE
            P = AA(KA,KA-1)
            Q = AA(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AA(KA+2,KA-1)
            X = ABS(P) + ABS(Q) + ABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( SQRT( P*P + Q*Q + R*R ),P )
            AA(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
C                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AA(KA,J) + Q * AA(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AA(KA+2,J)
               AA(KA+2,J) = AA(KA+2,J) - P * Z
            END IF
            AA(KA+1,J) = AA(KA+1,J) - P * Y
            AA(KA,J)   = AA(KA,J)   - P * X
 490     CONTINUE
C                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AA(II,KA) + Y * AA(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AA(II,KA+2)
               AA(II,KA+2) = AA(II,KA+2) - P * R
            END IF
            AA(II,KA+1) = AA(II,KA+1) - P * Q
            AA(II,KA)   = AA(II,KA) - P
 500     CONTINUE
C                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVEC(II,KA) + Y * EVEC(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVEC(II,KA+2)
               EVEC(II,KA+2) = EVEC(II,KA+2) - P * R
            END IF
            EVEC(II,KA+1) = EVEC(II,KA+1) - P * Q
            EVEC(II,KA)   = EVEC(II,KA) - P
 510     CONTINUE
 520  CONTINUE
      GO TO 390
C                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
 530  CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AA(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AA(I,I) - EVAL(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AA(I,N)
               DO 540 J = N2, N-1
                  R = R + AA(I,J) * AA(J,N)
 540           CONTINUE
               AA(I,N) = -R / W
               N2 = I
 550        CONTINUE
 560     CONTINUE
C                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVEC(I,J) = AA(I,J)
 570           CONTINUE
            END IF
 580     CONTINUE
C                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 610  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVEC(I,N) * AA(N,J)
 590              CONTINUE
                  EVEC(I,J) = Z
 600           CONTINUE
 610        CONTINUE
         END IF
         
      END IF

      DO 620 I = L, K
         DO 620 J = 1, M
            EVEC(I,J) = EVEC(I,J) * WK(I)
 620  CONTINUE
C                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = WK(I)
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL      = EVEC(I,N)
               EVEC(I,N) = EVEC(J,N)
               EVEC(J,N) = REPL
 630        CONTINUE
         END IF
 640  CONTINUE
      
      DO 660 I = K+1, M
         J = WK(I)
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL      = EVEC(I,N)
               EVEC(I,N) = EVEC(J,N)
               EVEC(J,N) = REPL
 650        CONTINUE
         END IF
 660  CONTINUE
      
      RETURN
      END
