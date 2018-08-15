c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /home/users/buras/cvs2svn-2.3.0/cvs/libRadtran/libsrc_f/BDREF.f,v 1.16 2011-02-24 19:14:22 bernhard Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      REAL FUNCTION RPV_REFLECTION (RHO0, K, THETA, SIGMA, T1, T2,
     &     SCALE, MU1EXT, MU2EXT, PHIEXT, BADMU)
C       Computes the Rahman, Pinty, Verstraete BRDF.  The incident
C     and outgoing cosine zenith angles are MU1 and MU2, respectively,
C     and the relative azimuthal angle is PHI.  In this case the incident
C     direction is where the radiation is coming from, so MU1>0 and 
C     the hot spot is MU2=MU1 and PHI=180 (the azimuth convention is
C     different from the original Frank Evans code). 
C     The reference is:
C       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
C       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
C       With NOAA Advanced Very High Resolution Radiometer Data,
C       J. Geophys. Res., 98, 20791-20801.
      IMPLICIT  NONE

      REAL RHO0, K, THETA, SIGMA, T1, T2, SCALE, MU1EXT, MU2EXT
      REAL PHIEXT, BADMU

      REAL M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG
      REAL HSPOT, T, G
      REAL MU1, MU2, PHI

C     Work on copies
      MU1 = MU1EXT
      MU2 = MU2EXT
      PHI = PHIEXT

C     Azimuth convention different from Frank Evans:
C     Here PHI=0 means the backward direction while 
C     while in DISORT PHI=0 means forward.
      PHI = 3.1415926 - PHI

C     Don't allow mu's smaller than BADMU because 
C     the albedo is larger than 1 for those 
      IF (BADMU .GT. 0) THEN
         IF (MU1 .LT. BADMU) THEN
            MU1 = BADMU
         ENDIF
         
         IF (MU2 .LT. BADMU) THEN
            MU2 = BADMU
         ENDIF
      ENDIF

C     Special treatment of hot spot
      IF ( PHI .EQ. 0.0 .AND. MU1 .EQ. MU2) THEN
         HSPOT = RHO0 * ((2*MU1*MU1*MU1)**(K-1) * 
     &        (1-THETA)/(1+THETA)/(1+THETA) * (2-RHO0) 
     &        + SIGMA/MU1) * (T1*EXP(3.1415926*T2)+1.0)
         
         RPV_REFLECTION = HSPOT * SCALE
         RETURN
      ENDIF

      M = (MU1 * MU2 * (MU1 + MU2))**(K-1)
      COSPHI = COS(PHI)
      SIN1 = SQRT(1.0-MU1**2)
      SIN2 = SQRT(1.0-MU2**2)
      COSG = MU1*MU2 + SIN1*SIN2*COSPHI
      G = ACOS(COSG)
      F = (1-THETA**2) / (1 + 2*THETA*COSG + THETA**2)**1.5

      TAN1 = SIN1/MU1
      TAN2 = SIN2/MU2
      CAPG = SQRT( TAN1**2 + TAN2**2 - 2*TAN1*TAN2*COSPHI )
      H = 1 + (1-RHO0)/(1+CAPG)
      T = 1 + T1*EXP(T2*(3.1415926-G))

      RPV_REFLECTION = RHO0 * (M * F * H + SIGMA/MU1) * T * SCALE
      
      IF (RPV_REFLECTION .LE. 0.0) THEN
         RPV_REFLECTION = 0.0
      ENDIF

      RETURN
      END 


      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI,
     &                      type, 
     &                      rho0, k, theta, sigma, t1, t2, scale,
     &                      iso, vol, geo,
     &                      u10, pcl, xsal,
     &                      callnum, FAST_SURFACE )
c
c      Supplies surface bi-directional reflectivity.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c
c   Called by- DREF, SURFAC

c +-------------------------------------------------------------------+
c
      IMPLICIT  NONE

c     .. Scalar Arguments ..
      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
      INTEGER   TYPE
      REAL      RHO0, K, THETA, SIGMA, T1, T2, SCALE
      REAL      ISO, VOL, GEO
      REAL      U10, PCL, XSAL
      INTEGER   CALLNUM
c     .. Internal variables
c     REAL CR, CI, 
      LOGICAL INPERR, FAST_SURFACE
      INTEGER IRMU
      REAL WAVLEN, TEMP, FLXALB, RMU, BADMU
      REAL      SWVNMHI, SWVNMLO, SRHO0, SK, STHETA
      REAL      SSIGMA, ST1, ST2, SSCALE
      REAL      SISO, SVOL, SGEO
      REAL      SMU, SMUP, SDPHI
      
c     .. Function
      REAL RPV_REFLECTION, DREF2
      Logical WRTBAD2

      DATA BADMU    /-999.0/
      SAVE BADMU
      SAVE SWVNMHI, SWVNMLO, SRHO0, SK, STHETA
      SAVE SSIGMA, ST1, ST2, SSCALE
      SAVE SISO, SVOL, SGEO

c     ??? need to adapt FIRSTRPV as in ocean.c
c     ??? it's not just the first call - it's any call after a change
c
c     wavelength = 1 / average wavenumber
      WAVLEN = 2.E4 / (WVNMLO + WVNMHI)
      TEMP = 273.

      if (TYPE .EQ. 1) THEN

c     Apply DISORT2 test; use exactly the same angular grids etc.
c     Determine the smallest angle for which the albedo is lower 
c     than 1 and use this angle as lower limit in the following 
c     BDREF calls.
         if ( (SWVNMLO .NE. WVNMLO) .OR. 
     &        (SWVNMHI .NE. WVNMHI) .OR. 
     &        (SRHO0   .NE. RHO0)   .OR.
     &        (SK      .NE. K)      .OR.
     &        (STHETA  .NE. THETA)  .OR.
     &        (SSIGMA  .NE. SIGMA)  .OR.
     &        (ST1     .NE. T1)     .OR.
     &        (ST2     .NE. T2)     .OR.
     &        (SSCALE  .NE. SCALE)) THEN

            SWVNMLO = WVNMLO
            SWVNMHI = WVNMHI
            SRHO0   = RHO0
            SK      = K
            STHETA  = THETA
            SSIGMA  = SIGMA
            ST1     = T1
            ST2     = T2
            SSCALE  = SCALE

            BADMU = 0.0

            DO IRMU = 100,0,-1

               RMU  = IRMU*0.01
               
               FLXALB = DREF2( WVNMLO, WVNMHI, RMU, 
     &              TYPE, RHO0, K, THETA, 
     &              SIGMA, T1, T2,
     &              SCALE, ISO, VOL, GEO, U10, PCL, XSAL, CALLNUM,
     &              FAST_SURFACE )
               
               IF( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 ) THEN
                  BADMU = (IRMU+1)*0.01
                  IF (BADMU .GT. 1.0)  BADMU = 1.0
                  write (0,*) '*** Using',BADMU,
     &                 ' as limiting mu in AMBRALS'
                  GOTO 80
               ENDIF
               
            ENDDO
            
 80         CONTINUE

         ENDIF

         BDREF = RPV_REFLECTION (RHO0, K, THETA, SIGMA, T1, T2, 
     &        SCALE, MUP, MU, DPHI, BADMU)
      ENDIF

         


      if (TYPE .EQ. 2) THEN

c    call C wrapper function
         call oceanfort( WVNMLO, WVNMHI, MU, MUP, DPHI, 
     &        U10, PCL, XSAL, CALLNUM, BDREF )

c     Remove BRDFs smaller than 0
         if (BDREF .LT. 0)  BDREF = 0

c     Check for NaN
         if (BDREF .NE. BDREF) THEN
            write (0,*) 'BDREF', WVNMLO, WVNMHI, MU, MUP, DPHI, BDREF 
            BDREF = 1
         ENDIF
      ENDIF


      if (TYPE .EQ. 3) THEN

c     mu=0 or dmu=0 cause problems
         if ( (SISO .NE. ISO) .OR. 
     &        (SVOL .NE. VOL) .OR. 
     &        (SGEO .NE. GEO)) THEN

            SISO = ISO
            SVOL = VOL
            SGEO = GEO

            BADMU = 0.0

            DO IRMU = 100,0,-1

               RMU  = IRMU*0.01
               
               FLXALB = DREF2( WVNMLO, WVNMHI, RMU, 
     &              TYPE, RHO0, K, THETA, 
     &              SIGMA, T1, T2,
     &              SCALE, ISO, VOL, GEO, U10, PCL, XSAL, CALLNUM,
     &              FAST_SURFACE )
               
               IF( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 ) THEN
                  BADMU = (IRMU+1)*0.01
                  IF (BADMU .GT. 1.0)  BADMU = 1.0
                  write (0,*) '*** Using',BADMU,
     &                 ' as limiting mu in AMBRALS'
                  GOTO 90
               ENDIF
               
            ENDDO
            
 90         CONTINUE

         ENDIF

c     convert phi to degrees
         SDPHI = DPHI
         SMUP = MUP
         SMU = MU

         dphi = dphi / 3.1415927 * 180.0

         IF (BADMU .GT. 0) THEN
            IF (MU .LT. BADMU) THEN
               MU = BADMU
            ENDIF
            
            IF (MUP .LT. BADMU) THEN
               MUP = BADMU
            ENDIF
         ENDIF

         call ambralsfort( iso, vol, geo, mu, mup, dphi, BDREF )
      
         DPHI = SDPHI
         MUP = SMUP
         MU = SMU

c     Check for NaN
         if (BDREF .NE. BDREF) THEN
            write (0,*) 'BDREF', iso, vol, geo, MU, MUP, DPHI, BDREF 
            BDREF = 1
         ENDIF
      ENDIF

      RETURN
      END
