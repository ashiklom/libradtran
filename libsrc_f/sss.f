      SUBROUTINE SSS (NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
     &                WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
     &                PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
     &                MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
     &                FLUP, DFDT, UAVG, UU, ALBMED, TRNMED,
     &                uavgdn, uavgso, uavgup, 
     &                btype, brho0, bk, btheta, bu10, bpcl, bsal )

      IMPLICIT NONE

      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI,
     &          MXSQT, MXSZA
      INCLUDE  'DISORT.MXD'
      PARAMETER ( MI = MXCMU / 2, MI9M2 = 9*MI - 2,
     &          NNLYRI = MXCMU*MXCLY, MXSQT = 1000 )
c     ..
c     .. Scalar Arguments ..

      CHARACTER HEADER*127
      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, NLYR,
     &          NMOM, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
      INTEGER   btype
      REAL      brho0, bk, btheta, bu10, bpcl, bsal
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( 5 )
      REAL      ALBMED( MAXUMU ), DFDT( MAXULV ), DTAUC( MAXCLY ),
     &          FLUP( MAXULV ), PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ),
     &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( MAXCLY ),
     &          TEMPER( 0:MAXCLY ), TRNMED( MAXUMU ), UAVG( MAXULV ),
     &          UMU( MAXUMU ), UTAU( MAXULV ),
     &          UU( MAXUMU, MAXULV, MAXPHI )
      REAL      UAVGDN( MAXULV ), UAVGUP( MAXULV ), UAVGSO( MAXULV )


c     .. local variables
      INTEGER LC, IU, J, K
      DOUBLE PRECISION UTOA, G, ALB, MUREL
      DOUBLE PRECISION PL, PLM1, PLM2
      DOUBLE PRECISION P(MXCLY), TAUC( 0:MXCLY )

      TAUC(0) = 0.0

c     .. currently only calculation for TOA up-welling is 
c     .. implemented; other altitudes and down-welling is 
c     .. straightforward to implement but takes some time
      IF (NTAU .GT. 1) THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ',
     $        'only top of atmosphere calculation implemented'
         STOP
      ENDIF

      IF (UTAU(1) .GT. 0.0) THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ',
     $        'only top of atmosphere calculation implemented'
         WRITE ( *, * )  ' UTAU(1) = ', UTAU(1)
         STOP
      ENDIF


c     .. loop over user zenith
      DO IU = 1, NUMU
         IF (UMU (IU) .LT. 0.0) THEN
            WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ',
     $           'UMU smaller than 0 (down-welling) not implemented'
            STOP
         END IF

c     .. loop over user azimuth
         DO J = 1, NPHI

c     .. calculate scattering angle
            MUREL = COS((PHI(J) - PHI0)/180.0*3.1415926 )
     $           * SQRT((1.0 - UMU(IU)*UMU(IU))*(1.0 - UMU0*UMU0))
     $           - UMU(IU) * UMU0
            
            G = 1.0 / UMU0 + 1.0 / UMU (IU)
            
c     .. calculate scattering phase function at scattering angle
            DO LC = 1, NLYR
               PLM2 = 1.0D0
               PLM1 = MUREL
               
               P(LC) = PMOM(0,LC) * PLM2
               
               IF (NMOM .GT. 0)  P(LC) = P(LC) + PMOM(1,LC)*PLM1*3.0D0
               
               DO K = 2, NMOM
                  PL = (DBLE(2*K-1) * MUREL * PLM1 - 
     $                 DBLE(K-1) * PLM2) / DBLE(K)
                  
                  P(LC) = P(LC) + PMOM(K,LC)*PL*DBLE(2*K+1)
                  
                  PLM2 = PLM1
                  PLM1 = PL
               ENDDO
            ENDDO
            
            UTOA = 0.0
c     .. calculate atmospheric contribution
            DO LC = 1, NLYR
               
               TAUC( LC ) = TAUC( LC - 1 ) + DTAUC( LC )
               
               ALB = SSALB(LC) * P(LC) / 4.0D0 / 3.1415926D0 / G
     $              * (1.0 - EXP (-DTAUC(LC) * G)) / UMU (IU)
               
               UTOA = UTOA + EXP(-TAUC( LC-1 ) * G ) * ALB
            ENDDO
            
c     .. add contribution of the surface 
            UTOA = UTOA + UMU0 / 3.1415926D0 * 
     $           EXP(-TAUC(NLYR) * G ) * ALBEDO
            
c     .. multiply with direct beam source
            UU( IU, 1, J ) = UTOA * FBEAM
         ENDDO
      ENDDO

      RETURN
      END
