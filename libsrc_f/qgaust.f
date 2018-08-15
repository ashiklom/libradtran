      SUBROUTINE QGAUST( NSTR, SZA, UMU0, RESULT)

c       Test if the solar zenith angle coincides with one of
c       the computational angles

      INCLUDE  'DISORT.MXD'
      INTEGER NSTR, NN, RESULT
      REAL UMU0, UMU0S, SZA, PI, CMU (MXCMU), CWT (MXCMU)

      PI = 2.*ASIN( 1.0 )

      RESULT=0

      NN = NSTR / 2

      CALL QGAUSN2( NN, CMU, CWT )

      DO IQ=1, NN
         IF( ABS( (UMU0 - CMU (IQ)) / UMU0 ) .LT.1.E-4 ) THEN
            UMU0S = UMU0
            IF (UMU0 .LT. CMU (IQ)) THEN
               UMU0  = CMU (IQ) * (1. - 1.1E-4)
            ELSE
               UMU0  = CMU (IQ) * (1. + 1.1E-4)
            ENDIF

            SZA   = 180./PI*ACOS (UMU0)
            WRITE (0,*) 
     &           '******* WARNING >>>>>> ',
     &           'SETDIS--beam angle=computational angle;',
     &           '******* changing solar zenith angle umu0 from', 
     &           UMU0S, ' to', UMU0
            RESULT=-1
         ENDIF
      ENDDO
      
      RETURN
      END
