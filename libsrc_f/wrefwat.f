c  Wrapper function for REFWAT, to make it accessible from C
c
      SUBROUTINE WREFWAT ( WAVLEN, TEMP, RE, IM )

      IMPLICIT  NONE


      REAL WAVLEN, TEMP, RE, IM
      COMPLEX REF

      COMPLEX REFWAT
      EXTERNAL REFWAT
      
      REF = REFWAT( WAVLEN, TEMP )

      RE = REAL(REAL(REF))
      IM = REAL(AIMAG(REF))

      RETURN
      END
