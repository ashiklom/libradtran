c  Wrapper function for REFICE, to make it accessible from C
c
      SUBROUTINE WREFICE ( WAVLEN, TEMP, RE, IM )

      IMPLICIT  NONE


      REAL WAVLEN, TEMP, RE, IM
      COMPLEX REF

      COMPLEX REFICE
      EXTERNAL REFICE

      REF = REFICE( WAVLEN, TEMP )

      RE = REAL(REAL(REF))
      IM = REAL(AIMAG(REF))

      RETURN
      END
