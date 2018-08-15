* Front end for plkavg because C doesn't like Fortran functions
      SUBROUTINE CPLKAVG( WNUMLO, WNUMHI, T, RAD )

      RAD = PLKAVG( WNUMLO, WNUMHI, T )
      RETURN
      END
