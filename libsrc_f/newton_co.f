      SUBROUTINE newton_co( np, x, y, a )
*
* Newton interpolation
*   
*   Reference: Cheney and Kincaid, Numerical Mathematics and 
*              Computing, 1980, Brooks/Cole Publishing Company
*              Monterey, California. */
*
*
      INTEGER np
      REAL x(*), y(*), a(*)
*
      INTEGER i, j
*
      DO i = 1, np
         a(i) = y(i)
      ENDDO
*
      DO j=1, np-1
         DO i=np,j+1,-1
            a(i) = (a(i) - a(i-1)) / (x(i) - x(i-j))
         ENDDO
      ENDDO
      RETURN
      END
