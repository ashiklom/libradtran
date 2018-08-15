      REAL FUNCTION eval ( n, x, a, t)
*
* Newton interpolation
*   
*   Reference: Cheney and Kincaid, Numerical Mathematics and 
*              Computing, 1980, Brooks/Cole Publishing Company
*              Monterey, California. */
*
      INTEGER n
      REAL x(*), a(*), t
*
      eval = a(n)
      DO i = n-1, 1, -1
         eval = eval * ( t - x(i) ) + a(i)
      ENDDO
*      
      RETURN
      END
