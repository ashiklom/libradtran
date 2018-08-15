*
      SUBROUTINE linpol( x1, y1, x2, y2, x, y )
*--------------------------------------------------------------------
* Linearly interpolate between (x1,y1) and (x2,y2), the new value
* at x is y.
*
      REAL a, b, x1, y1, x2, y2, x, y
      a = ( y2-y1 ) / ( x2 -x1 )
      b = y1 - a*x1
      y = a*x + b
      RETURN
      END
