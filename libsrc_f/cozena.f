      REAL FUNCTION  cozena( day, hour, dlat, dlon )
C+---------------------------------------------------------------------+
C!       Reference:    Woolf, H.M., NASA TM X-1646.                    !
C!                     On the computation of solar elevation           !
C!                     angles and the determination of sunrise         !
C!                     and sunset times.                               !
C+---------------------------------------------------------------------+
C!       D a t a               v a r i a b l e s:                      !
C+---------------------------------------------------------------------+
C!       epsiln   eccentricity of earth orbit                          !
C!       sinob    sine of obliquity of ecliptic                        !
C!       dpy      days per year (365.242)                              !
C!       dph      degrees per hour (360./24)                           !
C+---------------------------------------------------------------------+
C!       I n p u t             v a r i a b l e s:                      !
C+---------------------------------------------------------------------+
C!       day      days of year (jan. 1 is day 1)                       !
C!       hour     hours of greenwich time                              !
C!       dlon     degree longitude of grid point,                      !
C!                counted positive west of greenwich                   !
C!       dlat     degree latitude  of grid point,                      !
C!                counted positive northern hemisphere                 !
C+---------------------------------------------------------------------+
C!       O u t p u t           v a r i a b l e s:                      !
C+---------------------------------------------------------------------+
C!       cozena   cosine of zenith angle   ( reference -- eq. 1.1 )    !
C+---------------------------------------------------------------------+
C!       I n t e r n a l       v a r i a b l e s:                      !
C+---------------------------------------------------------------------+
C!       dpr      degree/radian (*57.29578)                            !
C!       rpd      radian/degree (*.01745329)                           !
C!       dang     angle measured from perihelion, in radians           !
C!                which is taken as midnight jan. 1                    !
C!       homp     hours of meridian passage or true solar noon         !
C!                ( reference -- eq. 1.6 )                             !
C!       hang     hour angle, a measure of the longitudinal distance   !
C!                to the sun from the point calculated  ( eq. 1.5 )    !
C!       sindlt   sine of declination angle  ( eq. 1.2 )               !
C!       cosdlt   cosine of declination angle                          !
C!       sigma    reference -- eq. 1.3a                                !
C!       ang      reference -- eq. 1.3b                                !
C+---------------------------------------------------------------------+
      REAL ang, cosdlt, dang, day, dlat,
     $     dlon, dph, dpy, epsiln, 
     $     hang, homp, hour, pi, rpd, dpr, sigma, sindlt, sinob
      DATA   epsiln/.016733/, sinob/.3978/, dpy/365.242/,
     $       dph/15.0/, pi/3.1415926535898/
      rpd    = pi/180.0
      dpr    = 1.0/rpd
      dang   = 2.0*pi*(day - 1.0)/dpy
      homp   = 12.0 + 0.123570*SIN( dang )
     $     - 0.004289*COS( dang )
     $     + 0.153809*SIN( 2.0*dang )
     $     + 0.060783*COS( 2.0*dang )
      hang   = dph*( hour-homp ) - dlon
      ang    = 279.9348*rpd + dang
      sigma  = ( ang*dpr + 0.4087*SIN( ang )
     $     + 1.8724*COS( ang )
     $     - 0.0182*SIN( 2.0*ang )
     $     + 0.0083*COS( 2.0*ang ) )*rpd
      sindlt = sinob*SIN( sigma )
      cosdlt = SQRT( 1.0-sindlt**2 )
      cozena = sindlt*SIN( rpd*dlat ) +
     $         cosdlt*COS( rpd*dlat )*COS( rpd*hang )
      RETURN
      END
