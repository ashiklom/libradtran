       SUBROUTINE inter( dim, npoints, itype, arg,
     $     xarr, yarr, ynew, hh)
*-------------------------------------------------------------------
* Copyright (C) 1994 Arve Kylling
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 1, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* To obtain a copy of the GNU General Public License write to the
* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
* USA.
*-------------------------------------------------------------------
*
*     Interpolates at the x-point arg from x-value array xarr and
*     y-value array yarr. xarr and yarr are expected to have
*     descending arguments, i.e. for atmospheric applications
*     xarr typically holds the altitude and inter expects
*     xarr(1) = altitude at top of atmosphere.
*
*     Input variables:
*     dim       Array dimension of xarr and yarr
*     npoints   No. points in arrays xarr and yarr
*     itype     Interpolation type
*     arg       Interpolation argument
*     xarr      array of x values
*     yarr      array of y values
*
*     Output variables:
*     ynew      Interpolated function value at arg
*     hh        gradient or scale height value  
*
      INTEGER dim
      REAL  arg, hh, xarr(dim), yarr(dim), ynew
*
       IF ( arg .LE. xarr(1) .AND. arg .GE. xarr(npoints)) THEN
         DO 10 iq = 1 , npoints-1
           IF ( arg .LE. xarr(iq) .AND. arg .GT. xarr(iq+1)) ip = iq
   10    CONTINUE
         IF ( arg .EQ. xarr(npoints)) ip = npoints - 1
       ELSEIF ( arg .GT. xarr(1)) THEN
         ip = 1
       ELSEIF ( arg .LT. xarr(npoints)) THEN
         ip = npoints - 1
       ENDIF
*
*     Interpolate function value at arg from data points ip to ip+1
*
*
*     exponential interpolation
*
       IF ( itype .EQ. 1 ) THEN
          IF ( yarr(ip+1) .EQ. yarr(ip) ) THEN
             hh   = 0.0
             ynew = yarr(ip)
          ELSE
             hh   = -( xarr(ip+1) - xarr(ip) ) / 
     $            LOG( yarr(ip+1) / yarr(ip))
             ynew = yarr(ip) * EXP(- ( arg - xarr(ip) ) / hh )
          ENDIF
*
*     linear interpolation
*
       ELSEIF ( itype .EQ. 2 ) THEN
         hh = ( yarr(ip+1) - yarr(ip) ) / ( xarr(ip+1) - xarr(ip) )
         ynew = yarr(ip) + hh*( arg - xarr(ip) )
*     
      ENDIF
      RETURN
      END       
