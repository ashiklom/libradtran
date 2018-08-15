      SUBROUTINE aerprof( nlyr, z, aerprf, seasn, vulcan, visib )
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
* Gets the aerosol profile depending on season, aerosol conditions
* and tropospheric visibility.
*
* Input/output variables described in phodis.f and uvspec.f
*
      IMPLICIT LOGICAL (A-Z)    ! To nearly get strong typing
      INTEGER ilyr, imlc, itype, jlyr, lc,
     $     m10lyr, m2lyr, m30lyr, mxmlyr
      PARAMETER ( mxmlyr = 34-1 )
*
      INTEGER nlyr, seasn, vulcan
      REAL z(0:*), aerprf(0:*), visib
      REAL aerpfm(mxmlyr), vs(5), zaerdm(mxmlyr), aertmp(mxmlyr)
      REAL const, h
*
      INCLUDE "prfinc.f"
*
      DATA m10lyr / 11 /, m2lyr / 3 /, m30lyr / 27 / 
      DATA vs / 50., 23., 10., 5., 2. /
*
*     First test some input variables
*
      IF ( seasn .NE. 1 .AND. seasn .NE. 2 ) CALL errmsg(
     $     'seasn not set correctly', .true.)
      IF ( vulcan.LT.1 .OR. vulcan.GT. 4 ) CALL errmsg(
     $     'vulcan wrong',.TRUE.)
*
*   Load aerosol profile according to season and aerosol type
*
*  0-2km                           ! Independent of season
*
      DO ilyr = 2, 5
         IF ( visib .GE. vs(ilyr) ) GOTO 10 
      ENDDO
 10   CONTINUE
      jlyr = ilyr
      IF ( jlyr .GT. 5 ) jlyr = 5
      const = 1. / ( 1./vs(jlyr) - 1./vs(jlyr-1) )
      DO imlc = 1, m2lyr
         aerpfm(imlc) = const * (
     $        (hzviz(imlc,jlyr)-hzviz(imlc,jlyr-1))/visib
     $        + hzviz(imlc,jlyr-1)/vs(jlyr)
     $        - hzviz(imlc,jlyr)/vs(jlyr-1)  )
      ENDDO
*
*  SPRING -- SUMMER
*
      IF ( seasn .EQ. 1 ) THEN              !  Spring-summer
*
*  >2-10km
*
         IF ( visib .LE. 23 ) THEN
            DO imlc = m2lyr+1, m10lyr
               IF ( zaerd(imlc) .LE. 4. ) THEN
                  aerpfm(imlc) = spsu23(imlc)
               ELSE
                  aerpfm(imlc) = spsu50(imlc)
               ENDIF
            ENDDO
         ELSE
            const = 1. / ( 1./23. - 1./50. )
            DO imlc = m2lyr+1, m10lyr
               IF ( zaerd(imlc) .LE. 4. ) THEN
                  aerpfm(imlc) = const * (
     $                 ( spsu23(imlc)-spsu50(imlc))/
     $                 visib+spsu50(imlc)/23.-spsu23(imlc)/50.)
               ELSE
                  aerpfm(imlc) = spsu50(imlc)
               ENDIF
            ENDDO
         ENDIF
*
*  >10-100km
*
         IF (vulcan .EQ. 1 ) THEN           !  Background 
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = bastss(imlc) 
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = upnatm(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 2 ) THEN     !  Moderate vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = vumoss(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 3 ) THEN     ! High vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = hivuss(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 4 ) THEN     ! Extreme vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = exvuss(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ENDIF
*
*   FALL -- WINTER
*

      ELSE IF ( seasn .EQ. 2 ) THEN !  Fall-winter
*
*  >2-10km
*
         IF ( visib .LE. 23 ) THEN
            DO imlc = m2lyr+1, m10lyr
               IF ( zaerd(imlc) .LE. 4. ) THEN
                  aerpfm(imlc) = fawi23(imlc)
               ELSE
                  aerpfm(imlc) = fawi50(imlc)
               ENDIF
            ENDDO
         ELSE
            const = 1. / ( 1./23. - 1./50. )
            DO imlc = m2lyr+1, m10lyr
               IF ( zaerd(imlc) .LE. 4. ) THEN
                  aerpfm(imlc) = const * (
     $                 ( fawi23(imlc)-fawi50(imlc))/
     $                 visib+fawi50(imlc)/23.-fawi23(imlc)/50.)
               ELSE
                  aerpfm(imlc) = fawi50(imlc)
               ENDIF
            ENDDO
         ENDIF
*
*  >10-100km
*
         IF (vulcan .EQ. 1 ) THEN           !  Background 
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = bastfw(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = upnatm(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 2 ) THEN     !  Moderate vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = vumofw(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 3 ) THEN     ! High vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = hivufw(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ELSE IF ( vulcan .EQ. 4 ) THEN     ! Extreme vulcanic
            DO imlc = m10lyr+1, m30lyr
               aerpfm(imlc) = exvufw(imlc)
            ENDDO
            DO imlc = m30lyr+1, mxmlyr
               aerpfm(imlc) = vutono(imlc)
            ENDDO
         ENDIF
      ENDIF
*
*     Invert profile so it agrees with phodis layering
*
      DO imlc = 1, mxmlyr
         aertmp(imlc) = aerpfm(imlc)
      ENDDO
      DO imlc = 1, mxmlyr
         aerpfm(imlc) = aertmp(mxmlyr-imlc+1)
         zaerdm(imlc) = zaerd(mxmlyr-imlc+1)
      ENDDO
*
*     Interpolate model atmosphere to user grid defined by z
*
      DO lc = 0, nlyr
*
*     Linear interpolation for aerosol.
*
         itype = 2
         CALL inter( maxcly, mxmlyr, itype, 
     $        z(lc), zaerdm, aerpfm, aerprf(lc), h)
      ENDDO
*
      RETURN
      END
