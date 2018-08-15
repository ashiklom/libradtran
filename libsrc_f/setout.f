      SUBROUTINE setout( dtauc, nlyr, ntau, utau, z, zout ) 
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
* Linearly interpolate to get approximate tau corresponding to 
* altitude zout
*
* Input/output variables described in phodis.f
*
      INCLUDE 'DISORT.MXD'
      INTEGER itau
      INTEGER nlyr, ntau
      REAL dtauc(*), hh, utau(*),
     $     z(*), zout(*), utauout
      REAL tauint(mxcly+1)

*     
      IF( mxcly.LT.nlyr ) THEN
         WRITE ( 0,* )  ' **** Error during setout (setout.f)'
         WRITE ( 0,* )  ' ****  Symbolic dimension mxcly',
     &        ' should be increased to at least ', NLYR,
     &        ' ****  in libsrc_f/DISORT.MXD'
         STOP
      ENDIF

      tauint(1) = 0.
      DO  lc = 1, nlyr
         tauint(lc+1) = tauint(lc) + dtauc(lc)
      ENDDO
*
      itype = 2
      DO itau=1, ntau
         CALL inter( mxcly+1, nlyr+1, itype, zout(itau),
     $        z, tauint, utau(itau), hh)
*     
      ENDDO

*     
      RETURN
      END
