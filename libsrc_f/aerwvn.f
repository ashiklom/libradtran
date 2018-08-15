      SUBROUTINE aerwvn( aerwl, aerabs, aersym, aerext, aer_pro, 
     $     aer_dtau, aer_gg, aer_ssa, nlyr, lambda, zd, nlambda,
     $     maxawvn)
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
* Calculates absorption and scattering optical depth and asymmetry
* factors due aerosols.
*     
* Input/output variables described in phodis.f
*     
      IMPLICIT LOGICAL (A-Z)    ! To nearly get strong typing
      INTEGER iv, iva, lc, maxawvn, nlambda, nlyr
      REAL aerwl(*), aerabs(maxawvn,0:*),
     $     aerext(maxawvn,0:*), aer_pro(0:*), aersym(maxawvn,0:*),
     $     aer_dtau(nlyr,*), aer_gg(nlyr,*), aer_ssa(nlyr,*),
     $     lambda(*), zd(0:*) 
      REAL aabs0, aabs1, aersca0, aersca1,
     $     aersym0, aersym1, asca0, asca1, asym0, asym1, 
     $     babsaer, basya, bscaaer, daerwl, deltazkm, dwvn, wvna
*
      DO iv = 1, nlambda
         DO lc = 1, nlyr
            deltazkm = ( zd(lc-1) - zd(lc) ) ! Must be in kilometers for aerosols
            wvna = lambda(iv)
            DO iva = 1, maxawvn-1 !Find right wavelength interval
               IF ((wvna .GE. aerwl(iva)) .AND.
     $              (wvna .LT. aerwl(iva+1)))        GOTO 95
            ENDDO
 95         CONTINUE
            IF ( iva+1 .GT. maxawvn ) iva = maxawvn-1
            IF (wvna .LT. aerwl(1))   iva = 1
            dwvn = aerwl(iva+1)-wvna
            daerwl = aerwl(iva+1) - aerwl(iva)
*     
*     Scattering coefficient for aerosols
*     
            aersca0 = aerext(iva,lc-1)-aerabs(iva,lc-1)
            aersca1 = aerext(iva+1,lc-1)-aerabs(iva+1,lc-1)
            asca0=(aersca0-aersca1)*dwvn/daerwl+aersca1
            aersca0 = aerext(iva,lc)-aerabs(iva,lc)
            aersca1 = aerext(iva+1,lc)-aerabs(iva+1,lc)
            asca1=(aersca0-aersca1)*dwvn/daerwl+aersca1
            bscaaer = 0.5*deltazkm*(asca0*aer_pro(lc-1)
     $           + asca1*aer_pro(lc))
            IF ( bscaaer .LT. 0.0 ) bscaaer = 0.0
*     
*     Absorption coefficient for aerosols
*     
            aabs0=(aerabs(iva,lc-1)-aerabs(iva+1,lc-1))
     $           *dwvn/daerwl+aerabs(iva+1,lc-1)
            aabs1=(aerabs(iva,lc)-aerabs(iva+1,lc))*
     $           dwvn/daerwl+aerabs(iva+1,lc)
            babsaer = 0.5*deltazkm*(aabs0*aer_pro(lc-1)
     $           + aabs1*aer_pro(lc))
            IF ( babsaer .LT. 0.0 ) babsaer = 0.0
*     
*     Now lets do the asymmetry factor for the aerosols
*     
            aersym0 = aersym(iva,lc-1)
            aersym1 = aersym(iva+1,lc-1)
            asym0=((aersym0-aersym1)*dwvn/daerwl+aersym1)
            aersym0 = aersym(iva,lc)
            aersym1 = aersym(iva+1,lc)
            asym1=((aersym0-aersym1)*dwvn/daerwl+aersym1)
            basya =  0.5 * (asym0 + asym1)
*
            aer_dtau(lc,iv) = bscaaer + babsaer
            aer_gg(lc,iv) = basya
            IF ( (bscaaer + babsaer) .LE. 0.0) THEN
               aer_ssa(lc,iv) = 0.0
            ELSE
               aer_ssa(lc,iv) = bscaaer / (bscaaer + babsaer)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
