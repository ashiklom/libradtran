*
* @20i@
* Function: wcloud
* Description:
* @code{wcloud} calculates the asymmetry factor, the single scattering
* albedo and the extinction optical depth of water clouds for a single
* wavelength. The parameterization due to Hu and Stamnes (1993) is used. 
* Parameters:
*    REAL lambda: The wavelength in nanometers. (input)
*    LOGICAL newsiz: Set @code{newsiz=.TRUE.} whenever new water cloud liquid 
*                    water content and/or effective radius are given. If 
*                    @code{wccon} and @code{wceffr} are the same as in the 
*                    previous call set @code{newsiz = .FALSE.} to save
*                    computer time. Must be equal to @code{.TRUE.} on the 
*                    first call. (input)
*    INTEGER nlyr:   Number of atmospheric layers. (input)
*    CHARACTER*(*) path: The filepath to water cloud parameterization files. 
*                    (input)
*    REAL(*) wccon:  The liquid water content of each layer in grams per cubic 
*                    meter. nlyr+1 elements.(input)
*    REAL(*) wceffr: The water droplet effective radius in microns. nlyr+1 elements.
*                    (input)
*    REAL(*) wc_dtau: The water cloud optical depth of each layer. nlyr elements. 
*                    (output)
*    REAL(*) wc_gg:  The water cloud asymmetry factor of each layer. nlyr elements. 
*                    (output)
*    REAL(*) wc_ssa: The water cloud single scattering of each layer. nlyr elements.
*                    (output)
*    REAL(*) zd:     The altitude of each level in km. @code{zd(nlyr)} is the 
*                    bottom of the atmosphere. nlyr+1 elements. (input)
*    INTEGER wclyr:  If wclyr is .EQ. 0, the cloud properties are defined 
*                    per level, otherwise per layer. (input)                    
*
* Return value:
* Example:
* See @file{uvspec.c}.
* Files: 
* @file{DISORT.MXD} with parameter @code{mxcly} must be present in the 
* same directory as @file{wcloud.f}.
* Known bugs:
* Author:  Arve Kylling    
* Important changes: RB 08.12.2010: replaced mxcly by 2^15-1 to allow
*                                   more atmospheric layers
*     
* @i20@
*    
      SUBROUTINE wcloud ( wavlen, newsiz, nlyr, path, nstring,
     $     wccon, wceffr, wc_dtau, wc_gg, wc_ssa, zd, wclyr)
*-------------------------------------------------------------------
* Copyright (C) 1997 Arve Kylling
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
      IMPLICIT LOGICAL (A-Z)    ! To nearly get strong typing

      INTEGER MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MXSZA
      INCLUDE 'DISORT.MXD'
*
**mxwvwc     Maximum number wavelength points for water cloud.
**mxsiz      Maximum number sizes for water clouds.
**mxcof      Number of coefficients for water cloud parameterization.
*
*
*IO-parameters
*
      INTEGER iunit
      PARAMETER ( iunit = 1 )
      INTEGER mxcof, mxsiz, mxwvwc
      PARAMETER (  mxwvwc = 74, mxsiz = 3, mxcof = 3 ) 
*
      REAL r1mach
      CHARACTER*(*) path
      INTEGER icof, isize, ivwc, iwcsiz(0:32767), lc, nlyr, nstring
      INTEGER wclyr
      LOGICAL first, first_line, newsiz
      REAL cl_asy(mxsiz,mxwvwc,mxcof), cl_ext(mxsiz,mxwvwc,mxcof),
     $     cl_ssa(mxsiz,mxwvwc,mxcof), wcwvn(mxwvwc), wccon(0:*), 
     $     wceffr(0:*), wc_dtau(*), wc_gg(*), wc_ssa(*), 
     $     wvnwcz(mxsiz+1), zd(0:*)
      CHARACTER*1  char
      CHARACTER*9 filenm    ! No. of characters must equal
                             ! length of filename.
      CHARACTER*100 fullname, line
      REAL basyd, bextd, bssad,
     $     deltaz, ex_v1, ex_v2, ex_v1_l1, ex_v1_l2,
     $     ex_v2_l1, ex_v2_l2, gg_v1, gg_v2, gg_v1_l1, gg_v1_l2,
     $     gg_v2_l1, gg_v2_l2, ss_v1, ss_v2, ss_v1_l1, ss_v1_l2,
     $     ss_v2_l1, ss_v2_l2, sprec, ssalbd, lambda, wavlen
*
      REAL wc_avg, reff_avg
*
      SAVE sprec, iwcsiz, cl_asy, cl_ext, cl_ssa, wcwvn, wvnwcz
      DATA first / .TRUE. /
*
      IF ( first ) THEN
         first = .FALSE.
         sprec = r1mach(1)
         IF( 32767.LT.nlyr ) THEN
            WRITE ( 0,* )  ' ****  wcloud can only handle',
     &           ' 32767 levels, change inside wcloud.f! ',
     &           ' **** '
            STOP
         ENDIF
*
*     Read in water cloud optical parameters.
*     
         filenm = 'wc/wc.asy'
         fullname = path(1:nstring)//filenm
         OPEN(UNIT=iunit,FILE=fullname,FORM='FORMATTED',STATUS='OLD',
     $        ERR=99  )
         DO 100 isize = 1, mxsiz
            READ(iunit,*) 
            READ(iunit,*) 
            READ(iunit,*) wvnwcz(isize), wvnwcz(isize)
            READ(iunit,*) 
            DO 101 ivwc = 1, mxwvwc
               READ(iunit,*) wcwvn(ivwc),
     $              (cl_asy(isize, ivwc,icof),icof=1,mxcof)
 101        CONTINUE
 100     CONTINUE
         CLOSE(iunit)  
*     
         filenm = 'wc/wc.ext'
         fullname = path(1:nstring)//filenm
         OPEN(UNIT=iunit,FILE=fullname,FORM='FORMATTED',STATUS='OLD',
     $        ERR=99  )
         DO 110 isize = 1, mxsiz
            READ(iunit,*) 
            READ(iunit,*) 
            READ(iunit,*) wvnwcz(isize), wvnwcz(isize)
            READ(iunit,*) 
            DO 111 ivwc = 1, mxwvwc
               READ(iunit,*) wcwvn(ivwc),
     $              (cl_ext(isize, ivwc,icof),icof=1,mxcof)
 111        CONTINUE
 110     CONTINUE
         CLOSE(iunit)  
*     
         filenm = 'wc/wc.ssa'
         fullname = path(1:nstring)//filenm
         OPEN(UNIT=iunit,FILE=fullname,FORM='FORMATTED',STATUS='OLD',
     $        ERR=99  )
         DO 120 isize = 1, mxsiz
            READ(iunit,*) 
            READ(iunit,*) 
            READ(iunit,*) wvnwcz(isize), wvnwcz(isize)
            READ(iunit,*) 
            DO 121 ivwc = 1, mxwvwc
               READ(iunit,*) wcwvn(ivwc),
     $              (cl_ssa(isize, ivwc,icof),icof=1,mxcof)
 121        CONTINUE
 120     CONTINUE
         CLOSE(iunit)  
         DO ivwc = 1, mxwvwc         ! Change from microns to nm
            wcwvn(ivwc) = wcwvn(ivwc) * 1.0E+03 
         ENDDO
      ENDIF
*     
*     Find the size-interval to be used for water clouds.
*     
      IF ( newsiz ) THEN
         DO 130 lc = 1, nlyr

            wc_avg = 0.0
            reff_avg = 0.0
            IF ( wclyr .EQ. 0 .AND. 
     $           wccon(lc-1).GT.0.0 .AND. wccon(lc).GT.0.0 ) THEN
               reff_avg = 0.5*(wceffr(lc-1)+wceffr(lc))
               wc_avg   = 0.5*(wccon(lc-1)+wccon(lc))
            ELSE
               reff_avg = wceffr(lc)
               wc_avg   = wccon(lc)
            ENDIF
            IF ( wc_avg .GT. 0.0) THEN
               IF ( reff_avg .LT. 2.5 .OR. reff_avg .GE. 60.0 ) THEN
                  GOTO 88
               ENDIF
            ENDIF
            IF ( wc_avg .GT. 0.0 .AND. reff_avg.GT.0.0 ) THEN
               isize = 1
               IF ( reff_avg .GE. 2.5  .AND.
     $              reff_avg .LT. wvnwcz(isize) ) THEN
                  iwcsiz(lc) = isize 
               ELSE 
                  DO 140 isize = 2, mxsiz
                     IF ( reff_avg .GE. wvnwcz(isize-1)  .AND.
     $                    reff_avg .LT. wvnwcz(isize) ) THEN
                        iwcsiz(lc) = isize
                     ENDIF
 140              CONTINUE
               ENDIF
            ENDIF
 130     CONTINUE
      ENDIF
*     
*     Mie scattering from water drops (water clouds)
*     

*     If wavelength outside the range, set it to the boundary
*     to make sure that an at least somewhat reasonable value is
*     returned, BM, April 6, 2001
*
      lambda = wavlen
      IF ( lambda.GE.wcwvn(1) ) THEN
         lambda = wcwvn(1)-1
      ENDIF
      
      IF ( lambda.LE.wcwvn(mxwvwc) ) THEN
         lambda = wcwvn(mxwvwc)+0.01
      ENDIF

      DO  200 lc = 1, nlyr
         deltaz = ( zd(lc-1) - zd(lc) ) ! Must be in kilometers for
         bextd  = 0.0
         ssalbd = 0.0
         basyd  = 0.0
*                                     ! water clouds
         IF ( wclyr .EQ. 0 ) THEN
*     optical properties defined per level
*     need to interpolate between adjacent levels
            IF (wccon(lc-1).GT.0.0 .AND. wccon(lc).GT.0.0) THEN
               DO 210 ivwc = 1, mxwvwc-1
                  IF ( lambda.LE.wcwvn(ivwc) .AND.
     $                 lambda.GT.wcwvn(ivwc+1)) THEN
*     
*     Linearly interpolate coefficients to lambda
*
                     reff_avg = 0.5*(wceffr(lc-1)+wceffr(lc))
                     wc_avg   = 0.5*(wccon(lc-1)+wccon(lc))
*
*     Extinction coefficient
*     
                     ex_v1    = ( cl_ext(iwcsiz(lc),ivwc,1) * (
     $                    reff_avg ** cl_ext(iwcsiz(lc),ivwc,2) )
     $                    + cl_ext(iwcsiz(lc),ivwc,3) ) * deltaz
                     ex_v2    = ( cl_ext(iwcsiz(lc),ivwc+1,1) * (
     $                    reff_avg ** cl_ext(iwcsiz(lc),ivwc+1,2) )
     $                    + cl_ext(iwcsiz(lc),ivwc+1,3) ) * deltaz
                     CALL linpol(wcwvn(ivwc), ex_v1, wcwvn(ivwc+1),
     $                    ex_v2, lambda, bextd )
                     bextd = bextd*wc_avg
*     
*     Asymmetry factor
*     
                     gg_v1 = cl_asy(iwcsiz(lc-1),ivwc,1) * (
     $                    reff_avg ** cl_asy(iwcsiz(lc-1),ivwc,2) )
     $                    + cl_asy(iwcsiz(lc-1),ivwc,3)
                     gg_v2 = cl_asy(iwcsiz(lc-1),ivwc+1,1) * (
     $                    reff_avg ** cl_asy(iwcsiz(lc-1),ivwc+1,2) )
     $                    + cl_asy(iwcsiz(lc-1),ivwc+1,3)
                     CALL linpol(wcwvn(ivwc), gg_v1, wcwvn(ivwc+1),
     $                    gg_v2, lambda, basyd )
*     
*     Single scattering albedo
*     
                     ss_v1 = 1. - (cl_ssa(iwcsiz(lc),ivwc,1) * (
     $                    reff_avg ** cl_ssa(iwcsiz(lc),ivwc,2) )
     $                    + cl_ssa(iwcsiz(lc),ivwc,3))
                     ss_v2 = 1. - (cl_ssa(iwcsiz(lc),ivwc+1,1) * (
     $                    reff_avg ** cl_ssa(iwcsiz(lc),ivwc+1,2) )
     $                    + cl_ssa(iwcsiz(lc),ivwc+1,3))
                     CALL linpol(wcwvn(ivwc), ss_v1, wcwvn(ivwc+1),
     $                    ss_v2, lambda, bssad )
                  ENDIF
 210           CONTINUE
               ssalbd = bssad
               IF ( bextd  .LT. sprec ) bextd  = sprec
               IF ( ssalbd .GT.  1.0  ) ssalbd =  1.0
               IF ( basyd  .LT.  0.0  ) basyd  =  0.0
               wc_gg(lc)      = basyd
               wc_dtau(lc)    = bextd
               wc_ssa(lc)     = ssalbd
            ENDIF
         ELSE
*     optical properties defined per layer
*     need to interpolate between adjacent levels
            IF (wccon(lc).GT.0.0) THEN
               DO 211 ivwc = 1, mxwvwc-1
                  IF ( lambda.LE.wcwvn(ivwc) .AND.
     $                 lambda.GT.wcwvn(ivwc+1)) THEN
*     
*     Linearly interpolate coefficients to lambda
*     
*     Extinction coefficient
*     
                     ex_v1 = ( cl_ext(iwcsiz(lc),ivwc,1) * (
     $                    wceffr(lc) ** cl_ext(iwcsiz(lc),ivwc,2) )
     $                    + cl_ext(iwcsiz(lc),ivwc,3) ) * deltaz
                     ex_v2 = ( cl_ext(iwcsiz(lc),ivwc+1,1) * (
     $                    wceffr(lc) ** cl_ext(iwcsiz(lc),ivwc+1,2) )
     $                    + cl_ext(iwcsiz(lc),ivwc+1,3) ) * deltaz
                     ex_v1 = ex_v1*wccon(lc)
                     ex_v2 = ex_v2*wccon(lc)
                     CALL linpol(wcwvn(ivwc), ex_v1, wcwvn(ivwc+1),
     $                    ex_v2, lambda, bextd )
*     
*     Asymmetry factor
*     
                     gg_v1 = cl_asy(iwcsiz(lc),ivwc,1) * (
     $                    wceffr(lc) ** cl_asy(iwcsiz(lc),ivwc,2) )
     $                    + cl_asy(iwcsiz(lc),ivwc,3)
                     gg_v2 = cl_asy(iwcsiz(lc),ivwc+1,1) * (
     $                    wceffr(lc) ** cl_asy(iwcsiz(lc),ivwc+1,2) )
     $                    + cl_asy(iwcsiz(lc),ivwc+1,3)
                     CALL linpol(wcwvn(ivwc), gg_v1, wcwvn(ivwc+1),
     $                    gg_v2, lambda, basyd )
*     
*     Single scattering albedo
*     
                     ss_v1 = 1. - (cl_ssa(iwcsiz(lc),ivwc,1) * (
     $                    wceffr(lc) ** cl_ssa(iwcsiz(lc),ivwc,2) )
     $                    + cl_ssa(iwcsiz(lc),ivwc,3))
                     ss_v2 = 1. - (cl_ssa(iwcsiz(lc),ivwc+1,1) * (
     $                    wceffr(lc) ** cl_ssa(iwcsiz(lc),ivwc+1,2) )
     $                    + cl_ssa(iwcsiz(lc),ivwc+1,3))
                     CALL linpol(wcwvn(ivwc), ss_v1, wcwvn(ivwc+1),
     $                    ss_v2, lambda, bssad )
                  ENDIF
 211           CONTINUE
               ssalbd = bssad
               IF ( bextd  .LT. sprec ) bextd  = sprec
               IF ( ssalbd .GT.  1.0  ) ssalbd =  1.0
               IF ( basyd  .LT.  0.0  ) basyd  =  0.0
               wc_gg(lc)      = basyd
               wc_dtau(lc)    = bextd
               wc_ssa(lc)     = ssalbd
            ENDIF
         ENDIF
 200  CONTINUE
*     
      RETURN
*     
 88   CONTINUE
      WRITE(*,*)'Error, in water cloud parameterization Hu and Stamnes'
      WRITE(*,*)'       droplet size', wceffr(lc-1),reff_avg, wceffr(lc)
      WRITE(*,*)'       < than 2.5 micron or larger than 60 micron for '
      WRITE(*,*)'       layer ', lc, '.'
*
      STOP
 99   CONTINUE
      WRITE(*,*)'Error during read of file:',fullname
*
      STOP
      END
