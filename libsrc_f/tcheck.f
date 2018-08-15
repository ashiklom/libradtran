      INTEGER FUNCTION TCHECK(NLYR,NTAU,NSTR,NUMU,NPHI,OPTIMIZE,DELTA)

      INTEGER NLYR, NTAU, NSTR, NUMU, NPHI, OPTIMIZE, DELTA

      INTEGER   MXCLY, MXULV
      INCLUDE  'DISORT.MXD'
      INTEGER   OPTIMIZE_ANY
      INTEGER   OPTIMIZE_MXCLY 
      INTEGER   OPTIMIZE_MXULV

      TCHECK = 0  

      IF( MXCLY.LT.NLYR ) THEN
         write (0,*) '* Error, too many atmospheric layers.' 
         write (0,*) '* Increase MXCLY in libsrc_f/DISORT.MXD'
         write (0,*) '* to at least',NLYR,' and recompile libRadtran!'
         TCHECK=-1
      ENDIF

      IF( MXULV.LT.NTAU ) THEN
         write (0,*) '* Error, too many output altitudes.' 
         write (0,*) '* Increase MXULV in libsrc_f/DISORT.MXD'
         write (0,*) '* to at least',NTAU,' and recompile libRadtran!'
         TCHECK=-1
      ENDIF

      OPTIMIZE_ANY = 0
      OPTIMIZE_MXCLY = 0
      OPTIMIZE_MXULV = 0

c     sometimes layer number are changing by one or two
c     this is taken care of by delta

      IF( NLYR.LT.MXCLY-DELTA .AND. NTAU.LT.2*MXCLY-DELTA 
     &                                         .AND. MXCLY.GT.1) THEN
         OPTIMIZE_MXCLY=MAX0(NLYR,NTAU/2,1)+DELTA
         OPTIMIZE_ANY = OPTIMIZE_ANY + 1
      ENDIF

c      IF( NTAU.LT.MXULV ) THEN
c         OPTIMIZE_MXULV=1
c         OPTIMIZE_ANY = 1
c      ENDIF

      IF ( OPTIMIZE_ANY .GT. 0 .AND. OPTIMIZE .EQ. 1 ) THEN
        write (0,*) ''
        write (0,*) ' ***  Optimisation for FORTRAN ',
     &              'array-dimensions possible:' 
        IF (OPTIMIZE_MXCLY .NE. 0) THEN 
         write (0,*)'      decrease MXCLY (atmospheric layers) to ',
     &       OPTIMIZE_MXCLY, ' (',MXCLY,')'
        ENDIF
c        IF (OPTIMIZE_MXULV .NE. 0) THEN 
c         write (0,*) '     decrease MXULV (output altitudes) to ',NTAU,
c     &        '(',MXULV,')'
c        ENDIF
        write (0,*)  '      in libsrc_f/DISORT.MXD, '

        IF (OPTIMIZE_MXCLY .NE. 0) THEN 
         write (0,*)'      and'
         write (0,*)'      decrease nvx (atmospheric layers) to ',
     &       OPTIMIZE_MXCLY, ' (',MXCLY,')'
         write (0,*)  '      in src/fl_radparams.inc, and '
         write (0,*)  '      in src/fl_radparams.cinc'
        ENDIF
        write (0,*)  '      and recompile libRadtran!'
        write (0,*) ''
      ENDIF

      RETURN
      END
