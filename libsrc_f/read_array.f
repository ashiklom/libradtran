      SUBROUTINE read_array( array, error, max_x, nx, ny, filenm )
*
*BOD
*DESCRIPTION
* Read two-dimensional array. Number of elements nx in x-direction
* (number of columns)  must be specified. Number of rows, ny, is
* calculated and returned.
*USAGE
**FORTRAN
*     CALL read_array( array, error, max_x, nx, ny, filenm )
**PERL
*INPUT
**filenm
*     character
*     File name, must contain NO BLANKS.
**max_x
*     integer
*     First dimension of array. To read a one dimensional array, set
*     max_x = nx = 1.
**nx
*     integer
*     Number of columns.
*OUTPUT
**array
*     real array(nx,ny)
*     The array read.
**error
*     integer
*     Error flag. Error checking performed only if first is true.
***0     Everything ok.
***1     Error occurred during reading of data file.
**ny
*     integer
*     Number of rows.
*ERROR_MESSAGES
*ROUTINES_CALLED
*     o3_molina
*FILES
*     consts.inc               Include file
*REFERENCES
*EOD
*
      IMPLICIT LOGICAL (A-Z)  ! To nearly get strong typing
      CHARACTER*(*) filenm 
      INTEGER  error, ios, iunit, ix, iy, max_x, nx, ny
      LOGICAL first
*
      REAL array(max_x, *)
*
* Internal variables
*     
      CHARACTER*100 line
*
      INCLUDE 'consts.inc'
*
      error = 0
*
* Read in all data on first entry only
*
      OPEN(UNIT=iunit,FILE=filenm,STATUS='OLD',
     $     FORM='FORMATTED',ERR=99, IOSTAT=ios)
      GOTO 991
 99   CONTINUE
      WRITE(*,'(A,A)') 'Read_array: error during read of file: ',
     $     filenm
      WRITE(*,'(A, i4)') 'IOSTAT: ', ios
      error = 1
      RETURN
 991  CONTINUE
*     
      iy = 1
      DO WHILE ( .TRUE. ) 
         READ(iunit,'(A)',END=98 ) line
*
* Read all the no-comment lines
*
         IF ( line(:1) .NE. comcha ) THEN
            READ(line,*) (array(ix, iy),ix=1, nx)
            iy = iy + 1
         ENDIF
      ENDDO
 98   CONTINUE
      ny = iy - 1
      CLOSE(iunit)
*
      RETURN
      END
      
