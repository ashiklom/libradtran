c=======================================================================
c SBDART gaseous absorption module
c
c BM extracted code from SBDART version 2.2, March 31, 2003
c Changes:
c  - moved all constants to sbtaugas.inc; added "INCLUDE 'sbtaugas.inc'
c    to every subroutine
c=======================================================================
      subroutine taucor(gwk, tau, amu, utau, cf)
      IMPLICIT NONE
c
c purpose: find a correction factor, cf, which corrects the
c  k-distribution optical depths so that 
c
c   sum(gwk(1:3)*exp(-cf*tau(1:3)/amu))=exp(-utau)
c
c input:
c 
c   gwk         k-distribution weights
c   tau         k-distribution optical depths
c   amu         cosine of solar zenith
c   utau        optical depth along slant path due to line absorption
c
c output:
c   cf          correction factor,  tau_corrected=cf*tau_uncorrected
c
c  
      INCLUDE 'sbtaugas.inc'
      integer i, j1
      double precision tau(3), gwk(3), amu, utau, cf, ff, f, fp, dcf
      double precision d1, d2
 
      cf = 1.
 
      if (utau .gt. 12.0) return 
 
      do i = 1, 20
         d1 = 0
         do j1 = 1, 3
            d1 = d1 + gwk(j1)*exp((-cf*tau(j1)/amu))
         end do
         ff = d1
         f = log(ff) + utau
         if (abs(f) .lt. 0.000001) return 
         d2 = 0
         do j1 = 1, 3
            d2 = d2 + gwk(j1)*exp((-cf*tau(j1)/amu))*tau(j1)
         end do
         fp = -d2/(ff*amu)                       ! df/dcf
         dcf = -f/fp
         cf = cf + dcf
      end do
 
      write (*, '(8a12)') '------------', 'gwk  ', '------------', 
     1   '------------', 'tau  ', '------------', 'amu', 'utau'
      write (*, '(8(1p,e12.4,0p))') gwk, tau, amu, utau
      stop 'TAUCOR: iteration did not converge'
 
      end 
