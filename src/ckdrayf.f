c CKDRAYF calculates the Rayleigh optical thickness profile 
c for band ib, subband ig
c ATTENTION: ckdrayf() requires that ckdfuf() has been called before
c            because all the initializations are done there!
c changed error message: increase nvx

      subroutine ckdrayf (ib, ig, u0, nlev, tau)
      
      IMPLICIT NONE

      INCLUDE 'fl_radparams.inc'
      INCLUDE 'fl_special.inc'

c Function parameters
      INTEGER IB, IG
      REAL u0
      INTEGER nlev
      REAL tau(nvx)

c Internal variables 
      INTEGER I

c Common blocks

c Rayleigh optical thickness
      REAL tr,wr,wwr
      common /ray/ tr(nvx), wr(nvx), wwr(nvx,4)

      LOGICAL lchou,lband6a,lpar,lray
      common /chou/ lchou,lband6a,lpar,lray
      


c initialize boundaries
      nv    = nlev-1
      nv1   = nv+1

      ndfs  = nv
      mdfs  = nv1
      ndfs4 = 4*ndfs 
      ndfs2 = 2*ndfs
      mb    = 18
      mbs   = 6
      mbir  = 12 
      nc    = 8 
      
      if(nlev.GT.nv1x) then
        write (*,*) 'error: increase nvx to at least', nlev-1
        write (*,*) '       in fl_radparams.inc and fl_radparams.cinc'
        write (*,*) '       and recompile the package (ckdrayf.f)'
        STOP
      end if

c initialize Rayleigh cross sections
      call rayle ( ib, u0, ig, lray)
      
      DO I = 1, nv
         tau(i) = tr(i)
      ENDDO
      
      END
