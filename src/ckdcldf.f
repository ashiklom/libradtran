c CKDCLDF calculates the water cloud optical properties profile 
c for a given effective radius and liquid water content
c ATTENTION: ckdcldf() requires that ckdfuf() has been called before
c            because all the initializations are done there!
c changed error message: increase nvx

      subroutine ckdcldf (z, reff, lwc, nlev, tau, gg, ssa)
      
      IMPLICIT NONE

      INCLUDE 'fl_radparams.inc'
      INCLUDE 'fl_special.inc'


c Function parameters
      REAL z(nvx), reff(nvx), lwc(nvx)
      INTEGER nlev
      REAL tau(mbx, nvx), ssa(mbx, nvx), gg(4,mbx,nvx)

c Internal variables 
      INTEGER I, J, IB

c Common blocks
      
c cloud physical properties 
      real pre,plwc,pde,piwc
      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)

c cloud optical properties
      real tw,ww,www
      common /wat/ tw(nvx), ww(nvx), www(nvx,4)

c dz, layer thickness; computation is a little bit different from
c the Fu and Liou code: 
c  here,        dz = z(i) - z(i+1)
c  Fu and Liou, dz calculated from pressure difference
      REAL dz
      common /thick/ dz(nvx)


c initialize boundaries
      nv    = nlev-1
      nv1   = nv+1
      
      ndfs  = nv
      mdfs  = nv1
      ndfs4 = 4*ndfs 
      ndfs2 =  2*ndfs
      mb    = 18
      mbs   = 6
      mbir  = 12 
      nc    = 8 
      
      if(nlev.GT.nv1x) then
        write (*,*) 'error: increase nvx to at least', nlev-1
        write (*,*) '       in fl_radparams.inc and fl_radparams.cinc'
        write (*,*) '       and recompile the package (ckdcldf.f)'
        STOP
      end if



c copy effective radius and liquid water content to common block
c need to shift by 1 to get libradtran layering convention
c      DO i=1, nv
c         pre(i)  = reff(i+1)
c         plwc(i) = lwc(i+1)
c      ENDDO
      DO i=1, nv
         pre(i)  = reff(i)
         plwc(i) = lwc(i)
      ENDDO

c initialize layer widths (replacement for thicks())
      DO I = 1, nv
         dz(i) = z(i) - z(i+1)
      ENDDO 

c Loop over bands
      DO ib = 1, mb

c Use Yong Hu's Water cloud optical properties for SW bands

         if (ib .LE. 6 ) then
            call water_hu (ib) 
         else
            call water ( ib )
         endif
         
         DO i = 1, nv1
            tau(ib,i) = tw(i)
            ssa(ib,i) = ww(i)
            DO j=1, 4
               gg(j,ib,i) = www(i,j)
            ENDDO
         ENDDO
         
      ENDDO
      END
