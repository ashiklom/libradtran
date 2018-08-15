c CKDICEF calculates the ice cloud optical properties profile 
c for a given effective size and liquid water content
c ATTENTION: ckdicef() requires that ckdfuf() has been called before
c            because all the initializations are done there!
c changed error message: increase nvx

      subroutine ckdicef (z, reff, lwc, nlev, tau, gg, ssa, unscaled)

      IMPLICIT NONE

      INCLUDE 'fl_radparams.inc'
      INCLUDE 'fl_special.inc'


c Function parameters
      REAL z(nvx), reff(nvx), lwc(nvx)
      INTEGER nlev
      REAL tau(mbx, nvx), ssa(mbx, nvx), gg(4,mbx,nvx)
      INTEGER unscaled

c Internal variables 
      INTEGER I, J, IB

c Common blocks

c cloud physical properties 
      real pre,plwc,pde,piwc
      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)

c cloud optical properties
      real ti,wi,wwi
      common /ic/ ti(nvx), wi(nvx), wwi(nvx,4)

c dz, layer thickness; computation is a little bit different from
c the Fu and Liou code: 
c  here,        dz = z(i) - z(i+1)
c  Fu and Liou, dz calculated from pressure difference
      REAL dz
      common /thick/ dz(nvx)


c fl93i, FALSE = new 1998 ice cld,  TRUE = old 1993 ice cld
      logical fl93i
      data fl93i /.FALSE./ 
      

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
        write (*,*) '       and recompile the package (ckdicef.f)'
        STOP
      end if

c copy effective radius and liquid water content to common block
c need to shift by 1 to get libradtran layering convention
c      DO i=1, nv
c         pde(i)  = reff(i+1)
c         piwc(i) = lwc(i+1)
c      ENDDO
      DO i=1, nv
         pde(i)  = reff(i)
         piwc(i) = lwc(i)
      ENDDO

c initialize layer widths (replacement for thicks())
      DO I = 1, nv
         dz(i) = z(i) - z(i+1)
      ENDDO 



c Loop over bands

      DO ib = 1, mb
         
         if ( fl93i ) then
            call ice ( ib )
         else
            call icenew ( ib, unscaled )
         endif

         DO i = 1, nv1
            tau(ib,i) = ti(i)
            ssa(ib,i) = wi(i)
            DO j=1, 4
               gg(j,ib,i) = wwi(i,j)
            ENDDO
         ENDDO
         
         
      ENDDO
      END


