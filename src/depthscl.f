c=======================================================================
c SBDART gaseous absorption module
c
c BM extracted code from SBDART version 1.50, September 12, 2002
c The code in depthscl.f required heavy changes to work with libRadtran;
c in particular, only the gaseous absorption things were kept while all 
c aerosol, cloud, and Rayleigh scattering stuff was removed.
c=======================================================================
      subroutine depthscl(kdist, kd, nk, nz, wl, dtausc,
     &     gwk, dtauk, dtaugc, 
     &     wt, dtau)
      IMPLICIT NONE

c input:
c   kdist   k-distribution mode index
c           -1  read k-distribution dtauk, gwk from LOUTF
c            0  tau set to -log(transmission_vertical_path)
c               slant optical depth used when tau_scat < 1 and wl < 4
c            1  no renormalization of lowtran k-distribution tau
c            2  corrected to solar beam slant path transmission
c            3  corrected to solar beam slant path transmission when
c               tau_scat < 1 and wl < 4
c
c  kd      k-distribution index
c  nk      number of k-distribution passes
c  nz      number of vertical layers
c  wl      wavelength
c  dtausc  scattering optical depth increment
c  gwk     k-distribution weights
c  dtauk   k-distribution optical depths for absorption by molecular lines
c  dtaugc  gas continuum optical depth
c
c output:
c  dtau   atmospheric optical depth increment
 
      INCLUDE 'sbtaugas.inc'
      integer kd, nk, nz, kdist, i, iw
 
      integer mk                        ! number of lowtran k-dist terms
      parameter (mk = 3)
 
      double precision dtausc(mxly),
     &   dtau(mxly), wt, dtaug(mxly), afac, tglv, tgls, tsc, 
     &   ramp, wl, gwk(mxkd), dtauk(mxly,mxkd), dtaugc(mxly)
 
      doubleprecision d1, d2, d3, d4, d5
      integer j1
                                ! get gas optical depth
 
      wt = gwk(kd)
      if (kdist .eq. (-1)) then                  ! tau from LOUTF
         do j1 = 1, nz
            dtaug(j1) = dtauk(j1,kd)
         end do
      else if (kdist.eq.0 .or. nk.eq.1) then    ! tau=-log(transmission)
         wt = 1.
         tsc = 0.
         tglv = 0.
         tgls = 0.
         do i = 1, nz
            tglv = tglv + dtauk(i,1)             ! vertical path
            tgls = tgls + dtauk(i,1+3)           ! amu0 * slant path
            tsc = tsc + dtausc(i)
            afac = 1.
            if (tglv .gt. .001) afac = tgls/tglv
            call depthscl_rolloff (wl, tsc, ramp)
            afac = afac*ramp + 1. - ramp
            dtaug(i) = dtaugc(i) + dtauk(i,1)*afac
         end do
      else if (kdist .eq. 1) then                ! uncorrected
         do j1 = 1, nz
            dtaug(j1) = dtaugc(j1) + dtauk(j1,kd)
         end do
      else if (kdist .eq. 2) then ! corrected to slant path transmission
         do j1 = 1, nz
            dtaug(j1) = dtaugc(j1) + dtauk(j1,kd+3)
         end do
      else                                       ! ramp to uncorrected
         tsc = 0.
         do i = 1, nz
            tsc = tsc + dtausc(i)
            call depthscl_rolloff (wl, tsc, ramp)
            dtaug(i)=dtaugc(i)+dtauk(i,kd)*(1.-ramp)+dtauk(i,kd+3)*ramp
         end do
      endif
                                ! zero out optical depth in sub-surface
      do i = 1, nz
         dtau(i) = dtaug(i) + dtausc(i)
      end do
 
      return 
      end 
c-----------------------------------------------------------------------
      subroutine depthscl_rolloff(wl, tsc, ramp)
c purpose: compute roll-off factor for transition from full slant-path
c          correction to uncorrected vertical optical depths. 
c
c input:
c  wl      wavelength
c  tsc     scattering optical depth
c output:
c  ramp    roll-off factor (drops to zero when either the scattering
c          optical depth is large, or the wavelength exceeds wlhi)

      INCLUDE 'sbtaugas.inc'

      double precision wl, tsc, ramp
      double precision wllo, wlhi      ! wavelength roll-off upper limit
      parameter (wllo = 3.9, wlhi = 4.1)
      doubleprecision d1, d2, d3, d4, d5
      integer j1
 
      ramp = (wlhi - wl)/(wlhi - wllo)
      ramp = max(min(one,ramp),zero)
      ramp = ramp*exp(1. - max(tsc,one))
      return 
 
      end 
