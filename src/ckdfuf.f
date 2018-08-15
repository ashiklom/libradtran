c CKDFUF does the initialization (set the variables, call ckd1_init(), etc.)
c and calculates the molecular absorption cross sections ck(ib,i,ig)
c changed error message: increase nvx 

      subroutine ckdfuf (co2, ch4, n2o, f11, f12, f22,
     $     rhoair, rhoh2o, rhoo3, 
     $     pres, temp, davg, z, nlev, h2o_cont, kg, delg, ck)

      IMPLICIT NONE

      INCLUDE 'fl_radparams.inc'
      INCLUDE 'fl_special.inc'


c Function parameters
      REAL co2, ch4, n2o, f11, f12, f22
      REAL rhoair(nv1x), rhoh2o(nv1x), rhoo3(nv1x)
      REAL pres(nv1x), temp(nv1x), davg(nv1x), z(nv1x)
      LOGICAL h2o_cont
      INTEGER kg(mbx)
      REAL delg(mbx, mg)
      REAL ck(mbx,nvx,mg)
      INTEGER nlev


c Internal variables 
      INTEGER I, IB, IG

c Common blocks

c Atmospheric parameters
      REAL    pp, pt, pa, ph, po
      common /atmos/ pp(nv1x), pt(nv1x), pa(nv1x), ph(nv1x), po(nv1x)

c Line absorption optical thickness
      REAL    tg
      common /gas/ tg(nvx)

c Continuum absorption optical thickness
      REAL    tgm
      common /con/ tgm(nvx)

c H20 Continuum absorbtion now for entire LW (5-2200cm-1)
      INTEGER iwtas
      common /cont_tas/ iwtas

      logical lchou,lband6a,lpar,lray
      common /chou/ lchou,lband6a,lpar,lray

c select solar spectrum
      integer isolar_spectrum
      common /select_solar_spectra/ isolar_spectrum


c dz, layer thickness; computation is a little bit different from
c the Fu and Liou code: 
c  here,        dz = z(i) - z(i+1)
c  Fu and Liou, dz calculated from pressure difference
      REAL dz
      common /thick/ dz(nvx)

c mixing ratios of the well-mixed gases, CO2, CH4, and N2O
      REAL umco2, umch4, umn2o
      common /umcon/ umco2, umch4, umn2o


c CFCs are only considered if idkfr is set to 1
c (and this is a very expensive option, see 
c  option 'if (idkfr .eq. 1) then' )
      REAL cfc_conc(3)
      common /cfcs/ cfc_conc


c band weights
      real hk
      real hk1,  fk1o3,  sol_spect, fk1h2o
      real hk2,  c2h2o
      real hk3,  c3h2o
      real hk4,  c4h2o
      real hk5,  c5h2o
      real hk6,  c6h2o
      real hk7,  c7h2o
      real hk8,  c8h2o
      real hk9,  c9h2o
      real hk10, c10h2o, c10ch4, c10n2o
      real hk11, c11h2o, c11ch4, c11n2o
      real hk12, c12o3,  c12h2o
      real hk13, c13h2o
      real hk14, c14hca, c14hcb
      real hk15, c15hca, c15hcb
      real hk16, c16h2o
      real hk17, c17h2o
      real hk18, c18h2o

      common /band1/  hk1(10),fk1o3(10),sol_spect(0:7),fk1h2o(10)
      common /band2/  hk2(8),c2h2o(3,11,8)
      common /band3/  hk3(12),c3h2o(3,11,12)
      common /band4/  hk4(7),c4h2o(3,11,7)
      common /band5/  hk5(12),c5h2o(3,11,12)
      common /band6/  hk6(5),c6h2o(3,11,5)
      common /band7/  hk7(2),c7h2o(3,19,2)
      common /band8/  hk8(3),c8h2o(3,19,3)
      common /band9/  hk9(4),c9h2o(3,19,4)
      common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
      common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
      common /band12/ hk12(5),c12o3(3,19,5),c12h2o(3,19)
      common /band13/ hk13(2),c13h2o(3,19,2)
      common /band14/ hk14(10),c14hca(3,19,10),c14hcb(3,19,10)
      common /band15/ hk15(12),c15hca(3,19,12),c15hcb(3,19,12)
      common /band16/ hk16(7),c16h2o(3,19,7)
      common /band17/ hk17(7),c17h2o(3,19,7)
      common /band18/ hk18(8),c18h2o(3,19,8)



c Window K's option
c If idkfr is changed, the number of subbands in bands 11, 12, and 13
c is changed (see below); this implies changes at other places in 
c libRadtran, e.g. the band weights data/correlated_k/fu/quad.dat 
c and the photon fraction files, data/correlated_k/fu/x_solar.dat, etc.
      INTEGER idkfr
      common /dkfrwn/ idkfr
      data idkfr /0/

c H20 Continuum absorbtion now for entire LW (5-2200cm-1)
      data iwtas /3/

c DO NOT CHANGE lchou and lpar because these options require 
c action within the solver, too (see rad_0100.f)
      data lchou   /.FALSE./
      data lband6a /.FALSE./

c lpar; PAR water absorption in band 1
      data lpar    /.TRUE./

c lray; if TRUE, the Rayleigh optical depth changes with subband
c which is neccessary because with changing ozone absorption in band 1 
c the effective wavelength for Rayleigh changes.
      data lray    /.TRUE./

c select solar spectrum
c be sure that this option equals the data stored in 
c data/solar_flux/fu because the flux is read from there!
c For consistency with the Fu and Liou code, the extraterrestrial 
c irradiance in data/solar_flux/fu has been scaled to a solar 
c constant of 1368 W/m2; 
      data isolar_spectrum /1/



      kg(1)  = 10
      kg(2)  = 8
      kg(3)  = 12
      kg(4)  = 7
      kg(5)  = 12
      kg(6)  = 5
      kg(7)  = 2
      kg(8)  = 3
      kg(9)  = 4
      kg(10) = 4
      kg(11) = 3
      kg(12) = 5
      kg(13) = 2
      kg(14) = 10
      kg(15) = 12
      kg(16) = 7
      kg(17) = 7
      kg(18) = 8

c  Depending on Window K's option, the number of subbands needs changes...

      if (idkfr .eq. 0) then
         kg(11) = 3
         kg(12) = 5
         kg(13) = 2
      endif

      if (idkfr .eq. 1) then
         kg(11) = 80
         kg(12) = 40
         kg(13) = 10
      endif

      if (idkfr .eq. 2) then
         kg(11) = 5
         kg(12) = 5
         kg(13) = 2
      endif

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
        write (*,*) '       and recompile the package (ckdfuf.f)'
        STOP
      end if

c copy atmospheric data to common blocks
      DO I = 1, nv1
        PP(I) = PRES(I)
        PT(I) = TEMP(I)
        PA(I) = DAVG(I)
        PH(I) = RHOH2O(I)/RHOAIR(I)
        PO(I) = RHOO3(I)/RHOAIR(I)
      ENDDO

      UMCO2 = CO2
      UMCH4 = CH4
      UMN2O = N2O

      cfc_conc(1) = F11
      cfc_conc(2) = F12
      cfc_conc(3) = F22


c initialize layer widths (replacement for thicks())
      DO I = 1, nv
         dz(i) = z(i) - z(i+1)
      ENDDO 

c initialize Rayleigh cross sections
c DO NOT REMOVE THIS! It is required in later calls to rayle()
      call rayle2
           
c initialization; very important; sets e.g. the band weights      
      call ckd1_init (isolar_spectrum, lray)


c Loop over bands
      DO ib = 1, mb

c Get the water vapor continuum for band ib
         if (h2o_cont) then
            call gascon_ckd_parm(ib, h2o_cont)
         endif

c Loop over the g's for this band
         DO IG = 1, kg(ib)
c Get the layer optical depths and compute the absorption coefficients
            CALL GASES (ib, IG, HK)
            DO I = 1, nv
               
               if (h2o_cont) then
                  ck(ib,i,ig) = (tg(I) + tgm(I)) / dz(I)
               else
                  ck(ib,i,ig) = tg(I) / dz(I)
               endif
               
               IF (ck(ib,i,ig)  .LT. 0.0) THEN
                  WRITE (0,'(1X,A,A,E10.3)') 'Something is wrong,',
     $                 ' negative absorption: ',ck(ib,i,ig) 
                  WRITE (0,'(1X,A,I2,A,I3,A,I2)')
     $                 'Band: ',ib, '   Level:',I, '   k:',IG
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

c Copy band weights to final destination      
      DO IG = 1, kg(1)
         delg(1,ig) = hk1(ig)
      ENDDO

      DO IG = 1, kg(2)
         delg(2,ig) = hk2(ig)
      ENDDO

      DO IG = 1, kg(3)
         delg(3,ig) = hk3(ig)
      ENDDO

      DO IG = 1, kg(4)
         delg(4,ig) = hk4(ig)
      ENDDO

      DO IG = 1, kg(5)
         delg(5,ig) = hk5(ig)
      ENDDO

      DO IG = 1, kg(6)
         delg(6,ig) = hk6(ig)
      ENDDO

      DO IG = 1, kg(7)
         delg(7,ig) = hk7(ig)
      ENDDO

      DO IG = 1, kg(8)
         delg(8,ig) = hk8(ig)
      ENDDO

      DO IG = 1, kg(9)
         delg(9,ig) = hk9(ig)
      ENDDO

      DO IG = 1, kg(10)
         delg(10,ig) = hk10(ig)
      ENDDO

      DO IG = 1, kg(11)
         delg(11,ig) = hk11(ig)
      ENDDO
      
      DO IG = 1, kg(12)
         delg(12,ig) = hk12(ig)
      ENDDO
      
      DO IG = 1, kg(13)
         delg(13,ig) = hk13(ig)
      ENDDO

      DO IG = 1, kg(14)
         delg(14,ig) = hk14(ig)
      ENDDO

      DO IG = 1, kg(15)
         delg(15,ig) = hk15(ig)
      ENDDO

      DO IG = 1, kg(16)
         delg(16,ig) = hk16(ig)
      ENDDO

      DO IG = 1, kg(17)
         delg(17,ig) = hk17(ig)
      ENDDO

      DO IG = 1, kg(18)
         delg(18,ig) = hk18(ig)
      ENDDO


      END
