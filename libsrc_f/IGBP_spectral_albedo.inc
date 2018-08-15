C *** MODULE SPECTRAL_DAT ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 16:25:20 9/11/06
C...Switches:                     
!
!
      INTEGER NIGBP
      PARAMETER (NIGBP = 20)
      INTEGER NBAND
      PARAMETER (NBAND = 15)
!
! *** Spectral Reflectances
!
      REAL SPECALBIGBP(15,20)
!
!
! *** Scene dependent solar zenith adjustment factor
!
!
      REAL D(20)                                 ! (20) SEA ICE
!
! *** Longwave Emmissivities in 12 Fu bands
!
      REAL EMISSTBL(12,20)
!
! Broadband and Window (8-12microns) Emissivities
! Integrated values of spectral emissivities, not
! currently used by Fu & Liou code.
!
!     Data e_bb/0.996,  &! ( 1) EVERGREEN NEEDLE FOR
!               0.996,  &! ( 2) EVERGREEN BROAD FOR
!               0.990,  &! ( 3) DECIDUOUS NEEDLE FOR
!               0.990,  &! ( 4) DECIDUOUS BROAD FOR
!               0.993,  &! ( 5) MIXED FOREST
!               0.984,  &! ( 6) CLOSED SHRUBS
!               0.954,  &! ( 7) OPEN/SHRUBS(DESERT)
!               0.993,  &! ( 8) WOODY SAVANNA
!               0.993,  &! ( 9) SAVANNA
!               0.993,  &! (10) GRASSLAND
!               0.992,  &! (11) WETLAND
!               0.981,  &! (12) CROPLAND
!               1.000,  &! (13) URBAN
!               0.983,  &! (14) CROP MOSAIC
!               1.000,  &! (15) ANTARCTIC SNOW
!               0.941,  &! (16) BARREN/DESERT
!               0.991,  &! (17) OCEAN WATER
!               0.992,  &! (18) TUNDRA
!               0.988,  &! (19) FRESH SNOW
!               0.979/   ! (20) SEA ICE
 
!     Data e_wn/0.990,  &! ( 1) EVERGREEN NEEDLE FOR
!               0.990,  &! ( 2) EVERGREEN BROAD FOR
!               0.978,  &! ( 3) DECIDUOUS NEEDLE FOR
!               0.978,  &! ( 4) DECIDUOUS BROAD FOR
!               0.984,  &! ( 5) MIXED FOREST
!               0.955,  &! ( 6) CLOSED SHRUBS
!               0.897,  &! ( 7) OPEN/SHRUBS(DESERT)
!               0.982,  &! ( 8) WOODY SAVANNA
!               0.982,  &! ( 9) SAVANNA
!               0.982,  &! (10) GRASSLAND
!               0.984,  &! (11) WETLAND
!               0.982,  &! (12) CROPLAND
!               1.000,  &! (13) URBAN
!               0.983,  &! (14) CROP MOSAIC
!               1.000,  &! (15) ANTARCTIC SNOW
!               0.869,  &! (16) BARREN/DESERT
!               0.986,  &! (17) OCEAN WATER
!               0.981,  &! (18) TUNDRA
!               0.988,  &! (19) FRESH SNOW
!               0.979/   ! (20) SEA ICE
!
      REAL BB(20)
      COMMON /SPECTRAL_DAT/SPECALBIGBP, EMISSTBL, BB
C *** END MODULE SPECTRAL_DAT ***
