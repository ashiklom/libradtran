/*--------------------------------------------------------------------
 * $Id: uvspec.h 3311 2017-12-07 16:19:50Z bernhard.mayer $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#ifndef __uvspec_h
#define __uvspec_h

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mystic.h"
#if HAVE_LIDAR
#include "lidar.h"
#endif
#include <time.h>

/* moved to Makeconf.in and lex_starter.l */
/*#define VERSION "1.6-beta"*/

/* Indices for the reading of the atmosphere files */
#define ATMFILE_ZD       0
#define ATMFILE_PRESS    1
#define ATMFILE_TEMPER   2
#define ATMFILE_AIR      3
#define ATMFILE_O3       4
#define ATMFILE_O2       5
#define ATMFILE_H2O      6
#define ATMFILE_CO2      7
#define ATMFILE_NO2      8
#define ATMFILE_SO2      9

/* Atmospheric molecules that may be included in spectral calculations */
#define MOL_NN       15
#define MOL_AIR       0
#define MOL_O3        1
#define MOL_O2        2
#define MOL_H2O       3
#define MOL_CO2       4
#define MOL_NO2       5
#define MOL_BRO       6
#define MOL_OCLO      7
#define MOL_HCHO      8
#define MOL_O4        9
#define MOL_SO2      10
#define MOL_CH4      11
#define MOL_N2O      12
#define MOL_CO       13
#define MOL_N2       14

/* Possible units of total molecule column */
#define MOL_UNIT_CM_2 1
#define MOL_UNIT_DU   2
#define MOL_UNIT_PPM  3
#define MOL_UNIT_MM   4

/* Possible units of density profiles */
#define CM_3          1    /* per cubic centimeter */
#define M_3           2    /* per cubic meter      */
#define MMR           3    /* mass mixing ratio    */
#define VMR           4    /* volumn mixing ratio  */
#define RH            5    /* relative humidity    */

#define RS_NO_DATA    0   /* no radiosonde data */
#define RS_FROM_FILE  1   /* radiosonde data source */
#define RS_FROM_ECMWF 2   /* ECMWF data converted with CDO */
#define RS_FROM_ECHAM 3   /* ECHAM data*/


#define FN_NN                 40  /* Number of general files, increase 
				     when adding new general files */

#define FN_PATH                0
#define FN_ALBEDO              1
#define FN_ALBEDO_LIB_PATH     2
#define FN_ALBEDO_LIB_WVL      3
#define FN_ALBEDO_MAP          4
#define FN_SURFACE_TYPE_MAP    5
#define FN_SURFACE_TEMP_MAP    6
#define FN_RPV                 7
#define FN_RPV_LIB_PATH        8
#define FN_ALTITUDE_MAP        9
#define FN_ATMOSPHERE         10
#define FN_RADIOSONDE         11
#define FN_ECMWF              12
#define FN_ECMWF_WIND_MAP     13
#define FN_U10_MAP            14
#define FN_SALINITY_MAP       15
#define FN_PIGMENTS_MAP       16
#define FN_ECHAM              17
#define FN_EXTRATERRESTRIAL   18
#define FN_SZA                19
#define FN_SATGEOMETRY        20
#define FN_WLTRANS            21
#define FN_WLBANDS            22
#define FN_IC_REFF            23
#define FN_CLOUD_FRACTION     24
#define FN_CLOUD_FRACTION_MAP 25
#define FN_SLITFUNCTION       26
#define FN_FILTERFUNCTION     27
#define FN_SPLINE             28
#define FN_MOL_TAU_SCA        29
#define FN_MOL_TAU_ABS        30
#define FN_REFIND             31
 //#define FN_CK_GENERIC         32 /* filename for Input.ck_scheme = CK_GENERIC; included Input.ck_scheme_filename instead */
#define FN_OUTFILE            33
#define FN_RANDOMSTATUS       34
#define FN_FLUORESCENCE       35
#define FN_REPTRAN            36
#define FN_ROSSLI             37
#define FN_HAPKE              38
#define FN_RAYTRACING         39

#define FN_AER_NN            9  /* Number of aerosol related files, increase 
				   when adding new aerosol files */
#define FN_AER_TAU           0
#define FN_AER_GG            1
#define FN_AER_SSA           2
#define FN_AER_MOMENTS       3
#define FN_AER_REF           4
#define FN_AER_SIZ           5
#define FN_AER_EXPLICIT      6
#define FN_AER_SPECIES       7
#define FN_AER_SPECIES_LIB   8

#define FN_MC_NN              16 /* Number of mystic related files, increase
			            when adding new mystic files */
#define FN_MC_PHOTONS         0
#define FN_MC_BASENAME        1
#define FN_MC_ALBEDO          2
#define FN_MC_ALBEDO_TYPE     3
#define FN_MC_ALBEDO_SPECTRAL 4
#define FN_MC_TEMPERATURE     5
#define FN_MC_ELEVATION       6
#define FN_MC_RPV_TYPE        7
#define FN_MC_RPV_SPECTRAL    8
#define FN_MC_ROSSLI          9
#define FN_MC_LIDAR          10
#define FN_MC_UMU            11
#define FN_MC_SUNSHAPE_FILE  12
#define FN_MC_AERIS          13
#define FN_MC_AIRMASS        14
#define FN_MC_AMBRALS        15

#define CK_CRS             0
#define CK_KATO            1
#define CK_KATO2           2
#define CK_KATO2_96        3
#define CK_FU              4
#define CK_AVHRR_KRATZ     5
#define CK_FILE	           6   /* CK_GENERIC before; renamed to CK_FILE for new input option name */
#define CK_LOWTRAN         7
#define CK_SBDART          7   /* Same as LOWTRAN; needed for new input option name */
#define CK_RAMAN           8   /* This is not a CK-distribution, but a way to do many RTE solutions for one wavelength, AK */
#define CK_REPTRAN         9
#define CK_REPTRAN_CHANNEL 10
#define CK_KATO2ANDWANDJI  12

#define CK_ABS_NN	   7   /* implemented for ck_lowtran_absorption (old: absorption); Merge with molecules for atmospheric molecules for spectral calculations (MOL_NN) */
#define CK_ABS_O4	   0
#define CK_ABS_N2	   1
#define CK_ABS_CO	   2
#define CK_ABS_SO2	   3
#define CK_ABS_NH3	   4
#define CK_ABS_NO	   5
#define CK_ABS_HNO3	   6

#define CK_H2OCONT_OFF     0
#define CK_H2OCONT_ON      1
#define CK_H2OCONT_V2_1    1
#define CK_H2OCONT_V2_4    2

#define REPTRAN_OPTION_NONE    0
#define REPTRAN_OPTION_FINE    1
#define REPTRAN_OPTION_MEDIUM  2
#define REPTRAN_OPTION_COARSE  3

#define WLGRID_NONE    0
#define WLGRID_CK      1
#define WLGRID_USER    2
#define WLGRID_BANDS   3
#define WLGRID_MOLABS  4
#define WLGRID_UVSPEC  5

#define SRC_NONE    0
#define SRC_SOLAR   1
#define SRC_THERMAL 2
#define SRC_BLITZ   3
#define SRC_LIDAR   4

#define CRS_MOL_NN       3
#define CRS_MOL_RAYLEIGH 0
#define CRS_MOL_O3       1
#define CRS_MOL_NO2      2

#define CRS_MODEL_BODHAINE	1
#define CRS_MODEL_NICOLET	2
#define CRS_MODEL_PENNDORF	3
#define CRS_MODEL_BODHAINE29	4
#define CRS_MODEL_BASS_AND_PAUR	5
#define CRS_MODEL_MOLINA	6
#define CRS_MODEL_DAUMONT	7
#define CRS_MODEL_BOGUMIL	8
#define CRS_MODEL_BURROWS	9


#define CO2_YOSHINO  1

#define O2_DEFAULT   1

#define SO2_BREMEN   1

#define BRO_IASB     1

#define OCLO_WAHNER  1

#define HCHO_IASB  1

#define O4_GREENBLATT  1

#define RAMAN_NONE         0
#define RAMAN_CALC         1   /* Raman scattering included */
#define RAMAN_DEFAULT      1   /* The default Raman scattering cross section */

#define RAYLEIGH_NONE      0
#define RAYLEIGH_CALC      1
#define RAYLEIGH_FILE_MONO 2
#define RAYLEIGH_FILE_SPEC 3

#define MOLABS_NONE      0
#define MOLABS_CALC      1
#define MOLABS_FILE_MONO 2
#define MOLABS_FILE_SPEC 3
#define MOLABS_LOOKUP    4

#define MOLABS_SRC_NONE   0
#define MOLABS_SRC_CALC   1
#define MOLABS_SRC_LOOKUP 2

#define MX_NN  9  
#define MX_O2  0
#define MX_H2O 1
#define MX_CO2 2
#define MX_NO2 3
#define MX_CH4 4
#define MX_N2O 5
#define MX_F11 6
#define MX_F12 7
#define MX_F22 8

#define CLD_ZD     0    
#define CLD_LWC    1
#define CLD_EFFR   2

/*merging brdf Hapke options*/
#define BRDF_HAPKE_NN	3
#define BRDF_HAPKE_B0	0
#define BRDF_HAPKE_H	1
#define BRDF_HAPKE_W	2
/*merging brdf Ross-Li options*/
#define BRDF_ROSSLI_NN	3
#define BRDF_ROSSLI_ISO	0
#define BRDF_ROSSLI_VOL	1
#define BRDF_ROSSLI_GEO	2
/*merging brdf Ross-Li options*/
#define BRDF_ROSSLI_HOTSPOT_OFF 0
#define BRDF_ROSSLI_HOTSPOT_ON  1
/*merging brdf rpv options*/
#define BRDF_RPV_NN	7
#define BRDF_RPV_K	0
#define BRDF_RPV_RHO0	1
#define BRDF_RPV_THETA	2
#define BRDF_RPV_SIGMA	3
#define BRDF_RPV_T1	4
#define BRDF_RPV_T2	5
#define BRDF_RPV_SCALE	6
/* moved to cdisort.c */
//#define BRDF_NONE   0
//#define BRDF_RPV    1  /* don't change these numbers as they are */
//#define BRDF_CM     2  /* used by Fortran code which of course   */
//#define BRDF_AMB    3  /* has no access to this header file      */

#define BPDF_NONE   0
#define BPDF_TSANG  1

#define RPV_CONSTANT       0
#define RPV_FROM_RPV_FILE  1
#define RPV_IGBP_LIBRARY   2
#define RPV_USER_LIBRARY   3

#define ROSSLI_CONSTANT          0
#define ROSSLI_AMBRALS_CONSTANT  1
#define ROSSLI_FROM_ROSSLI_FILE  2
#define ROSSLI_FROM_AMBRALS_FILE 3

#define HAPKE_CONSTANT       0
#define HAPKE_FROM_HAPKE_FILE  1

#define SOLVER_FDISORT1    1
#define SOLVER_SDISORT     2
#define SOLVER_FTWOSTR     3
//#define TWOSTRPP         4 /*bettina: obsolete, == CTWOSTR */
#define SOLVER_SOS         5
#define SOLVER_MONTECARLO  6
#define SOLVER_MYSTIC      6
#define SOLVER_POLRADTRAN  7
#define SOLVER_FDISORT2    8
#define SOLVER_SPSDISORT   9
//#define QDISORT         10 /* arve: made obsolete, cdisort has general source */
#define SOLVER_TZS        11
#define SOLVER_SSS        12
#define SOLVER_SSSI       13
#define SOLVER_NULL       16
#define SOLVER_RODENTS    17 /* ulrike: Robert Buras' twostream */
#define SOLVER_TWOSTREBE  18 /* ulrike: Bernhard Mayers twostream */
#define SOLVER_TWOMAXRND  19 /* ulrike: Bernhard Mayers twostream */
#define SOLVER_TWOSTR     20
#define SOLVER_DISORT     21
#define SOLVER_SSLIDAR    22

#define SSLIDAR_NN	  5 /*following SSLIDAR for new option names sslidar ... */
#define SSLIDAR_AREA	  0
#define SSLIDAR_E0	  1
#define SSLIDAR_EFF	  2
#define SSLIDAR_POSITION  3
#define SSLIDAR_RANGE	  4

#define DISORT_ICM_OFF     0
#define DISORT_ICM_PHASE   1
#define DISORT_ICM_MOMENTS 2

#define MC_SENSORPOSITION_NONE      0
#define MC_SENSORPOSITION_CARTESIAN 1
#define MC_SENSORPOSITION_SPHERICAL 2

#define MC_RIS_NN		2	/*new input option mc_ris factor/optical_depth value (old: mc_ris_factor,mc_ris_optical_depth) */
#define MC_RIS_FACTOR		0
#define MC_RIS_OPTICAL_DEPTH	1

#define POLRADTRAN_NN	        3
#define POLRADTRAN_AZIORDER     0
#define POLRADTRAN_NSTOKES      1
#define POLRADTRAN_SRC_CODE     2

#define SDISORT_NN	   3
#define SDISORT_ICHAPMAN   0
#define SDISORT_NSCAT	   1
#define SDISORT_NREFRAC	   2

#define PAN_MODE_NONE      0
#define PAN_MODE_CAMERA    1
#define PAN_MODE_SATELLITE 2

#define PAN_FLAG_NN                        6
#define PAN_FLAG_DISTR_PHOTONS_OVER_PIXEL  0
#define PAN_FLAG_NO_PIXEL                  1
#define PAN_FLAG_QUICKLOOK                 2
#define PAN_FLAG_WITH_DIRECT_RAD           3
#define PAN_FLAG_WEIGHT_WITH_COS           4
#define PAN_FLAG_CIRCUMSOLAR_VAR_RED       5

#define PAN_ALIGNMENT_NONE   0
#define PAN_ALIGNMENT_SUN    1
#define PAN_ALIGNMENT_MU     2
#define PAN_ALIGNMENT_ZENITH 3

#define IC_FU_REFF_DEF          0  /* new variable for option ic_fu reff_def */
#define IC_FU_DELTASCALING      1  /* new variable for option ic_fu deltascaling */

#define MODIFY_TYPE_NN          2  /* number of mod_id2 (modify set/scale )*/
#define MODIFY_VAR_NN           4  /* number of mod_id1 (modify gg/ssa/tau/tau550)*/
#define MODIFY_TYPE_SET         0  /* modify set id (modify profile optical properties)*/
#define MODIFY_TYPE_SCALE       1  /* modify set id (modify profile optical properties)*/
#define MODIFY_VAR_GG           0  /* modify gg  id (modify profile optical properties)*/
#define MODIFY_VAR_SSA          1  /* modify ssa id (modify profile optical properties)*/
#define MODIFY_VAR_TAU          2  /* modify tau id (modify profile optical properties)*/
#define MODIFY_VAR_TAU550       3  /* modify tau550 id (modify profile optical properties)*/

#define PROP_NONE         0   /* no caoth (= clouds and other things)      */
#define PROP_HU           1   /* "wc_file" and Hu and Stamnes [1993]       */
#define PROP_ECHAM4       2   /* "wc_file" and ECHAM4                      */
#define PROP_MIE          3   /* "wc_file" and "wc_properties mie"         */
#define PROP_FILE         4   /* "wc_file" and "wc_properties <filename>"  */
#define PROP_EXPLICIT     5   /* "wc_files"                                */
#define PROP_FU           6   /* "ic_file" and Fu et al. [1996/98]         */
#define PROP_KEY          7   /* "ic_file" and Key et al. [2002]           */
#define PROP_YANG         8   /* "ic_file" and Yang/Key/Mayer              */
#define PROP_IC_MIE       9   /* "ic_file" and "ic_properties mie"         */
#define PROP_BAUM        10   /* "ic_properties baum"                      */
#define PROP_BAUM_V36    11   /* "ic_properties baum_v36  "                */
#define PROP_HEY         12   /* "ic_properties hey"                       */
#define PROP_YANG2013    13   /* "ic_properties yang2013"                  */
#define PROP_RAYTRACING  14   /* "ic_properties raytracing"                */

#define REFF_FIXED       0   /* fixed r_eff in all heights                */
#define REFF_OU          1   /* Ou and Liou, 1995, Atm Res. 35, 127-138   */
#define REFF_FILE        2   /* read r_eff from file                      */

#define PROCESS_NONE        0
#define PROCESS_SUM         1
#define PROCESS_RGB         2
#define PROCESS_RGBNORM     3
#define PROCESS_INT         4
#define PROCESS_RAMAN       5

#define UNIT_NOT_DEFINED      -999
#define UNIT_PER_NM              0
#define UNIT_PER_CM_1            1
#define UNIT_PER_BAND         2

#define HEAT_NONE     0   /* no heating rate calculation                                               */
#define HEAT_LOCAL    1   /* heating rates from (k_abs * actinic fluxes)                               */
#define HEAT_LAYER_FD 2   /* heating rates calculated with Forward  Difference of the flux over layers */
#define HEAT_LAYER_CD 3   /* heating rates calculated with Centered Difference of the flux over layers */

/* main switch for interpolation onto zout and altitude levels */
#define Z_INTERPOLATE_LOGLIN  0      /* no interpolation                   */
#define Z_INTERPOLATE_LINMIX  1      /* linear mixing ratio interpolation  */
#define Z_INTERPOLATE_SPLINE  2      /* spline interpolation               */

/* assumptions for interpolation of profiles,   */
/* calculation of effective densities for       */
/* optical properties and columns               */
#define INTERP_METHOD_SPLINE     0  /* spline interpolation (columns not consistent)               */
#define INTERP_METHOD_LINEAR     1  /* linear profile                                              */
#define INTERP_METHOD_LOG        2  /* logarithmic profile                                         */
#define INTERP_METHOD_LINMIX     3  /* linear mixing ratio (only for gas densities other than air) */
#define INTERP_METHOD_LOG_SPLINE 4  /* first take logarithm, than spline (columns not consistent)  */

#define NO_ZOUT_INTERPOLATE    0  /* used for zout_interpolate */ 
#define ZOUT_INTERPOLATE       1  /* used for zout_interpolate */

#define TIME_NEAREST_DATE      0  /* used for time_interpolate */ 
#define TIME_INTERPOLATION     1  /* used for time_interpolate */

#define TIMEZONE_UTC        0
#define TIMEZONE_LAT        1

#define ZOUT_MYSTIC  -999
#define ZOUT_SURFACE -666
#define ZOUT_TOA      666
#define ZOUT_CPT      699

#define OUTPUT_USER_NONE             0
#define OUTPUT_USER_WAVE             1
#define OUTPUT_USER_WAVE_MAX         2
#define OUTPUT_USER_WAVE_MIN         3
#define OUTPUT_USER_WAVENUMBER       4
#define OUTPUT_USER_WAVENUMBER_MAX   5
#define OUTPUT_USER_WAVENUMBER_MIN   6
#define OUTPUT_USER_ZOUT_SUR         7
#define OUTPUT_USER_ZOUT_SEA         8
#define OUTPUT_USER_Z_SUR            9
#define OUTPUT_USER_P               10
#define OUTPUT_USER_T               11
#define OUTPUT_USER_T_D             12
#define OUTPUT_USER_T_SUR           13
#define OUTPUT_USER_THETA           14
#define OUTPUT_USER_THETA_E         15
#define OUTPUT_USER_EDIR            16
#define OUTPUT_USER_EGLO            17
#define OUTPUT_USER_EDN             18
#define OUTPUT_USER_EUP             19
#define OUTPUT_USER_ENET            20
#define OUTPUT_USER_ESUM            21
#define OUTPUT_USER_FDIR            22
#define OUTPUT_USER_FGLO            23
#define OUTPUT_USER_FDN             24
#define OUTPUT_USER_FUP             25
#define OUTPUT_USER_F               26
#define OUTPUT_USER_UU              27
#define OUTPUT_USER_UDIR            28
#define OUTPUT_USER_UGLO            29
#define OUTPUT_USER_UDN             30
#define OUTPUT_USER_UUP             31
#define OUTPUT_USER_U               32
#define OUTPUT_USER_ALB             33
#define OUTPUT_USER_SZA             34
#define OUTPUT_USER_C_P             35
#define OUTPUT_USER_HEAT            36
#define OUTPUT_USER_EMIS            37 
#define OUTPUT_USER_ABS             38 
#define OUTPUT_USER_W_RAD           39
#define OUTPUT_USER_M_RAD           40
#define OUTPUT_USER_N               41
#define OUTPUT_USER_RHO             42
#define OUTPUT_USER_VMR             43
#define OUTPUT_USER_MMR             44
#define OUTPUT_USER_RH              45
#define OUTPUT_USER_RH_ICE          46
#define OUTPUT_USER_CLWC            47
#define OUTPUT_USER_CLWD            48
#define OUTPUT_USER_TCLW            49
#define OUTPUT_USER_REFF_WAT        50     
#define OUTPUT_USER_CIWC            51
#define OUTPUT_USER_CIWD            52
#define OUTPUT_USER_TCIW            53        
#define OUTPUT_USER_REFF_ICE        54
#define OUTPUT_USER_TCW             55       
#define OUTPUT_USER_CC              56
#define OUTPUT_USER_TCC             57 
#define OUTPUT_USER_CLOUDS          58 
#define OUTPUT_USER_WIND_U          59
#define OUTPUT_USER_WIND_V          60
#define OUTPUT_USER_WIND_W          61
#define OUTPUT_USER_DTDX            62
#define OUTPUT_USER_DTDY            63
#define OUTPUT_USER_DTDZ            64
#define OUTPUT_USER_HEAT_AD_X       65
#define OUTPUT_USER_HEAT_AD_Y       66
#define OUTPUT_USER_HEAT_AD_Z       67
#define OUTPUT_USER_HEAT_AD         68
#define OUTPUT_USER_SSLIDAR_NPHOT   69
#define OUTPUT_USER_SSLIDAR_NPHOT_Q 70
#define OUTPUT_USER_SSLIDAR_RATIO   71
#define OUTPUT_USER_SPHER_ALB       72

#define OUTPUT_FORMAT_ASCII          0
#define OUTPUT_FORMAT_MATRIX         1
#define OUTPUT_FORMAT_FLEXSTOR       2
#define OUTPUT_FORMAT_NETCDF         3
#define OUTPUT_FORMAT_SAT_PICTURE    4

#define OUTLEVEL_ZOUT_ABOVE_SUR           0
#define OUTLEVEL_ZOUT_ABOVE_SEA           1
#define OUTLEVEL_PRESS                    2
#define OUTLEVEL_ATM_LEVELS               3
#define OUTLEVEL_ALL_LEVELS               4
#define OUTLEVEL_MODEL_LEVELS             5
#define OUTLEVEL_MODEL_LAYERS             6
#define OUTLEVEL_MODEL_LEVELS_AND_LAYERS  7

#define OUTCAL_ABSOLUTE        0
#define OUTCAL_THERMAL         1
#define OUTCAL_TRANSMITTANCE   2
#define OUTCAL_REFLECTIVITY    3
#define OUTCAL_BRIGHTNESS      4

#define ISCCP_WATER     0
#define ISCCP_ICE       1

#define MAXAWVN        47

#define FALSE           0
#define NO              0
#define TRUE            1
#define YES             1
#define SWITCH_OFF     0  /* new variable for option on/off */
#define SWITCH_ON      1  /* new variable for option on/off */

#define ALBEDO_CONSTANT             0 
#define ALBEDO_FROM_ALBEDO_FILE     1
#define ALBEDO_FROM_ALBEDO_MAP      2
#define ALBEDO_FROM_EMISSIVITY_MAP  3
#define ALBEDO_IGBP_LIBRARY         4
#define ALBEDO_USER_LIBRARY         5

#define FLUORESCENCE_CONSTANT                 0 
#define FLUORESCENCE_FROM_FLUORESCENCE_FILE   1

#define ALT_NOT_DEFINED  0
#define ALT_DIRECT_INPUT 1
#define ALT_FROM_MAP     2
#define ALT_FROM_ATM     3

#define NO_CAOTH		0
#define CAOTH_FROM_1D		1
#define CAOTH_FROM_IPA_FILES	2
#define CAOTH_FROM_3D		3
#define CAOTH_FROM_MOMENTS	4
#define CAOTH_FROM_ECMWF        5
#define CAOTH_FROM_ECHAM        6

#define DIM_NN              2
#define DIM_1D              0
#define DIM_3D              1

#define TIPA_NONE           0
#define TIPA_DIRDIFF        1   /* ulrike to choose between TICA where the whole 3d-cloud-matrix is tilted and the new input */
#define TIPA_DIR            2   /*    and 3dIPA where only the direct radiation is calculated with "tilted clouds" */
#define TIPA_DIR3D          3   /*    special version, MC photon is 3d for direct and ipa afterwards */

/* cloud overlap */
#define CLOUD_OVERLAP_RAND     0
#define CLOUD_OVERLAP_MAX      1
#define CLOUD_OVERLAP_MAXRAND  2
#define CLOUD_OVERLAP_OFF      3

#define SZA_BY_TIME_AND_LOCATION  0
#define SZA_DIRECT_INPUT          1
#define SZA_FROM_SZA_FILE         2
#define SZA_ECHAM                 3

#define LATITUDE_SIGNUM_N	  1
#define LATITUDE_SIGNUM_S	 -1
#define LONGITUDE_SIGNUM_E	  1
#define LONGITUDE_SIGNUM_W	 -1

#define NOT_DEFINED_DOUBLE  -999.0
#define NOT_DEFINED_FLOAT   -999.0
#define NOT_DEFINED_INTEGER -999

/* in the style of netCDF types NC_xxx */
#define TYPE_CHAR    2
#define TYPE_SHORT   3
#define TYPE_INT     4
#define TYPE_FLOAT   5
#define TYPE_DOUBLE  6

/* additional error numbers for netCDF output */
#define ERROR_READ_MISSING_VALUE  -501

/* define ascending or descending axis */
#define ASCENDING   0
#define DESCENDING  1

/* minimum number of photons for each band*/
#define MIN_MCPHOTONS 1000   

/* physical constants */
#define GAMMA 6.673E-11        /* +/- 0.010E-11  Newton's gravitational constant in N m2 kg-2                              */ 
                               /* accuracy of gamma is quite poor, as it is hard to measure                                */
                               /* recommended value from Fattori, M., Towards an atom interferometric determination of the */
                               /* Newtonian graviational constant, 2003, Physics Letters A 318 (2003), 184-191             */
#define BOLTZMANN  1.3806504e-23 /* boltzmann constant k, form wikipedia */
#define AVOGADRO   6.022142e23   /* 1 mol particles */
 

/* data of the planet */
/* CE: r_earth is now also an input option! Check whether this should be calculated*/
/* accurately based on r_earth input !*/
#define R_EARTH   6378000.0                       /*radius of the earth in m */
#define M_EARTH   5.97412E24                      /* mass of the earth in kg  */
#define G_SURFACE GAMMA*M_EARTH/(R_EARTH*R_EARTH) /* gravitational acceleration */

/* averaged data of "air" */
#define C_P_DRY_STD  1003.8   /* specific heat of dry air at constant pressure (J kg-1 K-1), standard T = 275 K */
#define C_V_DRY_STD   716.87  /* specific heat of dry air at constant volumn   (J kg-1 K-1), standard T = 275 K */ 
#define MOL_MASS_AIR 28.977   /* mass per mol particle in g / mol (from function ck.c) */
#define MOL_MASS_WV  18.015   /* mass per mol particle in g / mol (from function ck.c) */
#define R_AIR  BOLTZMANN*AVOGADRO*1000.0/MOL_MASS_AIR  /* gas constant for dry air */

/**********************************************/
/*            structure definitions           */
/**********************************************/

/**************************/
/* Monte Carlo STRUCTURES */
/**************************/

/* Monte Carlo backward structure */

typedef struct {
  int yes;  /* flag if backward is desired */

  int islower;   /* x-index of lower left  sample pixel */
  int jslower;   /* y-index of lower left  sample pixel */
  int isupper;   /* x-index of upper right sample pixel */
  int jsupper;   /* y-index of upper right sample pixel */

  int writeallpixels;
  int writeback;

  int output;      /* which quantity: edir, edn, eup, ...    */

  int absorption;  /* flag: calculation of absorption yes/no */

  int thermal_heating_method; /* flag: how to calculate thermal heating rates */
} mc_backward_struct;

/* Monte Carlo input structure */

typedef struct {
  long int photons;

  int     Nx_sample;
  int     Ny_sample;

  double  dx_sample;
  double  dy_sample;

  float   alpha;  /*opening angle for radiance sampling */

  int     Nr;   /* number of radii, for I3RC case 7 */
  int     Nt;   /* number of times, for I3RC case 7 */
  float   dt;   /* time interval,   for I3RC case 7 */

  int     tenstream;      // Tenstream solver by Fabian Jakub
  char   *tenstream_options;      //  number of streams, ... flag

  int     ipa;            /* MYSTIC in IPA mode                         */
  int     tipa;           /* MYSTIC in TIPA mode (ulrike)               */
  int     escape;         /* escape probabilities for radiances         */

  int     locest;         /* local estimator                            */
  int     jacobian;
  int     jacobian_std;
  int     reference_to_NN;         /* referencing to NN                 */

  int     ncirc;
  double  dcirc;

  int     vroom;          /* Variance Reduction Optimal Options Method  */
  int     maxscatters;    /* photon killed if number of scatters larger */
  int     minscatters;    /* local estimate only done after this        */
  int     minphotons;     /* minimum number of photons to simulate      */

  double *blitz_position; /* start and end position of linear blitz     */

  float   delta_scaling_mucut;
  int     delta_scaling_start;

  float   truncate;       /* truncate phase function                    */

  int     spectral_is;    /* spectral calculation with importance sampling  (ALIS)*/
  int     spectral_is_nwvl;
  float   *spectral_is_wvl; /* calculation wavelength for ALIS           */

  int     concentration_is; /* concentration importance sampling */

  int     boxairmass;     /* boxairmass factors                         */

  int    *spherical;      /* spherical flag                             */
  int     spherical3D_scene;
                          /* flag for restricted scene in spherical 3D  */
  double  spherical3D_scene_lon_min;
                          /* left longitude scene in spherical 3D       */
  double  spherical3D_scene_lon_max;
                          /* right longitude scene in spherical 3D      */
  double  spherical3D_scene_lat_min;
                          /* lower latitude scene in spherical 3D       */
  double  spherical3D_scene_lat_max;
                          /* upper latitude scene in spherical 3D       */

  int     DoLE;           /* Specular reflection for Ocean              */

  int     ixmin;          /* ranges of netcdf cloud to be used          */
  int     ixmax;
  int     iymin;
  int     iymax;

  int     polarisation;   /* polarisation flag                          */
  int     nstokes;        /* number of Stokes components                */
  int     polarisation_state;  /* initial polarisation state            */
 
  int     coherent_backscatter; /* coherent backscatter flag CB         */
  int     coherent_backscatter_lmode; /* cb lidar mode selection        */
  int     coherent_backscatter_off; /* use CB mode without CB result    */
  double  *ris;           /* ris [factor/optical_depth]                 */

  int     sample_cldprp;   /* FE: sampling of cloud properties flag     */

  int     bcond;          /* boundary conditions                        */

  int     reflectalways;  /* reflect each photon and weight with albedo */
  int     refraction;     /* enable refraction                          */

  int     sensorposition;   /* user-defined sensor position for         */
  double  sensorposition_x; /* sensors at point locations               */ 
  double  sensorposition_y;
  double  sensorposition_z;
  double  sensorposition_lon;
  double  sensorposition_lat;
  double  sensorposition_alt;
  int     panorama;         /* panorama flag (1:zenith mode,
			                      2:satellite mode)         */
  int     pan_alignment;    /* panorama flag, camera alignment          */
  int    *pan;              /* panorama flag (distr_photons_over_pixel,np_pixel,quicklook,with_direct_rad,weight_with_cos) */


  double *pan_angles;       /* panorama boundaries                      */

  int     panorama_forward; /* allsky (4pi) panorama calculated by forward MYSTIC */
  int     panorama_forward_Ntheta; /* number of theta angles */
  int     panorama_forward_Nphi;   /* number of phi angles */
  
  double  sun_radius;       /* angular radius of the sun disc in deg    */

  int     surfaceparallel;  /* irradiance parallel to the surface       */
                            /* instead of horizontal                    */

  int     sensordirection;    /* user-defined sampling direction for    */
  double  sensordirection_dx; /* specifically oriented sensors          */ 
  double  sensordirection_dy; 
  double  sensordirection_dz; 

  int     absorption;     /* calculate absorption for each box          */ 
  int     abs_unit;       /* output units for absorbed radiance         */

  char    **filename;

  int     azimuth;        /* new or old azimuth convention              */

  int     std;            /* calculate standard deviation               */
  int     randomseed;     /* seed for the random number generator       */
  int     readrandomseed;   /* read seed for the random number generator  */
  int     readrandomstatus; /* read status of the random number generator */

  int     spectral;

  int     allocate_umu_and_phi; /*Quix fix for mc_bw_umu_file, mc_panorama, mc_blitz_position; needs to be cleaned up;*/
  mc_backward_struct backward;  /* information for backward MYSTIC      */

  int visualize;          /* OpenGL visualization                       */
} mc_inp_struct;

/* Monte Carlo output structure */
/* (MYSTIC uses profiles in reverted order (the lowest level is 0)) */

typedef struct {
  float       *dt;
  float       *om;

  float       *g1;
  float       *g2;
  float       *ff;

  float       *ds;

  float       *re;
} mc_caoth_struct;


typedef struct {
  
  float       *dens;  /* air density needed to calculate spectral rayleigh coefficients*/
  
  int         nphamataer;
  int         *nmomaer;
  float       ***momaer;
  int         **nthetaaer;
  float       ***thetaaer;
  double      ***muaer;
  float       ***phaseaer;
  
  float       *z;
  float       ***temper;

  float       *refind;

  /*  mc_caoth_struct *caoth; */
  float       **dt;
  float       **om;

  float       **g1;
  float       **g2;
  float       **ff;

  float       **ds;

  float       **re;

  /* 3D molecular atmosphere */
  double       *****kabs3D; 
  double       *****ksca3D;
  
  alis_struct alis;          /* importance sampling    */
  sample_struct sample;      /* 3d sampling            */
  elevation_struct elev;     /* 2d elevation           */
  surftemp_struct surftemp;  /* 2d surface temperature */
} mc_out_struct;


/***********************************************/
/* Raditive Transfer Equation STRUCTURE        */
/* See disort.doc for description of variables */
/***********************************************/

typedef struct {
  int     deltam;
  int     fisot;
  int     nstr;  
  int     maxphi;
  int     maxumu;
  int     nphi;
  int     numu;
  int     disort_icm;
  int     ibcnd; 
  int     pseudospherical;
  int     solver;
  int     nprndis;
  int     reverse;
  int     sos_nscat;
  int*    prndis;
  int*    sdisort;
  int*    polradtran;
  char*   pol_quad_type;
  char    pol_stokes[4];
  double  pol_max_delta_tau;
  float*  phi;
  float*  umu;
  int*    cmuind; /* Indices for the cmu angles, used with Raman option */
  int*    umuind; /* Indices for the umu angles, used with Raman option */
  mc_inp_struct mc;
} rte_struct;


/* Caoth output STRUCTURE */

typedef struct {
  /* quantities as defined by the user */
  float *lwc;        /* z-profile of the liquid water content */
  float *lwc_zout;   /* z-profile of the liquid water content at zout levels */
  float tcw;         /* total column water */
  float *effr;       /* z-profile of the effective radius     */
  float *effr_zout;  /* z-profile of the effective radius at zout levels     */

  /* layer quantities */
  float *lwc_layer;  /* z-profile of the liquid water content */
  float *effr_layer; /* z-profile of the effective radius     */

  float effrmin;     /* minimum effective radius */
  float effrmax;     /* maximum effective radius */
  
  int   alloc;       /* flag, if memory is allocated          */
} caoth_microphys_struct;


/********************************************************/
/* optical properties of aerosol, water, and ice clouds */
/********************************************************/

typedef struct {
  float **dtau;     /* layer optical thickness [wavelength][layer]         */

  float **ff;       /* f  for double HG [wavelength][layer]                */
  float **g1;       /* g1 for double HG [wavelength][layer]                */
  float **g2;       /* g2 for double HG [wavelength][layer]                */

  float **ssa;      /* layer single scattering albedo [wavelength][layer]  */
  float **f;        /* delta-m scaling factor                              */

  int   nlam;       /* number of wavelengths                               */
	
  float ****moment; /* moments of the phase function for each              */
                    /* wavelength and layer and polarisation matrix element*/
  int   nlev;       /* number of levels                                    */

  int   **nmom;     /* number of moments for each wavelength and layer     */
  int   nphamat;    /* number of phase matrix elements                     */

  int ***ntheta;
  float ****theta;
  double ****mu;
  float ****phase;

  float **dscale;   /* layer delta-scaling factor                     TZds */

  int alloc;        /* flag, if memory is allocated;                       */
} optprop_struct;


/* TIPA structure (ulrike)*/

typedef struct {
  /* for tipa dir (tilting) */
  double **zsorted;/* [tilting-level][level=0..totlev-1]; contains z-levels */
  double *level;   /* tilting-levels */
  int ***msorted;  /* [tilting-level][level=0..totlev-1][0..2 i.e. nx,ny or nz]; contains pixel-info */
  int *totlev;     /* [tilting-level] each tilting-level has a different amount of intersections (beam<->wall) */
  int nztilt;      /* at how many levels tilting has to be performed = amount of "tilting-levels" */
  int tipawc;      /* = 1 in case there is a 3D-watercloud-file */
  int tipaic;      /* = 1 in case there is a 3D-icecloud-file */
  
  /* INPUT for tipa_xyzintersec (dir and dirdiff) STEPs 0.3, 1, 2, 3 */
  double szarad;
  double alpha;
  double dx;
  double dy;
  int maxzpixel;
  int minzpixel;
  int kxmax;
  int kxmin;
  int kxadd;
  int pymax;
  int pymin;
  int pyadd;
  int ireflevel;   /* index of level at which tilting has to be carried out */
  double reflevel; /* this is tipa->level[iz] for tipa dir and zout[0] in case of make_tipa (tipa dirdiff)*/
    
  /* OUTPUT of tipa_xyzintersec (dir and dirdiff) STEPs 0.3, 1, 2, 3 */
  int **zmatrix;
  int **xmatrix;
  int **ymatrix;
  double *zz;
  double *zx;
  double *zy;
/*   double ****taudircld;  taudircld[wavelength][tilting-level][iipa][jipa] */
  double **taudircld; /* taudircld[wavelength][tilting-level] NOTE eventually wavelength is not necessary! */
  
} tipa_struct; /* ulrike: output->wc.tipa or output->ic.tipa */

/* Caoth output STRUCTURE */
typedef struct {
  char *name;
  char *name1;
  char *name2;
  char *fullname;

  float *zd;          /* level altitudes                            */
  int   nlev;         /* number of levels                           */

  caoth_microphys_struct microphys;  /* microphysical properties      */
  optprop_struct         optprop;    /* optical properties            */
  ssprop_struct          ssprop;     /* single scattering properties  */

  crystal_prop_struct *raytracing_prop;   /* should be moved into microphys or optprop */
  int     n_raytracing_prop;              /* should be moved into microphys or optprop */
  
  /* the following two are used by MYSTIC only                      */
  int cldproperties; /* defined optical or microphysical properties */
  int nonHG;         /* Henyey-Greenstein or not                    */

  int newsiz;      /* flag for the Fortran function wcloud() */

  float *z3D;      /* checking of the 3D grid */ 
  int   Nz3D;
  
  tipa_struct tipa; /* ulrike */

  float cloudcover;
  
} caoth_out_struct;

/* Caothoff input STRUCTURE quickfix for new option names no_absorption and no_scattering*/

typedef struct {
  char   *name;
  int    no_absorption;
  int    no_scattering;
} caothoff_inp_struct;


/* Caoth input STRUCTURE */

typedef struct {
  char   *name;
  char   *name1;
  char   *name2;
  char   *fullname;

  int     source;
  int     layer;

  int     ipa;

  char   *filename;
  char   *properties_filename;

  float   cloudcover;

  float  **modify;

  float   reff; 

  int     no_scattering;

  int     properties;
  int     reff_prop;
  int     unscaled;
  int     fu2yang;
  int     *ic_fu;
  int     interpolate;
  int     habit;
  int     roughness;
} caoth_inp_struct;

/**********************/
/* Aerosol structures */
/**********************/

/* Aerosol input STRUCTURE */

typedef struct {
  int     standard;  /* Shettle standard aerosol defined    */
  int     spec;      /* specific aerosol properties defined */

  char **filename;   /* array of filenames */
  int haze;
  int vulcan;
  int seasn;
  float alpha;
  float beta;
  float alpha_0;
  float alpha_1;
  float alpha_2;

  float visibility;
  
  int   n_species;        /* number of aerosol species */
  char  **species_names;  /* aerosol species name [i_aer] */
  char  *mixture_name;  /* aerosol_species_file mixture_name for default aerosol mixtures */

  float re;          /* refractive index            */
  float im;          /* (independent of wavelength) */

  float tau_wvl_lambda;
  float tau_wvl_tau;


  float  **modify;

  int   no_scattering;
  int   profile_modtran;    /* squeeze aerosol profile as in MODTRAN */
} aer_inp_struct;

/* Aerosol output STRUCTURE */

typedef struct {
  float   visibility;

  float *zd;       /* level altitudes                  */
  int   nlev;      /* number of levels                 */

  float *zd_alt;       /* level altitudes, starting at 'altitude'   */
  int   nlev_alt;      /* number of levels, starting at 'altitude'  */

  /* optical properties profile */
  optprop_struct optprop;
  ssprop_struct ssprop; 
  optprop_struct *optprop_ipa; /* FIXMECE This should be removed. */
  optprop_struct *optprop_opac;

  int   n_species;        /* number of aerosol species */
  char  **species_names;  /* aerosol species name [i_aer] */
  float **massdens;       /* aerosol mass density profiles [i_aer][lev] */
  float **massdens_zout;  /* aerosol mass density at zout levels */
  char  *aerosol_library; /* library for aerosol optical properties */

  double   *xsiz;    /* size distribution              */
  double   *ysiz;
  int      nsiz;

  float   *re_r;     /* real and imaginary part of the */
  float   *im_r;     /* refractive index               */
} aer_out_struct;

/***************************/
/* Fluorescence structures */
/***************************/

/* Fluorescence input structure */

typedef struct {
  int    source;          /* where to get the flourescence from (constant, file, library) */
  double fluorescence;    /* constant input fluorescence value */
} flu_inp_struct;

/* Fluorescence output structure */

typedef struct {

  int    source;         /* where to get the fluorescence from (constant, file, library) */
  double fluorescence;   /* constant input fluorescence value */
  int    file_nlambda_fluorescence;  /* number of elements in the fluorescence file   */
  double *file_fluorescence;          /* fluorescence on the spectral grid of the file */
  double *file_lambda_fluorescence;   /* wavelength grid of the fluorescence file      */

  double *fluorescence_r;             /* fluorescence on the internal wavelength grid  */
  int    nlambda_fluorescence;

} flu_out_struct;

/*********************/
/* Albedo structures */
/*********************/

/* albedo input structure */

typedef struct {
  int    source;          /* where to get the albedo from (constant, file, library) */
  int    surface;         /* reference number for spectral albedo type inside a library */
  float  albedo;          /* constant input albedo value */
  char  *netCDF_alb_name; /* name of the albedo variable in the netCDF file */
  char  *netCDF_surf_name; /* name of the albedo variable in the netCDF file */
  char  *library;	  /* albedo library*/
  float  scale_factor;    /* in order to convert data in "per cent" or other ... */
  int    surf_type_map;   /* specifies, if the surface type is picked up in a map */
} alb_inp_struct;

/* albedo output structure */

typedef struct {

  int    source;         /* where to get the albedo from (constant, file, library) */
  int    surface;        /* reference number for spectral albedo type inside a library */
  float  albedo;         /* constant input albedo value */
  int    surf_type_map;  /* specifies, if the surface type is picked up in a map */

  int    file_nlambda_albedo;  /* number of elements in the albedo file   */
  float *file_albedo;          /* albedo on the spectral grid of the file */
  float *file_lambda_albedo;   /* wavelength grid of the albedo file      */

  float *albedo_r;             /* albedo on the actual internal wavelength grid  */
  int    nlambda_albedo;

} alb_out_struct;

/* RPV (Rahman, Pinty, Verstraete) structure */

typedef struct {
  float *rpv; /* sigma, t1, t2: Deguenther and Meerkoetter [2000]; t1: snow; includes forward peak; scale: optional scaling factor for the BRDF, to adjust albedo */
  int   source;
  char  *library; /* rpv_library*/
} rpv_inp_struct;


typedef struct {
  float *hapke;
  int source;
} hapke_inp_struct;

typedef struct {
  float *rossli;
  int source;
  int hotspot;
} rossli_inp_struct;


/* Cox and Munk structure */

typedef struct {
  char *pcl_netCDF_name;
  float pcl_scale_factor;
  float *param;
  int   solar_wind;
  char *sal_netCDF_name;
  float sal_scale_factor;
} cm_inp_struct;

/* Input structure for polarized reflectance function*/
/* So far only water BPDF by Tsang (Mishchenko) is included */
/* which has only wind speed as parameter */
typedef struct {
  int type;
  float u10;
} bpdf_inp_struct;


/* RossLi (Ambrals) structure */
/* to store wavelength dependent BRDF        */
typedef struct {
  int   n;

  float *lambda;
  float *iso;     /* input wavelength grid    */
  float *vol;
  float *geo;

  float *iso_r;   /* internal wavelength grid */
  float *vol_r;
  float *geo_r;
} rossli_out_struct;

/* Hapke structure */
/* to store wavelength dependent BRDF        */
typedef struct {
  int   n;

  float *lambda;
  float *w;     /* input wavelength grid    */
  float *b0;
  float *h;

  float *w_r;   /* internal wavelength grid */
  float *b0_r;
  float *h_r;
} hapke_out_struct;

/* RPV (Rahman, Pinty, Verstraete) structure */
/* to store wavelength dependent BRDF        */
typedef struct {
  int   n;

  float *lambda;
  float *rho0;     /* input wavelength grid    */
  float *k;
  float *theta;
  float *scale;
  float *sigma;
  float *t1;
  float *t2;

  float *rho0_r;   /* internal wavelength grid */
  float *k_r;
  float *theta_r;
  float *scale_r;
  float *sigma_r;
  float *t1_r;
  float *t2_r;

} rpv_out_struct;

/* Array of RPV (Rahman, Pinty, Verstraete) structures */
/* for different surfaces, addressed by a label        */
typedef struct {
  int   n;

  char  **label;
  char  **filename;

  float **albedo_r;   /* albedo on internal wavelength grid */
  rpv_out_struct *rpv;
} surfaces_out_struct;


typedef struct {
  float tautot;
  int lctop;
  int type;
  float ref;
} sssi_out_struct;


/*************************/
/* Wavelength structures */
/*************************/

/* Wavelength input grid STRUCTURE */

typedef struct {
  float start;
  float end;

  int start_index;
  int end_index;
} wl_inp_struct;

/* Wavelength output grid STRUCTURE */

typedef struct {
  float start;
  float end;

  int   type;

  int   nlambda_h;
  int   nlambda_s;
  int   nlambda_t;
  int   nlambda_r;

  int   nlambda_h_print3D;

  float *lambda_h;
  float *lambda_s;
  float *lambda_t;
  float *lambda_r;

  int    use_reptran;             /* Indicates whether representative wavelengths are used. */

  int    *reptran_band_t;         /* Assignes lambda_t to representative wavelengths bands. */

  int    *iwvl_in_reptran_file_r; /* Indices of lambda_r in reptran_file. */

  float  *extra_reptran_r;        /* Extraterrestrial flux on lambda_r grid. */

  int    *nlambda_in_reptran_band;/* Number of representative wavelengths per band. */
  int    **reptran_band;          /* Representative wavelengths for each band (refers to index in lambda_r grid). */
  double **weight_reptran_band;   /* Weights of representative wavelengths. */
  double *width_of_reptran_band;  /* Width of band (in nm for solar source and cm^-1 for thermal source). */

  int   *map_e2h;

  int   ignore_solar_file;

  float *fbeam;
  float *filter;

  float *wvnmlo_r; 
  float *wvnmhi_r; 
  
  float *wvnmlo_t; 
  float *wvnmhi_t; 

  float *wvnmlo_h;
  float *wvnmhi_h;

  float delta_wvl_raman; /* The +- wavelengths to be considered for Raman shifts. */
  float delta_wvl_raman_lower; /* The - wavelengths to be considered for Raman shifts. */
  float delta_wvl_raman_upper; /* The + wavelengths to be considered for Raman shifts. */
  float delta_wvl_raman_extra; /* Add a little extra to lower and upper to be on the safe side. */
  int   raman_start_id;
  int   raman_end_id;

  int  nlambda_rte_lower;
  int  nlambda_rte_upper;

  float *pfraction;
} wl_out_struct;


/**********************/
/* Altitude structure */
/**********************/

typedef struct {
  int    source;          /* where to get the altitude from (direct input or map)                        */
  float  altitude;        /* surface altitude over sea level                                             */
  float  scale_factor;    /* scale factor, e.g. 9.8 in order to convert from geopotential meter to meter */
  char  *netCDF_alt_name; /* name of the altitude variable in the netCDF file                            */
  float  altitude_dz;     /* fakes an increased vertical resolution, profiles are interpolated           */
} alt_inp_struct;

typedef struct {
  int    source;        /* where to get the altitude from (direct input, map, or atmosphere) */  
  float  altitude;      /* Same as altitude from input structure is set.
			   Otherwise obtained from atmosphere input file */
} alt_out_struct;

/*************************/
/* Atmosphere structures */
/*************************/

/* Atmosphere input data STRUCTURE */

typedef struct {
  char **filename;      /* Filename of density profile, input using dens_file  */
  int    rs_source;     /* where the radiosonde data comes from (file or ECMWF file) */
  int    n_rs_gas;      /* number of gases, which are given by radiosonde file */
  int   *rs_gas;        /* MOL-number for gases given by radiosonde data */
  int   *rs_unit;       /* specifies unit, in which radiosonde data is given in */  
  int    rs_add_upper_levels; /* marker in order to show, if levels above rs-levels are added or not */
  float *column;        /* Total column of a profile, input using dens_column  */
/*   float *scale_factor; /\* Factor to scale profile with, calculated           *\/ */
  int   *unit_column;   /* Unit of the column value, input using dens_column   */
  int   *unit_profile;  /* Unit of the profile, input using dens_file          */
  float *mol_mass;      /* molecular mass of one molucule in u                 */

  float *z_atm_forced_sea;  /* forced z_grid of the atmosphere relative to the sea surface */
  int    nz_atm_forced_sea; /* number of forced atmosphere levels                          */
  float *zout_sur;      /* output altitudes, relative to the surface           */
  float *zout_sea;      /* output altitudes, relative to sea level             */
  float *press_zout;    /* output levels in pressure  hPa                      */
  int    nzout;         /* number of output altitudes                          */
  int    zout_source;   /* input switch */

  float  sza;
  int    sza_source;    /* where to get sza from (direct input, file, or lat/lon/time)*/
  float  phi0;
  float  sza_spher;
  float  phi0_spher;
  int    zout_interpolate; /* switch for zout-level interpolation */
  int    time_interpolate; /* switch for time interpolation (of netCDF data input) */ 
  int    z_interpolate;    /* group specification of the interpolation methods     */
                           /* of the atmospheric profiles                          */
  int    interpol_method_press;  /* individual z_interpolation method  */
  int    interpol_method_temper; /* individual z_interpolation method  */
  int    interpol_method_refind; /* individual z_interpolation method  */
  int   *interpol_method_gas;    /* individual z_interpolation methods */
  int   *well_mixed_gas;         /* specifies if a gas is well mixed   */

  int    ECMWF_new_format;
  int    ECMWF_ozone_climatology;

} atm_inp_struct;

/* Atmosphere output data STRUCTURE */

///3DA
typedef struct {

  float ***press;
  float **refind;        /* wavelength dependent refractive index refind[iv][lc] */
  float ***temper;
  float ***temper_avg;     /* average temperature temper_avg[lc] in the layer from z[lc] to z[lc+1]  */ 

  float ****dens;          /* density profiles */
  float ****dens_avg;      /* layer average density dens_avg[lc][ix][iy] in the layer from z[lc] to z[lc+1]  */
  float ***denstab;      /* density profiles as function of altitude and sza */
  float ***denstab_avg;  /* layer average density profiles as function of altitude and sza */
  float **denstab_amf;   /* density profile as function of altitude and sza, for sdisort */
  float **dens_zout;     /* gas number densities  on the zout - grid [MOL_GAS][lc] */
  float *press_zout;     /* pressure              on the zout - grid */
  float *temper_zout;    /* temperature           on the zout - grid */
  float *temper_d_zout;  /* dewpoint temperature  on the zout - grid */
  float *dtheta_dx;     /* derivative of temperature with respect to x, this only works with atm_maps */
  float *dtheta_dy;     /* derivative of temperature with respect to y, this only works with atm_maps */
  float *dtheta_dz;      /* derivative of potential temperature with respect to z */
  float *c_p;            /* specific heat capacity at constant pressure of moist air */
  float *theta;          /* potential temperature at atmospheric grid */
  float *theta_zout;     /* potential temperature on the zout - grid  */
  float *theta_e_zout;   /* equivalent potential temperature on the zout - grid */

  float *sza_denstab;    /* Ary for the sza in denstab */
  
  float z_cpt_sea;       /* height of the cold point tropopause above sea level */
  float temper_cpt;      /* temperature at the cold point tropopause */
  float press_cpt;       /* pressure at the cold point tropopause */

  int   denstab_id;      /* Density id for denstab, MOL_* */
  int   nsza_denstab;    /* Number of sza in denstab */

} atm_microphys_struct;


//3DAbs
typedef struct {
  
  double *****tau_rayleigh_r;    /* wavelength-dependent Rayleigh scattering optical depth   */
  double *****tau_molabs_r;      /* wavelength-dependent molecular absorption optical depth  */
  double *****tau_molabs_md_r;   /* wavelength-dependent molecular absorption optical depth  */
			       /* with density-matrix absorption subtracted. For AMF       */
			       /* calculations                                             */
                               /* format: tau_molabs_r[lc][iv][iq]                         */
                               /* ??? this is not consistent with tau_molabs_user_r[][]    */
                               /* ??? and should be changed accordingly                    */

  /* The following are explicitely defined by the user and do not depend on subband iq     */
  /* they are copied to the above structures which are used as input to the solvers        */
  /* ATTENTION: These structures are not redistributed and must therefore not be used      */
  /* after the call to setup_redistribute(); optical_properties relies only on the above   */
  /* quantities.                                                                           */

  double  *tau_rayleigh_user;   /* monochromatic Rayleigh scattering optical depth         */
  double **tau_rayleigh_user_r; /* wavelength-dependend Rayleigh scattering optical depth  */
                                /* format: tau_rayleigh_user_r[iv][lc]                     */
  double  *tau_molabs_user;     /* monochromatic molecular absorption optical depth        */
  double **tau_molabs_user_r;   /* wavelength-dependend molecular absorption optical depth */
                                /* format: tau_molabs_user_r[iv][lc]                       */
} atm_optprop_struct;


typedef struct {
  char* pol_scat_files;
  float ***cdento;                /* average air density over a layer assuming  dens = exp (-c * z) */

  float *temper_out;

  atm_microphys_struct microphys;
  atm_optprop_struct   optprop;

  int   molabs;                 /* molecular absorption flag         */
  int   rayleigh;               /* Rayleigh scattering flag          */

  float *zd;
  int   nlev;
  int   nlyr;

  float *zout_sur;     /* output altitudes, relative to the surface    */
  float *zout_sur_org; /* first zout grid, which is in each case there */
  float *zout_sea;     /* output altitudes, relative to sea level      */
  float *zout_sea_org; /* first zout grid, which is in each case there */
  float *press_zout_org;    /* first pressure_out grid */
  int    nzout;        /* number of output altitudes                   */
  int    nzout_user;   /* number of output altitudes, same as nzout,   */
                       /* except when internally more zout are needed, */
                       /*eg. for raman scattering                      */
  int   *zout_user_index;/* indecies for the user output altitudes     */
  int   *zout_comp_index;/* indecies for the internal computational altitudes     */
  int    nzout_org;    /* number of output altitudes, of org grids     */
  int   *zout_index;   /* relationship between the zout indecies       */

  float *zd_common;
  int    nlev_common;
  int    nlyr_common;

  float *z_model_level;        /* altitudes from model data (ECMWF or ECHAM) */
  int   nz_model_level;       /* number of altitudes levels of the model data (ECMWF or ECHAM) */

  float *z_model_layer;        /* altitudes from model data (ECMWF or ECHAM) */
  int   nz_model_layer;       /* number of altitudes levels of the model data (ECMWF or ECHAM) */

  double **wght_r;        /* ck quadrature weight                   */
  double ***ipa_wght_r;   /* independent pixel ck quadrature weight */
  int    *nq_r;           /* number of ck quadrature weights        */

  float *file_phi0;
  float *file_sza;
  float *file_lambda_sza;
  float *sza_r;
  float *sza_t;
  float *phi0_r;

  int    file_nlambda_sza;
  int    nlambda_sza;

  int    nmom;  

  int nthreed;   /* number of 3D layers   */
  int *threed;   /* flag if a layer is 3D */
  
  /* 3D grids must correspond to caoths included with wc/ic/profile_file */
  int Nxcld;     /* 3D horizontal caoth grid */
  int Nycld;     /* 3D horizontal caoth grid */
  int Nzcld;     /* 3D horizontal caoth grid */
  double dxcld;
  double dycld;

  float *z_3d_layer;        /* altitudes from 3d atmosphere data*/
  
  /* 3D grids for 3D atmosphere file */
  int Nxatm;
  int Nyatm;
  int Nzatm; 
  double dxatm;
  double dyatm;
  
} atm_out_struct;

/*******************************************************************/
/* The ck structure holds all the information required for         */
/* a correlated-k approximation.                                   */
/* Components:                                                     */
/*    scheme               CK_KATO, CK_FU                          */
/*    wvlc  [1 ... n_wvl]  central wavelength [nm]                 */                   
/*    wvnlo [1 ... n_wvl]  lower wavenumber [cm-1] for each band   */    
/*    wvnhi [1 ... n_wvl]  upper wavenumber [cm-1] for each band   */    
/*    n_wvl                number of wavelength bands              */
/*    comp  [1 ... n_comp][1 ... n_wvl]                            */
/*                         monochromatic cross sections bands      */
/*                         which don't require c-k but only a      */
/*                         monochromatic calculation               */
/*                         (e.g. ozone in the UV/VIS)              */
/*    wght  [1 ... n_comp][1 ... n_wvl]                            */
/*                         number of quadrature points             */
/*******************************************************************/

typedef struct {   /* correlated-k structure */
  int scheme;

  double *wvlc;
  double *wvnlo;
  double *wvnhi;
  int n_wvl;

  double **comp;
  int    **wght;

  void *x_thermal;  /* photon fraction; x_thermal[iv][iq1][iq2][...] */
  void *x_solar;    /* photon fraction; x_solar  [iv][iq1][iq2][...] */

  int n_comp;
} ck_struct;

/* Cross section output data STRUCTURE */

typedef struct {
  char  **filename;  /*  Filename of cross section versus wavelength, input using crs_file  */
  float *****crs;      /* 5D ary, for cross section[ix, iy, temperature (or layer), type, wavelength] */
  float **crs_amf;   /* 2D ary for amf calculations, needed by sdisort */
  int number_of_ramanshifts;      /* Total number of Raman transitions considered */
  int number_of_ramanwavelengths; /* Total number of Raman transitions considered plus the wanted wavelength*/
  double *wvl_of_ramanshifts; /* Wavelengths of Raman transitions considered for each user wavelength*/
  double **crs_raman_RL;   /* 2D ary for  Raman scattering cross section [layer][Raman shifted wavelengths]              */
                          /* weighted by N2 and O2 volume mixing ratio. Raman loss terms, shifted from wanted wavelength*/  
  double **crs_raman_RG;   /* 2D ary for  Raman scattering cross section [layer][Raman shifted wavelengths]              */
                          /* weighted by N2 and O2 volume mixing ratio. Raman gain terms, shifted to wanted wavelength*/  
} crs_out_struct;

typedef struct {
  float** kabs; /* 2D ary, for absorption coefficient[temperature (or layer), wavelength] */
} kabs_out_struct;

typedef struct {
  double ****crs;      /* cross section [ix, iy, layer, quadrature point] */
  double *weight;    /* weight [quadrature point]               */
  int ngauss;        /* number of quadrature points             */
  int nlev;
} ck_profile;

typedef struct {
  ck_profile **profile;      /* 2D ary, for ck_profile [type, band] */
  ck_profile ***ipa_profile; /* same as profile for different independent columns [ipa][type][band] */
} crs_ck_out_struct;


typedef struct {
  float *cf;
  float cf_total;
  float *cf_zout;
  float *zd;
  int nlev;
} cf_out_struct;


/* Wind STRUCTURE */

typedef struct {
  /* quantities as defined by the user */
  float *u;   /* z-profile of the wind component in m/s */
  float *v;   /* z-profile of the wind component in m/s */
  float *w;   /* z-profile of the wind component in m/s */
  float *z;   /* z-levels in km */
  size_t nlev; /* number of levels from the arrays above */

  float *u_zout;   /* wind component in m/s at zout levels */
  float *v_zout;   /* wind component in m/s at zout levels */
  float *w_zout;   /* wind component in m/s at zout levels */
} wind_out_struct;

			    
/* INPUT STRUCTURE */

typedef struct {
  int quiet;
  int verbose;
  int test_optical_properties; 

  int          spectrum_unit; /* kind of units: W/(m**2 nm) or W/(m**2 cm-1) or W/(m**2 band) */

  char         **filename;    /* array of filenames */

  char	*atmosphere_filename; /*Can be either whole path of atmosphere file or 'subarctic_winter'/'_summer', 'midlatitude_winter'/'_summer', 'tropics', 'US_standard' */
  char         *atmosphere3d_filename;/* path of 3D atmosphere_ file */
  int          atmosphere3d;
  
  float        r_earth;     

  int          cloud_overlap; 
  
  int          ipa;
  int          ipa3d;         /* independent pixel for 3D-clouds (ulrike) */
  int          tipa;          /* tilted independent pixel/column approximation (ulrike) */

  int          filter_function_normalize;

  int          convolve;
  int          spline;
  int          header;
  int          output_format;
  int          source;
  double       bandwidth;
  int          bandwidth_unit;
  int          raman;                  /* Flag to turn Raman option on                          */
  int          raman_original;         /* Flag to turn original slow Raman option on            */
  int          raman_fast;             /* Flag to turn fast Raman option on, default on         */
  int          n_raman_transitions_N2; /* Number of transitions considered for N2               */
  int          n_raman_transitions_O2; /* Number of transitions considered for O2               */
  float        delta_wvl_raman;        /* The +- wavelengths to be considered for Raman shifts. */

  int	      *crs_model; /* Old: o3_spectrum, no2_spectrum, rayleigh_crs */

  int          rayleigh;
  float        rayleigh_depol;
  int          molabs;
  int          absorption;
  int          scattering;
  float        spline_lambda_0;
  float        spline_lambda_1;
  float        spline_lambda_step;
  struct tm    UTC;  /* time as Universal Time Coordinated*/
  struct tm    LAT;  /* time as Local Apparent Time */
  int          time_zone;
  float        delta_time_start;   /* in seconds, for an sza averaged over an time interval */
  float        delta_time_end;     /* in seconds, for an sza averaged over an time interval */
  float        latitude;     
  float        lat_degrees;     
  float        lat_minutes;     
  float        lat_seconds;     
  int 	       lat_signum;     
  float        longitude;   
  float        lon_degrees;     
  float        lon_minutes;     
  float        lon_seconds;     
  int 	       lon_signum;     
  int          sat_pixel_x;
  int          sat_pixel_y;
  int          co2_spectrum;
  int          o2_spectrum;
  float        no2_column;
  int          so2_spectrum;
  int          bro_spectrum;
  int          oclo_spectrum;
  int          hcho_spectrum;
  int          o4_spectrum;
  /* float        precip_water; */
  float        pressure;
  float        surface_temperature;
  char        *netCDF_name_surf_T;
  float       *mixing_ratio;
  int          calibration;
  int          optimize_fortran;
  int          optimize_delta;

  int write_ext_to_file; /* needed by ARLEM. This is a quick fix!!! This should be somewhere else BCA */
  int write_output_as_netcdf; /* option whether to write lidar output into ASCII or netcdf */

  int          disort2_brdf;   /* BRDF type for disort2 - 2D BRDFs for MYSTIC handled separately */
  /* float  u10; */  /* wind speed for Cox and Munk parameterization, probably not used */

  rte_struct     rte;
  atm_inp_struct atm;
  alt_inp_struct alt;
  alb_inp_struct alb;
  aer_inp_struct aer;
  flu_inp_struct flu;

  rpv_inp_struct rpv;
  hapke_inp_struct hapke;
  rossli_inp_struct rossli;
  cm_inp_struct cm;
  bpdf_inp_struct bpdf; 

  wl_inp_struct  wl;

  int            n_caoth;
  int            n_caothoff;
  caoth_inp_struct *caoth;
  caothoff_inp_struct *caothoff;
  int            i_wc;
  int            i_ic;
  int            ck_scheme;
  char          *ck_scheme_filename; /*initializing ck_scheme = CK_GENERIC*/
  int            ck_h2ocont;

  char          *ck_reptran_arg;   /*input argument for initializing ck_reptran_channel or ck_reptran_option*/
  int            ck_reptran_option;
  char          *ck_reptran_channel;

  int		*ck_abs;

  char		*output_unit_processing	;/*Input argument to initialize either processing or output_unit */
  int            processing;
  int            output_unit;
  int            heating;

  int    output_user_flag;
  int    n_output_user;
  int   *output_user;     /* output_user    [column of output] = kind of output */
  int   *output_user_gas; /* output_user_gas[column of output] = MOL_GAS        */

  double *sslidar;
  int sslidar_nranges;
  int sslidar_polarisation;
} input_struct;


/* OUTPUT STRUCTURE */

/* elements whose memory requirements have been free'ed are listed behind */
/* each declaration. Once they are all done rm this comment.              */

typedef struct {
  float       eccentricity_factor;
  float       sunshine_fraction;     /* part of the time interval, where sun is above horizon */
  int         spectrum_unit;         /* kind of units: W/(m**2 nm) or W/(m**2 cm-1) or W/(m**2 band) */
  int         bandwidth_unit;        /* kind of units: nm or cm-1 */

  float       *column;               /* Total column of a profile, input using dens_column */
  float       *scale_factor;         /* Factor to scale profile with, calculated           */
  float      **column_denstab;       /* Total column of a profile, input using dens_column [gas][sza] */
  float      **scale_factor_denstab; /* Factor to scale profile with, calculated           [gas][sza] */

  float       precip_water;
  float       water_scale_factor;

  float       *mixing_ratio;

  float       *rayleigh_depol;  /* wavelength-dependent depolarization */

  long int    mc_photons;
  double      **mc_photons_r;   /* photon fraction for each wavelength and quadrature point */

  float       *tausol;          /* =dtau_wc+dtau_ic+dtau_mol+dtau_aer; see ancillary.c (ulrike), for tipa dir */
  float       *dtauc;           /* free'ed */
  float       *dtauc_clr;       /* = dtauc - dtau[i_wc] = optical thickness excluding wc; for twomaxrnd, BM 13.10.2016 */
  float       *dtauc_md;        /* free'ed */
  float       *ssalb;
  float       *ssalbR;
  float       *ssalb_clr;       /* = ssalb excluding wc; for twomaxrnd, BM 13.10.2016 */
  float       *ssalbRL;
  float       *utau;
  int         nphamat;
  float       ***pmom;
  float       *pmom01_clr;      /* first moment (= asymmetry parameter) excluding wc; for twomaxrnd, BM 13.10.2016 */
  int         **ntheta;
  float       ***theta;
  double      ***mu;
  float       ***phase;

  float       surface_temperature;

  /* flux, actinic flux, and intensities on the output grid */ 
  float       **rfldir;     /* direct flux                    */
  float       **rfldn;      /* diffuse downward flux          */
  float       **flup;       /* diffuse upward flux            */ 

  float       **uavgso;     /* direct actinic flux            */
  float       **uavgdn;	    /* diffuse downward actinic flux  */
  float       **uavgup;	    /* diffuse upward actinic flux    */

  float       **uavg;       /* average intensity              */
  float       ***u0u;       /* azimutally averaged intensites */
  float       ****uu;       /* intensities                    */
  float       **heat;       /* heating rates,                   f(zout,lambda) */
  float       **emis;       /* thermal emitted energy in K/day, f(zout,lambda) */
  float       **w_zout;     /* vertical wind caused by radiative heating on the zout - grid */

  float       **sslidar_nphot;   /* nphoton count for sslidar */
  float       **sslidar_nphot_q; /* Q part of stokes of nphot */
  float       **sslidar_ratio;   /* lidar ratio for sslidar   */

  float       **albmed;     /* albedo of complete atmosphere */
  float       **trnmed;     /* transmittance of complete atmosphere */

  /* flux, actinic flux, and intensities on the internal wavelength grid */ 
  float       **rfldir_r;  /* direct flux                    */
  float       **rfldn_r;   /* diffuse downward flux          */
  float       **flup_r;	   /* diffuse upward flux            */
			                                       
  float       **uavgso_r;  /* direct actinic flux            */
  float       **uavgdn_r;  /* diffuse downward actinic flux  */
  float       **uavgup_r;  /* diffuse upward actinic flux    */
			                                       
  float       **uavg_r;	   /* average intensity              */
  float       ***u0u_r;	   /* azimutally averaged intensites */
  float       ****uu_r;	   /* intensities                    */
  float       **heat_r;    /* heating rates;                             depends on (zout,lambda) */
  float       **emis_r;    /* thermal emitted energy in K/day;           depends on (zout,lambda) */
  float       **w_zout_r;  /* vertical wind caused by radiative heating; depends on (zout,lambda) */

  float       **sslidar_nphot_r;   /* nphoton count for sslidar */
  float       **sslidar_nphot_q_r; /* Q part of stokes of nphot */
  float       **sslidar_ratio_r;   /* lidar ratio for sslidar   */

  float       **albmed_r;     /* albedo of complete atmosphere */
  float       **trnmed_r;     /* transmittance of complete atmosphere */
  
  /* flux, actinic flux, and intensities on the transmittance grid */ 
  float       **rfldir_t;  /* direct flux                    */
  float       **rfldn_t;   /* diffuse downward flux          */
  float       **flup_t;	   /* diffuse upward flux            */
			                                       
  float       **uavgso_t;  /* direct actinic flux            */
  float       **uavgdn_t;  /* diffuse downward actinic flux  */
  float       **uavgup_t;  /* diffuse upward actinic flux    */
			                                       
  float       **uavg_t;	   /* average intensity              */
  float       ***u0u_t;	   /* azimutally averaged intensites */
  float       ****uu_t;	   /* intensities                    */
  float       **heat_t;    /* heating rates;                             depends on (zout,lambda) */
  float       **emis_t;    /* thermal emitted energy in K/day;           depends on (zout,lambda) */
  float       **w_zout_t;  /* vertical wind caused by radiative heating; depends on (zout,lambda) */

  float       **sslidar_nphot_t;   /* nphoton count for sslidar */
  float       **sslidar_nphot_q_t; /* Q part of stokes of nphot */
  float       **sslidar_ratio_t;   /* lidar ratio for sslidar   */

  float       **albmed_t;     /* albedo of complete atmosphere */
  float       **trnmed_t;     /* transmittance of complete atmosphere */

  /* wavelength integrated flux, actinic flux, and intensities on the output grid */
  double      *rfldir_int;  /* direct flux                    */
  double      *rfldn_int;   /* diffuse downward flux          */
  double      *flup_int;    /* diffuse upward flux            */ 

  double      *uavgso_int;  /* direct actinic flux            */
  double      *uavgdn_int;  /* diffuse downward actinic flux  */
  double      *uavgup_int;  /* diffuse upward actinic flux    */

  double      *uavg_int;    /* average intensity              */
  double      **u0u_int;    /* azimutally averaged intensites */
  double      ***uu_int;    /* intensities                    */  

  double      *heat_int;    /* heating rates f(z) */
  double      *emis_int;    /* heating rates f(z) */
  double      *w_zout_int;  /* vertical wind caused by radiative heating; f(z) */

  float       **sslidar_nphot_int;   /* nphoton count for sslidar */
  float       **sslidar_nphot_q_int; /* Q part of stokes of nphot */
  float       **sslidar_ratio_int;   /* lidar ratio for sslidar   */

  double       *albmed_int;     /* albedo of complete atmosphere */
  double       *trnmed_int;     /* transmittance of complete atmosphere */

  /* integrated properties */
  double      incident;


  /* MYSTIC 3D fields */
  int islower;  /* to save memory, only those pixels are allocated */
  int isupper;  /* which are actually sampled; in partcular for    */
  int jslower;  /* backward Monte Carlo only a small fraction of   */
  int jsupper;  /* the pixels is usually sampled .                 */

  /* MYSTIC 3D fields on the output grid */  
  float       ****rfldir3d; 
  float       ****rfldn3d;
  float       ****flup3d;
  float       ****uavgso3d; 
  float       ****uavgdn3d; 
  float       ****uavgup3d; 
  float       ******radiance3d;
  float       ****abs3d;
  float       ****absback3d;

  /* corresponding variances */
  float       ****rfldir3d_var;
  float       ****rfldn3d_var;
  float       ****flup3d_var;
  float       ****uavgso3d_var; 
  float       ****uavgdn3d_var; 
  float       ****uavgup3d_var; 
  float       *****radiance3d_var;
  float       ****abs3d_var;
  float       ****absback3d_var;

  /* MYSTIC 3D fields on the internal wavelength grid */  
  float       ****rfldir3d_r;
  float       ****rfldn3d_r;
  float       ****flup3d_r;
  float       ****uavgso3d_r;
  float       ****uavgdn3d_r;
  float       ****uavgup3d_r;
  float       ******radiance3d_r;
  float       ****abs3d_r;
  float       ****absback3d_r;

  /* corresponding variances */
  float       ****rfldir3d_var_r; 
  float       ****rfldn3d_var_r;
  float       ****flup3d_var_r;
  float       ****uavgso3d_var_r;
  float       ****uavgdn3d_var_r;
  float       ****uavgup3d_var_r;
  float       *****radiance3d_var_r;
  float       ****abs3d_var_r;
  float       ****absback3d_var_r;

  /* MYSTIC 3D fields on the transmittance grid */  
  float       ****rfldir3d_t;
  float       ****rfldn3d_t;
  float       ****flup3d_t;
  float       ****uavgso3d_t;
  float       ****uavgdn3d_t;
  float       ****uavgup3d_t;
  float       ******radiance3d_t;
  float       ****abs3d_t;
  float       ****absback3d_t;

  /* corresponding variances */
  float       ****rfldir3d_var_t; 
  float       ****rfldn3d_var_t;
  float       ****flup3d_var_t;
  float       ****uavgso3d_var_t;
  float       ****uavgdn3d_var_t;
  float       ****uavgup3d_var_t;
  float       *****radiance3d_var_t;
  float       ****abs3d_var_t;
  float       ****absback3d_var_t;

  /* PolRadtran flux and intensities on the output grid */ 
  float       ***down_flux;     /* free'ed */
  float       ***up_flux;       /* free'ed */
  float       *****down_rad;    /* free'ed */
  float       *****up_rad;      /* free'ed */

  /* wavelength integrated PolRadtran flux and intensities on the output grid */ 
  double      **down_flux_int;     /* free'ed */
  double      **up_flux_int;       /* free'ed */
  double      ****down_rad_int;    /* free'ed */
  double      ****up_rad_int;      /* free'ed */

  /* PolRadtran flux and intensities on the internal wavelength grid */ 
  float       ***down_flux_r;   /* free'ed */
  float       ***up_flux_r;     /* free'ed */
  float       *****down_rad_r;  /* free'ed */
  float       *****up_rad_r;    /* free'ed */

  /* PolRadtran flux and intensities on the transmittance grid */ 
  float       ***down_flux_t;   /* free'ed */
  float       ***up_flux_t;     /* free'ed */
  float       *****down_rad_t;  /* free'ed */
  float       *****up_rad_t;    /* free'ed */

  float       *mu_values;       /* free'ed */
  float       *sza_h;           /* free'ed */
  float       *sza_s;           /* free'ed */

  atm_out_struct      atm;
  alt_out_struct      alt;
  alb_out_struct      alb;
  flu_out_struct      flu;
  rossli_out_struct   rossli;
  hapke_out_struct    hapke;
  rpv_out_struct      rpv;
  surfaces_out_struct rpv_surf;
  wind_out_struct     wind;

  crs_out_struct      crs;
  crs_ck_out_struct   crs_ck;
  kabs_out_struct     kabs;

  aer_out_struct      aer;

  caoth_out_struct *caoth;
  caoth_out_struct **caoth_ipa;
  caoth3d_out_struct *caoth3d;
  int Nx_3d;
  int Ny_3d;
  int delX_3d;
  int delY_3d;

  cf_out_struct     cf;
  cf_out_struct     *cfipa;

  mc_out_struct     mc;
  sssi_out_struct   sssi;

  float             *ipaweight;
  int               nipa;
  int               niipa;
  int               njipa;

  int               print_phi;
  int               molecular3d; 
  
  wl_out_struct     wl;
  ck_struct         ck;
} output_struct;

int uvspec(input_struct input, output_struct *output);
int uvspec_check(input_struct input);
int caoth_check(input_struct input);

int caoth_set(input_struct *input, output_struct *output);

int get_caothoff_index ( caothoff_inp_struct **caothoff,
		      int               *n_caothoff,
		      char              *name );
int get_caoth_index ( caoth_inp_struct **caoth,
		      int               *n_caoth,
		      char              *name,
		      int                old_input );

void pmesg (char *message, int verbose);


/* include all uvspec header files */
#include "albedo.h"
#include "aerosol.h"
#include "ancillary.h"
#include "atmosphere.h"
#include "ck.h"
#include "cloud.h"
#include "cloud3d.h"
#include "extraterrestrial.h"
#include "fluorescence.h"
#include "ipa.h"
#include "molecular.h"
#include "molecular3d.h"
#include "redistribute.h"
#include "solve_rte.h"
#include "sza.h"
#include "yang56.h"
#include "specrend_uvspec.h"

#endif
