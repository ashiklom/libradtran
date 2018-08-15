/************************************************************************
 * $Id: mystic.h 3311 2017-12-07 16:19:50Z bernhard.mayer $
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

#ifndef __mystic_h
#define __mystic_h

#if defined (__cplusplus)
extern "C" {
#endif

/* for debugging. This mode writes the randomstatus file, with */
/* help of it one can start shortly before a buggy photon      */
/* instead of having to simulate many good photons before      */
/* reaching the problematic one                                */
#undef WRITERANDOMSTATUS

#undef NEW_REFLECT 

#undef MUCHOUT
#define MUCHOUTPHOTON 99
#define CLDPRP

#include <stdio.h>  /* required for FILENAME_MAX */
#include "phasetable.h"
#include "cdisort.h"

#define NOT_A_NUMBER              0.0/0.0

/* radiation source */
#define MCSRC_NONE                0
#define MCSRC_SOLAR               1
#define MCSRC_THERMAL_SURFACE     2
#define MCSRC_THERMAL_ATMOSPHERE  3
#define MCSRC_THERMAL_BACKWARD    4/*TZ bt*/
#define MCSRC_LIDAR               5
#define MCSRC_BLITZ               6

/* albedo method */
#define MCALB_NONE                0
#define MCALB_LAM                 1
#define MCALB_LAM2D               2
#define MCALB_LAM2D_SPECTRAL      3
#define MCALB_RPV                 4
#define MCALB_RPV2D_SPECTRAL      5
#define MCALB_COXANDMUNK          6
#define MCALB_HAPKE               7
#define MCALB_ROSSLI              8
#define MCALB_ROSSLI2D            9
#define MCALB_TSANG               10

#define MCPAN_ALIGNMENT_NONE 0
#define MCPAN_ALIGNMENT_SUN  1
#define MCPAN_ALIGNMENT_MU   2

/* boundary conditions */
#define MCBCOND_NONE              0
#define MCBCOND_PERIODIC          1
#define MCBCOND_MIRROR            2
#define MCBCOND_ABSORB            3

/* calculate 3D absorption field? */
#define MCFORWARD_ABS_NONE        0    /* no absorption                     */
#define MCFORWARD_ABS_ABSORPTION  1    /* calculate absorption              */
#define MCFORWARD_ABS_EMISSION    2    /* calculate emission                */
#define MCFORWARD_ABS_HEATING     3    /* calculate (absorption + emission) */
#define MCFORWARD_ABS_ACTINIC     4    /* calculate actinic flux            */

#define MCABS_UNIT_W_PER_M2_AND_DZ  0  /* default absorbed energy per box */
#define MCABS_UNIT_W_PER_M3         1  /* normalised with height          */
#define MCABS_UNIT_K_PER_DAY        2  /* converted to K / day            */

#define MCBACKWARD_NONE           0 /* forward calculation           */
#define MCBACKWARD_RADIANCE       1 /* backward radiance             */
#define MCBACKWARD_EDIR           2 /* direct irradiance             */
#define MCBACKWARD_EDN            3 /* downward diffuse irradiance   */
#define MCBACKWARD_EUP            4 /* upward diffuse irradiance     */ 
#define MCBACKWARD_FDIR           5 /* direct irradiance             */
#define MCBACKWARD_FDN            6 /* downward diffuse irradiance   */
#define MCBACKWARD_FUP            7 /* upward diffuse irradiance     */ 
#define MCBACKWARD_ABS            8 /* absorption                    */
#define MCBACKWARD_ACT            9 /* actinic flux                  */
#define MCBACKWARD_EMIS           10 /* emission                      */
#define MCBACKWARD_HEAT           11 /* heating rate                  */
#define MCBACKWARD_CLDPRP         12 /* cloud properties              */
#define MCBACKWARD_EGLOB          13 /* direct + downward diffuse irradiance */
#define MCBACKWARD_EXP            14 /* diffuse irradiance in positive x direction  */
#define MCBACKWARD_EXN            15 /* diffuse irradiance in negative x direction  */
#define MCBACKWARD_EYP            16 /* diffuse irradiance in positive y direction  */
#define MCBACKWARD_EYN            17 /* diffuse irradiance in negative y direction  */


#define MCBACKWARD_HEAT_EMABS    1
#define MCBACKWARD_HEAT_EMABSOPT 2
#define MCBACKWARD_HEAT_DENET    3
#define MCBACKWARD_HEAT_HYBRID   4

/* flag what is defined in the wcloud3D.dat, 1,2,3 */
#define CLD_OPTPROP               1 /* explicit optical properties: ext, asy, ssa */
#define CLD_EXTREFF               2 /* extinction + effective radius              */
#define CLD_LWCREFF               3 /* liquid water content + effective radius    */


/* azimuth conventions */
#define MCAZIMUTH_NEW             0
#define MCAZIMUTH_OLD             1

/* small number for floating point comparisons */
#define MC_EPSILON                1E-6

/* scatter types */
#define MCCAOTH_NONE              -1
#define MCCAOTH_TOT               0
#define MCCAOTH_MOL               1   /* this has to be CAOTH_MOL+1 */
#define MCCAOTH_AER               2   /* this has to be CAOTH_AER+1 */
#define MCCAOTH_FIR               3   /* this has to be CAOTH_FIR+1 */

#define CAOTH_AER                 1   /* this has to be MCCAOTH_AER-1 */

#define MCSCAT_MOL                0
#define MCSCAT_AER                1
#define MCSCAT_HG1                2
#define MCSCAT_HG2                3
#define MCSCAT_PFT                4
#define MCSCAT_SKIP               5

/* status of photon */
#define MCSTATUS_DEFAULT          0
#define MCSTATUS_TRAVEL           1
#define MCSTATUS_SCATTER          2 
#define MCSTATUS_ABSORB           3 
#define MCSTATUS_SURFACE          4 
#define MCSTATUS_OUTOFDOMAIN      5 
#define MCSTATUS_BOUNDARY_LOWER   6 
#define MCSTATUS_BOUNDARY_UPPER   7 
#define MCSTATUS_HIT_LIDAR        8 
#define MCSTATUS_HIT_VIS_FOD_CONE 9 
#define MCSTATUS_HIT_DET_PLANE    10
#define MCSTATUS_VIRTUALSCATTER   11
#define MCSTATUS_VISQIDD_ACTION   12
#define MCSTATUS_RISQIDD_ACTION   13
#define MCSTATUS_INCREMENT        14
#define MCSTATUS_CLONE            15
#define MCSTATUS_SPLIT            16
#define MCSTATUS_PURGE            17

/* mode of Scattering, e.g. Delta Scaling */
#define MCSC_MODE_NORMAL          0
#define MCSC_MODE_DELTA_SCALE     1

/* mode of Real Importance Sampling */
#define MCRIS_MODE_NORMAL         0
#define MCRIS_MODE_MAS            1

/* mode of Virtual Importance Sampling */
#define MCVIS_MODE_NORMAL         0
#define MCVIS_MODE_QIDD           1

/* mode of Detector Directional Importance Sampling */
#define MCDDIS_NONE               0
#define MCDDIS_UPF                1
#define MCDDIS_UDA                2

/* version of jacobian */
#define MC_JAC_VAR_BETA           0
#define MC_JAC_VAR_ABS            1

/* number of phase matrix elements */
#define NPHAMAT                   6

typedef enum {NONE, UPPER, LOWER} BOUNDARY;

typedef struct {
  double a;
  double b;
  double c;
  double d;
  double upper;
} surface;

/******************************************************************************************/
/* Description of profile structures:                                                     */
/* =================================  						          */
/*  											  */
/* The atmospheric profiles have become a much more compact and flexible structure. 	  */
/*  										          */
/* A typical call to a profile data point will look like this: 			          */
/*  											  */
/* atmos->ksca   [p->SC_mode][p->RIS_mode] ->prof [ist][p->kc] 			          */
/* atmos->ksca3D [p->SC_mode][p->RIS_mode] ->prof [ist][p->kc][p->ic][p->jc] 		  */
/*  											  */
/* Here, ksca and ksca3D stand for the scattering rates in 1D and 3D, resp., 		  */
/* p->SC_mode is the mode of scattering, e.g. NORMAL or DELTA_SCALE, 			  */
/* p->RIS_mode is the mode of real importance sampling, e.g. NORMAL or 		          */
/*                                 enhancement of molecular scattering, 		  */
/* prof[ist] is the profile structure of particle type ist, 				  */
/* ist is the particle type, e.g. WC, IC, MOL, AER, TOT (=total) 			  */
/*  											  */
/* A new variable is kext: 								  */
/*  											  */
/* atmos->kext   [p->SC_mode][p->RIS_mode][p->VIS_mode] ->prof [ist][p->kc] 		  */
/* atmos->kext3D [p->SC_mode][p->RIS_mode][p->VIS_mode] ->prof [ist][p->kc][p->ic][p->jc] */
/*  											  */
/* Here, only ist=TOT exists. 							          */
/* p->VIS_mode is the mode of virtual importance sample, e.g. NORMAL or 		  */
/*                                    enhancement inverse to detector distance 	          */
/*  											  */
/* kext... is the interaction rate that is used to reduce tau in travel_tau 		  */
/*  											  */
/* ksca...[RIS_mode != NORMAL] is the interaction rate used to choose scattering 	  */
/*        ( if ksca < kext, virtual scatter can happen. These do local and 		  */
/*          directional estimates, but no change of photon direction       ) 		  */
/*  											  */
/* ksca...[p->RIS_mode == NORMAL] is the true scatter rate 				  */
/*        ( if ksca[NORMAL] < ksca[RIS_mode] the photon weight changes: 		  */
/*          the ratio times exponential of the difference times pathlength )  	          */
/*  											  */
/******************************************************************************************/

typedef struct {
  int n_caoth;
  int n;
  double **prof;
} profile;

 
/* profile3D is required to store the profiles of  */
/* absorption and scattering coefficients; float   */
/* should be accurate enough                       */

typedef struct {
  int n_caoth;
  int Nz;
  int Nx;
  int Ny;
  int **threed;    /* flag if a cloud/layer is 3D    */
  int *nthreed;     /* number of 3D layers per cloud */
  float ****prof;
  int *tocalloc;
} profile3D;


/* ddprofile3D is required to store the profile of */
/* the absorbed energy; as all output quantities,  */
/* this parameter is also stored in double rather  */
/* than float                                      */

typedef struct {
  int Nz;
  int Nx;
  int Ny;
  int *threed;
  int nthreed;
  /* double ***aer; */
  /* double ***ozo; */
  /* double ***ice; */
  double ***cld;
  double ***tot;
} ddprofile3D;


typedef struct {
  double mol;
  double aer;
  double cld;
  double ice;
  double tot;
} optical_depth;


/* FE - struct to store results of cloud properties sampling */
#ifdef CLDPRP
  typedef struct {
    int sample_cldprp;
    double **reff_wc;
    double **reff_ic;
    double **rhit_wc;
    double **rhit_ic;
    double **tau_wc;
    double **tau_ic;
    double **dxlwc;
    double **dylwc;
    double **dzlwc;
    double **dxiwc;
    double **dyiwc;
    double **dziwc;
    double **totweights;
    double **wc_weights;
    double **ic_weights;
  } result_cldprp_struct;
#endif

  
typedef struct {
  int    **ndir;    /* number of photons     */
  int    **ndn;
  int    **nup;

  double **edir;     /* irradiance           */
  double **edn; 
  double **eup;  

  double **fdir;     /* actinic flux         */
  double **fdn; 
  double **fup;  

  /* CE included additional dimension for all radiances to account for polarisation */
  double ****raddir;  /* radiance             */
  double ****raddif;
  double ****radesc;

  /* introduced for I3RC case7: radiance binned in radial and pathlength bins */
  double ****radpat; /* cone sampling  */
  double ****radpes; /* local estimate */

  /* standard deviations */
  int std_alloc;     /* flag if allocated */
  double **edir2;    /* irradiance        */
  double **edn2; 
  double **eup2;  
  
  double **fdir2;    /* actinic flux      */
  double **fdn2; 
  double **fup2;  

  double ****raddir2; /* radiance          */
  double ****raddif2;
  double ****radesc2;
  double ****radpat2; 
  double ****radpes2; 

  int ****nraddir;    /* number of photons */
  int ****nraddif;    /* (radiance)        */
  int ****nradesc;
  int ****nradpat;
  int ****nradpes;

#ifdef CLDPRP
  result_cldprp_struct cldprp;
#endif
  
} radiation_field;

  
typedef struct {
  double k;
  double rho0;
  double theta;
  double scale;
  double sigma;
  double t1;
  double t2;
} rpv_struct;


/* Replaced by structures in cdisort.h
typedef struct {
  double iso;
  double vol;
  double geo;
  int    hotspot;
  int isCaM;
  } rossli_struct;

typedef struct {
  double h;
  double b0;
  double w;
} hapke_struct;
*/


typedef struct {
  int method;
  char filename[FILENAME_MAX];
  char spectralfilename[FILENAME_MAX];

  double **albedo2D;
  double albedo;

  double *cos_lat;
  double *X;
  double *Y;
  double xmax;
  double ymax;

  double delX;
  double delY;
  int Nx;
  int Ny;
  double Wsurf;
  float bplkavg;

  int reflectalways;
  int spherical3D_scene;

  char rpv_filename [FILENAME_MAX];
  rpv_struct *rpv;   /* rpv structure                              */
  int *rpv_count;    /* number of occurrence of each RPV index     */
  unsigned char **rpv_index;  /* xy-field of RPV indices                    */

  double *alb_type;

  /* Cox and Munk */
  float u10;
  float pcl;
  float xsal;
  float uphi;
  int   solar_wind;

  char rossli_filename [FILENAME_MAX];
  rossli_brdf_spec **rossli2D;  /* xy-field of Rossli BRDF structure */
  rossli_brdf_spec *rossli;     /* homogeneous Rossli BRDF */
  hapke_brdf_spec *hapke;       /* homogeneous Hapke BRDF */

  double *spectral_albedo;      /* ALIS spectral albedo [iv]           */
  double **spectral_alb_type;   /* ALIS spectral albedo types [il][iv] */

} albedo_struct;


/* CB - struct to store various stuff needed for coherent backscattering */
typedef struct{
  double start[3];       /* coordinates for start, end and the two first and two */
  double sca1[3];        /* last scatters of the normal path */
  double sca2[3];
  double scam[3];
  double scan[3];
  double end[3];

  int Zfirst;              /* indicates if first scatter has happend */
  double phamfirst;
  double Zlast[4][4];      /* saves last phase matrix */
  double Zfirstalt[4][4];  /* saves last phase matrix of the reverse path */
  double Zmiddlealt[4][4]; /* saves the reverse product of the middle phase matrices */
  double stokes_alt[4];    /* normal stokes vector of the reverse path */
  double stokes_cb[4];     /* cb-part of the stokes vector for both paths */
  double ddx;              /* path difference of different paths times wavenumber */
  double pld;              /* path difference of different paths */
  int cb;                  /* indicates that cb SVs have been calculated succesfully */

  double tauext_mid;       /* tauext and tausca needed to find tauabs for reverse path */
  double tauext_sta;
  double tauext_las;
  double tauext_sta_alt;
  double tauext_end_alt;

  double tausca_sta_alt;

  double tausca;           /* weights needed besides the phase matrix */
  double tauext;
  double weight;
  double pdir;

  double weight_hit;
  double cosalpha;

  double phamat0[6];
  double phamat_save0[4][4];

  double mu_ini;

  double start2[3];
  double end2[3];

  double Zlastalt[4][4];
} cohebasca_struct;

/* CB - struct to save certain int values from photon_struct needed for the */
/* computation of the phase matrices for the reverse path */
typedef struct{
  int kc_f;         /* grid box indices of the first scatter of the original path */
  int ic_f;
  int jc_f;
  int kc_s;
  int ic_s;
  int jc_s;
  int SC_mode_f;    /* scatter mode of the first scatter of the original path */
  int scaco_f;      /* scattercounter of the first scatter of the original path */
  double weight_f;  /* save first p->weight */
  double weight_sl; /* save second-to-last p->weight */
  double tausca_f;
} p_save_struct;


/* FE - struct to sample cloud properties during photon tracing */
#ifdef CLDPRP
typedef struct {	
  double reff_wc;
  double reff_ic;
  double rhit_wc;
  double rhit_ic;
  double tau_wc;
  double tau_ic;
  double dxlwc;
  double dylwc;
  double dzlwc;
  double dxiwc;
  double dyiwc;
  double dziwc;
} sample_cldprp_struct;  
#endif
  
typedef struct photon_struct *photon_sub_struct;
                    /* this trick is needed in order to define a */
		    /* photon_struct within locest_struct within */
                    /* a photon_struct                           */


typedef struct {
  double dx[3];
  int hit[3];
  double cotheta; 
} direction;


typedef struct {
  direction dir;           /* direction                                  */

  double cosalpha;         /* cosine of detector surface angle           */
  double pdir;             /* probability                                */
  double pdir_iso;         /* probability for isoenergetic scattering    */
  double **pdir_sct;

  double dist;             /* distance to hit point (detector/beam end)  */
  double distinv;          /* inverse distance to hit point              */

  double r_det;
  double z_det;
  double x_cc[3];
  double t_det;

  double x_hit[3];
  double weight_hit;

  int in_cone;
  int will_hit_cone;
  int behind_detector;
  int will_hit_det_plane;
  int lidar_outside_grid;

  double hit_det_plane_step; /* distance till photon hits detector plane    */
  double vis_fod_step;       /* distance till something happens for VIS-FOD */
  double vis_fod_step2;      /* distance till something happens for VIS-FOD */
  double vis_fod_kext;       /* interaction rate inside VIS-FOD cone        */

  double hitpoint[3];        /* the detector hit point to be hit in         */
                             /* local estimator mode                        */

} locest_struct;

/* structure to store the number of local estimates of certain weight "add" */
struct tnode_LE {
  int add;	       /* approximate weight of bin */
  int nevents;         /* sum of events within bin  */
  struct tnode_LE *left;  /* left child                */
  struct tnode_LE *right; /* right child               */
};

typedef struct {
  int n_caoth;
  int ntypes;
  int norders;
  int Nz_jac;
  int nz;
  int nt;
  int np;
  int no;
  int na;
  double ******radloc;     /* radiance for lidar                      */
  double ******radloc2;    /* standard deviation                      */
  double *******radloc_t;  /* radiance for lidar                      */
  double *******radloc_t2; /* standard deviation                      */
  double ******jacobian_t;
  double ******jacobian_t2;
  double *****jacobian_t2act;
  int *actis;
  double **photoncounter;  /* number of photons that are simulated    */
  struct tnode_LE ****dtree; /* local estimator size count tree    */

  long int** localEstimates;  /* number of local estimates per range gate for radar */
  double** weightedPathlengths;  /* cumulative weighted pathlengths for radar */
  double** weights;      /* cumulative weights for radar */
} lidar_result_field;

typedef struct {
  int n_caoth;
  int Nx;
  int Ny;
  int Nz_jac;
  double *****jacobian_t;
  double *****jacobian_t2;
  double ***jacobian_tact;
//  double **photoncounter;  /* number of photons that are simulated    */
//  struct tnode_LE ****dtree; /* local estimator size count tree    */
} jacobian_result_field;

typedef struct {
  int Nc;
  int Nv; 
  int Nx; 
  int Ny; 
  int Np;
  double *****rad_t;        /* spectral concentration dependent radiation field 
			       [Nsc,Nv,Nx,Ny,Np] */
} radiation_field_t;

typedef struct {
  double R1t[6];
  double RMt[6];
  double M_Elements[4][5];
  double M_Elements_S0[4][5];
  double M_Elements_SS0[4][5];
  double M_Elements_SS0_t[4][5];                                                                                                                                      
  long int R1C;                                                                                                                                                       
  long int RMC;                                                                                                                                                       
} mishchenko_cb_struct;                                                                                                                                               
                                                                                                                                                                      
typedef struct {                                                                                                                                                      
  double mu_ave;                                                                                                                                                      
  long int mu_counter;                                                                                                                                                
  double path_ave_ris;                                                                                                                                                
  long int path_counter;                                                                                                                                              
  double meanPLD;                                                                                                                                                     
  long int cPLD;   
} lidar_cb_struct;     

typedef struct {
  radiation_field *surf;
  radiation_field **alt;
  lidar_result_field *lidar;
  jacobian_result_field *jacobian;
  radiation_field_t *surf_t;
  radiation_field_t **rad_t;   /* importance sampling (spectral+concentration)*/
  mishchenko_cb_struct* mish;
  lidar_cb_struct* lidcb;
    
  double ***back;
  double ***back2;     /* BCA: Maybe a comment on what is back2? */
  double ***back_dEnet; /* **CK thermal backward heating rates*/
  double ***back_dEnet2; /* **CK thermal backward heating rates; standard deviation*/
  double **backemis;
  double *****back_t;   /* importance sampling (spectral+concentration) backward */
    
  ddprofile3D *absorption3D;    /* average            */
  ddprofile3D *absorption3D2;   /* standard deviation */
    
#ifdef CLDPRP
  int sample_cldprp;
  float wavelength;
#endif
    
  double **circcontr;
  
  double* pathlength_per_layer_tot; 

} result_struct;


/* Structure to store a direction and opening angle */
/* for radiance calculations; convention:           */ 
/*    phi=0,   South                                */
/*    phi=90,  West                                 */
/*    phi=180, North                                */
/*    phi=270, East                                 */
/*                                                  */
/*    theta > 0 upward                              */
/*    theta < 0 downward                            */

typedef struct {
  double theta;         /* polar angle                      */
  double phi;           /* azimuth, internal                */
  double externalphi;   /* azimuth, as defined by the user  */
  double alpha; 
  double cosalpha; 
  double omega; 
  direction dir;
} radang;

typedef struct {
  double theta;         /* polar angle                      */
  double phi;           /* azimuth, internal                */
  double externalphi;   /* azimuth, as defined by the user  */
  double alpha; 
  double cosalpha; 
  double omega; 
  direction dir;

  double x[3];          /* position                           */
  double radius;        /* radial size of laser/detector area */
  double area;          /* laser/detector area                */
  int divergence;

} lidar_laser;

typedef struct {
  double theta;         /* polar angle                      */
  double phi;           /* azimuth, internal                */
  double externalphi;   /* azimuth, as defined by the user  */
  double alpha; 
  double *cosalpha; 
  double *omega;
  double *azphi;
  direction dir;
  double dx_perp1[3];
  double dx_perp2[3];

  double x[3];          /* position                           */
  double radius;        /* radial size of laser/detector area */
  double area;          /* laser/detector area                */

  double s_det;
  double c_det;
  double t_det;
  double z_det;
  double x_cc[3];

} lidar_detector;


/* structure to store the sum of photon weights for an individual photon */
struct tnode {
  double *loc;	  /* pointer to storage    */
  double weight;  /* sum of photon weights */
  struct tnode *left;       /* left child  */
  struct tnode *right;      /* right child */
};

/* structure to store the sum of photon weights for an individual photon */
struct anode {
  double *loc;	  /* pointer to storage    */
  double *array;  /* sum of photon weights */
  struct anode *left;       /* left child  */
  struct anode *right;      /* right child */
};


/* linked list element to store a photon path */
struct photon_path {
  double x[3];
  struct photon_path *next;
};


typedef struct photon_struct {
  /* current status of photon */
  double x[3];               /* photon location                              */
  direction dir;             /* photon direction                             */
  direction dir0;            /* initial photon direction                     */
  direction dir00;            /* initial photon direction                    */
  direction dir_ort;         /* orthogonal direction to scattering plane     */

  double lest_phi; 	/* Local estimate direction in regard to center of sun */
  double lest_alpha; 	/* Local estimate direction in regard to center of sun */
  double dstep_cone;
  
  int ic;                    /* box index, x                                 */
  int jc;                    /* box index, y                                 */
  int kc;                    /* box index, z                                 */
  int jcg;                   /* box index, y, for spherical 3D               */
  int iv_alis;

  int photon_status;         /* flag: current status of photon               */
  int reallyhit[3];          /* flag: crossed box boundary                   */
  int muchoutcounter;

  int SC_mode;               /* actual mode for Scattering                   */
  int DDIS_SC_mode;          /* actual mode for DDIS Scattering              */
  int RIS_mode;              /* actual mode for Real Importance Sampling     */
  int VIS_mode;              /* actual mode for Virtual Importance Sampling  */

  /* counters */
  int direct;                /* direct or diffuse                            */

  int scattercounter;        /* number of scattering + reflection events     */ 
  int reflectcounter;        /* number of reflection events     */ 
  int photoncounter;         /* number of current photon                     */
  int backward_is;           /* number of backward index                     */
  int backward_js;           /* number of backward index                     */
  int rayleighcounter;       /* count number of rayleigh scatters            */
  int isclone;
  int clonescattercounter;
  int escapescattercounter;
  int spikewarningcounter;

  int update_atmos_coord;    /* whether to update atmos_coord in step3d (spherical3d) */

  double pathlength;         /* photon pathlength                            */
  optical_depth tauabs;      /* absortion optical depth along the path       */
  double tauris;
  double p_norm;
  double *pathlength_per_layer; /* photon pathlength per layer               */       

  /* weights */
  double weight;              /* general weight, contains what is not below   */
  double weight1;             /* general weight, contains what is not below   */
  double weight_emis;         /* weight for emission (thermal backward) **CK 2013 */
  double phi0;                /* initial azimuthal direction, needed for pol  */
  double fw_phi0;             /* sun azimuthal direction, needed for pol      */
  double fw_phi;              /* detector azimuthal direction, needed for pol */
  double *stokes0;            /* initial stokes vector for polarization       */
  double *stokes;             /* stokes vector for polarization               */
  double **phamat;            /* phase matrix of photon                       */
  double q_isoene;            /* weight for isoenergetic scattering           */
  double ***q_jacobian;       /* weight for jacobian                          */
  double *r_jacobian;         /* pathlength for jacobian                      */
  int Nz_alis;                /* number of atmosph. layers, needed by ALIS or boxairmass */
  int Nz_jac;                 /* number of atmosph. layers, needed by jacobian*/
  int n_caoth;                /* number of caoth, needed by jacob.            */
  /* CE replaced the following by pathlengths in each layer */
  /* double *dtauabs_spectral;    weight vector for spectral absorption       */ 
  int nlambda;                /* number of wavelengths included in q_spectral */
  double *q_spectral;         /* weight for spectral calculation, absorption and scattering coefficient  */
  double *q2_spectral;        /* weigh for spectral variation of phase matrix */
  double *q_albedo_spectral;  /* weight for spectral albedo                   */
  int    Nc;                  /* number of scaling factors                    */ 
  double *q_concentration;    /* weight for concentration importance sampling */
  double *q2_concentration;    
  double special_weight;

  /* stuff for local estimator / escape radiances */
  double *pdir;              /* probability for direction radang[]           */
  locest_struct lest;        /* variables needed for local estimator         */

  /* honestly, these should not be in this structure */
  int ipa;                   /* independent pixel approximation              */
  double maxpathlength;      /* evtl only allow this as maximal pathlength   */

  /* miscellaneous */
  struct tnode *wtree;       /* photon weight tree                           */
  struct anode *wtree2;       /* photon weight tree                          */
  struct photon_path *path;  /* photon path linked list                      */

  double vis_beta;
  double risqidd_beta;

  double tauext_tot;
  double pathlength2;
  double pathlength3;

  photon_sub_struct parent_photon;

  /* coherent backscatter structures */
  cohebasca_struct cb;  
  p_save_struct pss;
  
  /* structure for sampling of cloudproperties by photons */
#ifdef CLDPRP
  sample_cldprp_struct cldprp;
#endif
 
  int dEnet_component;  /* photon component for thermal backward heating rates - indicating the starting direction/area  **CK */
  double n_dEnet;  /* denet method - number of photons (depending on face area for std calculation  **CK */

  int doling;

  int ncircos;
  double **tocirco;

} photon_struct;

  
typedef struct {
  int user_defined_btemp;      /* flag if user defined surface_temperature  */
  float btemp;                 /* constant surface temperature              */

  int surf2D; /* 2D temperature field */
  double **temp2D;
  double **plkavg2D;
  double delX;
  double delY;
  int Nx;
  int Ny;
} surftemp_struct;

typedef struct {
  int elev2D;

  surface **surf;
  double surfmax;
  double delX;
  double delY;
  int Nx;
  int Ny;

} elevation_struct;

typedef struct {
  char   *name;              
  double fraction;            /* fraction of this crystal type                           */
  double oriented_fraction;   /* fraction of oriented crystals                           */
  double angdist_width;       /* width of the angular distribution of oriented particles */
  int orientation_dof;        /* degrees of freedom of oriented particles                */
} crystal_prop_struct;        /* crystal properties for raytracing                       */


typedef struct {
  double delX;
  double delY;
  double *X;
  double *Y;
  float  *Z;   /* all altitudes are floats to avoid roundoff problems */
  double *Yg;
  double *sin2_lat;
  double *tan2_lat;
  int *realY;
  int Nx;
  int Ny;
  int Nz;
  int Nyg;
  int jc_eq;
  double xmax;
  double ymax;
  int **threed;
  int *nthreed;
  int nscaDS;   /* number of scatter modes, e.g. delta-scaling */
  int nscaRIS;  /* number of RIS scatter modes */
  int nscaVIS;  /* number of VIS scatter modes */
  double *ddis_eps;
  double r_earth;
  int spherical3D_scene;
  int n_caoth;
  char **caoth_name;
  int *scatter_type;

  profile ***ksca;
  profile *kabs;
  profile ****kext;
  profile *g1;
  profile *g2;
  profile *ff;
  profile *reff;
  profile3D *Bplanck;
  
  double ***kabs_spectral; /* spectral absorption coefficients */ 
  double ***ksca_spectral;
  /* double **kabs_aer_spectral; /\* spectral absorption coefficients aer *\/  */
  /*   double **ksca_aer_spectral; */
  /*   double **kabs_wc_spectral; /\* spectral absorption coefficients wc *\/  */
  /*   double **ksca_wc_spectral; */
  /*   double **kabs_ic_spectral; /\* spectral absorption coefficients ic*\/  */
  /*   double **ksca_ic_spectral; */

  int nlambda_abs;
  float *lambda;
  double **Bplanck_spectral;
  int nlambda_ref;
  int *ilambda_ref;
  
  double **kabs_scaled; /* concentration importance sampling */ 
  double **ksca_scaled;
  int Nc;

  profile3D ***ksca3D;
  profile3D *kabs3D;
  profile3D ****kext3D;
  int ***spiky_box;

  profile3D *g1_3D;
  profile3D *g2_3D;
  profile3D *ff_3D;
  profile3D *reff_3D;

  int *i_wc;     // indices where profiles for water are stored
  int *i_ic;     // indices where profiles for ice  are stored
  
#ifdef CLDPRP

  int sample_cldprp;  
  profile3D *dxlwc_3D;
  profile3D *dylwc_3D;
  profile3D *dzlwc_3D;

#endif
  
  crystal_prop_struct *raytracing_prop;
  int n_raytracing_prop;
  
  double rayleigh_depol; /* Rayleigh depolarisation */

  phase_function_table **phase;/* complete cloud single scattering properties */
  pft *phase_aer;             /* aerosol layer single scattering properties */

  double *tabstot;
  double *tscatot;

  double Watm;
  double maxemis;

  int nztilt;
  float *tipa_transmission;
  float ***tipa_transmission3D;
  float *tipa_P_start_level;
  float ***tipa_P_start_level3D;

  double ris_factor;

} atmosphere_struct;


typedef struct {
  int *sample;
  int nzout;
  radang *rad;
  direction rad_dir0; /* backup of rad.dir */
  double **umu2D;
  double **phi2D;
  int Nd;         /* number of angular directions */
  int Nr;         /* number of radii              */
  int Nt;         /* number of time intervals     */
  float dr;       /* radius interval              */
  float dt;       /* time interval                */
  double delX;
  double delY;
  int Nx;
  int Ny;
  int surface;
  int escape;
  int delta_scaling;
  int tipadir;
 
  int coherent_backscatter;        /* CB */
  int coherent_backscatter_lmode;  /* CB lidar mode */
  double ris_factor;
  int cldprp;
  double ris_optical_depth;
 
  double escape_eps_ddis_upf; /* DDISsing probability */
  int LidarLocEst;        /* local estimator for Lidar */
  int ili;                /* active index of lidar detectors    */
  int reference_to_NN;        /* referencing to NN */

#if HAVE_VROOM
  int vroom;
  int LLE_RIS_MAS;        /* Real Importance Sampling for Molec. and Aerosol Scat. */
  int LLE_D_DIS;          /* Detector Directional Imp. Sampl.                      */
  int LLE_VIS_FOD;        /* Virt. Imp. Sampl. in FOV Of Detector                  */
  int LLE_VIS_QIDD;       /* Virtual IS according to Quadr. Inv. Dist. to Det.     */
  int LLE_sponti;         /* Spontaneous cloning                                   */
  int LLE_turnmax;        /* Turn off when converged                               */
  int LLE_moon;        /* Replace laser by moon or sun                           */
  int RIS_MS;
  int VIS_CS;
  int splitter;
  int ntupelLE;
  int startCP;
  int LEperCP;
  int use_p_norm;
  double n_split_max;
  double n_split_min;
  int LLE_jacobian;
  int abs_jacobian;
  int LLE_jacobian_std;
  int abs_jacobian_std;
  int vroomreflectalways;

  int LLE_channels;
  int LLE_stokes;
  int LLE_pol_rings;

  int LLE_RIS_MAS_start;
  int LLE_RIS_MAS_end;

  double LE_taucrit;      /* Iwabuchi critical tau                                 */

  int Nli;                /* number of lidar detectors          */
  int LLE_Nt;             /* number of time bins for lidar      */
  int LLE_No;             /* number of opening angles for lidar */
  int LLE_Na;             /* number of ring segments for lidar  */
  int LLE_islinear;       /* whether or not opening angles increase linearly */
  int LLE_polarisation;   /* polarisation of the laser beam     */
  int LLE_scatterout;
  double LLE_afac;        /* factor of opening angles for lidar */
  double LLE_pulse;       /* pulse length of the laser beam     */
  lidar_laser *laser;     /* emitter for local estimator        */
  lidar_detector *lidar;  /* detector for local estimator       */
  double *lidar_t;        /* time bins for differential loc.est.*/
  double *lidar_dt;       /* width of time bins                 */
  int n_phase_max;
  pft **phase_max;

  double LLE_eps_ddis_upf;     /* DDISsing probability */
  double LLE_eps_ddis_uda;     /* DDISsing probability */
  double LLE_eps_fod_dis_phi;
  double LLE_taumax;
  double jac_eps;

  double visqidd_betamax;
  double visqidd_rmax;
  double visqidd_facs;
  double risqidd_betamax;
  double risqidd_facs;

  double MPdynDDIS;
  double split_max;
  double split_min;
#endif

  int    maxscatters;        /* evtl only allow this as maximum number of scatters */
  int    minscatters;        /* evtl only allow this as maximum number of scatters */
  double *blitz_position;
  double maxpathlength;

  int passback3D;
  int std;        /* calculate standard deviation */

  int spectral_is;  /*spectral calculation with importance sampling */
  int concentration_is;  /* compute several concentrations with importance sampling */
  int boxairmass;

  int spherical;
  int spherical3D;
  int spherical3D_scene;
  double spherical3D_scene_lon_min;
  double spherical3D_scene_lon_max;
  double spherical3D_scene_lat_min;
  double spherical3D_scene_lat_max;
  int refraction;

  int ixmin;
  int ixmax;
  int iymin;
  int iymax;

  int polarisation; 
  int polarisation_state; 
  int nstokes;         /* Number of Stokes components       */   
  int nphamat;         /* Number of Phase matrix components */   

  int bcond;      /* boundary conditions */

  int backward; 
  int backward_is;
  int backward_js;

  int backward_islower;     /* sample pixels to be calculated; */
  int backward_jslower;     /* lower left - upper right        */
  int backward_isupper;
  int backward_jsupper;

  int backward_islower_org; /* same as above, to store original */
  int backward_jslower_org; /* data defined in the input file   */
  int backward_isupper_org; /* (they may be changed e.g. by     */
  int backward_jsupper_org; /*  forward2backward())             */

  int backward_writeallpixels;
  int backward_writeback;

  int backemis_kc;

  /*** forward tracing                             ***/
  double forward_sza;    /* solar zenith angle       */
  double forward_phi0;   /* solar azimuth angle      */
  double forward_vza;    /* viewing zenith angle     */
  double forward_phi;    /* viewing azimuth angle    */     

  float forward_zout;    /* output altitude          */
  float forward_zstart;  /* start altitude           */

  /*** backward tracing                            ***/
  double backward_sza;   /* solar zenith angle       */
  double backward_phi0;  /* solar azimuth angle      */
  double backward_vza;   /* viewing zenith angle     */
  double backward_phi;   /* viewing azimuth angle    */     
  double phi_target;     /* target azimuth angle, needed for polarisation */

  float backward_zout;    /* output altitude         */
  float backward_zstart;  /* start altitude          */
  int sample_forward_kc;  /* forward sampling layer  */
  int sample_backward_kc; /* backward sampling layer */

  float* sample_backward_sunshape_p; /* sunshape distribution function for direct rad*/
  float* sample_backward_sunshape_p2; /* sunshape distribution function for diffuse rad*/
  double* sample_backward_sunshape_F; /* cumulative distribution function of sunshape_p2 */
  double* sample_backward_sunshape_alpha; /* Angle alpha from center of sun in radians */
  double* sample_backward_sunshape_slope; /* Slope of sun shape */
  int    sample_backward_sunshape_n_lines; /* Number of lines in sun shape file */
  int    sample_backward_sunshape; /* Switch */
  double sun_radius; /* Sun radius in degrees */

  /* actual angles */
  double sza;
  double phi0;
  double vza;
  double phi;

  float zstart;                 /* photon start position         */

  char backward_altstr[255];    /* store the altitude level for the output filename */

  int panorama;            /*  panorama flag                */
  int panorama_forward;    /*  panorama flag                */
  int pan_alignment;       /*  panorama flag                */
  int pan_no_pixel;        /*  panorama flag                */
  int pan_quicklook;       /*  panorama flag                */
  int pan_distr_photons_over_pixel; /*   panorama flag      */
  int pan_circumsolar_var_red;       /*   panorama flag               */ 
  int pan_weight_with_cos; /*   panorama flag               */
  int pan_with_direct_rad; /*   panorama flag               */
  double pan_umu_min;      /* panorama boundaries           */
  double pan_umu_max;      /* panorama boundaries           */
  double pan_theta_min;    /* panorama boundaries           */
  double pan_theta_max;    /* panorama boundaries           */
  double pan_phi_min;      /* panorama boundaries           */
  double pan_phi_max;      /* panorama boundaries           */
  double pan_dphi; 	/* panorama pixel width */
  double pan_dtheta; 	/* panorama pixel heigth */
  double pan_FOV;          /* Field of view of the actual panorama pixel in sr */
  double pan_picture;
  double align_theta;
  double align_phi;
  double add_phi_when_aligning;

  int Nz_jac;

  int sensorposition;           /* user-defined sensor position  */
  double sensorposition_x[3];   /* only works in backward mode   */

  int sensordirection;          /* user-defined sensor-direction */
  double sensordirection_dx[3]; /* only works in backward mode   */

  int surfaceparallel;          /* sample surface-parallel       */

  double **surface_area;        /* surface area of sample pixels (2D elevation, surface-parallel, forward) */
  double **surface_area_std;    /* standard deviation                                                      */
  int **surface_area_counter;   /* counter for averaging                                                   */
  /*  int lchatta; */
  
  int DoLE;

  int ncirc;
  double dcirc;

  /* **CK */
  int n_xy;  /* number of photons distributed on xy-plane of a volume pixel for thermal backward heating rates */
  int n_xz;  /* number of photons distributed on xz-plane of a volume pixel for thermal backward heating rates */
  int n_yz;  /* number of photons distributed on yz-plane of a volume pixel for thermal backward heating rates */
  double weight_heat ;  /* number of photons distributed on yz-plane of a volume pixel for thermal backward heating rates */
  int **heat_flag;  /* flag for heating rate method  **CK */

 
  int planar_switch;
  int mish_switch;
  int spec_otp_switch;
  int nocb_switch;
  int use_phase_matrix_total;
  int det_eq_tra;                                                                                                                                                     
  double wavelength;                                                                                                                                                  
  double wavenumber;                                                                                                                                                  
  long int nphotons;  

  double pathlength;
  
} sample_struct;

typedef struct {
  int                  first;
  /* spectral importance sampling (molecules) */
  float                *lambda;
  int                   nlambda_abs;
  int nlambda_ref;
  int *ilambda_ref;
  /* int                   rayleigh_crs; */
/*   float                 rayleigh_depol_user; */
/*   float                *dens; */
/*   float                 mixing_ratio_co2;  */
  double             ***dt;
  double             ***om;
  double               *albedo;      /* spectral albedo albedo [iv]      */
  double               **alb_type;   /* 2D albedo field alb_type[il][iv] */
  /* concentration importance sampling (aerosols) */
  int                   Nc;
  double               *aer_scaling_factors;
} alis_struct;



/* 3d caoth structure; included in uvspec caoth_out_struct; */
/* used by MYSTIC and IPA                                   */
typedef struct {
  char *name;
  char *fullname;

  int Nx;            /* number of grid points in x-direction */
  int Ny;            /* number of grid points in y-direction */
  int nlyr;          /* number of layers                     */
  
  double delX;       /* grid size in x-direction    [m] */
  double delY;       /* grid size in y-direction    [m] */
  
  float *zd;         /* level altitudes, 0 ... nlyr [m] */
                     /* all altitudes are float to      */
                     /* avoid roundoff problems         */

  float ***lwc;      /* liquid water content    [g/m3]   */
  float ***reff;     /* effective radius        [micron] */

  float ***ext;      /* extinction coefficient  [1/m]       */
  float ***ssa;      /* single scattering albedo            */
  float ***g1;       /* asymmetry parameter 1 for double-HG */
  float ***g2;       /* asymmetry parameter 2 for double-HG */
  float ***ff;       /* forward fraction for  double-HG     */
  float ***f;        /* delta scaling factor                */
  float ***dscale;   /* delta-scaling factor Tobias         */ 
  float reffmin;     /* minimum effective radius in caoth   */
  float reffmax;     /* maximum effective radius in caoth   */
  
  float ***tausol3d;  /* contains tau_wc or tau_ic at [nz][nx][ny], necessary for mc_tipa dir (ulrike) */
  int nztilt;         /* at how many levels tilting has to be performed = amount of "tilting-levels" */

  phase_function_table phase;  /* complete single scattering properties: */
                               /* tables of phase function and           */
                               /* cumulative phase function              */            

  int cldproperties; /* as defined in the wcloud3D header                  */
  int nonHG;         /* flag: Henyey-Greenstein or explicit phase function */

  crystal_prop_struct *raytracing_prop;   /* should be moved into microphys or optprop */
  int     n_raytracing_prop;              /* should be moved into microphys or optprop */
  
  int *threed;       /* flag: 3D layer or 1D layer */
  int nthreed;       /* number of 3D layers        */

} caoth3d_out_struct;

/* prototypes */
int mystic (int                  *nlyr,
	    int                   n_caoth,
            float               **dt_s,
	    float               **om_s,
	    float               **g1_s,
	    float               **g2_s,
	    float               **f_s,
            float               **ds_s,
            alis_struct          *alis,
            float                *refind,
            float                 r_earth,
            float                *rayleigh_depol,
            caoth3d_out_struct   *caoth3d,
            float               **re_s,
            float              ***temper,
	    int                   temper3d, 
            float                *zprof,
            float              ***momaer,
	    int                  *nmomaer,
	    int                 **nthetaaer,
	    float              ***thetaaer,
	    double             ***muaer,
	    float              ***phaseaer,
	    int                   nphamataer,
	    int                  *spectral,
            float                *alb,
            float                *alb_type,
	    float                 sza,
	    float                 phi0,
	    float                 sza_spher,
	    float                 phi0_spher,
	    long int              nphotons,
            int                  *source,
	    float                *wvnmlo,
	    float                *wvnmhi,
	    float                *wavelength,
            float                *zoutorg,
	    int                  *nzoutorg,
            float                *rpv_rho0,
	    float                *rpv_k,
	    float                *rpv_theta, 
            float                *rpv_scale,
	    float                *rpv_sigma,
	    float                *rpv_t1,
	    float                *rpv_t2,
	    float                *hapke_h,
	    float                *hapke_b0,
	    float                *hapke_w,
	    float                *rossli_iso,
	    float                *rossli_vol,
	    float                *rossli_geo,
	    int                   rossli_hotspot,
            float                *u10,
	    float                *pcl,
	    float                *xsal,
	    float                *uphi,
	    int                  *solar_wind,
            float                *bpdf_u10, 
            int                  *tenstream,
	    char                 *tenstream_options,
            int                  *ipa,
	    int                  *absorption,
	    int                  *thermal_heating_method,
	    int                  *loaddata,
            sample_struct        *sample,
            elevation_struct     *elev,
            surftemp_struct      *surftemp,
            char                 *basename, 
            char                 *umufilename, 
            char                 *sunshape_filename, 
            char                 *albedo_filename,
            char                 *alb_type_filename,
            char                 *rpvfilename,
	    char                **rpv_labels,
	    int                  *rpv_nlabels,
            char                 *ambralsfilename,
            char                 *rosslifilename,
            float                *delta_scaling_mucut,
            float                *truncate,
            int                  *reflectalways,
            int                   quiet,
            float                *rfldir,
	    float                *rfldn,
	    float                *flup,
            float                *uavgso,
	    float                *uavgdn,
	    float                *uavgup,
            float              ***rfldir3d,
	    float              ***rfldn3d,
	    float              ***flup3d, 
            float              ***uavgso3d,
	    float              ***uavgdn3d,
	    float              ***uavgup3d, 
            float             ****radiance3d,
	    float              ***absback3d,
            float              ***abs3d, 
            float           ******radiance3d_is,
            float              ***rfldir3d_var,
	    float              ***rfldn3d_var,
	    float              ***flup3d_var, 
            float              ***uavgso3d_var,
	    float              ***uavgdn3d_var,
	    float              ***uavgup3d_var, 
            float             ****radiance3d_var,
	    float              ***absback3d_var,
            float              ***abs3d_var, 
	    int                   nxcld,
	    int                   nycld,
	    int                   nzcld,
	    double                dxcld,
	    double                dycld,
            char                 *datapath,
            char                 *randomstatusfilename,
	    int                   readrandomstatus,
	    int                   write_output_as_netcdf,
            int                   write_files,
	    int                   visualize );
  
double sc_Rayleigh_mu();
double sc_Rayleigh_mu_depol (double depol);
double sc_Isotropic_mu();
double sc_Lambertian_mu();
double sc_HG_mu(double g);
double sc_HG2_mu(double g1, double g2, double ff);
double random_tau();
  

double rpv_brdf (double rho0, double k, double theta, 
		 double scale, double sigma, double t1, double t2,
		 double mu1, double mu2, double phi);

double rpv_albedo (double rho0, double k, double theta, double mu1);

int random_rpv_direction (double rho0, double k, double theta, double mu1,
                          double *mu2, double *phi);
  
photon_struct *calloc_photon ( sample_struct *sample,
			       int            Nz,
			       int            source,
			       int            nlambda_abs,
                               int            Nsc,
                               int            n_caoth );

void init_direction (double sinsza, double cossza, double sinphi,
                     double cosphi, 
                     direction *dir);

void gen_default_sunshape(float wvnmlo, float wvnmhi, int N, float *pd, double *alpha); 

int limb_dark_lest_dir (direction *dir_sun_center, sample_struct *sample,
			double *phi, double *alpha);

int random_rossli_direction (double iso, double vol, double geo, 
			     double mu1, double *mu2, double *phi);

void new_direction (double mu, double phi, direction *dir, double phi_inc);

void elev_coord (photon_struct *p, elevation_struct *elev, int *ie, int *je);

#if HAVE_LIDAR
int read_lidar (char *filename,  sample_struct *sample,
		int quiet);
#endif

#if defined (__cplusplus)
}
#endif

/* prototypes of internal functions */

double sc_mu (pft *iphase, int newrand, int scaled, double F_max, double F_min);

void cp_direction (direction *dest, direction *source);


int get_phase_matrix_total ( atmosphere_struct *atmos,
			     photon_struct     *p,
			     double             mu,
			     int                np,
			     int                vismode,
			     int                spectral_is,
			     int                concentration_is,
                             double             ris_factor,
			     int                release,
			     /* Output */
			     double            *phase_matrix,
			     double           **P_caoth,
			     double            *P_spec,
			     double            *weight );

double get_ksca (atmosphere_struct *atmos, photon_struct *p, int isp);
double get_kscaIS (atmosphere_struct *atmos, photon_struct *p, int isp);
double get_kabs (atmosphere_struct *atmos, photon_struct *p, int ist);

int err_out (char *output, int status);
int mem_err_out (char *output, int line, char *fff, char *file);

double sc_Isotropic_phi();

void cp_photon_struct ( photon_struct *dest,
			photon_struct *source, 
			sample_struct *sample,
			int            n_caoth );

void mat_v_mult(double **a, double *b, double *c);
void v_mult  (double *a, double *b, double *c);
void v_sub   (double *a, double *b, double *c);
void v_add   (double *a, double *b, double *c);
void vn_mult (double *a, double *b, double *c);
void v_cross_product (double *a, double *b, double *c);
int v_cross_product_norm (double *a, double *b, double *c);
void v_mult_mu (double *a, double *b, double *c);
void v_fabs (double *a, double *c);
void v_neg (double *a);

double cosd ( double angle );
double sind ( double angle );
double tand ( double angle );
double acosd ( double rad );
double asind ( double rad );
double atand ( double rad );

void destroy_photon (photon_struct *p, int n_caoth);

int direct_radiation (sample_struct     *sample,
		      atmosphere_struct *atmos,
		      photon_struct     *p,
		      result_struct     *result,
		      elevation_struct  *elev,
		      int               *source, 
		      float             *refind);

int photon_journey ( photon_struct     *p,
		     atmosphere_struct *atmos,
		     sample_struct     *sample,
		     result_struct     *result,
		     elevation_struct  *elev,
		     albedo_struct     *albedo,
		     surftemp_struct   *surftemp,
		     int               *absorption,
		     int               *source, 
		     int                plotpath,
		     int                visualize,
		     float             *wvnmlo,
		     float             *wvnmhi,
		     float             *refind,
		     int                quiet );

void random_Lambertian_normal (direction *dir, double *norm);

int mu_scatter ( atmosphere_struct    *atmos,
		 photon_struct        *p,
		 int                   isp,
		 double               *mu );

int sample_cldprp (photon_struct *p, atmosphere_struct *atmos, double step);

int travel_tau (photon_struct *p, atmosphere_struct *atmos,
		sample_struct *sample, result_struct *result,
		elevation_struct *elev, double tau, 
		int absorption, int source, int plotpath, 
		int visualize, float *refind, int quiet);

int step1D (photon_struct *p, atmosphere_struct *atmos, double step, 
	    int bcond, int photonpath, int visualize);

void intersection1D (photon_struct *p,
		     atmosphere_struct *atmos,
		     double tau, double tausca,
		     double *step, int source);
int intersection1D_spherical (photon_struct *p,
			      atmosphere_struct *atmos,
			      int refraction, 
			      float *refind, double *step);

int cross_surface (photon_struct *p, double step,
		   elevation_struct *elev,
		   int bcond,
		   double *gamma);


void phase_matrix_mult ( double          **phase_matrix,
			 double           *P,
			 double           *pdir_sca,
			 double           *pdir_inc,
			 double            phi_sun,
			 double            phi_det,
			 int               escape,
			 int               scattercounter,
			 int               backward );

/* CE: function is not needed in mystic.c anymore, only in lidar.c (coherent backscattering), */
/* should be removed if possible. */
double derive_deltaphi_horz ( double *dx,
			      double   phi_esc,
			      int      scattercounter,
			      int      LidarLocEst,
			      double   phi_init );

struct tnode *addtree_stddev (struct tnode *p, double *loc, double weight);
void treeprint_stddev (struct tnode *p);

void random_Isotropic_normal (direction *dir, double *norm);
double random_cone_mu (double *cotheta);
double random_cone_mu_gauss (double *cotheta);

void atmos_coord  (photon_struct *p, atmosphere_struct *atmos, int *ic, int *jc);

void hitflag (direction *dir);

int set_photon_z (float z, atmosphere_struct *atmos, photon_struct *p);

double get_F (pft *iphase, double mu, int scaled);

double calc_phi_horz (double *dx, double *n);

int find_tau_ext ( photon_struct     *photon, 
		   atmosphere_struct *atmos,
		   sample_struct     *sample,
		   elevation_struct  *elev,
		   int                ic,
		   int                jc,
		   int                kc,
		   float             *refind,
		   double            *dx, 
		   double            *start,
		   double             path,
		   double            *tau_save );

int find_ris_factor (photon_struct *photon,
                     atmosphere_struct *atmos,
                     sample_struct *sample,
                     elevation_struct *elev,
                     float *refind );

int summarize_result_mishchenko_cb (mishchenko_cb_struct* mish, char* basename );
void save_phase_matrix_elements ( mishchenko_cb_struct *mish, double totweight, int scattercounter, double *stokes0 );
void save_phase_matrix_elements_temp ( mishchenko_cb_struct *mish, double **phamat, int scattercounter, double* phaseMatrix );
double kahanSum(double *f,int N);
double fiveKahanSum(double a, double b, double c, double d, double e);
double fourKahanSum(double a, double b, double c, double d);
double threeKahanSum(double a, double b, double c);

void print_sr_error_message(char* histfilename);

#endif
