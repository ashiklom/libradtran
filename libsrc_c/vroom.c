#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include "mystic.h"
#include "vroom.h"
#include "uvspecrandom.h"
#include <time.h>
#include "errors.h"
#if HAVE_LIDAR
#include "lidar.h"
#endif
#ifndef PI
#define PI 3.14159265358979323846264338327
#endif
static inline void q40(scadis_struct*q28,int q41);static int q42(pft**
phase_max,int n_phase_max,scadis_struct q28,int q27,int q43,double*mu);static 
inline int q44(double q45,double*mu);static double q46(double q34,double q23,
int q27,int behind_detector,double q47,pft**phase_max,int n_phase_max,double 
q48,double q38,double q49,scadis_struct q28,int q43);int set_vroom_settings(
int vroom,sample_struct*sample,int q6){switch(vroom){case 1:sample->vroom=1;
sample->escape_eps_ddis_upf=0.1;sample->ntupelLE=22;sample->startCP=4;sample->
LEperCP=11;sample->RIS_MS=0;sample->splitter=1;sample->use_p_norm=1;sample->
split_max=3.0;sample->split_min=0.3;sample->n_split_max=5000.0;sample->
n_split_min=0.2;sample->LE_taucrit=3.0;sample->MPdynDDIS=0.1;sample->VIS_CS=1;
sample->vroomreflectalways=0;if(!q6){fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\124\165\162\156\151\156\147\40\157\156\40\126\122\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\115\56\56\56\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\50\126\141\162\151\141\156\143\145\40\122\145\144\165\143\164\151\157\156\40\117\160\164\151\155\141\154\40\117\160\164\151\157\156\163\40\115\145\164\150\157\144\51\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\111\146\40\171\157\165\40\141\162\145\40\165\163\151\156\147\40\166\162\157\157\155\54\40\160\154\145\141\163\145\40\143\151\164\145\72\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\122\56\40\102\165\162\141\163\40\141\156\144\40\102\56\40\115\141\171\145\162\40\50\62\60\61\61\51\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\n");fprintf(stderr,"\52\52\52\40\126\122\117\117\115\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);}break;case 0:
sample->vroom=0;sample->escape_eps_ddis_upf=0.;sample->ntupelLE=0;sample->
startCP=0;sample->LEperCP=0;sample->RIS_MS=0;sample->splitter=0;sample->
use_p_norm=0;sample->split_max=1e14;sample->split_min=0.;sample->n_split_max=
0.;sample->n_split_min=0.;sample->LE_taucrit=9999.;sample->MPdynDDIS=0.;sample
->VIS_CS=0;sample->vroomreflectalways=0;if(!q6)fprintf(stderr,"\166\162\157\157\155\40\157\146\146\41\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);break;default:
fprintf(stderr,"\123\124\117\120\41\40\131\157\165\40\141\162\145\40\165\163\151\156\147\40\115\131\123\124\111\103\163\40\154\157\143\141\154\40\145\163\164\151\155\141\164\145\40\164\145\143\150\156\151\161\165\145\54\40\141\156\144\40\171\157\165\162\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\141\162\145\40\163\160\151\153\171\54\n\40\142\165\164\40\171\157\165\40\150\141\166\145\40\156\157\164\40\163\160\145\143\151\146\151\145\144\40\167\150\145\164\150\145\162\40\171\157\165\40\167\141\156\164\40\164\157\40\165\163\145\40\164\150\145\40\166\141\162\151\141\156\143\145\40\162\145\144\165\143\164\151\157\156\40\155\145\164\150\157\144\40\126\122\117\117\115\56\n\40\120\154\145\141\163\145\40\163\160\145\143\151\146\171\40\145\151\164\150\145\162\40\47\155\143\137\166\162\157\157\155\40\157\156\47\40\157\162\40\47\155\143\137\166\162\157\157\155\40\157\146\146\47\40\151\156\40\171\157\165\162\40\151\156\160\165\164\40\146\151\154\145\56\n\40\47\157\156\47\40\151\163\40\162\145\143\157\155\155\145\156\144\145\144\40\146\157\162\40\171\157\165\162\40\143\165\162\162\145\156\164\40\141\160\160\154\151\143\141\164\151\157\156\56\n\105\170\151\164\151\156\147\56\56\56"
);exit(0);}return 0;}int mc_vroom_check_and_verbose(sample_struct*sample,int 
q6,int q7){
#if HAVE_LIDAR
int q50=0;
#endif
if(!(q7)){sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q51){q50=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q50);}
#endif
}if(sample->ntupelLE&&!q6)fprintf(stderr,
"\45\144\55\164\165\160\145\154\40\114\105\n",sample->ntupelLE);if(sample->
splitter){if(!(sample->escape_eps_ddis_upf||sample->LLE_D_DIS)){fprintf(stderr
,"\127\141\162\156\151\156\147\41\40\143\141\156\47\164\40\165\163\145\40\163\160\154\151\164\164\145\162\40\151\146\40\164\150\145\162\145\40\151\163\40\156\157\40\104\111\123\41\n"
);fprintf(stderr,"\56\56\56\56\56\56\56\56\164\165\162\156\151\156\147\40\157\146\146\40\163\160\154\151\164\164\145\162\56\56\56\56\56\56\56\56\56\56\n"
);sample->splitter=0;}else if(!q6)fprintf(stderr,
"\163\160\154\151\164\164\151\156\147\40\141\142\157\166\145\40\45\145\40\n",
sample->split_max);if(!q6)fprintf(stderr,"\155\141\170\40\163\160\154\151\164\164\151\156\147\40\142\145\154\157\167\40\45\145\40\n"
,sample->n_split_max);}if(sample->LidarLocEst){if(sample->escape){fprintf(
stderr,"\41\41\41\40\123\164\157\160\41\41\41\41\40\131\157\165\40\167\141\156\164\40\155\145\40\164\157\40\145\163\164\151\155\141\164\145\40\145\163\143\141\160\145\40\162\141\144\151\141\156\143\145\163\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\141\156\144\40\154\157\143\141\154\40\145\163\164\151\155\141\164\157\162\40\141\164\40\164\150\145\40\163\141\155\145\40\164\151\155\145\41\41\41\40\124\150\151\163\40\151\163\40\40\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\156\157\164\40\147\157\151\156\147\40\164\157\40\167\157\162\153\41\41\41\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\41\41\41\n"
);return-1;}}if(!sample->escape&&sample->vroom){fprintf(stderr,"\41\41\41\40\105\162\162\157\162\41\41\41\40\126\122\117\117\115\40\151\163\40\157\156\40\142\165\164\40\105\123\103\101\120\105\40\151\163\40\157\146\146\41\40\n"
);fprintf(stderr,"\41\41\41\40\123\157\155\145\164\150\151\156\147\40\151\163\40\167\162\157\156\147\41\40\103\157\156\164\141\143\164\40\164\150\145\40\144\145\166\145\154\157\160\145\162\163\40\50\143\157\144\145\40\122\102\51\41\n"
);return-1;}if(!q6&&sample->LidarLocEst){fprintf(stderr,"\40\56\56\56\40\162\165\156\156\151\156\147\40\122\125\114\105\123\40\50\114\151\144\141\162\40\105\155\165\154\141\164\157\162\40\45\144\51\40\n"
,sample->LidarLocEst);if(sample->LLE_D_DIS){fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\104\145\164\145\143\164\157\162\40\104\151\162\145\143\164\151\157\156\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\n"
);fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40");if(
sample->LLE_eps_ddis_upf)fprintf(stderr,"\40\165\163\151\156\147\40\120\150\141\163\145\40\106\165\156\143\164\151\157\156\40\50\45\145\51"
,sample->LLE_eps_ddis_upf);if(sample->LLE_eps_ddis_uda)fprintf(stderr,"\40\165\163\151\156\147\40\104\145\164\145\143\164\157\162\40\101\162\145\141\40\50\45\145\51"
,sample->LLE_eps_ddis_uda);if(sample->LLE_eps_fod_dis_phi)fprintf(stderr,"\40\165\163\151\156\147\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\157\146\40\104\145\164\145\143\164\157\162\40\50\45\145\51"
,sample->LLE_eps_fod_dis_phi);fprintf(stderr,"\n");if(sample->LLE_VIS_FOD)
fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\117\146\40\104\145\164\145\143\164\157\162\n"
);}if(sample->LLE_VIS_QIDD)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\121\165\141\144\162\141\164\151\143\40\111\156\166\145\162\163\145\40\104\145\164\145\143\164\157\162\40\104\151\163\164\141\156\143\145\n"
);if(sample->LLE_RIS_MAS)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\122\145\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\115\157\154\145\143\165\154\141\162\40\141\156\144\40\101\145\162\157\163\157\154\40\123\143\141\164\164\145\162\151\156\147\n"
);if(sample->LLE_sponti)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\123\160\157\156\164\151\163\160\154\151\164\40\166\145\162\163\151\157\156\40\45\144\n"
,sample->LLE_sponti);
#if HAVE_LIDAR
if(sample->LLE_channels==LIDAR_CHANNEL_RAMAN)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\122\141\155\141\156\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);if(sample->LLE_channels==LIDAR_CHANNEL_HSRL)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\110\123\122\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);
#endif
if(sample->LLE_turnmax)fprintf(stderr,
"\40\40\40\40\40\124\165\162\156\155\141\170\151\156\147\56\56\56\n");if(
sample->LLE_taumax)fprintf(stderr,
"\40\40\40\40\40\124\141\165\155\141\170\151\156\147\40\45\145\56\56\56\n",
sample->LLE_taumax);if(sample->LE_taucrit)fprintf(stderr,
"	\40\40\40\111\167\141\142\165\143\150\151\156\147\40\45\145\56\56\56\n",
sample->LE_taucrit);fprintf(stderr,"\n");}return 0;}int mc_vroom_prepare(
sample_struct*sample,atmosphere_struct*q8,float q9,int q10,int q11,int q6){int
 q52=0,q53=0,q54=0,q55=0,kc=0,q56=0,q57=0,q58=0,q59=0;double q60=0.0,q61=0.0,
q62=0.0;double q63=0.0;int q50=0;int nphamat=1;pft**q64=NULL;int*q65=NULL;int*
*q66=NULL;float***q67=NULL,***q68=NULL;double***q69=NULL;double*q70=NULL;float
**q71=NULL,**q72=NULL;double**q73=NULL;int*q74=NULL;int q29=0,q75=0,q76=0;
double q77=0.0;double*F=NULL;double q78=0.0;if(!q6)fprintf(stderr,"\52\52\104\104\111\123\72\40\104\145\146\151\156\151\156\147\40\157\160\164\151\155\141\154\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\40\56\56\56\n"
);sample->n_phase_max=0;for(q29=1;q29<=q8->n_caoth;q29++){switch(q8->
scatter_type[q29]){case MCSCAT_MOL:case MCSCAT_SKIP:break;case MCSCAT_AER:
sample->n_phase_max+=q10;break;case MCSCAT_HG1:case MCSCAT_HG2:(sample->
n_phase_max)++;break;case MCSCAT_PFT:sample->n_phase_max+=q8->phase[q29]->n;
break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q8->scatter_type[q29]);return-1;}}if(sample->n_phase_max!=0){sample->
n_phase_max+=3;q64=calloc(sample->n_phase_max,sizeof(pft*));q65=calloc(sample
->n_phase_max,sizeof(int));q70=calloc(sample->n_phase_max,sizeof(double));q56=
0;for(q29=1;q29<=q8->n_caoth;q29++){switch(q8->scatter_type[q29]){case 
MCSCAT_MOL:q77=0.0;q64[q56]=calloc(1,sizeof(pft));q65[q56]=1;q50=
create_iphase_from_HG(q64[q56++],q77,q6);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q50);if(!q6)fprintf(stderr,"\52\52\104\104\111\123\40\155\157\154\145\143\165\154\141\162\40\144\165\155\155\171\72\40\165\163\151\156\147\40\147\75\60\40\141\163\40\151\163\157\164\162\157\160\151\143\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\n"
);break;case MCSCAT_AER:for(q58=0;q58<q10;q58++)if(q8->phase_aer[q58].nphamat
!=0)q64[q56++]=&(q8->phase_aer[q58]);if(!q6)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\117\120\101\103\40\141\145\162\157\163\157\154\40\154\141\171\145\162\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q8->caoth_name[q29]);break;case MCSCAT_HG1:case MCSCAT_HG2:q77=0.0;for(q58=0;
q58<q10;q58++){if(q8->threed[q29][q58]>=1)for(q75=0;q75<q8->Nx;q75++)for(q76=0
;q76<q8->Ny;q76++){if(q77<q8->g1_3D->prof[q29][q58][q75][q76])q77=q8->g1_3D->
prof[q29][q58][q75][q76];if(q77<q8->g2_3D->prof[q29][q58][q75][q76])q77=q8->
g2_3D->prof[q29][q58][q75][q76];}else{if(q77<q8->g1->prof[q29][q58])q77=q8->g1
->prof[q29][q58];if(q77<q8->g2->prof[q29][q58])q77=q8->g2->prof[q29][q58];}}
q64[q56]=calloc(1,sizeof(pft));q65[q56]=1;q50=create_iphase_from_HG(q64[q56++]
,q77,q6);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q50);if(!q6)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\155\141\170\40\162\137\145\146\146\40\110\107\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q8->caoth_name[q29]);break;case MCSCAT_PFT:for(q58=0;q58<q8->phase[q29]->n;
q58++)q64[q56++]=q8->phase[q29]->iphase[q58];if(!q6)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\162\137\145\146\146\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q8->caoth_name[q29]);break;case MCSCAT_SKIP:if(!q6)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\163\153\151\160\160\151\156\147\40\50\144\165\155\155\171\51\56\56\56\n"
,q8->caoth_name[q29]);break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q8->scatter_type[q29]);return-1;}}if(q56>sample->n_phase_max){fprintf(stderr,
"\45\163\45\144\45\163\45\144\45\163","\105\162\162\157\162\41\40\124\150\145\40\163\145\164\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,sample->n_phase_max,"\51\40\163\155\141\154\154\145\162\n\40\164\150\141\156\40\164\150\145\40\156\165\155\142\145\162\40\157\146\40\144\145\146\151\156\145\144\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,q56,"\51\41\n\40\103\157\156\164\141\143\164\40\122\157\142\145\162\164\40\140\115\145\163\163\171\140\40\102\165\162\141\163\40\146\157\162\40\143\157\155\160\154\141\151\156\164\56\40\105\170\151\164\151\156\147\56\56\56\n"
);return-1;}sample->n_phase_max=q56;if(!q6)fprintf(stderr,"\52\52\104\104\111\123\72\40\146\151\156\141\154\154\171\40\165\163\145\144\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\72\40\45\144\n"
,sample->n_phase_max);q74=calloc(nphamat,sizeof(int));q72=calloc(nphamat,
sizeof(float*));q73=calloc(nphamat,sizeof(double*));q71=calloc(nphamat,sizeof(
float*));q66=calloc(nphamat,sizeof(int*));q67=calloc(nphamat,sizeof(float**));
q69=calloc(nphamat,sizeof(double**));q68=calloc(nphamat,sizeof(float**));for(
q59=0;q59<nphamat;q59++){q66[q59]=calloc(sample->n_phase_max,sizeof(int));q67[
q59]=calloc(sample->n_phase_max,sizeof(float*));q69[q59]=calloc(sample->
n_phase_max,sizeof(double*));q68[q59]=calloc(sample->n_phase_max,sizeof(float*
));for(q56=0;q56<sample->n_phase_max;q56++){q66[q59][q56]=q64[q56]->n[q59];q69
[q59][q56]=q64[q56]->mu[q59];q68[q59][q56]=calloc(q66[q59][q56],sizeof(float))
;for(q58=0;q58<q66[q59][q56];q58++)q68[q59][q56][q58]=(float)q64[q56]->p[
MCSC_MODE_NORMAL][q59][q58];q67[q59][q56]=calloc(q66[q59][q56],sizeof(double))
;for(q58=0;q58<q66[q59][q56];q58++)q67[q59][q56][q58]=acos(q69[q59][q56][q58])
;}}q50=sort_and_add_weighted_phase(sample->n_phase_max,q70,q66,q67,q69,q68,&(
q74),&(q72),&(q73),&(q71),nphamat,1,q6);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q50);sample->phase_max=calloc(1,sizeof(pft*));sample->phase_max[0]=calloc(1,
sizeof(pft));q50=calc_cumulative_table(q73,q71,q74,nphamat,-1.0,sample->
phase_max[0],q9,q6);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\143\165\155\165\154\141\164\151\166\145\137\164\141\142\154\145\50\51\n"
,q50);if(q8->nscaDS>1){for(q56=0;q56<sample->n_phase_max;q56++)for(q58=0;q58<
q66[0][q56];q58++)q68[0][q56][q58]=(float)q64[q56]->p[MCSC_MODE_DELTA_SCALE*(
q64[q56]->nscales>1)][0][q58];q50=sort_and_add_weighted_phase(sample->
n_phase_max,q70,q66,q67,q69,q68,&(q74),&(q72),&(q73),&(q71),1,1,q6);if(q50!=0)
return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q50);for(q57=0;q57<q74[0];q57++)if(q73[0][q57]>q9)break;for(q58=q57+1;q58<q74
[0];q58++)if(q71[0][q58]<q71[0][q57])q71[0][q58]=q71[0][q57];F=calloc(q74[0],
sizeof(double));normalize_phase(q73,q71,F,q74,nphamat,!q6);for(q56=0;q56<
sample->phase_max[0]->n[0];q56++){sample->phase_max[0]->p[
MCSC_MODE_DELTA_SCALE][0][q56]=q71[0][q56];sample->phase_max[0]->F[
MCSC_MODE_DELTA_SCALE][q56]=F[q57];}calc_iphase_coeffs(sample->phase_max[0],
MCSC_MODE_DELTA_SCALE);sample->phase_max[0]->dscale=-999.0;free(F);}for(q56=0;
q56<sample->n_phase_max;q56++)if(q65[q56]){free_iphase(q64[q56]);free(q64[q56]
);}free(q64);free(q65);free(q70);for(q57=0;q57<nphamat;q57++){free(q72[q57]);
free(q73[q57]);free(q71[q57]);}free(q74);free(q72);free(q73);free(q71);for(q59
=0;q59<nphamat;q59++){for(q58=0;q58<sample->n_phase_max;q58++){free(q68[q59][
q58]);free(q67[q59][q58]);}free(q66[q59]);free(q67[q59]);free(q69[q59]);free(
q68[q59]);}free(q66);free(q67);free(q69);free(q68);sample->n_phase_max=1;q78=
0.0;for(q56=0;q56<sample->phase_max[0]->n[0];q56++){if(q78<sample->phase_max[0
]->p[MCSC_MODE_NORMAL][0][q56])q78=sample->phase_max[0]->p[MCSC_MODE_NORMAL][0
][q56];if(q78<1./sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q56])q78=1./
sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q56];}if(q78<q1){sample->
n_phase_max=0;free(sample->phase_max[0]);free(sample->phase_max);if(!q6)
fprintf(stderr,"\52\52\116\157\40\156\145\145\144\40\146\157\162\40\104\104\111\123\54\40\164\165\162\156\151\156\147\40\157\146\146\40\126\122\117\117\115\40\50\151\156\40\143\141\163\145\40\151\164\40\167\141\163\40\157\156\51\n"
);}}if(sample->n_phase_max!=0){q50=calloc_hybrid3D_field(&(q8->spiky_box),q8);
if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\154\157\143\137\150\171\142\162\151\144\63\104\137\146\151\145\154\144\50\51\n"
,q50);for(kc=0;kc<q8->Nz;kc++){if(q8->threed[MCCAOTH_TOT][kc]>=1){for(q29=0;
q29<q8->n_caoth;q29++){if(q8->threed[q29][kc]>=1){for(q52=0;q52<q8->Nx;q52++)
for(q53=0;q53<q8->Ny;q53++)if(q8->ksca3D[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]
->prof[q29][kc][q52][q53]>0.0)q8->spiky_box[kc][q52][q53]=1;}else if(q8->ksca[
MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[q29][kc]>0.0)for(q52=0;q52<q8->Nx;
q52++)for(q53=0;q53<q8->Ny;q53++)q8->spiky_box[kc][q52][q53]=1;}}else{for(q29=
0;q29<q8->n_caoth;q29++){if(q8->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->
prof[q29][kc]>0.0)q8->spiky_box[kc][0][0]=1;}}}}else{sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q51){q50=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q50);}
#endif
}
#ifdef NEWVISQIDD
if(sample->LLE_VIS_QIDD){q60=10.;sample->visqidd_betamax=-q60/sample->lidar[
sample->ili].z_det;sample->visqidd_rmax=-sample->lidar[sample->ili].z_det;
sample->visqidd_facs=sqrt(sample->visqidd_betamax*sample->visqidd_rmax*sample
->visqidd_rmax);fprintf(stderr,"\126\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->visqidd_betamax,sample->visqidd_rmax,sample->visqidd_facs);
#ifdef NEWRISQIDD
q60=1.;sample->risqidd_betamax=-q60/sample->lidar[sample->ili].z_det;sample->
risqidd_facs=sample->risqidd_betamax*sample->visqidd_rmax*sample->visqidd_rmax
;fprintf(stderr,"\122\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->risqidd_betamax,sample->visqidd_rmax,sample->risqidd_facs);
#endif
}
#else
if(sample->LLE_VIS_QIDD){q60=0000.0;q61=-q60/sample->lidar[sample->ili].z_det;
for(kc=0;kc<q8->Nz;kc++){q62=(0.5*(q8->Z[kc]+q8->Z[kc+1])-sample->lidar[0].x[2
])/(-sample->lidar[0].dir.dx[2]*sample->lidar[sample->ili].z_det);if(q62>0.0){
if(q62<1.0)q63=q61;else q63=q61/(q62*q62);if(q8->threed[MCCAOTH_TOT][kc]>=1){
for(q52=0;q52<q8->Nx;q52++)for(q53=0;q53<q8->Ny;q53++)for(q54=0;q54<q8->nscaDS
;q54++)for(q55=0;q55<q8->nscaRIS;q55++)if(q63>q8->kext3D[q54][q55][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc][q52][q53])q8->kext3D[q54][q55][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc][q52][q53]=q63;}else{for(q54=0;q54<q8->
nscaDS;q54++)for(q55=0;q55<q8->nscaRIS;q55++)if(q63>q8->kext[q54][q55][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc])q8->kext[q54][q55][MCVIS_MODE_QIDD]->
prof[MCCAOTH_TOT][kc]=q63;}}}}
#endif
return 0;}int mc_vroom_cloning(photon_struct*p,atmosphere_struct*q8,
sample_struct*sample,result_struct*q12,elevation_struct*q13,albedo_struct*q14,
surftemp_struct*q15,int*q16,int*q17,int q18,int q19,float*q20,float*q21,float*
q22,int q6){int q50=0;photon_struct*q79=NULL;if(sample->ntupelLE>1&&!p->
isclone&&p->scattercounter>=sample->startCP&&p->escapescattercounter+sample->
LEperCP<p->scattercounter+sample->ntupelLE){q79=calloc_photon(sample,q8->Nz,*
q17,q8->nlambda_abs,q8->Nc,q8->n_caoth);
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\154\157\156\151\156\147\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode);
#endif
cp_photon_struct(q79,p,sample,q8->n_caoth);q79->isclone=1;q79->wtree=p->wtree;
q79->photon_status=MCSTATUS_SPLIT;q50=photon_journey(q79,q8,sample,q12,q13,q14
,q15,q16,q17,q18,q19,q20,q21,q22,q6);if(q50<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\143\154\157\156\145\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q50,p->photoncounter);return-1;}destroy_photon(q79,q8->n_caoth);p->
escapescattercounter=p->scattercounter+sample->ntupelLE-1;}return 0;}int 
mc_vroom_splitting_and_rr(photon_struct*p,atmosphere_struct*q8,sample_struct*
sample,result_struct*q12,elevation_struct*q13,albedo_struct*q14,
surftemp_struct*q15,int*q16,int*q17,int q18,int q19,float*q20,float*q21,float*
q22,int q6){double q23=0.0,q82=0.0,q83=0.0;int q50=0,q84=0,q85=0;
#if HAVE_LIDAR
locest_struct lest;
#endif
photon_struct*q86=NULL;double n_split_max=0.0,n_split_min=0.0;int q87=0;if(
sample->LLE_sponti){if(p->scattercounter<7&&sample->LLE_sponti==1)p->
special_weight*=1.5;else p->special_weight*=1.0+0.1/(p->scattercounter+1);}if(
sample->splitter){if(sample->escape)v_mult_mu(p->dir.dx,sample->rad[0].dir.dx,
&q23);
#if HAVE_LIDAR
if(sample->LidarLocEst){q50=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x,&lest);if(q50<0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q50);v_mult_mu(p->dir.dx,lest.dir.dx,&q23);}
#endif
q82=get_phase_max(sample->phase_max,sample->n_phase_max,q23,p->DDIS_SC_mode);
if(sample->use_p_norm){q82/=p->p_norm;if(p->p_norm>1.0){if(q82<1.0)p->p_norm*=
q82;if(p->p_norm<1.0)p->p_norm=1.0;}}q83=p->special_weight*p->weight*p->stokes
[0]*exp(-p->tauris)*q82;
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\162\145\163\164\40\45\145\40\167\145\151\147\150\164\40\45\145\40\160\155\141\170\40\45\145\40\160\156\157\162\155\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->special_weight*p->stokes[0]*exp(-p->tauris),p->weight,q82*p->p_norm,p->
p_norm);
#endif
n_split_max=sample->n_split_max;n_split_min=sample->n_split_min;
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q83,n_split_max);
#endif
if(q83>sample->split_max&&!(sample->ntupelLE&&!p->isclone&&p->scattercounter>=
sample->startCP)){
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q83,n_split_max);
#endif
q83=q83>n_split_max?n_split_max:q83;q84=(int)q83;p->weight/=(double)q84;q86=
calloc_photon(sample,q8->Nz,*q17,q8->nlambda_abs,q8->Nc,q8->n_caoth);
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\145\163\164\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q83,q84);
#endif
if(q84>1000){p->spikewarningcounter++;if(p->spikewarningcounter>1){fprintf(
stderr,"\127\141\162\156\151\156\147\41\40\123\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\40\141\164\40\143\154\157\156\145\40\163\143\141\164\164\145\162\40\157\162\144\145\162\40\45\144\40\163\143\141\164\164\145\162\40\45\144\40\160\150\157\164\157\156\40\45\144\n"
,p->spikewarningcounter,p->clonescattercounter,p->scattercounter,p->
photoncounter);fprintf(stderr,"\40\161\137\163\160\40\45\145\40\167\40\45\145\40\111\60\40\45\145\40\145\170\160\50\55\164\141\165\51\40\45\145\40\120\40\45\145\n"
,p->special_weight,p->weight*(double)q84,p->stokes[0],exp(-p->tauris),q82);}
q87=1;}for(q85=0;q85<q84-1;q85++){cp_photon_struct(q86,p,sample,q8->n_caoth);
q86->wtree=p->wtree;p->photon_status=MCSTATUS_TRAVEL;
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\157\160\171\156\165\155\142\145\162\72\40\156\137\163\160\154\151\164\40\45\144\40\43\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q84,q85);
#endif
q50=photon_journey(q86,q8,sample,q12,q13,q14,q15,q16,q17,q18,q19,q20,q21,q22,
q6);if(q50<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\163\160\154\151\164\40\45\144\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q50,p->photoncounter,q85);return-1;}}if(q87){if(p->spikewarningcounter>1)
fprintf(stderr,"\114\145\141\166\151\156\147\40\163\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\n"
,p->spikewarningcounter);p->spikewarningcounter--;}destroy_photon(q86,q8->
n_caoth);}if(q83<p->weight)q83=p->weight;if(q83<sample->split_min){q83=q83<
n_split_min?n_split_min:q83;if(uvspec_random()>q83)return MCSTATUS_PURGE;p->
weight/=q83;}}return MCSTATUS_DEFAULT;}double get_phase_max(pft**phase_max,int
 n_phase_max,double q23,int SC_mode){int q85=0,q50=0;double q82=0.0,phase=0.0;
for(q85=0;q85<n_phase_max;q85++){q50=get_phase_matrix_pft(phase_max[q85],q23,
SC_mode,1,&phase);if(q50){fct_err_out(q50,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\160\146\164"
,ERROR_POSITION);return-1.0;}q82+=phase;}return q82/((double)n_phase_max);}
void cp_locest(locest_struct*q24,locest_struct*q17,sample_struct*sample,int 
n_caoth){int q85=0,q29=0,q59=0;cp_direction(&(q24->dir),&(q17->dir));q24->
cosalpha=q17->cosalpha;q24->pdir=q17->pdir;q24->pdir_iso=q17->pdir_iso;if(q24
->pdir_sct!=NULL)for(q29=0;q29<n_caoth+1;q29++)for(q59=0;q59<sample->nstokes;
q59++)q24->pdir_sct[q29][q59]=q17->pdir_sct[q29][q59];q24->dist=q17->dist;q24
->distinv=q17->distinv;q24->r_det=q17->r_det;q24->z_det=q17->z_det;for(q85=0;
q85<3;q85++)q24->x_cc[q85]=q17->x_cc[q85];q24->t_det=q17->t_det;for(q85=0;q85<
3;q85++)q24->x_hit[q85]=q17->x_hit[q85];q24->weight_hit=q17->weight_hit;q24->
in_cone=q17->in_cone;q24->will_hit_cone=q17->will_hit_cone;q24->
behind_detector=q17->behind_detector;q24->will_hit_det_plane=q17->
will_hit_det_plane;q24->lidar_outside_grid=q17->lidar_outside_grid;q24->
hit_det_plane_step=q17->hit_det_plane_step;q24->vis_fod_step=q17->vis_fod_step
;q24->vis_fod_step2=q17->vis_fod_step2;q24->vis_fod_kext=q17->vis_fod_kext;for
(q85=0;q85<3;q85++)q24->hitpoint[q85]=q17->hitpoint[q85];}static inline void 
q40(scadis_struct*q28,int q41){int q85=0;q28->q2=0.0;q28->q3=0.0;q28->q4=2.0;
for(q85=0;q85<3;q85++)q28->dirold_dx[q85]=0.0;q28->d_phi=0.0;q28->epsfac=1.0;
q28->mu_max=1.0;q28->mu_min=-1.0;for(q85=0;q85<q41;q85++)q28->q5[q85]=2.0;for(
q85=0;q85<q41;q85++)q28->F_min[q85]=0.0;}int mc_vroom_prep_DDIS(sample_struct*
sample,photon_struct*p,atmosphere_struct*q8,int*q25,int*q26,int*q27,
locest_struct*lest,scadis_struct*q28){int q85=0;
#if HAVE_LIDAR
int q50=0;double q88=0.0,q89=0.0,q90=0.0;double q91=0.0,q92=0.0;double q93=0.0
;double q94=0.0;
#endif
q40(q28,sample->n_phase_max);if(sample->escape_eps_ddis_upf){if(sample->
ntupelLE&&!p->isclone&&q28->epsfac){if(p->scattercounter<sample->startCP)q28->
epsfac=sample->MPdynDDIS/sample->escape_eps_ddis_upf;else q28->epsfac=0.0;}if(
uvspec_random()<sample->escape_eps_ddis_upf*q28->epsfac)*q25=MCDDIS_UPF;for(
q85=0;q85<3;q85++)lest->dir.dx[q85]=sample->rad[0].dir.dx[q85];}
#if HAVE_LIDAR
if(sample->LLE_D_DIS){q50=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x_cc,lest);if(q50!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q50);q28->q2=-lest->dist*lest->cosalpha;q28->q3=lest->dist*sqrt(1.-lest->
cosalpha*lest->cosalpha);q88=1./(q28->q2-sample->lidar[sample->ili].z_det);if(
!p->lest.behind_detector){q91=-(q28->q3+sample->lidar[sample->ili].q95)*q88;
q92=acos(lest->cosalpha);q28->q4=cos(atan(q91)-q92);}if(q8->nthreed==0){q93=(p
->x[2]-q8->Z[p->kc])/(q8->Z[p->kc+1]-q8->Z[p->kc]);q28->epsfac*=q93*q8->q96[p
->kc+1]+(1.-q93)*q8->q96[p->kc];if(q28->epsfac<0.)q28->epsfac=0.0;}if(sample->
lidar[sample->ili].cosalpha[0]<lest->cosalpha){if(q28->q2<sample->lidar[sample
->ili].z_det){if(p->lest.behind_detector){fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\142\145\150\151\156\144\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\151\156\40\146\162\157\156\164\41\n"
);fprintf(stderr,
"\154\157\143\141\164\151\157\156\40\45\145\40\45\145\40\45\145\40\n",p->x[0],
p->x[1],p->x[2]);return-1;}q94=uvspec_random();if(q94<(sample->
LLE_eps_ddis_upf+sample->LLE_eps_ddis_uda)*q28->epsfac)*q25=MCDDIS_UPF;if(q94<
sample->LLE_eps_ddis_uda*q28->epsfac)*q25=MCDDIS_UDA;}else{if(!p->lest.
behind_detector)fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\151\156\40\146\162\157\156\164\40\157\146\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\142\145\150\151\156\144\41\n"
);*q27=-1;}}else{if(sample->lidar[sample->ili].cosalpha[0]<-lest->cosalpha){*
q27=-1;}else{*q27=1;q94=uvspec_random();if(q94<(sample->LLE_eps_ddis_upf+
sample->LLE_eps_ddis_uda)*q28->epsfac)*q25=MCDDIS_UPF;if(!p->lest.
behind_detector&&q94<sample->LLE_eps_ddis_uda*q28->epsfac)*q25=MCDDIS_UDA;if(*
q25&&uvspec_random()<sample->LLE_eps_fod_dis_phi)*q26=1;q90=-q28->q3/q28->q2;
if(q88>0.)q89=-(q28->q3+sample->lidar[sample->ili].q95)*q88;else q89=-(q28->q3
-sample->lidar[sample->ili].q95)*q88;q93=atan((q89-q90)/(1.0+q89*q90));q28->
mu_max=cos(q93);if(q93<0.)q28->mu_max=-q28->mu_max;if(q93<0.&&q93>-1e-10){q93=
0.0;q28->mu_max=1.0;}q28->mu_min=(q28->q2*sample->lidar[sample->ili].q97-q28->
q3*sample->lidar[sample->ili].q98)*lest->distinv;if(q28->mu_min<-1.0){if(q28->
mu_min<-1.0-q99){fprintf(stderr,"\105\162\162\157\162\41\40\155\165\137\155\151\156\40\75\40\45\145\40\151\156\40\163\143\141\164\164\145\162\151\156\147\50\51\40\151\163\40\165\156\160\150\171\163\151\143\141\154\41\n"
,q28->mu_min);return-1;}q28->mu_min=-1.0;}
#ifdef NOFODDIS
q28->mu_max=1.0;q28->mu_min=-1.0;
#endif
for(q85=0;q85<sample->n_phase_max;q85++){q28->q5[q85]=q100(sample->phase_max[
q85],q28->mu_max,p->DDIS_SC_mode);q28->F_min[q85]=q100(sample->phase_max[q85],
q28->mu_min,p->DDIS_SC_mode);}}}}if(*q27==1){q28->d_phi=sample->lidar[sample->
ili].q98/q28->q3*lest->dist;if(q28->d_phi>1.0)q28->d_phi=1.0;q28->d_phi=q101(
q28->d_phi);}else q28->d_phi=180.;
#endif
return 0;}int mu_scatter_special(atmosphere_struct*q8,photon_struct*p,pft**
phase_max,int n_phase_max,int q29,double*mu,scadis_struct q28,int q25,int q27)
{int q50=0;switch(q25){case MCDDIS_NONE:q50=mu_scatter(q8,p,q29,mu);if(q50!=0)
return err_out("\105\162\162\157\162\54\40\155\165\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);break;case MCDDIS_UPF:q50=q42(phase_max,n_phase_max,q28,q27,p->
DDIS_SC_mode,mu);if(q50!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);break;case MCDDIS_UDA:q50=q44(q28.q4,mu);if(q50!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q25);return-1;}return 0;}int random_reflection_special(sample_struct*sample,
atmosphere_struct*q8,elevation_struct*q13,photon_struct*p,pft**phase_max,int 
n_phase_max,double*mu,double*phi,double*q30,scadis_struct*q28,locest_struct 
lest,int q31,int q25,int q26,int q27){int q50=0,q85=0;double q102=0.0;double 
q103=0.0;for(q85=0;q85<3;q85++)q28->dirold_dx[q85]=p->dir.dx[q85];switch(q25){
case MCDDIS_NONE:if(q31==1)random_Lambertian_normal(&(p->dir),q30);else 
random_Isotropic_normal(&(p->dir),q30);v_mult_mu(p->dir.dx,q28->dirold_dx,mu);
break;case MCDDIS_UPF:q50=q42(phase_max,n_phase_max,*q28,q27,p->DDIS_SC_mode,
mu);if(q50!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);break;case MCDDIS_UDA:q50=q44(q28->q4,mu);if(q50!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q25);return-1;}switch(q25){case MCDDIS_NONE:break;case MCDDIS_UPF:case 
MCDDIS_UDA:if(q26)*phi=q28->d_phi*(2.0*uvspec_random()-1.0);else*phi=
sc_Isotropic_phi();if(p->scattercounter==0)q102=p->phi0;else q102=0.0;for(q85=
0;q85<3;q85++)p->dir.dx[q85]=lest.dir.dx[q85];new_direction(*mu,*phi-90.,&(p->
dir),q102);v_mult_mu(p->dir.dx,q30,&q103);if(q103<0.0){p->weight=0.0;}break;
default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q25);return-1;}return 0;}static int q42(pft**phase_max,int n_phase_max,
scadis_struct q28,int q27,int q43,double*mu){int iphase=0;iphase=(int)(
uvspec_random()*((double)n_phase_max-1e-11));switch(q28.mu_min!=-1.0||q28.
mu_max!=1.0){case 0:*mu=sc_mu(phase_max[iphase],1,q43,0.,0.);break;case 1:*mu=
sc_mu(phase_max[iphase],1,q43,q28.q5[iphase],q28.F_min[iphase]);break;default:
fprintf(stderr,"\105\122\122\117\122\41\40\163\157\155\145\164\150\151\156\147\40\166\145\162\171\40\163\164\162\141\156\147\145\40\150\141\160\160\145\156\145\144\40\151\156\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\41\n"
);return-1;}return 0;}static int inline q44(double q45,double*mu){double q104=
0.0,q105=0.0,q94=0.0;q94=2.*uvspec_random();q104=q45*q45;q105=2./(1.-q104);if(
q94<q105*(q45-q104)){if(q94==0.)*mu=0.;else*mu=1./(1.+(1.-q45)*(1.-q45)*q105/
q94);}else*mu=q104+q94/q105;return 0;}int mc_vroom_set_mus_and_phis(
sample_struct*sample,photon_struct*p,int q25,locest_struct lest,scadis_struct 
q28,double*mu,double*q23,double*q37,double phi,double*q38){
#if HAVE_LIDAR
int q50=0;
#endif
switch(q25){case MCDDIS_NONE:v_mult_mu(lest.dir.dx,p->dir.dx,q23);
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){*q37=sqrt(1.0-*q23**q23);q50=q106(p->dir.dx,lest.dir.
dx,*q37,q38);if(q50)return err_out("\105\122\122\117\122\41\40\144\145\162\151\166\145\137\143\160\150\151\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);}
#endif
break;case MCDDIS_UPF:case MCDDIS_UDA:*q23=*mu;v_mult_mu(q28.dirold_dx,p->dir.
dx,mu);if(sample->LLE_VIS_FOD){*q37=sqrt(1.0-*q23**q23);*q38=cosd(phi);}break;
default:fprintf(stderr,"\105\122\122\117\122\41\40\104\104\111\123\163\151\156\147\40\147\151\166\145\163\40\163\157\155\145\164\150\151\156\147\40\163\164\162\141\156\147\145\41\40\45\144\n"
,q25);return-1;}return 0;}int mc_vroom_scattering_calc_phases_and_jacobians(
sample_struct*sample,photon_struct*p,atmosphere_struct*q8,double q32,int q33,
double*q34,double*q35,double*q36){static double**q107=NULL;int q29=0;int q50=0
;if(q33==1){if(q107!=NULL){for(q29=0;q29<=q8->n_caoth;q29++)free(q107[q29]);
free(q107);q107=NULL;}return 0.0;}if(sample->LLE_jacobian||sample->
abs_jacobian)if(q107==NULL){q107=calloc((size_t)q8->n_caoth+1,sizeof(double*))
;for(q29=0;q29<=q8->n_caoth;q29++)q107[q29]=calloc(1,sizeof(double));}q50=
get_phase_matrix_total(q8,p,q32,1,0,sample->spectral_is,sample->
concentration_is,q8->ris_factor,0,q34,q107,q35,&(p->weight));if(q50)return 
fct_err_out(q50,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\164\157\164\141\154"
,ERROR_POSITION);if(sample->LLE_jacobian)*q36=*q34-q107[MCCAOTH_MOL][0];if(*
q34<=0.0||*q35<=0.0||*q36<0.0){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\156\157\162\155\54\40\120\137\163\160\145\143\54\40\141\156\144\40\120\137\151\163\157\145\156\145\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\156\157\162\155\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\151\163\157\145\156\145\40\45\145\40\n"
,*q34,*q35,*q36);return-1;}if(sample->LLE_jacobian){if((*q34)!=0.0)for(q29=1;
q29<q8->n_caoth+1;q29++)p->q_jacobian[0][q29-1][p->kc]+=q107[q29][0]/(*q34);if
((*q36)!=0.0)for(q29=1;q29<q8->n_caoth+1;q29++)p->q_jacobian[1][q29-1][p->kc]
+=q107[q29][0]/(*q36);}if(sample->abs_jacobian){for(q29=1;q29<q8->n_caoth+1;
q29++)if(*q34!=0.0)p->q_jacobian[0][q29-1][p->kc]+=q107[q29][0]/(*q34);}return
 0;}int mc_vroom_DDIS_weight_and_prep_stuff(sample_struct*sample,
atmosphere_struct*q8,double q23,double q37,double q38,int q25,int q27,double 
q34,double q35,double q36,locest_struct lest,scadis_struct q28,photon_struct*p
){double cosalpha=0.0;double q108=q35;
#if HAVE_LIDAR
int q50=0;
#endif
if((p->RIS_mode!=MCRIS_MODE_NORMAL)||sample->LLE_D_DIS||sample->
escape_eps_ddis_upf||sample->LLE_channels||q8->ris_factor!=1.){if(sample->
escape_eps_ddis_upf)q108=q46(q35,q23,0,0,sample->escape_eps_ddis_upf*q28.
epsfac,sample->phase_max,sample->n_phase_max,0,0,0,q28,p->DDIS_SC_mode);
#if HAVE_LIDAR
if(sample->LLE_D_DIS)q108=q46(q35,q23,q27,p->lest.behind_detector,sample->
LLE_eps_ddis_upf*q28.epsfac,sample->phase_max,sample->n_phase_max,sample->
LLE_eps_fod_dis_phi,q38,sample->LLE_eps_ddis_uda*q28.epsfac,q28,p->
DDIS_SC_mode);
#endif
if(!(q108>0.0)){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\163\160\145\143\137\104\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\163\160\145\143\137\104\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\156\157\162\155\40\45\145\40\n"
,q108,q35,q34);return-1;}p->weight*=q34/q108;
#ifdef MUCHOUT
if(p->q80==q81)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\156\145\167\40\167\145\151\147\150\164\40\45\145\40\146\162\157\155\40\45\145\40\57\40\45\145\40\167\151\164\150\40\155\165\62\40\45\145\40\141\156\144\40\120\163\160\145\143\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->weight,q34,q108,q23,q35);
#endif
p->q_isoene*=q36/q34;}if(sample->LLE_VIS_FOD||sample->LLE_eps_ddis_uda){
v_mult_mu(p->dir.dx,sample->lidar[sample->ili].dir.dx,&cosalpha);p->lest.
hit_det_plane_step=(q28.q2-sample->lidar[sample->ili].z_det)/cosalpha;p->lest.
will_hit_det_plane=(p->lest.hit_det_plane_step>0.0);}
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){q50=q109(q23,q37,q38,q28.q3,q28.q2,lest.distinv,sample
->lidar[sample->ili].t_det,&p->lest);if(q50!=0)return err_out("\105\162\162\157\162\54\40\143\141\154\143\137\144\151\163\164\141\156\143\145\137\164\157\137\143\157\156\145\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q50);}
#endif
return 0;}static double q46(double q34,double q23,int q27,int behind_detector,
double q47,pft**phase_max,int n_phase_max,double q48,double q38,double q49,
scadis_struct q28,int q43){double q35=0.0;double q110=0.0;double q111=0.0;
double q112=0.0;int q85=0;if(behind_detector){q47+=q49;q49=0.0;}if(q49>0.0)if(
q23>=0.){if(q23>q28.q4)q111=2./(1.-q28.q4*q28.q4);else{q111=(1.-q28.q4)/(1.-
q23);q111*=2./(1.-q28.q4*q28.q4)*q111;}}if(q27==-1)return q34;if(q28.mu_min==-
1.0&&q28.mu_max==1.0){if(q47>0.0)q110=get_phase_max(phase_max,n_phase_max,q23,
q43);q112=1.0;}else{if(q23<=q28.mu_max&&q23>=q28.mu_min)q110=get_phase_max(
phase_max,n_phase_max,q23,q43)*2.0/(q28.q5[q85]-q28.F_min[q85]);if(q38>cosd(
q28.d_phi))q112=1.0-q48*(1.0-180.0/q28.d_phi);else q112=1.0-q48;}q35=(1.0-q47-
q49)*q34+q112*(q47*q110+q49*q111);if(!(q35>0.0)){fprintf(stderr,"\105\162\162\157\162\40\151\156\40\120\137\163\160\145\143\72\40\45\145\40\72\40\145\160\163\137\144\144\151\163\137\165\160\146\40\45\145\40\145\160\163\137\144\144\151\163\137\165\144\141\40\45\145\40\120\137\156\157\162\155\40\45\145\40\146\141\143\137\160\150\151\40\45\145\40\120\137\144\144\151\163\137\165\160\146\40\45\145\40\120\137\144\144\151\163\137\165\144\141\40\45\145\40\155\165\62\40\45\145\40\155\165\155\141\170\40\45\145\40\155\165\155\151\156\40\45\145\n"
,q35,q47,q49,q34,q112,q110,q111,q23,q28.mu_max,q28.mu_min);fprintf(stderr,"\156\137\160\150\141\163\145\137\155\141\170\40\45\144\40\106\137\155\141\170\40\45\145\40\106\137\155\151\156\40\45\145\40\155\165\137\165\144\141\137\144\142\40\45\145\n"
,n_phase_max,q28.q5[0],q28.F_min[0],q28.q4);}
#ifdef MUCHOUT
if(q81==110)fprintf(stderr,"\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\40\143\141\154\143\137\160\137\163\160\145\143\72\40\45\145\40\45\145\40\45\144\40\n"
,q28.mu_min,q28.mu_max,q23<q28.mu_max);
#endif
return q35;}int calloc_hybrid3D_field(int****q39,atmosphere_struct*q8){int kc=
0,q52=0;*q39=calloc((size_t)q8->Nz,sizeof(int**));if(*q39==NULL){fprintf(
stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(kc=0;kc<q8->Nz;kc++){if(q8->threed[MCCAOTH_TOT][kc]>=1){(*q39)
[kc]=calloc((size_t)q8->Nx,sizeof(int*));if((*q39)[kc]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(q52=0;q52<q8->Nx;q52++){(*q39)[kc][q52]=calloc((size_t)q8->Ny,
sizeof(int));if((*q39)[kc][q52]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}else{(*q39)[kc]=calloc((size_t)1,sizeof(int*));if((*q39)[kc]==
NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}q52=0;(*q39)[kc][q52]=calloc((size_t)1,sizeof(int));if((*q39)[kc][
q52]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}return 0;}void free_hybrid3D_field(int****q39,atmosphere_struct*
q8){int kc=0,q52=0;for(kc=0;kc<q8->Nz;kc++){if(q8->threed[MCCAOTH_TOT][kc]>=1)
{for(q52=0;q52<q8->Nx;q52++)free((*q39)[kc][q52]);free((*q39)[kc]);}else{free(
(*q39)[kc][0]);free((*q39)[kc]);}}free(*q39);}
