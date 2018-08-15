/*--------------------------------------------------------------------
 * $Id: mie_lex.l 3194 2015-11-25 18:19:46Z svn-kylling $
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

%{

#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "./mie.h"

/* following line needed to suppress compiler warning about unused function "input" */
#define YY_NO_INPUT

int mie_check (input_struct input);

float  yytext2float(char* yytext, int which_token, int max_token);
float* yytext2floats(char* yytext, int nary);
int    yytext2nfloats(char* yytext);
int    yytext2int(char* yytext, int which_token, int max_token);
int*   yytext2integers(char* yytext, int nary);
int    yytext2nintegers(char* yytext);
char*  yytext2string(char* yytext, int which_token, int max_token);
int    yytext2ntokens(char* yytext);
int    yytextCstring(char *dest, char* yytext, int which_token, int max_token);
double yytext2double(char* yytext, int which_token, int max_token);


char* s, *t;
float* tmp;
int nary;
int ierror=0;
int line_number;
input_struct Input;
output_struct Output;
int ntokens=0;
int i=0;
int tmpint;

%}

%option noyywrap
%option nounput

DIGIT    [0-9]
SPACE    [ \t]
CCHAR    [#]
COMMENT  ({CCHAR}.*)
INTEGER  {DIGIT}+
SPACE_INTEGER {SPACE}+{INTEGER}
INTEGERS {SPACE_INTEGER}+ 
EXPON   ([DdEe][+-]?{DIGIT}+)
FLOAT    ((({DIGIT}+"."{DIGIT}*)|({DIGIT}+"."{DIGIT}*{EXPON}))|({INTEGER}))
SPACE_FLOAT {SPACE}+{FLOAT}
FLOATS   {SPACE_FLOAT}+ 
SIGN     [-+]
SIGNED_FLOAT {SIGN}?{FLOAT}
SPACE_SIGNED_FLOAT {SPACE}+{SIGNED_FLOAT}
SIGNED_FLOATS {SPACE_SIGNED_FLOAT}+
SIGNED_INTEGER {SIGN}?{INTEGER}
WORD     [_\.\/A-Za-z0-9][_\.\/A-Za-z0-9\-]*
SPACE_WORD {SPACE}+{WORD}
WORDS {SPACE_WORD}*


ZERO_OR_TWO_SIGNED_FLOATS ({SPACE}?{SIGNED_FLOAT}+{SPACE}+{SIGNED_FLOAT})?
ZERO_OR_TWO_OR_WORD {ZERO_OR_TWO_SIGNED_FLOATS}|{WORD}
FLOATS_OR_WORD ({FLOATS}|{SPACE_WORD})

%%

  /* These are the allowed input lines in the mie input file 
     The input line specifications must come before the eating of
     blanks, otherwise things get messed up 
     Obviously this whole lot could have been handled more properly
     using bison etc., but than rather few would know what was going 
     on??? 
   */

accuracy_phase{SPACE}+{SIGNED_FLOAT}  {Input.accuracy = yytext2float(yytext,2,2); }

aerosol_type{SPACE}+{WORD}          {s =  yytext2string(yytext,2,2);
                                     Input.medium = AEROSOL;
                                     if (Input.distribution == DIST_NONE)
                                         Input.distribution = DIST_AER;
                                     if  (strcmp(s,"inso")==0) {Input.aerosol_type = TYPE_INSO;}
				     else if (strcmp(s,"waso")==0) {Input.aerosol_type = TYPE_WASO;}
				     else if (strcmp(s,"soot") ==0) {Input.aerosol_type = TYPE_SOOT;}
				     else if (strcmp(s,"ssam") ==0) {Input.aerosol_type = TYPE_SSAM;}
                                     else if (strcmp(s,"sscm") ==0) {Input.aerosol_type = TYPE_SSCM;}
                                     else if (strcmp(s,"minm") ==0) {Input.aerosol_type = TYPE_MINM;}
                                     else if (strcmp(s,"miam") ==0) {Input.aerosol_type = TYPE_MIAM;}
                                     else if (strcmp(s,"micm") ==0) {Input.aerosol_type = TYPE_MICM;}
                                     else if (strcmp(s,"mitr") ==0) {Input.aerosol_type = TYPE_MITR;}
                                     else if (strcmp(s,"suso") ==0) {Input.aerosol_type = TYPE_SUSO;}
                                     else {fprintf(stderr,"Unknown 'aerosol type %s' on line %d\n", s, line_number+1);
				       ierror++;}
                                     free(s);
                                     }

basename{SPACE}+{WORD}              {
                                     Input.basename=(char *) calloc (strlen(yytext)+1,sizeof(char));
                                     Input.basename =yytext2string(yytext,2,2);
                                    }  

data_files_path{SPACE}+{WORD}       {yytextCstring (Input.data_files_path, yytext, 2, 2); 
                                     /* append slash */
                                     if (Input.data_files_path[strlen(Input.data_files_path)-1] != '/')
				       strcat (Input.data_files_path, "/");
                                    }

distribution{SPACE}+{WORD}{SPACE}+{FLOAT} {ntokens = yytext2ntokens(yytext);
                                     if (ntokens != 2) {
                                       fprintf(stderr, "Wrong number of arguments for 'distribution %s' on line %d\n", t, line_number+1);
				       ierror++;}
                                     s =  yytext2string(yytext,2,ntokens+1);
				     if      (strncmp(s,"none", 4)     ==0||strncmp(s,"NONE", 4)     ==0) {Input.distribution = DIST_NONE;}
				     else if (strncmp(s,"gamma",5)     ==0||strncmp(s,"GAMMA",5)     ==0) {Input.distribution = DIST_GAMMA;}
				     else if (strncmp(s,"lognormal",9) ==0||strncmp(s,"LOGNORMAL",9) ==0) {Input.distribution = DIST_LOGNORMAL;}
				     else if (strncmp(s,"file",4)      ==0||strncmp(s,"FILE",4)      ==0) {Input.distribution = DIST_FILE;}
                                     else {fprintf(stderr,"Unknown 'distribution %s' on line %d\n", s, line_number+1);
				       ierror++;}
                                     free(s);
                                     if (Input.distribution == DIST_GAMMA || Input.distribution == DIST_LOGNORMAL)
                                       Input.alpha = yytext2float(yytext, 3, ntokens+1);
                                     else 
                                       Input.sd_filename = yytext2string(yytext,3,ntokens+1);
                                    }

dx_max{SPACE}+{SIGNED_FLOAT}  {Input.dx_max = yytext2float(yytext,2,2); }

mass_density{SPACE}+{SIGNED_FLOAT}  {Input.rho_user = yytext2float(yytext,2,2); }



mie_program{SPACE}+{WORD}           {s =  yytext2string(yytext,2,2);
				     if      (strncmp("MIEV0",s,5)==0) {Input.program = MIEV0;}
				     else if (strncmp(s,"BH",2)   ==0) {Input.program = BHMIE;}
				     else {fprintf(stderr,"Unknown 'mie_program %s' on line %d\n", t, line_number+1);
				       ierror++;}
                                     free(s);
                                    }
mimcut{SPACE}+{SIGNED_FLOAT}        {Input.mie.mimcut = yytext2float(yytext,2,2); }
nmom{SPACE}+{INTEGER}               {Input.mie.nmom   = yytext2int(yytext,2,2); 
                                     Input.mie.momdim = Input.mie.nmom+1;}

nmom_netcdf{SPACE}+{INTEGER}        Input.nmom_netcdf = yytext2float(yytext,2,2);
n_r_max{SPACE}+{SIGNED_FLOAT}       {Input.n_r_max = yytext2float(yytext,2,2); }
nthetamax{SPACE}+{INTEGER}          {Input.nthetamax = yytext2int(yytext,2,2); } 
nstokes{SPACE}+{INTEGER}            {Input.nstokes   = yytext2int(yytext,2,2); } 
                                     

output_user{WORDS}                  {Input.n_output_user = yytext2ntokens(yytext);
				     Input.output_user = calloc (Input.n_output_user, sizeof (int));
				     for (i=2;i<=Input.n_output_user+1; i++) {
				       t = yytext2string(yytext, i, Input.n_output_user+1);  
                                       if      (strcmp("none",t)      ==0)  {Input.output_user[i-2] = OUTPUT_USER_NONE;}
				       else if (strcmp("lambda",t)==0 || strcmp("wavelength",t)==0)
                                                                            {Input.output_user[i-2] = OUTPUT_USER_WAVE;}
				       else if (strcmp("wavenumber",t) ==0) {Input.output_user[i-2] = OUTPUT_USER_WAVENUMBER;}
				       else if (strcmp("refrac_real",t)==0) {Input.output_user[i-2] = OUTPUT_USER_REFRAC_REAL;}
				       else if (strcmp("refrac_imag",t)==0) {Input.output_user[i-2] = OUTPUT_USER_REFRAC_IMAG;}
                                       else if (strcmp("r_eff"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_R_EFF;}
				       else if (strcmp("qext"       ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_QEXT;}
                                       else if (strcmp("qsca"       ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_QSCA;}
				       else if (strcmp("qback"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_QBACK;}
				       else if (strcmp("omega"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_OMEGA;}
				       else if (strcmp("gg"         ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_GG;}
				       else if (strcmp("spike"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_SPIKE;}
				       else if (strcmp("pmom"       ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_PMOM;}
				       else if (strcmp("sforw"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_SFORW;}
				       else if (strcmp("sback"      ,t)==0) {Input.output_user[i-2] = OUTPUT_USER_SBACK;}
				       else if (strcmp("cloudprp"   ,t)==0) {Input.output_format = OUTPUT_FORMAT_CLOUDPROP;}
				       else if (strcmp("aerosolprp" ,t)==0) {Input.output_format = OUTPUT_FORMAT_AERPROP;}
				       else if (strcmp("netcdf"     ,t)==0) {Input.output_format = OUTPUT_FORMAT_NETCDF;}
				       else {fprintf(stderr,"Unknown 'output_user %s' on line %d\n", t, line_number+1);
				       ierror++;}
				       free(t);
				     }
                                    }
r_eff{SIGNED_FLOATS}+{SPACE_WORD}?|r_mean{SIGNED_FLOATS}+{SPACE_WORD}? {ntokens = yytext2ntokens(yytext);
                                     Input.r_eff_min    = yytext2double(yytext, 2, ntokens+1);
                                     if (ntokens>1)
                                       Input.r_eff_max  = yytext2double(yytext, 3, ntokens+1);
                                     if (ntokens>2)
                                       Input.r_eff_step = yytext2double(yytext, 4, ntokens+1);
                                     if (ntokens>3) {
				       t = yytext2string(yytext, 5, ntokens+1);  
                                       if (strcasecmp ("log",t)==0)  {Input.r_eff_log  = 1;}
                                     }
                                    }

refrac{SPACE}+{WORD}{SPACE}?{ZERO_OR_TWO_OR_WORD}  {
                                       s = (char *) calloc (strlen(yytext)+1,sizeof(char));
                                       strcpy (s, yytext); 
				       t =  yytext2string(s,2,4);
				       if (strncmp("ice",t,3)==0) 
					 Input.medium = ICE;
				       else  {
					 if (strncmp("water",t,5)==0) 
					   Input.medium = WATER;
					 else  {
					   if (strncmp("user",t,4)==0) {
					     Input.medium = USER;
					     
					     strcpy (s, yytext); 
					     nary = yytext2nfloats(s);
					     
					     strcpy (s, yytext);
					     tmp = yytext2floats(s, nary);
					     
					     Input.crefin.re = tmp[1];
					     Input.crefin.im = tmp[2];
					   }
					   else {
					     if (strncmp("file",t,4)==0) {
					       Input.medium = USER;
					       Input.ri_filename = yytext2string(yytext,3,3); 
					     }
                                             else {
                                               fprintf(stderr,"Unknown command on line %d\n", line_number+1);
                                               fprintf(stderr,"%s\n", yytext);
                                               ierror++;
                                             }
					   }
					 }
				       }
				       free(s);}
size_distribution_file{SPACE}+{WORD} {Input.sd_filename = yytext2string(yytext,2,2);
                                      Input.distribution = DIST_FILE; }
temperature{SPACE}+{SIGNED_FLOAT}    {Input.temperature   = yytext2float(yytext,2,2); }

verbose                             {Input.verbose = 1;}

wavelength{FLOATS_OR_WORD}|wvn{FLOATS_OR_WORD} {
                                       ntokens = yytext2ntokens(yytext);
                                       if (ntokens == 1) {
                                         s = yytext2string(yytext,2,ntokens+1);
				         if      (strlen(s)==8 && strncasecmp("kato2.96",s,8)    ==0) {Input.ck_scheme = CK_KATO2_96;}
				         else if (strlen(s)==5 && strncasecmp("kato2",s,5)       ==0) {Input.ck_scheme = CK_KATO2;}
				         else if (strlen(s)==5 && strncasecmp("kato2andwandji",s,14) ==0) {Input.ck_scheme = CK_KATO2ANDWANDJI;}
				         else if (strlen(s)==4 && strncasecmp("kato",s,4)        ==0) {Input.ck_scheme = CK_KATO;}
				         else if (                strncasecmp("fu",s,2)          ==0) {Input.ck_scheme = CK_FU;}
				         else if (                strncasecmp("avhrr_kratz",s,11)==0) {Input.ck_scheme = CK_AVHRR_KRATZ;}
				         else if (                strncasecmp("sbdart",s,6)      ==0 || 
                                                                  strncasecmp("lowtran",s,7)     ==0) {Input.ck_scheme = CK_LOWTRAN;}
				         else                                                         {Input.ck_scheme = CK_GENERIC;
                                           Input.filename_ck_generic = (char *) calloc (strlen(s)+1,sizeof(char));
				           strcpy (Input.filename_ck_generic, s);}
				         free(s);
                                       }
                                       else if (ntokens == 2) {
				         Input.wl.wl_start  = yytext2double(yytext,2,ntokens+1); 
                                         Input.wl.wl_end    = yytext2double(yytext,3,ntokens+1);
                                       }
                                       else {
                                         fprintf(stderr,"Unknown command on line %d\n", line_number+1);
                                         fprintf(stderr,"%s\n", yytext);
                                         ierror++;
                                       } 
				     }
wavelength_step{SPACE}+{SIGNED_FLOAT}|wvn_step{SPACE}+{SIGNED_FLOAT}       {Input.wl.wl_step   = yytext2double(yytext,2,2); }
wavelength_index{SPACE_INTEGER}+{SPACE_INTEGER}|wvn_index{SPACE_INTEGER}+{SPACE_INTEGER}  {s = (char *) calloc (strlen(yytext)+1,sizeof(char));
                                     strcpy (s, yytext);
				     Input.lambda_i_start = yytext2int(s,2,3); 
                                     Input.lambda_i_end   = yytext2int(yytext,3,3);
				     if (Input.lambda_i_start > Input.lambda_i_end) {
				       tmpint = Input.lambda_i_start;
				       Input.lambda_i_start = Input.lambda_i_end;
				       Input.lambda_i_end   = tmpint;
				     }
                                     free(s);}


{SPACE}+       ;         /* eat up blanks */ 

{COMMENT}*     ;         /* eat up comments */

\n                        { ++line_number;};

{WORD}                    {fprintf(stderr,"Unknown command on line %d\n", line_number+1);
                           fprintf(stderr,"%s\n", yytext);
                           ierror++;}

<<EOF>> { 
  if (ierror > 0) 
    exit(1);
  
  yyterminate();
}      

%%

/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void usage()
{
  fprintf (stderr, "\nCalculates scattering and absorption properties using Mie theory.\n");
  fprintf (stderr, "See libRadtran documentation for more information.\n\n");
  fprintf (stderr, "No error checking yet, so be careful.\n\n");

  fprintf (stderr, "Usage: mie [-h] < input_file > Output_file\n");
  fprintf (stderr, " -h    Print this message.\n");
  fprintf (stderr, " -f    Input filename. \n");  
  fprintf (stderr, "\n");  
}

int main(int argc, char **argv) 
{
  
  /*  @12c@                

<lpdoc>

\section{Tool for Mie calculations - \code{mie}}
\label{sec:mie}

The various input parameters of the \code{mie} tool are described in the following.

\begin{description}

\option{accuracy_phase}
Accuracy of the phase function values in the netcdf output. Default value is $0.01$ (1 percent).
\option{aerosol\_type}
With this option Mie calculations are performed for the specified aerosol type. 
\begin{Verbatim}
aerosol_type type
\end{Verbatim}
The 
aerosol properties (refractive index, size distribution, density, humidity) are taken 
from the OPAC database \citep{hess98:_optic_proper_aeros_cloud}
Possible values for \code{type} are
\begin{description}
\parameter{inso }
Water insoluble aerosol consists mostly of soil particles with a certain amount of organic material. 
\parameter{waso}
Water soluble aerosol originates from gas to particle conversion and consists of various types of sulfates, nitrates, and other, also organic water-soluble substances.
\parameter{soot}
Soot is absorbing black carbon, which is not soluble in water. In reality soot particles have a chain-like character, which of course is not accounted for in Mie calculations of optical properties. The optical properties are calculated assuming many very small spherical particles.  
\parameter{ssam }
Sea salt particles consist of the various kinds of salt contained in seawater. The different modes are given to allow for a different wind-speed-dependant increase of particle number for particles of different size.
This aerosol type represents the accumulation mode.
\parameter{sscm }
Sea salt particles (coarse mode).
\parameter{minm  }
Mineral aerosol or desert dust is produced in arid regions. It consists of a mixture of quartz and clay minerals and is modeled with three modes to allow to consider increasing relative amount of large particles for increased turbidity.
This aerosol type represents the nucleation mode.
\parameter{miam}
Mineral aerosol (accumulation mode). 
\parameter{micm}
Mineral aerosol (coarse mode). 
\parameter{mitr  }
Mineral transported is used to describe desert dust that is transported over long distances with a reduced amount of large particles.
\parameter{suso}
The sulfate component is used to describe the amount of sulfate found in the Antarctic aerosol. This component is not suited to describe antropogenic sulfate aerosols that are included in the water-soluble component.  
\end{description}

\option{basename}
Filename for output of Mie program. 
\begin{Verbatim}
basename filebase
\end{Verbatim}
This option is only used in combination 
with \code{output\_user netcdf}. The default is wc. for water, ic. for
ice, or waso., inso. etc. for OPAC aerosols.

\option{distribution}
If specified the effective radius is converted into a size distribution of droplets.
\begin{Verbatim}
distribution distribution_type distribution_parameter
\end{Verbatim}
where distribution type is one of the two following:
\begin{description}
\parameter{GAMMA}
The  Gamma distribuition of cloud droplet sizes is 
\begin{equation}
n(r) = a r^{\alpha} \exp(-b r), 
\end{equation}
where $\alpha$ is the distribution parameter given as second argument. ($a$ and $b$
are determined automatically.) The effective radius of the distribution is 
$r_{\text{eff}} = (\alpha+3)/b$.
A typical value for water clouds is $\alpha=7$. For ice clouds a typical value is $\alpha=1$.
A large value of $\alpha$ gives close to a monodisperse distribution.
\parameter{LOGNORMAL}
The lognormal distribution of cloud droplet sizes is 
\begin{equation}
n(r) = \frac{a}{r} \exp\left(-\frac{1}{2}\left(\frac{\ln(r)-\ln(r_0)}{ln(\sigma)}\right)^2\right),
\end{equation}
where $r_0$ is the logarithmic mode of the distribution (provided by option \code{r\_eff}) 
and $\sigma$ is the standard deviation, which is given by the second argument.
\end{description}

\option{dx\_max}
This option makes sense in combination with \code{distribution GAMMA} or \code{distribution LOGNORMAL}.
It can be used to specify the maximum widths of the size distribution
bins, which are sampled on a size parameter ($\frac{2\pi r}{\lambda}$) grid. The default value 
is 0.03 which is not very accurate for small wavelengths. In order to get accurate phase matrices 
this value should be decreased. 

\option{mass\_density}
Specifies the mass density of the medium.
\begin{Verbatim}
mass_density value
\end{Verbatim}
Useful in combinition with \code{refrac user} and \code{output\_user cloudprop},
as the format of \code{cloudprop} specifies the extinction coefficient per 
mass and not per volume as usual in this mie program.

\option{mie\_program}
Specify which Mie program to use:
\begin{Verbatim}
mie_program type
\end{Verbatim}
where \code{type} is one of
\begin{description}
\parameter{BH}
The Mie scattering program by Bohren and Hoffmann, 
ftp://ftp.astro.princeton.edu/draine/scat/bhmie/bhmie.f
\parameter{MIEV0}
The Mie scattering program by W. Wiscombe. For documentation see 
\file{libsrc\_f/MIEV.doc} and the NCAR Mie report at 
ftp://climate1.gsfc.nasa.gov/wiscombe/Single\_Scatt/Homogen\_Sphere/Exact\_Mie/
\end{description}

\option{mimcut}
(positive) value below which imaginary refractive
index is regarded as zero (computation proceeds
faster for zero imaginary index). Only used by \code{mie\_program} \code{MIEV0}.
\begin{Verbatim}
mimcut value
\end{Verbatim}

\option{nmom }
Number of moments of the phase function to be calculated (default: 0). 
\begin{Verbatim}
nmom value
\end{Verbatim}
Only possible with \code{mie\_program} \code{MIEV0}.

\option{nmom\_netcdf}
Specify the number of Legendre polynomials that are written to the
netcdf file. 
\begin{Verbatim}
nmom_netcdf value
\end{Verbatim}
This option only makes sense if \code{output\_user netcdf}
is specified. If not specified, all polynomials are written.
For the calculation of the phase function all polynomials are of
course considered. 

\option{n\_r\_max}
This option makes sense in combination with \code{distribution GAMMA} or \code{distribution LOGNORMAL}. 
It defines the upper cutoff value for the size distribution in terms of effective radius $r_{\text{eff}}$. The default is 5, 
which means that the size distribution is cut off at a value of  $5\cdot r_{\text{eff}}$. This value should 
be increased if only small $r_{\text{eff}}$ are calculated.

\option{nstokes }
Number of Stokes parameters (default: 1). 
\begin{Verbatim}
nstokes value
\end{Verbatim}
For \code{nstokes=1} the
Legendre polynomials of the phase function will be calculated. To
calculate all phase matrix elements required for polarized radiative
transfer, set \code{nstokes=4}.      

\option{nthetamax }
Specify the maximum number of scattering angles to be used to sample
the phase matrix. 
\begin{Verbatim}
nthetamax value
\end{Verbatim}
The default value is 1000. If the accuracy of the
phase function is less than 1\% for \code{nthetamax} angles a warning
is printed to the screen. The option is only meaningful in combination
with \code{output\_user netcdf}, otherwise phase functions are not
computed. 

\option{output\_user}
The mie output is one line of output quantities to standard output (stdout) for each 
wavelength and each particle radius.
With this option the user may specify the columns desired for output:
\begin{verbatim}
output_user output_1 output_2 ... output_n
\end{verbatim}
where output\_i is one of following arguments:
\begin{description}
\parameter{lambda}
Wavelength in nm.
\parameter{wavenumber}
Wave number in cm-1.
\parameter{r\_eff}
particle radius in micro meter.
\parameter{refrac\_real}
The real part of the refractive index.
\parameter{refrac\_imag}
The imaginary part of the refractive index.
\parameter{qext}
The extinction efficiency factor, if \code{r\_eff} is specified, or
the extinction coefficient [km-1] per unit concentration [cm$^3$/m$^3$],
if a \code{size\_distribution\_file} is specified. If the medium is 
liquid water, 1 cm$^3$/m$^3$ equals a liquid water content of 1g/m$^3$ because 
the density of water is close to 1 g/cm$^3$. For ice and other substances,
the density has to be considered (0.917 g/cm$^3$ for ice at 273K).
\parameter{qsca}
The scattering efficiency factor, if \code{r\_eff} is specified, or
the scattering coefficient [km-1] per unit concentration [cm$^3$/m$^3$],
if a \code{size\_distribution\_file} is specified.
\parameter{qback}
The backscattering efficiency factor, if \code{r\_eff} is specified, or
the backscattering coefficient [km-1] per unit concentration [cm$^3$/m$^3$],
if a \code{size\_distribution\_file} is specified. Works only with \code{mie_program BH}. 
\parameter{omega}
The single scattering albedo.
\parameter{gg}
The asymmetry parameter.
\parameter{sforw}
(Complex) forward-scattering amplitude S1 at 0 degrees.
\parameter{sback}
(Complex) back-scattering amplitude S1 at 180 degrees. 
\parameter{spike}
To quote from Wiscombe's \file{MIEV0.doc}:
     \begin{quotation}
(REAL) magnitude of the smallest denominator of
            either Mie coefficient (a-sub-n or b-sub-n),
            taken over all terms in the Mie series past
            N = size parameter XX.  Values of SPIKE below
            about 0.3 signify a ripple spike, since these
            spikes are produced by abnormally small denominators
            in the Mie coefficients (normal denominators are of
            order unity or higher).  Defaults to 1.0 when not
            on a spike.  Does not identify all resonances
            (we are still working on that).
     \end{quotation}
Meaningless if a \code{size\_distribution\_file} was specified.
\parameter{pmom}
The nmom+1 moments (from 0 to nmom, see option \code{nmom}) of the phase function. 
The phase function 
 $p(\mu)$ is \\*
\begin{equation}
   p (\mu) = \sum_{m=0}^{\infty} (2m+1) \cdot k_m \cdot P_m (\mu) 
\end{equation}
  where $k_m$ is the m'th moment and $P_m (\mu)$ is the m'th Legendre polynomial.
\parameter{cloudprp}
This is a special option which, if specified, must be the only option of \code{output\_user},
as \code{cloudprp} specifies a whole format of the output. In particular this option is 
useful when a correlated-k wavelength grid is specified with \code{wavelength}. 
If specified, the output will be written in a format, which can be directly used by 
{\sl libRadtran}. See also \code{output netcdf}.
See \code{uvspec}-options \code{ic\_properties} and \code{wc\_properties} and 
there the items \code{mie} and \code{filename}.
\parameter{aerosolprp}
This option is similar to the \code{cloudprp} option. The only difference is that 
the effective radius dimension is replaced by humidity values of the aerosol. 
\parameter{netcdf }
This option writes the output to a netCDF 
file which can be used by \code{uvspec} using the options \code{ic\_properties} 
and \code{wc\_properties}.
\end{description}The default output is: 
\begin{verbatim}
lambda refrac_real refrac_imag qext omega gg spike pmom
\end{verbatim}

\option{r\_eff}
The radius [micron] of the particle to calculate single scattering properties of.
Used together with the wavelength information to calculate the 
Mie size parameter.
\begin{verbatim}
r_eff radius
\end{verbatim}
The user can optionally specify a 2nd and 3rd argument to make a loop over several radii:
\begin{verbatim}
r_eff radius_min radius_max radius_step 
\end{verbatim}
Alternatively, a logarithmic grid may be defined by specifying a multiplicative increment
(indicated by the keyword \code{log}):
\begin{verbatim}
r_eff radius_min radius_max radius_increment log
\end{verbatim}
In this case the radius is increased by a factor \code{radius_increment} for each step.

\option{refrac}
Specify which refractive index to use. 
\begin{Verbatim}
refrac type
\end{Verbatim}
The following choices for \code{type} are valid: 
\begin{description}
\parameter{ice}
The complex refractive index is taken from the REFICE function of W. Wiscombe. 
\parameter{water}
The complex refractive index is taken from the REFWAT function of W. Wiscombe. 
\parameter{user }$<$re$>$ $<$im$>$
A user defined refractive index. re and im are the real and imaginary parts
(both positive numbers).
\parameter{file }$<$filename$>$
Read refractive index from a three-column file containing wavelength [nm],
and the real and imaginary parts of the refractive index (both positive numbers).
The Mie calculation is done for each wavelength defined here.
\end{description}

\option{size\_distribution\_file}
Specify a two column file, r [micron], dn(r)/dr, which describes a size distribution 
of droplets. 
\begin{Verbatim}
size_distribution_file file
\end{Verbatim}
The Mie calculation is repeated for each value of r found in the 
size distribution file, and the final result is a weighted average of these 
values. The user himself has to choose a set of r's suited for 
his specific purpose.

\option{temperature}
Ambient temperature, used to calculate the refractive indices of water and ice.
\begin{Verbatim}
temperature value
\end{Verbatim}
Temperature dependence is only considered above 10 micron (water) and 167 micron (ice), 
respectively. Default: 300K.

\option{verbose }
If specified abundances of informative messages are output to stderr. To make
use of this information, you may want to write the standard mie output to 
one file and the diagnostic messages to another. To do so, try
\code{(mie < mie.inp > mie.out) >\& verbose.txt} (depending on your shell you 
might need a slightly different syntax). 

\option{wavelength}
Sets the wavelength range, in nm. 
\begin{verbatim}
wavelength lambda_min lambda_max
\end{verbatim}
The wavelength step is specified by \code{wavelength\_step}.
For unregular wavelength grid it is also possible to specify a file, where the 
wavelength grid is stored.
\begin{verbatim}
wavelength wvl_filename
\end{verbatim}
where \code{wvl\_filename} is the path and name of the file, which contains 
the wavelength grid. It is expected that the wavelength values in nm are stored 
in the second column.
For the \code{correlated\_k} schemes implemented in \code{uvspec} you can use following abbrevations
instead of a filename (in this case \code{data\_files\_path} must be specified also): 
\code{kato}, \code{kato2}, \code{kato2andwandji}, \code{kato2.96}, \code{fu}, and \code{avhrr\_kratz}.
This option is ignored if \code{refrac file} is specified.
\option{wavelength\_step}
The wavelength step, in nm. Ignored if \code{refrac file} is specified.
\begin{verbatim}
wavelength_step value
\end{verbatim}

\option{wavelength\_index}
Set the wavelengths to be selected. This might be the normal wavelength grid defined by 
\code{wavelength} and \code{wavelength\_step} or a \code{correlated\_k} wavelength grid.
E.g., in combination with
\code{wavelength AVHRR\_KRATZ}, \code{wavelength\_index 15 15} will select wavelength
index 15 which corresponds to channel 4, or \code{wavelength\_index 10 14} will select 
those bands required for channel 3. Indices start from 1.

\end{description} 

</lpdoc>

      @c12_03@   */

  /***************** End of documentation *********************/


  /* prototype of main function */
  int mie(input_struct input, output_struct *output);

  int c;
  /* int iu, iv, j, lev; */
  int status=0;

  char infilename[FILENAME_MAX]="";
  
  /*  Input_struct Input;
   Output_struct Output;*/

  /* Initialization */

 

  Input.verbose         = 0;

  Input.medium          =  WATER;
  Input.temperature     =  300.0;
  Input.rho_user        = -999.0;
  Input.r_eff_min       = -999.0;
  Input.r_eff_max       = -999.0;
  Input.distribution    =  DIST_NONE;
  Input.alpha           = -999.0;
  Input.sd_filename     =  "";
  Input.ri_filename     =  "";
  Input.basename        =  ""; 
  Input.program         =  MIEV0;
  Input.nstokes         = 1;
  Input.nthetamax       = 1000;
  Input.n_r_max         = 5;
  Input.dx_max          = 0.03;

  Input.wl.wl_start     = 0.0;
  Input.wl.wl_end       = 0.0;
  Input.mie.anyang      = 1;
  Input.mie.ipolzn      = 0;
  Input.mie.momdim      = 0;
  Input.mie.nmom        = 0; 
  Input.mie.numang      = 0;
  Input.mie.mimcut      = 0;
  Input.mie.perfct      = 0;
  Input.mie.prnt[0]     = 0;
  Input.mie.prnt[1]     = 0;
  Input.mie.xmu         = (float *) calloc (Input.mie.numang+1, sizeof(float));
  Input.mie.s1          = (mie_complex *) calloc (Input.mie.numang+1, sizeof(mie_complex));
  Input.mie.s2          = (mie_complex *) calloc (Input.mie.numang+1, sizeof(mie_complex));
  Input.n_output_user   = 0;
  Input.ck_scheme       = CK_CRS;
  Input.lambda_i_start  = -999; 
  Input.lambda_i_end    = -999;
  Input.data_files_path = (char *) calloc (FILENAME_MAX,sizeof(char));
  Input.nmom_netcdf     = 0;
  Input.accuracy        = 0.01;

  
  while ((c=getopt (argc, argv, "hf:")) != EOF)  {
    switch(c)  {
    case 'h': 
      usage();
      return (-1);
      break;
    case 'f':
      strcpy (infilename, optarg);
      break;
    default:
      usage();
      return (-1);
    }
  }

  /* Read the input file */
  if (strlen(infilename)>0)
     yyin = fopen (infilename, "r");
  yylex();

  if (strlen(infilename)>0)
    fclose (yyin);
   
  /* initialisation dependent on input data */

  if (Input.r_eff_max == -999.)
    Input.r_eff_max = Input.r_eff_min;
  
  switch (Input.medium) {
    case WATER:
      Output.rho_medium = RHO_WATER;
      break;
    case ICE:
      Output.rho_medium = RHO_ICE;
      break;
    case USER:
      if (Input.rho_user > 0.0)
        Output.rho_medium = Input.rho_user;
      else 
        Output.rho_medium = 1.0;
      break;
    case AEROSOL: 
       fprintf(stderr, "Mie calculation for aerosol. Size distributions, refractive indices\n");
       fprintf(stderr, "and the aerosol density are taken from the OPAC database.\n");
       break; 
    default:
      fprintf (stderr, "Error, unknown medium %d\n", Input.medium);
      break; 
  }

   if ((Input.output_format == OUTPUT_FORMAT_AERPROP || Input.output_format == OUTPUT_FORMAT_CLOUDPROP ||  Input.output_format == OUTPUT_FORMAT_NETCDF) && 
        Input.mie.nmom == 0){
       fprintf(stderr, "Warning: You have chosen *output_user* aerosolprp, cloudprp, or netcdf.\n");
       fprintf(stderr, "         This option requires that Legendre polynomials are calculated.\n");
       fprintf(stderr, "         nmom is now set to a default value of 1000, please check whether this number is ok. \n");
       Input.mie.nmom   = 1000; 
       Input.mie.momdim = Input.mie.nmom+1;
    }

 
  /* checks */

  status = mie_check(Input);
  if (status!=0) {
    fprintf (stderr, "%d error(s) in mie input-file, aborting\n", status);
    return status;
  }

  /* Do the Mie calculation */

  status = mie(Input, &Output);
  if (status != 0.0) {
    fprintf (stderr, "Error during execution of mie\n");
    return status;
  }

  return status;

}

float yytext2float(char* yytext, int which_token, int max_token) {
  double dnumber;
  char token[] = " ";
  char *p, *s=NULL;
  int count;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  p = strtok(yytext2,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }
  if ((count != which_token) && (count>max_token)) {
    fprintf(stderr, "Wrong number of arguments (%d) on line %d\n", count, line_number+1);
  }
  dnumber = atof(s);
  /*   printf("dnumner %f\n",dnumber);*/
  free(yytext2);
  
  return (float)dnumber;
}

float* yytext2floats(char* yytext, int nary) {
  int i;
  char token[] = " ";
  char *p;
  char *yytext2=NULL;

  /* Fill up ary */
  float *ary    = (float *) calloc (nary, sizeof(float));

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  p = strtok(yytext2,token);
  i = 0;
  while ((p = strtok(NULL,token)) != NULL)
    ary[i++] = atof(p);

  return ary;
}

int yytext2nfloats(char* yytext) {
  int nary;
  char token[] = " ";
  char *p;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  /* Get number of ary elements */
  nary = 0;
  p = strtok(yytext2,token);
  while ((p = strtok(NULL,token)) != NULL) 
    nary++;

  return nary;
}

int* yytext2integers(char* yytext, int nary) {
  int i;
  char token[] = " ";
  char *p;
  char *yytext2=NULL;

  /* Fill up ary */
  int *ary    = (int *) calloc (nary, sizeof(int));

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  p = strtok(yytext2,token);
  i = 0;
  while ((p = strtok(NULL,token)) != NULL)
    ary[i++] = atoi(p);

  return ary;
}

int yytext2nintegers(char* yytext) {
  int nary;
  char token[] = " ";
  char *p;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  /* Get number of ary elements */
  nary = 0;
  p = strtok(yytext2,token);
  while ((p = strtok(NULL,token)) != NULL)
    nary++;

  return nary;
}

int yytext2int(char* yytext, int which_token, int max_token) {
  char token[] = " ";
  char *p, *s=NULL;
  int count;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  p = strtok(yytext2,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }

  if ((count != which_token) && (count>max_token))
    fprintf(stderr, "Wrong number of arguments on line %d\n", line_number+1);

  return atoi(s);
}

char* yytext2string(char* yytext, int which_token, int max_token) {
  char token[] = " ";
  char *p, *s=NULL, *string=NULL;
  int count;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  p = strtok(yytext2,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }
  if ((count != which_token) && (count>max_token)) {
    fprintf(stderr, "Wrong number of arguments (%d) on line %d\n", count, line_number+1);
  }
  string = (char *) calloc (strlen(s)+1, sizeof (char));
  strcpy (string, s);
   
  free(yytext2);   
  return string;
}

int yytext2ntokens(char* yytext) {
  int nary=0;
  char token[]=" ";
  char *p=NULL;
  char *yytext2=NULL;

  yytext2 = (char *) calloc (strlen(yytext)+1,sizeof(char));
  strcpy (yytext2, yytext);

  /* Get number of ary elements */
  nary = 0;
  p = strtok(yytext2,token);

  while ((p = strtok(NULL,token)) != NULL) 
    nary++;

  free(yytext2);
  return nary;
}

/* like yytext2string, but string is copied to a pre-allocated string dest */
int yytextCstring(char *dest, char* yytext, int which_token, int max_token) {
  char token[]=" ";
  char *p=NULL, *s=NULL, *tmpstring=NULL;
  int count=0;

  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }

  if ((count != which_token) && (count>max_token))
    fprintf(stderr, "Wrong number of arguments on line %d\n", line_number+1);

  strcpy (dest, s);

  free(tmpstring);

  return 0;
}

double yytext2double(char* yytext, int which_token, int max_token) {
  double dnumber=0;
  char token[]=" ";
  char *p=NULL, *s=NULL;
  int count=0;
  char *tmpstring=NULL;

  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }

  if ((count != which_token) && (count>max_token))
    fprintf(stderr,"Wrong number of arguments on line %d\n", line_number+1);

  dnumber = atof(s);

  free(tmpstring);
  return dnumber;
}




/**************************************************/
/* Check the Mie input structure for plausibility */
/**************************************************/

int mie_check (input_struct input)
{
  int status=0;
  status =0;

  if (!( input.nstokes==1 || input.nstokes==4)){
    fprintf(stderr, "Error: nstokes must be either 1 or 4 \n");
    status = -1;
   }

  if (input.nstokes >1 && input.program != MIEV0){
    fprintf(stderr, "Error: To calclate the phase matrix (*nstokes*=4) you must use the mie_program MIEV0.\n");
    status = -1;
   }
  
   
  if (input.mie.nmom >0 && input.program != MIEV0){
    fprintf(stderr, "Error: To calclate the Legendre polynomials you must use the mie_program MIEV0.\n");
    status = -1;
   } 

 if (Input.aerosol_type!=0 && Input.medium != AEROSOL){
     fprintf(stderr, "*** Warning: You have chosen an aerosol type and provided a refractive index.\n");
     fprintf(stderr, "*** The refractive index of the aerosol type will be used for the calculation,\n");
     fprintf(stderr, "*** not the one that you have provided. \n");
   }

   if (Input.aerosol_type!=0 && Input.distribution != DIST_AER){
        fprintf(stderr, "Error: You have chosen an aerosol type and provided a size distribution.\n");
        fprintf(stderr, "       This is not possible since the OPAC size distributions are \n");
        fprintf(stderr, "       taken automatically. Please remove the option *distribution*.\n");
        status = -1;
   }
    
   if (Input.medium != AEROSOL &&  Input.r_eff_min==-999 && Input.distribution == DIST_NONE ){
       fprintf(stderr, "Error: You have to provide a particle size or a particle size distribution \n"); 
       fprintf(stderr, "       to perform a Mie calculation. \n");
       fprintf(stderr, "       Please specify *r_eff* or *size_distribution_file*.\n"); 
       status = -1; 
    }

   if (Input.medium != AEROSOL && Input.output_format == OUTPUT_FORMAT_AERPROP){
       fprintf(stderr, "Error: You have chosen *output_user* aerosolprp but you want to do a calculation \n");
       fprintf(stderr, "       for clouds. Please use cloudprp instead. \n");
       status = -1;
    }

    if (Input.medium == AEROSOL && Input.output_format == OUTPUT_FORMAT_CLOUDPROP){
       fprintf(stderr, "Error: You have chosen *output_user* cloudprp but you want to do a calculation \n");
       fprintf(stderr, "       for aerosols. Please use aerosolprp instead. \n");
       status = -1;
     }

    if (Input.output_format == OUTPUT_FORMAT_CLOUDPROP && Input.distribution == DIST_NONE){
       fprintf(stderr, "Error: You have chosen *output_user* cloudprp. For this option you have \n");
       fprintf(stderr, "       to specifiy a size distribution using the option *distribution*. \n");
       status = -1; 
     }

   return status;
}
