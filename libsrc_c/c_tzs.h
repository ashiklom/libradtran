
int c_tzs( int nlyr, float *dtauc, int nlev_c, float *zd_c,
	   int nzout, float *zout, 
	   float *ssalb, float *temper, 
	   float wvnmlo, float wvnmlhi, 
	   int usrtau, int ntau, float *utau, 
	   int usrang, int numu, float *umu, int nphi, float *phi, 
	   float albedo, double btemp, float ttemp, 
	   float temis, int planck, 
	   int *prnt, char *header,                //  int prndis[7]; char header[127]
	   float *rfldir, float *rfldn, float *flup,
	   float *dfdt, float *uavg,
	   float ***uu, int quiet );

void errmsg_tzs(char *messag, int   type);

void prtinp_tzs( int nlyr, int maxnlyr, float *dtauc, float *ssalb, float *temper,
		 float wvnmlo, float wvnmhi, int ntau, float *utau, int numu, float *umu,
		 int nphi, float *phi, float albedo, double btemp, float ttemp, float temis,
		 double *tauc, int nlev_c, float *zd_c, int nzout, float *zout );

void chekin_tzs( int nlyr, int maxnlyr, float *dtauc, int nlev_c, float *zd_c,
		 int nzout, float *zout, 
		 float *ssalb, float *temper, 
		 float wvnmlo, float wvnmlhi, 
		 int usrtau, int ntau, float *utau, 
		 int usrang, int numu, float *umu, int nphi, float *phi, 
		 float albedo, double btemp, float ttemp, 
		 float temis, int planck,
		 double *tauc, int quiet );

#define TZS_WARNING 0
#define TZS_ERROR   1
