*DECK SWDE                                                                 SWDE.1     
      SUBROUTINE SWDE (PGG,PREF,PRMUZ,PTO1,PW                              SWDE.2     
     S     ,PRE1,PRE2,PTR1,PTR2         )                                  SWDE.3     
C                                                                          SWDE.4     
C**** *SWDE* - DELTA-EDDINGTON IN A CLOUDY LAYER                           SWDE.5     
C                                                                          SWDE.6     
C     PURPOSE.                                                             SWDE.7     
C     --------                                                             SWDE.8     
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY       SWDE.9     
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.                     SWDE.10    
C                                                                          SWDE.11    
C**   INTERFACE.                                                           SWDE.12    
C     ----------                                                           SWDE.13    
C          *SWDE* IS CALLED BY *SWR*, *SW2S*                               SWDE.14    
C                                                                          SWDE.15    
C     *CALL*     SWDE (KDLON,PGG,PREF,PRMUZ,PTO1,PW                        SWDE.16    
C    S                ,      PRE1,PRE2,PTR1,PTR2         )                 SWDE.17    
C                                                                          SWDE.18    
C        EXPLICIT ARGUMENTS :                                              SWDE.19    
C        --------------------                                              SWDE.20    
C PGG    : (NDLON)             ; ASSYMETRY FACTOR                          SWDE.21    
C PREF   : (NDLON)             ; REFLECTIVITY OF THE UNDERLYING LAYER      SWDE.22    
C PRMUZ  : (NDLON)             ; COSINE OF SOLAR ZENITH ANGLE              SWDE.23    
C PTO1   : (NDLON)             ; OPTICAL THICKNESS                         SWDE.24    
C PW     : (NDLON)             ; SINGLE SCATTERING ALBEDO                  SWDE.25    
C     ==== OUTPUTS ===                                                     SWDE.26    
C PRE1   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING NO            SWDE.27    
C                              ; REFLECTION FROM UNDERLYING LAYER          SWDE.28    
C PTR1   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING NO          SWDE.29    
C                              ; REFLECTION FROM UNDERLYING LAYER          SWDE.30    
C PRE2   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING               SWDE.31    
C                              ; REFLECTION FROM UNDERLYING LAYER          SWDE.32    
C PTR2   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING             SWDE.33    
C                              ; REFLECTION FROM UNDERLYING LAYER          SWDE.34    
C                                                                          SWDE.35    
C        IMPLICIT ARGUMENTS :   NONE                                       SWDE.36    
C        --------------------                                              SWDE.37    
C                                                                          SWDE.38    
C     METHOD.                                                              SWDE.39    
C     -------                                                              SWDE.40    
C                                                                          SWDE.41    
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.                    SWDE.42    
C                                                                          SWDE.43    
C     EXTERNALS.                                                           SWDE.44    
C     ----------                                                           SWDE.45    
C                                                                          SWDE.46    
C          NONE                                                            SWDE.47    
C                                                                          SWDE.48    
C     REFERENCE.                                                           SWDE.49    
C     ----------                                                           SWDE.50    
C                                                                          SWDE.51    
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND             SWDE.52    
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IN CORE MODEL"    SWDE.53    
C                                                                          SWDE.54    
C     AUTHOR.                                                              SWDE.55    
C     -------                                                              SWDE.56    
C        JEAN-JACQUES MORCRETTE  *ECMWF*                                   SWDE.57    
C                                                                          SWDE.58    
C     MODIFICATIONS.                                                       SWDE.59    
C     --------------                                                       SWDE.60    
C        ORIGINAL : 88-12-15                                               SWDE.61
C
C    Removed the longitude index, because in uvspec swde is only 
c    called for one location at a time. (CE, 2006-12-05)
c      
C     ------------------------------------------------------------------   SWDE.62    
C                                                                          SWDE.63    
C*       0.1   ARGUMENTS                                                   SWDE.64    
C              ---------                                                   SWDE.65    
C                                                                          SWDE.66    
c     REAL PGG(KDLON),PREF(KDLON),PRMUZ(KDLON),PTO1(KDLON),PW(KDLON)       SWDE.67    
c     REAL PRE1(KDLON),PRE2(KDLON),PTR1(KDLON),PTR2(KDLON)                 SWDE.68    

      REAL PGG, PREF, PRMUZ, PTO1, PW                                      CE.01
      REAL PRE1, PRE2,PTR1,PTR2                                            CE.02  
C                                                                          SWDE.69    
C     ------------------------------------------------------------------   SWDE.70    
C                                                                          SWDE.71    
C*       0.2   FUNCTIONS                                                   SWDE.72    
C              ---------                                                   SWDE.73    
C                                                                          SWDE.74    
CDIR$ VFUNCTION EXPHF                                                      SWDE.75    
*IF DEF,NEC                                                                US101298.3     
      EXPHF(X)=EXP(X)                                                      US101298.4     
*ENDIF                                                                     US101298.5     
C     ------------------------------------------------------------------   SWDE.76    
C                                                                          SWDE.77    
C*         1.      DELTA-EDDINGTON CALCULATIONS                            SWDE.78    
C                                                                          SWDE.79    
 100  CONTINUE                                                             SWDE.80    
C                                                                          SWDE.81    
c     DO 131 JL   =   1 , KDLON                                            SWDE.82
         
c      DO 131 JL   =   1 , 1                                                CE.03
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS                    SWDE.84    
C                                                                          SWDE.85    
 110  CONTINUE                                                             SWDE.86    

c      write(*,*) "pgg, pref, prmuz",  PGG, PREF, PRMUZ

      ZFF = PGG*PGG
      ZGP = PGG/(1.+PGG)          
      ZTOP = (1.- PW * ZFF) * PTO1
      ZWCP = (1-ZFF)* PW /(1.- PW * ZFF)                   
      ZDT = 2./3.                                                          SWDE.92    
      ZX1 = 1.-ZWCP*ZGP                                                    SWDE.93    
      ZWM = 1.-ZWCP                                                        SWDE.94    
      ZRM2 =  PRMUZ * PRMUZ                                
      ZRK = SQRT(3.*ZWM*ZX1)                                               SWDE.96    
      ZX2 = 4.*(1.-ZRK*ZRK*ZRM2)                                           SWDE.97    
      ZRP=ZRK/ZX1                                                          JM920928.1     
      ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2                               SWDE.99    
      ZBETA = 3.*ZWCP* PRMUZ *(1.+3.*ZGP*ZRM2*ZWM)/ZX2     
      ZARG=MIN(ZTOP/PRMUZ,87.)                                            JM920101.1     
      ZEXMU0=EXPHF(-ZARG)                                                  US070593.27    
      ZARG2=MIN(ZRK*ZTOP,87.)                                              JM920101.3     
      ZEXKP=EXPHF(ZARG2)                                                   US070593.28    
      ZEXKM = 1./ZEXKP 
      ZXP2P = 1.+ZDT*ZRP                                                   SWDE.104   
      ZXM2P = 1.-ZDT*ZRP                                                   SWDE.105   
      ZAP2B = ZALPHA+ZDT*ZBETA                                             SWDE.106   
      ZAM2B = ZALPHA-ZDT*ZBETA                                             SWDE.107   
C                                                                          SWDE.108   
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER            SWDE.109   
C                                                                          SWDE.110   
 120  CONTINUE                                                             SWDE.111   
C                                                                          SWDE.112   
      ZA11 = ZXP2P                                                         SWDE.113   
      ZA12 = ZXM2P                                                         SWDE.114   
      ZA13 = ZAP2B                                                         SWDE.115   
      ZA22 = ZXP2P*ZEXKP                                                   SWDE.116   
      ZA21 = ZXM2P*ZEXKM                                                   SWDE.117   
      ZA23 = ZAM2B*ZEXMU0                                                  SWDE.118   
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12                                    SWDE.119
c      write(*,*) "zdena",  ZDENA
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA                                   SWDE.120   
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA                                   SWDE.121   
      ZRI0A = ZC1A+ZC2A-ZALPHA                                             SWDE.122   
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA                                        SWDE.123   
      PRE1 = (ZRI0A-ZDT*ZRI1A)/ PRMUZ                      
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0                          SWDE.125   
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0                     SWDE.126   
      PTR1 = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ               
C                                                                          SWDE.128   
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER               SWDE.129   
C                                                                          SWDE.130   
 130  CONTINUE                                                             SWDE.131   
C                                                                          SWDE.132   
      ZB21 = ZA21- PREF *ZXP2P*ZEXKM                      
      ZB22 = ZA22- PREF *ZXM2P*ZEXKP                      
      ZB23 = ZA23- PREF *ZEXMU0*(ZAP2B - PRMUZ )          
      ZDENB = ZA11 * ZB22 - ZB21 * ZA12                                    SWDE.136   
c      write(*,*) "zdenB",  ZDENB, ZDENB- ZDENA
      ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB                                   SWDE.137   
      ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB                                   SWDE.138   
      ZRI0C = ZC1B+ZC2B-ZALPHA                                             SWDE.139   
      ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA                                        SWDE.140   
      PRE2 = (ZRI0C-ZDT*ZRI1C) / PRMUZ                    
      ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0                      SWDE.142   
      ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0               SWDE.143   
      PTR2 = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ         
C                                                                          SWDE.145   
c 131  CONTINUE                                                             SWDE.146   
      RETURN                                                               SWDE.147   
      END                                                                  SWDE.148   
