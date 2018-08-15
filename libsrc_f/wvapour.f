C$$$ Algorithms, Comparisons and Source References by Schlatter and Baker
C 
C From schlatter@profsc.fsl.noaa.gov Wed Jun 12 16:25:46 1991
C 
C General information:
C --------------------
C This is an index of the thermodynamic subprograms located in THERMO.OLB.  They
C are listed according to the general type of calculation that is desired,
C divided into the following categories: 
C 
C         (a) moisture parameters
C         (b) latent heat
C         (c) pressure
C         (d) temperature
C         (e) thickness
C 
C These algorithms were collected, edited, commented, and tested by Thomas W.
C Schlatter and Donald V. Baker from August to October 1981 in the PROFS Program
C Office, NOAA Environmental Research Laboratories, Boulder, Colorado.  Where
C possible, credit has been given to the original author of the algorithm and a
C reference provided.
C 
C 
C The input/output units are as follows:
C 
C Temperature                     Celsius
C Pressure                        millibars
C Relative humidity               percent
C Saturation specific humidity    g vapor/kg moist air
C Mixing ratio                    g vapor/kg dry air
C Thickness                       meters
C Precipitable water              centimeters
C Latent heat                     joules/kg
C 
C 
C The following symbols are used in subprogram calls:
C 
C EW      water vapor pressure
C EQPT    equivalent potential temperature
C P       pressure
C PC      pressure at the convective condensation level
C PS      surface pressure
C RH      relative humidity
C T       temperature
C TC      temperature at the lifting condensation level
C TD      dew point
C TH      potential temperature (a dry adiabat)
C THW     wet bulb potential temperature (a moist adiabat)
C W       mixing ratio
C WBAR    mean mixing ratio from surface to pressure pm
C 
C Note:  all routines are function subprograms except for "PTLCL", which is a
C        subroutine.

C ###############################################################################
C
C TEST OF FUNCTIONS DPT,DEWPT
C PURPOSE: CALCULATE DEW POINT GIVEN SATURATION VAPOR PRESSURE.
C
C SAT.VP.    DPT    DEWPT    SMITH.(351)
C -------  ------  ------    ------
C 0.1891  -40.002  -40.025    -40
C 1.2540  -20.000  -20.032    -20
C 2.8627  -10.000  -10.022    -10
C 6.1078    0.000   -0.010      0
C 12.272   10.000   10.000     10
C 23.373   20.000   20.003     20
C
C EFFICIENCY TEST: 3000 CALCULATIONS AT SVP=.1891MB.
C
C           FUNC.  T(SEC)
C           -----  ------
C           DPT     1.92
C           DEWPT   0.59
C

C $$$-------------------------------------------------------------------------
C $$$-------------------------------------------------------------------------



C$$$  SUBPROGRAM DOCUMENTATION BLOCK 
C                .      .    .                                       . 
C SUBPROGRAM:    DPT       COMPUTE DEWPT GIVEN WATER VAPOR PRESSURE
C       Baker, Schlatter  17-MAY-1982     Original version.
C
C ABSTRACT: THIS FUNCTION RETURNS THE DEW POINT DPT (CELSIUS), GIVEN THE
C   WATER VAPOR PRESSURE EW (MILLIBARS). 
C
C USAGE:    DPT(EW)
C
C   INPUT ARGUMENT LIST: 
C     EW       - REAL    WATER VAPOR PRESSURE (MB)
C 
C   OUTPUT ARGUMENT LIST:   
C     DPT    - REAL    (CELSIUS)
C 
C REMARKS: 
C   APPROXIMATE DEW POINT BY MEANS OF TETEN'S FORMULA.
C   THE FORMULA APPEARS AS EQ.(8) IN BOLTON, DAVID, 1980:
C   "THE COMPUTATION OF EQUIVALENT POTENTIAL TEMPERATURE,"
C   MONTHLY WEATHER REVIEW, VOL 108, NO. 7 (JULY), P.1047.
C   THE FORMULA IS EW(T) = ES0*10**(7.5*T/(T+237.3))
C            OR    EW(T) = ES0*EXP(17.269388*T/(T+237.3))
C   THE INVERSE FORMULA IS USED BELOW.

        FUNCTION DPT(EW)

        IMPLICIT NONE

        REAL EW
        REAL ES0
        REAL DPT
        REAL X
        REAL DNM
        REAL T
        REAL FAC
        REAL EDP
        REAL ESW
        REAL DTDEW
        REAL DT
        DATA ES0/6.1078/

C   ES0 = SATURATION VAPOR PRESSURE (MB) OVER WATER AT 0C
C   RETURN A FLAG VALUE IF THE VAPOR PRESSURE IS OUT OF RANGE.

        IF (EW.GT..06.AND.EW.LT.1013.) GO TO 5
        DPT = 9999.
        RETURN
    5   CONTINUE

        X = ALOG(EW/ES0)
        DNM = 17.269388-X
        T = 237.3*X/DNM
        FAC = 1./(EW*DNM)

C   LOOP FOR ITERATIVE IMPROVEMENT OF THE ESTIMATE OF DEW POINT

   10   CONTINUE

C   GET THE PRECISE VAPOR PRESSURE CORRESPONDING TO T.

        EDP = ESW(T)

C   ESTIMATE THE CHANGE IN TEMPERATURE CORRESPONDING TO (EW-EDP)
C   ASSUME THAT THE DERIVATIVE OF TEMPERATURE WITH RESPECT TO 
C   VAPOR PRESSURE (DTDEW) IS GIVEN BY THE DERIVATIVE OF THE
C   INVERSE TETEN FORMULA.

        DTDEW = (T+237.3)*FAC
        DT = DTDEW*(EW-EDP)
        T = T+DT
        IF (ABS(DT).GT.1.E-04) GO TO 10
        DPT = T
        RETURN
        END


C$$$  SUBPROGRAM DOCUMENTATION BLOCK 
C                .      .    .                                       . 
C SUBPROGRAM:    DEWPT       COMPUTE DEWPT GIVEN WATER VAPOR PRESSURE
C   PRGMMR: STAN BENJAMIN    ORG: FSL/PROFS  DATE: 90-07-19 
C 
C ABSTRACT: C   THIS FUNCTION YIELDS THE DEW POINT DEWPT (CELSIUS), GIVEN THE
C   WATER VAPOR PRESSURE EW (MILLIBARS).
C 
C PROGRAM HISTORY LOG: 
C   82        DON BAKER         ORIGINAL VERSION
C 
C USAGE:    DEWPT(EW)
C
C   INPUT ARGUMENT LIST: 
C     EW       - REAL    WATER VAPOR PRESSURE (MB)
C 
C   OUTPUT ARGUMENT LIST:   
C     DEWPT    - REAL    (CELSIUS)
C 
C REMARKS: 
C   THE EMPIRICAL FORMULA APPEARS IN BOLTON, DAVID, 1980:
C   "THE COMPUTATION OF EQUIVALENT POTENTIAL TEMPERATURE,"
C   MONTHLY WEATHER REVIEW, VOL. 108, NO. 7 (JULY), P. 1047, EQ.(11).
C   THE QUOTED ACCURACY IS 0.03C OR LESS FOR -35 < DEWPT < 35C.
C 
C ATTRIBUTES: 
C   LANGUAGE: FORTRAN-77
C   MACHINE:  NAS-9000 
C$$$

        FUNCTION DEWPT(EW)

        IMPLICIT   NONE
        REAL DEWPT,EW,ENL

        ENL = ALOG(EW)
        DEWPT = (243.5*ENL-440.8)/(19.48-ENL) 
        RETURN
        END




C$$$  SUBPROGRAM DOCUMENTATION BLOCK 
C                .      .    .                                       . 
C SUBPROGRAM:    DEWPNT      COMPUTE DEWPNT GIVEN TEMP AND RH
C   PRGMMR: STAN BENJAMIN    ORG: FSL/PROFS  DATE: 90-07-19 
C 
C ABSTRACT: FUNCTION TO RETURN THE DEW POINT GIVEN TEMPERATURE
C   AND RELATIVE HUMIDITY.
C 
C PROGRAM HISTORY LOG: 
C   82        DON BAKER         ORIGINAL VERSION
C 
C USAGE:    DEWPNT(TK,RH)
C
C   INPUT ARGUMENT LIST: 
C     TK       - REAL    TEMPERATURE (KELVIN)
C     RH       - REAL    RELATIVE HUMIDITY (PERCENT)
C 
C   OUTPUT ARGUMENT LIST:   
C     DEWPNT   - REAL    (KELVIN)
C 
C REMARKS: 
C   THIS FUNCTION RETURNS THE DEW POINT (KELVIN) GIVEN THE TEMPERATURE
C   (KELVIN) AND RELATIVE HUMIDITY (%). THE FORMULA IS USED IN THE
C   PROCESSING OF U.S. RAWINSONDE DATA AND IS REFERENCED IN PARRY, H.
C   DEAN, 1969: "THE SEMIAUTOMATIC COMPUTATION OF RAWINSONDES,"
C   TECHNICAL MEMORANDUM WBTM EDL 10, U.S. DEPARTMENT OF COMMERCE,
C   ENVIRONMENTAL SCIENCE SERVICES ADMINISTRATION, WEATHER BUREAU,
C   OFFICE OF SYSTEMS DEVELOPMENT, EQUIPMENT DEVELOPMENT LABORATORY,
C   SILVER SPRING, MD (OCTOBER), PAGE 9 AND PAGE II-4, LINE 460.
C 
C ATTRIBUTES: 
C   LANGUAGE: FORTRAN-77
C   MACHINE:  NAS-9000 
C$$$
C*******************************************************
C
      FUNCTION DEWPNT(TK,RH)

      IMPLICIT NONE

      REAL TK
      REAL RH
      REAL T
      REAL X
      REAL DPD
      REAL DWPTC
      REAL DEWPNT

      IF(TK.GT.150.) THEN
       T = TK -273.15
      ELSE
       T = TK
       ENDIF
      X = 1.-0.01*RH

C   COMPUTE DEW POINT DEPRESSION.

      DPD =(14.55+0.114*T)*X+((2.5+0.007*T)*X)**3+(15.9+0.117*T)*X**14
      DWPTC = T-DPD
      DEWPNT= DWPTC + 273.15
      RETURN
       END



C$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
C                .      .    .                                       .  
C SUBPROGRAM:    TCON        COMPUTES TEMPERATURE AT THE LCL            
C   PRGMMR: PATTY MILLER     ORG: FSL/PROFS  DATE: 90-06-15             
C                                                                       
C ABSTRACT: THIS FUNCTION RETURNS THE TEMPERATURE TCON (CELSIUS) AT     
C   THE LIFTING CONDENSATION LEVEL, GIVEN THE TEMPERATURE T (CELSIUS)   
C   AND THE DEW POINT D (CELSIUS).                                      
C                                                                       
C PROGRAM HISTORY LOG:                                                  
C   82-05-17  D. BAKER, T. SCHLATTER   ORIGINAL VERSION                 
C                                                                       
C USAGE:    TCON(T,D)                                                   
C   INPUT ARGUMENT LIST:                                                
C     T        - REAL TEMPERATURE (C)                                   
C     D        - REAL DEWPOINT TEMPERATURE (C)                          
C                                                                       
C   OUTPUT ARGUMENT LIST:                                               
C                                                                       
C REMARKS: NONE                                                         
C                                                                       
C ATTRIBUTES:                                                           
C   LANGUAGE: FORTRAN-77                                                
C   MACHINE:  NAS-9000                                                  
C$$$                                                                    
        FUNCTION TCON(T,D)                                              
        IMPLICIT NONE

        REAL T
        REAL D
        REAL S
        REAL DLT
        REAL TCON
C                                                                       
C   COMPUTE THE DEW POINT DEPRESSION S.                                 
        S = T-D                                                         
C   THE APPROXIMATION BELOW, A THIRD ORDER POLYNOMIAL IN S AND T,       
C   IS DUE TO HERMAN WOBUS. THE SOURCE OF DATA FOR FITTING THE          
C   POLYNOMIAL IS UNKNOWN.                                              
                                                                        
        DLT = S*(1.2185+1.278E-03*T+                                    
     1        S*(-2.19E-03+1.173E-05*S-5.2E-06*T))                      
        TCON = T-DLT                                                   
        RETURN                                                         
        END                                                             

C$$$   SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    TEMP       CALCULATES TEMPERATURE FROM POTENTIAL TEMP AND PRES
C
C   PRGMMR:  BENJAMIN, STAN ORG: ERL/PROFS      DATE: 93-01-12
C
C ABSTRACT:  CALCULATES TEMPERATURE BY SOLVING POISSON EQUATION
C
C PROGRAM HISTORY LOG:
C
C USAGE:           TEMP = TEMP(TP, P)
C
C   INPUT ARGUMENT LIST:
C   TP          - REAL           POTENTIAL TEMPERATURE IN KELVIN
C   P           - REAL           PRESSURE IN MB
C
C   OUTPUT ARGUMENT LIST:
C   TEMP        - REAL           TEMPERATURE IN KELVIN
C
C REMARKS: NONE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN-77 + EXTENSIONS
C   MACHINE: DEC - VAX, VMS
C
C$$$
C
C**********************************************************************


C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    ESW         COMPUTE THE SATURATION VAPOR PRESSURE
C   PRGMMR: STAN BENJAMIN    ORG: FSL/PROFS  DATE: 90-07-19
C
C ABSTRACT: FUNCTION TO RETURNS THE SATURATION VAPOR PRESSURE
C   ESW (MILLIBARS) OVER LIQUID WATER GIVEN THE TEMPERATURE
C   T (CELSIUS OR KELVIN).
C
C PROGRAM HISTORY LOG:
C   82-05-17  TOM SCHLATTER   ORIGINAL VERSION
C
C USAGE:    ESW(T)
C   INPUT ARGUMENT LIST:
C     T        - REAL    TEMPERATURE (CELSIUS OR KELVIN)
C
C   OUTPUT ARGUMENT LIST:
C
C REMARKS: THE POLYNOMIAL APPROXIMATION BELOW IS DUE TO HERMAN WOBUS,
C   A MATHEMATICIAN WHO WORKED AT THE NAVY WEATHER RESEARCH FACILITY,
C   NORFOLK, VIRGINIA. THE COEFFICIENTS OF
C   THE POLYNOMIAL WERE CHOSEN TO FIT THE VALUES IN TABLE 94 ON
C   PP. 351-353 OF THE SMITHSONIAN METEOROLOGICAL TABLES BY ROLAND
C   LIST (6TH EDITION). THE APPROXIMATION IS VALID FOR
C   -50 < T < 100C.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN-77
C   MACHINE:  NAS-9000
C$$$
        FUNCTION ESW(T1)

        IMPLICIT NONE

        REAL T1
        REAL T
        REAL POL
        REAL ESW
        REAL ES0
C
C   THIS FUNCTION RETURNS THE SATURATION VAPOR PRESSURE ESW (MILLIBARS)
C   OVER LIQUID WATER GIVEN THE TEMPERATURE T (CELSIUS).
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C
C   THE POLYNOMIAL APPROXIMATION BELOW IS DUE TO HERMAN WOBUS, A 
C   MATHEMATICIAN WHO WORKED AT THE NAVY WEATHER RESEARCH FACILITY, 
C   NORFOLK, VIRGINIA. THE COEFFICIENTS OF THE 
C   POLYNOMIAL WERE CHOSEN TO FIT THE VALUES IN TABLE 94 ON PP. 351-353
C 
C   OF THE SMITHSONIAN METEOROLOGICAL TABLES BY ROLAND LIST (6TH 
C   EDITION). THE APPROXIMATION IS VALID FOR -50 < T < 100C.
C
C   ES0 = SATURATION VAPOR RESSURE OVER LIQUID WATER AT 0C
        DATA ES0/6.1078/
        T = T1
        IF(T.GT.150.) T=T-273.15
        POL = 0.99999683       + T*(-0.90826951E-02 +
     1     T*(0.78736169E-04   + T*(-0.61117958E-06 +
     2     T*(0.43884187E-08   + T*(-0.29883885E-10 +
     3     T*(0.21874425E-12   + T*(-0.17892321E-14 +
     4     T*(0.11112018E-16   + T*(-0.30994571E-19)))))))))
        ESW = ES0/POL**8
        RETURN
         END



C$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
C                .      .    .                                       .  
C SUBPROGRAM:    WMR         APPROXIMATES MIXING RATIO                  
C   PRGMMR: PATTY MILLER     ORG: FSL/PROFS  DATE: 90-06-15             
C                                                                       
C ABSTRACT: THIS FUNCTION APPROXIMATES THE MIXING RATIO WMR             
C   (GRAMS OF WATER VAPOR PER KILOGRAM OF DRY AIR) GIVEN THE            
C   PRESSURE P (MB) AND THE TEMPERATURE T (CELSIUS). THE FORMULA        
C   USED IS GIVEN ON P. 302 OF THE SMITHSONIAN METEOROLOGICAL TABLES    
C   BY ROLAND LIST (6TH EDITION).                                       
C                                                                       
C PROGRAM HISTORY LOG:                                                  
C   82-05-17  D. BAKER, T. SCHLATTER   ORIGINAL VERSION                 
C                                                                       
C USAGE:    WMR(P,T)                                                    
C   INPUT ARGUMENT LIST:                                                
C     P        - REAL PRESSURE (MB)                                     
C     T        - REAL TEMPERATURE (C)                                   
C                                                                       
C   OUTPUT ARGUMENT LIST:                                               
C                                                                       
C REMARKS: NONE                                                         
C                                                                       
C ATTRIBUTES:                                                           
C   LANGUAGE: FORTRAN-77                                                
C   MACHINE:  NAS-9000                                                  
C$$$    
                                                                
C        FUNCTION WMR(P,T)                                               
C                                                                       
C   EPS = RATIO OF THE MEAN MOLECULAR WEIGHT OF WATER (18.016 G/MOLE)   
C         TO THAT OF DRY AIR (28.966 G/MOLE)                            
C                                                                       
C        DATA EPS/0.62197/                                               
C                                                                       
C   THE NEXT TWO LINES CONTAIN A FORMULA BY HERMAN WOBUS FOR THE        
C   CORRECTION FACTOR WFW FOR THE DEPARTURE OF THE MIXTURE OF AIR       
C   AND WATER VAPOR FROM THE IDEAL GAS LAW. THE FORMULA FITS VALUES     
C   IN TABLE 89, P. 340 OF THE SMITHSONIAN METEOROLOGICAL TABLES,       
C   BUT ONLY FOR TEMPERATURES AND PRESSURES NORMALLY ENCOUNTERED IN     
C   IN THE ATMOSPHERE.                                                  
                                                                        
C        X = 0.02*(T-12.5+7500./P)                                       
C        WFW = 1.+4.5E-06*P+1.4E-03*X*X                                  
C        FWESW = WFW*ESW(T)                                             
C        R = EPS*FWESW/(P-FWESW)                                         
                                                                        
C   CONVERT R FROM A DIMENSIONLESS RATIO TO GRAMS/KILOGRAM.             
                                                                        
C        WMR = 1000.*R                                                   
C        RETURN                                                          
C        END                                                             


C$$$  SUBPROGRAM DOCUMENTATION BLOCK                                    
C                .      .    .                                       .  
C SUBPROGRAM:    EPT         COMPUTES EQUIV. POT. TEMP.                
C   PRGMMR: PATTY MILLER     ORG: FSL/PROFS  DATE: 90-06-15             
C                                                                      
C ABSTRACT:  THIS FUNCTION RETURNS THE EQUIVALENT                       
C   POTENTIAL TEMPERATURE EPT                                          
C   (CELSIUS) FOR A PARCEL OF AIR INITIALLY AT TEMPERATURE T (CELSIUS), 
C   DEW POINT TD (CELSIUS) AND PRESSURE P (MILLIBARS). THE FORMULA USED 
C   IS EQ.(43) IN BOLTON, DAVID, 1980: "THE COMPUTATION OF EQUIVALENT   
C   POTENTIAL TEMPERATURE," MONTHLY WEATHER REVIEW, VOL. 108, NO. 7  
C   (JULY), PP. 1046-1053. THE MAXIMUM ERROR IN EPT IN 0.3C.  IN MOST   
C   CASES THE ERROR IS LESS THAN 0.1C.                                  
C                                                                       
C PROGRAM HISTORY LOG:                                                  
C   82-05-17  T. SCHLATTER           ORIGINAL VERSION                  
C                                                                       
C USAGE:    EPT(T,TD,P)                                                
C   INPUT ARGUMENT LIST:                                                
C     P        - REAL PRESSURE (MB)                                     
C     TD       - REAL DEWPOINT TEMPERATURE (C)                          
C     T        - REAL TEMPERATURE (C)                                   
C                                                                      
C   OUTPUT ARGUMENT LIST:                                              
C                                                                       
C REMARKS: NONE                                                        
C                                                                      
C ATTRIBUTES:                                                           
C   LANGUAGE: FORTRAN-77                                               
C   MACHINE:  NAS-9000                                                 
C$$$                      


C###############################################################################
C
C TEST OF FUNCTIONS OE,EPT
C PURPOSE: CALCULATE EQUIVALENT POTENTIAL TEMPERATURE GIVEN TEMP-
C          ERATURE, DEW POINT, AND PRESSURE.
C
C TEMP  DWPT  PRES    OE     EPT
C ----  ----  ----   -----  -----
C  30    15   1000   62.24  62.49
C   0   -20    700   33.10  33.01
C  10     0    850   36.97  37.01
C -15   -25    500   45.04  45.02
C  25    20    900   85.02  84.60
C
C EFFICIENCY TEST: 1000 CALCULATIONS AT T=20C, TD=10C, P=1000MB.
C
C           FUNC.  T(SEC)
C           -----  ------
C           OE      18.50
C           EPT      0.80
C                                            
C        FUNCTION EPT(T,TD,P)                                           
C                                                                      
C   COMPUTE THE MIXING RATIO (GRAMS OF WATER VAPOR PER KILOGRAM OF      
C   DRY AIR).                                                           
                                                                        
C        W = WMR(P,TD)  
C        write(*,*) "hallo world"                                                      
                                                                        
C   COMPUTE THE TEMPERATURE (CELSIUS) AT THE LIFTING CONDENSATION LEVEL.
                                                                        
C        TLCL = TCON(T,TD)                                               
C        TK = T+273.15                                                   
C        TL = TLCL+273.15                                                
C        PT = TK*(1000./P)**(0.2854*(1.-0.00028*W))                      
C        EPTK = PT*EXP((3.376/TL-0.00254)*W*(1.+0.00081*W))              
C        EPT= EPTK-273.15                                                
C        RETURN                                                         
C        END                                                            


        FUNCTION EPT(T,TD,P,N_AIR,N_H2O)                                           
C                                                                      
C   REPLACED THE WMR (COMPUTE THE WATER MIXING RATIO) FUNCTION AS 
C   IT IS NOT VALID FOR ALL WANTED TEMPERATURES AND PRESSURES !!!            
        IMPLICIT NONE   

        REAL T
        REAL TD
        REAL P
        REAL N_AIR
        REAL N_H2O
        REAL CP   
        REAL CV
        REAL kappa
        REAL EPTK
        REAL EPT
        REAL TLCL
        REAL PT
        REAL TL
        REAL TK
        REAL W
        REAL TCON


        CP = 1003.8  
        CV = 716.87
        kappa = (CP-CV) / CP

        W = N_H2O * 18.015 /  (N_AIR * 28.966) * 1000.      
C   1000 == kg(water)/kg(dry air) -> g(water)/kg(dry air)                                                  
                                                      
C   COMPUTE THE TEMPERATURE (CELSIUS) AT THE LIFTING CONDENSATION LEVEL.
                                                                        
        TLCL = TCON(T,TD)                                               
        TK = T+273.15                                                   
        TL = TLCL+273.15                                                
        PT = TK*(1000./P)**(kappa*(1.-0.00028*W))                      
        EPTK = PT*EXP((3.376/TL-0.00254)*W*(1.+0.00081*W))              
        EPT= EPTK-273.15                                                
        RETURN                                                         
        END              
