      SUBROUTINE GEOFACPR
     $       (brosza, sza_bro, fac, nfac, szaloc,
     $                fac_2, NFAC_2, szaloc_2,
     $            z_lay, nlyr, zd, zenang, Re, VN)
*     
*     
*     
*     ********************************************************************
*
*     Calculates the geometric correction factor needed for spherical
*     geometry.
*     INCLUDING REFRACTION EFFECT.  FOR THIS EFFECT, GEOFACPR IS
*     SLOWER BUT MORE PRECISE THAN GEOFAST.    
*
*     INPUTS:
*            brosza: if =1, BrO density as function of alt and SZA
*            z_lay: Where in the layer the Chapman function is to be
*                   computed. E.g. 0.0, bottom of layer, 0.5, middle
*                   of layer and 1.0 top of layer.
*            nlyr:  Number of layers in atmospheric model
*            zd(lc):lc = 0, nlyr. zd(lc) is distance from bottom 
*                   surface to top of layer lc. zd(nlyr) = 0.0 km
*            zenang:Solar zenith angle as seen from bottom surface
*            Re:    Radius of earth. NOTE: Use the same dimension as zd,
*                   in km.
*            VN: refraction indexes of the layers
*
*     OUTPUTS:
*            nfac:
*            fac:   geometrical correction factor ds/dh
*            szaloc: local solar zenithal angles
*            fac_2: ds/dh beyond the tangent point of the path
*            szaloc_2: local sza beyond the tangent point
*         
*     
*     ********************************************************************  
*
*
      
      INCLUDE 'DISORT.MXD'
      integer nlyr, nfac(mxcly), NFAC_2(mxcly), brosza 
      real*8  z_lay, zd(0:*), zenang, Re, 
     $        VN(0:*), 
     $        SZALOC(mxcly,0:mxcly),
     $        SZALOC_2(mxcly,0:mxcly), 
     $        fac(mxcly,mxcly), 
     $        fac_2(mxcly,mxcly),sza_bro(mxsza)
     
*
*      
*     VZ: altitudes in cm   
*
      real*8 VZ(0:mxcly)
*      
      do i=0,nlyr
	VZ(i)=zd(i)*1.d+5
      enddo
*
*     Earth radius: km --> cm
*     
      Re=Re*1.d+5
*
      do i=1,nlyr
        nfac(i)  =0
	NFAC_2(i)=0
	do j = 1, nlyr
	  fac(i,j)  =0.D0
	  fac_2(i,j)=0.D0
	enddo
      enddo      
      do i=1, nlyr
        do j=0,nlyr
	  SZALOC(i,j)  = 0.D0
	  SZALOC_2(i,j)= 0.D0
	enddo
      enddo   
*            
      nbmax=nlyr
     

      CALL opathpr(zenang, VZ, VN, Re,
     $  	 nlyr, z_lay, FAC,
     $  	 SZALOC, FAC_2, SZALOC_2,
     $  	 NFAC, NFAC_2, nbmax)
     
      
      
      if (brosza.eq.1) then	       
        do i=1,nlyr
          do j=0,nlyr    
*         If local sza is lower than the minimum sza of BrO matrix,
*         put its value to that minimum	   
            if(szaloc(i,j).lt.sza_bro(1))
     $                     szaloc(i,j)  = sza_bro(1)
            if(szaloc_2(i,j).lt.sza_bro(1))
     $                     szaloc_2(i,j)= sza_bro(1)
          enddo
        enddo 
      endif
       
      end
      
      
      
****************************************************************************      
      
      
      

      SUBROUTINE opathpr (zenang, VZ, VN, Re,
     $                  nlyr, z_lay, FAC,
     $                  SZALOC, FAC_2, SZALOC_2,
     $                  NFAC, NFAC_2, nbmax)
*
*
*    ***********************************************************   
*     Compute optical path and local SZA in successive 
*     atmospheric layers.  
*     Refraction is taken into account.
*
*     INPUTS:
*       zenang	: SZA out of atmosphere in degrees
*       VZ	: Altitude vector [CM]  
*       VN	: refraction index vector
*       nlyr    : number of layers of atmosphere
*       Re:       radius of earth [CM]
*
*     OUTPUTS:
*       fac:      ds/dh geometrical factors before (including)
*                 tangent point
*       SZALOC:   local SZAs before (including) tangent point
*       fac_2:    ds/dh geometrical factors after tangent point
*       SZALOC_2: local SZAs after tangent point 
*  
*    ***********************************************************
*     
* 
      INCLUDE 'DISORT.MXD'
      integer nfac(mxcly),nfac_2(mxcly)
      real*8 zenang, Re, z_lay
      real*8 VN(0:*), VZ(0:mxcly),
     $       FAC(mxcly,mxcly),
     $       FAC_2(mxcly,mxcly),
     $       SZALOC(mxcly,0:mxcly),
     $       SZALOC_2(mxcly,0:mxcly)
    
      
      real*8 ds, deltaz, DegToRad, RadToDeg, 
     $       Rgoal, Zi, r0, r1, r2, BETAGOAL,   
     $       ca, theta0, theta1, s0, s1, beta0, beta1, 
     $       hh, Ratm0, ratm1, n0, n1, Rreach,
     $       dsdx0, dsdx1, dbetadx0, dbetadx1, Rlev, sign, 
     $       x0, x1, dx, c, rt,PI,rad(mxcly+1),n_ref(mxcly+1),
     $       DR, DR0, DR1, DZ0, DZ1, DERIVDR, deltas1,
     $       theta2, precis
  
  
  
      PI=DACOS(-1.D0)  
      DegToRad=pi/180.d0
      RadToDeg=1/DegToRad
      
*     *********************************************
*     STEP SIZE OF THE PATHS THROUGH THE ATMOSPHERE
*     *********************************************      
         
      ds=VZ(0)/dfloat(nbmax)

     
*     *********************************************
*     PRECISION TO BE REACHED ON THE LOCAL VERTICAL
*     (1 meter, except for 90<SZA<90.6: 1 cm)
*     *********************************************

      precis = 1.d2
      if(zenang.gt.90.d0.and.
     $   zenang.lt.90.6) precis= 1.d0 
     
*     ************************************
*     REFRACTION INDEXES OF THE ATMOSPHERE
*     ************************************


      do i=1,nlyr+1
        rad(i)  = Re+VZ(i-1)
	n_ref(i)= VN(i-1)
      enddo      
      

*     *************************************
*     AIMED POINTS ON THE LOCAL VERTICAL,  
*     FROM TOP OF ATMOSPHERE DOWN TO GROUND
*     *************************************

    
      DO i=1,nlyr
          
	
	  kount=0
	  kount2=0
	       
          deltaz=0.D0
*           
*         ALTITUDE AND RADIUS OF i_TH AIMED POINT
*         	  
          Zi = VZ(i) + z_lay*(VZ(i-1)-VZ(i))
	  Rgoal=Re+Zi


*         ********************************************************		
*         LOOP (PATHS plottingS) UNTIL REACHED POINT AND AIMED 
*         POINT ON THE LOCAL VERTICAL ARE CLOSE ENOUGH (0.1 METER)
*         ********************************************************
 
 
          ibissec=0
	  isecant=0
*         	  
	  IMETHOD=0    
*===        	
  2       continue  
*===        
         
	  
          nfac(i)  =0
	  nfac_2(i)=0
         
          do ij=1,nlyr
	    fac(i,ij)  =0.D0
	    fac_2(i,ij)=0.D0
          enddo
          
          do ij=0,nlyr
	    SZALOC(i,ij)  = 0.D0
	    SZALOC_2(i,ij)= 0.D0
  	  enddo
       
	  
*         PROPERTIES OF INITIAL POINT OF PATH (mainly fct of deltaz)
          r0       = Re + VZ(0)	  
	  call POINTS_OF_INCIDENCE_PR
     $         (VZ,Zi+deltaz,zenang,nlyr,Re,ca)
	  theta0   = zenang*DegToRad - ca
	  x0       = -r0*dcos(theta0)
   	  s0       = 0.d0
	  beta0    = 0.d0
	  call indrefpr (nlyr+1, r0, Re,
     $                 rad, n_ref, n0, hh, Ratm0)


          dsdx0    = 1/(1-Ratm0*(DSIN(theta0)**2))
	  dbetadx0 = dsdx0*DSIN(theta0)/r0 
	
*         BETA ANGLE TO BE REACHED (OVERTAKEN) BY PATH
          BETAGOAL = ca
	
*         OPTICAL CONSTANT OF PATH	
	  c        = n0*r0*dsin(theta0)

*        
*         RADIUS OF TANGENT POINT OF PATH, rt (<0 IF NON EXISTENT)
	  
	  call tangradpr (nlyr+1,rad,n_ref,c,RE,rt)
	  if(rt.lt.0.d0.and.zenang.gt.90.d0)then
*
*           *****************************************************              	  
*           light cannot reach targeted level (nor inferior ones)
*           => return	
*           *****************************************************

	    do ik=i,nlyr
	      NFAC(IK)=-1
	    enddo
	    
*           >>>>>>>	    
	    goto 99
*	    >>>>>>>

	  endif
*         
*	  
*         ******************************************
*         PATH plotting   
*         FROM TOP OF ATMOSPHERE DOWN TO AIMED POINT
*         ******************************************        
	 
	  j = 1
	  
*         sign=1.: before tangent point; -1.:after.	  
          sign  = 1.d0
	  itgpt = 0
	  
*         first atmospheric level to be reached	  
	  k     = 1
	  klev  = k
	  Rlev  = Re+VZ(k)
	  ilev  = 0
	  
	  SZALOC(i,0) = theta0*RadToDeg
	
*===         	
  3       continue
*===
            if(r0.ne.rt)then                
              r1 = r0 - ds*dcos(theta0)    
	    else
	      r1 = r2
	    endif    	         
	         
	    if (itgpt.eq.0 .and. r1.LT.rt) then
     
*             ************************    
*             TANGENT POINT IS REACHED
*             ************************

	      r1    = rt
	      theta1= 90.d0*DegtoRad
	      itgpt = 1
	      ktgpt = k
	      k     = k-1
	      if(k.le.nlyr)Rlev = Re+VZ(k)
	      sign  = -1.d0
	      ilev  = 0
	    endif
	    	
	    if( (itgpt.eq.0 .and. r1.le.Rlev).or.
     $          (itgpt.gt.0 .and. r1.ge.Rlev)) then

*              ******************************     
*              ATMOSPHERIC LEVEL K IS REACHED
*              ******************************
          
               ilev = 1
	       klev = k
	       
	       if(itgpt.eq.0)then
	         k = k+1
	       else
	         k=k-1
		 if(k.lt.0)then
	           write(*,*)'out of atm',i
		   
*                  no hope for i
*                  >>>>>>>		   
		   go to 5
*                  >>>>>>>
		   
		 endif
	       endif
	       
*              radius at reached level 	       
	       r1=Rlev
	       
*              next level to be reached       
	       if(k.le.nlyr) then
	          Rlev = Re+VZ(k)
	       endif
	       
	    endif
	   
	    	  
	    call indrefpr(nlyr+1, r1, Re,
     $              rad, n_ref, n1, hh, ratm1)
                       
	 
	    if(r1.ne.rt)then
	      if( 1.d0-(c/(r1*n1))**2 .LT.0.d0)then
	        write(*,*)'problem with DSQRT in geofac4new'
		write(*,*) i, r1*1.d-5, re*1.d-5
	        go to 5
	      endif 
	      theta1 = dacos( sign*DSQRT( 1.d0-(c/(r1*n1))**2 ) )
	    endif
	    
	      
	    x1       = -1.d0*r1*dcos(theta1)
	    dx       = x1-x0
	    dsdx1    = 1.d0/( 1.d0-ratm1*(dsin(theta1))**2 )
	    s1       = s0 + 0.5d0*( dsdx0+dsdx1 )*dx
            dbetadx1 = dsdx1* dsin(theta1)/r1
	    beta1    = beta0 
     $                 + 0.5d0*( dbetadx0+dbetadx1 )*dx



	    if(ilev.eq.1)then
	    
*             *****************************************************
*             LEVEL REACHED => RECORD GEOMETRY FACTOR AND LOCAL SZA
*             *****************************************************

	        if(itgpt.eq.0)then

*                 ********************		
*                 BEFORE TANGENT POINT
*                 ********************
		
	          nfac(i) = nfac(i)+1
	          if(klev.le.0.or.klev.gt.nlyr)then
	            write(*,*)'level problem!'
		    stop	
	          endif  
	          FAC(i,klev)    = s1/(VZ(klev-1)-VZ(klev))
		  SZALOC(i,klev) = theta1*RadToDeg
		 
		else 
		
*                 *******************		
*                 AFTER TANGENT POINT
*                 *******************
	         
		  itgpt           = itgpt+1
	          nfac_2(i)       = nfac_2(i)+1
	          FAC_2(i,klev+1) = s1/(VZ(klev)-VZ(klev+1))
		  SZALOC_2(i,klev)= theta1*RadToDeg
		  
	        endif
		
	        ilev= 0
	        s1  = 0.d0 
               
	    endif
	    
	    
	               
            if(beta1.ge.BETAGOAL)then
	    
*             *************************************
*             LOCAL VERTICAL IS REACHED (OVERTAKEN)
*             *************************************	          

*             Reached point on the local vertical
                         
	   
              Rreach = 
     $          r0*dsin(theta0)/dsin(theta0+BETAGOAL-beta0)
              if(Rreach.lt.0.D0)write(*,*)'R reached < 0 !?!?'
	     
	      
              deltas1= r1*(dsin(beta1-BETAGOAL)/
     $                       dsin(theta0+BETAGOAL-beta0) )
             
              theta2 = (theta1+theta0)*RadToDeg/2.d0
                  
	      if(dabs(Rgoal-Rreach).lt.precis)then	      	    
*          	              
*               ***********************************************
*               REACHED POINT ON THE VERTICAL IS OK, PATH IS OK
*               (distance from aimed point < precis)
*               ***********************************************
* 
*               A SZA>90 must be treated like a SZA<90 (isza90=1) 
*               IF:
*               There is no tangent point (itgpt=0)
*               OR
*               The tangent point is near the local vertical
*               (itgpt=1)

                isza90=0
                IF (itgpt.le.1)isza90=1 

*              
*               =============================================
*               tangent factor: from "fac_2" to "fac" factors 
*               =============================================        		
*                   
		if (nfac_2(i).gt.0.and.ktgpt.le.nlyr
     $               	.and. isza90.eq.0 ) then	  
    
		   fac(i,ktgpt)   = fac_2(i,ktgpt)
		   fac_2(i,ktgpt) = 0.d0
		   nfac(i)        = nfac(i)+1
		   nfac_2(i)      = nfac_2(i)-1
		   if(ktgpt.ge.1)szaloc(i,ktgpt) =  
     $  				  szaloc_2(i,ktgpt-1)
     
		endif

*               ===========================
*               factors near local vertical 
*               ===========================
		
                if (z_lay.ne.0.d0) then
                
*                 --------------------------		
*                 for midlayer-related paths
*                 --------------------------
		
		  if( zenang.le.90.d0 
     $                      .or.isza90.eq.1)then                        
		    nfac(i)        = nfac(i)+1
             	    fac(i,i)       = (s1-deltas1)/(VZ(i-1)-VZ(i))
		    SZALOC(i,i)    = theta2
		  else if(s1.ne.0.d0)then
		    if(fac_2(i,i).eq.0.d0)nfac_2(i)=nfac_2(i)+1
		    fac_2(i,i)     = (s1-deltas1)/(VZ(i-1)-VZ(i))
		    szaloc_2(i,i-1)= theta2
		  endif
		  
                else 
		
*                 -----------------------		
*                 for level-related paths
*                 ----------------------- 
                		
  		  if(zenang.le.90.d0 .or. isza90.eq.1)then
                    if(nfac(i).eq.i-1)then
                      nfac(i)      = nfac(i)+1
		      fac(i,i)     = (s1-deltas1)/(VZ(i-1)-VZ(i))
		      SZALOC(i,i)  = theta2
		    else if(nfac(i).ne.i)then
		      write(*,*)'problem'
		      stop 	  
		    endif
		  else
		    if(fac_2(i,i+1).eq.0.d0.AND.s1.NE.0.D0)then
                      nfac_2(i)    = nfac_2(i)+1
		      fac_2(i,i+1) = (s1-deltas1)/(VZ(i)-VZ(i+1))
		      szaloc_2(i,i)= theta2
		    endif    
		  endif 
		   	 
		endif
		
*               >>>>>>	    
	        goto 5
*               >>>>>>

              else
	      	
		DR = Rgoal-Rreach
			
		if(kount.eq.0)then
		  DR0 = DR
		  DZ0 = deltaz
		endif
			
*               >>>>>>
                goto 4
*               >>>>>>
              	   
	      endif
	      	      
	    endif
	    
	    
	               
	    
	    if(k.eq.nlyr+1)then

*             *****************	    
*             GROUND IS REACHED
*             ***************** 

              Rreach =  r1*dsin(theta1)/
     $                  dsin(theta1+betagoal-beta1)
             
	      
	      if( dabs(Rgoal-Rreach).gt.precis)then

 	      
*               DR (deltaz)	      
	        DR = Rgoal-Rreach
               			      
*               >>>>>>   
                goto 4
*               >>>>>>

              else

*                 REACHED POINT AND AIMED POINT ARE CLOSE ENOUGH
*                                   
*                          ===> NEXT AIMED POINT
	        
*               >>>>>>	      
	        goto 5
*               >>>>>>
		
	      endif
	      
	    endif
	  
	  
*           ***************************************
*           PROPERTIES OF CURRENT POINT OF THE PATH
*           ***************************************
            
	    r2       = r0	   	                 
	    r0       = r1
	    theta0   = theta1
	    s0       = s1
	    beta0    = beta1
	    ratm0    = ratm1
	    n0       = n1   
	    dsdx0    = dsdx1
	    dbetadx0 = dbetadx1
	    x0       = x1
	    

*           *****************************************
*           CONTINUE THE PATH (IF IT IS NOT TOO LONG)
*           *****************************************
	                  
	    j = j+1
	  
            if (j.le.nbmax*nbmax) then
 	     
*         <<<<<<	  
	  goto 3
*         <<<<<<
	    
	    else
	    
	      write(*,*)'too many iterations (geofac4)'         
	      stop
	        
	    endif
	  
 
   4      continue
   	  	  
*         **********************************************************	      
*         REACHED POINT ON THE VERTICAL IS TOO FAR FROM AIMED POINT!
*         **********************************************************

          kount = kount+1

*         *************************************
*         KOUNT NOT TOO LARGE => TRY A NEW PATH
*         *************************************
	  
	  if(kount.lt.80) then

              
              IF (IMETHOD.EQ.0) THEN
	      
	        deltaz = deltaz+DR
	       	
	      ELSE IF (IMETHOD.EQ.1) THEN

*               ************************************
*               SECANT METHOD [harsh Newton-Raphson]
*               for a root of function DR(deltaz) at
*               precision = 0.1 [m].
*               ************************************
 	
		if(isecant.eq.0)then
*                 test if secant method can start		
		  if(DR.ne.DR0)isecant=1
		endif
		
		if(isecant.eq.0)then
*                 secant method not possible yet		
	          deltaz   = deltaz+DR
	        else
*                 secant method
		  DR1      = DR
		  DZ1      = deltaz
		  DERIVDR  = (DR1-DR0)/(DZ1-DZ0)
		  if(DERIVDR.eq.0.d0)then
*                   secant method failed
                    kount  = 1000
		    deltaz = 0.d0
		    write(*,*)'secant method failed'		  
		  else
		    deltaz= deltaz-(DR1/DERIVDR)
		  endif
		  DR0      = DR1
		  DZ0      = DZ1
	        endif 
		
	      ELSE IF (IMETHOD.EQ.2) THEN
	      
*               ********************************	      
*               BISSECTION METHOD for DR(deltaz)
*               ********************************
	      
                if(ibissec.eq.0)then
*                 test if bissection method can start		
	          if(DR*DR0.lt.0.d0)ibissec = 1
	        endif
	      
	        if(ibissec.eq.0)then
*                 bissection method not possible yet		
	          deltaz   = deltaz+DR
	        else
*                 bissection method		
	          if(ibissec.eq.1)then
*                   starting bissection method		  
		    DR1    = DR
		    DZ1    = deltaz
		    ibissec= 2
		  else
		    if(DR*DR0.lt.0.d0)then
		      DR1  = DR
		      DZ1  = deltaz
		    else
		      DR0  = DR
		      DZ0  = deltaz
		    endif
		  endif
*                 Bissection		  
		  deltaz   = (DZ0+DZ1)/2.d0
*                 Regula-falsi	
*                 deltaz   = DZ1-DR1*((DZ1-DZ0)/(DR1-DR0))
	        endif
	        
	      ENDIF  
		 
*         <<<<<<<		  	  
          go to 2
*         <<<<<<<
		 
	  else if (IMETHOD.le.1)then
	  
	      IMETHOD = IMETHOD+1
	      write(*,*)'changing method for refraction effect...'
	      kount   = 0
	      deltaz  = 0.d0
	    
*         <<<<<<<		  	  
          go to 2
*         <<<<<<<
  
	      	      
	       
	  else    
	      
              write(*,*)' '
    	      write(*,*)
     $     'for chapman factors: numerical problem at ALTITUDE=',
     $                (RGOAL-re)*1.d-5, 'km'
              write(*,*)'reached pt is too far from aimed one: '
       	      write(*,*)'        ',DABS(rgoal-rreach)*1.d-2, ' m'
              write(*,*)' '
	      stop
	          
          endif
	  
  5       continue
          
	  		
      ENDDO
      
 99   return            
      end




                                           
******************************************************************************    
    
    
    
    
      subroutine POINTS_OF_INCIDENCE_PR 
     $           ( Z, Z0, zenang, nlyr, Re, ca )

*     points_of_incidence : CALCULATE POINTS OF INCIDENCE OF THE RAYS 
*                           INTO THE HIGHEST LAYER

      real*8 Z(0:*), Z0, zenang, ca
      real*8 hp       
*            horizontal path between point of incidence of a ray 
*            into the atmosphere and its height in the zenith    
      real*8 b, c, p1, p2, p3, Re, REA, PI
*      
      PI  = DACOS(-1.D0)     
            
*     Top of atmosphere                    
      REA = Re + Z(0)                
        
*     Elevation angle               
      c   = 90.d0 - zenang
      
      if ( c.NE.90.d0 ) then
     	b  = dtan(c*PI/180.d0)
        p1 = 1.d0 + b*b
        p2 = Re + Z0
        p3 = REA*REA*p1-p2*p2
	if(p3.LT.0.D0)then
	  write(*,*)'p3<0 in POINTS OF INCIDENCE2'
	  stop
	endif  
        if (p3 .GT. 0.d0) then
            hp = (-b*p2 + dsqrt(p3))/p1
	    if(dabs(hp/REA).GT.1.D0)then
	      write(*,*)'problem with hp in POINTS OF INC2'
	      stop
	    endif
            ca = dasin(hp/REA)
        else 
	    ca = -1.d0
        endif	
      else 
        ca = 0.d0
      endif	


      return
      end





************************************************************************
      SUBROUTINE tangradpr (ndim,rad,n_ref,c,RE,rt)
      
      real*8 rad(ndim), n_ref(ndim), c, rt
      real*8 limval, rad1, rad2, f1, f2, rarg,
     >       RE, ynew, hh, ratm, fnew
      
      limval = 1.D-30
      
      ind = 1
   4  f1 = (c/rad(ind))-n_ref(ind)
      f2 = (c/rad(ndim))-n_ref(ndim)
      rad1=rad(ind)
      rad2=rad(ndim)
      
      if ( f1*f2 .lt. 0.d0 ) then
     
*       Search for tangent radius by dichotomy

   5	rarg=(rad1+rad2)/2.d0
        call indrefpr(ndim, rarg, RE,
     >     rad, n_ref, ynew, hh, ratm)
        fnew=(c/rarg)-ynew
        if(dabs(fnew).LT.limval)then
	  rt=rarg
	else
	  if(fnew*f1.LT.0)then
	    rad2=rarg
	    f2=fnew
	  else
	    rad1=rarg
	    f1=fnew
	  endif
*         dichotomy method	  
	  goto 5
	endif
	        
      else
      
        ind=ind+1
	if(ind.lt.ndim)then
	  goto 4
	else
*         no tangent point      
          rt=-1.d0
	endif
	  
      endif
      
           
      RETURN
      END

***********************************************************************


      SUBROUTINE indrefpr( ndim, rarg, RE,
     $     rad, n_ref, ynew, hh, ratm)
*-------------------------------------------------------------------
*
*      Variation of inter for exponential interpolation of 
*      refraction index n and atmospheric parameter 
*                                             R_atm=-r/(n/(dn/dr)).
*   
*-------------------------------------------------------------------
*
*     Interpolates at the x-point rarg from x-value array rad and
*     y-value array n_ref. rad and n_ref are expected to have
*     descending arguments, i.e. for atmospheric applications
*     rad typically holds the altitude and inter expects
*     rad(1) = altitude at top of atmosphere.
*
*     Input variables:
*     dim       Array dimension of rad and n_ref
*     rarg      Interpolation argument
*     rad       array of radius values
*     n_ref     array of refraction index values
*     RE: radius of earth
*
*     Output variables:
*     ynew      Interpolated function value at rarg
*     ratm  
*     hh        gradient or scale height value  
*


      INTEGER ndim
      REAL*8  rarg, hh, rad(ndim), n_ref(ndim), 
     $        ynew, ratm, RE
*


      IF ( rarg .LE. rad(1) .AND. rarg .GE. rad(ndim)) THEN
         DO 10 iq = 1 , ndim-1
           IF ( rarg .LE. rad(iq) .AND. rarg .GT. rad(iq+1)) ip = iq
  10     CONTINUE
         IF ( rarg .EQ. rad(ndim)) ip = ndim - 1
      ELSEIF ( rarg .GT. rad(1)) THEN
         ip = 1
	 
      ELSEIF ( rarg .LT. rad(ndim)) THEN
         ip = ndim - 1
	
      ENDIF
*
*     Interpolate function value at rarg from data points ip to ip+1
*     (exponential interpolation)
*      
      IF ( n_ref(ip+1) .EQ. n_ref(ip) ) THEN
           
             hh   = 0.d0
             ynew = n_ref(ip)
             write(*,*)'problem with refraction index (indrefpr)', 
     $            ip, n_ref(ip), n_ref(ip+1)
*	     stop
	  	    
      ELSE          
             hh   = -( rad(ip+1) - rad(ip) ) / 
     $            DLOG( n_ref(ip+1) / n_ref(ip))
             ynew = n_ref(ip) * DEXP(- ( rarg - rad(ip) ) / hh)
	     ratm=rarg/hh
      ENDIF
*
     

      RETURN
      END       

      SUBROUTINE GEOFAST
     $            (brosza,sza_bro,fac,nfac,szaloc,
     $                   fac_2,NFAC_2,szaloc_2,
     $                 z_lay,nlyr,zd,zenang,r,VN )
*     
*     ********************************************************************
*
*     Calculates the geometric correction factor needed for spherical
*     geometry. 
*     WITH OR WITHOUT REFRACTION EFFECT. IF REFRACTION EFFECT, GEOFAST IS
*     LESS PRECISE BUT FASTER THAN GEOFACPR.    
*
*     INPUTS:
*            brosza: if =1, BrO density is function of alt and SZA.
*            z_lay: Where in the layer the Chapman function is to be
*                   computed. E.g. 0.0, bottom of layer, 0.5, middle
*                   of layer and 1.0 top of layer.
*            nlyr:  Number of layers in atmospheric model
*            zd(lc):lc = 0, nlyr. zd(lc) is distance from bottom 
*                   surface to top of layer lc. zd(nlyr) = 0.0 km
*            zenang:Solar zenith angle as seen from bottom surface
*            r:     Radius of earth. NOTE: Use the same dimension as zd,
*                   in km.
*            VN: refraction indexes of the layers
*
*     OUTPUTS:
*            nfac:
*            fac:   geometrical correction factor ds/dh
*            szaloc: local solar zenithal angles
*            altloc: local altitudes
*            fac_2: ds/dh beyond the tangent point of the path
*            szaloc_2: local sza beyond the tangent point
*            altloc_2: local altitudes beyond the tangent point
*           
*     
*     ********************************************************************  
*
*
      INCLUDE 'DISORT.MXD'
      integer brosza, nlyr, nfac(1:mxcly), NFAC_2(1:mxcly) 
      real*8  z_lay, zd(0:*), zenang, r, 
     $        VN(0:mxcly), VNi(0:mxcly), 
     $        VN2(0:mxcly), VZ(0:mxcly),
     $        dsdh1(0:mxcly-1,0:mxcly-1),
     $        dsdh2(0:mxcly-1,0:mxcly-1),  
     $        SZALOC1(0:mxcly-1,0:mxcly), 
     $        SZALOC2(0:mxcly-1,0:mxcly), 
     $        SZALOC(1:mxcly,0:mxcly), 
     $        SZALOC_2(1:mxcly,0:mxcly), 
     $        fac(1:mxcly,1:mxcly), fac_2(1:mxcly,1:mxcly),
     $        sza_bro(mxsza)
     
     
*
*     VZ: altitudes in cm and in increasing order
*     VNi: refraction indexes corresponding to VZ    
*      
     
      
      do i=0,nlyr
	VZ(i)=zd(nlyr-i)*1.d+5
	VNi(i)=VN(nlyr-i)
      enddo
      
*
*     Earth radius: km --> cm
*     
      r=r*1.d+5
*
      do i=1,nlyr
        nfac(i)=0
	NFAC_2(i)=0
	do j=1,nlyr
	  fac(i,j)=0.D0
	  fac_2(i,j)=0.D0
	enddo
      enddo      
      
     
      
      do i=0, nlyr-1
        do j=0,nlyr-1
	  dsdh1(i,j)=0.D0
	  dsdh2(i,j)=0.D0
	  SZALOC1(i,j)=0.D0
	  SZALOC2(i,j)=0.D0
	enddo
	SZALOC1(i,nlyr)=0.D0
	SZALOC2(i,nlyr)=0.D0
      enddo   
*     
           
      CALL opathfast(zenang,VZ,VNi,VN2,r,nlyr,z_lay,dsdh1,
     $           SZALOC1,dsdh2,SZALOC2, NFAC,NFAC_2)
      
     
*                                                      
* |-------> FAC, FAC_2, SZALOC, SZALOC_2, ALTLOC, ALTLOC_2  
*  
      do i=1, nlyr
        do j=1, nlyr
	  FAC(i,j)=dsdh1(nlyr-i,nlyr-j)
	  FAC_2(i,j)=dsdh2(nlyr-i,nlyr-j)
	 
	
	 
	  szaloc(i,j)=szaloc1(nlyr-i,nlyr-j)
	  szaloc_2(i,j)=szaloc2(nlyr-i,nlyr-j)
*                                          
	  if(j.eq.nfac(i)-nfac_2(i)-1)then 
	    szaloc_2(i,j)=zenang
	  endif
*          	  
	  if(brosza.eq.1)then
*           If local sza is lower than the minimum sza of BrO matrix,
*           put its value to that minimum	  
	    if(szaloc(i,j).lt.sza_bro(1))
     $	            szaloc(i,j)=sza_bro(1)
	    if(szaloc_2(i,j).lt.sza_bro(1))
     $	            szaloc_2(i,j)=sza_bro(1)
	  endif
	enddo
*		
	szaloc(i,0)   = szaloc1(nlyr-i,nlyr)
	szaloc_2(i,0) = szaloc2(nlyr-i,nlyr)
*           
        if(brosza.eq.1)then
          if(szaloc(i,0).lt.sza_bro(1))
     $      szaloc(i,0)=sza_bro(1)
	  if(szaloc_2(i,0).lt.sza_bro(1))
     $      szaloc_2(i,0)=sza_bro(1)                            	
	endif
      enddo
      
     
      return 	   
      end
      
      
      
   
      

      SUBROUTINE opathfast 
     $     (khi,VZ,VN,VN2,Re,nlyr,z_lay,dsdh1,SZALOC1,
     $      dsdh2,SZALOC2,NFAC,NFAC_2)
*
*
*    ***********************************************************   
*     Compute optical path and local SZA in successive 
*     atmospheric layers Refraction is taken into account.
*
*     INPUTS:
*       khi	: SZA out of atmosphere in degrees
*       VZ	: Altitude vector [CM]  
*       VN	: refraction index vector
*       nlyr    : number of layers of atmosphere
*       Re:       radius of earth [CM]
*
*     OUTPUTS:
*       dsdh: ds/dh geometrical factors
*       SZALOC: local SZAs
*       ALTLOC: local altitudes
*          index 1 : before tangent point (included)
*          index 2 : after tangent point
*    ***********************************************************
*     
*      
      PARAMETER ( NMAX = 5)
      INCLUDE 'DISORT.MXD'
*     Variables XX and YY below must be XX(0:NMAX),YY(0:NMAX)      
      INTEGER  NL, ind, nlyr, cptr, flag_secant, scratch,
     $         NFAC(1:mxcly), NFAC_2(1:mxcly)
      REAL*8   PI, DegToRad, khi, deltaz_top, tmp_top, 
     $         A, B, Zi, VZ(0:mxcly),  deltaz_bot, ca, ia,
     $         VN(0:mxcly), n1, n2, ra, ztop, zbot, r1, r2, Re, g,
     $         local_sza_bot, gamma1, path,tanheight,
     $         f, gamma2, new_bot, old_bot, old_top, new_top,
     $         Sx,Sx2,Sxy,Sy,Delta, local_sza_top,
     $         XX(0:5), yy(0:5),VN2(0:mxcly), z_lay, yp1,ypn,
     $         dsdh1(0:mxcly-1,0:mxcly-1), SZALOC1(0:mxcly-1,0:mxcly),
     $         dsdh2(0:mxcly-1,0:mxcly-1), SZALOC2(0:mxcly-1,0:mxcly)
*
*       
                                           
      PI=DACOS(-1.D0)
      DegToRad=PI/180.D0
      RadToDeg= 1/DegToRad
*     
*     ray tracing calculated only for slant path geometry
      if ( khi.eq.0.D0 ) then
        write(*,*)'no ray tracing for vertical path'
	return
      endif	    

*     Initialisations:                        

      deltaz_top = 0.d0
      old_top = 0.d0
      tmp_top = 0.d0
      NL = 0
      ind = 0  
      A = 0.D0
      B = 0.D0
*
*     ************************************************************
*     ************************************************************
*     Scan scattering altitudes from top down to bottom atmosphere
*     ************************************************************
*     ************************************************************
*   
      DO i=nlyr-1,0,-1
*            
	Zi=VZ(i)+z_lay*(VZ(i+1)-VZ(i))
*	
        if (NL.gt.2) deltaz_top = -dexp( A + B*Zi )
*        
        cptr = 0 
        deltaz_bot = 0.0 
        flag_secant = 0

*       ----------------------------------------------------------
*       loop until final point is close enough (less or equal 1 m)
*       from the aimed point
*       ----------------------------------------------------------
       
  2     nfac(nlyr-i)=0
        NFAC_2(nlyr-i)=0
        cptr=cptr+1
        if ( cptr.gt.15 )  go to 3

	call points_of_incidence_fast
     $       (VZ,Zi-deltaz_top,khi,nlyr,Re,ca)
	
	ia = khi*DegToRad - ca
	local_sza_top = ia*RadToDeg
	
* |--->	szaloc1
	SZALOC1(i,nlyr)=local_sza_top

	
	n2 = VN(nlyr)

	
*       	
*       *****************************************	
*       ray traycing : scans atmosphere from top 
*       down to tangent point or scattering level
*       *****************************************
*      
      
 
        DO j=nlyr-1,-1,-1
*    
	    if(j.lt.0) then
	      n1= VN(j+1) 
	    else
	      n1= VN(j)
	    endif
	   
     
            if(dabs(n2/n1*dsin(ia)).le.1.D0)then
              ra = dasin( n2/n1 * dsin(ia) )
            else if(dabs((n2/n1*dsin(ia))-1.d0).lt.1.d-15)then
              ra=dasin(1.D0)
	    else
	    
	      go to 25
            endif
	    	   
            ztop = VZ(j+1)
	   
	    
*	    add a 2 km thick layer below the atmosphere
            
            if(j.lt.0) then
	      zbot=(ztop-2.d5)
	    else
	      zbot=VZ(j)  
	    endif
	    
            r1 = Re + ztop
            r2 = Re + zbot
            g = r1/r2 * dsin(ra)
*           Precision correction	    
	    if(dabs(g-1.d0).lt.1.d-15)g=1.D0
* 
          
            if ( g .lt. 1.d0 ) then
             
*             ===================	      
*   	      above tangent point
*             ===================

              ia = dasin(g)
              local_sza_bot = ia*RadToDeg
              gamma1 = ia - ra
*            
                  
              call optical_path_fast
     $             (re,ra,ca,gamma1,r1,r2,local_sza_bot,path)
              
*                             
              if ( path.eq. -1.d0) then 
*	         scattering layer is reached
                 go to 21 
              endif      
*	                 
              if ( j.ge.0 ) then	   

* |----------------> I. dsdh1              
                     dsdh1(i,j)=path/(VZ(j+1)-VZ(j))
		     SZALOC1(i,j)=local_sza_bot
		     nfac(nlyr-i)=nfac(nlyr-i)+1
              endif              
            else

*             ========================	            
*             tangent point is reached
*             ========================
              
              scratch = 0 
*
*             
	      tanheight = r1*dsin(ra)-Re
*	      
	      if(khi.gt.90.AND.tanheight.lt.0.d0) then
		
*                tangent point under earth surface, ground has been 
*                reached, light cannot reach targeted level (inferior ones
*                neither)
*              
*           		       
		 do ik=nlyr-i,nlyr
		   nfac(ik)=-1
		 enddo
*
*                exit loop i	       
	         GOTO 27
	      endif
*             	      
	                  	
              if ( j.ge.0 ) then
	                        
*               avoid discontinuity at tangent height 
   23		tanheight = r1*dsin(ra)-Re
              
		if ( tanheight.eq.r2-Re )  goto 231
	       
	        yp1=1.D+31
	        ypn=yp1    
		
	 	call spline2 (VZ,VN,nlyr+1,yp1,ypn,VN2)
		call splint2 (VZ,VN,VN2,nlyr+1,tanheight,n1)
		
		
		ra = dasin(n2/n1*dsin(ia))
		scratch=scratch+1
		if( scratch.lt.6 ) then
		  goto 23
	        endif
              endif
  231         r2 = r1
              	
              local_sza_bot = (PI-ra)*RadToDeg
              gamma1 = PI - 2.D0*ra
           
	        	      	
              call optical_path_fast
     $             (re,ra,ca,gamma1,r1,r2,local_sza_bot,path)
                             
*             ==========================
*             Layer of the tangent point
*             ==========================
	      			          
              if ( j.ge.0 ) then

* |-----------------> II. dsdh1
                      dsdh1(i,j)=path/(VZ(j+1)-VZ(j))

		      
* |------>	      szaloc:	      
		      szaloc1(i,j)=local_sza_bot
		    
		      nfac(nlyr-i)=nfac(nlyr-i)+1            
              endif 
              	              
              if ( ca.gt.1.D-10 )then
*                  ca > 0 with a certain precision
                  
*        	scans atmosphere up to scattering level, 
*               starting from first layer above tangent point.

                ia = ra
              	     
              	do k=j+1,nlyr
		  if(k.eq.nlyr)then
		    n2 = VN(k-1)
		  else   
		    n2 = VN(k)
		  endif  
              	  f = n1/n2*dsin(ia)
              	  
		 
              	  if ( f .le. 1.D0 )then
                    ra = dasin(f)
		    if(k.eq.nlyr) then
*		      add a 3 km thick layer above atmosphere
		      ztop=VZ(k)+3.d5
		    else
		      ztop=VZ(k+1)
		    endif
                    zbot = VZ(k)
                    r1 = Re + zbot
                    r2 = Re + ztop
                    ia = dasin(r1/r2*dsin(ra))
                    gamma2 = ra - ia
                    ra = PI - ra
                    local_sza_top = (PI-ia)*RadToDeg  
                   
                    call optical_path_fast
     $                 (re,ra,ca,gamma2,r1,r2,local_sza_top,path)
                  
*                   
*                   ====================
*                   Beyond tangent point
*                   ====================

                    if ( k.lt.nlyr )then
* |-----------------> III. dsdh2
                      dsdh2(i,k)=path/(VZ(k+1)-VZ(k))
		      szaloc2(i,k)=local_sza_bot
		      szaloc2(i,k+1)=local_sza_top
		      if(path.gt.0.d0) NFAC_2(nlyr-i)=NFAC_2(nlyr-i)+1 				 
		    endif
		    
		  endif	    
                            
                  n1 = n2
                  local_sza_bot = local_sza_top
                         
                  if (ca .lt. 0.) then
		    goto 24
		  else 
*    		    The ray leaves the atmosphere before reaching 
*                   direction of looking  ???
		  endif 	  
                enddo 
  24		continue     
                 
              endif   
            	  
              GO TO 21  
            endif

    
    
            n2 = n1
            local_sza_top = local_sza_bot
	    

*
*           ================================
*           exit test ======================   
            if ( ca .lt. 1.d-7 ) GOTO 21
*           ==================== exit loop j
*           ================================	    

   25     continue  
        ENDDO
	
*      	
*       *************	          
*       End of loop j        
*       *************
*     
	
	  
   21   continue
   
        if ( ca .gt. 0.5*gamma1 )then
                   
*   
*            ======================================   
*            the ray hit the ground before reaching  
*            direction of looking !!!!!!!!!!!!!!
*            ======================================
*             
*                
            
	  go to 3
         
        endif
	
	
        new_bot = r2-Re-Zi

	
        old_top = tmp_top
        old_bot = deltaz_bot
	 
        deltaz_bot = new_bot
        
*       if bad convergence => turn to secant method  
           
        if ( (cptr.gt.1) .and. 
     $   (dabs(deltaz_bot) .gt. dabs(old_bot)*0.5D0) )
     $	                     flag_secant = 1
        
        if ( flag_secant.eq.1 )then    
*	  (x1.y2 - y1.x2)/(y2-y1)
          new_top = (old_top*deltaz_bot-old_bot*deltaz_top) 
     $                  / (deltaz_bot-old_bot)
        else 
          new_top = deltaz_top + deltaz_bot
        endif  
       
          
        tmp_top = deltaz_top
        deltaz_top = new_top
	       
        
*       --------------------------------------	  
        IF (dabs(deltaz_bot).gt. 1.d2 ) goto 2
*       --------------------------------------	       

   3    continue
       
        if ( deltaz_top .lt. -10.d0 ) then
	
	 
*          Linear regression on NMAX last points to extrapolate
*          deltaz_top at next layer position
          
                    
           NL=NL+1  
	   IF(NL.gt.NMAX)then
	     N=NMAX
	   else
	     N=NL
	   endif		       
	    
	   if(deltaz_top.ge.0)then
	     write(*,*)'error in geofac4 - opathfast'
	     stop
	   endif
	   
	   XX(ind) = Zi 
	   yy(ind) = dlog(- deltaz_top)   		    
          
	   Sx = 0.d0 
	   Sx2 = 0.d0
	   Sy = 0.d0
	   Sxy =  0.d0
	    
	   do k=0, N-1
		Sx  = Sx+XX(k)
		Sx2 = Sx2 + XX(k)*XX(k)
		Sy  = Sy + yy(k)
		Sxy = Sxy + XX(k)*yy(k)
	   enddo
	   
	   Delta = N*Sx2-Sx*Sx 
	     
	   if (Delta.ne.0.D0)then
		A = (Sx2*Sy - Sx*Sxy)/Delta
		B = (N*Sxy-Sx*Sy)/Delta
	 
	   endif  
	   
	   if(ind.eq.NMAX-1)then
	     ind=0
	   else
	     ind=ind+1
	   endif
         
        endif
  26    continue	
        
      ENDDO

*     **************************       
*     ************************** 
*     End of loop i=nlyr-1,-1,-1          
*     **************************
*     **************************
      
  27  continue  
      
      return
      end





      subroutine points_of_incidence_fast 
     $           ( Z, Z0, khi, nlyr, Re, ca )

*     points_of_incidence_fast : CALCULATE POINTS OF INCIDENCE OF THE RAYS 
*                           INTO THE HIGHEST LAYER

      real*8 Z(0:*), Z0, khi, ca
      real*8 hp       
*            horizontal path between point of incidence of a ray 
*            into the atmosphere and its height in the zenith    
      real*8 b, c, p1, p2, p3, Re, REA, PI
*      
      PI=DACOS(-1.D0)     
            
*     Top of atmosphere                    
      REA = Re + Z(nlyr)                
        
*     Elevation angle               
      c = 90.d0 - khi
      
      if ( c.NE.90.d0 ) then
     	b  = dtan(c*PI/180.d0)
        p1 = 1.d0 + b*b
        p2 = Re + Z0
        p3 = REA*REA*p1-p2*p2  
        if (p3 .GT. 0.d0) then
            hp = (-b*p2 + dsqrt(p3))/p1
            ca = dasin(hp/REA)
        else 
	  ca = -1.d0
        endif	
      else 
        ca = 0.d0
      endif	


      return
      end



 

      subroutine optical_path_fast
     $           (re,ra,ca,gamma1,r1,r2,local_sza,path)
      
*     optical_path_fast : Calculates optical path in meters in each layer crossed. 
*                    A ray reaching the direction of looking will be cut.     
 
      real*8 ra, ca, gamma1, r1, r2, local_sza, path
      real*8 rd             
*              Distance between the point of looking and 
*              the point where the ray enters the layer  
      real*8 rl            
*              Distance between the point of looking and the   
*              point where the ray meets the looking direction 
      real*8 ra1           
*              Angle between rd and r1 
      real*8 p1,p2,p3
      real*8 LA        
      real*8 Z0, SCA, a, PI, RadToDeg, re
*      
*
*     RE: Earth radius in cm

      
      PI=DACOS(-1.D0)
      RadToDeg=180.D0/PI
      LA=0.D0
      Z0=0.D0
      
      if ( ca .LT. 0.d0 ) then
        path=-1.d0            
        return
      endif
      
      ca = ca-gamma1
      
      if ( ca .LT. 0.d0 .AND. DABS(ca).gt.1.D-10) then 
*                             (avoid false ca<0)           
*                     If the ray beames over the direction of looking,
*                     it is cut there
     
        ca = ca + gamma1
	               
        SCA = PI + LA - ca - ra         
*            scattering angle = PI-apparent_SZA       
        local_sza = (PI - SCA)*RadToDeg
	
        p1  = (Re+Z0)*(Re+Z0) + r1*r1 - 2.d0*(Re+Z0)*r1*dcos(ca)
        r2 = r1*dsin(ra)/dsin(PI-ra-ca)


        if ( r2 .GT. (Re+Z0-1.d2) ) then    
*                         allow 1 m tolerance  
       	  if (p1.GT.0.d0) then
             	rd  = dsqrt(p1)
             	a = ((Re+Z0)/rd)*dsin(ca)   
             	if ( dabs(a).GT. 1.D0 ) then  
		  ca = ca-gamma1  
		  path=0.d0
		  go to 99 
             	endif
                ra1 = dasin(a)
                rl  = rd*dsin(ra-ra1)/dsin(SCA)
                p2 = rl*rl + rd*rd - 2.d0*rl*rd*dcos(ca+ra1-LA)
                if (p2.GT.0.0) then
		  path = dsqrt(p2)
*		  CUT PATH THROUGH LAST LAYER		    
                else 
		  path = 0.d0
		endif  
                p3 = rl*rl+(Re+Z0)*(Re+Z0)-2.d0*rl*(Re+Z0)*dcos(PI-LA)
                if (p3.GT.0.d0) r2 = dsqrt(p3)
                ca = ca-gamma1
          else  
                ca = -1.d0
                r2 = Re + Z0
                path = 0.d0
          endif		
            
        else                       
*         if scattering height is too low (even with 1m tolerance!) 
          ca = -1.d0
          path = 0.d0
        endif
     
      else             
*                     
*       FULL PATH THROUGH A LAYER
* 
        p1 = r1*r1 + r2*r2 - 2.d0*r1*r2*dcos(gamma1)
        if (p1.GT.0.d0) then
	  path = dsqrt(p1)
        else 
	  path = 0.d0
	endif
	
      endif
 99   continue  
      return
      end 
   

      subroutine chpman2(NLYR,zenang,DTAUC,
     $                    NFAC,FAC,
     $                    NFAC_2,FAC_2,
     $                    CHP, mxcly)
      
*     *****************************************************
*
*     CHAPMAN FUNCTION, with refraction taken into account
*
*     INPUTS: 
*       nlyr: number of atmospheric layers.
*       dtauc: optical depth of layers. 
*       fac,nfac,fac_2,NFAC_2: outputs of GEOFAC4.
*   
*     OUTPUT:
*       CHP: Chapman factors
*     
*     ***************************************************** 

      INTEGER nlyr, nfac(1:mxcly), NFAC_2(1:mxcly) 
      REAL*8  fac(1:mxcly,1:mxcly), fac_2(1:mxcly,1:mxcly), 
     $        dtauc(*), zenang

      REAL*8  CHP(1:mxcly)
  
     

      do lc=1, nlyr
        
	CHP(lc)=0.D0
	     
*       *****************************************	
*       before tangent point (included if exists)
*       *****************************************	

	if(nfac(lc).gt.0)then                 	
	  do j=1,nfac(lc)
	    CHP(lc)=CHP(lc)+FAC(lc,j)*
     $                      (dtauc(j))
	  enddo
	else
	  if(zenang.gt.90.D0) CHP(lc)=1.D+20  
	endif
	
*       *******************	
*       after tangent point
*       *******************
      
	if(NFAC_2(lc).gt.0)then
	  do j=1,NFAC_2(lc)
	  
	    CHP(lc)=CHP(lc)+FAC_2(lc,nfac(lc)-j)*
     $                         (dtauc(nfac(lc)-j))
	  enddo
	endif
	
      enddo    
     
      
      return
      end
      
      



      SUBROUTINE DTBRO (zd,zenang,sigbro,densbro,sza_bro,
     $                  NFAC, NFAC_2, 
     $                  SZALOC, SZALOC_2, 
     $     nlyr, ndenssza,
     $                  DTAUBRO, DTAUBRO_2,DBRO,DBRO_2,
     $                  VDBR)
     
*
*     ***************************************************
*
*     Calculates the BrO optical depth of layers, for  
*     SZA-dependent BrO densities.  This in order to 
*     calculate the Chapman function. 
*
* 
*     OUTPUTS:
*         DTAUBRO: BrO optical depth of layers before the
*                  tangent point of the rays.
*         DTAUBRO_2: BrO optical depth of layers after 
*                    the tangent point of the rays.            
*
*     *************************************************** 
*    
      INCLUDE 'DISORT.MXD'
      integer nlyr, nfac(1:mxcly), NFAC_2(1:mxcly), ndenssza 
      REAL*8 zd(0:*), sigbro(0:*),DENSBRO(mxsza,0:*),
     $       SZALOC(1:mxcly,0:mxcly),
     $       SZALOC_2(1:mxcly,0:mxcly),
     $       sza_bro(mxsza), dbro(1:mxcly,0:mxcly),
     $       dbro_2(1:mxcly,0:mxcly),zenang,
     $       hh
      REAL*8 DTAUBRO(mxcly,*), DTAUBRO_2(mxcly,*)      
      REAL*8 sord(MXsza),sabs(MXsza)
      real*8 pi,VDBR(mxcly)
      logical expon
* 
    
      
*     choice of exponential interpolation      
      expon=.true.
*      
      pi = 2.D0 * DASIN(1.D0)
* 
*     ************************     
*     BEFORE THE TANGENT POINT (maybe included)
*     ************************
*

     
      do lc=1, nlyr
       
        if(nfac(lc).gt.0)then
	  
	  do j=1, nfac(lc)
	    
*           interpolation of BrO densities
            do m=1,ndenssza
	      sabs(m)=sza_bro(ndenssza-m+1)
	      sord(m)=densbro(ndenssza-m+1,j-1)
	    enddo
	    	   
            CALL dpinter(mxsza,ndenssza,1,szaloc(lc,j-1),
     $                sabs,sord,dbro(lc,j-1),hh)
	    		    
	    
              
     
                                  
*           interpolation of BrO densities
            do m=1,ndenssza
	        sord(m)=densbro(ndenssza-m+1,j)
		if(j.EQ.nfac(lc).AND.zenang.gt.90.D0)
     $                 sord(m)=densbro(ndenssza-m+1,j-1)
	    enddo 
	     	    
	    CALL dpinter(mxsza,ndenssza,1,szaloc(lc,j),
     $                sabs,sord,dbro(lc,j),hh)    
	      
                  
        
	 linear =0
            IF(expon .EQV. .false.  .or.
     $	       SIGBRO(j-1)*dbro(lc,j-1).EQ.0.D0 .or.
     $         SIGBRO(j)*dbro(lc,j).EQ.0.D0.or.
     $         SIGBRO(j-1)*dbro(lc,j-1).EQ. 
     $         SIGBRO(j)*dbro(lc,j)) THEN
C             LINEAR 	    
               linear=1
	      dtaubro(lc,j)=(zd(j-1)-zd(j))*1.D+5
     $                    *0.5D0*( SIGBRO(j-1)*dbro(lc,j-1)
     $	                          +  SIGBRO(j)*dbro(lc,j) )			  
            ELSE
C             EXPON	    
	      dtaubro(lc,j)=(zd(j-1)-zd(j))*1.D+5
     $                    *( SIGBRO(j-1)*dbro(lc,j-1)
     $	                  -  SIGBRO(j)*dbro(lc,j) )/
     $                    DLOG( (SIGBRO(j-1)*dbro(lc,j-1))/
     $                          (SIGBRO(j)*dbro(lc,j))  )
            ENDIF
c	    endif 
	    
	  enddo
	endif
       enddo 
 
 
       
          
*     
*     *********************** 
*     AFTER THE TANGENT POINT 
*     ***********************
*
      do lc=1, nlyr
      
        if(NFAC_2(lc).gt.0)then
	
	    	
	  do j=1, NFAC_2(lc)	    
*           interpolation of BrO densities
            do m=1,ndenssza
	      sabs(m)=sza_bro(ndenssza-m+1)
	      sord(m)=densbro(ndenssza-m+1,nfac(lc)-j)
	    enddo 	
	    CALL dpinter(mxsza,ndenssza,1,szaloc_2(lc,nfac(lc)-j),
     $                sabs,sord,dbro_2(lc,nfac(lc)-j),hh)    
	      
	
              
*             interpolation of BrO densities
              do m=1,ndenssza
	        sord(m)=densbro(ndenssza-m+1,nfac(lc)-j-1)
	      enddo 
	      CALL dpinter(mxsza,ndenssza,1,szaloc_2(lc,nfac(lc)-j-1),
     $                sabs,sord,dbro_2(lc,nfac(lc)-j-1),hh)        
	   
               
	   
  

            IF(expon .EQV. .false. .or.
     $       SIGBRO(nfac(lc)-j-1)*dbro_2(lc,nfac(lc)-j-1).EQ.0.D0
     $       .OR. SIGBRO(nfac(lc)-j)*dbro_2(lc,nfac(lc)-j).EQ.0.D0
     $       .or. SIGBRO(nfac(lc)-j-1)*dbro_2(lc,nfac(lc)-j-1).EQ. 
     $         SIGBRO(nfac(lc)-j)*dbro_2(lc,nfac(lc)-j) ) THEN
C             LINEAR 	    
	      dtaubro_2(lc,nfac(lc)-j)=(zd(nfac(lc)-j-1)
     $                               -zd(nfac(lc)-j))*1.D+5
     $    *0.5D0*(SIGBRO(nfac(lc)-j-1)*dbro_2(lc,nfac(lc)-j-1)
     $	      +  SIGBRO(nfac(lc)-j)*dbro_2(lc,nfac(lc)-j) )			  
            ELSE
C             EXPONENTIAL	    
	      dtaubro_2(lc,nfac(lc)-j)=(zd(nfac(lc)-j-1)
     $                               -zd(nfac(lc)-j))*1.D+5*
     $        ( SIGBRO(nfac(lc)-j-1)*dbro_2(lc,nfac(lc)-j-1)
     $	      -  SIGBRO(nfac(lc)-j)*dbro_2(lc,nfac(lc)-j) )
     $       /DLOG( (SIGBRO(nfac(lc)-j-1)*dbro_2(lc,nfac(lc)-j-1))
     $             /(SIGBRO(nfac(lc)-j)*dbro_2(lc,nfac(lc)-j)) )
            ENDIF
	    
	     		  

	    
	  enddo
	endif
	
      enddo           
     
     

  
      return      
      
      END 


       SUBROUTINE dpinter( dim, npoints, itype, arg,
     $     xarr, yarr, ynew, hh)
*-------------------------------------------------------------------
* Copyright (C) 1994 Arve Kylling
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 1, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* To obtain a copy of the GNU General Public License write to the
* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
* USA.
*-------------------------------------------------------------------
*
*     Interpolates at the x-point arg from x-value array xarr and
*     y-value array yarr. xarr and yarr are expected to have
*     descending arguments, i.e. for atmospheric applications
*     xarr typically holds the altitude and inter expects
*     xarr(1) = altitude at top of atmosphere.
*
*     Input variables:
*     dim       Array dimension of xarr and yarr
*     npoints   No. points in arrays xarr and yarr
*     itype     Interpolation type
*     arg       Interpolation argument
*     xarr      array of x values
*     yarr      array of y values
*
*     Output variables:
*     ynew      Interpolated function value at arg
*     hh        gradient or scale height value  
*
      INTEGER dim
      REAL*8  arg, hh, xarr(dim), yarr(dim), ynew
*
       IF ( arg .LE. xarr(1) .AND. arg .GE. xarr(npoints)) THEN
         DO 10 iq = 1 , npoints-1
           IF ( arg .LE. xarr(iq) .AND. arg .GT. xarr(iq+1)) ip = iq
   10    CONTINUE
         IF ( arg .EQ. xarr(npoints)) ip = npoints - 1
       ELSEIF ( arg .GT. xarr(1)) THEN
         ip = 1
       ELSEIF ( arg .LT. xarr(npoints)) THEN
         ip = npoints - 1
       ENDIF
*
*     Interpolate function value at arg from data points ip to ip+1
*
*
*     exponential interpolation
*
       IF ( itype .EQ. 1 ) THEN
          IF ( yarr(ip+1) .EQ. yarr(ip) ) THEN
             hh   = 0.0
             ynew = yarr(ip)
          ELSE
             hh   = -( xarr(ip+1) - xarr(ip) ) / 
     $            LOG( yarr(ip+1) / yarr(ip))
             ynew = yarr(ip) * EXP(- ( arg - xarr(ip) ) / hh )
          ENDIF
*
*     linear interpolation
*
       ELSEIF ( itype .EQ. 2 ) THEN
         hh = ( yarr(ip+1) - yarr(ip) ) / ( xarr(ip+1) - xarr(ip) )
         ynew = yarr(ip) + hh*( arg - xarr(ip) )
*     
      ENDIF
      RETURN
      END       


      subroutine chapman4(NLYR,zenang,
     $                    NFAC,FAC,DTAUC_MB,DTAUBRO,
     $                    NFAC_2,FAC_2,DTAUBRO_2, mxcly,
     $                    CHP)
      
*     *****************************************************
*
*     CHAPMAN FUNCTION, with refraction taken into account
*
*     INPUTS: 
*       nlyr: number of atmospheric layers.
*       dtauc_mb: output of OPTICP, optical depth for
*                 atmosphere without BrO. 
*       fac,nfac,fac_2,NFAC_2: outputs of GEOFAC4.
*       dtaubro, dtaubro_2: outputs of DTBRO, BrO term of 
*                           optical depths.
*     OUTPUT:
*       CHP: Chapman factors
*     
*     ***************************************************** 

      INTEGER nlyr, nfac(1:mxcly), NFAC_2(1:mxcly) 
      REAL*8  DTAUBRO(mxcly,*), DTAUBRO_2(mxcly,*),
     $        fac(1:mxcly,1:mxcly), fac_2(1:mxcly,1:mxcly), 
     $        dtauc_mb(*), zenang

      REAL*8  CHP(1:mxcly)

     

      do lc=1, nlyr
        
	CHP(lc)=0.D0
	     
*       *****************************************	
*       before tangent point (included if exists)
*       *****************************************	

	if(nfac(lc).gt.0)then
                 	
	  do j=1,nfac(lc)
	    CHP(lc)=CHP(lc)+FAC(lc,j)*
     $              (dtauc_mb(j) + dtaubro(lc,j))
	  enddo
         
	else
	  if(zenang.gt.90.D0) CHP(lc)=1.D+20  

	endif
	
*       *******************	
*       after tangent point
*       *******************
      
	if(NFAC_2(lc).gt.0)then
	  do j=1,NFAC_2(lc)
	   
	    CHP(lc)=CHP(lc)+FAC_2(lc,nfac(lc)-j)*
     $      (dtauc_mb(nfac(lc)-j) + dtaubro_2(lc,nfac(lc)-j))
           
	  enddo

	endif
	
      enddo    
     
      
      return
      end
      
      
      SUBROUTINE  dpSGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )
C
C         FACTORS A REAL*8 BAND MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.
C
C     INPUT:
C
C        ABD     REAL*8(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL*8
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL*8(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     EXAMPLE:  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
C
C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C
C     ROUTINES CALLED:  FROM LINPACK: SGBFA
C                       FROM BLAS:    SAXPY, DPSDOT, SSCAL, DPSASUM
C                       FROM FORTRAN: ABS, dmax1, MAX0, MIN0, SIGN
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER  LDA, N, ML, MU, IPVT(*)
      REAL*8     ABD(LDA,*), Z(*)
      REAL*8    RCOND
C
      REAL*8     DPSDOT, EK, T, WK, WKM
      REAL*8     ANORM, S, DPSASUM, SM, YNORM
      INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
C
C
C                       ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM, DPSASUM(L,ABD(IS,J), 1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C                                               ** FACTOR
      CALL dpSGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                     ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K))/ABS(EK-Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (ABD(M,K) .NE. 0.0E0) THEN
            WK  = WK /ABD(M,K)
            WKM = WKM/ABD(M,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
         MM = M
         IF (KP1 .LE. JU) THEN
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
C
C                         ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML, N-K)
         IF (K .LT. N) Z(K) = Z(K) + DPSDOT(LM,ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
C
      YNORM = 1.0E0
C                         ** SOLVE L*V = Y
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML, N-K)
         IF (K .LT. N) CALL dpSAXPY(LM, T, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0/DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C                           ** SOLVE  U*Z = W
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K)) / ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
         LM = MIN0(K, M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL dpSAXPY(LM, T, ABD(LA,K), 1, Z(LZ), 1)
  160 CONTINUE
C                              ** MAKE ZNORM = 1.0
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE  dpSGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )
C
C         FACTORS A REAL*8 BAND MATRIX BY ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     INPUT:  SAME AS 'SGBCO'
C
C     ON RETURN:
C
C        ABD,IPVT    SAME AS 'SGBCO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C                       FROM FORTRAN: MAX0, MIN0
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER  LDA, N, ML, MU, IPVT(*), INFO
      REAL*8     ABD(LDA,*)
C
      REAL*8     T
      INTEGER  I,DPISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C
      M = ML + MU + 1
      INFO = 0
C                        ** ZERO INITIAL FILL-IN COLUMNS
      J0 = MU + 2
      J1 = MIN0(N, M) - 1
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
      JZ = J1
      JU = 0
C
C                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1 = N - 1
      DO 120 K = 1, NM1
         KP1 = K + 1
C                                  ** ZERO NEXT FILL-IN COLUMN
         JZ = JZ + 1
         IF (JZ .LE. N) THEN
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
         ENDIF
C                                  ** FIND L = PIVOT INDEX
         LM = MIN0(ML, N-K)
         L = DPISAMAX(LM+1, ABD(M,K), 1) + M - 1
         IPVT(K) = L + K - M
C
         IF (ABD(L,K) .EQ. 0.0E0) THEN
C                                      ** ZERO PIVOT IMPLIES THIS COLUMN
C                                      ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                ** INTERCHANGE IF NECESSARY
            IF (L .NE. M) THEN
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
            ENDIF
C                                   ** COMPUTE MULTIPLIERS
            T = -1.0E0 / ABD(M,K)
            CALL dpSSCAL(LM, T, ABD(M+1,K), 1)
C
C                               ** ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
            MM = M
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .NE. MM) THEN
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
               ENDIF
               CALL dpSAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE
C
         ENDIF
C
  120 CONTINUE
C
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE  dpSGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )
C
C         SOLVES THE REAL*8 BAND SYSTEM
C            A * X = B  OR  TRANSPOSE(A) * X = B
C         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     INPUT:
C
C        ABD     REAL*8(LDA, N)
C                THE OUTPUT FROM SBGCO OR SGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.
C
C        B       REAL*8(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, DPSDOT
C                       FROM FORTRAN: MIN0
C
      implicit real*8 (a-h, o-z)
      INTEGER  LDA, N, ML, MU, IPVT(*), JOB
      REAL*8     ABD(LDA,*), B(*)
C
      REAL*8     DPSDOT,T
      INTEGER  K,KB,L,LA,LB,LM,M,NM1
C
C
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                               ** JOB = 0 , SOLVE  A * X = B
C                               ** FIRST SOLVE L*Y = B
         IF (ML .NE. 0) THEN
            DO 20 K = 1, NM1
               LM = MIN0(ML, N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .NE. K) THEN
                  B(L) = B(K)
                  B(K) = T
               ENDIF
               CALL dpSAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
         ENDIF
C                           ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / ABD(M,K)
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL dpSAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE
C
      ELSE
C                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                  ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = DPSDOT(LM, ABD(LA,K), 1, B(LB), 1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C                                  ** NOW SOLVE TRANS(L)*X = Y
         IF (ML .NE. 0) THEN
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML, N-K)
               B(K) = B(K) + DPSDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
               L = IPVT(K)
               IF (L .NE. K) THEN
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
               ENDIF
   80       CONTINUE
         ENDIF
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE  dpSGECO( A, LDA, N,IPVT, RCOND, Z )
C
C         FACTORS A REAL*8 MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C
C     ON ENTRY
C
C        A       REAL*8(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL*8
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL*8(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     ROUTINES CALLED:  FROM LINPACK: SGEFA
C                       FROM BLAS:    SAXPY, DPSDOT, SSCAL, DPSASUM
C                       FROM FORTRAN: ABS, dmax1, SIGN
C
      implicit real*8 (a-h,o-z)
      INTEGER  LDA, N, IPVT(*)
      REAL*8     A(LDA,*), Z(*)
      REAL*8     RCOND
C
      REAL*8     DPSDOT,EK,T,WK,WKM
      REAL*8     ANORM,S,DPSASUM,SM,YNORM
      INTEGER  INFO,J,K,KB,KP1,L
C
C
C                        ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = DMAX1( ANORM, DPSASUM(N,A(1,J),1) )
   10 CONTINUE
C                                      ** FACTOR
     
      CALL dpSGEFA(A,LDA,N,IPVT,INFO)
c      if (dabs(a(12,12)).lt.1.D-5)
c     $  write(*,*)a(12,12),'sgefa!'
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C                        ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
C
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K)) / ABS(EK-Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .NE. 0.0E0) THEN
            WK  = WK  / A(K,K)
            WKM = WKM / A(K,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE
C
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
C                                ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DPSDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
C                                 ** SOLVE L*V = Y
      YNORM = 1.0E0
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL dpSAXPY(N-K, T, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE
C
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
C                                  ** SOLVE  U*Z = V
      YNORM = S*YNORM
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K))/ABS(Z(K))
            CALL dpSSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL dpSAXPY(K-1, T, A(1,K), 1, Z(1), 1)
  160 CONTINUE
C                                   ** MAKE ZNORM = 1.0
      S = 1.0E0 / DPSASUM(N, Z, 1)
      CALL dpSSCAL(N, S, Z, 1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE  dpSGEFA( A, LDA, N, IPVT, INFO )
C
C         FACTORS A REAL*8 MATRIX BY GAUSSIAN ELIMINATION.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     INPUT:  SAME AS 'SGECO'
C
C     ON RETURN:
C
C        A,IPVT  SAME AS 'SGECO'
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C
      implicit real*8 (a-h,o-z)
      INTEGER  LDA, N, IPVT(*), INFO
      REAL*8     A(LDA,*)
C
      REAL*8     T
      INTEGER  DPISAMAX,J,K,KP1,L,NM1
C
C
C                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      INFO = 0
      NM1 = N - 1
      DO 60 K = 1, NM1
         KP1 = K + 1
C                                            ** FIND L = PIVOT INDEX
         L = DPISAMAX( N-K+1, A(K,K), 1) + K-1
         IPVT(K) = L
C
         IF (A(L,K) .EQ. 0.0E0) THEN
C                                     ** ZERO PIVOT IMPLIES THIS COLUMN
C                                     ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
C                                     ** INTERCHANGE IF NECESSARY
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            ENDIF
	    
C                                     ** COMPUTE MULTIPLIERS
            T = -1.0E0 / A(K,K)
            CALL dpSSCAL( N-K, T, A(K+1,K), 1 )
C
C                              ** ROW ELIMINATION WITH COLUMN INDEXING
          
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               ENDIF
	      
               CALL dpSAXPY( N-K, T, A(K+1,K), 1, A(K+1,J), 1 )
c>>>>>>>>>	       if(dabs(a(k+1,j)).lt.1.D-5)write(*,*)K,J,KP1,N,a(k+1,j)
   30       CONTINUE

C
         ENDIF
C
   60 CONTINUE
C
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE  dpSGESL( A, LDA, N,IPVT, B, JOB )
C
C         SOLVES THE REAL*8 SYSTEM
C            A * X = B  OR  TRANS(A) * X = B
C         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
C
C     ON ENTRY
C
C        A       REAL*8(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL*8(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C
C     ROUTINES CALLED:  FROM BLAS:    SAXPY, DPSDOT
C
      implicit real*8 (a-h,o-z)
      INTEGER  LDA, N, IPVT(*), JOB
      REAL*8     A(LDA,*), B(*)
C
      REAL*8     DPSDOT,T
      INTEGER  K,KB,L,NM1
C
C
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
C                                 ** JOB = 0 , SOLVE  A * X = B
C                                     ** FIRST SOLVE  L*Y = B
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .NE. K) THEN
               B(L) = B(K)
               B(K) = T
            ENDIF
            CALL dpSAXPY( N-K, T, A(K+1,K), 1, B(K+1), 1 )
   20    CONTINUE
C                                    ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / A(K,K)
            T = -B(K)
            CALL dpSAXPY( K-1, T, A(1,K), 1, B(1), 1 )
   40    CONTINUE
C
      ELSE
C                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                    ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            T = DPSDOT( K-1, A(1,K), 1, B(1), 1 )
            B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
C                                    ** NOW SOLVE  TRANS(L)*X = Y
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DPSDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
            L = IPVT(K)
            IF (L .NE. K) THEN
               T = B(L)
               B(L) = B(K)
               B(K) = T
            ENDIF
   80    CONTINUE
C
      ENDIF
C
      RETURN
      END
      REAL*8 FUNCTION  dpSASUM( N, SX, INCX )
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- DPSASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
C
      implicit real*8 (a-h,o-z)
      REAL*8 SX(*)
C
C
      DPSASUM = 0.0
      IF( N.LE.0 )  RETURN
      IF( INCX.NE.1 ) THEN
C                                          ** NON-UNIT INCREMENTS
          DO 10 I = 1, 1+(N-1)*INCX, INCX
             DPSASUM = DPSASUM + ABS(SX(I))
   10     CONTINUE
      ELSE
C                                          ** UNIT INCREMENTS
         M = MOD(N,6)
         IF( M.NE.0 ) THEN
C                             ** CLEAN-UP LOOP SO REMAINING VECTOR
C                             ** LENGTH IS A MULTIPLE OF 6.
            DO 30  I = 1, M
              DPSASUM = DPSASUM + ABS(SX(I))
   30       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 6
           DPSASUM = DPSASUM + ABS(SX(I)) + ABS(SX(I+1)) + ABS(SX(I+2))
     $                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE     dpSAXPY( N, SA, SX, INCX, SY, INCY )
C
C          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH
C                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
C            WHERE LX = 1          IF INCX .GE. 0,
C                     = (-INCX)*N  IF INCX .LT. 0
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      implicit real*8 (a-h,o-z)
      REAL*8 SX(*), SY(*), SA
C
C
cc      if(dabs(sy(12)).ne.0.AND.dabs(sy(12)).lt.1.D-10)
cc     $      write(*,*)sy(12)

      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,4)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
              SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 4
            SY(I)   = SY(I)   + SA * SX(I)
            SY(I+1) = SY(I+1) + SA * SX(I+1)
            SY(I+2) = SY(I+2) + SA * SX(I+2)
            SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SY(IY) = SY(IY) + SA*SX(IX)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C     
      ENDIF
     
      RETURN
      END
      REAL*8 FUNCTION  DPSDOT( N, SX, INCX, SY, INCY )
C
C          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C     DPSDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C            WHERE  LX = 1          IF INCX .GE. 0,
C                      = (-INCX)*N  IF INCX .LT. 0,
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      implicit real*8 (a-h,o-z)
      REAL*8 SX(*), SY(*)
C
C
      DPSDOT = 0.0
      IF( N.LE.0 )  RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             DPSDOT = DPSDOT + SX(I) * SY(I)
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,5)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
               DPSDOT = DPSDOT + SX(I) * SY(I)
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 5
            DPSDOT = DPSDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)
     $                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     $                  + SX(I+4)*SY(I+4)
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            DPSDOT = DPSDOT + SX(IX) * SY(IY)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE     dpSSCAL( N, SA, SX, INCX )
C
C         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
C            SA  SINGLE PRECISION SCALE FACTOR
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX)
C                FOR I = 0 TO N-1
C
      implicit real*8 (a-h,o-z)
      REAL*8 SA, SX(*)
C
C
      IF( N.LE.0 ) RETURN
C
      IF( INCX.NE.1 ) THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SX(I) = SA * SX(I)
   10     CONTINUE
C
      ELSE
C
         M = MOD(N,5)
         IF( M.NE.0 ) THEN
C                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                           ** IS A MULTIPLE OF 5.
            DO 30  I = 1, M
               SX(I) = SA * SX(I)
   30       CONTINUE
         ENDIF
C                             ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 5
            SX(I)   = SA * SX(I)
            SX(I+1) = SA * SX(I+1)
            SX(I+2) = SA * SX(I+2)
            SX(I+3) = SA * SX(I+3)
            SX(I+4) = SA * SX(I+4)
   50    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE     dpSSWAP( N, SX, INCX, SY, INCY )
C
C          INTERCHANGE S.P VECTORS  X  AND  Y
C
C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
C
C --OUTPUT--
C       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
C       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)
C
C     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
C     WHERE LX = 1          IF INCX .GE. 0,
C              = (-INCX)*N  IF INCX .LT. 0
C     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      implicit real*8 (a-h,o-z)
      REAL*8 SX(*), SY(*), STEMP1, STEMP2, STEMP3
C
C
      IF( N.LE.0 ) RETURN
C
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
C
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             STEMP1 = SX(I)
             SX(I) = SY(I)
             SY(I) = STEMP1
   10     CONTINUE
C
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
C
C                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,3)
         IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 3.
            DO 20  I = 1, M
               STEMP1 = SX(I)
               SX(I) = SY(I)
               SY(I) = STEMP1
   20       CONTINUE
         ENDIF
C                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 3
            STEMP1  = SX(I)
            STEMP2  = SX(I+1)
            STEMP3  = SX(I+2)
            SX(I)   = SY(I)
            SX(I+1) = SY(I+1)
            SX(I+2) = SY(I+2)
            SY(I)   = STEMP1
            SY(I+1) = STEMP2
            SY(I+2) = STEMP3
   30    CONTINUE
C
      ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            STEMP1 = SX(IX)
            SX(IX) = SY(IY)
            SY(IY) = STEMP1
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
C
      ENDIF
C
      RETURN
      END
      INTEGER FUNCTION  dpISAMAX( N, SX, INCX )
C
C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
C
C --OUTPUT-- DPISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
C                         ABS(SX(1+(I-1)*INCX))
C
      implicit real*8 (a-h,o-z)
      REAL*8 SX(*), SMAX, XMAG
C
C
      IF( N.LE.0 ) THEN
         DPISAMAX = 0
      ELSE IF( N.EQ.1 ) THEN
         DPISAMAX = 1
      ELSE
         SMAX = 0.0
         II = 1
         DO 20  I = 1, 1+(N-1)*INCX, INCX
            XMAG = ABS(SX(I))
            IF( SMAX.LT.XMAG ) THEN
               SMAX = XMAG
               DPISAMAX = II
            ENDIF
            II = II + 1
   20    CONTINUE
      ENDIF
C
      RETURN
      END


      SUBROUTINE spline2(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(0:n-1),y(0:n-1),y2(0:n-1)
      PARAMETER (NMAX=1400)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(0:NMAX-1)
      
      if (yp1.gt.0.99d30) then
        y2(0)=0.d0
        u(0)=0.d0
      else
        y2(0)=-0.5d0
        u(0)=(3.d0/(x(1)-x(0)))*((y(1)-y(0))/(x(1)-x(0))-yp1)
      endif
      
      do 11 i=1,n-2
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p      
11    continue

      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/
     $	   (x(n-1)-x(n-2)))*(ypn-(y(n-1)-y(n-2))/(x(n-1)-x(n-2)))
      endif
      y2(n-1)=(un-qn*u(n-2))/(qn*y2(n-2)+1.D0)
      do 12 k=n-2,0,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software )?".

      SUBROUTINE splint2(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(0:n-1),y2a(0:n-1),ya(0:n-1)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=0
      khi=n-1
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) write (*,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)
     *+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software )?".
