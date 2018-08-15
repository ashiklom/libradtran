C
C PAGE 277-278: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: SPL3.FOR
C
C NATURAL CUBIC SPLINE FUNCTION FOR SIN(X) AT EQUIDISTANT POINTS (SPL3,ZSPL3)
C
  
      SUBROUTINE ZSPL3(N,T,Y,H,B,U,V,Z) 
      DIMENSION  T(N),Y(N),H(N),B(N),U(N),V(N),Z(N)       
      DO 2 I = 1,N-1
        H(I) = T(I+1) - T(I)
        B(I) = (Y(I+1) -Y(I))/H(I)    
   2  CONTINUE    
      U(2) = 2.0*(H(1) + H(2))
      V(2) = 6.0*(B(2) - B(1))
      DO 3 I = 3,N-1
        U(I) = 2.0*(H(I) + H(I-1)) - H(I-1)**2/U(I-1)     
        V(I) = 6.0*(B(I) - B(I-1)) - H(I-1)*V(I-1)/U(I-1) 
   3  CONTINUE    
      Z(N) = 0.0  
      DO 4 I = N-1,2,-1     
        Z(I) = (V(I) - H(I)*Z(I+1))/U(I)
   4  CONTINUE    
      Z(1) = 0.0
      RETURN
      END 
  
      FUNCTION SPL3(N,T,Y,Z,X)
      DIMENSION  T(N),Y(N),Z(N)       
      DO 2 I = N-1,1,-1     
        DIFF = X - T(I)     
        IF(DIFF .GE. 0.0)  GO TO 3    
   2  CONTINUE    
      I = 1 
   3  H = T(I+1) - T(I)     
      B = (Y(I+1) - Y(I))/H - H*(Z(I+1) + 2.0*Z(I))/6.0 
      P = 0.5*Z(I) + DIFF*(Z(I+1) - Z(I))/(6.0*H) 
      P = B + DIFF*P
      SPL3 = Y(I) + DIFF*P  
      RETURN      
      END 

      SUBROUTINE TRI(N,A,D,C,B,X)     
      DIMENSION  A(N),D(N),C(N),B(N),X(N) 
      DO 2 I = 2,N
        XMULT = A(I-1)/D(I-1) 
        D(I) = D(I) - XMULT*C(I-1)    
        B(I) = B(I) - XMULT*B(I-1)    
   2  CONTINUE    
      X(N) = B(N)/D(N)      
      DO 3 I = N-1,1,-1     
        X(I) = (B(I) - C(I)*X(I+1))/D(I)
   3  CONTINUE    
      RETURN      
      END 
  
