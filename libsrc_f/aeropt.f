      SUBROUTINE aeropt( vulcan, aerext, aerabs, aersym, aerwvn, rh,
     $     nlyr, zd, haze, maxawvn )
*
* Gets the appropriate asymmetry, absorption and extinction
* coefficients for the different aerosol models and cirrus
* clouds. In part adopted from MODTRAN, hence see MODTRAN
* documentation for details.
*
* Input/output variables described in phodis.f
*   
      IMPLICIT LOGICAL (A-Z)    ! To nearly get strong typing
      INTEGER haze, maxawvn, nlyr, vulcan
      REAL aerabs(maxawvn,0:*), aerext(maxawvn,0:*), 
     $     aersym(maxawvn,0:*), aerwvn(*), rh(0:*), zd(0:*) 
*
      INTEGER irh, iv, jrh1, jrh2, lc
      REAL a, a1, a2, x, x1, x2, y, y1, y2, z1, z2, zk 
*     
      REAL rhzone(4)
*
      INCLUDE "extinc.f"
*
      DATA rhzone/0.,70.,80.,99./
*     
      DO lc = 0, nlyr       
*     
*     Get right relative humidity interval for 0-10km
*
            IF ( zd(lc) .GE. 0.0 .AND. zd(lc) .LE. 10.0 ) THEN
               DO irh = 2, 4
                  IF ( rh(lc) .LT. rhzone(irh) ) GOTO 10 
               ENDDO
 10            CONTINUE
               jrh1 = irh
               IF ( jrh1 .GT. 4 ) jrh1 = 4
               jrh2 = jrh1 - 1
               IF ( rh(lc) .GT. 0.0 .AND. rh(lc) .LT. 99. ) 
     $              x = alog(100.0-rh(lc))
               x1 = alog(100.0-rhzone(jrh2))
               x2 = alog(100.0-rhzone(jrh1))
               IF ( rh(lc) .GE. 99.0 ) x = x2
               IF ( rh(lc) .LE.  0.0 ) x = x1
            ENDIF
*     
*     0-2km
*     
            IF ( zd(lc) .GE. 0.0 .AND. zd(lc) .LE. 2.0 ) THEN
               DO iv = 1, maxawvn
                  IF ( haze .EQ. 1 ) THEN 
                     IF ( rurext(iv,jrh1).GT.0.0 .AND.
     $                    rurext(iv,jrh2).GT.0.0 ) THEN
                        y2=alog(rurext(iv,jrh1))
                        y1=alog(rurext(iv,jrh2))
                        y=y1+(y2-y1)*(x-x1)/(x2-x1)
                        aerext(iv,lc)=exp(y)
                     ELSE
                        aerext(iv,lc)=0.0
                     ENDIF
                     IF ( rurabs(iv,jrh1).GT.0.0 .AND.
     $                    rurabs(iv,jrh2).GT.0.0 ) THEN
                        z2=alog(rurabs(iv,jrh1))
                        z1=alog(rurabs(iv,jrh2))
                        zk=z1+(z2-z1)*(x-x1)/(x2-x1)
                        aerabs(iv,lc)=exp(zk)
                     ELSE
                        aerabs(iv,lc)=0.0
                     ENDIF
                     IF ( rursym(iv,jrh1).GT.0.0 .AND.
     $                    rursym(iv,jrh2).GT.0.0 ) THEN
                        a2=alog(rursym(iv,jrh1))
                        a1=alog(rursym(iv,jrh2))
                        a=a1+(a2-a1)*(x-x1)/(x2-x1)
                        aersym(iv,lc)=exp(a)
                     ELSE
                        aersym(iv,lc)=0.0
                     ENDIF
                  ELSE IF ( haze .EQ. 4 ) THEN 
                     IF ( ocnext(iv,jrh1).GT.0.0 .AND.
     $                    ocnext(iv,jrh2).GT.0.0 ) THEN
                        y2=alog(ocnext(iv,jrh1))
                        y1=alog(ocnext(iv,jrh2))
                        y=y1+(y2-y1)*(x-x1)/(x2-x1)
                        aerext(iv,lc)=exp(y)
                     ELSE
                        aerext(iv,lc)=0.0
                     ENDIF
                     IF ( ocnabs(iv,jrh1).GT.0.0 .AND.
     $                    ocnabs(iv,jrh2).GT.0.0 ) THEN
                        z2=alog(ocnabs(iv,jrh1))
                        z1=alog(ocnabs(iv,jrh2))
                        zk=z1+(z2-z1)*(x-x1)/(x2-x1)
                        aerabs(iv,lc)=exp(zk)
                     ELSE
                        aerabs(iv,lc)=0.0
                     ENDIF
                     IF ( ocnsym(iv,jrh1).GT.0.0 .AND.
     $                    ocnsym(iv,jrh2).GT.0.0 ) THEN
                        a2=alog(ocnsym(iv,jrh1))
                        a1=alog(ocnsym(iv,jrh2))
                        a=a1+(a2-a1)*(x-x1)/(x2-x1)
                        aersym(iv,lc)=exp(a)
                     ELSE
                        aersym(iv,lc)=0.0
                     ENDIF
                  ELSE IF ( haze .EQ. 5 ) THEN 
                     IF ( urbext(iv,jrh1).GT.0.0 .AND.
     $                    urbext(iv,jrh2).GT.0.0 ) THEN
                        y2=alog(urbext(iv,jrh1))
                        y1=alog(urbext(iv,jrh2))
                        y=y1+(y2-y1)*(x-x1)/(x2-x1)
                        aerext(iv,lc)=exp(y)
                     ELSE
                        aerext(iv,lc)=0.0
                     ENDIF
                     IF ( urbabs(iv,jrh1).GT.0.0 .AND.
     $                    urbabs(iv,jrh2).GT.0.0 ) THEN
                        z2=alog(urbabs(iv,jrh1))
                        z1=alog(urbabs(iv,jrh2))
                        zk=z1+(z2-z1)*(x-x1)/(x2-x1)
                        aerabs(iv,lc)=exp(zk)
                     ELSE
                        aerabs(iv,lc)=0.0
                     ENDIF
                     IF ( urbsym(iv,jrh1).GT.0.0 .AND.
     $                    urbsym(iv,jrh2).GT.0.0 ) THEN
                        a2=alog(urbsym(iv,jrh1))
                        a1=alog(urbsym(iv,jrh2))
                        a=a1+(a2-a1)*(x-x1)/(x2-x1)
                        aersym(iv,lc)=exp(a)
                     ELSE
                        aersym(iv,lc)=0.0
                     ENDIF
                  ELSE IF ( haze .EQ. 6 ) THEN 
                     IF ( troext(iv,jrh1).GT.0.0 .AND.
     $                    troext(iv,jrh2).GT.0.0 ) THEN
                        y2=alog(troext(iv,jrh1))
                        y1=alog(troext(iv,jrh2))
                        y=y1+(y2-y1)*(x-x1)/(x2-x1)
                        aerext(iv,lc)=exp(y)
                     ELSE
                        aerext(iv,lc)=0.0
                     ENDIF
                     IF ( troabs(iv,jrh1).GT.0.0 .AND.
     $                    troabs(iv,jrh2).GT.0.0 ) THEN
                        z2=alog(troabs(iv,jrh1))
                        z1=alog(troabs(iv,jrh2))
                        zk=z1+(z2-z1)*(x-x1)/(x2-x1)
                        aerabs(iv,lc)=exp(zk)
                     ELSE
                        aerabs(iv,lc)=0.0
                     ENDIF
                     IF ( trosym(iv,jrh1).GT.0.0 .AND.
     $                    trosym(iv,jrh2).GT.0.0 ) THEN
                        a2=alog(trosym(iv,jrh1))
                        a1=alog(trosym(iv,jrh2))
                        a=a1+(a2-a1)*(x-x1)/(x2-x1)
                        aersym(iv,lc)=exp(a)
                     ELSE
                        aersym(iv,lc)=0.0
                     ENDIF
                  ENDIF
               ENDDO
*

*     >2-10km
*     
            ELSE IF ( zd(lc) .GT. 2.0 .AND. zd(lc) .LE. 10.0 ) THEN
               DO iv = 1, maxawvn
                  IF ( troext(iv,jrh1).GT.0.0 .AND.
     $                 troext(iv,jrh2).GT.0.0 ) THEN
                     y2=alog(troext(iv,jrh1))
                     y1=alog(troext(iv,jrh2))
                     y=y1+(y2-y1)*(x-x1)/(x2-x1)
                     aerext(iv,lc)=exp(y)
                  ELSE
                     aerext(iv,lc)=0.0
                  ENDIF
                  IF ( troabs(iv,jrh1).GT.0.0 .AND.
     $                 troabs(iv,jrh2).GT.0.0 ) THEN
                     z2=alog(troabs(iv,jrh1))
                     z1=alog(troabs(iv,jrh2))
                     zk=z1+(z2-z1)*(x-x1)/(x2-x1)
                     aerabs(iv,lc)=exp(zk)
                  ELSE
                     aerabs(iv,lc)=0.0
                  ENDIF
                  IF ( trosym(iv,jrh1).GT.0.0 .AND.
     $                 trosym(iv,jrh2).GT.0.0 ) THEN
                     a2=alog(trosym(iv,jrh1))
                     a1=alog(trosym(iv,jrh2))
                     a=a1+(a2-a1)*(x-x1)/(x2-x1)
                     aersym(iv,lc)=exp(a)
                  ELSE
                     aersym(iv,lc)=0.0
                  ENDIF
               ENDDO
*     
*     >10-30km
*     
            ELSE IF ( zd(lc) .GT. 10.0 .AND. zd(lc) .LE. 30.0 ) THEN
               IF ( vulcan .EQ. 1 .OR. vulcan .EQ. 2 ) THEN
                  DO iv = 1, maxawvn
                     aerabs(iv,lc) = bstabs(iv)
                     aerext(iv,lc) = bstext(iv)
                     aersym(iv,lc) = bstsym(iv)
                  ENDDO
               ELSE IF (vulcan .EQ. 3 ) THEN
                  DO iv = 1, maxawvn
                     aerabs(iv,lc) = avoabs(iv)
                     aerext(iv,lc) = avoext(iv)
                     aersym(iv,lc) = avosym(iv)
                  ENDDO
               ELSE IF (vulcan .EQ. 4 ) THEN
                  DO iv = 1, maxawvn
                     aerabs(iv,lc) = fvoabs(iv)
                     aerext(iv,lc) = fvoext(iv)
                     aersym(iv,lc) = fvosym(iv)
                  ENDDO
               ENDIF
*     
*     >30-100km
*     
            ELSE IF ( zd(lc) .GT. 30.0 .AND. zd(lc) .LE. 100.0 ) THEN
               DO iv = 1, maxawvn
                  aerabs(iv,lc) = dmeabs(iv)
                  aerext(iv,lc) = dmeext(iv)
                  aersym(iv,lc) = dmesym(iv)
               ENDDO
            ENDIF
*     
         ENDDO
*     
         DO iv = 1, maxawvn       !Convert from micro-meter to nanometers
            aerwvn(iv) = vx2(iv)*1.E+03 
         ENDDO
*     
      RETURN
      END
