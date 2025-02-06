     FUNCTION znplus2(n11,n22)
         use input
         use hmtconstants
         use constants
         
         IMPLICIT NONE
       
       
         COMPLEX*16 :: znplus2
         REAL*8     :: factorial, dfactorial
         REAL*8     :: n11, n22
         
         call centerinput
       
         znplus2 = (0.00d+00,0.00d+00)
         
         mx = (n11 + n22)
         
         do i = 0,int(mx),2
           
           id = dfloat(i)
         
           ii = max(0,int(id-n22))
           jj = min(int(id),int(n11))
         
           f = (0.00d+00, 0.00d+00)
         
           do j = ii,jj
             jd = dfloat(j)
         
             f = f + ( (factorial(n11)*factorial(n22)*((PAz)**(n11-jd))*((PBz)**(n22-id+jd)))/&
                       (factorial(jd)*factorial(n11-jd)*factorial(id-jd)*factorial(n22-id+jd)) )
           enddo
         
           znplus2 = znplus2 + f * ( (dfactorial(id + 1.00d+00)*sqrt(pi))/&
                    ( (2.00d+00)**((id + 2.00d+00)/2.00d+00) * gmma**((id + 3.00d+00)/2.00d+00) ) )
         enddo
          
     END FUNCTION

! -------------------------------------------------- !

     FUNCTION znplus1(n11,n22)
         use input
         use hmtconstants
         use constants
         
         IMPLICIT NONE
       
         COMPLEX*16 :: znplus1
         REAL*8     :: factorial, dfactorial
         REAL*8     :: n11, n22
       
         call centerinput
       
         znplus1 = (0.00d+00, 0.00d+00)
         
         mx = (n11 + n22)
       
         if (mx .eq. 0.00d+00) then
           znplus1 = (0.00d+00,0.00d+00) 
             
          
         else 
           do i = 1,int(mx),2
             
             id = dfloat(i)
             
             ii = max(0,int(id-n22))
             jj = min(int(id),int(n11))
             
             f = 0.00d+00
             
             do j = ii,jj
               jd = dfloat(j)
              
               f = f + ( (factorial(n11)*factorial(n22)*((PAz)**(n11-jd))*((PBz)**(n22-id+jd)))/&
                         (factorial(jd)*factorial(n11-jd)*factorial(id-jd)*factorial(n22-id+jd)) )
             enddo
              
             znplus1 = znplus1 + f * ( (dfactorial(id)*sqrt(pi))/&
                      ( (2.00d+00)**((id + 1.00d+00)/2.00d+00) * gmma**((id + 2.00d+00)/2.00d+00) ) )
           enddo
          
         endif
       
     END FUNCTION 

!--------------------------------------

     FUNCTION xterm(l11, l22, kk)
         use input
         use constants
         use hmtconstants
       
         IMPLICIT NONE
       
         REAL*8     :: l11, l22, kk
         COMPLEX*16 :: xterm
         REAL*8     :: factorial, dfactorial
       
!-------------------------------------

         call centerinput
       
         xterm = (0.0000d+00, 0.0000d+00)
       
         mx = l11 + l22
       
         do i = 0,int(mx)
           id = dfloat(i)
       
           ii = max(0,i-int(l22))
           jj = min(i,int(l11))
       
           f = (0.000d+00, 0.000d+00)
       
           do j = ii,jj
             jd = dfloat(j)
       
             f = f + ( (factorial(l11)*factorial(l22)*((PAx)**(l11-jd))*((PBx)**(l22-id+jd)))/&
                       (factorial(jd)*factorial(l11-jd)*factorial(id-jd)*factorial(l22-id+jd)) )
           enddo
       
       
           intt = int(i/2)
       
           eikx = (0.00d+00,0.00d+00)
           do t = 0,intt
             td = dfloat(t)
       
             eikx  = eikx + ((-1.00d+00)**td)*( (-kk/sqrt(gmma))**(id-2.00d+00*td) )/(factorial(td)*factorial(id-2.00d+00*td))
           enddo
       
           xterm = xterm + f*(1.00d+00/(4.00d+00*gmma)**(id/2.00d+00))*(-1.00d+00)**(id)*(io)**(id)*factorial(id)*eikx
         enddo
       
       
         xterm = (pi/gmma)**(1.00d+00/2.00d+00)*xterm*cdexp(-(kk**2.00d+00)/(4.00d+00*gmma) + io*kk*px)
       
     END FUNCTION

!---------------------------------------

     FUNCTION ovlpart(lone,ltwo,ppa,ppb)
     
         use input
         use constants
         use hmtconstants
       
         IMPLICIT NONE
       
         REAL*8     :: lone, ltwo
         COMPLEX*16 :: ovlpart
         REAL*8     :: factorial, dfactorial
         REAL*8     :: ppa, ppb
       
         call centerinput
       
         ovlpart = (0.00d+00,0.00d+00)
       
         mx = (lone + ltwo)/2.00d+00
       
         do i = 0,int(mx)
           id = dfloat(i)
       
           ii = max(0,int(2.00d+00*id-ltwo))
           jj = min(int(2.00d+00*id),int(lone))
       
           f = (0.00d+00,0.00d+00)
       
           do j = ii,jj
             jd = dfloat(j)
       
             f = f + ( (factorial(lone)*factorial(ltwo)*((ppa)**(lone-jd))*((ppb)**(ltwo-2.00d+00*id+jd)))/&
                (factorial(jd)*factorial(lone-jd)*factorial(2.00d+00*id-jd)*factorial(ltwo-2.00d+00*id+jd)) )
           enddo
       
       
           ovlpart = ovlpart + f * (dfactorial(2.00d+00*id - 1.00d+00))/(2.00d+00*gmma)**id
         enddo
       
         ovlpart = dsqrt(pi/gmma)*ovlpart                                                                                                                
     END FUNCTION
     
     
     !-----------------------
     
     SUBROUTINE centerinput
     
         use input
         use constants
         use hmtconstants
       
         IMPLICIT NONE
       
       
         gmma = alpha1 + alpha2
       
         px = (alpha1*ax + alpha2*bx)/(alpha1 + alpha2)
         py = (alpha1*ay + alpha2*by)/(alpha1 + alpha2)
         pz = (alpha1*az + alpha2*bz)/(alpha1 + alpha2)
       
         PAx = px - ax
         PBx = px - bx
       
         PAy = py - ay
         PBy = py - by
       
         PAz = pz - az
         PBz = pz - bz
       
         ABx = ax - bx
         ABy = ay - by
         ABz = az - bz
       
         PCx = px - cx
         PCy = py - cy
         PCz = pz - cz
       
     END SUBROUTINE
     
     !--------------------
     
     FUNCTION expp(ABx, ABy, ABz, alp1, alp2)
     
         use input
         implicit none
         
         real*8 :: ABx, ABy, ABz, ABsqr
         real*8 :: expp
         real*8 :: alp1, alp2, gmma
       
         gmma = alp1 + alp2
       
         ABsqr= (ABx**2.00d+00) + (ABy**2.00d+00) + (ABz**2.00d+00)
       
         expp = dexp((-alp1*alp2*ABsqr)/gmma)
       
     END FUNCTION
     
     
     
     FUNCTION norm(l11,m11,n11,l22,m22,n22,alp1,alp2)
     
         implicit none
       
         real*8 :: norm 
         real*8, parameter :: pi = 3.14159265358979d+00
         real*8 :: l11, m11, n11, l22, m22, n22
         real*8 :: dfactorial
         real*8 :: alp1, alp2
          
         norm = ( (2.00d+00*alp1/pi)**(3.00d+00/4.00d+00)*( ((4.00d+00*alp1)**(l11 + m11 + n11))/&
                  (dfactorial(2.00d+00*l11-1.00d+00)*dfactorial(2.00d+00*m11-1.00d+00)*&
                   dfactorial(2.00d+00*n11-1.00d+00)) )**(1.00d+00/2.00d+00) )* &
                ( (2.00d+00*alp2/pi)**(3.00d+00/4.00d+00)*( ((4.00d+00*alp2)**(l22 + m22 + n22))/&
                  (dfactorial(2.00d+00*l22-1.00d+00)*dfactorial(2.00d+00*m22-1.00d+00)*&
                   dfactorial(2.00d+00*n22-1.00d+00)) )**(1.00d+00/2.00d+00) )
         
     END FUNCTION
     
     
     
     FUNCTION factorial(n)
         IMPLICIT NONE
       
         real*8 :: factorial, n
       
         if (n == 0.00d+00) then
           factorial = 1.00d+00
         elseif (n == -1.00d+00) then
           factorial = 1.00d+00 
         elseif (n == 1.00d+00) then
           factorial = 1.00d+00
         elseif (n == 2.00d+00) then
           factorial = 2.00d+00
         elseif (n == 3.00d+00) then
           factorial = 6.00d+00
         elseif (n == 4.00d+00) then
           factorial = 24.00d+00
         elseif (n == 5.00d+00) then
           factorial = 120.00d+00
         elseif (n == 6.00d+00) then
           factorial = 720.00d+00
         elseif (n == 7.00d+00) then
           factorial = 5040.00d+00
         elseif (n == 8.00d+00) then
           factorial = 40320.00d+00
         elseif (n == 9.00d+00) then
           factorial = 362880.00d+00
         elseif (n == 10.00d+00) then
           factorial = 3628800.00d+00
         elseif (n == 11.00d+00) then
           factorial = 39916800.00d+00
         elseif (n == 12.00d+00) then
           factorial = 479001600.00d+00
         endif
     END FUNCTION
     
     
     
     FUNCTION dfactorial(n)
         IMPLICIT NONE
       
         real*8 :: dfactorial, n
       
       
         if (n == 0.00d+00) then
           dfactorial = 1.00d+00
         elseif (n == -1.00d+00) then
           dfactorial = 1.00d+00
         elseif (n == 1.00d+00) then
           dfactorial = 1.00d+00
         elseif (n == 2.00d+00) then
           dfactorial = 2.00d+00
         elseif (n == 3.00d+00) then
           dfactorial = 3.00d+00
         elseif (n == 4.00d+00) then
           dfactorial = 8.00d+00
         elseif (n == 5.00d+00) then
           dfactorial = 15.00d+00
         elseif (n == 6.00d+00) then
           dfactorial = 48.00d+00
         elseif (n == 7.00d+00) then
           dfactorial = 105.00d+00
         elseif (n == 8.00d+00) then
           dfactorial = 384.00d+00
         elseif (n == 9.00d+00) then
           dfactorial = 945.00d+00
         elseif (n == 10.00d+00) then
           dfactorial = 3840.00d+00
         elseif (n == 11.00d+00) then
           dfactorial = 10395.00d+00
         elseif (n == 12.00d+00) then
           dfactorial = 46080.00d+00
         endif
     END FUNCTION
     

