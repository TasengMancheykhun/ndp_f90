!-------------------------------------------------------------------!
      subroutine dpint 
!-------------------------------------------------------------------!

      use grid
      use atom  
      use basis  
      use mos 
      use gpt  
      use constants

      implicit none

      integer(kind=4)          :: p,q,i,j,ij,ik,jk,il,jl,kl 
      integer(kind=4)          :: na,la,ma,nb,lb,mb,in1,en1,in2,en2
      integer(kind=4)          :: nc,lc,mc
      real(kind=8)             :: totm,ccx,ccy,ccz,px,py,pz
      real(kind=8)             :: gmma,a,b,ax,ay,az,bx,by,bz 
      real(kind=8)             :: ab2,pax,pay,paz,pbx,pby,pbz
      real(kind=8)             :: val1,val2,mmx,mmy,mmz,pcx,pcy,pcz
      real(kind=8)             :: sumx,sumy,sumz,norma,normb 
      real(kind=8)             :: sumall,num1,c1,c2,sx,sy,sz
      real(kind=8)             :: sxm,sym,szm,mux,muy,muz
      real(kind=8)             :: mux1,mux2,muy1,muy2,muz1,muz2
      real(kind=8)             :: mux3,muy3,muz3,mux4,muy4,muz4
      real(kind=8)             :: summx,summy,summz,sumdp
      real(kind=8),allocatable :: npgxm(:,:),npgym(:,:),npgzm(:,:)
      real(kind=8),allocatable :: nbasxm(:,:),nbasym(:,:),nbaszm(:,:)
      real(kind=8)             :: con=2.541766 !(1.0d0/0.3934560d0)

      ALLOCATE(npgxm(npg, npg), npgym(npg, npg), npgzm(npg, npg))
      ALLOCATE(nbasxm(nbas, nbas), nbasym(nbas, nbas), nbaszm(nbas, nbas))

      npgxm = 0.0d0; npgym = 0.0d0; npgzm = 0.0d0
      nbasxm = 0.0d0; nbasym = 0.0d0; nbaszm = 0.0d0  

      ! --- center of mass ----------- !      
      totm = 0.0d0; ccx = 0.0d0; ccy = 0.0d0; ccz = 0.0d0

      DO i = 1, natom  
        totm = totm + ch(i)
        ccx = ccx + ch(i)*xa(i)
        ccy = ccy + ch(i)*ya(i)
        ccz = ccz + ch(i)*za(i)
        WRITE(*,*) xa(i), ya(i), za(i), ch(i)
      ENDDO 

      ccx = ccx/totm  
      ccy = ccy/totm 
      ccz = ccz/totm 

      WRITE(*,*)
      WRITE(*,*) ccx, ccy, ccz

      !--- specify the specific point for integral to be calculated 

!      ccx = 0.000000; ccy = 0.000000; ccz = -0.123509

      !--- when calcuating the dipole moment

      ! ccx = 0.000000; ccy = 0.000000; ccz = 0.000000
      nc = 0; lc = 0; mc = 1

      do 1201 p = 1, npg       
        a=npgalp(p); norma=npgnorm(p)          
        na=pn(p); la=pl(p); ma=pm(p)
        ax=pxa(p); ay=pya(p); az=pza(p)
        do 1202 q = 1,npg
          b=npgalp(q); normb=npgnorm(q) 
          nb=pn(q); lb=pl(q); mb=pm(q)
          bx=pxa(q); by=pya(q); bz=pza(q)

          gmma=a+b
          px=(a*ax + b*bx)/(gmma)
          py=(a*ay + b*by)/(gmma)
          pz=(a*az + b*bz)/(gmma)
          ab2=((ax-bx)**(two))+((ay-by)**(two))+((az-bz)**(two))
          val1=(pi/gmma)**(thby2)
          val2=dexp(-(a*b*ab2)/(gmma))
          pax=px-ax; pay=py-ay; paz=pz-az
          pbx=px-bx; pby=py-by; pbz=pz-bz 
          pcx=px-ccx; pcy=py-ccy; pcz=pz-ccz
      
          sx=sumx(na,nb,a,b,pax,pbx)
          sy=sumy(la,lb,a,b,pay,pby)
          sz=sumz(ma,mb,a,b,paz,pbz)
          sxm=summx(na,nb,nc,a,b,pax,pbx,pcx)
          sym=summy(la,lb,lc,a,b,pay,pby,pcy) 
          szm=summz(ma,mb,mc,a,b,paz,pbz,pcz)
     
          mmx=norma*normb*val1*val2*sxm*sym*szm
          !mmy=norma*normb*val1*val2*sx*sym*sz
          !mmz=norma*normb*val1*val2*sx*sy*szm

          npgxm(p,q)=mmx
          !npgym(p,q)=mmy
          !npgzm(p,q)=mmz

1202    enddo
1201  enddo

      !--- contraction level normalization

      !--- X INTEGRAL 
        sumall=0.0d0; num1=0.0d0; nbasxm=0.0d0
        do 1203 p = 1,nbas
          do 1204 q = 1,nbas          
            in1=inde1(p); en1=inde2(p)
            in2=inde1(q); en2=inde2(q)
            sumall=0.0d0; num1=0.0d0
            do 1205 i = in1,en1
              do 1206 j = in2,en2
               c1=npgcff(i); c2=npgcff(j)
               num1=(c1*c2*npgxm(i,j))
               sumall=sumall+num1
1206          enddo
1205        enddo
              nbasxm(p,q)=sumall
1204      enddo
1203    enddo

        !--- Y INTEGRAL      
!        sumall=0.0d0; num1=0.0d0; nbasym=0.0d0
!        do 1207 p = 1,nbas
!          do 1208 q = 1,nbas
!            in1=inde1(p); en1=inde2(p)
!            in2=inde1(q); en2=inde2(q)
!            sumall=0.0d0; num1=0.0d0
!            do 1209 i = in1,en1
!              do 1210 j = in2,en2
!               c1=npgcff(i); c2=npgcff(j)
!               num1=(c1*c2*npgym(i,j))
!               sumall=sumall+num1
!1210          enddo
!1209        enddo
!              nbasym(p,q)=sumall
!1208      enddo
!1207    enddo

        !--- Z INTEGRAL      
!        sumall=0.0d0; num1=0.0d0; nbaszm=0.0d0
!        do 1211 p = 1,nbas
!          do 1212 q = 1,nbas
!            in1=inde1(p); en1=inde2(p)
!            in2=inde1(q); en2=inde2(q)
!            sumall=0.0d0; num1=0.0d0
!            do 1213 i = in1,en1
!              do 1214 j = in2,en2
!               c1=npgcff(i); c2=npgcff(j)
!               num1=(c1*c2*npgzm(i,j))
!               sumall=sumall+num1
!1214          enddo
!1213        enddo
!              nbaszm(p,q)=sumall
!1212      enddo
!1211    enddo

      !--- end of moment integral calculations

      !--- Dipole Moment calculation
        
!      mux=0.0d0; mux1=0.0d0; mux2=0.0d0; mux3=0.0d0; mux4=0.0d0
!      muy=0.0d0; muy1=0.0d0; muy2=0.0d0; muy3=0.0d0; muy4=0.0d0
!      muz=0.0d0; muz1=0.0d0; muz2=0.0d0; muz3=0.0d0; muz4=0.0d0

     !--- term1
!      do 1216 i = 1,nbas  
!        do 1217 j = 1,nbas
!          mux1=mux1+((pmatrix(i,j))*(nbasxm(i,j)))
!          muy1=muy1+((pmatrix(i,j))*(nbasym(i,j)))
!          muz1=muz1+((pmatrix(i,j))*(nbaszm(i,j)))
!1217    enddo
!1216  enddo      

     !--- term2 
!      do 1218 i = 1,natom
!        mux2=mux2+((ch(i))*(xa(i)))
!        muy2=muy2+((ch(i))*(ya(i)))
!        muz2=muz2+((ch(i))*(za(i)))
!1218  enddo

!        mux=(-mux1+mux2)*con
!        muy=(-muy1+muy2)*con
!        muz=(-muz1+muz2)*con

!      write(*,'(a14,3F14.7)') "dipole term1 ",mux,muy,muz

      !--- printing the moments

        OPEN(12, file = 'moments.txt')
        DO i = 1, nbas
          DO j = 1, i
            write(12,*) cmplx(i,j), nbasxm(i,j)
          ENDDO
        ENDDO
        CLOSE(12)
       
!-------------------------------------------------------------------!      
      end subroutine
!-------------------------------------------------------------------!


!-------------------------------------------------------------------!
        FUNCTION sumx(l1,l2,alp1,alp2,pax,pbx)
!-------------------------------------------------------------------!
        
          use constants

          implicit none
          integer(kind=4)         :: ix,l1,l2
          real(kind=8)            :: sumx,alp1,alp2
          real(kind=8)            :: pax,pbx,gma,dfac,coeff

          gma=alp1+alp2
          sumx=0.0d0

          do 1203 ix=0,((l1+l2)/2)
            sumx=sumx+coeff((2*ix),l1,l2,pax,pbx)* &
                 & (dfac(2*(ix)-1)/((two*gma)**(ix)))   
1203      enddo
        
!-------------------------------------------------------------------!
        END
!-------------------------------------------------------------------!


!-------------------------------------------------------------------!
        function summx(l1,l2,l3,alp1,alp2,pax,pbx,pcx)
!-------------------------------------------------------------------!
        
        use constants

        implicit none
        integer(kind=4)         :: ix,ix1,ix2,l1,l2,l3,ii,ii1
        integer(kind=4)         :: rr,nn,fact,ncr
        real(kind=8)            :: summx,alp1,alp2
        real(kind=8)            :: pax,pbx,pcx,gma,dfac,coeff
        real(kind=8)            :: sumx1,sumx2,prd

        gma=alp1+alp2
        sumx1=0.0d0; sumx2=0.0d0; summx=0.0d0

        if (l3 .ge. 1) then 
          nn=l3  

          do 1203 ix1=1,(l3+1)  
            rr=ix1-1 
            ii1=l3-rr
            ncr=((fact(nn))/((fact(nn-rr))*(fact(rr))))
            prd=ncr*(pcx**rr)
            sumx1=0.0d0

            do 1204 ix2=0,(l1+l2)
              ii=ix2+ii1
              if (mod(ii,2) .eq. 0) then
                sumx1=sumx1+coeff(ix2,l1,l2,pax,pbx)* &
                     & (dfac(ii-1)/((two*gma)**(ii/2)))      
              endif        
1204        enddo
            sumx2=sumx2+(prd*sumx1)
1203      enddo

        summx = sumx2

        elseif (l3 .eq. 0) then
                
          do 1205 ix=0,((l1+l2)/2)
          summx=summx+coeff((2*ix),l1,l2,pax,pbx)* &
               & (dfac(2*(ix)-1)/((two*gma)**(ix)))
1205      enddo

        endif

!-------------------------------------------------------------------!
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
        function sumy(l1,l2,alp1,alp2,pax,pbx)
!-------------------------------------------------------------------!
        
        use constants

        implicit none
        integer(kind=4)         :: ix,l1,l2
        real(kind=8)            :: sumy,alp1,alp2
        real(kind=8)            :: pax,pbx,gma,dfac,coeff

        gma=alp1+alp2
        sumy=0.0d0

        do 1203 ix=0,((l1+l2)/2)
          sumy=sumy+coeff((2*ix),l1,l2,pax,pbx)* &
               & (dfac(2*(ix)-1)/((two*gma)**(ix)))   
1203    enddo
        
!-------------------------------------------------------------------!
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
        function summy(l1,l2,l3,alp1,alp2,pax,pbx,pcx)
!-------------------------------------------------------------------!
        
        use constants

        implicit none
        integer(kind=4)         :: ix,ix1,ix2,l1,l2,l3,ii,ii1
        integer(kind=4)         :: rr,nn,fact,ncr
        real(kind=8)            :: summy,alp1,alp2
        real(kind=8)            :: pax,pbx,pcx,gma,dfac,coeff
        real(kind=8)            :: sumy1,sumy2,prd

        gma=alp1+alp2
        sumy1=0.0d0; sumy2=0.0d0; summy=0.0d0

        if (l3 .ge. 1) then 
          nn=l3  

          do 1203 ix1=1,(l3+1)  
            rr=ix1-1 
            ii1=l3-rr
            ncr=((fact(nn))/((fact(nn-rr))*(fact(rr))))
            prd=ncr*(pcx**rr)
            sumy1=0.0d0

            do 1204 ix2=0,(l1+l2)
              ii=ix2+ii1
              if (mod(ii,2) .eq. 0) then
                sumy1=sumy1+coeff(ix2,l1,l2,pax,pbx)* &
                     & (dfac(ii-1)/((two*gma)**(ii/2)))      
              endif        
1204        enddo
              sumy2=sumy2+(prd*sumy1)
1203      enddo

        summy = sumy2

        elseif (l3 .eq. 0) then
                
          do 1205 ix=0,((l1+l2)/2)
          summy=summy+coeff((2*ix),l1,l2,pax,pbx)* &
               & (dfac(2*(ix)-1)/((two*gma)**(ix)))
1205      enddo

        endif

!-------------------------------------------------------------------!
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
        function sumz(l1,l2,alp1,alp2,pax,pbx)
!-------------------------------------------------------------------!
        
        use constants

        implicit none
        integer(kind=4)         :: ix,l1,l2
        real(kind=8)            :: sumz,alp1,alp2
        real(kind=8)            :: pax,pbx,gma,dfac,coeff

        gma=alp1+alp2
        sumz=0.0d0

        do 1203 ix=0,((l1+l2)/2)
          sumz=sumz+coeff((2*ix),l1,l2,pax,pbx)* &
               & (dfac(2*(ix)-1)/((two*gma)**(ix)))   
1203    enddo
        
!-------------------------------------------------------------------!
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
        function summz(l1,l2,l3,alp1,alp2,pax,pbx,pcx)
!-------------------------------------------------------------------!
        
        use constants

        implicit none
        integer(kind=4)         :: ix,ix1,ix2,l1,l2,l3,ii,ii1
        integer(kind=4)         :: rr,nn,fact,ncr
        real(kind=8)            :: summz,alp1,alp2,prd
        real(kind=8)            :: pax,pbx,pcx,gma,dfac,coeff
        real(kind=8)            :: sumz1,sumz2

        gma=alp1+alp2
        sumz1=0.0d0; sumz2=0.0d0; summz=0.0d0

        if (l3 .ge. 1) then 
          nn=l3  

          do 1203 ix1=1,(l3+1)  
            rr=ix1-1 
            ii1=l3-rr
            ncr=((fact(nn))/((fact(nn-rr))*(fact(rr))))
            prd=ncr*(pcx**rr)
            sumz1=0.0d0

            do 1204 ix2=0,(l1+l2)
              ii=ix2+ii1
              if (mod(ii,2) .eq. 0) then
                sumz1=sumz1+coeff(ix2,l1,l2,pax,pbx)* &
                     & (dfac(ii-1)/((two*gma)**(ii/2)))      
              endif        
1204        enddo
              sumz2=sumz2+(prd*sumz1)
1203      enddo

        summz = sumz2

        elseif (l3 .eq. 0) then
                
          do 1205 ix=0,((l1+l2)/2)
          summz=summz+coeff((2*ix),l1,l2,pax,pbx)* &
               & (dfac(2*(ix)-1)/((two*gma)**(ix)))
1205      enddo

        endif

!-------------------------------------------------------------------!
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!        
        function coeff(i,l1,l2,a,b)

        use constants        
        implicit none
        real(kind=8)    :: coeff,a,b
        integer(kind=4) :: i,l1,l2,ia,ib,fact,ix

        ia=max(0,i-l2)
        ib=min(i,l1)        

        coeff=zero

        do ix = ia,ib
          coeff=coeff+((fact(l1)*fact(l2)*(a**(l1-ix)))* &
               & (b**(l2-i+ix))) / &
               & (fact(ix)*fact(l1-ix)* &
               & fact(i-ix)*fact(l2-i+ix))    

        enddo
        
        end
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!        
        function fact(num)

        integer(kind=4) :: num,fact,ix

        fact=1

        if (num .le. 1) then
          fact=1
        else
          do ix = 1,num
            fact=fact*ix
          enddo      
        endif        

        end
!-------------------------------------------------------------------!        

