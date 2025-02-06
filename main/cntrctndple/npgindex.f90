!--------------------------------------------------------------------!
      subroutine npgindex
!--------------------------------------------------------------------!
      use constants 
      use atom
      use basis
      use gpt
      use pdensitygrid 
      use mos
      implicit none
      integer(kind=4)   :: ia,ib,ic,ii1,ii2,ii,jj,in1,in2,en1,en2
      real(kind=8)      :: vv1

! ALLOCATING THE INDEXES OF npg's WHICH FORM THE BASIS FUNCTIONS 
      allocate(inde1(nbas),inde2(nbas),indat(nbas))
!---------------------------------------------------------------

!--------------------------------Mishu----------------
 allocate(pnpgnorm(npg),cffwnorm(npg),alpinv(npg),sqap(npg),sqapcub(npg),hfap(npg),osqapcub(npg))
!-----------------------------------------------------


! INDEXING OF BASIS FUNCTIONS TO npg's      
      ii=0; jj=0; inde1=0; inde2=0; indat=0 
      do 701 ia = 1,natom
        in1=basat1(atid(ia)); en1=basat2(atid(ia))
        do 702 ib=in1,en1
          in2=func1(ib); en2=func2(ib)
          jj=jj + 1
          ii1=ii + 1
          do 703 ic=in2,en2
            ii=ii+1
703       enddo
          ii2=ii
          inde1(jj)=ii1
          inde2(jj)=ii2
          indat(jj) = ia 
702     enddo
701   enddo          

! ALLOCATING THE ARRAYS TO STORE ALL THE PRIMITIVES EXPONENT
! ,THEIR CONTRACTION COEFFICIENTS,COORDINATES,N,L,M VALUES          
      allocate(npgalp(npg),npgnorm(npg),npgcff(npg))
      allocate(pxa(npg),pya(npg),pza(npg))
      allocate(pn(npg),pl(npg),pm(npg),spd(npg))

! ASSIGNING THE VALUES TO THE REQUIRED ARRAYS 
      npgalp=0.0d0; npgnorm=0.0d0; npgcff=0.0d0
      pxa=0.0d0; pya=0.0d0; pza=0.0d0; pn=0; pl=0; pm=0
      ii=0     

      do 801 ia = 1,natom
        in1=basat1(atid(ia)); en1=basat2(atid(ia))
        do 802 ib = in1,en1
          in2=func1(ib); en2=func2(ib)
          do 803 ic = in2,en2
            ii = ii + 1
            npgalp(ii)=pgalp(ic)
            npgnorm(ii)=normnpg(ic)
            npgcff(ii)=npgcont(ic)
!--------------------------Mishu------------------------------
           pnpgnorm(ii) = pnormnpg(ic)
           cffwnorm(ii) = npgcff(ii)*pnpgnorm(ii)
           alpinv(ii)   = 1.0d0/npgalp(ii) 
             sqap(ii)   = dsqrt(alpinv(ii))
           sqapcub(ii)  = (sqap(ii))**1.50d0
          osqapcub(ii)  = onesix*sqapcub(ii)
              hfap(ii)  = half*alpinv(ii)
!-------------------------------------------------------------
            pxa(ii)=xa(ia)
            pya(ii)=ya(ia)
            pza(ii)=za(ia)
            if (spdfpg(ic) .eq. "S") then
              pn(ii)=0; pl(ii)=0; pm(ii)=0
              spd(ii)="S"
            elseif (spdfpg(ic) .eq. "PX") then
              pn(ii)=1; pl(ii)=0; pm(ii)=0
              spd(ii)="PX"
            elseif (spdfpg(ic) .eq. "PY") then
              pn(ii)=0; pl(ii)=1; pm(ii)=0
             spd(ii)="PY"
           elseif (spdfpg(ic) .eq. "PZ") then
             pn(ii)=0; pl(ii)=0; pm(ii)=1
             spd(ii)="PZ"
           elseif (spdfpg(ic) .eq. "DXX") then
             pn(ii)=2; pl(ii)=0; pm(ii)=0
             spd(ii)="DXX"
           elseif (spdfpg(ic) .eq. "DYY") then
             pn(ii)=0; pl(ii)=2; pm(ii)=0
             spd(ii)="DYY"
           elseif (spdfpg(ic) .eq. "DZZ") then
             pn(ii)=0; pl(ii)=0; pm(ii)=2
             spd(ii)="DZZ"
           elseif (spdfpg(ic) .eq. "DXY") then
             pn(ii)=1; pl(ii)=1; pm(ii)=0
             spd(ii)="DXY"
           elseif (spdfpg(ic) .eq. "DXZ") then
             pn(ii)=1; pl(ii)=0; pm(ii)=1
             spd(ii)="DXZ"
           elseif (spdfpg(ic) .eq. "DYZ") then
             pn(ii)=0; pl(ii)=1; pm(ii)=1
             spd(ii)="DYZ"
           elseif (spdfpg(ic) .eq. "FXXX") then
             pn(ii)=3; pl(ii)=0; pm(ii)=0
             spd(ii)="FXXX"
           elseif (spdfpg(ic) .eq. "FYYY") then
             pn(ii)=0; pl(ii)=3; pm(ii)=0
             spd(ii)="FYYY"
           elseif (spdfpg(ic) .eq. "FZZZ") then
             pn(ii)=0; pl(ii)=0; pm(ii)=3
             spd(ii)="FZZZ"
           elseif (spdfpg(ic) .eq. "FXXY") then
             pn(ii)=2; pl(ii)=1; pm(ii)=0
             spd(ii)="FXXY"
           elseif (spdfpg(ic) .eq. "FXXZ") then
             pn(ii)=2; pl(ii)=0; pm(ii)=1
             spd(ii)="FXXZ"
           elseif (spdfpg(ic) .eq. "FYYX") then
             pn(ii)=1; pl(ii)=2; pm(ii)=0
             spd(ii)="FYYX"
           elseif (spdfpg(ic) .eq. "FYYZ") then
             pn(ii)=0; pl(ii)=2; pm(ii)=1
             spd(ii)="FYYZ"
           elseif (spdfpg(ic) .eq. "FZZX") then
             pn(ii)=1; pl(ii)=0; pm(ii)=2
             spd(ii)="FZZX"
           elseif (spdfpg(ic) .eq. "FZZY") then
             pn(ii)=0; pl(ii)=1; pm(ii)=2
             spd(ii)="FZZY"
           elseif (spdfpg(ic) .eq. "FXYZ") then
             pn(ii)=1; pl(ii)=1; pm(ii)=1
             spd(ii)="FXYZ"
           endif
803       enddo
802     enddo
801   enddo
 
      write(2,'(a32)') "READING OF ATOMS, BASIS COMPLETE"     
      write(2,'(a10)') "    "  
      
!      call chgcalcr 

      !--- Contraction coefficients normalized
!      ii=0; jj=0
!      do 804 ia = 1,natom
!        in1=basat1(atid(ia)); en1=basat2(atid(ia))
!        do 805 ib = in1,en1
!          in2=func1(ib); en2=func2(ib)
!          jj = jj + 1
!          do 806 ic = in2,en2
!            ii = ii + 1
!            npgcff(ii)=npgcff(ii)/dsqrt(ovpmat(jj,jj))
!806       enddo  
!805     enddo   
!804   enddo        

        !do ia = 1,natom
        !  write(*,*) xa(ia),ya(ia),za(ia)
        !enddo
!-------------------------------------------------------------------!
      end subroutine         
!-------------------------------------------------------------------!              
