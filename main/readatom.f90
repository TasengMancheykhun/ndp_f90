  module constants  
  implicit none
  real(kind=8)                 :: pi=3.14159265358979d0,two=2.0d0
  real(kind=8)                 :: thfr=0.750d0,half=0.50d0,four=4.0d0
  end module        


  module atom
  implicit none
  integer(kind=4)              :: natom 
  real(kind=8),allocatable     :: ch(:),xxa(:),yya(:),zza(:)
  integer(kind=4),allocatable  :: atid(:),nid(:)
  character(len=4),allocatable :: nmid(:)
  end module atom


  module basis
  implicit none
  integer(kind=4)              :: nbas,npg,totbas,shl
  integer(kind=4)              :: n,l,m,nlmsum    
  integer(kind=4),allocatable  :: is1(:),is2(:),atsh(:)
  integer(kind=4),allocatable  :: ib1(:),ib2(:),atpg(:)
  real(kind=8),allocatable     :: basalp(:),bascff(:)
  real(kind=8),allocatable     :: normnpg(:),npgcont(:)
  real(kind=8),allocatable     :: normcon(:),npgalp(:) 
  real(kind=8),allocatable     :: ovpnpg(:,:)
  character(len=4),allocatable :: atnm(:),spdf(:)
  character(len=4),allocatable :: spdfpg(:),spdfbas(:)
  integer(kind=4),allocatable  :: atsh1(:),atsh2(:)
  integer(kind=4),allocatable  :: sh1(:),sh2(:)
  real(kind=8),allocatable     :: alp(:),cff(:)
  integer(kind=4),allocatable  :: basat1(:),basat2(:)
  integer(kind=4),allocatable  :: func1(:),func2(:)
  real(kind=8),allocatable     :: pgalp(:),pgcff(:)
  integer(kind=4),allocatable  :: inde1(:),inde2(:),indat(:)
  end module basis


  module gpt  
  implicit none
  real(kind=8),allocatable     :: pxa(:),pya(:),pza(:)
  integer(kind=4),allocatable  :: pn(:),pl(:),pm(:)     
  real(kind=8),allocatable     :: gpt1(:),gp(:,:),gpt3(:)
  real(kind=8),allocatable     :: npgnorm(:),npgcff(:),prdg(:)
  character(len=4),allocatable :: spd(:) 
  end module 

   subroutine gptvalues
      use atom
      use basis
      use gpt
      implicit none  
      integer(kind=4)   :: ia,ib,ic,id,ix,iy,ii
      integer(kind=4)   :: in1,en1,in2,en2,in3,en3,in4,en4
      real(kind=8)      :: g1,gamm,ab2  

      !--- allocating the arrays 
      allocate(npgalp(npg),npgnorm(npg),npgcff(npg))
      allocate(pxa(npg),pya(npg),pza(npg))  
      allocate(pn(npg),pl(npg),pm(npg),spd(npg))

      npgalp=0.00d+00; npgnorm=0.00d+00; npgcff=0.00d+0
      pxa=0.00d+00; pya=0.00d+00; pza=0.00d+00; pn=0; pl=0; pm=0

     ! Assigning the values to the arrays required 
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
            pxa(ii)=xxa(ia)
            pya(ii)=yya(ia)
            pza(ii)=zza(ia)
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


 
!-------------------------------------------------------------------!
      deallocate(npgnorm,npgcff)
      deallocate(spd)
!   
      end subroutine gptvalues
!-------------------------------------------------------------------!      




!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
! Author - Nitin Kumar Singh || Date - 26/06/2020 ||
! Affiliation - Department of Chemical Sciences, IISER Mohali      
!      
! A subroutine which reads the atom name, charge, coordinates  
! from a input file called "input.dat" provided by the user.
!
! The "input.dat" file should be provided in a specific format.
! This format is exactly like that of a gamess .dat file. 
! To make input.dat file ---> run a gamess file and store its
! filename.dat file.       
!
! After calling, the subroutine provides the number of atoms(natom),      
! the number of primitive gaussians(npg),the number of cartesian
! gaussian functions(nbas),total number of functions(totbas) and 
! total number of basis set shell(shl). Also, the exponents(basalp)
! and contraction coefficients(bascff) are stored along with the 
! index maps in 1d arrays of each shell (is1 and is2)      
!      
!-------------------------------------------------------------------!
      subroutine readatombasis
!-------------------------------------------------------------------!              
      use constants  
      use atom
      use basis
 
      implicit none
      integer(kind=4)    :: ia,ib,ic,c,cn,iat,ii,ik,ij
      character(len=15)  :: dum,space
      character(len=4)   :: at
      real(kind=8)       :: a,chg,x,y,z,expo,coef   
      real(kind=8)       :: ch1,ch2,ch3,ch4,pnorm,dfac
      integer(kind=4)    :: mval,coun,total,in1,en1,in2,en2
      integer(kind=4)    :: in3,en3,in4,en4,itot,atbas,atnpg  
      integer(kind=4)    :: ibas,ipg1,ipg2
      space = "        "


      natom = 0
      !---Section to read count natom,npg,nbas,totbas,shell----------!
      open(1,file='input.dat')
      read(1,'(a6)') dum(1:6)
      if (dum(1:6) .eq. " $DATA") then
        read(1,*) dum
        read(1,*) dum
      else 
        write(*,'("Incorrect file formating, Check")')        
      stop
      endif
      npg = 0; nbas = 0; totbas = 0; shl = 0
      read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
      do 101 while(dum(1:10) .ne. " $END     ")
        do 102 while(dum(1:8) .ne. space)
          read(1,'(a5,i7)') dum,c
          do 10 ia = 1,c
            read(1,*) cn,expo,coef
10        enddo 
          totbas = totbas + c
          if (dum(1:4) .eq. "   S") then 
            nbas = nbas + 1
            npg = npg + (c*1)
            shl = shl + 1
          elseif (dum(1:4) .eq. "   P") then 
            nbas = nbas + 3
            npg = npg + (c*3)
            shl = shl + 1
          elseif (dum(1:4) .eq. "   D") then 
            nbas = nbas + 6
            npg = npg + (c*6)
            shl = shl + 1
          elseif (dum(1:4) .eq. "   F") then 
            nbas = nbas + 10
            npg = npg + (c*10)
            shl = shl + 1
          endif
102     enddo
        natom = natom + 1
        read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
101   enddo
      close(1)
      !---Display to check if the count is correct-------------------!
      !write(*,*) "Number of atoms ", natom
      !write(*,*) "Number of primitives", npg
      !write(*,*) "Number of basis functions", nbas
      !write(*,*) "Number of functions", totbas
      !write(*,*) "Number of basis set shells", shl
      !---Allocate dimensions to arrays -----------------------------!
      allocate(ch(natom),xxa(natom),yya(natom),zza(natom),atnm(natom))
      !---Initialization of arrays-----------------------------------!
      open(1,file='input.dat')
      read(1,'(a6)') dum(1:6)
      if (dum(1:6) .eq. " $DATA") then
      read(1,*) dum
      read(1,*) dum
      else
      write(*,'("Incorrect file formating, Check")')
      stop
      endif
      iat=0  
      read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z  
      do 104 while(dum(1:10) .ne. " $END     ")
        iat = iat + 1
        atnm(iat)=dum(1:4); ch(iat)=chg
        xxa(iat)=x; yya(iat)=y; zza(iat)=z
        do 105 while(dum(1:8) .ne. space)
        read(1,'(a5,i7)') dum,c
          do 11 ia = 1,c
            read(1,*) cn,expo,coef
11        enddo
105     enddo
        read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z   
104   enddo        
      close(1)
      !--- Checking the type of atoms and their numbers
      allocate(atid(natom)) 
      atid=0; ii=0; ik=0;ij=0
      do 106 ia=1,natom
        ii=1 
        do 107 ib=ia+1,natom
          if (atnm(ia) .eq. atnm(ib)) then
            ii=ii + 1                   
          endif
107     enddo
        if ((ii .gt. 1) .and. (atid(ia) .eq. 0)) then   
          ik=ik+1; ij=ij+1       
          do 12 ib=ia+1,natom
            if (atnm(ia) .eq. atnm(ib)) then 
              atid(ia)=ik
              atid(ib)=ik
            endif
12        enddo
        elseif ((ii .eq. 1) .and. (atid(ia) .eq. 0)) then   
          ij=ij+1      
          ik=ik+1
          atid(ia)=ik
        endif
106   enddo
      !--- number of times a atom appears (stored in nid array)
      mval=maxval(atid,natom) 
      allocate(nid(mval),nmid(mval))
      atid=0; ii=0; ik=0;ij=0
      do 108 ia=1,natom
        ii=1 
        do 109 ib=ia+1,natom
          if (atnm(ia) .eq. atnm(ib)) then
            ii=ii + 1                   
          endif
109     enddo
        if ((ii .gt. 1) .and. (atid(ia) .eq. 0)) then   
          ik=ik+1; ij=ij+1       
          do 13 ib=ia+1,natom
            if (atnm(ia) .eq. atnm(ib)) then 
              atid(ia)=ik
              atid(ib)=ik
            endif
13        enddo
          nid(ij)=ii  
          nmid(ij)=atnm(ia) 
          !write(*,*) nmid(ij),nid(ij)  
        elseif ((ii .eq. 1) .and. (atid(ia) .eq. 0)) then   
          ij=ij+1      
          ik=ik+1
          atid(ia)=ik
          nid(ij)=ii 
          nmid(ij)=atnm(ia) 
          !write(*,*) nmid(ij),nid(ij)  
        endif
108   enddo
      !--- storing the exponents and contraction coeff only different
      !--- atom types
      iat=0;  
      total=0; in2=0; in1=0;en1=0
      in3=0; en3=0; atbas=0; atnpg=0 
      do ik=1,mval
        in1=total+1
        at=nmid(ik)
        !write(*,*) "finding  " ,at 
        coun=0; iat=0
        open(1,file='input.dat')
        read(1,'(a6)') dum(1:6)
        if (dum(1:6) .eq. " $DATA") then
          read(1,*) dum
          read(1,*) dum
        else 
          write(*,'("Incorrect file formating, Check")')      
        endif        
        !ibas =0; iat = 0; ish = 0; total = 0
        read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
        do 110 while(dum(1:10) .ne. " $END     ")
          iat = iat + 1
          atnm(iat)=dum(1:5); ch(iat)=chg
          xxa(iat)=x; yya(iat)=y; zza(iat)=z
          !write(*,*) atnm(iat),dum(1:5),at
          if (at .eq. atnm(iat)) then
            coun=coun+1      
          endif        
          do 111 while(dum(1:8) .ne. space)
            if ((at .eq. atnm(iat)) .and. (coun .eq. 1)) then
              read(1,'(a5,i7)') dum,c
              do 14 ia = 1,c
                read(1,*) cn,expo,coef
                !write(*,*) cn,expo,coef
14            enddo
              !write(*,*) dum,c
              if (c .ne. 0) then
              in2=in2 + 1 
              total=total+c      
              en3=total; in3=en3-(c-1)
              !in1=total-c+1; en1=total
              !write(*,*) in2,in3,en3
              endif
              if (dum(1:4) .eq. "   S") then
                atbas = atbas + 1
                atnpg = atnpg + (c*1)
              elseif (dum(1:4) .eq. "   P") then
                atbas = atbas + 3
                atnpg = atnpg + (c*3)
              elseif (dum(1:4) .eq. "   D") then
                atbas = atbas + 6
                atnpg = atnpg + (c*6)
              elseif (dum(1:4) .eq. "   F") then
                atbas = atbas + 10
                atnpg = atnpg + (c*10)
              endif
            else  
              read(1,'(a5,i7)') dum,c
              do 15 ia = 1,c
                read(1,*) cn,expo,coef
15            enddo
            endif    
111       enddo
          read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
110     enddo
        close(1)
        en1=total
      enddo
      !write(*,*) "writing ", atbas,atnpg 
      !-------------------------------------------------------------!
      allocate(atsh1(mval),atsh2(mval),sh1(in2),sh2(in2))
      allocate(alp(en1),cff(en1),spdf(in2))
      allocate(basat1(mval),basat2(mval),func1(atbas),func2(atbas))
      allocate(pgalp(atnpg),pgcff(atnpg))
      allocate(spdfpg(atnpg),spdfbas(atbas))
      !-------------------------------------------------------------!
      !-------------------------------------------------------------!
      alp=0.0d0; cff=0.0d0  
      !-------------------------------------------------------------!  
      iat=0; total=0; in2=0; en2=0; in1=0; en1=0
      in3=0; en3=0; itot=0; atbas=0; atnpg=0; ibas=0
      do 112 ik=1,mval
        in1=total+1
        at=nmid(ik)
        en2=in2+1
        in4 = atbas+1
        !write(*,*) "finding  " ,at(1:4) 
        coun=0; iat=0
        open(1,file='input.dat')
        read(1,'(a6)') dum(1:6)
        if (dum(1:6) .eq. " $DATA") then
          read(1,*) dum
          read(1,*) dum
        else 
          write(*,'("Incorrect file formating, Check")')      
        endif        
        !ibas =0; iat = 0; ish = 0; total = 0
        read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
        do 113 while(dum(1:10) .ne. " $END     ")
          iat = iat + 1
          atnm(iat)=dum(1:4); ch(iat)=chg
          xxa(iat)=x; yya(iat)=y; zza(iat)=z
          !write(*,*) atnm(iat)
          if (at .eq. atnm(iat)) then
            coun=coun+1      
          endif        
          do 114 while(dum(1:8) .ne. space)
            if ((at .eq. atnm(iat)) .and. (coun .eq. 1)) then
              read(1,'(a5,i7)') dum,c
              do 16 ia = 1,c
                itot=itot+1
                read(1,*) cn,alp(itot),cff(itot)
16            enddo
              if (c .ne. 0) then
              in2=in2 + 1 
              total=total+c
              en3=total; in3=en3-(c-1)
              sh1(in2)=in3; sh2(in2)=en3
              spdf(in2)=dum(4:4)
              !write(*,*) in2,spdf(in2),sh1(in2),sh2(in2)
              ibas = atbas+1
              if (dum(1:4) .eq. "   S") then
                ! S Function
                atbas = atbas + 1
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="S"
                enddo
                spdfbas(atbas)="S"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
              elseif (dum(1:4) .eq. "   P") then
                ! PX Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="PX"
                enddo
                spdfbas(atbas)="PX"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! PY Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="PY"
                enddo
                spdfbas(atbas)="PY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! PZ function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="PZ"
                enddo  
                spdfbas(atbas)="PZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
              elseif (dum(1:4) .eq. "   D") then
                ! DXX Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DXX"
                enddo  
                spdfbas(atbas)="DXX"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! DYY Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DYY"
                enddo  
                spdfbas(atbas)="DYY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! DZZ Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DZZ"
                enddo  
                spdfbas(atbas)="DZZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! DXY Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DXY"
                enddo  
                spdfbas(atbas)="DXY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! DXZ Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DXZ"
                enddo  
                spdfbas(atbas)="DXZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! DYZ Function      
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="DYZ"
                enddo  
                spdfbas(atbas)="DYZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
              elseif (dum(1:4) .eq. "   F") then
                !atbas = atbas + 10
                !atnpg = atnpg + (c*10)
                ! FXXX Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FXXX"
                enddo  
                spdfbas(atbas)="FXXX"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FYYY Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FYYY"
                enddo  
                spdfbas(atbas)="FYYY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FZZZ Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FZZZ"
                enddo  
                spdfbas(atbas)="FZZZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FXXY Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FXXY"
                enddo  
                spdfbas(atbas)="FXXY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FXXZ Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FXXZ"
                enddo  
                spdfbas(atbas)="FXXZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FYYX Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FYYX"
                enddo  
                spdfbas(atbas)="FYYX"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FYYZ Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FYYZ"
                enddo  
                spdfbas(atbas)="FYYZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FZZX Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FZZX"
                enddo  
                spdfbas(atbas)="FZZX"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FZZY Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FZZY"
                enddo  
                spdfbas(atbas)="FZZY"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
                ! FXYZ Function
                atbas = atbas + 1 
                do ii=1,c
                  atnpg=atnpg+1
                  pgalp(atnpg)=alp(itot-c+ii)
                  pgcff(atnpg)=cff(itot-c+ii)
                  spdfpg(atnpg)="FXYZ"
                enddo  
                spdfbas(atbas)="FXYZ"
                ipg1=atnpg-c+1; ipg2=atnpg
                func1(atbas)=ipg1; func2(atbas)=ipg2
                !write(*,*) atbas,ipg1,ipg2
              endif
              en4=atbas
              !func1(in2)=ibas; func2(in2)=atbas
              !write(*,*) func1(in2),func2(in2)
              endif
            else  
              read(1,'(a5,i7)') dum,c
              do 17 ia = 1,c
                read(1,*) cn,expo,coef
17            enddo
            endif    
114       enddo
          read(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z
113     enddo
        close(1)
        en1=total
        atsh1(ik)=en2; atsh2(ik)=in2
        basat1(ik)=in4; basat2(ik)=en4
        !write(*,*) "final",en2,in2
        !write(*,*) "final functions",basat1(ik),basat2(ik)
112    enddo 
!-----------------------------------------------------------------------!
!--- NORMALIZATION CONSTANT CALCULATION AND CONTRACTION COEFFICIENTS----!

      allocate(npgcont(atnpg),normnpg(atnpg))  
      in1=0; in2=0; in3=0;  
      en1=0; en2=0; en3=0;
      do 115 ia = 1,mval
        in1 = basat1(ia); en1 = basat2(ia)
        !write(*,*) nmid(ia),in1,en1
        do 116 ib = in1,en1
          in2 = func1(ib); en2=func2(ib)
          !write(*,*) in2,en2
          do 117 ic = in2,en2
            if (spdfpg(ic) .eq. "S") then
              a=pgalp(ic)
              n=0; l=0; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "PX") then   
              a=pgalp(ic)
              n=1; l=0; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "PY") then   
              a=pgalp(ic)
              n=0; l=1; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "PZ") then   
              a=pgalp(ic)
              n=0; l=0; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DXX") then   
              a=pgalp(ic)
              n=2; l=0; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DYY") then   
              a=pgalp(ic)
              n=0; l=2; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DZZ") then   
              a=pgalp(ic)
              n=0; l=0; m=2
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DXY") then   
              a=pgalp(ic)
              n=1; l=1; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DXZ") then   
              a=pgalp(ic)
              n=1; l=0; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "DYZ") then   
              a=pgalp(ic)
              n=0; l=1; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FXXX") then   
              a=pgalp(ic)
              n=3; l=0; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FYYY") then   
              a=pgalp(ic)
              n=0; l=3; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FZZZ") then   
              a=pgalp(ic)
              n=0; l=0; m=3
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FXXY") then   
              a=pgalp(ic)
              n=2; l=1; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FXXZ") then   
              a=pgalp(ic)
              n=2; l=0; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FYYX") then   
              a=pgalp(ic)
              n=1; l=2; m=0
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FYYZ") then   
              a=pgalp(ic)
              n=0; l=2; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FZZX") then   
              a=pgalp(ic)
              n=1; l=0; m=2
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FZZY") then   
              a=pgalp(ic)
              n=0; l=1; m=2
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            elseif (spdfpg(ic) .eq. "FXYZ") then   
              a=pgalp(ic)
              n=1; l=1; m=1
              call normalize(a,n,l,m,pnorm)
              normnpg(ic)=pnorm
              npgcont(ic)=pgcff(ic)    
            endif               
117       enddo
116     enddo
115   enddo 

      xxa = xxa*1.8897261300d+00
      yya = yya*1.8897261300d+00
      zza = zza*1.8897261300d+00


!-------------------------------------------------------------------!


!      deallocate(atid) 
!      deallocate(nid,nmid)
!
!      deallocate(atsh1,atsh2,sh1,sh2)
!      deallocate(alp,cff,spdf)
!      deallocate(basat1,basat2,func1,func2)
!      deallocate(pgalp,pgcff)
!      deallocate(spdfpg,spdfbas)
! 
!      deallocate(ch,xa,ya,za,atnm)
!      deallocate(npgcont,normnpg)
   
      end subroutine readatombasis      
!-------------------------------------------------------------------!             



!-------------------------------------------------------------------!
      function dfac(num)
!-------------------------------------------------------------------!
      implicit none
      real(kind=8)  :: dfac,prod
      integer(kind=4)  :: i,n,num

      n = num
      prod = 1.0d0
      do 117 while (n .gt. 0)
        prod = prod*n
        n = n-2
117   enddo
      dfac = prod

!-------------------------------------------------------------------!
      end function
!-------------------------------------------------------------------!




!-------------------------------------------------------------------!
      subroutine normalize(a0,n0,l0,m0,pnorm0)
!-------------------------------------------------------------------!
      use constants
      implicit none
      integer(kind=4)   :: n0,l0,m0,nlmsum
      real(kind=8)      :: a0,pnorm0
      real(kind=8)      :: frac1,frac2,numo,demo,prod,dfac

      frac1=1.0d0; frac2=1.0d0; prod=1.0d0
      numo=1.0d0; demo=1.0d0; pnorm0=1.0d0

      frac1=((two*a0)/pi)
      frac1 = frac1**thfr
      nlmsum=n0+l0+m0
      numo=(four*a0)**(nlmsum)
      demo=(dfac(2*n0-1))*(dfac(2*l0-1))*(dfac(2*m0-1))
      frac2=(numo/demo)**(half)
      prod=frac1*frac2
      pnorm0=prod

!-------------------------------------------------------------------!
      end subroutine
!-------------------------------------------------------------------!
 
