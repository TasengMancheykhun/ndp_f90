        SUBROUTINE gptvalues      

          USE atom
          USE basis
          USE gpt

          IMPLICIT NONE 

          INTEGER*4 :: ia, ib, ic, id, ix, iy, ii
          INTEGER*4 :: in1, en1, in2, en2, in3, en3, in4, en4
          REAL*8    :: g1, gamm, ab2  

          !--- allocating the arrays

          ALLOCATE(npgalp(npg), npgnorm(npg), npgcff(npg))
          ALLOCATE(pxa(npg), pya(npg), pza(npg))  
          ALLOCATE(pn(npg), pl(npg), pm(npg), spd(npg))

          npgalp = 0.00d+00; npgnorm = 0.00d+00; npgcff = 0.00d+0
          pxa = 0.00d+00; pya = 0.00d+00; pza = 0.00d+00; pn = 0; pl = 0; pm = 0

          ! Assigning the values to the arrays required 

          ii = 0

          do 801 ia = 1,natom
            in1=basat1(atid(ia)); en1=basat2(atid(ia))
            do 802 ib = in1,en1
              in2=func1(ib); en2=func2(ib)
              do 803 ic = in2,en2
                ii = ii + 1
                npgalp(ii)=pgalp(ic)
                npgnorm(ii)=normnpg(ic)
                npgcff(ii)=npgcont(ic)
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
803           enddo
802         enddo
801       enddo 


 

!          DEALLOCATE(npgnorm, npgcff)
!          DEALLOCATE(spd)
   
        END SUBROUTINE gptvalues

        ! ------------------------------------------------------------- !      




        ! ------------------------------------------------------------- !
        !
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
        ! The subroutine provides : 
        ! the number of atoms (natom),      
        ! the number of primitive gaussians (npg), 
        ! the number of cartesian gaussian functions (nbas), 
        ! total number of functions (totbas) and 
        ! total number of basis set shell(shl). 
        !
        ! Also, the exponents(basalp) and contraction coefficients(bascff) 
        ! are stored along with the 
        ! index maps in 1d arrays of each shell (is1 and is2)      
        !      
        !---------------------------------------------------------------!

        SUBROUTINE readatombasis
           
          USE constants  
          USE atom
          USE basis

          IMPLICIT NONE

          INTEGER*4         :: ia, ib, ic, c, cn, iat, ii, ik, ij
          CHARACTER(LEN=15) :: dum, space
          CHARACTER(LEN=4)  :: at
          REAL*8            :: a, chg, x, y, z, expo, coef   
          REAL*8            :: ch1, ch2, ch3, ch4, pnorm, dfac
          INTEGER*4         :: mval, coun, total, in1, en1, in2, en2
          INTEGER*4         :: in3, en3, in4, en4, itot, atbas, atnpg  
          INTEGER*4         :: ibas, ipg1, ipg2

          space = "        "

          natom = 0

          !---Section to read count natom,npg,nbas,totbas,shell----------!

          OPEN(1,FILE='input.dat')
          READ(1,'(a6)') dum(1:6)
          IF (dum(1:6) .eq. " $DATA") THEN
            READ(1,*) dum
            READ(1,*) dum
          ELSE 
            WRITE(*,'("Incorrect file formating, Check")')        
          STOP
          ENDIF


          npg = 0; nbas = 0; totbas = 0; shl = 0

          READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z

          DO 101 WHILE(dum(1:10) .ne. " $END     ")          
            DO 102 WHILE(dum(1:8) .ne. space)
              READ(1,'(a5,i7)') dum, c

              DO 10 ia = 1, c
                READ(1,*) cn, expo, coef
10            ENDDO

              totbas = totbas + c

              IF (dum(1:4) .eq. "   S") THEN
                nbas = nbas + 1
                npg = npg + (c*1)
                shl = shl + 1

              ELSEIF (dum(1:4) .eq. "   P") THEN
                nbas = nbas + 3
                npg = npg + (c*3)
                shl = shl + 1

              ELSEIF (dum(1:4) .eq. "   D") THEN 
                nbas = nbas + 6
                npg = npg + (c*6)
                shl = shl + 1

              ELSEIF (dum(1:4) .eq. "   F") THEN 
                nbas = nbas + 10
                npg = npg + (c*10)
                shl = shl + 1

              ENDIF
             
102         ENDDO

            natom = natom + 1

            READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z

101       ENDDO
          CLOSE(1)


          ! -------- Print to check if the count is correct -------------!

          ! write(*,*) "Number of atoms ", natom
          ! write(*,*) "Number of primitives", npg
          ! write(*,*) "Number of basis functions", nbas
          ! write(*,*) "Number of functions", totbas
          ! write(*,*) "Number of basis set shells", shl

          !---------- Allocate dimensions to arrays ---------------------!

          ALLOCATE(ch(natom), xa(natom), ya(natom), za(natom), atnm(natom))

          !----------- Initialization of arrays -------------------------!

          OPEN(1,file = 'input.dat')
          READ(1,'(a6)') dum(1:6)

          IF (dum(1:6) .eq. " $DATA") THEN
            READ(1,*) dum
            READ(1,*) dum
          else
            WRITE(*,'("Incorrect file formating, Check")')
          STOP
          ENDIF

          iat = 0  
          READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z 

          DO 104 WHILE(dum(1:10) .ne. " $END     ")

            iat = iat + 1

            atnm(iat) = dum(1:4); ch(iat) = chg
            xa(iat) = x; ya(iat) = y; za(iat) = z

            DO 105 WHILE(dum(1:8) .ne. space)
            READ(1,'(a5,i7)') dum, c
              DO 11 ia = 1, c
                READ(1,*) cn, expo, coef
11            ENDDO
105         ENDDO

            READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z 
           
104       ENDDO       
          CLOSE(1)





           
          !--- Checking the type of atoms and their numbers

          ! atid = Array that stores ID of a particular atom.
          !    For eg. : ID is given as 1 for one type of atom
          !              and 2 for another type of atom, so on.    

          ALLOCATE(atid(natom))

          atid = 0; ii = 0; ik = 0

          DO 106 ia = 1, natom
            ii = 1 
            DO 107 ib = ia + 1, natom
              IF ( atnm(ia) .eq. atnm(ib) ) THEN
                ii = ii + 1                   
              ENDIF
107         ENDDO

            IF ( (ii .gt. 1) .and. (atid(ia) .eq. 0) ) THEN
              ik = ik + 1 

              DO 12 ib = ia + 1 , natom
                IF ( atnm(ia) .eq. atnm(ib) ) THEN
                  atid(ia) = ik
                  atid(ib) = ik
                ENDIF
12            ENDDO

            ELSEIF ( (ii .eq. 1) .and. (atid(ia) .eq. 0) ) THEN         
              ik = ik + 1

              atid(ia) = ik
            ENDIF
106       ENDDO


          !---- Number of times an atom appears (stored in nid array) ----!

          !  nid = Array that stores no. of counts of a particular type of atom
          ! nmid = Array that stores name of all type of atoms

          mval = maxval(atid, natom) 

          ALLOCATE( nid(mval), nmid(mval) )

          atid = 0; ii = 0; ik = 0; ij = 0

          DO 108 ia = 1, natom
            ii = 1 
            DO 109 ib = ia + 1, natom
              IF ( atnm(ia) .eq. atnm(ib) ) THEN
                ii = ii + 1                   
              ENDIF
109         ENDDO
          
            IF ( (ii .gt. 1) .and. (atid(ia) .eq. 0) ) THEN  
              ik = ik + 1; ij = ij + 1       

              DO 13 ib = ia + 1, natom
                IF (atnm(ia) .eq. atnm(ib)) THEN 
                  atid(ia) = ik
                  atid(ib) = ik
                ENDIF
13            ENDDO

              nid(ij) = ii  
              nmid(ij) = atnm(ia) 
              !write(*,*) nmid(ij),nid(ij)  

            ELSEIF ((ii .eq. 1) .and. (atid(ia) .eq. 0)) THEN

              ij = ij + 1      
              ik = ik + 1

              atid(ia) = ik
              nid(ij) = ii 
              nmid(ij) = atnm(ia) 
              !write(*,*) nmid(ij),nid(ij) 
              
            ENDIF
108       ENDDO



          !--- Storing the exponents and contraction coeff 
          !--- of only different atom types

          iat = 0 
          total = 0; in2 = 0; in1 = 0; en1 = 0
          in3 = 0; en3 = 0; atbas = 0; atnpg = 0 

          DO ik = 1, mval
            in1 = total + 1
             at = nmid(ik)
            ! write(*,*) "finding  " ,at 

            coun = 0; iat = 0

            OPEN(1, file = 'input.dat')
            READ(1,'(a6)') dum(1:6)

            IF (dum(1:6) .eq. " $DATA") THEN
              READ(1,*) dum
              READ(1,*) dum            
            ELSE
              WRITE(*,'("Incorrect file formating, Check")')  
              STOP    
            ENDIF      


            ! ibas = 0; iat = 0; ish = 0; total = 0
            READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z

            DO 110 WHILE(dum(1:10) .ne. " $END     ")
              iat = iat + 1

              atnm(iat) = dum(1:5); ch(iat) = chg
              xa(iat) = x; ya(iat) = y; za(iat) = z
              ! write(*,*) atnm(iat), dum(1:5), at

              IF ( at .eq. atnm(iat) ) THEN
                coun = coun + 1      
              ENDIF 

              DO 111 WHILE( dum(1:8) .ne. space )

                IF ( (at .eq. atnm(iat)) .and. (coun .eq. 1) ) THEN
                  READ(1,'(a5,i7)') dum, c

                  DO 14 ia = 1, c                  
                    READ(1,*) cn, expo, coef
                    ! write(*,*) cn, expo, coef
14                ENDDO

                  ! write(*,*) dum, c
                  IF (c .ne. 0) THEN

                    in2 = in2 + 1 
                    total = total + c      
                    en3 = total; in3 = en3 - (c - 1)
                    ! in1 = total - c + 1; en1 = total
                    ! write(*,*) in2, in3, en3
                  ENDIF

                  IF (dum(1:4) .eq. "   S") THEN
                    atbas = atbas + 1
                    atnpg = atnpg + (c*1)

                  ELSEIF (dum(1:4) .eq. "   P") THEN
                    atbas = atbas + 3
                    atnpg = atnpg + (c*3)

                  ELSEIF (dum(1:4) .eq. "   D") THEN
                    atbas = atbas + 6
                    atnpg = atnpg + (c*6)

                  ELSEIF (dum(1:4) .eq. "   F") THEN
                    atbas = atbas + 10
                    atnpg = atnpg + (c*10)

                  ENDIF

                ELSE                        
                  READ(1,'(a5,i7)') dum, c

                  DO 15 ia = 1, c
                    read(1,*) cn, expo, coef
15                ENDDO

                ENDIF    
111           ENDDO

              READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z
110         ENDDO
            CLOSE(1)
            
            en1 = total
          ENDDO
          !write(*,*) "writing ", atbas,atnpg 



          ! ------------------------------------------------------------------ !

          ALLOCATE( atsh1(mval), atsh2(mval), sh1(in2), sh2(in2) )
          ALLOCATE( alp(en1), cff(en1), spdf(in2) )
          ALLOCATE( basat1(mval), basat2(mval), func1(atbas), func2(atbas) )
          ALLOCATE( pgalp(atnpg), pgcff(atnpg) )
          ALLOCATE( spdfpg(atnpg), spdfbas(atbas) )

          ! ------------------------------------------------------------------ !
    
          alp = 0.0d0; cff = 0.0d0  

          iat = 0
          total = 0; in2 = 0; en2 = 0; in1 = 0; en1 = 0
          in3 = 0; en3 = 0; itot = 0; atbas = 0; atnpg = 0; ibas = 0

          DO 112 ik = 1, mval

            in1 = total + 1
             at = nmid(ik)
            en2 = in2 + 1
            in4 = atbas + 1
            ! write(*,*) "finding  " ,at(1:4) 

            coun = 0; iat = 0

            OPEN(1, file = 'input.dat')
            READ(1,'(a6)') dum(1:6)

            IF (dum(1:6) .eq. " $DATA") THEN
              READ(1,*) dum
              READ(1,*) dum
            ELSE 
              WRITE(*,'("Incorrect file formating, Check")')      
              STOP 
            ENDIF 


            ! ibas = 0; iat = 0; ish = 0; total = 0
            READ(1,'(a11,f7.1,3f18.10)') dum, chg, x, y, z

            DO 113 WHILE(dum(1:10) .ne. " $END     ")
              iat = iat + 1

              atnm(iat) = dum(1:4); ch(iat) = chg              
              xa(iat) = x; ya(iat) = y; za(iat) = z
              !write(*,*) atnm(iat)

              IF ( at .eq. atnm(iat) ) THEN
                coun = coun + 1      
              ENDIF        

              DO 114 WHILE( dum(1:8) .ne. space )
              
                IF ( (at .eq. atnm(iat)) .and. (coun .eq. 1) ) THEN
                  READ(1,'(a5,i7)') dum, c

                  DO 16 ia = 1, c
                    itot = itot + 1

                    READ(1,*) cn, alp(itot), cff(itot)
16                ENDDO

                  IF (c .ne. 0) THEN

                    in2 = in2 + 1 
                    total = total + c
                    en3 = total; in3 = en3 - (c - 1)

                    sh1(in2) = in3; sh2(in2) = en3
                    spdf(in2) = dum(4:4)
                    ! write(*,*) in2, spdf(in2), sh1(in2), sh2(in2)

                    ibas = atbas + 1

                    IF (dum(1:4) .eq. "   S") THEN
                      ! S Function
                      atbas = atbas + 1                   
                      DO ii = 1, c                   
                         atnpg = atnpg + 1                   
                         pgalp(atnpg) = alp(itot - c + ii)
                         pgcff(atnpg) = cff(itot - c + ii)
                        spdfpg(atnpg) = "S"                   
                      ENDDO
                     
                      spdfbas(atbas) = "S"                 
                      ipg1 = atnpg - c + 1 ; ipg2 = atnpg
                      func1(atbas) = ipg1 ; func2(atbas) = ipg2
                      !write(*,*) atbas,ipg1,ipg2


                    ELSEIF (dum(1:4) .eq. "   P") THEN
                      ! PX Function      
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "PX"
                      ENDDO

                      spdfbas(atbas) = "PX"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas,ipg1,ipg2
                      ! PY Function      
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "PY"
                      ENDDO

                      spdfbas(atbas) = "PY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2

                      ! write(*,*) atbas,ipg1,ipg2
                      ! PZ function   

                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "PZ"

                      ENDDO

                      spdfbas(atbas) = "PZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas,ipg1,ipg2

                    ELSEIF (dum(1:4) .eq. "   D") THEN
                      ! DXX Function      

                      atbas = atbas + 1 
                    
                      DO ii = 1, c                  
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                        spdfpg(atnpg) = "DXX"
                      ENDDO

                      spdfbas(atbas) = "DXX"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas,ipg1,ipg2
                      ! DYY Function      

                      atbas = atbas + 1

                      DO ii = 1, c                  
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "DYY"
                      ENDDO

                      spdfbas(atbas) = "DYY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas,ipg1,ipg2
                      ! DZZ Function      

                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "DZZ"
                      ENDDO

                      spdfbas(atbas) = "DZZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! DXY Function      
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "DXY"
                      ENDDO  

                      spdfbas(atbas) = "DXY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! DXZ Function      
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg+1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "DXZ"
                      ENDDO  

                      spdfbas(atbas) = "DXZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! DYZ Function      
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "DYZ"
                      ENDDO  

                      spdfbas(atbas) = "DYZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                    ELSEIF (dum(1:4) .eq. "   F") THEN
                      ! atbas = atbas + 10
                      ! atnpg = atnpg + (c*10)
                      ! FXXX Function
                      atbas = atbas + 1

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FXXX"
                      ENDDO 

                      spdfbas(atbas) = "FXXX"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FYYY Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FYYY"
                      ENDDO 

                      spdfbas(atbas) = "FYYY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FZZZ Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FZZZ"
                      ENDDO

                      spdfbas(atbas) = "FZZZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FXXY Function
                      atbas = atbas + 1

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FXXY"
                      ENDDO

                      spdfbas(atbas) = "FXXY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas,ipg1,ipg2

                      ! FXXZ Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FXXZ"
                      ENDDO  

                      spdfbas(atbas) = "FXXZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FYYX Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FYYX"
                      ENDDO  

                      spdfbas(atbas) = "FYYX"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FYYZ Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FYYZ"
                      ENDDO

                      spdfbas(atbas) = "FYYZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FZZX Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FZZX"
                      ENDDO  

                      spdfbas(atbas) = "FZZX"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      !write(*,*) atbas, ipg1, ipg2

                      ! FZZY Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FZZY"
                      ENDDO  

                      spdfbas(atbas) = "FZZY"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                      ! FXYZ Function
                      atbas = atbas + 1 

                      DO ii = 1, c
                        atnpg = atnpg + 1
                        pgalp(atnpg) = alp(itot - c + ii)
                        pgcff(atnpg) = cff(itot - c + ii)
                       spdfpg(atnpg) = "FXYZ"
                      ENDDO  

                      spdfbas(atbas) = "FXYZ"
                      ipg1 = atnpg - c + 1; ipg2 = atnpg
                      func1(atbas) = ipg1; func2(atbas) = ipg2
                      ! write(*,*) atbas, ipg1, ipg2

                    ENDIF
                    en4 = atbas
                    ! func1(in2) = ibas; func2(in2) = atbas
                    ! write(*,*) func1(in2), func2(in2)

                  ENDIF
                 
                ELSE  

                  READ(1,'(a5,i7)') dum, c

                  DO 17 ia = 1, c
                    READ(1,*) cn, expo, coef
17                ENDDO

                ENDIF
114           ENDDO

            READ(1,'(a11,f7.1,3f18.10)') dum,chg,x,y,z

113       ENDDO
          CLOSE(1)

          en1 = total
          atsh1(ik) = en2; atsh2(ik) = in2
          basat1(ik) = in4; basat2(ik) = en4

          ! write(*,*) "final", en2, in2
          ! write(*,*) "final functions", basat1(ik), basat2(ik)

112       ENDDO 

          !-------------------------------------------------------------------------!
          !---- NORMALIZATION CONSTANT CALCULATION AND CONTRACTION COEFFICIENTS ----!

          ALLOCATE(npgcont(atnpg), normnpg(atnpg))  
          in1 = 0; in2 = 0; in3 = 0;  
          en1 = 0; en2 = 0; en3 = 0;
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
117           enddo
116         enddo
115       enddo 

      xa = xa*1.8897261300d+00
      ya = ya*1.8897261300d+00
      za = za*1.8897261300d+00


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
   
      END SUBROUTINE readatombasis      
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
 
