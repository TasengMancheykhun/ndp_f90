!     --------------------------------------------------------------------------
      module inpar
!     --------------------------------------------------------------------------
        implicit none
        
        real*8, parameter ::  zero = 0.0000000d+00,  one = 1.00000000d+00
        real*8, parameter :: quatr = 0.2500000D+00, half = 0.50000000D+00
        real*8, parameter ::   two = 2.0000000D+00, thre = 3.00000000D+00
        real*8, parameter ::  four = 4.0000000D+00
        real*8, parameter ::  pi = 3.1415926535897932D+00, twpi = two*pi
       
        real*8, parameter :: ev = 27.2113846900D+00
       
        complex*16, parameter ::   io = (zero, one)
        complex*16, parameter :: aval = (one, zero), bval = (zero, zero)
      end
!     -------------------------------------------------------------------------


!     -------------------------------------------------------------------------
      program main
!     -------------------------------------------------------------------------
        use inpar
        implicit none
        
        integer(kind=4), parameter :: n = 69
        real(kind=8)               :: smat(n,n), xmat(n,n), xmti(n,n),T
        real(kind=8)               :: coef(n,n), coftr(n,n), cofti(n,n)
        complex(kind=16)           :: S 

        real(kind=8)               :: dzr(n,n), dzi(n,n)
        integer(kind=4)            :: i1,i2,i3,j1
       

        complex*16 :: coft(n,n), tcoft(n,n)

        open(15,file = 'onelec.txt')
        open(16,file = './mat/hmt145br.txt')
        open(17,file = './mat/hmt145bi.txt')
        do i1 = 1, n
          do i2=1,i1 
            read(15,*) smat(i1,i2),T,T
                       smat(i2,i1) = smat(i1,i2)

            read(16,*) dzr(i1,i2)
                       dzr(i2,i1) = dzr(i1,i2)

            read(17,*) dzi(i1,i2)
                       dzi(i2,i1) = dzi(i1,i2)
          enddo
        enddo
        close(16)
        close(15)
        close(17)

        call XMAT_CAL_ (smat,n,xmat,xmti)
       
        open(11,file="coef.txt")
          do i1 = 1, n
            read(11,*) (coef(i1,i2),i2 = 1,n)
          enddo
        close(11)
       
        xmat = zero; xmat = matmul(xmti, coef)
!       coft = matmul(transpose(xmat),xmat)
       
        coftr = matmul(transpose(xmat),matmul(dzr,xmat))
        cofti = matmul(transpose(xmat),matmul(dzi,xmat))
        
!        open(32,file="mohmt.txt")  
!        do i1 = 1, n
!          write(32,*) cmplx(i1,i1), coftr(i1,i1) + io*cofti(i1,i1)
!          write(32,'(*(f8.3))') (coft(i1,i2),i2=1,n)
!        enddo
!        close(32)
        
! ------------- Taseng ------------!        
        do i1 = 1, n
          do j1 = 1, n
            if (coftr(i1,j1) .lt. 1E-10) then
              coftr(i1,j1) = 0.00d+00
            endif
          enddo
        enddo
        
        do i1 = 1, n
          do j1 = 1, n
            if (cofti(i1,j1) .lt. 1E-10) then
              cofti(i1,j1) = 0.00d+00
            endif
          enddo
        enddo
!-----------------------------------!
     
      coft = coftr + io*cofti

      tcoft= transpose(coft)
    
      open(122,file = './result/diagcoft.txt')
      do i1 = 1, n
        write(122,*) cmplx(i1,i1), real(coft(i1,i1)), aimag(coft(i1,i1))
      enddo
      close(122)



      open(122,file='./result/coft.txt')
      open(123,file="./result/tcoft.txt")
      do i1 = 1, n
        do j1 = 1, n
          write(122,*) cmplx(i1,j1), real(coft(i1,j1)), aimag(coft(i1,j1))
          write(123,*) cmplx(i1,j1), real(tcoft(i1,j1)), aimag(tcoft(i1,j1))
        enddo
      enddo
      close(122)
      close(123)



      end
!     ---------------------------------------------------------------------!




!     ---------------------------------------------------------------------!
      SUBROUTINE XMAT_CAL_ (smat, norb, xmat, xmti)
!     ---------------------------------------------------------------------!

      USE inpar
      IMPLICIT NONE

       REAL*8, allocatable :: sevl(:), exsvl(:,:), umat(:,:)

       real*8, intent(in)  :: smat(norb, norb)
       real*8, intent(out) :: xmat(norb, norb), xmti(norb, norb)

       INTEGER*4  :: norb, IX, IY, nelec

      allocate(umat(norb, norb), sevl(norb))

      sevl = zero; umat = zero; umat = smat
      call CALL_DSYEV_ (umat, norb, sevl)

      allocate(exsvl(norb, norb))

      exsvl = zero
        DO ix = 1, norb
          exsvl(ix,ix) = (sevl(ix))**(-half)
        ENDDO

         xmat = zero
       do ix = 1,  norb
        do iy = 1,  norb
         xmat(ix,iy) = exsvl(ix,ix)*umat(iy,ix)
        enddo
       enddo

        exsvl = zero; exsvl = xmat; xmat = zero
        call dgemm('N', 'N',  norb,  norb,  norb, one, umat,  norb, &
                     exsvl,  norb, zero, xmat,  norb)

         exsvl = zero
        DO ix = 1, norb
          exsvl(ix,ix) = (sevl(ix))**half
        ENDDO

         xmti = zero
       do ix = 1,  norb
        do iy = 1,  norb
         xmti(ix,iy) = exsvl(ix,ix)*umat(iy,ix)
        enddo
       enddo

        exsvl = zero; exsvl = xmti; xmti = zero
        call dgemm('N', 'N',  norb,  norb,  norb, one, umat,  norb, &
                     exsvl,  norb, zero, xmti,  norb)

      deallocate(umat, sevl, exsvl)

      END
!     --------------------------------------------------------------------!


!     --------------------------------------------------------------------!
      SUBROUTINE CALL_DSYEV_ (H, NDIM, eigvals)
!     --------------------------------------------------------------------!
       IMPLICIT NONE

       CHARACTER (LEN = 1), PARAMETER  :: JOBZ="V"
       CHARACTER (LEN = 1), PARAMETER  :: UPLO="L"

       REAL    (KIND = 8), ALLOCATABLE :: WORK(:)

       INTEGER (KIND = 4) :: NDIM,INFO
       INTEGER (KIND = 4) :: LDA, LDWORK

       REAL    (KIND = 8) :: H(NDIM,NDIM)
       REAL    (KIND = 8) :: EIGVALS(NDIM)

        LDA = NDIM
        LDWORK = 3*NDIM - 1

        ALLOCATE(WORK(LDWORK))
        CALL DSYEV (JOBZ, UPLO, NDIM, H, LDA, &
                          EIGVALS, WORK, LDWORK, INFO)
        DEALLOCATE(WORK)

      END


