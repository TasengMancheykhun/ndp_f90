      module constants
        implicit none
        real(kind=8)                 :: pi=3.14159265358979d0
        real(kind=8)                 :: two=2.0d0, zero=0.00d+00
        real(kind=8)                 :: thfr=0.750d0, half=0.50d0, four=4.0d0
        real(kind=8)                 :: bohr=1.889726130d0, thby2 = 3.00d+00/2.00d+00
        integer(kind=4)              :: fv=5,fo=4,one=1
      end module

      module filenames
        implicit none        
        character(len=40)            :: flinp,flout,flatombas   
        character(len=40)            :: flmor,flmocr,flmoci   
        character(len=40)            :: flrdencub,flcrdencub,flcidencub
        character(len=40)            :: mrc,denyn,grdyn,hessyn
      end module filenames

      module grid
        implicit none        
        real(kind=8)                 :: xmin,xmax,ymin,ymax,zmin,zmax
        real(kind=8)                 :: dx,dy,dz
        integer(kind=4)              :: nx,ny,nz,npt,ix,iy,iz
        real(kind=8),allocatable     :: den(:),gdxyz(:,:)
        real(kind=8),allocatable     :: lap(:),hxx(:),hyy(:),hzz(:)
        real(kind=8),allocatable     :: hxy(:),hyz(:),hxz(:)
        complex(kind=16),allocatable :: cden(:),cgdxyz(:,:)
        complex(kind=16),allocatable :: clap(:),chxx(:),chyy(:),chzz(:)
        complex(kind=16),allocatable :: chxy(:),chyz(:),chxz(:)
      end module       

      module atom
        implicit none     
        integer(kind=4)              :: natom
        integer(kind=4),allocatable  :: atid(:),nid(:)
        real(kind=8),allocatable     :: ch(:),xa(:),ya(:),za(:)      
        character(len=4),allocatable :: nmid(:)
      end module       

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
      end module        

      module gpt
        implicit none
        real(kind=8),allocatable     :: pxa(:),pya(:),pza(:)
        integer(kind=4),allocatable  :: pn(:),pl(:),pm(:)
        real(kind=8),allocatable     :: gpt1(:),gp(:,:),gpt3(:)
        real(kind=8),allocatable     :: npgnorm(:),npgcff(:),prdg(:)
        character(len=4),allocatable :: spd(:)
      end module

      module mos
        implicit none
        integer(kind=4)              :: imo,mu,nu,munu,ie,iele
        real(kind=8),allocatable     :: mocofr(:),mocofi(:)
        complex(kind=16),allocatable :: zmo(:),zpmat(:,:)
        real(kind=8),allocatable     :: pmatrix(:,:),orb(:),porb(:)
        real(kind=8),allocatable     :: dorbx(:),dorby(:),dorbz(:)
        real(kind=8),allocatable     :: dorbxx(:),dorbyy(:),dorbzz(:)
        real(kind=8),allocatable     :: dorbxy(:),dorbxz(:),dorbyz(:)
        real(kind=8),allocatable     :: bo(:,:),psi(:),coeff(:,:)
      end module mos

