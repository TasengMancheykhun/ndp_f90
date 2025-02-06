!--------Module-containing-constants-and-variables--------------------!
        module baspara

        implicit none
        real(kind=8)  :: thtw=1.50d0
        real(kind=8)  :: thr=3.0d0
        real(kind=8)  :: rthr=1.73205080756888d0
        real(kind=8)  :: a,b,b2,ab2,ax,bx,ay,by,az,bz
        real(kind=8)  :: ppx,ppy,ppz,ccx,ccy,ccz
        real(kind=8)  :: prodab,sumab,sumab2,sumab3
        real(kind=8)  :: diffx,diffy,diffz,na,nb
        real(kind=8)  :: piab,exppd,ovp,ssoo
        real(kind=8)  :: frac,frac1,frac2,frac3,frac4

        end module
!---------------------------------------------------------------------!

      module constants
      implicit none
      !real(kind=8)                 :: pi=3.141592653589793238460d0,two=2.0d0
      real(kind=8)                 :: pi=3.14159265358979323844D+00,two=2.0d0
      real(kind=8)                 :: thfr=0.750d0,half=0.50d0,four=4.0d0
      !real(kind=8)                 :: bohr=1.889726130d0,cf=2.871234000188190d0
      real(kind=8)                 :: bohr=(1.0d0/0.52917724924d0) !1.88972598857892d0
      real(kind=8)                 :: cf=2.871234000188190d0
      real(kind=8)                 :: onei=0.1250d0!,fvbyth=1.6666666666666667
      real(kind=8)                 :: epsl=2.87d-05
      real(kind=8)                 :: oner=1.0d0,fvbyth=(5.0d0/3.0d0)
      integer(kind=4)              :: fv=5,fo=4,one=1
      real(kind=8)                 :: quart=0.250d0,thby2=1.50d0,zero=0.0d0
      complex(kind=8)              :: zio =(0.0d0,1.0d0), zio3 = (0.0d0,-1.0d0)
      real(kind=8)                 :: onesix = 0.16666666666666666666666666666666666670d0
      end module

      module filenames
      implicit none
      character(len=40)            :: flinp,flout,flatombas
      character(len=40)            :: dummy,morc,flmor,flmocr,flmoci
      character(len=40)            :: prop,prop2,flelfcub
      character(len=40)            :: flrdencub,flcrdencub,flcidencub,flmo,flmo2
      character(len=4)             :: gpath
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
      real(kind=8),allocatable     :: pnormnpg(:)
      end module

      module gpt
      implicit none
      real(kind=8),allocatable     :: pxa(:),pya(:),pza(:)
      integer(kind=4),allocatable  :: pn(:),pl(:),pm(:),nlm(:)
      real(kind=8),allocatable     :: gpt1(:),gp(:,:),gpt3(:),pgcoff(:)
      real(kind=8),allocatable     :: npgnorm(:),npgcff(:),prdg(:)
      character(len=4),allocatable :: spd(:)
      real(kind=8),allocatable     :: pnpgnorm(:)
      real(kind=8),allocatable     :: momx(:,:),momy(:,:),momz(:,:)
      real(kind=8),allocatable     :: momxx(:,:),momyy(:,:),momzz(:,:)
      real(kind=8),allocatable     :: momxy(:,:),momyz(:,:),momxz(:,:)
      end module

      module mos
      implicit none
      integer(kind=4)              :: imo,mu,nu,munu,ie,iele,nocc
      real(kind=8),allocatable     :: mocofr(:),mocofi(:)
      complex(kind=8),allocatable :: zmo(:),zpmat(:,:) !,zmo2(:,:)
      real(kind=8),allocatable     :: pmatrix(:,:),orb(:),porb(:)
      real(kind=8),allocatable     :: ovpmat(:,:),povp(:,:)
      real(kind=8),allocatable     :: ps(:,:),sps(:,:),shalf(:,:)
      real(kind=8),allocatable     :: dorbx(:),dorby(:),dorbz(:)
      real(kind=8),allocatable     :: dorbxx(:),dorbyy(:),dorbzz(:)
      real(kind=8),allocatable     :: dorbxy(:),dorbxz(:),dorbyz(:)
      real(kind=8),allocatable     :: bo(:,:),psi(:),coeff(:,:),occ(:)
      real(kind=8),allocatable     :: gdcr(:,:),gdci(:,:),vd(:,:)
      complex(kind=8),allocatable  :: hcxx(:),hcyy(:),hczz(:)
      complex(kind=8),allocatable  :: hcxy(:),hcyz(:),hcxz(:)
      real(kind=8),allocatable     :: pobp(:)
      complex(kind=8),allocatable :: zpobp(:),zpox(:),zpoy(:),zpoz(:), zdxx(:),zdyy(:),zdzz(:),zdxy(:),zdyz(:),zdxz(:) 
      complex(kind=8),allocatable :: orbp(:),zorbx(:),zorby(:),zorbz(:),pox(:),poy(:),poz(:)
      complex(kind=8),allocatable :: zorbxx(:),zorbyy(:),zorbzz(:),zorbxy(:),zorbyz(:),zorbxz(:)  
      end module mos


      module grid
      implicit none
      real(kind=8)                 :: xmin,xmax,ymin,ymax,zmin,zmax
      real(kind=8)                 :: dx,dy,dz,stp
      integer(kind=4)              :: nx,ny,nz,npt,ix,iy,iz,npts,pts1
      integer(kind=4)              :: gcir,ncp,pts,cpm1,cpm3,cpp1,cpp3
      real(kind=8),allocatable     :: den(:),gdxyz(:,:),elf(:)
      real(kind=8),allocatable     :: lap(:),hxx(:),hyy(:),hzz(:)
      real(kind=8),allocatable     :: hxy(:),hyz(:),hxz(:)
      complex(kind=8),allocatable :: cden(:),cgdxyz(:,:)
      complex(kind=8),allocatable :: clap(:),chxx(:),chyy(:),chzz(:)
      complex(kind=8),allocatable :: chxy(:),chyz(:),chxz(:),lapc(:)
      end module

      module calctime
      implicit none
      real(kind=8)                 :: ts,tf     
      end module     








      module pfilenames
      implicit none        
      character(len=14)            :: flemdcub,flemdcubx,flemdcuby,flemdcubz,flemdcubxx,flemdcubyy,flemdcubzz,flcpemd
      character(len=40)            :: pdenyn,pgrdyn,phessyn
      end module pfilenames

      module pdensitygrid
      implicit none
      real(kind=8)                 :: px,py,pz,px2,py2,pz2,px3,py3,pz3,p3,pxpy,pypz,pxpz,pxpypz
      integer(kind=4)              :: npx,npy,npz,ipx,ipy,ipz,iorb,maxorb,kk!,p1,p2,npts
      real(kind=8)                 :: pxmin, pymin, pzmin, pxmax, pymax, pzmax  
      real(kind=8)                 :: dpx,dpy,dpz
      real(kind=8),allocatable     :: pden(:),pdenx(:),pdeny(:),pdenz(:),pdenxx(:),pdenyy(:),pdenzz(:),pgradsum(:),mop(:),plap(:)
      real(kind=8),allocatable     :: pdenxy(:),pdenyz(:),pdenxz(:),alpinv(:),cffwnorm(:),sqap(:),sqapcub(:)
      real(kind=8),allocatable     :: osqapcub(:),hfap(:),orbp2(:),orbpimag(:) 
      real(kind=8),allocatable     :: orsum(:,:),oisum(:,:),osum(:,:) 
      real(kind=8)      :: pxAx,pyAy,pzAz,pA,ap,pna,cf,Ax,Ay,Az
      real(kind=8)      :: sqrtap,pxsqap,pysqap,pzsqap,pnacf,sqap3,p2,palp,val3
      real(kind=8)      :: halfap,osqap3,sum2,expterm1,pxap,pyap,pzap,pxpyap,pypzap,pxpzap,expterm0
      complex(kind=8)  :: expterm,sum1,pe,pezio,expterm2, sum1x,sum1y,sum1z,difexpx,difexpy,difexpz 
      real(kind=8)      :: Cx2,Cy2,Cz2,Cx3,Cx4,Cy3,Cy4,Cz3,Cz4        
      real(kind=8)      :: recp
      complex(kind=8)  :: difexpx2,difexpy2,difexpz2,pehfap ,pe2,zexpterm,zexpterm2
      complex(kind=8)  :: difexpxy,difexpyz,difexpxz,sumxx,sumyy,sumzz,sumxy,sumyz,sumxz
      real(kind=8), allocatable :: sum4(:) 
      end module pdensitygrid

     module polargrid
      implicit none
      integer(kind=4)           :: nt,np,nr,nrtp
      Real(kind=8)              :: ar, rmin,tmin,pmin,rmax,pmax,tmax,rmax2
      Real(kind=8),allocatable  :: pdenp(:),pdent(:),pdenphi(:)   !wp(:),wt(:),wr(:)
      Real(kind=8)    :: dr,dt,dp,dummy,fin,start,c,prod,rr,tt,pp,p,cp,ct,sp,st,c2p,s2p,s2t,c2t
      Real(kind=8) :: ctct,cpcp,stst,spsp,ctst,ctsp,spct,ctcp,cpct,stsp,spst,stcp,cpst,stct,cpsp,spcp,pstst
      integer(kind=4)           :: ir,it,ip,ipts,ipp
      Real(kind=8) ::  d1p,d1phi,d1t,pxo,pyo,pzo,pxo2,pyo2,pzo2,pxo3,pyo3,pzo3,po,po2,pxopyo, pyopzo, pxopzo, pxopyopzo
      end module polargrid

!------------------------------------------------------------------------------------
      module newrap
      integer(kind=4)           :: maxit
      Real(kind=8)              :: eps, error,fx, fy, fz, fdx, fdy, fdz, maxm,minm
      Real(kind=8),allocatable  ::  gpx(:), gpy(:), gpz(:) 
      end module newrap

      module secder
      real(kind=8)  :: d1x,d1y,d1z,dxx,dyy,dzz,dxy,dyz,dxz,val,val1 
      end module secder   


      module para
      integer(kind=4)   :: pjj,blck,fin,ini,moini,mofin,monum,moind
      integer(kind=4)   :: ierr
      integer(kind=4)   :: rank
      integer(kind=4)   :: nprocs
      !integer(kind=4),dimension(MPI_STATUS_SIZE) :: status1
      integer(kind=4)   :: namesize
      integer(kind=4)   :: pn,pn1,st,ed,ii3,jj3,ia3
      character(len=20) :: file1,fln,file2, file3, file4, file5,file6,fname2
      end module para
      
