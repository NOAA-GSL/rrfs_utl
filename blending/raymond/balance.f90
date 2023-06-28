      subroutine balance(ml,nl,kk1, &
                       qw,qvapor,PH1,PHB,t,tw, &
                       MU1,psw,MUB_2d,HGT_2d, &
                       DNW_1d,ZNU_1d,RDNW_1d,RDN_1d,ptop, &
                       PH,MU_2d)

!=======================================================================
! Purpose:  balance WRF initial state
!
! The hydrostatic dry surface pressure (MU) and geopotential
!
! height (PH) are predicted variables in the WRF model,
!
! and should be diagnosed from T, QVAPOR, and PSFC (Hsiao et al. 2010)
!=======================================================================

  implicit none

  integer :: i,j,k,ml,nl,kk1
  real*8    :: sdmd,slmd,ppb,ttb,albn,aln,qvf1,qv1,qv1_b,qv2,qv2_b
  real*8    :: t_full,rho_full,PH_full,PH_hd,ptop


  real*8, dimension(ml,nl,kk1-1) :: qw,qvapor,t,tw
  real*8, dimension(ml,nl,kk1)   :: PH1,PHB
  real*8, dimension(ml,nl) :: MU1,psw,MUB_2d,HGT_2d
  real*8, dimension(kk1-1) :: DNW_1d,ZNU_1d,RDNW_1d,RDN_1d
  real*8, dimension(ml,nl,kk1-1) :: TP,Brho,PP,PPA
  real*8, dimension(ml,nl) :: MU2
  !real*8, intent(inout), dimension(ml,nl,kk1)   :: PH
  !real*8, intent(inout), dimension(ml,nl) :: MU_2d
  real*8, dimension(ml,nl,kk1)   :: PH
  real*8, dimension(ml,nl) :: MU_2d
!F2PY intent(in,out,copy) PH,MU_2d


  real*8, parameter    :: tis0=290.     ! Value used in WRF.
  real*8, parameter    :: ts0=300.      ! Value used in WRF.
  real*8, parameter    :: ps0=100000.   ! Value used in WRF.
!!  real*8, parameter    :: ptop=3000.    ! Value used in WRF.
  real*8, parameter    :: tlp=50.       ! Value used in WRF.
  real*8, parameter    :: gas_constant = 287.0      ! Value used in WRF.
  real*8, parameter    :: gas_constant_v = 461.6    ! Value used in WRF.
  real*8, parameter    :: cp = 7.0*gas_constant/2.0 ! Value used in WRF.
  real*8, parameter    :: kappa = gas_constant / cp
  real*8, parameter    :: cvpm = -(1.-gas_constant / cp)
  real*8, parameter    :: cpovcv = cp/(cp-gas_constant)
  real*8, parameter    :: rd_over_rv = gas_constant / gas_constant_v
  real*8, parameter    :: gravity = 9.81        ! m/s - value used in WRF.

!c=======================================================================
!cds  check MU & PH not change when increment = 0
!c=======================================================================
!psw(:,:)=0
!tw(:,:,:)=0
!qw(:,:,:)=0
!c=======================================================================
!cds  compute increments of dry column air mass (MU)
!c=======================================================================
   do j=1,nl
   do i=1,ml
     sdmd=0.0
     slmd=0.0
      do k=1,kk1-1
         sdmd=sdmd+qw(i,j,k)*DNW_1d(k)
         slmd=slmd+(1.0+qvapor(i,j,k))*DNW_1d(k)
      enddo
     MU2(i,j)=-(psw(i,j)+(MU1(i,j)+MUB_2d(i,j))*sdmd)/slmd
   enddo
   enddo
!c=======================================================================
!cds  compute increments of geopotential (PH)
!c=======================================================================
!compute P before relocation
   do j=1,nl
   do i=1,ml
   do k=1,kk1-1
     ppb=ZNU_1d(k)*MUB_2d(i,j)+ptop
     ttb=(tis0+tlp*log(ppb/ps0))*(ps0/ppb)**kappa
     albn=(gas_constant/ps0)*ttb*(ppb/ps0)**cvpm
     qvf1=1.+qvapor(i,j,k)/rd_over_rv
     aln=-1./(MUB_2d(i,j)+MU1(i,j))*(albn*MU1(i,j) &
          +RDNW_1d(k)*(PH1(i,j,k+1)-PH1(i,j,k)))
     TP(i,j,k)=ps0*((gas_constant*(ts0+t(i,j,k))*qvf1)/ &
          (ps0*(aln+albn)))**cpovcv
     Brho(i,j,k)=1.0/(albn+aln)
     PP(i,j,k)=TP(i,j,k)-ppb
   enddo
   enddo
   enddo
!compute P perturbation's increment (after relocation)
   do j=1,nl
   do i=1,ml
   k=kk1-1
     qv1=0.5*(qw(i,j,k)+qw(i,j,k))
     qv1_b=0.5*(qvapor(i,j,k)+qvapor(i,j,k))
     qv2=-qv1/((1.+qv1_b)*(1.+qv1_b))
     qv2_b=1./(1.+qv1_b)
     qv1=qv1*qv2_b+qv1_b*qv2
     qv1_b=qv1_b*qv2_b
     PPA(i,j,k)=(-0.5/RDNW_1d(k))* &
                ((MU2(i,j)+qv1*MUB_2d(i,j))/qv2_b &
                -(MU1(i,j)+qv1_b*MUB_2d(i,j))*qv2/(qv2_b*qv2_b))
   do k=kk1-2,1,-1
     qv1=0.5*(qw(i,j,k)+qw(i,j,k+1))
     qv1_b=0.5*(qvapor(i,j,k)+qvapor(i,j,k+1))
     qv2=-qv1/((1.+qv1_b)*(1.+qv1_b))
     qv2_b=1./(1.+qv1_b)
     qv1=qv1*qv2_b+qv1_b*qv2
     qv1_b=qv1_b*qv2_b
     PPA(i,j,k)=PPA(i,j,k+1) &
                -(1./RDN_1d(k+1))* &
                ((MU2(i,j)+qv1*MUB_2d(i,j))/qv2_b &
                -(MU1(i,j)+qv1_b*MUB_2d(i,j))*qv2/(qv2_b*qv2_b))
   enddo
   enddo
   enddo
!compute density, PH_full & background hydrostatic PH(PH_hd) after relocation
!and then obtain PH
   do j=1,nl
   do i=1,ml
     PH_full=HGT_2d(i,j)*gravity
     PH_hd=HGT_2d(i,j)*gravity
   do k=1,kk1-1
     t_full=(ts0+t(i,j,k)+tw(i,j,k))*((TP(i,j,k)+PPA(i,j,k))/ps0)**kappa  !theta to temp.
     rho_full=(TP(i,j,k)+PPA(i,j,k))/ &
              (gas_constant*t_full*(1.0+(qw(i,j,k)+qvapor(i,j,k)) &
              /rd_over_rv))
     PH_full=PH_full-DNW_1d(k)* &
             (MU1(i,j)+MUB_2d(i,j)+MU2(i,j))/rho_full
     PH_hd=PH_hd-DNW_1d(k)* &
             (MU1(i,j)+MUB_2d(i,j))/Brho(i,j,k)
     PH(i,j,k+1)=PH_full-PHB(i,j,k+1) &
                 +(PH1(i,j,k+1)+PHB(i,j,k+1)-PH_hd)
   enddo
     PH(i,j,1)=0.
   enddo
   enddo

   do j=1,nl
   do i=1,ml
     MU_2d(i,j)=MU2(i,j)+MU1(i,j)
   enddo
   enddo

      return
      end

