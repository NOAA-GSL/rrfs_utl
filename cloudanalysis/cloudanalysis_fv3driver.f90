program cloudanalysis
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  cloudanalysis      driver for generalized No Varational cloud/hydrometeor analysis
!
!   PRGMMR: Ming Hu          ORG: GSL/AVID        DATE: 2020-08-10
!
! ABSTRACT: 
!  This subroutine serves as a driver for generalized No Varational cloud/hydrometeor analysis
!
! PROGRAM HISTORY LOG:
!    2008-12-20  Hu  Add NCO document block
!
!   input argument list:
!     mype     - processor ID that does this IO
!
!   output argument list:
!
! USAGE:
!   INPUT FILES: 
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 
!   MACHINE:  Linux cluster (WJET) at NOAA/ESRL - Boulder, CO
!
!$$$
!
!_____________________________________________________________________
!
! 
  use mpi
  use kinds,   only: r_single,i_kind, r_kind
!  use wrf_mass_guess_mod, only: soil_temp_cld,isli_cld,ges_xlon,ges_xlat,ges_tten

  use constants, only: init_constants,init_constants_derived
  use constants, only: rd_over_cp, h1000
  use constants, only: zero,one,rad2deg,fv

  use rapidrefresh_cldsurf_mod, only: init_rapidrefresh_cldsurf
  use rapidrefresh_cldsurf_mod, only: dfi_radar_latent_heat_time_period,   &
                                      metar_impact_radius,                 &
                                      l_cleanSnow_WarmTs,l_conserve_thetaV,&
                                      r_cleanSnow_WarmTs_threshold,        &
                                      i_conserve_thetaV_iternum,           &
                                      l_cld_bld, cld_bld_hgt,              &
                                      build_cloud_frac_p, clear_cloud_frac_p, &
                                      nesdis_npts_rad, &
                                      iclean_hydro_withRef, iclean_hydro_withRef_allcol, &
                                      l_use_hydroretrieval_all, &
                                      i_lightpcp, l_numconc, qv_max_inc,ioption, &
                                      l_precip_clear_only,l_fog_off,cld_bld_coverage,cld_clr_coverage,&
                                      i_T_Q_adjust,l_saturate_bkCloud,i_precip_vertical_check,l_rtma3d, &
                                      l_qnr_from_qr, n0_rain, &
                                      r_cloudfrac_threshold

  use namelist_mod, only: load_namelist
  use namelist_mod, only: iyear,imonth,iday,ihour,iminute,isecond
  use namelist_mod, only: fv3_io_layout_y
  use namelist_mod, only: fv3sar_bg_opt

  use get_fv3sar_bk_parall_mod, only: lon2,lat2,nsig
  use get_fv3sar_bk_parall_mod, only: mype_istart,mype_jstart
  use get_fv3sar_bk_parall_mod, only: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  use get_fv3sar_bk_parall_mod, only: read_fv3sar_bk_full
  use get_fv3sar_bk_parall_mod, only: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg, &
                               ges_qnr,ges_qni,ges_qnc,ges_qcf
  use get_fv3sar_bk_parall_mod, only: xlon,xlat,xland,soiltbk
  use get_fv3sar_bk_parall_mod, only: release_mem_fv3sar
  use get_fv3sar_bk_parall_mod, only: update_fv3sar
  use get_fv3sar_bk_parall_mod, only: fv3_io_layout_end 
  use get_fv3sar_bk_parall_mod, only: read_fv3sar_init 
  use gsi_rfv3io_tten_mod, only: nlon_regional,nlat_regional,nsig_regional
!
  implicit none

! MPI variables
  integer :: npe, mype, ierror

! Declare passed variables
!
  integer :: regional_time(6)
! background
!
!  real(r_single),allocatable:: xlon(:,:)        ! 2D longitude in each grid
!  real(r_single),allocatable:: xlat(:,:)        ! 2D latitude in each grid
!  real(r_single),  allocatable:: xland(:,:)
!  real(r_single),allocatable:: soiltbk(:,:)
!
!  surface observation
!
  character(len=7) :: obstype
  character(len=20):: isis
  integer :: nreal,nchanl,ilat,ilon,ndatafv3
  integer :: istart,jstart
!
  integer(i_kind) :: nvarcld_p
  parameter (nvarcld_p=13)

  integer(i_kind)              :: numsao
  real(r_single), allocatable  :: oi(:)
  real(r_single), allocatable  :: oj(:)
  integer(i_kind),allocatable  :: ocld(:,:)
  character*10,   allocatable  :: owx(:)
  real(r_single), allocatable  :: oelvtn(:)
  real(r_single), allocatable  :: odist(:)
  character(8),   allocatable  :: cstation(:)
  real(r_single), allocatable  :: oistation(:)
  real(r_single), allocatable  :: ojstation(:)
  real(r_single), allocatable  :: wimaxstation(:)
!
  integer(i_kind),allocatable  :: osfc_station_map(:,:)
!
!  lightning observation: 2D field in RR grid
!
  real(r_single),allocatable  :: lightning(:,:)
!
!  GOES - NASA LaRc cloud products: several 2D fields in RR grid
!
  real(r_single),allocatable  :: nasalarc_cld(:,:,:)

!
!  radar observation : 3D reflectvity in RR grid
!
  real(r_single),allocatable :: ref_mos_3d(:,:,:)
  real(r_single),allocatable :: ref_mos_3d_tten(:,:,:)
  real(r_single),allocatable :: ref_mosaic31(:,:,:)
  integer(i_kind)          :: nmsclvl_radar 
!
!  GOES - NESDIS cloud products : 2d fields
!
  real(r_single), allocatable :: sat_ctp(:,:)
  real(r_single), allocatable :: sat_tem(:,:)
  real(r_single), allocatable :: w_frac(:,:)
  integer(i_kind),allocatable :: nlev_cld(:,:)
!
! cloud/hydrometeor analysis variables
!
!=========================================================================
!  cld_cover_3d in the Generalized Cloud/Hydrometeor Analysis code
!   Definition:  3-d gridded observation-based information
!      including 0-1 cloud-fraction (with zero value for clear)
!      and negative values indicating "unknown" status.
!   cld_cover_3d is initialized with negative values.
!   cld_type_3d, pcp_type_3d, wthr_type_2d - similar to cld_cover_3d
!=========================================================================

  real(r_single), allocatable :: cld_cover_3d(:,:,:)  ! cloud cover
  integer(i_kind),allocatable :: cld_type_3d(:,:,:)   ! cloud type
  integer(i_kind),allocatable :: pcp_type_3d(:,:,:)   ! precipitation type
  integer(i_kind),allocatable :: wthr_type_2d(:,:)    ! weather type type
  integer(i_kind),allocatable :: cloudlayers_i(:,:,:) ! 5 different layers 
                                                      ! 1= the number of layers
                                                      ! 2,4,... bottom
                                                      ! 3,5,... top
!
  real(r_single),allocatable :: cldwater_3d(:,:,:)    ! cloud water
  real(r_single),allocatable :: nwater_3d(:,:,:)      ! cloud water number concentration
  real(r_single),allocatable :: cldice_3d(:,:,:)      ! cloud ice
  real(r_single),allocatable :: nice_3d(:,:,:)        ! cloud ice number concentration
  real(r_single),allocatable :: rain_3d(:,:,:)        ! rain
  real(r_single),allocatable :: nrain_3d(:,:,:)       ! rain number concentration
  real(r_single),allocatable :: snow_3d(:,:,:)        ! snow
  real(r_single),allocatable :: graupel_3d(:,:,:)     ! graupel
  real(r_single),allocatable :: cldtmp_3d(:,:,:)      ! cloud temperature

  real(r_single),allocatable :: rain_1d_save(:)       ! rain
  real(r_single),allocatable :: nrain_1d_save(:)      ! rain number concentration    
  real(r_single),allocatable :: snow_1d_save(:)       ! snow
  real(r_single),allocatable :: vis2qc(:,:)           ! fog

  real(r_kind)    ::  thunderRadius=2.5_r_kind
  integer(i_kind) :: miss_obs_int
  real(r_kind)    :: miss_obs_real
  parameter ( miss_obs_int = -99999  )
  parameter ( miss_obs_real = -99999.0_r_kind )
  real(r_single)  ::  krad_bot          ! radar bottom level

!
! collect cloud
  real(r_kind)    :: cloud_def_p
  data  cloud_def_p       / 0.000001_r_kind/
  real(r_kind),allocatable :: sumqci(:,:,:)  ! total liquid water
  real(r_kind),allocatable :: watericemax(:,:)  ! max of total liquid water
  integer(i_kind),allocatable :: kwatericemax(:,:)  ! lowest level of total liquid water
  real(r_single),allocatable::temp1(:,:),tempa(:)
  real(r_single),allocatable::all_loc(:,:)
  real(r_single),allocatable::strp(:)
  integer(i_kind) :: im,jm
!
! option in namelist
!
  integer(i_kind) :: opt_cloudwaterice_retri  ! method for cloud water retrieval
  integer(i_kind) :: opt_hydrometeor_retri    ! method for precipitation retrieval
  integer(i_kind) :: opt_cloudtemperature     ! if open temperature adjustment scheme
  integer(i_kind) :: istat_surface,istat_nesdis,istat_radar    ! 1 has observation
  integer(i_kind) :: istat_nasalarc,istat_lightning            ! 0 no observation
  integer(i_kind) :: imerge_nesdis_nasalarc  !  =1 merge NASA LaRC with NESDIS
                                             !  =2 use NASA LaRC only
                                             !  = other, use NESDIS only
!
!
  real(r_kind), pointer :: ges_z (:,:  )=>NULL()  ! geopotential height
  real(r_kind), pointer :: ges_ps(:,:  )=>NULL()  ! surface pressure
  real(r_single), pointer :: ges_tv(:,:,:)=>NULL()  ! virtual temperature
  real(r_single), pointer :: ges_q (:,:,:)=>NULL()  ! specifici humidity
!
!  misc.
!
  integer(i_kind) :: ytotal,ybegin,yend
  integer(i_kind) :: i,j,k,kk
  integer(i_kind) :: iglobal,jglobal,ilocal,jlocal
  logical :: ifindomain
  integer(i_kind) :: imaxlvl_ref
  real(r_kind)    :: max_retrieved_qrqs,max_bk_qrqs,ratio_hyd_bk2obs
  real(r_kind)    :: qrqs_retrieved
  real(r_kind)    :: qrlimit,qrlimit_lightpcp
  real(r_kind)    :: qnr_limit
  real(r_kind)    :: dbz_clean_graupel
  integer(i_kind) :: ilat1s,ilon1s
  integer(i_kind) :: clean_count,build_count,part_count,miss_count
  integer :: sss,rrr

  real(r_kind)    :: refmax,snowtemp,raintemp,nraintemp,graupeltemp
  real(r_kind)    :: snowadd,ratio2
  integer(i_kind) :: imax, jmax, ista, iob, job
  real(r_kind)    :: dfi_lhtp, qmixr, tsfc
  real(r_kind)    :: Temp, watwgt
  real(r_kind)    :: cloudwater, cloudice

  real(r_kind),parameter    :: pi = 4._r_kind*atan(1._r_kind)
  real(r_kind),parameter    :: rho_w = 999.97_r_kind, rho_a = 1.2_r_kind
  real(r_kind),parameter    :: cldDiameter = 10.0E3_r_kind

  real(r_kind),parameter :: am_r = pi * 1000.0_r_kind / 6.0_r_kind
  real(r_kind)           :: lambda

! local variables used for adjustment of qr/qs for RTMA_3D to alleviate ghost reflectivity
  logical         :: print_verbose
  logical         :: verbose
  integer(i_kind) :: k_cap            ! highest level when adjument is done (used for adjust qr/qs for RTMA_3D)
  logical         :: fileexist
  character(len=80) :: obsfile
  integer         :: lunin
  integer         :: nsat1
  integer         :: istatus
!
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  write(obsfile,'(a,I4.4)') 'stdout_cloudanalysis.d',mype
  open(6, file=trim(obsfile),form='formatted',status='unknown')
  write(6,*) '===> cloud analysis over subdomain = ', mype
!
  clean_count=0
  build_count=0
  part_count=0
  miss_count=0
!!
  write(6,*) '========================================'
  write(6,*) 'gsdcloudanalysis: Start generalized cloud analysis', mype
  write(6,*) '========================================'
!
!
  call init_constants(.true.)
  call init_constants_derived

  call init_rapidrefresh_cldsurf
  call load_namelist(mype)
!
! get background ready
!
  call read_fv3sar_init(fv3sar_bg_opt,mype,npe)
  call mpi_barrier(MPI_COMM_WORLD,ierror)

  call read_fv3sar_bk_full(mype)
!
!
  krad_bot=7.0_r_single

  opt_hydrometeor_retri=3       ! 1=Kessler 2=Lin 3=Thompson
  opt_cloudtemperature=3        ! 3=latent heat, 4,5,6 = adiabat profile
  opt_cloudwaterice_retri=1     ! 1 = RUC layer saturation and autoconvert
                                ! 2 = convective 
  imerge_nesdis_nasalarc=2      !  =1 merge NASA LaRC with NESDIS
                                !  =2 use NASA LaRC only
                                !  =3 No Satellite cloud top used
                                !  = other, use NESDIS only

!
! initialize the observation flag  
!
  istat_surface=0
  istat_nesdis=0
  istat_radar=0
  istat_lightning=0
  istat_nasalarc=0

  print_verbose=.false.
  if (verbose) print_verbose=.true.

  write(6,*) "analysis time is=",iyear,imonth,iday,ihour,iminute,isecond
  regional_time(1)=iyear
  regional_time(2)=imonth
  regional_time(3)=iday
  regional_time(4)=ihour
  regional_time(5)=iminute
  regional_time(6)=isecond
  
  write(6,*) 'fv3_io_layout_y=',fv3_io_layout_y
  write(6,*) 'fv3_io_layout_end=',fv3_io_layout_end
!  call load_gsdpbl_hgt(mype)
!
!  check consistency of the options
!

! Now either stratiform or cumulus cloud is considered in the cloud
!  water calculation. This leads to a limitation for the temperature
!  adjustment when stratiform cloud is chosen because adiabat profile
!  scheme based on the convection. This limitation may change when 
!  stratiform and cumulus cloud are both considered at the same time in the future.

  if(opt_cloudwaterice_retri == 1 .and. opt_cloudtemperature >= 4) then
     write(6,*) 'gsdcloudanalysis: ',&
       'inconsistent option for opt_cloudwaterice_retri and opt_cloudtemperature'
     write(6,*) 'gsdcloudanalysis: ',&
       'opt_cloudtemperature must be set to 3 when opt_cloudwaterice_retri =1'
     call stop2(113)
  endif
!!
!!----------------------------------------------
!! 2. read observations                  
!!----------------------------------------------
!!
!! 1.1   allocate observation fields
!!
!
  ybegin=1
  yend=nlat_regional
  ytotal=nlat_regional
  if(fv3_io_layout_y > 1) then
     ybegin=1
     if(mype>0) ybegin=fv3_io_layout_end(mype)+1
     yend=fv3_io_layout_end(mype+1)
     ytotal=fv3_io_layout_end(fv3_io_layout_y)
  endif

  allocate(ref_mos_3d(lon2,lat2,nsig))
  allocate(ref_mos_3d_tten(lon2,lat2,nsig))
  ref_mos_3d=miss_obs_real
  ref_mos_3d_tten=miss_obs_real

  allocate(lightning(lon2,lat2))
  lightning=-9999.0_r_kind

  allocate(sat_ctp(lon2,lat2))
  allocate(sat_tem(lon2,lat2))
  allocate(w_frac(lon2,lat2))
  allocate(nlev_cld(lon2,lat2))
  sat_ctp=miss_obs_real
  sat_tem=miss_obs_real
  w_frac=miss_obs_real
  nlev_cld=miss_obs_int

  allocate(osfc_station_map(lon2,lat2))
  osfc_station_map=miss_obs_int
!
!!
!! 1.2 start to read observations                 
!!
  nmsclvl_radar = -999
  lunin=55

!!  1.2.2 read in surface observations
!!
  istart=mype_istart
  jstart=mype_jstart
  numsao=0
!
  fileexist=.false.
  obsfile='fv3_metarcloud.bin'
  inquire(file=trim(obsfile),exist=fileexist)
  if(fileexist) then
     open(lunin, file=trim(obsfile),form='unformatted')
     read(lunin) obstype,isis,nreal,nchanl,ilat,ilon,ndatafv3
     write(6,*) 'metar cloud=',obstype,isis,nreal,nchanl,ilat,ilon,ndatafv3
     numsao=ndatafv3
     allocate(oi(numsao))
     allocate(oj(numsao))
     allocate(ocld(nvarcld_p,numsao))
     allocate(owx(numsao))
     allocate(oelvtn(numsao))
     allocate(odist(numsao))
     allocate(cstation(numsao))
     allocate(oistation(numsao))
     allocate(ojstation(numsao))
     allocate(wimaxstation(numsao))
     call read_Surface(mype,lunin,istart,jstart,lon2,lat2, &
                       numsao,nvarcld_p,oi,oj,ocld,owx,oelvtn,&
                       odist,cstation,oistation,ojstation)
     if(mype == 0) write(6,*) 'gsdcloudanalysis: ',                                  &
                      'Surface cloud observations are read in successfully'
     istat_surface=1
     close(lunin)
  endif

!!  1.2.6 read in reflectivity mosaic
!!
  fileexist=.false.
  if(fv3_io_layout_y==1) then
    write(obsfile,'(a)') 'RefInGSI3D.dat_01'
  else
    write(obsfile,'(a,I4.4,a)') 'RefInGSI3D.dat.',mype,'_01'
  endif
  write(6,*)
  write(6,*) 'processing ',trim(obsfile)

  inquire(file=trim(obsfile),exist=fileexist)
  if(fileexist) then
     nsat1=0
     nmsclvl_radar=33
     open(lunin, file=trim(obsfile),form='unformatted')
     allocate( ref_mosaic31(lon2,lat2,nmsclvl_radar) )
     ref_mosaic31=-99999.0_r_single

     call read_radar_ref_bin(mype,lunin,istart,jstart,lon2,lat2,nmsclvl_radar,ref_mosaic31)
     write(6,*) 'gsdcloudanalysis: ',                         &
                   ' radar reflectivity is read in successfully'
     istat_radar=1
     close(lunin)
  endif
!
!  1.2.8 read in lightning
!
  fileexist=.false.
  obsfile='LightningInFV3LAM.dat'
  call read_Lightning2cld(obsfile,lon2,lat2,istart,jstart,lightning, &
                          istat_lightning)
  write(6,*) 'gsdcloudanalysis: Lightning is read in successfully'
!
!  1.2.9 read in NASA LaRC cloud products
!
  fileexist=.false.
  obsfile='NASALaRC_cloud4fv3.bin'
  inquire(file=trim(obsfile),exist=fileexist)
  if(fileexist) then
     open(lunin, file=trim(obsfile),form='unformatted')
     allocate(nasalarc_cld(lon2,lat2,5))
     nasalarc_cld=miss_obs_real

     call read_NASALaRC_fv3(mype,lunin,lon2,lat2,istart,jstart,nasalarc_cld)
     write(6,*) 'gsdcloudanalysis:',                       &
                  'NASA LaRC cloud products are read in successfully'
     istat_nasalarc = 1
     close(lunin)
  endif
!
!!
!!  1.4  if there are NASA LaRC cloud products, use them to replace NESDIS ones.
!!       So we use NASA LaRC data in the same way as NESDIS ones
!!
  if(imerge_nesdis_nasalarc == 1 ) then
     if(istat_nasalarc == 1 ) then
        do j=2,lat2-1
           do i=2,lon2-1
             if(sat_ctp(i,j) < -99990.0) then   ! missing value is -999999.0
                sat_ctp(i,j) = nasalarc_cld(i,j,1)
                sat_tem(i,j) = nasalarc_cld(i,j,2)
                w_frac(i,j)  = nasalarc_cld(i,j,3)
                nlev_cld(i,j)= int(nasalarc_cld(i,j,5))
                istat_nesdis =istat_nasalarc
             endif
           enddo
        enddo
     endif
  elseif ( imerge_nesdis_nasalarc == 2) then
     if(istat_nasalarc == 1 ) then
       sat_ctp(:,:) = nasalarc_cld(:,:,1)
       sat_tem(:,:) = nasalarc_cld(:,:,2)
       w_frac(:,:)  = nasalarc_cld(:,:,3)
       nlev_cld(:,:)= int(nasalarc_cld(:,:,5))
       istat_nesdis =istat_nasalarc
     endif
  elseif ( imerge_nesdis_nasalarc == 3) then
       istat_nesdis = 0
  endif
!!
!!
!!  1.6 collect the cloud information from whole domain
!!       and assign the background cloud to each METAR obs
!!
  allocate(sumqci(lon2,lat2,nsig))
  do k=1,nsig
     do j=1,lat2
        do i=1,lon2
           sumqci(i,j,k)= ges_ql(i,j,k) + ges_qi(i,j,k)
        enddo
     enddo
  enddo
!
  allocate(watericemax(lon2,lat2))
  allocate(kwatericemax(lon2,lat2))
  watericemax=0._r_kind
  wimaxstation=0.0_r_single
  kwatericemax=-1
  do j=1,lat2
     do i=1,lon2
       do k = 1,nsig
          watericemax(i,j) = max(watericemax(i,j),sumqci(i,j,k))
       end do
       do k=1,nsig
          if (sumqci(i,j,k) > cloud_def_p .and. kwatericemax(i,j) == -1) then
             kwatericemax(i,j) = k
          end if
       end do
     enddo
  enddo
!!
!  im=nlon_regional
!  jm=nlat_regional
!  allocate(all_loc(lat2,lon2))
!  allocate(strp(lat1*lon1))
!  allocate(tempa(itotsub))
!  allocate(temp1(im,jm))
!  do j=1,lat2
!     do i=1,lon2
!        all_loc(j,i) = watericemax(i,j)
!     enddo
!  enddo
!  call strip(all_loc,strp)
!  call mpi_allgatherv(strp,ijn(mype+1),mpi_real4, &
!            tempa,ijn,displs_g,mpi_real4,mpi_comm_world,ierror)
!  ierror=0
!  if(ierror /= 0 ) then
!     write(*,*) 'MPI error: cloud analysis=',mype
!  endif
!  temp1=0.0_r_single
!  call unfill_mass_grid2t(tempa,im,jm,temp1)
!
!  if(istat_surface==1) then
!     do ista=1,numsao
!        iob = min(max(int(oistation(ista)+0.5),1),im)
!        job = min(max(int(ojstation(ista)+0.5),1),jm)
!        wimaxstation(ista)=temp1(iob,job)
!        if(wimaxstation(ista) > 0._r_single) then
!            i=int(oi(ista))
!            j=int(oj(ista))
!        endif
!     enddo
!  endif
!  deallocate(all_loc,strp,tempa,temp1)
!
!! make a surface station map in grid coordinate
! doesn't work for FV3LAM
!  if(istat_surface==1) then
!     do ista=1,numsao
!        iob = int(oistation(ista)-jstart(mype+1)+2)
!        job = int(ojstation(ista)-istart(mype+1)+2)
!        if(iob >=1 .and. iob<=lon2-1 .and. job >=1 .and. job<=lat2-1) then
!           osfc_station_map(iob,job)=1
!           osfc_station_map(iob+1,job)=1
!           osfc_station_map(iob,job+1)=1
!           osfc_station_map(iob+1,job+1)=1
!        endif
!!     enddo
!  endif

!!
!!  1.8 check if data available: if no data in this subdomain, return. 
!!
  if( (istat_radar + istat_surface + istat_nesdis + istat_lightning ) == 0 ) then
     write(6,*) ' No cloud observations available, skip this domain ', mype
     if(allocated(ref_mos_3d))      deallocate(ref_mos_3d)
     if(allocated(ref_mos_3d_tten)) deallocate(ref_mos_3d_tten)
     if(allocated(lightning))       deallocate(lightning)
     if(allocated(sat_ctp))         deallocate(sat_ctp)
     if(allocated(sat_tem))         deallocate(sat_tem)
     if(allocated(w_frac))          deallocate(w_frac)
     if(allocated(nlev_cld))        deallocate(nlev_cld)
  else
!!
!!----------------------------------------------
!! 2. allocated background arrays and read background  
!!    further observation data process before cloud analysis
!!----------------------------------------------
!
!!
!! 2.2   allocate background and analysis fields
!!
  allocate(cldwater_3d(lon2,lat2,nsig))
  allocate(cldice_3d(lon2,lat2,nsig))
  allocate(rain_3d(lon2,lat2,nsig))
  allocate(nrain_3d(lon2,lat2,nsig))
  allocate(snow_3d(lon2,lat2,nsig))
  allocate(graupel_3d(lon2,lat2,nsig))
  allocate(cldtmp_3d(lon2,lat2,nsig))
  allocate(vis2qc(lon2,lat2))
  cldwater_3d=miss_obs_real
  cldice_3d=miss_obs_real
  rain_3d=miss_obs_real
  nrain_3d=miss_obs_real
  snow_3d=miss_obs_real
  graupel_3d=miss_obs_real
  cldtmp_3d=miss_obs_real
  vis2qc=miss_obs_real
  allocate(rain_1d_save(nsig))
  allocate(nrain_1d_save(nsig))
  allocate(snow_1d_save(nsig))
  rain_1d_save=miss_obs_real
  nrain_1d_save=miss_obs_real
  snow_1d_save=miss_obs_real
  allocate(nice_3d(lon2,lat2,nsig))
  allocate(nwater_3d(lon2,lat2,nsig))
  nice_3d=miss_obs_real
  nwater_3d=miss_obs_real
!!          
!! 2.4 read in background fields
!!          
!   call read_fv3sar_fix
!        zh(i,j)     !  terrain in meter
!        ps_bk(i,j)  !  surace pressure in mb
!        xland(i,j)  !  0=water, 1=land, 2=ice
!        soiltbk(i,j)!  soil temperature
!        xlon(i,j)   !  longitude back to degree
!        xlat(i,j)   !  latitude  back to degree
!        q_bk(i,j,k) ! specific humidity
!        qmixr = q_bk(i,j,k)/(one - q_bk(i,j,k))     ! covert from specific humidity to mixing ratio
!        t_bk(i,j,k)=ges_tv(j,i,k)/                                  &
!                     (one+fv*q_bk(i,j,k))   ! virtual temp to temp
!  call BackgroundCld(mype,lon2,lat2,nsig,t_bk,p_bk,ps_bk,q_bk,h_bk,    &
!             zh,pt_ll,eta1_ll,aeta1_ll,eta2_ll,aeta2_ll,regional,wrf_mass_regional)
!
!! 
!!  2.6 vertical interpolation of radar reflectivity
!!
!  call MPI_Barrier(mpi_comm_world, ierror)
  if(istat_radar ==  1 ) then
     call vinterp_radar_ref(mype,lon2,lat2,nsig,nmsclvl_radar, &
                          ref_mos_3d,ref_mosaic31,h_bk,zh)
     deallocate( ref_mosaic31 )
     ref_mos_3d_tten=ref_mos_3d
     call build_missing_REFcone(mype,lon2,lat2,nsig,krad_bot,ref_mos_3d_tten,h_bk,pblh)
  endif
!!
!!  2.8 convert lightning to reflectivity 
!!  
  if(istat_lightning ==  1 ) then
     call convert_lghtn2ref(mype,lon2,lat2,nsig,ref_mos_3d_tten,lightning,h_bk)
  endif
!!
!!
!!----------------------------------------------
!! 3.  Calculate 3-d cloud cover obs-information field (cld_cover_3d), 
!!               cloud type, precipitation type 
!!----------------------------------------------
!!
  allocate(cld_cover_3d(lon2,lat2,nsig))
  allocate(cld_type_3d(lon2,lat2,nsig))
  allocate(wthr_type_2d(lon2,lat2))
  allocate(pcp_type_3d(lon2,lat2,nsig))
  allocate(cloudlayers_i(lon2,lat2,21))
  cld_cover_3d=miss_obs_real
  cld_type_3d =miss_obs_int
  wthr_type_2d=miss_obs_int
  pcp_type_3d =miss_obs_int
!!
!!
!  call MPI_Barrier(mpi_comm_world, ierror)
  if(istat_surface ==  1) then
     call cloudCover_surface(mype,lat2,lon2,nsig,thunderRadius,           &
              cld_bld_hgt,t_bk,p_bk,q_bk,h_bk,zh,                         &
              numsao,nvarcld_p,numsao,oi,oj,ocld,owx,oelvtn,odist,        &
              cld_cover_3d,cld_type_3d,wthr_type_2d,pcp_type_3d,          &
              wimaxstation, kwatericemax,vis2qc)
     write(6,*) 'gsdcloudanalysis:',                        &  
                   'success in cloud cover analysis using surface data'
  endif
!
  if(istat_nesdis == 1 ) then
     call cloudCover_NESDIS(mype,regional_time,lat2,lon2,nsig,            &
                         xlon,xlat,t_bk,p_bk,h_bk,xland,                  &
                         soiltbk,sat_ctp,sat_tem,w_frac,                  &
                         l_cld_bld,cld_bld_hgt,                           &
                         build_cloud_frac_p,clear_cloud_frac_p,nlev_cld,  &
                         cld_cover_3d,cld_type_3d,wthr_type_2d,osfc_station_map)
     write(6,*) 'gsdcloudanalysis:',                        & 
                   ' success in cloud cover analysis using NESDIS data'
  endif
!  call MPI_Barrier(mpi_comm_world, ierror)
!
!! for Rapid Refresh application, turn off the radar reflectivity impact 
!! on cloud distribution  (Oct. 14, 2010)
!!  if(istat_radar == 1 .or. istat_lightning == 1 ) then
!!     call cloudCover_radar(mype,lat2,lon2,nsig,h_bk,ref_mos_3d,  &
!!                           cld_cover_3d,wthr_type_2d)
!!     if(mype == 0) write(6,*) 'gsdcloudanalysis: ',                 & 
!!                   ' success in cloud cover analysis using radar data'
!!  endif
!
!!
!!----------------------------------------------
!! 4.  Calculate 3-d cloud water and ice  
!!     Calculate 3-d hydrometeors 
!!     Calculate radar temperature tendency
!!     Calculate in cloud temperature
!!     moisture adjustment (to do)
!!----------------------------------------------
!!
!! 4.2 find the cloud layers
!!

  call cloudLayers(lat2,lon2,nsig,h_bk,zh,cld_cover_3d,               &
                   cld_type_3d,cloudlayers_i)
  write(6,*) 'gsdcloudanalysis: success in finding cloud layers'
!!
!! 4.4 decide the cloud type
!!
  call cloudType(lat2,lon2,nsig,h_bk,t_bk,p_bk,ref_mos_3d,            &
                 cld_cover_3d,cld_type_3d,wthr_type_2d,cloudlayers_i)
  write(6,*) 'gsdcloudanalysis: success in deciding cloud types'
!!
!! 4.6 calculate liquid water content
!!
  if(opt_cloudwaterice_retri == 1 ) then
     call cloudLWC_stratiform(mype,lat2,lon2,nsig,q_bk,t_bk,p_bk,      &
                  cld_cover_3d,cld_type_3d,wthr_type_2d,cloudlayers_i, &
                  cldwater_3d,cldice_3d)
     write(6,*) 'gsdcloudanalysis: ',                      &
                 'success in modifying hydrometeors for stratiform clouds '

  elseif (opt_cloudwaterice_retri == 2) then
     call cloudLWC_Cumulus(lat2,lon2,nsig,h_bk,t_bk,p_bk,              &
                  cld_cover_3d,cld_type_3d,wthr_type_2d,cloudlayers_i, &
                  cldwater_3d,cldice_3d,cldtmp_3d)
     write(6,*) 'gsdcloudanalysis: ',                      &
                 ' success in modifying hydrometeors from radar reflectivity'
  else
     write(6,*)'gsdcloudanalysis: ',                                   &
      ' Invalid cloud water calculation option, check opt_cloudwaterice_retri'
     stop 113
  endif
!!
!! 4.8 calculate hydrometeors
!!
!
!  call MPI_Barrier(mpi_comm_world, ierror)
  if(istat_radar == 1 .or. istat_lightning == 1) then
     call PrecipType(lat2,lon2,nsig,t_bk,p_bk,q_bk,ref_mos_3d,         &
                    wthr_type_2d,pcp_type_3d)
     write(6,*) 'gsdcloudanalysis: ',                      &
                 ' success in deciding precipitation type'

!     call PrecipMxR_radar(mype,lat2,lon2,nsig,      &
!                  t_bk,p_bk,ref_mos_3d, &
!                  pcp_type_3d,rain_3d,nrain_3d,snow_3d,graupel_3d,opt_hydrometeor_retri) 
     call hydro_mxr_thompson (lon2,lat2,nsig, t_bk,p_bk,ref_mos_3d, &
                   rain_3d,nrain_3d,snow_3d,istatus, mype )
     write(6,*) 'gsdcloudanalysis: ',                               &
                 ' success in determining hydrometeor types from radar refl'
     if(opt_hydrometeor_retri.ne.3) then
        do k=1,nsig
           do j=1,lat2
              do i=1,lon2
                 nrain_3d(i,j,k)= ges_qnr(i,j,k)
              enddo
           enddo
        enddo

     endif
  endif
!!
!! 4.10 radar temperature tendency for DFI
!!
!  dfi_lhtp=dfi_radar_latent_heat_time_period
!  if (istat_NESDIS /= 1) sat_ctp=miss_obs_real
!  call radar_ref2tten(mype,istat_radar,istat_lightning,lon2,lat2,nsig,ref_mos_3d_tten, &
!                       cld_cover_3d,p_bk,t_bk,ges_tten(:,:,:,1),dfi_lhtp,krad_bot,pblh,sat_ctp)
!
!!
!! 4.12  temperature adjustment
!!
!!  call TempAdjust(mype,lat2,lon2,nsig,opt_cloudtemperature, t_bk, p_bk, w_bk, q_bk, &
!!                   cldwater_3d,cldice_3d,cldtmp_3d)
!!
!!----------------------------------------------
!! 5.  the final analysis or update background
!!----------------------------------------------
!!
!!  the final analysis of cloud 
!!
!
!  call MPI_Barrier(mpi_comm_world, ierror)
  do k=1,nsig
     ! The phy_data file uses vertical levels from bottom to top, opposite of
     ! the other intermediate netcdf files and the expected direction in the
     ! Fortran code. I use the kk index to reverse the direction. CSH
     kk=nsig+1-k
     do j=1,lat2
        do i=1,lon2
           ! clean  cloud
           if( cld_cover_3d(i,j,k) > -0.001_r_kind .and. cld_cover_3d(i,j,k) <= cld_clr_coverage) then 
              cldwater_3d(i,j,k) = zero
              cldice_3d(i,j,k)   = zero
              nice_3d(i,j,k)     = zero
              nwater_3d(i,j,k)   = zero
              clean_count        = clean_count+1
           ! build cloud
           elseif( cld_cover_3d(i,j,k) > cld_bld_coverage .and. cld_cover_3d(i,j,k) < 2.0_r_kind .and. ges_qcf(i,j,kk) < r_cloudfrac_threshold ) then
              cloudwater         =0.001_r_kind*cldwater_3d(i,j,k)
              cloudice           =0.001_r_kind*cldice_3d(i,j,k)
              cldwater_3d(i,j,k) = max(cloudwater,ges_ql(i,j,k))
              cldice_3d(i,j,k)   = max(cloudice,ges_qi(i,j,k))
              ! mhu: Feb2017: set qnc=1e8 and qni=1e6 when build cloud
              if(cloudwater > 1.0e-7_r_kind .and. cloudwater >= ges_ql(i,j,k)) then
                 nwater_3d(i,j,k) = 1.0E8_r_single
              else
                 nwater_3d(i,j,k) = ges_qnc(i,j,k)
              endif
              if(cloudice > 1.0e-7_r_kind .and. cloudice >= ges_qi(i,j,k)) then
                 nice_3d(i,j,k) = 1.0E6_r_single
              else
                 nice_3d(i,j,k) = ges_qni(i,j,k)
              endif
              build_count=build_count+1
           ! unknown or partial cloud, using background values
           else  
              cldwater_3d(i,j,k) = ges_ql(i,j,k)
              cldice_3d(i,j,k)   = ges_qi(i,j,k)
              nice_3d(i,j,k)     = ges_qni(i,j,k)
              nwater_3d(i,j,k)   = ges_qnc(i,j,k)
              if( cld_cover_3d(i,j,k) > cld_clr_coverage ) then
                 part_count=part_count+1
              else
                 miss_count=miss_count+1
              endif
           endif
        end do
     end do
  end do
!
!  the final analysis of precipitation
!
!  move surface Ts check here.  Feb.6 2013
!     l_cleanSnow_WarmTs - if clean snow retrieval when Ts > 5C
!     r_cleanSnow_WarmTs_threshold - threshold for the Ts used in cleaning snow
! If the 1st level temperature is less than 5 degree, then keep 
! snow. Otherwise, use the hybrometeors in the maximum reflectivity level to
!  tune the background hydrometeors.
!
!  NOTE:  l_cleanSnow_WarmTs is alwasy true, it is not an option now. (Feb.4
!  2013)
!

  if(l_use_hydroretrieval_all .or. l_rtma3d) then !RTMA
     qrlimit=15.0_r_kind*0.001_r_kind
     do k=1,nsig
        do j=1,lat2
        do i=1,lon2
           snowtemp=snow_3d(i,j,k)
           raintemp=rain_3d(i,j,k)
           nraintemp=nrain_3d(i,j,k)
           rain_3d(i,j,k) = ges_qr(i,j,k)
           nrain_3d(i,j,k)= ges_qnr(i,j,k)
           snow_3d(i,j,k) = ges_qs(i,j,k)
           graupel_3d(i,j,k) = ges_qg(i,j,k)
           if(ref_mos_3d(i,j,k) > zero ) then
!             snow_3d(i,j,k) = MIN(max(max(snowtemp,zero)*0.001_r_kind,ges_qs(i,j,k)),qrlimit)
              snow_3d(i,j,k) = MIN(    max(snowtemp,zero)*0.001_r_kind               ,qrlimit)
              raintemp = max(raintemp,zero)*0.001_r_kind  
              if(raintemp <= qrlimit) then
                 rain_3d(i,j,k) = raintemp
                 nrain_3d(i,j,k)= nraintemp
              else
                 rain_3d(i,j,k) = qrlimit
                 nrain_3d(i,j,k)= nraintemp*(qrlimit/raintemp)
              endif
           elseif( ref_mos_3d(i,j,k) <= zero .and. & 
                   ref_mos_3d(i,j,k) > -100.0_r_kind ) then
              rain_3d(i,j,k) = zero
              nrain_3d(i,j,k) = zero
              snow_3d(i,j,k) = zero
              graupel_3d(i,j,k) = zero
           else
              rain_3d(i,j,k) = ges_qr(i,j,k)
              nrain_3d(i,j,k)= ges_qnr(i,j,k)
              snow_3d(i,j,k) = ges_qs(i,j,k)
              graupel_3d(i,j,k) = ges_qg(i,j,k)
           endif
        end do
        end do
     end do
 
!    ---- verical check and adjustment to the analysis of precipitation
!         in order to remove/reduce the backround reflectivity "ghost" in
!         analysis.
!         Note: here rain_3d, snow_3d have been already changed into unit of kg/kg.
!     if(l_precip_vertical_check) then
     if(i_precip_vertical_check > 0) then

        if(print_verbose) then
           write(6,'(1x,A,I4.4,A)')"SUB gsdcloudanalysis: precip_vertical_check start... (for pe=",mype,")."
        else
           if(mype == 0) then
              write(6,'(1x,A,I4.4,A)')"SUB gsdcloudanalysis: precip_vertical_check start ... (only print for pe=",mype,")."
           end if
        end if

        qnr_limit=200000_r_kind
        dbz_clean_graupel=35.0_r_kind

        do j=1,lat2
           do i=1,lon2

!          1. search the max reflectivity obs for veritcal profile at each grid
!             point (same code used in hydrometeor anlysis for RAP forecast)
              refmax=-999.0_r_kind
              imaxlvl_ref=0
              do k=1,nsig
                 if(ref_mos_3d(i,j,k) > refmax) then
                    imaxlvl_ref=k
                    refmax=ref_mos_3d(i,j,k)
                 endif
              enddo
!          2. check and adjustment along the profile at each grid point
              if( refmax > 0 .and. (imaxlvl_ref > 0 .and. imaxlvl_ref < nsig ) ) then
                 ! cleaning the Graupel, if refmax <= dbz_clean_graupel (35dbz)
                 ! because graupel is copied from background, not retrieved in cloud analysis.
                 ! (as seen above, graupel_3d(i,j,k) = ges_qg(j,i,k) )
                 if( refmax <= dbz_clean_graupel ) graupel_3d(i,j,:) = zero

                 ! adjusting hydrometeors based on maximum reflectivity level
                 select case (i_precip_vertical_check)
                    case(2)    ! adjust each level along the profile (1:nsig)
                       max_retrieved_qrqs=snow_3d(i,j,imaxlvl_ref)+rain_3d(i,j,imaxlvl_ref)
                       do k=1,nsig
                          qrqs_retrieved=snow_3d(i,j,k)+rain_3d(i,j,k)
                          if(qrqs_retrieved > max_retrieved_qrqs .and. qrqs_retrieved > 0.0001_r_kind &
                             .and. max_retrieved_qrqs > 0.0000001_r_kind) then
                             ratio_hyd_bk2obs=max(min(max_retrieved_qrqs/qrqs_retrieved,1.0_r_kind),0.0_r_kind)
                             if(rain_3d(i,j,k) > zero) then
                                rain_3d(i,j,k) = rain_3d(i,j,k)*ratio_hyd_bk2obs
                                nrain_3d(i,j,k)= min(nrain_3d(i,j,k)/ratio_hyd_bk2obs*2.5_r_kind,qnr_limit)
                             endif
                             if(snow_3d(i,j,k) > zero) then
                                snow_3d(i,j,k) = snow_3d(i,j,k)*ratio_hyd_bk2obs
                             end if
                          end if
                       end do
                    case(3)    ! adjust the dbz-obs-missed levels below max-dbz layer (1:kcap)
                               ! based on the qr+qs on max-refl level
                               ! keep the retrieved cloud analysis as much as possible
                       max_retrieved_qrqs=snow_3d(i,j,imaxlvl_ref)+rain_3d(i,j,imaxlvl_ref)
                       k_cap=min(imaxlvl_ref,nsig)
                       do k=k_cap,1,-1
                          if( ref_mos_3d(i,j,k) <= -100.0_r_kind ) then   !  dbz-obs-missing level
                             qrqs_retrieved=snow_3d(i,j,k)+rain_3d(i,j,k)
                             if(qrqs_retrieved > max_retrieved_qrqs .and. qrqs_retrieved > 0.0001_r_kind &
                                .and. max_retrieved_qrqs > 0.0000001_r_kind) then
                                ratio_hyd_bk2obs=max(min(max_retrieved_qrqs/qrqs_retrieved,1.0_r_kind),0.0_r_kind)
                                if(rain_3d(i,j,k) > zero) then
                                   rain_3d(i,j,k) = rain_3d(i,j,k)*ratio_hyd_bk2obs
                                   nrain_3d(i,j,k)= min(nrain_3d(i,j,k)/ratio_hyd_bk2obs*2.5_r_kind,qnr_limit)    ! 2.5(old) or 1.0(new4) or 1.5(new5/6) 2.5(old, new7) 2.0(new8)
                                endif
                                if(snow_3d(i,j,k) > zero) then
                                   snow_3d(i,j,k) = snow_3d(i,j,k)*ratio_hyd_bk2obs
                                end if
                             end if
                          end if
                       end do
                    case default
                       rain_3d(i,j,k) = rain_3d(i,j,k)
                       nrain_3d(i,j,k)= nrain_3d(i,j,k)
                       snow_3d(i,j,k) = snow_3d(i,j,k)
                 end select

              end if

           end do
        end do

        if(print_verbose) then
           write(6,'(1x,A,I4.4,A)')"SUB gsdcloudanalysis: precip_vertical_check is done ... (for pe=",mype,")."
        else
           if(mype == 0) then
              write(6,'(1x,A,I4.4,A)')"SUB gsdcloudanalysis: precip_vertical_check is done ... (only print for pe=",mype,")."
           end if
        end if

     end if 

  elseif(l_precip_clear_only) then !only clear for HRRRE
     do k=1,nsig
        do j=1,lat2
           do i=1,lon2
              if( ref_mos_3d(i,j,k) <= zero .and. ref_mos_3d(i,j,k) > -100.0_r_kind ) then
                 rain_3d(i,j,k) = zero
                 nrain_3d(i,j,k) = zero
                 snow_3d(i,j,k) = zero
                 graupel_3d(i,j,k) = zero
              else 
                 rain_3d(i,j,k) = ges_qr(i,j,k)
                 nrain_3d(i,j,k)= ges_qnr(i,j,k)
                 snow_3d(i,j,k) = ges_qs(i,j,k)
                 graupel_3d(i,j,k) = ges_qg(i,j,k)
              endif
           enddo
        enddo
     enddo
  else  ! hydrometeor anlysis for RAP forecast
     qrlimit=3.0_r_kind*0.001_r_kind
     qrlimit_lightpcp=1.0_r_kind*0.001_r_kind
     do j=1,lat2
        do i=1,lon2
           refmax=-999.0_r_kind
           imaxlvl_ref=0
           do k=1,nsig
              if(ref_mos_3d(i,j,k) > refmax) then
                 imaxlvl_ref=k
                 refmax=ref_mos_3d(i,j,k)
              endif
              rain_3d(i,j,k)=max(rain_3d(i,j,k)*0.001_r_kind,zero)
              snow_3d(i,j,k)=max(snow_3d(i,j,k)*0.001_r_kind,zero)
              rain_1d_save(k)=rain_3d(i,j,k)
              snow_1d_save(k)=snow_3d(i,j,k)
              nrain_1d_save(k)=nrain_3d(i,j,k)
!              ges_qnr(i,j,k)=max(ges_qnr(i,j,k),zero)
           enddo
           if( refmax > 0 .and. (imaxlvl_ref > 0 .and. imaxlvl_ref < nsig ) ) then       ! use retrieval hybrometeors
              tsfc=t_bk(i,j,1)*(p_bk(i,j,1)/h1000)**rd_over_cp - 273.15_r_kind
              if(tsfc  < r_cleanSnow_WarmTs_threshold) then    ! add snow on cold sfc   
                 do k=1,nsig
                    snowtemp=snow_3d(i,j,k) 
                    rain_3d(i,j,k) = ges_qr(i,j,k)
                    nrain_3d(i,j,k)= ges_qnr(i,j,k)
                    snow_3d(i,j,k) = ges_qs(i,j,k)
                    graupel_3d(i,j,k) = ges_qg(i,j,k)
                    if(ref_mos_3d(i,j,k) > zero ) then
                       snowtemp = MIN(max(snowtemp,ges_qs(i,j,k)),qrlimit)
                       snowadd = max(snowtemp - snow_3d(i,j,k),zero)
                       snow_3d(i,j,k) = snowtemp
                       raintemp=rain_3d(i,j,k) + graupel_3d(i,j,k)
                       if(raintemp > snowadd ) then
                          if(raintemp > 1.0e-6_r_kind) then
                             ratio2=1.0_r_kind - snowadd/raintemp
                             rain_3d(i,j,k) = rain_3d(i,j,k) * ratio2
                             graupel_3d(i,j,k) = graupel_3d(i,j,k) * ratio2
                          endif
                       else
                          rain_3d(i,j,k) = 0.0_r_kind
                          graupel_3d(i,j,k) = 0.0_r_kind
                       endif
                    endif
                 end do
              else    !  adjust hydrometeors based on maximum reflectivity level
                 max_retrieved_qrqs=snow_3d(i,j,imaxlvl_ref)+rain_3d(i,j,imaxlvl_ref)
                 max_bk_qrqs=-999.0_r_kind
                 do k=1,nsig
                    if(ges_qr(i,j,k)+ges_qs(i,j,k) > max_bk_qrqs) then
                        max_bk_qrqs = ges_qr(i,j,k)+ges_qs(i,j,k)
                    endif
                 enddo
                 if( max_bk_qrqs > max_retrieved_qrqs) then ! tune background hyhro
                    ratio_hyd_bk2obs=max(min(max_retrieved_qrqs/max_bk_qrqs,1.0_r_kind),0.0_r_kind)
                    do k=1,nsig
                       graupel_3d(i,j,k) = ges_qg(i,j,k)
                       rain_3d(i,j,k) = ges_qr(i,j,k)
                       nrain_3d(i,j,k)= ges_qnr(i,j,k)
                       snow_3d(i,j,k) = ges_qs(i,j,k)
                       if(ges_qr(i,j,k) > zero) then
                          rain_3d(i,j,k) = ges_qr(i,j,k)*ratio_hyd_bk2obs
                          nrain_3d(i,j,k)= ges_qnr(i,j,k)*ratio_hyd_bk2obs
                       endif
                       if(ges_qs(i,j,k) > zero) &
                          snow_3d(i,j,k) = ges_qs(i,j,k)*ratio_hyd_bk2obs
                    enddo
                 else      !  use hydro in max refl level
                    do k=1,nsig
                       graupel_3d(i,j,k) = ges_qg(i,j,k)
                       if(k==imaxlvl_ref) then
                          snow_3d(i,j,k) = MIN(snow_3d(i,j,k),qrlimit)
                          rain_3d(i,j,k) = MIN(rain_3d(i,j,k),qrlimit)  ! do we need qrlimit?              
                          nrain_3d(i,j,k) = nrain_3d(i,j,k)
                       else
                          rain_3d(i,j,k) = ges_qr(i,j,k)
                          snow_3d(i,j,k) = ges_qs(i,j,k)
                          nrain_3d(i,j,k) = ges_qnr(i,j,k)
                       endif
                    end do
                 endif
                 if(i_lightpcp == 1) then
! keep light precipitation between 28-15 dBZ
                    do k=1,nsig
                       if(ref_mos_3d(i,j,k) >=15.0_r_single .and. &
                          ref_mos_3d(i,j,k) <=28.0_r_single ) then
                          rain_3d(i,j,k) = max(min(rain_1d_save(k),qrlimit_lightpcp),rain_3d(i,j,k))
                          snow_3d(i,j,k) = max(min(snow_1d_save(k),qrlimit_lightpcp),snow_3d(i,j,k)) 
                          nrain_3d(i,j,k)= max(nrain_1d_save(k),nrain_3d(i,j,k))
                       endif
                    enddo  ! light pcp
                 endif
              endif
           else        ! clean if ref=0 or use background hydrometeors
              do k=1,nsig
                 rain_3d(i,j,k) = ges_qr(i,j,k)
                 nrain_3d(i,j,k)= ges_qnr(i,j,k)
                 snow_3d(i,j,k) = ges_qs(i,j,k)
                 graupel_3d(i,j,k) = ges_qg(i,j,k)
                 if((iclean_hydro_withRef==1)) then
                    if( iclean_hydro_withRef_allcol==1 .and. &
                       (refmax <= zero .and. refmax >= -100_r_kind) .and. &
                       (sat_ctp(i,j) >=1010.0_r_kind .and. sat_ctp(i,j) <1050._r_kind)) then     
                       rain_3d(i,j,k) = zero
                       nrain_3d(i,j,k)= zero
                       snow_3d(i,j,k) = zero
                       graupel_3d(i,j,k) = zero
                    else
                       if((ref_mos_3d(i,j,k) <= zero .and.       &
                           ref_mos_3d(i,j,k) > -100.0_r_kind)) then
                          rain_3d(i,j,k) = zero
                          nrain_3d(i,j,k)= zero
                          snow_3d(i,j,k) = zero
                          graupel_3d(i,j,k) = zero
                       endif
                    endif
                 endif
              end do
           endif
        end do
     end do
  endif

!
! If requested, compute rain number concentration from rain mixing
! ratio by assuming an exponential distribution.  The method is a
! simplified version of the make_RainNumber function in the UFS
! module_mp_thompson_make_number_concentrations.F90.
!
  if (l_qnr_from_qr) then
    write(6,*) 'compute rain number concentration from rain mixing ratio'
    write(6,*) 'n0_rain = ', n0_rain
    do k=1,nsig
       do j=1,lat2
          do i=1,lon2
             if (rain_3d(i,j,k) < 0.000000000001_r_kind) then
                rain_3d(i,j,k)  = zero
                nrain_3d(i,j,k) = zero
             else
                lambda = sqrt(sqrt(n0_rain*am_r*6.0_r_kind/rain_3d(i,j,k)))
                nrain_3d(i,j,k) = rain_3d(i,j,k) / 6.0_r_kind &
                                  * lambda*lambda*lambda / am_r
             endif
          end do
       end do
    end do
  endif


!
!  remove any negative hydrometeor mixing ratio or number concentration values
!
  do k=1,nsig
     do j=1,lat2
        do i=1,lon2
           cldwater_3d(i,j,k)= max(0.0_r_single,cldwater_3d(i,j,k))
           cldice_3d(i,j,k)  = max(0.0_r_single,cldice_3d(i,j,k))
           rain_3d(i,j,k)    = max(0.0_r_single,rain_3d(i,j,k))
           nrain_3d(i,j,k)   = max(0.0_r_single,nrain_3d(i,j,k))
           snow_3d(i,j,k)    = max(0.0_r_single,snow_3d(i,j,k))
           graupel_3d(i,j,k) = max(0.0_r_single,graupel_3d(i,j,k))
           nice_3d(i,j,k)    = max(0.0_r_single,nice_3d(i,j,k))
           nwater_3d(i,j,k)  = max(0.0_r_single,nwater_3d(i,j,k))
        end do
     end do
  end do
 
!!
!! move clean process up.   Feb. 6 , 2013
!!  clean the hydrmeteors on grid that:
!!       1)   convective suppress map shows 0 (no convection)
!!       2)   the whole column has no grid whose echo is larger than 0
!!       3)   reflectivity observation show no echo at this grid
!!
!!  if(iclean_hydro_withRef==1) then
!!     do j=2,lat2-1
!!        do i=2,lon2-1
!!           if( abs(ges_tten(j,i,nsig,1)) < 1.0e-5_r_single ) then
!!              inumlvl_ref=0
!!!              do k=1,nsig
!!                if(ref_mos_3d(i,j,k) > zero) then
!!                  inumlvl_ref=inumlvl_ref+1
!!                endif
!!              enddo
!!              if(inumlvl_ref==0) then
!!                 do k=1,nsig
!!                    if(ref_mos_3d(i,j,k) <= zero .and. ref_mos_3d(i,j,k) > -100.0_r_kind ) then
!!                       rain_3d(i,j,k)    = 0.0_r_single
!!!                       snow_3d(i,j,k)    = 0.0_r_single
!!                       graupel_3d(i,j,k) = 0.0_r_single
!!                    endif
!!                 end do
!!              endif
!!           endif
!!        end do
!!!     end do
!!  endif
!!
!!
  if(i_T_Q_adjust==2 .or. i_T_Q_adjust==0) then
     allocate(ges_tv(lon2,lat2,nsig))
     allocate(ges_q(lon2,lat2,nsig))
     ges_tv=t_bk
     ges_q=q_bk
  endif

  call cloud_saturation(mype,l_conserve_thetaV,i_conserve_thetaV_iternum,  &
                 lat2,lon2,nsig,q_bk,t_bk,p_bk,      &
                 cld_cover_3d,wthr_type_2d,cldwater_3d,cldice_3d,sumqci,qv_max_inc) !, l_saturate_bkCloud)

!
!  add fog  (12/08/2015)
!
  if (.not. l_fog_off) then
     do j=2,lat2-1
        do i=2,lon2-1
           if( vis2qc(i,j) > zero ) then
              do k=1,2
                 Temp = t_bk(i,j,k)*(p_bk(i,j,k)/h1000)**rd_over_cp
                 watwgt = max(0._r_kind,min(1._r_kind,(Temp-263.15_r_kind)/&
                                     (268.15_r_kind - 263.15_r_kind)))
                 cldwater_3d(i,j,k) = max(watwgt*vis2qc(i,j),cldwater_3d(i,j,k))
                 cldice_3d(i,j,k)   = max((1.0_r_single-watwgt)*vis2qc(i,j),cldice_3d(i,j,k))
              enddo
           endif
        enddo
     enddo
  endif

!
!  call check_cloud(mype,lat2,lon2,nsig,q_bk,rain_3d,snow_3d,graupel_3d, &
!             cldwater_3d,cldice_3d,t_bk,p_bk,h_bk,                      &
!             numsao,nvarcld_p,numsao,oi,oj,ocld,owx,oelvtn,cstation,    &
!             sat_ctp,cld_cover_3d,xland)
!----------------------------------------------
! 6.  save the analysis results
!----------------------------------------------
!
! for Rapid Refresh application, turn off the hydrometeors 
! (Oct. 14, 2010)
!
! T/Q update
  do k=1,nsig
     do j=1,lat2
        do i=1,lon2
           ! =0 no T/Q adjustment
           if(i_T_Q_adjust==0) then
              write(6,*) 'gsdcloudanalysis: no T/Q adjustment',mype
              t_bk=ges_tv
              q_bk=ges_q
           ! =1 default T/Q adjustment
           elseif(i_T_Q_adjust==1) then
              !t_bk(i,j,k)=t_bk(i,j,k)*(p_bk(i,j,k)/h1000)**rd_over_cp * (one+fv*q_bk(i,j,k))
              t_bk(i,j,k)=t_bk(i,j,k)*(p_bk(i,j,k)/h1000)**rd_over_cp
              ! Here q is mixing ratio kg/kg, need to convert to specific humidity
              q_bk(i,j,k)=q_bk(i,j,k)/(1+q_bk(i,j,k)) 
           ! =2 T/Q adjustment only for case of clearing
           elseif(i_T_Q_adjust==2) then
              !t_bk(i,j,k)=max(ges_tv(i,j,k),t_bk(i,j,k)*(p_bk(i,j,k)/h1000)**rd_over_cp * (one+fv*q_bk(i,j,k)))
              t_bk(i,j,k)=max(ges_tv(i,j,k),t_bk(i,j,k)*(p_bk(i,j,k)/h1000)**rd_over_cp)
              ! Here q is mixing ratio kg/kg, need to convert to specific humidity
              q_bk(i,j,k)=min(ges_q(i,j,k),q_bk(i,j,k)/(1+q_bk(i,j,k)))
           else
              write(6,*) 'gsdcloudanalysis: WARNING no T/Q adjustment, check i_T_Q_adjust value',mype
           endif
        enddo 
     enddo
  enddo
  if(i_T_Q_adjust==2 .or. i_T_Q_adjust==0) then
     deallocate(ges_tv)
     deallocate(ges_q)
  endif

  do k=1,nsig
     do j=1,lat2
        do i=1,lon2
           ! hydrometeor update
!if(mype==1)  then
!           ges_ql(i,j,k)=ges_ql(i,j,k)+0.001*float(j)
!endif
           ges_qr(i,j,k)=rain_3d(i,j,k)
           ges_qs(i,j,k)=snow_3d(i,j,k)
           ges_qg(i,j,k)=graupel_3d(i,j,k)
           ges_ql(i,j,k)=cldwater_3d(i,j,k)
           ges_qi(i,j,k)=cldice_3d(i,j,k)
           ges_qnr(i,j,k)=nrain_3d(i,j,k)
           ! cloud number concentration update
           if( l_numconc ) then
             ges_qni(i,j,k)=nice_3d(i,j,k)
             ges_qnc(i,j,k)=nwater_3d(i,j,k)
           endif
        enddo 
     enddo
  enddo
!
!  check observations
!
if(1==2) then

  do k=1,nsig
     do j=1,lat2
        do i=1,lon2
           ges_qg(i,j,k)=ref_mos_3d(i,j,k)
           ges_qs(i,j,k)=ref_mos_3d_tten(i,j,k)
           if(k<=33) then
!             ges_qni(i,j,k)=ref_mosaic31(i,j,k)
           endif
           if(k<=5 .and. istat_nasalarc == 1) then
             ges_qr(i,j,k)=nasalarc_cld(i,j,k)
           endif
           if(k==6) then
             ges_qr(i,j,k)=lightning(i,j)
           endif
        enddo
     enddo
  enddo
endif
if(1==2) then
 ges_qi=-99999
 loopstation: DO ista=1,numsao
     i = int(oi(ista)+0.0001_r_kind)
     j = int(oj(ista)+0.0001_r_kind)

     if ( ( i >=1 .and. i <= lon2 ) .and. &
          ( j >=1 .and. j <= lat2 )) then
        do k=1,nvarcld_p
           ges_qi(i,j,k)=ocld(k,ista)
        enddo
        ges_qi(i,j,nvarcld_p+1)=Odist(ista)
     end if
   enddo loopstation

endif
!
!  update background for analysis results
!
  call update_fv3sar(mype)
!!
!!----------------------------------------------
!! 7.  release space
!!----------------------------------------------
!!
  deallocate(cld_cover_3d,cld_type_3d,wthr_type_2d, &
             pcp_type_3d,cloudlayers_i)
  deallocate(sumqci)
  deallocate(cldwater_3d,cldice_3d,rain_3d,nrain_3d,snow_3d,graupel_3d,cldtmp_3d)
  deallocate(nice_3d,nwater_3d)
  deallocate(vis2qc)

  if(istat_surface ==  1 ) then
     deallocate(oi,oj,ocld,owx,oelvtn,odist,cstation,oistation,ojstation,wimaxstation)
     deallocate(watericemax,kwatericemax) 
  endif
  if(istat_nasalarc == 1 ) then
     deallocate(nasalarc_cld)
  endif

  deallocate(sat_ctp,sat_tem,w_frac,nlev_cld)
  deallocate(ref_mos_3d,ref_mos_3d_tten,lightning)

  endif
  call release_mem_fv3sar
!
  write(6,*) "CLDcount", r_cloudfrac_threshold,clean_count,build_count,part_count,miss_count
  write(6,*) '========================================'
  write(6,*) 'gsdcloudanalysis: generalized cloud analysis finished:',mype
  write(6,*) '========================================'
!
!
  call MPI_FINALIZE(ierror)
  stop
!


end program cloudanalysis
