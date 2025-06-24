program lake_surgery_hrrr
!
  use kinds, only: i_kind,r_kind,r_single,i_byte
  use module_ncio, only : ncio
  use module_map_utils, only : map_util
!  use module_surface, only : use_surface 
  use mpi
!
  implicit none
!
  type(ncio)     :: raphrrr,rrfs
  type(map_util) :: map
!  type(use_surface) :: sfc
!
!
!  namelist files
!
  character*180 :: rrfs_lam_source
  character*180 :: rrfs_lam_target
  namelist/setup/ rrfs_lam_source,rrfs_lam_target
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! define grid
  integer :: nx_rrfs,ny_rrfs
  integer :: nx_hrrr,ny_hrrr
  integer :: nx_rrfs4,ny_rrfs4
  integer :: nx_rrfs3,ny_rrfs3
  !
  character(len=5),parameter :: gridid='C3463'
  character(len=5) :: haloid
!
! define map
  real(r_single),allocatable,target :: rlon2d_rrfs(:,:),rlat2d_rrfs(:,:)
  real(r_single),allocatable,target :: rlon2d_hrrr(:,:),rlat2d_hrrr(:,:)
  real(r_single),allocatable :: lakemask_hrrr(:,:)
  real(r_single),allocatable :: lakemask_rrfs(:,:)
  real(r_single),allocatable :: lakemask_rrfshrrr(:,:)
  real(r_single),allocatable :: land_frac_rrfs(:,:)
! lake variables
  real(r_single),allocatable :: lake_depth_hrrr(:,:)
  real(r_single),allocatable :: lake_depth_rrfshrrr(:,:)
! C3463_oro_data.tile7.halo0.nc
! lake_depth
  real(r_single),allocatable :: lake_depth_source(:,:)
  real(r_single),allocatable :: lake_depth_target(:,:)
! slmsk
  real(r_single),allocatable :: slmsk_source(:,:)
  real(r_single),allocatable :: slmsk_target(:,:)
! land_frac
  real(r_single),allocatable :: land_frac_source(:,:)
  real(r_single),allocatable :: land_frac_target(:,:)
! vegetation_type_pct
  integer,parameter :: num_veg_cat=20
  real(r_single),allocatable :: vegetation_type_pct_source(:,:,:)
  real(r_single),allocatable :: vegetation_type_pct_target(:,:,:)
! soil_type_pct
  integer,parameter :: num_soil_cat=16
  real(r_single),allocatable :: soil_type_pct_source(:,:,:)
  real(r_single),allocatable :: soil_type_pct_target(:,:,:)
! C3463.soil_type.tile7.halo0.nc
! soil_type_pct
! soil_type
  real(r_single),allocatable :: soil_type_source(:,:)
  real(r_single),allocatable :: soil_type_target(:,:)
  real(r_single),allocatable :: soil_type_pct_source2(:,:,:)
  real(r_single),allocatable :: soil_type_pct_target2(:,:,:)
! C3463.substrate_temperature.tile7.halo0
! substrate_temperature
  real(r_single),allocatable :: substrate_temperature_source(:,:)
  real(r_single),allocatable :: substrate_temperature_target(:,:)
! C3463.vegetation_greenness.tile7.halo0.nc
! vegetation_greenness
  integer,parameter :: num_time=12
  real(r_single),allocatable :: vegetation_greenness_source(:,:,:)
  real(r_single),allocatable :: vegetation_greenness_target(:,:,:)
! C3463.vegetation_type.tile7.halo0.nc 
! vegetation_type_pct
! vegetation_type
  real(r_single),allocatable :: vegetation_type_source(:,:)
  real(r_single),allocatable :: vegetation_type_target(:,:)
  real(r_single),allocatable :: vegetation_type_pct_source2(:,:,:)
  real(r_single),allocatable :: vegetation_type_pct_target2(:,:,:)
! C3463.snowfree_albedo.tile7.halo0.nc
! visible_black_sky_albedo=1 visible_white_sky_albedo=2
! near_IR_black_sky_albedo=3 near_IR_white_sky_albedo=4
  real(r_single),allocatable :: sky_albedo_source(:,:,:,:)
  real(r_single),allocatable :: sky_albedo_target(:,:,:,:)
! C3463.facsf.tile7.halo0.nc 
! facsf
  real(r_single),allocatable :: facsf_source(:,:)
  real(r_single),allocatable :: facsf_target(:,:)
! C3463.maximum_snow_albedo.tile7.halo0.nc 
! maximum_snow_albedo
  real(r_single),allocatable :: maximum_snow_albedo_source(:,:)
  real(r_single),allocatable :: maximum_snow_albedo_target(:,:)
!
!  temperal arrary
  real(r_single),allocatable :: r4d2(:,:),r4d3(:,:,:), r4d4(:,:,:,:)
! 
  integer,parameter :: ncheck=10
  integer :: icheck(ncheck),jcheck(ncheck)
  data icheck /2162,1984,2715,2712,2714,3004,2585,2591,2508,2780/
  data jcheck / 244, 829,1162,1163,1163,1176,1225,1225,1261,1303/
!  data icheck /2598,2508,1766,1742,1741,1742,1744,1705,2239,2839,2842,2715,2466,2407,2581,2582,2584,2588,2589/
!  data jcheck / 552, 674, 801, 806, 808, 808, 808,1079,1096,1161,1161,1162,1171,1193,1220,1223,1224,1224,1224/
!
  integer :: i,j,k,n
  character*180 :: rrfsfile_source,rrfsfile_target,hrrrfile

  real :: xc,yc
  integer :: ii,jj,iii,jjj,ndist,mindist,nsearch
  integer :: n_lake2land,n_land2lake
  integer :: vegenum, albdonum, facnum
  logical :: l_lake_surgery, l_consist
  logical :: update_halo4,update_halo3

!**********************************************************************
!
!HRRR lake mask is from HRRR initial condition, for example:
!/mnt/lfs5/BMC/wrfruc/mhu/wcoss/nco/com/hrrr/prod/hrrr.20250529/conus/hrrr.t19z.wrf_inout
!LAKEMASK: 0 is land and 1 is water.

!RRFS lake mask is from the fix file:
!/mnt/lfs5/BMC/nrtrr/FIX_RRFS/lam/RRFS_NA_3km_C3463_Lake_fracSV
!file: C3463_oro_data.tile7.halo0.nc
!for lake_frac: any number larger than 0.75 is lake.

!If we change a land point to a lake point in RRFS, we will change lake_frac from 0 to 1.0 .
! (RRFS treats the whole grid point as lake if lake_frac>0.75)
!if we change a lake point to land, we change lake_frac to 0 

!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
  if(mype==0) then
!
     haloid='halo0'
     l_lake_surgery=.true.
     l_consist=.false.
     update_halo4=.true.
     update_halo3=.true.
!
     rrfs_lam_source='missing'
     rrfs_lam_target='missing'
     open(15, file='lake_surgery.namelist')
        read(15,setup)
     close(15)
     write(*,setup)

     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
     hrrrfile='hrrr.wrf_inout'

     write(*,*) 'open rrfs file ',trim(rrfsfile_source)
!
! read in rrfs dimension, latlon, and lake fraction
!
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_dim("lon",nx_rrfs)
     call rrfs%get_dim("lat",ny_rrfs)
     write(*,*) 'nx_rrfs,ny_rrfs=',nx_rrfs,ny_rrfs

     allocate(rlon2d_rrfs(nx_rrfs,ny_rrfs))
     allocate(rlat2d_rrfs(nx_rrfs,ny_rrfs))
     allocate(lakemask_rrfs(nx_rrfs,ny_rrfs))
     allocate(land_frac_rrfs(nx_rrfs,ny_rrfs))
     call rrfs%get_var("geolon",nx_rrfs,ny_rrfs,rlon2d_rrfs)
     call rrfs%get_var("geolat",nx_rrfs,ny_rrfs,rlat2d_rrfs)
     call rrfs%get_var("lake_frac",nx_rrfs,ny_rrfs,lakemask_rrfs)
     call rrfs%get_var("land_frac",nx_rrfs,ny_rrfs,land_frac_rrfs)
     call rrfs%close()
     write(*,*) 'rrfs lon=',maxval(rlon2d_rrfs),minval(rlon2d_rrfs)
     write(*,*) 'rrfs lat=',maxval(rlat2d_rrfs),minval(rlat2d_rrfs)
     write(*,*) 'rrfs lake fraction=',maxval(lakemask_rrfs),minval(lakemask_rrfs)
     write(*,*) '1811,783=',rlon2d_rrfs(1811,783),rlat2d_rrfs(1811,783)
!
! read in HRRR dimension, latlon, and lake mask
!
     call raphrrr%open(trim(hrrrfile),"r",200)
     call raphrrr%get_dim("west_east",nx_hrrr)
     call raphrrr%get_dim("south_north",ny_hrrr)
     write(*,*) 'nx_rap,ny_rap=',nx_hrrr,ny_hrrr

     allocate(rlon2d_hrrr(nx_hrrr,ny_hrrr))
     allocate(rlat2d_hrrr(nx_hrrr,ny_hrrr))
     allocate(lakemask_hrrr(nx_hrrr,ny_hrrr))
     allocate(lake_depth_hrrr(nx_hrrr,ny_hrrr))
     call raphrrr%get_var("XLONG",nx_hrrr,ny_hrrr,rlon2d_hrrr)
     call raphrrr%get_var("XLAT",nx_hrrr,ny_hrrr,rlat2d_hrrr)
     call raphrrr%get_var("LAKEMASK",nx_hrrr,ny_hrrr,lakemask_hrrr)
     call raphrrr%get_var("LAKEDEPTH2D",nx_hrrr,ny_hrrr,lake_depth_hrrr)
     call raphrrr%close()
     write(*,*) 'hrrr lon=',maxval(rlon2d_hrrr),minval(rlon2d_hrrr)
     write(*,*) 'hrrr lat=',maxval(rlat2d_hrrr),minval(rlat2d_hrrr)
     write(*,*) 'hrrr lake fraction=',maxval(lakemask_hrrr),minval(lakemask_hrrr)
     write(*,*) 'hrrr lake depth=',maxval(lake_depth_hrrr),minval(lake_depth_hrrr)

!
! now replace rrfs lake mask with HRRR lake mask
!
     call map%init_general_transform(nx_hrrr,ny_hrrr,rlat2d_hrrr,rlon2d_hrrr)
     deallocate(rlon2d_hrrr)
     deallocate(rlat2d_hrrr)
!
     allocate(lakemask_rrfshrrr(nx_rrfs,ny_rrfs))
     allocate(lake_depth_rrfshrrr(nx_rrfs,ny_rrfs))
     lakemask_rrfshrrr=lakemask_rrfs
     lake_depth_rrfshrrr=10.0
     write(*,*) 'generate the new lake mask'
     call update_lake_mask(map,nx_rrfs,ny_rrfs,nx_hrrr,ny_hrrr,rlon2d_rrfs,rlat2d_rrfs,&
                          lakemask_hrrr,lakemask_rrfshrrr, &
                          lake_depth_hrrr, lake_depth_rrfshrrr)
     call map%destory_general_transform()
!
     if(.not.l_lake_surgery) then
        lakemask_rrfshrrr=lakemask_rrfs
     endif
!     deallocate(rlon2d_rrfs)
!     deallocate(rlat2d_rrfs)
     deallocate(lakemask_hrrr)
     deallocate(lake_depth_hrrr)
!
!  read all lake variables
!
     allocate(lake_depth_source(nx_rrfs,ny_rrfs))
     allocate(lake_depth_target(nx_rrfs,ny_rrfs))
     allocate(slmsk_source(nx_rrfs,ny_rrfs))
     allocate(slmsk_target(nx_rrfs,ny_rrfs))
     allocate(land_frac_source(nx_rrfs,ny_rrfs))
     allocate(land_frac_target(nx_rrfs,ny_rrfs))
     allocate(vegetation_type_pct_source(nx_rrfs,ny_rrfs,num_veg_cat))
     allocate(vegetation_type_pct_target(nx_rrfs,ny_rrfs,num_veg_cat))
     allocate(soil_type_pct_source(nx_rrfs,ny_rrfs,num_soil_cat))
     allocate(soil_type_pct_target(nx_rrfs,ny_rrfs,num_soil_cat))

     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("lake_depth",nx_rrfs,ny_rrfs,lake_depth_source)
     lake_depth_target=lake_depth_source
     call rrfs%get_var("slmsk",nx_rrfs,ny_rrfs,slmsk_source)
     slmsk_target=slmsk_source
     call rrfs%get_var("land_frac",nx_rrfs,ny_rrfs,land_frac_source)
     land_frac_target=land_frac_source
     call rrfs%get_var("vegetation_type_pct",nx_rrfs,ny_rrfs,num_veg_cat,vegetation_type_pct_source)
     vegetation_type_pct_target=vegetation_type_pct_source
     call rrfs%get_var("soil_type_pct",nx_rrfs,ny_rrfs,num_soil_cat,soil_type_pct_source)
     soil_type_pct_target=soil_type_pct_source
     call rrfs%close()

     allocate(soil_type_source(nx_rrfs,ny_rrfs))
     allocate(soil_type_target(nx_rrfs,ny_rrfs))
     allocate(soil_type_pct_source2(nx_rrfs,ny_rrfs,num_soil_cat))
     allocate(soil_type_pct_target2(nx_rrfs,ny_rrfs,num_soil_cat))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.soil_type.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("soil_type",nx_rrfs,ny_rrfs,soil_type_source)
     soil_type_target=soil_type_source
     call rrfs%get_var("soil_type_pct",nx_rrfs,ny_rrfs,num_soil_cat,soil_type_pct_source2)
     soil_type_pct_target2=soil_type_pct_source2
     call rrfs%close()
!
     allocate(substrate_temperature_source(nx_rrfs,ny_rrfs))
     allocate(substrate_temperature_target(nx_rrfs,ny_rrfs))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.substrate_temperature.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("substrate_temperature",nx_rrfs,ny_rrfs,substrate_temperature_source)
     substrate_temperature_target=substrate_temperature_source
     call rrfs%close()
     !
     allocate(vegetation_greenness_source(nx_rrfs,ny_rrfs,num_time))
     allocate(vegetation_greenness_target(nx_rrfs,ny_rrfs,num_time))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.vegetation_greenness.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("vegetation_greenness",nx_rrfs,ny_rrfs,num_time,vegetation_greenness_source)
     vegetation_greenness_target=vegetation_greenness_source
     call rrfs%close()
     !
     allocate(vegetation_type_source(nx_rrfs,ny_rrfs))
     allocate(vegetation_type_target(nx_rrfs,ny_rrfs))
     allocate(vegetation_type_pct_source2(nx_rrfs,ny_rrfs,num_veg_cat))
     allocate(vegetation_type_pct_target2(nx_rrfs,ny_rrfs,num_veg_cat))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.vegetation_type.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("vegetation_type",nx_rrfs,ny_rrfs,vegetation_type_source)
     vegetation_type_target=vegetation_type_source
     call rrfs%get_var("vegetation_type_pct",nx_rrfs,ny_rrfs,num_veg_cat,vegetation_type_pct_source2)
     vegetation_type_pct_target2=vegetation_type_pct_source2
! vegetation_type_pct
     call rrfs%close()
     allocate(sky_albedo_source(nx_rrfs,ny_rrfs,num_time,4))
     allocate(sky_albedo_target(nx_rrfs,ny_rrfs,num_time,4))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.snowfree_albedo.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("visible_black_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_source(:,:,:,1))
     call rrfs%get_var("visible_white_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_source(:,:,:,2))
     call rrfs%get_var("near_IR_black_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_source(:,:,:,3))
     call rrfs%get_var("near_IR_white_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_source(:,:,:,4))
     sky_albedo_target=sky_albedo_source
     call rrfs%close()
     ! facsf_source
     allocate(facsf_source(nx_rrfs,ny_rrfs))
     allocate(facsf_target(nx_rrfs,ny_rrfs))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.facsf.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("facsf",nx_rrfs,ny_rrfs,facsf_source(:,:))
     facsf_target=facsf_source
     call rrfs%close()
     ! maximum_snow_albedo
     allocate(maximum_snow_albedo_source(nx_rrfs,ny_rrfs))
     allocate(maximum_snow_albedo_target(nx_rrfs,ny_rrfs))
     rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.maximum_snow_albedo.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_source),"r",200)
     call rrfs%get_var("maximum_snow_albedo",nx_rrfs,ny_rrfs,maximum_snow_albedo_source(:,:))
     maximum_snow_albedo_target=maximum_snow_albedo_source
     call rrfs%close()
!
!  convert 
!
     n_lake2land=0
     n_land2lake=0
     nsearch=10
     if(l_consist) then
! fix vegetation_greenness,substrate_temperature,sky_albedo
       do j=1,ny_rrfs
       do i=1,nx_rrfs
          if(abs(slmsk_source(i,j)-1.0) < 0.01) then
             vegenum=0
             albdonum=0
             facnum=0
             do k=1,num_time
                 if(vegetation_greenness_source(i,j,k) < -100.0) vegenum=vegenum+1
             enddo
             do k=1,num_time
                 if(sky_albedo_source(i,j,k,1) < -100.0) albdonum=albdonum+1
                 if(sky_albedo_source(i,j,k,2) < -100.0) albdonum=albdonum+1
                 if(sky_albedo_source(i,j,k,3) < -100.0) albdonum=albdonum+1
                 if(sky_albedo_source(i,j,k,4) < -100.0) albdonum=albdonum+1
             enddo
             if( vegenum > 0 .or. albdonum > 0 .or. substrate_temperature_source(i,j) < 100.0 ) then
                  !need to fix missing values
                  ! find the close land point
                  mindist=9999
                  do jj=max(j-nsearch,1),min(j+nsearch,ny_rrfs)
                  do ii=max(i-nsearch,1),min(i+nsearch,nx_rrfs)
                     ndist=(jj-j)*(jj-j)+(ii-i)*(ii-i)
                     if( abs(slmsk_source(ii,jj)-1.0) < 0.01 .and. substrate_temperature_source(ii,jj) > 100.0 ) then
                       vegenum=0
                       albdonum=0
                       do k=1,num_time
                          if(vegetation_greenness_source(ii,jj, k) > -100.0) vegenum=vegenum+1
                       enddo
                       do k=1,num_time
                          if(sky_albedo_source(ii,jj,k,1) > -100.0) albdonum=albdonum+1
                          if(sky_albedo_source(ii,jj,k,2) > -100.0) albdonum=albdonum+1
                          if(sky_albedo_source(ii,jj,k,3) > -100.0) albdonum=albdonum+1
                          if(sky_albedo_source(ii,jj,k,4) > -100.0) albdonum=albdonum+1
                       enddo
                       if( vegenum==num_time .and. albdonum==4*num_time ) then
                          if( mindist > ndist ) then
                             mindist=ndist
                             iii=ii
                             jjj=jj
                          endif
                       endif
                    endif
                  enddo
                  enddo
                  if(mindist<1000) then
                       write(*,'(a,2I5,a,2I5)') "fix missing value ",i,j, " using ",iii,jjj
                       substrate_temperature_target(i,j)=substrate_temperature_source(iii,jjj)
                       vegetation_greenness_target(i,j,:)=vegetation_greenness_source(iii,jjj,:)
                       sky_albedo_target(i,j,:,:)=sky_albedo_source(iii,jjj,:,:)
                  else
                      write(*,*) 'cannot find point to fill in missing point',i,j
                  endif
             endif
          endif
       enddo
       enddo
     !  substrate_temperature_source=substrate_temperature_target
     !  vegetation_greenness_source=vegetation_greenness_target
     !  sky_albedo_source=sky_albedo_target
     endif
       !
     if(l_lake_surgery) then
       write(*,*) 'converting ....'
       do j=1,ny_rrfs
       do i=1,nx_rrfs
         if( lakemask_rrfs(i,j) > 0.5 .and. lakemask_rrfshrrr(i,j) < 0.5 ) then 
           if( ( i>=2585-8 .and. i<=2585+12 .and. j>=1225-8 .and. j <=1225+5) .or. &
               ( i>=2714-5 .and. i<=2714+5  .and. j>=1163-5 .and. j <=1163+3) ) then
! covert lake to water point in great lake
!land_frac=0
!lake_frac = 0
!lake_depth = 0
!slmsk=0
              lake_depth_target(i,j)=0.0
              slmsk_target(i,j)=0.0
              land_frac_target(i,j)=0.0

              vegetation_type_pct_target(i,j,:)=0.0
              vegetation_type_pct_target(i,j,17)=1.0
              soil_type_pct_target(i,j,:)=0.0
              soil_type_pct_target(i,j,14)=1.0

              soil_type_pct_target2(i,j,:)=0.0
              soil_type_pct_target2(i,j,14)=1.0
              soil_type_target(i,j)=14

              substrate_temperature_target(i,j)=-999.0
              vegetation_greenness_target(i,j,:)=-999.0
              vegetation_type_pct_target2(i,j,:)=0.0
              vegetation_type_pct_target2(i,j,17)=1.0
              vegetation_type_target(i,j)=17
              sky_albedo_target(i,j,:,:)=-999.9
              facsf_target(i,j)=-999.0
              maximum_snow_albedo_target(i,j)=-999.0

!  convert lake point to land point
           else
             n_lake2land=n_lake2land+1
           ! find the close land point
             mindist=9999
             do jj=max(j-nsearch,1),min(j+nsearch,ny_rrfs)
             do ii=max(i-nsearch,1),min(i+nsearch,nx_rrfs)
                ndist=(jj-j)*(jj-j)+(ii-i)*(ii-i)
                if( lakemask_rrfs(ii,jj) < 0.5 .and. mindist > ndist ) then
                   if(substrate_temperature_source(ii,jj) > 100.0 ) then
                      vegenum=0
                      albdonum=0
                      facnum=0
                      if(facsf_target(ii,jj) > -100.0) facnum=facnum+1
                      if(maximum_snow_albedo_target(ii,jj) > -100.0) facnum=facnum+1
                      do k=1,num_time
                         if(vegetation_greenness_source(ii,jj, k) > -100.0) vegenum=vegenum+1
                      enddo
                      do k=1,num_time
                         if(sky_albedo_source(ii,jj,k,1) > -100.0) albdonum=albdonum+1
                         if(sky_albedo_source(ii,jj,k,2) > -100.0) albdonum=albdonum+1
                         if(sky_albedo_source(ii,jj,k,3) > -100.0) albdonum=albdonum+1
                         if(sky_albedo_source(ii,jj,k,4) > -100.0) albdonum=albdonum+1
                      enddo
                      if( vegenum==num_time .and. albdonum==4*num_time .and. facnum==2 ) then
                         mindist=ndist
                         iii=ii
                         jjj=jj
                      endif
                   endif
                endif
             enddo
             enddo
             if(abs(slmsk_source(iii,jjj)-1.0) > 0.01) then
                 write(*,*) 'make sure it is land ',iii,jjj,slmsk_source(iii,jjj),lakemask_rrfs(iii,jjj)
                 stop 999
             endif
             if(mindist<1000) then
                write(*,'(a,2I5,f5.1,a,2I5,f5.1)') "change lake ",i,j,lakemask_rrfs(i,j), &
                                                 " to land",iii,jjj,lakemask_rrfs(iii,jjj)

!   1) set lake_frac=0
!   2) lake_depth: 0
!   3) set slmsk=1 and land_frac=1
!   4) vegetation_type_pct: look for nearest neighbor land
!   5) soil_type_pct: nearest neighbor
!   6) soil_type_pct nearest neighbor land
!   7) soil_type dominant: nearest neighbor land
!   8) set substrate_temperature to nearest neighbor land
!   9) set vegetation_greenness= nearest neighbor land
!   10) vegetation_type_pct as item 4
!   11) vegetation_type=nearest neighbor land
!   12) set all 4 to nearest neighbor land
                lake_depth_target(i,j)=0.0
                slmsk_target(i,j)=1.0
                land_frac_target(i,j)=1.0
                vegetation_type_pct_target(i,j,:)=vegetation_type_pct_source(iii,jjj,:)
                soil_type_pct_target(i,j,:)=soil_type_pct_source(iii,jjj,:)

                soil_type_target(i,j)=soil_type_source(iii,jjj)
                soil_type_pct_target2(i,j,:)=soil_type_pct_source2(iii,jjj,:)
  
                substrate_temperature_target(i,j)=substrate_temperature_source(iii,jjj)
                vegetation_greenness_target(i,j,:)=vegetation_greenness_source(iii,jjj,:)
                vegetation_type_target(i,j)=vegetation_type_source(iii,jjj)
                vegetation_type_pct_target2(i,j,:)=vegetation_type_pct_source2(iii,jjj,:)
                sky_albedo_target(i,j,:,:)=sky_albedo_source(iii,jjj,:,:)
                facsf_target(i,j)=facsf_source(iii,jjj)
                maximum_snow_albedo_target(i,j)=maximum_snow_albedo_source(iii,jjj)
             else
                write(*,*) "cannot find nearest land for", i,j
                do jj=max(j-nsearch,1),min(j+nsearch,ny_rrfs)
                  write(*,'(40f4.1)') (lakemask_rrfs(ii,jj),ii=max(i-nsearch,1),min(i+nsearch,nx_rrfs))
                enddo
                stop 7777
             endif
           endif
!  convert land point to lake point
         elseif ( (lakemask_rrfs(i,j) < 0.5 .and. lakemask_rrfshrrr(i,j) > 0.5 ) ) then
           n_land2lake=n_land2lake+1
           
!  1) set lake_frac=1
!  2) lake_depth: use value from HRRR lake depth. But if HRRR lake depth is 50m, then set it to 10m.
!  3) set slmsk=0 and land_frac=0
!  4) vegetation_type_pct: set index 17 (num_veg_cat) to 1, others indexes are 0
!  5) soil_type_pct: set index 14 (num_soil_cat) to 1, other indexes are 0
!  6) soil_type_pct same as item 5.
!  7) soil_type: set it to 14.
!  8) set substrate_temperature to -999.0
!  9) set vegetation_greenness=0
!  10) vegetation_type_pct as item 4
!  11) vegetation_type=17
!  12) set all 4 to -999.9

           lake_depth_target(i,j)=lake_depth_rrfshrrr(i,j)
           if(abs(lake_depth_rrfshrrr(i,j)-50.0) < 0.1) lake_depth_target(i,j)=10.0
           slmsk_target(i,j)=0.0
           land_frac_target(i,j)=0.0
           vegetation_type_pct_target(i,j,:)=0.0
           vegetation_type_pct_target(i,j,17)=1.0
           soil_type_pct_target(i,j,:)=0.0
           soil_type_pct_target(i,j,14)=1.0

           soil_type_pct_target2(i,j,:)=0.0
           soil_type_pct_target2(i,j,14)=1.0
           soil_type_target(i,j)=14

           substrate_temperature_target(i,j)=-999.0
           vegetation_greenness_target(i,j,:)=-999.0
           vegetation_type_pct_target2(i,j,:)=0.0
           vegetation_type_pct_target2(i,j,17)=1.0
           vegetation_type_target(i,j)=17
           sky_albedo_target(i,j,:,:)=-999.9
           facsf_target(i,j)=-999.0
           maximum_snow_albedo_target(i,j)=-999.0
         endif
       enddo
       enddo
     endif
     write(*,*) 'number points changing from lake to land point=',n_lake2land
     write(*,*) 'number points changing from land to lake point=',n_land2lake
!
!  turn Pamlico Sound and Albemarle Sound in NC to lake
!
     if(l_lake_surgery) then
        write(*,*) 'turn Pamlico Sound and Albemarle Sound in NC to lake'
        call change_Outer_Banks(nx_rrfs,ny_rrfs,lakemask_rrfshrrr,land_frac_rrfs,lake_depth_target)
     endif

if(1==1) then
 do n=1,ncheck
     iii=icheck(n)
     jjj=jcheck(n)
     write(*,*)
     write(*,*) 'check ===',iii,jjj
     write(*,*) 'lon lat', rlon2d_rrfs(iii,jjj)-360.0,rlat2d_rrfs(iii,jjj)

     if(iii==2585 .and. jjj==1225) then
     write(*,*) 'lakemask new'
     do j=jjj-8,jjj+5
       write(*,'(25f5.2)') (lakemask_rrfshrrr(i,j),i=iii-8,iii+12)
     enddo
     write(*,*) 'lakemask old'
     do j=jjj-8,jjj+5
       write(*,'(25f5.2)') (lakemask_rrfs(i,j),i=iii-8,iii+12)
     enddo
     endif
     if(iii==2714 .and. jjj==1163) then
     write(*,*) 'lakemask new'
     do j=jjj-5,jjj+3
       write(*,'(25f5.2)') (lakemask_rrfshrrr(i,j),i=iii-5,iii+5)
     enddo
     write(*,*) 'lakemask old'
     do j=jjj-5,jjj+3
       write(*,'(25f5.2)') (lakemask_rrfs(i,j),i=iii-5,iii+5)
     enddo
     endif


     write(*,*) 'lakemask'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (lakemask_rrfshrrr(i,j),i=iii-4,iii+4), &
                            (lakemask_rrfs(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'lake depth'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (lake_depth_target(i,j),i=iii-4,iii+4), &
                            (lake_depth_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'slmsk'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (slmsk_target(i,j),i=iii-4,iii+4), &
                            (slmsk_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'land_frac'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (land_frac_target(i,j),i=iii-4,iii+4), &
                            (land_frac_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'vegetation_type_pct 17'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (vegetation_type_pct_target(i,j,17),i=iii-4,iii+4), &
                            (vegetation_type_pct_source(i,j,17),i=iii-4,iii+4)
     enddo
     write(*,*) 'soil_type_pct_target 14'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (soil_type_pct_target(i,j,14),i=iii-4,iii+4), &
                            (soil_type_pct_source(i,j,14),i=iii-4,iii+4)
     enddo
     write(*,*) 'soil_type'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (soil_type_target(i,j),i=iii-4,iii+4), &
                            (soil_type_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'soil_type_pct 2 - 14'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (soil_type_pct_target2(i,j,14),i=iii-4,iii+4), &
                            (soil_type_pct_source2(i,j,14),i=iii-4,iii+4)
     enddo
     write(*,*) 'substrate_temperature_target'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (substrate_temperature_target(i,j),i=iii-4,iii+4), &
                            (substrate_temperature_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'vegetation_greenness_source'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (vegetation_greenness_target(i,j,6),i=iii-4,iii+4), &
                            (vegetation_greenness_source(i,j,6),i=iii-4,iii+4)
     enddo
     write(*,*) 'vegetation_type_source'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (vegetation_type_target(i,j),i=iii-4,iii+4), &
                            (vegetation_type_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'vegetation_type_pct 2 17'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (vegetation_type_pct_target2(i,j,17),i=iii-4,iii+4), &
                            (vegetation_type_pct_source2(i,j,17),i=iii-4,iii+4)
     enddo
     write(*,*) 'sky_albedo_target'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (sky_albedo_target(i,j,6,1),i=iii-4,iii+4), &
                            (sky_albedo_source(i,j,6,1),i=iii-4,iii+4)
     enddo
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (sky_albedo_target(i,j,6,2),i=iii-4,iii+4), &
                            (sky_albedo_source(i,j,6,2),i=iii-4,iii+4)
     enddo
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (sky_albedo_target(i,j,6,3),i=iii-4,iii+4), &
                            (sky_albedo_source(i,j,6,3),i=iii-4,iii+4)
     enddo
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (sky_albedo_target(i,j,6,4),i=iii-4,iii+4), &
                            (sky_albedo_source(i,j,6,4),i=iii-4,iii+4)
     enddo
     write(*,*) 'facsf_target'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (facsf_target(i,j),i=iii-4,iii+4), &
                            (facsf_source(i,j),i=iii-4,iii+4)
     enddo
     write(*,*) 'maximum_snow_albedo_target'
     do j=jjj-4,jjj+4
       write(*,'(9f8.2,5x,9f8.2)') (maximum_snow_albedo_target(i,j),i=iii-4,iii+4), &
                            (maximum_snow_albedo_source(i,j),i=iii-4,iii+4)
     enddo
  enddo
endif


!
!  release memory
!
     deallocate(lakemask_rrfs)
     deallocate(lake_depth_source)
     deallocate(slmsk_source)
     deallocate(land_frac_source)
     deallocate(vegetation_type_pct_source)
     deallocate(vegetation_type_pct_source2)
     deallocate(soil_type_pct_source)
     deallocate(soil_type_pct_source2)
     deallocate(soil_type_source)
     deallocate(substrate_temperature_source)
     deallocate(vegetation_greenness_source)
     deallocate(vegetation_type_source)
     deallocate(sky_albedo_source)
     deallocate(facsf_source)
     deallocate(maximum_snow_albedo_source)

!
!  update lake variables
!
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
     write(*,'(a,a)') 'open to write:',trim(rrfsfile_target)
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("lake_frac",nx_rrfs,ny_rrfs,lakemask_rrfshrrr)
     call rrfs%replace_var("lake_depth",nx_rrfs,ny_rrfs,lake_depth_target)
     call rrfs%replace_var("slmsk",nx_rrfs,ny_rrfs,slmsk_target)
     call rrfs%replace_var("land_frac",nx_rrfs,ny_rrfs,land_frac_target)
     call rrfs%replace_var("vegetation_type_pct",nx_rrfs,ny_rrfs,num_veg_cat,vegetation_type_pct_target)
     call rrfs%replace_var("soil_type_pct",nx_rrfs,ny_rrfs,num_soil_cat,soil_type_pct_target)
     call rrfs%close()

     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.soil_type.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("soil_type",nx_rrfs,ny_rrfs,soil_type_target)
     call rrfs%replace_var("soil_type_pct",nx_rrfs,ny_rrfs,num_soil_cat,soil_type_pct_target2)
     call rrfs%close()

     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.substrate_temperature.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("substrate_temperature",nx_rrfs,ny_rrfs,substrate_temperature_target)
     call rrfs%close()
     !
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.vegetation_greenness.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("vegetation_greenness",nx_rrfs,ny_rrfs,num_time,vegetation_greenness_target)
     call rrfs%close()
     !
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.vegetation_type.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("vegetation_type",nx_rrfs,ny_rrfs,vegetation_type_target)
     call rrfs%replace_var("vegetation_type_pct",nx_rrfs,ny_rrfs,num_veg_cat,vegetation_type_pct_target2)
     call rrfs%close()
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.snowfree_albedo.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("visible_black_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_target(:,:,:,1))
     call rrfs%replace_var("visible_white_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_target(:,:,:,2))
     call rrfs%replace_var("near_IR_black_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_target(:,:,:,3))
     call rrfs%replace_var("near_IR_white_sky_albedo",nx_rrfs,ny_rrfs,num_time,sky_albedo_target(:,:,:,4))
     call rrfs%close()
     !
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.facsf.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("facsf",nx_rrfs,ny_rrfs,facsf_target)
     call rrfs%close()
     !
     rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.maximum_snow_albedo.tile7.'//haloid//'.nc'
     call rrfs%open(trim(rrfsfile_target),"w",200)
     call rrfs%replace_var("maximum_snow_albedo",nx_rrfs,ny_rrfs,maximum_snow_albedo_target)
     call rrfs%close()

!
!  update other halo
!
     if(update_halo4) then
        haloid='halo4'
        nx_rrfs4=nx_rrfs+8
        ny_rrfs4=ny_rrfs+8
        allocate(lakemask_rrfs(nx_rrfs4,ny_rrfs4))
        allocate(lake_depth_source(nx_rrfs4,ny_rrfs4))
        allocate(slmsk_source(nx_rrfs4,ny_rrfs4))
        allocate(land_frac_source(nx_rrfs4,ny_rrfs4))
        allocate(vegetation_type_pct_source(nx_rrfs4,ny_rrfs4,num_veg_cat))
        allocate(soil_type_pct_source(nx_rrfs4,ny_rrfs4,num_soil_cat))

        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs)
        call rrfs%get_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
        call rrfs%get_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
        call rrfs%get_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%get_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source)
        call rrfs%get_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source)
        call rrfs%close()
        lakemask_rrfs(5:nx_rrfs+4,5:ny_rrfs+4)=lakemask_rrfshrrr(1:nx_rrfs,1:ny_rrfs)
        lake_depth_source(5:nx_rrfs+4,5:ny_rrfs+4)=lake_depth_target(1:nx_rrfs,1:ny_rrfs)
        slmsk_source(5:nx_rrfs+4,5:ny_rrfs+4)=slmsk_target(1:nx_rrfs,1:ny_rrfs)
        land_frac_source(5:nx_rrfs+4,5:ny_rrfs+4)=land_frac_target(1:nx_rrfs,1:ny_rrfs)
        vegetation_type_pct_source(5:nx_rrfs+4,5:ny_rrfs+4,:)=vegetation_type_pct_target(1:nx_rrfs,1:ny_rrfs,:)
        soil_type_pct_source(5:nx_rrfs+4,5:ny_rrfs+4,:)=soil_type_pct_target(1:nx_rrfs,1:ny_rrfs,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs)
        call rrfs%replace_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
        call rrfs%replace_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
        call rrfs%replace_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%replace_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source)
        call rrfs%replace_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source)
        call rrfs%close()

        deallocate(lakemask_rrfs)
        deallocate(lake_depth_source)
        deallocate(slmsk_source)
        deallocate(land_frac_source)
        deallocate(vegetation_type_pct_source)
        deallocate(soil_type_pct_source)

        allocate(soil_type_source(nx_rrfs4,ny_rrfs4))
        allocate(soil_type_pct_source2(nx_rrfs4,ny_rrfs4,num_soil_cat))

        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.soil_type.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("soil_type",nx_rrfs4,ny_rrfs4,soil_type_source)
        call rrfs%get_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source2)
        call rrfs%close()

        soil_type_source(5:nx_rrfs+4,5:ny_rrfs+4)=soil_type_target(1:nx_rrfs,1:ny_rrfs)
        soil_type_pct_source2(5:nx_rrfs+4,5:ny_rrfs+4,:)=soil_type_pct_target2(1:nx_rrfs,1:ny_rrfs,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.soil_type.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("soil_type",nx_rrfs4,ny_rrfs4,soil_type_source)
        call rrfs%replace_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source2)
        call rrfs%close()

        deallocate(soil_type_source)
        deallocate(soil_type_pct_source2)

        allocate(substrate_temperature_source(nx_rrfs4,ny_rrfs4))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.substrate_temperature.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("substrate_temperature",nx_rrfs4,ny_rrfs4,substrate_temperature_source)
        call rrfs%close()

        substrate_temperature_source(5:nx_rrfs+4,5:ny_rrfs+4)=substrate_temperature_target(1:nx_rrfs,1:ny_rrfs)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.substrate_temperature.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("substrate_temperature",nx_rrfs4,ny_rrfs4,substrate_temperature_source)
        call rrfs%close()
        deallocate(substrate_temperature_source)
     !
        allocate(vegetation_greenness_source(nx_rrfs4,ny_rrfs4,num_time))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.vegetation_greenness.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("vegetation_greenness",nx_rrfs4,ny_rrfs4,num_time,vegetation_greenness_source)
        call rrfs%close()

        vegetation_greenness_source(5:nx_rrfs+4,5:ny_rrfs+4,:)=vegetation_greenness_target(1:nx_rrfs,1:ny_rrfs,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.vegetation_greenness.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("vegetation_greenness",nx_rrfs4,ny_rrfs4,num_time,vegetation_greenness_source)
        call rrfs%close()
        deallocate(vegetation_greenness_source)
     !
        allocate(vegetation_type_source(nx_rrfs4,ny_rrfs4))
        allocate(vegetation_type_pct_source2(nx_rrfs4,ny_rrfs4,num_veg_cat))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.vegetation_type.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("vegetation_type",nx_rrfs4,ny_rrfs4,vegetation_type_source)
        call rrfs%get_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source2)
        call rrfs%close()

        vegetation_type_source(5:nx_rrfs+4,5:ny_rrfs+4)=vegetation_type_target(1:nx_rrfs,1:ny_rrfs)
        vegetation_type_pct_source2(5:nx_rrfs+4,5:ny_rrfs+4,:)=vegetation_type_pct_target2(1:nx_rrfs,1:ny_rrfs,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.vegetation_type.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("vegetation_type",nx_rrfs4,ny_rrfs4,vegetation_type_source)
        call rrfs%replace_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source2)
        call rrfs%close()
        deallocate(vegetation_type_source)
        deallocate(vegetation_type_pct_source2)

        allocate(sky_albedo_source(nx_rrfs4,ny_rrfs4,num_time,4))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.snowfree_albedo.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("visible_black_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,1))
        call rrfs%get_var("visible_white_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,2))
        call rrfs%get_var("near_IR_black_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,3))
        call rrfs%get_var("near_IR_white_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,4))
        call rrfs%close()

        sky_albedo_source(5:nx_rrfs+4,5:ny_rrfs+4,:,:)=sky_albedo_target(1:nx_rrfs,1:ny_rrfs,:,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.snowfree_albedo.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("visible_black_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,1))
        call rrfs%replace_var("visible_white_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,2))
        call rrfs%replace_var("near_IR_black_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,3))
        call rrfs%replace_var("near_IR_white_sky_albedo",nx_rrfs4,ny_rrfs4,num_time,sky_albedo_source(:,:,:,4))
        call rrfs%close()
        deallocate(sky_albedo_source)
! facsf
        allocate(facsf_source(nx_rrfs4,ny_rrfs4))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.facsf.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("facsf",nx_rrfs4,ny_rrfs4,facsf_source)
        call rrfs%close()

        facsf_source(5:nx_rrfs+4,5:ny_rrfs+4)=facsf_target(1:nx_rrfs,1:ny_rrfs)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.facsf.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("facsf",nx_rrfs4,ny_rrfs4,facsf_source)
        call rrfs%close()
        deallocate(facsf_source)
! maximum_snow_albedo
        allocate(maximum_snow_albedo_source(nx_rrfs4,ny_rrfs4))
        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'.maximum_snow_albedo.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("maximum_snow_albedo",nx_rrfs4,ny_rrfs4,maximum_snow_albedo_source)
        call rrfs%close()

        maximum_snow_albedo_source(5:nx_rrfs+4,5:ny_rrfs+4)=maximum_snow_albedo_target(1:nx_rrfs,1:ny_rrfs)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'.maximum_snow_albedo.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("maximum_snow_albedo",nx_rrfs4,ny_rrfs4,maximum_snow_albedo_source)
        call rrfs%close()
        deallocate(maximum_snow_albedo_source)
     endif

     if(update_halo3) then
        haloid='halo3'
        nx_rrfs4=nx_rrfs+6
        ny_rrfs4=ny_rrfs+6
        allocate(lakemask_rrfs(nx_rrfs4,ny_rrfs4))
        allocate(lake_depth_source(nx_rrfs4,ny_rrfs4))
        allocate(slmsk_source(nx_rrfs4,ny_rrfs4))
        allocate(land_frac_source(nx_rrfs4,ny_rrfs4))

        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs)
        call rrfs%get_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
        call rrfs%get_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
        call rrfs%get_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%close()
        lakemask_rrfs(4:nx_rrfs+3,4:ny_rrfs+3)=lakemask_rrfshrrr(1:nx_rrfs,1:ny_rrfs)
        lake_depth_source(4:nx_rrfs+3,4:ny_rrfs+3)=lake_depth_target(1:nx_rrfs,1:ny_rrfs)
        slmsk_source(4:nx_rrfs+3,4:ny_rrfs+3)=slmsk_target(1:nx_rrfs,1:ny_rrfs)
        land_frac_source(4:nx_rrfs+3,4:ny_rrfs+3)=land_frac_target(1:nx_rrfs,1:ny_rrfs)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs)
        call rrfs%replace_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
        call rrfs%replace_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
        call rrfs%replace_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%close()

        deallocate(lakemask_rrfs)
        deallocate(lake_depth_source)
        deallocate(slmsk_source)
        deallocate(land_frac_source)
     endif

!
!  release memory
!
     deallocate(lakemask_rrfshrrr)
     deallocate(lake_depth_target)
     deallocate(slmsk_target)
     deallocate(land_frac_target)
     deallocate(vegetation_type_pct_target)
     deallocate(vegetation_type_pct_target2)
     deallocate(soil_type_pct_target)
     deallocate(soil_type_pct_target2)
     deallocate(soil_type_target)
     deallocate(substrate_temperature_target)
     deallocate(vegetation_greenness_target)
     deallocate(vegetation_type_target)
     deallocate(sky_albedo_target)
     deallocate(facsf_target)
     write(6,*) "=== LAKE SURGERY REPROCCESS SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!
end program lake_surgery_hrrr


subroutine update_lake_mask(map,nx_rrfs,ny_rrfs,nx_hrrr,ny_hrrr,rlon2d_rrfs,rlat2d_rrfs,&
                            lakemask_hrrr,lakemask_rrfshrrr, &
                            lake_depth_hrrr, lake_depth_rrfshrrr)
!                .      .    .                                       .
! subprogram:   update_lake_mask
!
    use kinds, only: r_single
    use module_map_utils, only : map_util
    implicit none

    type(map_util), intent(in) :: map
    integer, intent(in)        :: nx_rrfs,ny_rrfs,nx_hrrr,ny_hrrr
    real(r_single), intent(in) :: rlon2d_rrfs(nx_rrfs,ny_rrfs)
    real(r_single), intent(in) :: rlat2d_rrfs(nx_rrfs,ny_rrfs)
    real(r_single), intent(inout) :: lakemask_hrrr(nx_hrrr,ny_hrrr)
    real(r_single), intent(inout) :: lake_depth_hrrr(nx_hrrr,ny_hrrr)
    real(r_single), intent(inout) :: lakemask_rrfshrrr(nx_rrfs,ny_rrfs)
    real(r_single), intent(inout) :: lake_depth_rrfshrrr(nx_rrfs,ny_rrfs)
!
    real(r_single) :: xc,yc
    integer  :: ixc,jyc
    integer  :: i,j,n
    integer  :: ixc_walker,jyc_walker
    real(r_single) :: lon_lake,lat_lake
    integer :: nsize
!
! remove the following lakes
! 1. Honey Lake, CA 40.23206° N, 120.31139°
! 2. Carson Lake, NV 39.32246° N, 118.70291° W
! 3. Humboldt Lake, NV, 39.95264° N, 118.62924° W
! 4. Sevier Lake, UT, 38.96909° N, 113.14425° W - remove
!
     do n=1,4
        if(n==1) then
           lon_lake=-120.31139
           lat_lake=40.23206
           nsize=5
        elseif(n==2) then
           lon_lake=-118.70291
           lat_lake=39.32246
           nsize=5
        elseif(n==3) then
           lon_lake=-118.62924
           lat_lake=39.95264
           nsize=5
        elseif(n==4) then
           lon_lake=-113.14425
           lat_lake=38.96909
           nsize=5
        endif
        call map%tll2xy(lon_lake,lat_lake,xc,yc)
        ixc=int(xc+0.5)
        jyc=int(yc+0.5)
        if( (ixc > 0 .and. ixc < nx_hrrr + 1) .and. &
            (jyc > 0 .and. jyc < ny_hrrr + 1) ) then
           write(*,*) 'remove lake=',n,ixc,jyc,lakemask_hrrr(ixc,jyc)
           do j=max(1,jyc-nsize),min(jyc+nsize,ny_hrrr)
              write(*,'(100I3)') (int(lakemask_hrrr(i,j)+0.5),i=max(1,ixc-nsize),min(ixc+nsize,nx_hrrr))
              do i=max(1,ixc-nsize),min(ixc+nsize,nx_hrrr)
                  if(lakemask_hrrr(i,j) > 0.9) then
                    lakemask_hrrr(i,j)=0.0
                    lake_depth_hrrr(i,j)=-2
                  endif
              enddo
           enddo
        endif
     enddo
!
! replace rrfs lake mask with hrrr lake mask
!
! 0. Walker Lake, NV - keep the RRFS old 38.69482° N, 118.71828° W
     call map%tll2xy(-118.71828,38.69482,xc,yc)
     ixc_walker=int(xc+0.5)
     jyc_walker=int(yc+0.5)
     nsize=10
     if( (ixc_walker > 0 .and. ixc_walker < nx_hrrr + 1) .and. &
         (jyc_walker > 0 .and. jyc_walker < ny_hrrr + 1) ) then
         write(*,*) "keep Walker Lake",ixc_walker,jyc_walker,nsize
     endif

     do j=1,ny_rrfs
     do i=1,nx_rrfs
        call map%tll2xy(rlon2d_rrfs(i,j),rlat2d_rrfs(i,j),xc,yc)
        ixc=int(xc+0.5)
        jyc=int(yc+0.5)
        if( (ixc > 0 .and. ixc < nx_hrrr + 1) .and. &
            (jyc > 0 .and. jyc < ny_hrrr + 1) ) then
            if(abs(ixc-ixc_walker)<nsize .and. abs(jyc-jyc_walker)<nsize) then
               write(*,'(a10,2I5,2f6.2)') 'skip',i,j,lakemask_rrfshrrr(i,j),lakemask_hrrr(ixc,jyc)
            else
               lakemask_rrfshrrr(i,j)=lakemask_hrrr(ixc,jyc)
               lake_depth_rrfshrrr(i,j)=lake_depth_hrrr(ixc,jyc)
            endif
        endif
     enddo
     enddo
!
end subroutine update_lake_mask

subroutine change_Outer_Banks(nx_rrfs,ny_rrfs,lakemask_rrfshrrr,land_frac_rrfs,lake_depth_target)
      implicit none
! - Outer Banks - OBX - "Outer BanX" -
!      North Carolina offshore islands usually smaller than dx=3km

      integer, intent(in) :: nx_rrfs,ny_rrfs
      real, intent(in) :: land_frac_rrfs(nx_rrfs,ny_rrfs)
      real, intent(inout) :: lakemask_rrfshrrr(nx_rrfs,ny_rrfs)
      real, intent(inout) :: lake_depth_target(nx_rrfs,ny_rrfs)
      !
      integer,parameter :: n_obs_points=6
      integer :: iobx_points (n_obs_points)
      integer :: jobx_points (n_obs_points)
! --- from south to north, defining i/j points corresponding to the
! Outer Banks islands
      data iobx_points /  &
           3090, 3096, 3108, 3098, 3091, 3074   /
      data jobx_points /  &
            847,  871,  883,  902,  907,  925  /

      real    ::  slope(n_obs_points)   ! slope x/y for each line segment moving northward along Outer Banks
      integer ::  jlen(n_obs_points)
      integer :: i,iy,iybase,j

!  -- set up line segments defining Outer Banks
!      - Define as eastmost point of Pamlico Sound and
!   Albemarle Sound as a function of j starting at j=840
      integer i_east_limit (100),j_east_limit(100)   ! defined relative to j points moving northward from j=840
      real :: x1, x2, y1, y2
      integer :: iseg

      print *, 'n_obs_points =', n_obs_points

!  -- slope needed here is the del x /del y.   The code later will search upward.
!       do iseg = 1, n_obx_points-1
      do iseg = 1, n_obs_points-1
         x1 = float(iobx_points(iseg  ) )
         x2 = float(iobx_points(iseg+1) )
         y1 = float(jobx_points(iseg  ) )
         y2 = float(jobx_points(iseg+1) )
         print *, 'x1,x2,y1,y2 = ',x1,x2,y1,y2
         slope(iseg)  = (x2-x1) / (y2-y1)
         jlen (iseg)  = int( y2 - y1 )
         print *, 'slope, jlen, iseg = ', slope(iseg),jlen(iseg),iseg
      end do

! -- calculate i_obx_limit - eastern I value for Outer Banks islands.
! Function of j

!       i_east_limit(1) = obx_points(2,2)    !   Start with (3090,845)
!       iybase = 844
      iybase = 0
      do iseg = 1, n_obs_points-1
         do iy = 1,jlen(iseg)
            iybase = iybase + 1
            j_east_limit(iybase) = iybase+jobx_points(1)-1
            i_east_limit(iybase) = iobx_points(iseg)  &
                      + int(slope(iseg)*float(iy-1))
           write(*,'(a,5I10)') ' iseg, iy, iybase, yyy, i_east_limit =', &
                iseg, iy,iybase, j_east_limit(iybase), i_east_limit(iybase)
         end do
      end do
      iybase=iybase+1
      j_east_limit(iybase) = jobx_points(n_obs_points)
      i_east_limit(iybase) = iobx_points(n_obs_points)

      Do j = 1,iybase
        print *, j, j_east_limit(j),i_east_limit(j)
        do i=3050,i_east_limit(j)
           if(land_frac_rrfs(i,j_east_limit(j)) < 0.3) then
              lakemask_rrfshrrr(i,j_east_limit(j))=1.0
              lake_depth_target(i,j_east_limit(j))=5.0
           endif
        enddo
      end do

end subroutine 
