program fix_soil_vegetation_pct
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
!  type(use_surface) :: sfc
!
!
!  namelist files
!
  character*80 :: rrfs_lam_source
  character*80 :: rrfs_lam_target
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
  real(r_single),allocatable :: lakemask_rrfs_source(:,:)
  real(r_single),allocatable :: lakemask_rrfs(:,:)
! lake variables
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
  ! inland
  real(r_single),allocatable :: inland_source(:,:)
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
!
  integer :: i,j,k,n
  character*180 :: rrfsfile_source,rrfsfile_target,hrrrfile

  real :: xc,yc
  real :: rsum,rmax
  integer :: ii,jj,iii,jjj,ndist,mindist,nsearch
  integer :: n_lake2land,n_land2lake
  integer :: vegenum, albdonum
  logical :: l_fill_in, l_consist
  logical :: update_halo4,update_halo3
  integer :: ii17,ii14
  logical :: l_convert

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
     l_fill_in=.true.
     l_consist=.true.
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
     call rrfs%get_var("geolon",nx_rrfs,ny_rrfs,rlon2d_rrfs)
     call rrfs%get_var("geolat",nx_rrfs,ny_rrfs,rlat2d_rrfs)
     call rrfs%get_var("lake_frac",nx_rrfs,ny_rrfs,lakemask_rrfs)
     call rrfs%close()
     write(*,*) 'rrfs lon=',maxval(rlon2d_rrfs),minval(rlon2d_rrfs)
     write(*,*) 'rrfs lat=',maxval(rlat2d_rrfs),minval(rlat2d_rrfs)
     write(*,*) 'rrfs lake fraction=',maxval(lakemask_rrfs),minval(lakemask_rrfs)
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
     allocate(inland_source(nx_rrfs,ny_rrfs))

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
     call rrfs%get_var("inland",nx_rrfs,ny_rrfs,inland_source)
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
!   check consitent of lake mask and land mask
     do j=1,ny_rrfs
     do i=1,nx_rrfs
        if(inland_source(i,j) >0.99) then
           if(rlon2d_rrfs(i,j) > 267.6764 .and. rlon2d_rrfs(i,j) < 283.9097 .and. &
              rlat2d_rrfs(i,j) > 41.349211 .and. rlat2d_rrfs(i,j) <  49.0725 ) then
           else
             if( slmsk_source(i,j) > 0.99 .and. lakemask_rrfs(i,j) > 0.5) then
                write(*,*) 'land and lake mismatch =',i,j,rlon2d_rrfs(i,j),rlat2d_rrfs(i,j),slmsk_source(i,j),lakemask_rrfs(i,j)
             endif
             if( slmsk_source(i,j) < 0.09 .and. lakemask_rrfs(i,j) < 0.5)  then
                write(*,*) 'water and no lake mismatch=',i,j,rlon2d_rrfs(i,j),rlat2d_rrfs(i,j),slmsk_source(i,j),lakemask_rrfs(i,j)
                write(*,*) land_frac_source(i,j),lake_depth_source(i,j)
                write(*,'(30f5.2)') vegetation_type_pct_source(i,j,:)
                write(*,'(30f5.2)') soil_type_pct_source(i,j,:)
                ! change this part to lake. This dry lake will be fixed later in lake surgery
                lakemask_rrfs(i,j)=1.0
                lake_depth_target(i,j)=10.0
             endif
           endif
        endif
     enddo
     enddo
     n_lake2land=0
     n_land2lake=0
     nsearch=10
     if(l_consist) then
       vegetation_type_pct_target2=vegetation_type_pct_source
       soil_type_pct_target2=soil_type_pct_source
       do j=1,ny_rrfs
       do i=1,nx_rrfs
! consistent check:soil_type, vegetation_type
          l_convert=.false.
          rsum=0.0
          rmax=-1.0
          ii17=0
          ii14=0
          do k=1,num_veg_cat
             rsum=rsum+vegetation_type_pct_source(i,j,k)
             if(rmax < vegetation_type_pct_source(i,j,k)) then
               rmax=vegetation_type_pct_source(i,j,k)
               iii=k
             endif
          enddo
          if(abs(vegetation_type_target(i,j)-float(iii)) > 0.01) then
            write(*,'(a,2I5,4F8.4)') 'inconsist vege type ',i,j,vegetation_type_target(i,j),float(iii),slmsk_source(i,j),lakemask_rrfs(i,j)
            write(*,'(30f5.2)') vegetation_type_pct_source(i,j,:)
            write(*,'(30f5.2)') vegetation_type_pct_source2(i,j,:)
!            if(abs(vegetation_type_target(i,j)-17.0) < 0.01 .and. iii==15) then
!               write(*,*) "skip the ice point and keep it as water point"
!               vegetation_type_pct_target2(i,j,:)=vegetation_type_pct_source2(i,j,:)
!               soil_type_pct_target2(i,j,:)=soil_type_pct_source2(i,j,:)
!               slmsk_target(i,j)=0.0
!               land_frac_target(i,j)=0.0
!               ii17=17
!               ii14=14
!               write(*,'(30f5.2)')  vegetation_type_pct_target2(i,j,:)
!               write(*,'(30f5.2)')  soil_type_pct_target2(i,j,:)
!               write(*,*) vegetation_type_target(i,j),soil_type_target(i,j)
!               cycle
!            elseif ( iii==15) then
!                    write(*,*) 'grid point changes to ice'
!            else
               vegetation_type_target(i,j)=float(iii)
               ii17=iii
               l_convert=.true.
!            endif
          endif
          if(abs(rsum -1.0) > 0.01) then
             write(*,*) 'sum of veg type pct is not 1', i,j, rsum
          endif

          rsum=0.0
          rmax=-1.0
          do k=1,num_soil_cat
             rsum=rsum+soil_type_pct_source(i,j,k)
             if(rmax < soil_type_pct_source(i,j,k)) then
               rmax=soil_type_pct_source(i,j,k)
               iii=k
             endif
          enddo
          if(abs(soil_type_target(i,j)-float(iii)) > 0.01) then
            write(*,'(a,2I5,4F8.4)') 'inconsist soil type ',i,j,soil_type_target(i,j),float(iii),slmsk_source(i,j),lakemask_rrfs(i,j)
            write(*,'(30f5.2)') soil_type_pct_source(i,j,:)
            write(*,'(30f5.2)') soil_type_pct_source2(i,j,:)
            soil_type_target(i,j)=float(iii)
            ii14=iii
            l_convert=.true.
          endif
          if(abs(rsum -1.0) > 0.01) then
             write(*,*) 'sum of soil type pct is not 1', i,j, rsum
          endif

          if(l_fill_in .and. l_convert) then
! fix vegetation_greenness,substrate_temperature,sky_albedo (from water to land)
             if((abs(vegetation_type_source(i,j)-17.0) <0.01) .and. & 
                (abs(soil_type_source(i,j)-14.0) <0.01) .and. & 
                 abs(slmsk_source(i,j)-1.0) < 0.01) then
                n_lake2land=n_lake2land+1
           ! find the close land point
                mindist=9999
                do jj=max(j-nsearch,1),min(j+nsearch,ny_rrfs)
                do ii=max(i-nsearch,1),min(i+nsearch,nx_rrfs)
                  ndist=(jj-j)*(jj-j)+(ii-i)*(ii-i)
                  if( slmsk_source(ii,jj) > 0.9 .and. mindist > ndist ) then
                    if(substrate_temperature_source(ii,jj) > 100.0 ) then
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
                          mindist=ndist
                          iii=ii
                          jjj=jj
                       endif
                    endif
                 endif
               enddo
               enddo

               if(mindist<1000) then
                 write(*,'(a,2I5,f5.1,a,2I5,f5.1)') "change water ",i,j,slmsk_source(i,j), &
                                                 " to land",iii,jjj,slmsk_source(iii,jjj)
                 if(abs(lakemask_rrfs(i,j))>0.001 .or. abs(lake_depth_target(i,j)) >0.001) then
                         write(*,*) "warning lake,",i,j,lakemask_rrfs(i,j),lake_depth_target(i,j)
                 endif
                 if(abs(slmsk_target(i,j)-1.0)>0.001 .or. abs(land_frac_target(i,j)-1.0) >0.001) then
                         write(*,*) "warning small land fraction,",i,j,slmsk_target(i,j),land_frac_target(i,j)
                 endif
               !      lakemask_rrfs(i,j)=0.0
               !      lake_depth_target(i,j)=0.0
               !      slmsk_target(i,j)=1.0
               !      land_frac_target(i,j)=1.0
                     substrate_temperature_target(i,j)=substrate_temperature_source(iii,jjj)
                     vegetation_greenness_target(i,j,:)=vegetation_greenness_source(iii,jjj,:)
                     sky_albedo_target(i,j,:,:)=sky_albedo_source(iii,jjj,:,:)
                     facsf_target(i,j)=facsf_source(iii,jjj)
                     maximum_snow_albedo_target(i,j)=maximum_snow_albedo_source(iii,jjj)
               else
                  write(*,*) 'cannot find point to fill in missing point',i,j
                  stop 7777
               endif
             else
               if(abs(slmsk_source(i,j)) < 0.01) then
                  write(*,*) "convert land to water=",soil_type_target(i,j),vegetation_type_target(i,j),&
                                soil_type_source(i,j),vegetation_type_source(i,j)
                 n_land2lake=n_land2lake+1
           
                 if(lakemask_rrfs(i,j) <0.749) then
                         write(*,*) "warning no lake,",i,j,lakemask_rrfs(i,j),lake_depth_target(i,j)
                 endif
                 if(abs(slmsk_target(i,j))>0.001 .or. land_frac_target(i,j) >0.501) then
                        write(*,*) "warning large land fraction,",i,j,slmsk_target(i,j),land_frac_target(i,j)
                endif
               !  slmsk_target(i,j)=0.0
               !  land_frac_target(i,j)=0.0

                 substrate_temperature_target(i,j)=-999.0
                 vegetation_greenness_target(i,j,:)=-999.0
                 sky_albedo_target(i,j,:,:)=-999.9
                 facsf_target(i,j)=-999.0
                 maximum_snow_albedo_target(i,j)=-999.0
               else
                 write(*,*) 'No fill for this point',i,j,soil_type_target(i,j),vegetation_type_target(i,j),&
                                soil_type_source(i,j),vegetation_type_source(i,j)
               endif
             endif
          endif
       enddo
       enddo
     endif
     write(*,*) 'number points changing from water to land point=',n_lake2land
     write(*,*) 'number points changing from land to water point=',n_land2lake

     do j=1,ny_rrfs
     do i=1,nx_rrfs
        if( abs(rlat2d_rrfs(i,j)-60.64530)<0.01 .and. abs(rlon2d_rrfs(i,j)-285.0872) <0.01) then
              write(*,*) i,j,soil_type_target(i,j),vegetation_type_target(i,j)
              write(*,*) '1=',lakemask_rrfs(i,j),lake_depth_target(i,j)
              write(*,*) '2=',slmsk_target(i,j),land_frac_target(i,j)
              write(*,*) '3=',vegetation_type_pct_target2(i,j,:)
              write(*,*) '4=',soil_type_pct_target2(i,j,:)
              write(*,*) '5=',substrate_temperature_target(i,j)
              write(*,*) '6=',vegetation_greenness_target(i,j,6)
              write(*,*) '7=',sky_albedo_target(i,j,6,:)
              write(*,*) '8=',facsf_target(i,j),facsf_source(i,j)
              write(*,*) '9=',maximum_snow_albedo_target(i,j),maximum_snow_albedo_source(i,j)
        endif
     enddo
     enddo
!
!  release memory
!
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
     call rrfs%replace_var("lake_frac",nx_rrfs,ny_rrfs,lakemask_rrfs)
     call rrfs%replace_var("lake_depth",nx_rrfs,ny_rrfs,lake_depth_target)
!     call rrfs%replace_var("slmsk",nx_rrfs,ny_rrfs,slmsk_target)
!     call rrfs%replace_var("land_frac",nx_rrfs,ny_rrfs,land_frac_target)
!     call rrfs%replace_var("vegetation_type_pct",nx_rrfs,ny_rrfs,num_veg_cat,vegetation_type_pct_target)
!     call rrfs%replace_var("soil_type_pct",nx_rrfs,ny_rrfs,num_soil_cat,soil_type_pct_target)
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
        allocate(lakemask_rrfs_source(nx_rrfs4,ny_rrfs4))
        allocate(lake_depth_source(nx_rrfs4,ny_rrfs4))
        allocate(slmsk_source(nx_rrfs4,ny_rrfs4))
        allocate(land_frac_source(nx_rrfs4,ny_rrfs4))
        allocate(vegetation_type_pct_source(nx_rrfs4,ny_rrfs4,num_veg_cat))
        allocate(soil_type_pct_source(nx_rrfs4,ny_rrfs4,num_soil_cat))

        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs_source)
        call rrfs%get_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
!        call rrfs%get_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
!        call rrfs%get_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
!        call rrfs%get_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source)
!        call rrfs%get_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source)
        call rrfs%close()
        lakemask_rrfs_source(5:nx_rrfs+4,5:ny_rrfs+4)=lakemask_rrfs(1:nx_rrfs,1:ny_rrfs)
        lake_depth_source(5:nx_rrfs+4,5:ny_rrfs+4)=lake_depth_target(1:nx_rrfs,1:ny_rrfs)
!        slmsk_source(5:nx_rrfs+4,5:ny_rrfs+4)=slmsk_target(1:nx_rrfs,1:ny_rrfs)
!        land_frac_source(5:nx_rrfs+4,5:ny_rrfs+4)=land_frac_target(1:nx_rrfs,1:ny_rrfs)
!        vegetation_type_pct_source(5:nx_rrfs+4,5:ny_rrfs+4,:)=vegetation_type_pct_target(1:nx_rrfs,1:ny_rrfs,:)
!        soil_type_pct_source(5:nx_rrfs+4,5:ny_rrfs+4,:)=soil_type_pct_target(1:nx_rrfs,1:ny_rrfs,:)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs_source)
        call rrfs%replace_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
!        call rrfs%replace_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
!        call rrfs%replace_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
!        call rrfs%replace_var("vegetation_type_pct",nx_rrfs4,ny_rrfs4,num_veg_cat,vegetation_type_pct_source)
!        call rrfs%replace_var("soil_type_pct",nx_rrfs4,ny_rrfs4,num_soil_cat,soil_type_pct_source)
        call rrfs%close()

        deallocate(lakemask_rrfs_source)
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
        allocate(lakemask_rrfs_source(nx_rrfs4,ny_rrfs4))
        allocate(lake_depth_source(nx_rrfs4,ny_rrfs4))
        allocate(slmsk_source(nx_rrfs4,ny_rrfs4))
        allocate(land_frac_source(nx_rrfs4,ny_rrfs4))

        rrfsfile_source=trim(rrfs_lam_source)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_source),"r",200)
        call rrfs%get_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs_source)
        call rrfs%get_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
!        call rrfs%get_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
!        call rrfs%get_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%close()
        lakemask_rrfs_source(4:nx_rrfs+3,4:ny_rrfs+3)=lakemask_rrfs(1:nx_rrfs,1:ny_rrfs)
        lake_depth_source(4:nx_rrfs+3,4:ny_rrfs+3)=lake_depth_target(1:nx_rrfs,1:ny_rrfs)
!        slmsk_source(4:nx_rrfs+3,4:ny_rrfs+3)=slmsk_target(1:nx_rrfs,1:ny_rrfs)
!        land_frac_source(4:nx_rrfs+3,4:ny_rrfs+3)=land_frac_target(1:nx_rrfs,1:ny_rrfs)

        rrfsfile_target=trim(rrfs_lam_target)//"/"//gridid//'_oro_data.tile7.'//haloid//'.nc'
        call rrfs%open(trim(rrfsfile_target),"w",200)
        call rrfs%replace_var("lake_frac",nx_rrfs4,ny_rrfs4,lakemask_rrfs_source)
        call rrfs%replace_var("lake_depth",nx_rrfs4,ny_rrfs4,lake_depth_source)
!        call rrfs%replace_var("slmsk",nx_rrfs4,ny_rrfs4,slmsk_source)
!        call rrfs%replace_var("land_frac",nx_rrfs4,ny_rrfs4,land_frac_source)
        call rrfs%close()

        deallocate(lakemask_rrfs_source)
        deallocate(lake_depth_source)
        deallocate(slmsk_source)
        deallocate(land_frac_source)
     endif

!
!  release memory
!
     deallocate(lakemask_rrfs)
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
     deallocate(maximum_snow_albedo_target)
     write(6,*) "=== LAKE SURGERY REPROCCESS SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!
end program fix_soil_vegetation_pct
