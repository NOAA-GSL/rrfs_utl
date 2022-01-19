program  process_metar_cloud
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2009-09-04
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products and 
!     interpolate them into GSI mass grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use mpi
  use kinds, only: r_kind,i_kind,r_single
  use module_map_utils
!  use misc_definitions_module , only : PROJ_LC, PROJ_ROTLL
  use constants_module ,only : EARTH_RADIUS_M
  use constants, only: init_constants_derived, deg2rad
  use gridmod_gsimap ,only : nlon,nlat,init_general_transform,tll2xy,txy2ll

  use cld_parm_array_mod, only: region_dy,region_dx
  use cld_parm_array_mod, only: metar_impact_radius
  use cld_parm_array_mod, only: l_metar_impact_radius_change, &
                            metar_impact_radius_max,metar_impact_radius_min,&
                            metar_impact_radius_max_height,metar_impact_radius_min_height
  use cld_parm_array_mod, only: init_cld_parm

  use cld_parm_array_mod, only : obstype, sis, nchanl,nreal,ilat,ilon,ndata
  use cld_parm_array_mod, only : cdata_regular

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! MPI variables
  integer :: npe, mype,ierror
!SATID
  real     :: rad2deg = 180.0/3.1415926
!
!  grid
  integer(i_kind) :: nlonfv3,nlatfv3
  real,allocatable:: xlonfv3(:,:)    !
  real,allocatable:: ylatfv3(:,:)    !

  real,allocatable:: xlon(:,:)    !
  real,allocatable:: ylat(:,:)    !
  real(r_kind),allocatable:: rxlon(:,:)    !
  real(r_kind),allocatable:: rylat(:,:)    !

  real ::  userDX, userDY, CEN_LAT, CEN_LON
  real ::  userTRUELAT1,userTRUELAT2,MOAD_CEN_LAT,STAND_LON
  integer :: MAP_PROJ

  REAL :: truelat1, truelat2, stdlon, lat1, lon1, r_earth

!
  real(r_kind),allocatable,dimension(:,:):: cdata_fv3
  integer,allocatable,dimension(:,:):: index_regular
!

  CHARACTER*180   geofile
!
!  For NASA LaRC 
!
  CHARACTER*180   workPath
!     ****VARIABLES FOR THIS NETCDF FILE****
!
! namelist
!
  integer :: analysis_time
  integer :: analysis_minute
  character(len=100) :: prepbufrfile
  real(r_kind)       :: twindin
  namelist/setup/ analysis_time,analysis_minute,prepbufrfile,twindin,&
                  metar_impact_radius,l_metar_impact_radius_change, &
                  metar_impact_radius_max,metar_impact_radius_min, &
                  metar_impact_radius_max_height,metar_impact_radius_min_height
!
!  ** misc
      
  real(r_kind)        :: xc  ! x-grid coordinate (grid units)
  real(r_kind)        :: yc  ! y-grid coordinate (grid units)
  real(r_kind)        :: rlon  ! earth longitude (radians)
  real(r_kind)        :: rlat  ! earth latitude  (radians)

  logical     ::outside     ! .false., then point is inside x-y domain
                              ! .true.,  then point is outside x-y domain

  integer i,j,ii,jj,ij,ndatafv3,maxfv3

  integer :: NCID
!  integer :: status
  logical :: ifexist
  integer :: lunout

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) then


  lunout=15
  call init_cld_parm
!  get namelist
!
  analysis_time=2018051718
  prepbufrfile='prepbufr'
  analysis_minute=0
  twindin=0.5
 
  inquire(file='namelist_metarcld', EXIST=ifexist )
  if(ifexist) then
    open(10,file='namelist_metarcld',status='old')
       read(10,setup)
    close(10)
    write(*,*) 'Namelist setup are:'
    write(*,setup)
  else
    write(*,*) 'No namelist file exist, use default values'
    write(*,*) analysis_time
  endif
 

! set geogrid fle name
!

  call init_constants_derived

  workPath='./'
  write(geofile,'(a,a)') trim(workPath), 'geo_em.d01.nc'

  write(*,*) 'geofile', trim(geofile)
  call GET_DIM_ATT_geo(geofile,nlon,nlat)
  write(*,*) 'NLON,NLAT',nlon,nlat

  call GET_MAP_ATT_geo(geofile, userDX, userDY, CEN_LAT, CEN_LON, &
                userTRUELAT1,userTRUELAT2,MOAD_CEN_LAT,STAND_LON,MAP_PROJ)
  write(*,*) userDX, userDY, CEN_LAT, CEN_LON
  write(*,*) userTRUELAT1,userTRUELAT2,MOAD_CEN_LAT,STAND_LON,MAP_PROJ
  region_dy=userDX
  region_dx=userDY
!
!  get GSI horizontal grid in latitude and longitude
!
  allocate(xlon(nlon,nlat),rxlon(nlon,nlat))
  allocate(ylat(nlon,nlat),rylat(nlon,nlat))

  call OPEN_geo(geofile, NCID)
  call GET_geo_sngl_geo(NCID,Nlon,Nlat,ylat,xlon)
  call CLOSE_geo(NCID)

!
!  get fv3sar grid 
!
  write(geofile,'(a,a)') './', 'fv3sar_grid_spec.nc'
  call GET_DIM_ATT_fv3sar(geofile,NLONFV3,NLATFV3) 
  allocate(xlonfv3(nlonfv3,nlatfv3))
  allocate(ylatfv3(nlonfv3,nlatfv3))
  call OPEN_geo(geofile, NCID)
  call GET_geo_sngl_fv3sar(NCID,Nlonfv3,Nlatfv3,ylatfv3,xlonfv3)
  call CLOSE_geo(NCID)
  write(*,*) 'FV3SAR grid'
  write(*,*) 'nlonfv3,nlatfv3=', nlonfv3,nlatfv3
  write(*,*) 'max, min lon=', maxval(xlonfv3),minval(xlonfv3)
  write(*,*) 'max, min lat=', maxval(ylatfv3),minval(ylatfv3)
!
  mype=0
  rylat=ylat*deg2rad
  rxlon=xlon*deg2rad
  call init_general_transform(rylat,rxlon,mype)
!
  call read_prepbufr_metarcld(prepbufrfile,analysis_time,analysis_minute,twindin)

  write(*,*) 'number of cloud on regular grid=',ndata
  write(*,*) obstype,nreal,nchanl,ilat,ilon,sis
  open(lunout,file='regular_metarcloud.bin',form='unformatted')
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndata
  write(lunout) cdata_regular
  close(lunout)
!
! make a index for the location of the observation in regular grid
  allocate(index_regular(nlon,nlat))
  index_regular=0
  do i=1,ndata
     ii=int(cdata_regular(2,i))
     jj=int(cdata_regular(3,i))
     index_regular(ii,jj)=i
  enddo
!
! fill in the fv3 column with metar obs
!
  maxfv3=nlatfv3*nlonfv3
  allocate(cdata_fv3(nreal,maxfv3))
  cdata_fv3=-9999.0
  ndatafv3=0
  do j=1,nlatfv3
     do i=1,nlonfv3
        rlon=xlonfv3(i,j)
        rlat=ylatfv3(i,j)
        if(xlonfv3(i,j) > 180.0_r_kind) rlon=rlon-360.0_r_kind
        rlon=rlon*deg2rad
        rlat=rlat*deg2rad
        call tll2xy(rlon,rlat,xc,yc)
!
!! find the i,j of FV3 column in WRF grid (nearst point)
        call find_ij(xc,yc,nlon,nlat,ii,jj)

        if(ii > 0 .and. jj > 0 ) then
           ij=index_regular(ii,jj)
           if(ij > 0 .and. ij <= ndata) then
              ndatafv3=ndatafv3+1
              cdata_fv3(:,ndatafv3) = cdata_regular(:,ij)
              cdata_fv3(2,ndatafv3) = i
              cdata_fv3(3,ndatafv3) = j
           endif
        endif
     enddo
  enddo
  deallocate(cdata_regular)
!!
  write(*,*) 'observation number on fv3 grid is=',ndatafv3
! write results
  open(lunout,file='fv3_metarcloud.bin',form='unformatted')
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndatafv3
  write(lunout) cdata_fv3(1:nreal,1:ndatafv3)
  close(lunout)
  deallocate(cdata_fv3)
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

  endif ! if mype==0 

  call MPI_FINALIZE(ierror)

end program process_metar_cloud

subroutine find_ij(xc,yc,nlon,nlat,i,j)
!
! based on input xc, cy to find the nearest grid point
!
  use kinds, only: r_kind,i_kind,r_single

  implicit none

  real(r_kind),intent(in)     :: xc,yc
  integer, intent(in)         :: nlon,nlat
  integer, intent(out)        :: i,j

!
  integer :: ii,jj

! missing value
  i=0
  j=0
!
  ii=int(xc)
  jj=int(yc)
  if( (jj>=1 .and. jj<=nlat-1) .and. (ii>=1 .and. ii <= nlon-1) ) then
     i=ii
     j=jj
     if(xc-ii > 0.5_r_kind) i=ii+1
     if(yc-jj > 0.5_r_kind) j=jj+1
  elseif( ii == 0 .and. (jj>=1 .and. jj<=nlat-1) ) then
     if(xc-ii > 0.5_r_kind) then
       i=1
       j=jj
       if(yc-jj > 0.5_r_kind) j=jj+1
     endif
  elseif( ii == nlon .and. (jj>=1 .and. jj<=nlat-1) ) then
     if(xc-ii < 0.5_r_kind) then
       i=nlon
       j=jj
       if(yc-jj > 0.5_r_kind) j=jj+1
     endif
  elseif( jj == 0 .and. (ii>=1 .and. ii <= nlon-1) ) then
     if(yc-jj > 0.5_r_kind) then
       j=1
       i=ii
       if(xc-ii > 0.5_r_kind) i=ii+1
     endif
  elseif( jj == nlat .and. (ii>=1 .and. ii <= nlon-1) ) then
     if(yc-jj < 0.5_r_kind) then
       j=nlat
       i=ii
       if(xc-ii > 0.5_r_kind) i=ii+1
     endif
  elseif( ii == 0 .and. jj==0 ) then
     if( (xc-ii > 0.5_r_kind) .and. (yc-jj > 0.5_r_kind) ) then
        i=1
        j=1
     endif
  elseif( ii == 0 .and. jj==nlat ) then
     if( (xc-ii > 0.5_r_kind) .and. (yc-jj < 0.5_r_kind) ) then
        i=1
        j=nlat
     endif
  elseif( ii == nlon .and. jj==0 ) then
     if( (xc-ii < 0.5_r_kind) .and. (yc-jj > 0.5_r_kind) ) then
        i=nlon
        j=1
     endif
  elseif( ii == nlon .and. jj==nlat ) then
     if( (xc-ii < 0.5_r_kind) .and. (yc-jj < 0.5_r_kind) ) then
        i=nlon
        j=nlat
     endif
  endif
!  write(*,'(2f10.4,10I6)') xc,yc,ii,jj,i,j

end subroutine find_ij
