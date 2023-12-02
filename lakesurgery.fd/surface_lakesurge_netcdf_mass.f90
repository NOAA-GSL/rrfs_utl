program surface_lakesurge_netcdf_mass
!$$$  documentation block
!                .      .    .                                       .
! attributes:
!   language: f90
!
!$$$

  use mpi
  use kinds, only: r_single,i_kind
  implicit none

  INCLUDE 'netcdf.inc'

! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
  integer, parameter :: WRF_INTEGER = 106
!
! Declare local parameters

  character(len=120) :: flnm1
  character(len=120) :: flnm2
  character(len=19)  :: DateStr1
  integer(i_kind)    :: dh1, dh2
  
  integer(i_kind) :: i,j,k
  integer(i_kind) :: ndim1
  integer(i_kind) :: WrfType
  integer(i_kind), dimension(4)  :: start_index, end_index
  character (len= 4) :: staggering=' N/A'
  character (len= 3) :: ordering
  character (len=31) :: name,name1,name2,name3,name4,name5
  
  character (len=80), dimension(3)  ::  dimnames
  character (len=80) :: SysDepInfo
  
  integer(i_kind) :: l, n
  
  integer(i_kind) :: ierr, ier, Status, Status_next_time

! rmse stuff
  
  character (len=31) :: rmse_var
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
  integer(i_kind) nlon_regional,nlat_regional,nsig_regional

  integer(i_kind) wrf_real

!  surface variables
  real(r_single),allocatable::newdepth(:,:)
  real(r_single),allocatable::lakemask(:,:)
  real(r_single),allocatable::lu_index(:,:)
  real(r_single),allocatable::lake_depth(:,:)
  real(r_single),allocatable::lakedepth2d(:,:)
  real(r_single),allocatable::t_lake3d(:,:,:)
!
!
  integer :: iii,jjj,num,numh,i4,j4
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
  if(mype==0) then

  wrf_real=104

!   transfer code from diffwrf for converting netcdf wrf nmm restart file
!      to temporary binary format

  call ext_ncd_ioinit(sysdepinfo,status)
  
  flnm1='wrf_inout_newlakedepth' ! new lake depth from NCAR
  flnm2='wrf_inout_hrrr'         ! HRRR wrfinput_d01 with replaced lake depth


!-------------  get date info from flnm2
  call ext_ncd_open_for_read( trim(flnm2), 0, 0, "", dh2, Status)
  if ( Status /= 0 )then
     write(6,*)'CONVERT_NETCDF_MASS:  problem with flnm2 = ',&
          trim(flnm2),', Status = ', Status
     stop 74 
  endif

  call ext_ncd_get_next_time(dh2, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond

  rmse_var='T_LAKE3D'

  call ext_ncd_get_var_info (dh2,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )                !DEDE

  write(6,*)' d21  = ',dh2         !DEDE
  write(6,*)'file name = ',trim(flnm2)
  write(6,*)'rmse_var = ',trim(rmse_var)
  write(6,*)'ndim1 = ',ndim1
  write(6,*)'ordering = ',trim(ordering)
  write(6,*)'staggering = ',trim(staggering)
  write(6,*)'start_index = ',start_index
  write(6,*)'end_index = ',end_index
  write(6,*)'WrfType = ',WrfType
  write(6,*)'ierr  = ',ierr   !DEDE

  nlon_regional=end_index(1)
  nlat_regional=end_index(2)
  nsig_regional=end_index(3)
  write(6,*)' nlon,nlat,sig_regional=',nlon_regional,nlat_regional,nsig_regional

  allocate(newdepth(nlon_regional,nlat_regional))
  allocate(lake_depth(nlon_regional,nlat_regional))
  allocate(lakedepth2d(nlon_regional,nlat_regional))
  allocate(lu_index(nlon_regional,nlat_regional))
  allocate(lakemask(nlon_regional,nlat_regional))
  allocate(t_lake3d(nlon_regional,nlat_regional,nsig_regional))
!
  write(6,*) '================================================='
  rmse_var='LAKEMASK'
  call ext_ncd_get_var_info (dh2,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)'file name = ',trim(flnm2)
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh2,DateStr1,TRIM(rmse_var),              &
       lakemask,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min LAKEMASK=',maxval(lakemask),minval(lakemask)

  call ext_ncd_ioclose(dh2, Status)

! -- open the file with new lake depths

  call ext_ncd_open_for_read( trim(flnm1), 0, 0, "", dh1, Status)
  if ( Status /= 0 )then
     write(6,*)'CONVERT_NETCDF_MASS:  problem with flnm1 = ',&
          trim(flnm1),', Status = ', Status
     stop 74 
  endif
  call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond

  write(6,*) '================================================='
  rmse_var='LAKE_DEPTH'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' file name = ',trim(flnm1)
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering       
  write(6,*)' start_index=',start_index     
  end_index(3) =1
  write(6,*)' end_index=',end_index         
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       newdepth,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
     write(6,*)' max,min newdepth=',maxval(newdepth(:,:)),minval(newdepth(:,:))
!
  write(6,*) '================================================='
  rmse_var='LU_INDEX'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' file name = ',trim(flnm1)
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       lu_index,WRF_REAL,0,0,0,ordering,    &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min LU_INDEX=',maxval(lu_index),minval(lu_index)

  call ext_ncd_ioclose(dh1, Status)

!
!           update mass core netcdf file with new lake depth
!
  call ext_ncd_open_for_update( trim(flnm2), 0, 0, "", dh2, Status)
  if ( Status /= 0 )then
    write(6,*)'UPDATE_NETCDF_MASS:  problem with flnm2 = ',&
         trim(flnm2),', Status = ', Status
    stop 75
  endif
     
!-------------  get date info

  call ext_ncd_get_next_time(dh2, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*) '================================================='
  write(6,*) ' Update lake depths in HRRR wrfinput:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond

!-------------  get grid info
  rmse_var='T_LAKE3D'
  call ext_ncd_get_var_info (dh2,rmse_var,ndim1,ordering,staggering, &
                               start_index,end_index, WrfType, ierr    )
   
!-- read in old LAKE_DEPTH that will be updated
  write(6,*) '================================================='
  rmse_var='LAKE_DEPTH'
  call ext_ncd_get_var_info (dh2,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' file name = ',trim(flnm1)
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh2,DateStr1,TRIM(rmse_var),              &
       lake_depth,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
     write(6,*)' max,min Old lake_depth=',maxval(lake_depth(:,:)),minval(lake_depth(:,:))

! -- update lake_depth and lakedepth2d with newdepth
  DO J=1,nlat_regional
  DO I=1,nlon_regional
     if(lu_index(i,j) == 21)  then
       if(newdepth(i,j) > 0..and. newdepth(i,j) /= 10.) then
         lake_depth(i,j) = newdepth(i,j)
         lakedepth2d(i,j) = newdepth(i,j)
       else
         lakedepth2d(i,j) = 50.
       endif
     else
       lakedepth2d(i,j) = 0.
     endif

  ENDDO
  ENDDO
!

!-- write updated LAKE_DEPTH
  write(6,*) '================================================='
  rmse_var='LAKE_DEPTH'
  call ext_ncd_get_var_info (dh2,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  end_index(3) = 1
  write(6,*)' end_index=',end_index
  call ext_ncd_write_field(dh2,DateStr1,TRIM(rmse_var),              &
       lake_depth,WrfType,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min New lake_depth=',maxval(lake_depth(:,:)),minval(lake_depth(:,:))
  write(6,*) '================================================='
  rmse_var='LAKEDEPTH2D'
  call ext_ncd_get_var_info (dh2,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_write_field(dh2,DateStr1,TRIM(rmse_var),              &
       lakedepth2d,WrfType,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min New lakedepth2d=',maxval(lakedepth2d(:,:)),minval(lakedepth2d(:,:))
  write(6,*) '================================================='
  call ext_ncd_ioclose(dh2, Status)

  deallocate(t_lake3d)
  deallocate(lakemask)
  deallocate(lu_index)
  deallocate(newdepth)
  deallocate(lake_depth)
  deallocate(lakedepth2d)
  
  write(6,*) "=== HRRR LAKE DEPTH UPDATE SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)

end program surface_lakesurge_netcdf_mass

SUBROUTINE wrf_debug( level , str )
!  USE module_wrf_error
  IMPLICIT NONE
  CHARACTER*(*) str
  INTEGER , INTENT (IN) :: level
  INTEGER               :: debug_level
  CHARACTER (LEN=256) :: time_str
  CHARACTER (LEN=256) :: grid_str
  CHARACTER (LEN=512) :: out_str
!  CALL get_wrf_debug_level( debug_level )
  IF ( level .LE. debug_level ) THEN
    ! old behavior
!      CALL wrf_message( str )
  ENDIF
  write(*,*) 'wrf_debug called !'
  RETURN
END SUBROUTINE wrf_debug

