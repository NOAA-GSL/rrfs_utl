program surface_lakesurge_netcdf_mass
!$$$  documentation block
!                .      .    .                                       .
!   Security check for consistency of all land surface parameters on water/ice

!   prgmmr: Ming Hu                 date: 2019-11-18
!
! program history log:
!
!
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
  integer :: nlon, nlat
!
! Declare local parameters

  character(len=120) :: flnm1
  character(len=19)  :: DateStr1
  integer(i_kind)            :: dh1
  
  integer(i_kind) :: i,j,k
  integer(i_kind) :: ndim1
  integer(i_kind) :: WrfType
  integer(i_kind), dimension(4)  :: start_index, end_index
  character (len= 4) :: staggering=' N/A'
  character (len= 3) :: ordering
  character (len=31) :: name,name1,name2,name3,name4,name5
  
  character (len=80), dimension(3)  ::  dimnames
  character (len=80) :: SysDepInfo
  
  character (len=5) :: lutype
  integer(i_kind) :: l, n
  
  integer(i_kind) :: ierr, ier, Status, Status_next_time

! rmse stuff
  
  character (len=31) :: rmse_var
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
  integer(i_kind) nlon_regional,nlat_regional,nsig_regional

  integer(i_kind) wrf_real

!  surface parameters
  real(r_single),allocatable::lakemask(:,:)
  real(r_single),allocatable::tsk(:,:)
  real(r_single),allocatable:: t_lake3d(:,:,:)
!
  real(r_single),allocatable :: xlon(:,:)    !
  real(r_single),allocatable :: ylat(:,:)    !
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
  
  flnm1='wrf_inout'        ! for full cycle
  call ext_ncd_open_for_read( trim(flnm1), 0, 0, "", dh1, Status)
  if ( Status /= 0 )then
     write(6,*)'CONVERT_NETCDF_MASS:  problem with flnm1 = ',&
          trim(flnm1),', Status = ', Status
     stop 74 
  endif

!-------------  get date info

  call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*)' precipiation and snow data from background file at time:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond

  rmse_var='T_LAKE3D'

  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )                !DEDE

  write(6,*)' dh1  = ',dh1         !DEDE
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
  write(6,*)' nlon,lat,sig_regional=',nlon_regional,nlat_regional,nsig_regional

  allocate(xlon(nlon_regional,nlat_regional))
  allocate(ylat(nlon_regional,nlat_regional))
  allocate(tsk(nlon_regional,nlat_regional))
  allocate(lakemask(nlon_regional,nlat_regional))
  allocate(t_lake3d(nlon_regional,nlat_regional,nsig_regional))
!
  write(6,*) '================================================='
  rmse_var='XLONG'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       xlon,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min xlon=',maxval(xlon),minval(xlon)
  write(6,*) '================================================='
  rmse_var='XLAT'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       ylat,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min YLAT=',maxval(ylat),minval(ylat)
  write(6,*) '================================================='
  rmse_var='TSK'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       tsk,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min TSK=',maxval(tsk),minval(tsk)

  write(6,*) '================================================='
  rmse_var='T_LAKE3D'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering       
  write(6,*)' start_index=',start_index     
  write(6,*)' end_index=',end_index         
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       t_lake3d,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  do k=1,nsig_regional
     write(6,*)' max,min t_lake3d=',k,maxval(t_lake3d(:,:,k)),minval(t_lake3d(:,:,k))
  enddo 
  write(6,*) '================================================='
  rmse_var='LAKEMASK'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering  
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index         
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       lakemask,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min LAKEMASK=',maxval(lakemask),minval(lakemask)
  write(6,*) '================================================='

  call ext_ncd_ioclose(dh1, Status)

!
  DO J=1,nlat_regional
  DO I=1,nlon_regional
     if ((xlon(i,j).gt.-117.7 .and. xlon(i,j).lt.-111.5) .and.    &
         (ylat(i,j).gt.40.5   .and. ylat(i,j).lt.42.)    )     then
         if(abs(lakemask(i,j) -1.0) < 0.1 ) then
            if((xlon(i,j).gt.-112.2 .and. xlon(i,j).lt.-112.0) .and.    &
               (ylat(i,j).gt.41.2   .and. ylat(i,j).lt.41.5)    )     then
            else
               tsk(i,j)=max(275.0,tsk(i,j))
               do k=1,nsig_regional
                  t_lake3d(i,j,k)=max(t_lake3d(i,j,k),275.0)
               enddo
            endif
         endif
     endif
!     if ((xlon(i,j).gt.-117.7 .and. xlon(i,j).lt.-111.5) .and.    &
!              ! excludes Willard Bay
!              .not. (xlon(i,j).gt.-112.104 .and. xlon(i,j).lt.-112.100))then
!             if(ylat(i,j).gt.40.5 .and. ylat(i,j).lt.41.22) then
!                if(abs(lakemask(i,j) -1.0) < 0.1) tsk(i,j)=100
!             elseif(( ylat(i,j).ge.41.22 .and. ylat(i,j).lt.42.) .and. .not. &
!                    (ylat(i,j).gt.41.352 .and. ylat(i,j).lt.41.354)) then
!                if(abs(lakemask(i,j) -1.0) < 0.1) tsk(i,j)=100
!             endif ! xlat
!     endif ! xlong
  ENDDO
  ENDDO
!
!           update mass core netcdf file with snow,snowh,snowc
!
  write(6,*) ' ================== '
  write(6,*) ' check for salt lake temperature'
  write(6,*) ' ================== '
  flnm1='wrf_inout'
  call ext_ncd_open_for_update( trim(flnm1), 0, 0, "", dh1, Status)
  if ( Status /= 0 )then
    write(6,*)'UPDATE_NETCDF_MASS:  problem with flnm1 = ',&
         trim(flnm1),', Status = ', Status
    stop 75
  endif
     
!-------------  get date info

  call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*) ' Update sruface parameters in background at time:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond

!-------------  get grid info
  rmse_var='T_LAKE3D'
  call ext_ncd_get_var_info (dh1,rmse_var,ndim1,ordering,staggering, &
                               start_index,end_index, WrfType, ierr    )
   
  write(6,*) '================================================='
  rmse_var='TSK'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       tsk,WrfType,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  rmse_var='T_LAKE3D'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       t_lake3d,WrfType,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  call ext_ncd_ioclose(dh1, Status)
  deallocate(t_lake3d)
  deallocate(tsk)
  deallocate(lakemask)
  deallocate(xlon,ylat)
  
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

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

