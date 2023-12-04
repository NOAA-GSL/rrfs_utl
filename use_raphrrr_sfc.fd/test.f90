program testlib
!
  use kinds, only: r_kind,r_single,len_sta_name
  use module_ncio, only : ncio
  use module_map_utils, only : map_util
  use module_time, only : mtime
!
  implicit none

  type(ncio) :: bk
  type(map_util) :: map
  type(mtime) :: mt
  
!
  character(180) :: sfclatlonfile   ! surface station latlon file
  character(180) :: sndlatlonfile   ! sounding station latlon file
  character(180) :: bkfile
!
  real, allocatable :: fld(:,:,:)
  real(r_single),allocatable,target :: rlon2d(:,:),rlat2d(:,:)
  real(r_single) :: rlon,rlat
  real(r_single) :: xc,yc
  integer :: nx,ny,nz
  real(r_single) :: u,v,u0,v0
!
  integer :: idate,idate2
  integer :: JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR,hh,ff

  idate=mt%date2mins(2021,08,30,23,10)
  write(*,*) idate
  call mt%mins2date(idate,IYEAR,MONTH,IDAY,hh,ff)
  write(*,*) IYEAR,MONTH,IDAY,hh,ff
  idate=idate+60
  call mt%mins2date(idate,IYEAR,MONTH,IDAY,hh,ff)
  write(*,*) IYEAR,MONTH,IDAY,hh,ff

  write(*,*)  mt%elapseday(2018,01,30)
  JLDAYN=mt%juliandaynum(2017,11,20)
  write(*,*) JLDAYN
  call mt%julian2date(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
  write(*,*) IYEAR,MONTH,IDAY,IDAYWK,IDAYYR
!
  bkfile='/mnt/lfs4/BMC/rtwbl/mhu/test/process/sst/updatesst/wrf_inout'
  call bk%open(bkfile,"r",200)
  call bk%get_dim("west_east",nx)
  call bk%get_dim("south_north",ny)
  call bk%get_dim("bottom_top",nz)
  write(*,*) 'nx,ny,nz=',nx,ny,nz
!
  allocate(rlon2d(nx,ny))
  allocate(rlat2d(nx,ny))
!
! initialize map transfer parameters
!
  call bk%get_var("XLONG",nx,ny,rlon2d)
  call bk%get_var("XLAT",nx,ny,rlat2d)
  call map%init_general_transform(nx,ny,rlat2d,rlon2d)
!  call obsinfo3d%calculatxy(map)
!  call obsinfo3d%listall()
  call bk%close()
!
end program testlib
