
module module_esggrid_util

  use pesg, only : gtoxm_ak_dd,xmtog_ak_dd
  use pkind, only: dp
  use pietc, only: rtod,dtor

  implicit none

  public :: esggrid_util
  public :: edp
!
!-------------------------------------------------------------------------

! set default to private
  private
  integer,parameter  :: edp=dp
  type :: esggrid_util
! define esg grid
      real(dp)              :: a,k
      real(dp)              :: pdlat,pdlon,pdazi
      real(dp)              :: delx,dely
      real(dp)              :: dlat,dlon
      integer               :: lx,ly

      character(len=20) :: grid_type
      logical               :: if_initial=.false.
    contains
      procedure :: init
      procedure :: lltoxy
      procedure :: xytoll
      procedure :: close
  end type esggrid_util
!
! constants
!
contains

  subroutine lltoxy(this,dlon,dlat,x,y)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!         dlat,dlon
!   output argument list:
!         x,y
!
    implicit none

    class(esggrid_util) :: this
    real(dp),intent(in) :: dlat
    real(dp),intent(in) :: dlon
    real(dp),intent(out) :: x
    real(dp),intent(out) :: y

    real(dp),dimension(2) :: xm
    logical               :: ff
 
    if(this%if_initial) then
       x=-1
       y=-1
       call gtoxm_ak_dd(this%A,this%K, this%pdlat,this%pdlon,this%pdazi, &
                        this%delx,this%dely, dlat,dlon, xm, FF)
       if(FF) then
         write(*,*) "WARNING: return bad flag from gtoxm_ak_dd"
       else
         x=(xm(1)+this%lx)*0.5_dp
         y=(xm(2)+this%ly)*0.5_dp
       endif
    else
       write(*,*) "ERROR: lltoxy, not initilizated yet"
       stop 123
    endif

  end subroutine lltoxy

  subroutine xytoll(this,dlon,dlat,x,y)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!         x,y
!   output argument list:
!         dlat,dlon
!
    implicit none

    class(esggrid_util) :: this
    real(dp),intent(out) :: dlat
    real(dp),intent(out) :: dlon
    real(dp),intent(in) :: x
    real(dp),intent(in) :: y

    real(dp),dimension(2) :: xm
    logical               :: ff

    if(this%if_initial) then
       dlat=-1
       dlon=-1
       xm(1)=x*2.0_dp-this%lx
       xm(2)=y*2.0_dp-this%ly
       call xmtog_ak_dd(this%A,this%K, this%pdlat,this%pdlon,this%pdazi, &
                        this%delx,this%dely, xm, dlat,dlon, FF)
       if(FF) then
         write(*,*) "WARNING: return bad flag from gtoxm_ak_dd"
       endif
    else
       write(*,*) "ERROR: lltoxy, not initilizated yet"
       stop 123
    endif

  end subroutine xytoll

  subroutine init(this,grid_type)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!   output argument list:
!
    implicit none

    class(esggrid_util) :: this
    character(len=*),intent(in) :: grid_type
 
! define esg grid
    this%grid_type=trim(grid_type)
    if(trim(grid_type)=="CONUS_13km") then
       this%pdlat = 38.5_dp
       this%pdlon = -97.5_dp
       this%pdazi = 0.0_dp
       this%delx = 0.0580264051680313_dp*dtor
       this%dely = 0.0583034722837032_dp*dtor
       this%a = 0.112806865743845_dp
       this%k = -0.350522864992446_dp
       this%lx=397
       this%ly=233
       this%if_initial=.true.
    else
       write(*,*) 'ERROR: esggrid_util, unknow grid type ', trim(grid_type)
       stop 123
    endif

  end subroutine init

  subroutine close(this)
!                .      .    .                                       .
! subprogram:    close
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    implicit none

    class(esggrid_util) :: this

    this%grid_type=""
    this%pdlat = 0.0_dp
    this%pdlon = 0.0_dp
    this%pdazi = 0.0_dp
    this%delx = 0.0_dp
    this%dely = 0.0_dp
    this%a =0.0_dp
    this%k =0.0_dp
    this%lx=0
    this%ly=0
    this%if_initial=.false.

  end subroutine close

end module module_esggrid_util


