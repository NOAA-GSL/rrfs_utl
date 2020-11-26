!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE CONSTANTS_MODULE
!
! This module defines constants that are used by other modules 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cld_parm_array_mod


  use kinds, only: r_kind, i_kind
  implicit none

! set default to private
  private

  public :: metar_impact_radius
  public :: l_metar_impact_radius_change
  public :: metar_impact_radius_max
  public :: metar_impact_radius_min
  public :: metar_impact_radius_max_height
  public :: metar_impact_radius_min_height
  public :: region_dy,region_dx


  real(r_kind)  metar_impact_radius
  logical l_metar_impact_radius_change
  real(r_kind)  metar_impact_radius_max
  real(r_kind)  metar_impact_radius_min
  real(r_kind)  metar_impact_radius_max_height
  real(r_kind)  metar_impact_radius_min_height

  real :: region_dy,region_dx
   
  public :: init_cld_parm

  public :: obstype, sis, nchanl,nreal,ilat,ilon,ndata
  public :: cdata_regular

  character(len=7) :: obstype
  character(len=20):: sis
  integer(i_kind)  :: nchanl
  integer(i_kind)  :: nreal
  integer(i_kind)  :: ilon,ilat
  integer(i_kind)  :: ndata
  real(r_kind),allocatable,dimension(:,:):: cdata_regular

contains

subroutine init_cld_parm


  metar_impact_radius = 10.0_r_kind                 ! in grid
  l_metar_impact_radius_change = .false.            ! .true. =radius change vertically 
  metar_impact_radius_max        = 50000.0_r_kind   ! in meter
  metar_impact_radius_min        = 20000.0_r_kind   ! in meter
  metar_impact_radius_max_height = 3000.0_r_kind    ! in meter
  metar_impact_radius_min_height = 200.0_r_kind     ! in meter
  region_dy=3000.0
  region_dx=3000.0

end subroutine init_cld_parm


end module cld_parm_array_mod
