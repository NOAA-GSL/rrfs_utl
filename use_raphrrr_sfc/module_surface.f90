!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
module module_surface

  use kinds, only: r_kind,r_single,i_kind,i_byte,i_byte2
  implicit none

  public :: use_surface
!
!-------------------------------------------------------------------------

! set default to private
  private
  type :: use_surface
      integer :: halo
      integer :: nlat,nlon,nlev
      integer :: nlat_target,nlon_target
      integer(i_byte2), allocatable :: index_x(:,:)
      integer(i_byte2), allocatable :: index_y(:,:)
    contains
      procedure :: init
      procedure :: build_mapindex
      procedure :: use_sfc
      procedure :: close
  end type use_surface
!
! constants
!
contains
   
  subroutine use_sfc(this,rapfile,rrfsfile)
!                .      .    .                                       .
! subprogram:   build_mapindex
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
    use module_ncio, only : ncio
  
    implicit none

    character*80,intent(in) :: rapfile
    character*80,intent(in) :: rrfsfile

    class(use_surface) :: this
    type(ncio)     :: raphrrr,rrfs
!
    real(r_single),allocatable,target :: tmp2d4b(:,:)
    real(r_kind),  allocatable,target :: tmp2d8b(:,:)
!
    integer  :: nx_rap,ny_rap
    integer  :: nx_rrfs,ny_rrfs,nz_rrfs
    integer  :: i,j,ix,jx
!

    nx_rrfs=this%nlon
    ny_rrfs=this%nlat
    nz_rrfs=this%nlev
    nx_rap=this%nlon_target
    ny_rap=this%nlat_target

    call raphrrr%open(trim(rapfile),"r",200)
    call rrfs%open(trim(rrfsfile),"r",200)
!
!
    allocate(tmp2d4b(nx_rap,ny_rap))
    call raphrrr%get_var("LANDMASK",nx_rap,ny_rap,tmp2d4b)

    allocate(tmp2d8b(nx_rrfs,ny_rrfs))
    call rrfs%get_var("slmsk",nx_rrfs,ny_rrfs,tmp2d8b)
    tmp2d8b=2
    do j=1,ny_rrfs
      do i=1,nx_rrfs
         ix=this%index_x(i,j)
         jx=this%index_y(i,j)
         if( (ix >= 1 .and. ix <= nx_rap) .and. &
             (jx >= 1 .and. jx <= ny_rap) ) then
            tmp2d8b(i,j)=tmp2d4b(ix,jx)
         endif
      enddo
    enddo
    deallocate(tmp2d4b)
    call rrfs%close()

    call rrfs%open(trim(rrfsfile),"w",200)
    call rrfs%replace_var("slmsk",nx_rrfs,ny_rrfs,tmp2d8b)
    deallocate(tmp2d8b)
!
    call raphrrr%close()
    call rrfs%close()

  end subroutine use_sfc

  subroutine build_mapindex(this,map,rlon,rlat,landmask_this,landmask_target)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    map: grid convert class
!    landmask_this: land mask for this grid
!    landmask_target: land mask for target grid
!
!   output argument list:
!
    use module_map_utils, only : map_util
    implicit none

    class(use_surface) :: this
    type(map_util),  intent(in) :: map
    integer(i_byte), intent(in) :: landmask_this(this%nlon,this%nlat)
    integer(i_byte), intent(in) :: landmask_target(this%nlon_target,this%nlat_target)
    real(r_single),  intent(in) :: rlon(this%nlon,this%nlat)
    real(r_single),  intent(in) :: rlat(this%nlon,this%nlat)
!
    real(r_single) :: xc,yc
    integer  :: ixc,jyc
    integer  :: i,j,ix,jx
    integer  :: ii,jj,ixx,jxx
    real(r_single) :: dist,thisdist
!
    do j=1,this%nlat
    do i=1,this%nlon
       call map%tll2xy(rlon(i,j),rlat(i,j),xc,yc)
       ixc=int(xc+0.5)
       jyc=int(yc+0.5)
       if( (ixc > 0 .and. ixc < this%nlon_target + 1) .and. &
           (jyc > 0 .and. jyc < this%nlat_target + 1) ) then
          ix=min(max(ixc,1),this%nlon_target)
          jx=min(max(jyc,1),this%nlat_target)
          if( landmask_this(i,j)==landmask_target(ix,jx) ) then
             this%index_x(i,j)=ix
             this%index_y(i,j)=jx 
          else
             dist=99999.0
             do jj=max(1,jx-this%halo),min(jx+this%halo,this%nlat_target)
               do ii=max(1,ix-this%halo),min(ix+this%halo,this%nlon_target)
                 if(landmask_this(i,j)==landmask_target(ii,jj)) then
                    thisdist=(jj-jx)*(jj-jx)+(ii-ix)*(ii-ix)
                    thisdist=sqrt(thisdist)
                    if(thisdist <= dist) then
                       dist=thisdist
                       ixx=ii
                       jxx=jj
                    endif
                 endif
               enddo
             enddo
             if(dist < 88888.0) then
               this%index_x(i,j)=ixx
               this%index_y(i,j)=jxx
             else
               !write(*,*) i,j,this%index_x(i,j),this%index_y(i,j)
             endif
          endif
       endif
    enddo
    enddo
!     
!
  end subroutine build_mapindex

  subroutine init(this,nlon,nlat,nlev,nlon_target,nlat_target,halo)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    nlon,nlat - grid dimension
!
!   output argument list:
!
    implicit none

    class(use_surface) :: this
    integer,intent(in) :: halo
    integer,intent(in) :: nlon,nlat,nlev
    integer,intent(in) :: nlon_target,nlat_target

    this%halo=halo
    this%nlon=nlon
    this%nlat=nlat
    this%nlev=nlev
    this%nlon_target=nlon_target
    this%nlat_target=nlat_target
    allocate(this%index_x(this%nlon,this%nlat))
    allocate(this%index_y(this%nlon,this%nlat))
    this%index_x=-999
    this%index_y=-999

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

    class(use_surface) :: this

    this%nlon=0
    this%nlat=0
    this%nlon_target=0
    this%nlat_target=0
    deallocate(this%index_x)
    deallocate(this%index_y)

  end subroutine close

end module module_surface
