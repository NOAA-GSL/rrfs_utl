module module_io_fv3lam_bk
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-03-08
!
! ABSTRACT: 
!     This module read and write fv3lam fields for soil ajustment
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  imssnow
!   OUTPUT FILES: updated surface file
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

  use kinds, only: r_kind,r_single
  implicit none
!
! Rset default to private
!
  
  private

  public :: io_fv3lam_bk

  type :: io_fv3lam_bk

      character(len=80)   :: sfcfile
      character(len=80)   :: gridfile
      character(len=80)   :: tracerfile
      character(len=80)   :: dynfile

      integer             :: nlon
      integer             :: nlat
      integer             :: nlvl
      real                :: p_top
      integer             :: fv3_io_layout_y
      integer,allocatable :: fv3_layout_begin(:)
      integer,allocatable :: fv3_layout_end(:)

      real(r_kind),dimension(:,:),pointer  :: field2d
      real(r_kind),dimension(:,:,:),pointer  :: field3d

    contains
      procedure :: init
      procedure :: setup_grid
      procedure :: read_field
      procedure :: close
  end type io_fv3lam_bk
!
! constants
!
contains

  subroutine init(this,fv3_io_layout_y)
!                .      .    .                                       .
! subprogram: 
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

    integer,intent(in) :: fv3_io_layout_y
    class(io_fv3lam_bk) :: this

!
    this%fv3_io_layout_y=fv3_io_layout_y
!
    this%sfcfile='sfc_data.nc'
    this%gridfile='fv3_grid_spec'
    this%tracerfile='fv_tracer.res.tile1.nc'
    this%dynfile='fv_core.res.tile1.nc'
!
  end subroutine init

  subroutine close(this)
    implicit none
    class(io_fv3lam_bk) :: this

    this%fv3_io_layout_y=0
    if(allocated(this%fv3_layout_begin)) deallocate(this%fv3_layout_begin)
    if(allocated(this%fv3_layout_end)) deallocate(this%fv3_layout_end)

    if(associated(this%field2d)) deallocate(this%field2d)
    if(associated(this%field3d)) deallocate(this%field3d)

  end subroutine close

  subroutine setup_grid(this)
    use module_ncio, only: ncio
    implicit none
    class(io_fv3lam_bk) :: this
    type(ncio) :: fv3io
!
    character(len=80) :: thisfv3file
    integer :: id,iy
    integer :: nlon,nlat,nlvl
    real,allocatable :: r1d4b(:)
!
!
    allocate(this%fv3_layout_begin(this%fv3_io_layout_y))
    allocate(this%fv3_layout_end(this%fv3_io_layout_y))
    this%fv3_layout_begin=0
    this%fv3_layout_end=0
!
    iy=0
    do id=1,this%fv3_io_layout_y
       if(this%fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%dynfile),".",id-1
       else
          thisfv3file=trim(this%dynfile)
       endif
       call fv3io%open(trim(thisfv3file),'r',0)
       call fv3io%get_dim("yaxis_2",nlat)
       if(id==1) then
         call fv3io%get_dim("xaxis_1",nlon)
         call fv3io%get_dim("zaxis_1",nlvl)
       endif
       call fv3io%close
       this%fv3_layout_begin(id)=iy+1
       iy=iy+nlat
       this%fv3_layout_end(id)=iy
    enddo
    this%nlon=nlon
    this%nlvl=nlvl
    this%nlat=this%fv3_layout_end(this%fv3_io_layout_y)
    write(6,'(a,2I10)') " nlon,nlat,nlvl=",this%nlon,this%nlat,this%nlvl
    write(6,'(a20,20I6)') "fv3_layout_begin=",this%fv3_layout_begin
    write(6,'(a20,20I6)') "fv3_layout_end=",this%fv3_layout_end
!
! read model top pressure
    call fv3io%open('./fv_core.res.nc','r',200)
    call fv3io%get_dim("xaxis_1",nlvl)
    allocate(r1d4b(nlvl))
    call fv3io%get_var("ak",nlvl,r1d4b)
    this%p_top=r1d4b(1)
    call fv3io%close
    deallocate(r1d4b)
    write(6,*) 'model top pressure (pa)=',this%p_top


  end subroutine setup_grid

  subroutine read_field(this,varname_bdy)
!
    use module_ncio, only: ncio
    implicit none
    character(len=*),intent(in) :: varname_bdy
    class(io_fv3lam_bk) :: this
    type(ncio) :: fv3io
!
    character(len=80) :: thisfv3file
    character(len=80) :: thisfv3filebase
    integer :: id,fv3_io_layout_y
    integer :: nlon,nlat,nlat_local,nz
    integer :: i,j,k
    real,allocatable :: r3d4b(:,:,:)
    real,allocatable :: r2d4b(:,:)
!
    real,parameter :: grav=9.81
    real,parameter :: rgrav=1./grav
    character(len=20) :: varname
!
    if(trim(varname_bdy)=='ps') then
       varname="delp"
    elseif(trim(varname_bdy)=='w') then
       varname='W'
    elseif(trim(varname_bdy)=='t') then
       varname='T'
    elseif(trim(varname_bdy)=='delz') then
       varname='DZ'
    else
       varname=trim(varname_bdy)
    endif
!
    if(trim(varname)=='u' .or. trim(varname)=='v' .or.   &
       trim(varname)=='W' .or. trim(varname)=='DZ' .or.  &
       trim(varname)=='T' .or. trim(varname)=='delp') then
       thisfv3filebase=trim(this%dynfile)
    else
       thisfv3filebase=trim(this%tracerfile)
    endif
!
    fv3_io_layout_y=this%fv3_io_layout_y   
    nlon=this%nlon
    nlat=this%nlat
    nz=this%nlvl
!
    if(trim(varname)=='u') then
       nlat=nlat+1
    elseif(trim(varname)=='v') then
       nlon=nlon+1
    endif
!
    if(associated(this%field3d)) deallocate(this%field3d)
    allocate(this%field3d(nlon,nlat,nz))
    if(associated(this%field2d)) deallocate(this%field2d)
    allocate(this%field2d(nlon,nlat))

! read    
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(thisfv3filebase),".",id-1
       else
          thisfv3file=trim(thisfv3filebase)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1
       if(trim(varname)=='u') nlat_local=nlat_local+1

       call fv3io%open(trim(thisfv3file),'r',0)

       allocate(r3d4b(nlon,nlat_local,nz))
       call fv3io%get_var(trim(varname),nlon,nlat_local,nz,r3d4b)
       if(trim(varname)=='u') then
          this%field3d(:,this%fv3_layout_begin(id):this%fv3_layout_end(id)+1,:)=r3d4b(:,:,:)
       else
          this%field3d(:,this%fv3_layout_begin(id):this%fv3_layout_end(id),:)=r3d4b(:,:,:)
       endif
       deallocate(r3d4b)

!       if(trim(varname)=='DZ') then
!          allocate(r2d4b(nlon,nlat_local))
!          call fv3io%get_var('phis',nlon,nlat_local,r2d4b)
!          this%field2d(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d4b
!          deallocate(r2d4b)
!       endif

       call fv3io%close
    enddo
    write(6,*) 'read varname=',trim(varname),nlon,nlat,nz
    do k=1,nz
       write(6,*) 'field3d=',k,maxval(this%field3d(:,:,k)),minval(this%field3d(:,:,k))
    enddo
!
! get zh if DZ is read in
!    if(trim(varname)=='DZ') then
!       this%field2d=this%field2d*rgrav  ! terrain in meter
!
!       allocate(r3d4b(nlon,nlat,nz))
!       r3d4b=this%field3d
!       if(associated(this%field3d)) deallocate(this%field3d)
!       allocate(this%field3d(nlon,nlat,nz+1))
!       this%field3d(:,:,nz+1)=this%field2d
!       do k=nz,1,-1   
!          this%field3d(:,:,k)=this%field3d(:,:,k+1)-r3d4b(:,:,k)
!       enddo
!       deallocate(r3d4b)
!       do k=1,nz+1   
!          write(6,*) 'zh=',k,maxval(this%field3d(:,:,k)),minval(this%field3d(:,:,k))
!       enddo
!    endif
!
! get surface pressure if delp is read in
    if(trim(varname)=='ps') then
       if(associated(this%field2d)) deallocate(this%field2d)
       allocate(this%field2d(nlon,nlat))
       this%field2d=this%p_top
       do k=1,nz
          this%field2d=this%field2d+this%field3d(:,:,k)
       enddo
       write(6,*) 'surface pressure=',maxval(this%field2d),minval(this%field2d)

       if(associated(this%field3d)) deallocate(this%field3d)
       allocate(this%field3d(nlon,nlat,1))
       this%field3d(:,:,1)=this%field2d(:,:)
    endif
!
!
  end subroutine read_field

end module module_io_fv3lam_bk
