module module_io_fv3lam_bdy
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

  public :: io_fv3lam_bdy

  type :: io_fv3lam_bdy

      character(len=80)   :: bdyfile
      character(len=80)   :: bdyfile_new
      real(r_single)      :: p_top

      integer             :: ndims_bc         !<-- # of dimensions in the BC files
      integer             :: halo_integrate   !<-- Halo depth for FV3 integration
      integer,allocatable :: dimsize_bc(:)    !<-- Dimensions in the BC file.
      character(len=5),allocatable :: dimname_bc(:) ! <-- dimension names in the BC file 

! i/j distribution of the boundary 
      integer :: i_top_bottom_size
      integer :: j_top_bottom_size
      integer :: i_right_left_size
      integer :: j_right_left_size
      integer,allocatable :: i_bottom(:)
      integer,allocatable :: j_bottom(:)
      integer,allocatable :: i_top(:)
      integer,allocatable :: j_top(:)
      integer,allocatable :: i_right(:)
      integer,allocatable :: j_right(:)
      integer,allocatable :: i_left(:)
      integer,allocatable :: j_left(:)
! for u/v w
      integer :: iw_top_bottom_size
      integer :: jw_top_bottom_size
      integer :: iw_right_left_size
      integer :: jw_right_left_size
      integer,allocatable :: i_w_bottom(:)
      integer,allocatable :: j_w_bottom(:)
      integer,allocatable :: i_w_top(:)
      integer,allocatable :: j_w_top(:)
      integer,allocatable :: i_w_right(:)
      integer,allocatable :: j_w_right(:)
      integer,allocatable :: i_w_left(:)
      integer,allocatable :: j_w_left(:)
! for u/v s
      integer :: is_top_bottom_size
      integer :: js_top_bottom_size
      integer :: is_right_left_size
      integer :: js_right_left_size
      integer,allocatable :: i_s_bottom(:)
      integer,allocatable :: j_s_bottom(:)
      integer,allocatable :: i_s_top(:)
      integer,allocatable :: j_s_top(:)
      integer,allocatable :: i_s_right(:)
      integer,allocatable :: j_s_right(:)
      integer,allocatable :: i_s_left(:)
      integer,allocatable :: j_s_left(:)

      character(len=10) :: foursides(4)

      real,allocatable :: bdy_bottom(:,:,:)  ! boundary  
      real,allocatable :: bdy_top(:,:,:)  ! boundary  
      real,allocatable :: bdy_right(:,:,:)  ! boundary  
      real,allocatable :: bdy_left(:,:,:)  ! boundary  

    contains
      procedure :: init
      procedure :: read_bdy_ij
      procedure :: read_bdy
      procedure :: update_bdy
      procedure :: create_new_bdy
      procedure :: close
  end type io_fv3lam_bdy
!
! constants
!
contains

  subroutine init(this,bdyfile)
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
    use module_ncio, only: ncio
    implicit none

    character(len=*),intent(in) :: bdyfile

    class(io_fv3lam_bdy) :: this
    type(ncio) :: fv3io
!
    integer :: k, dimsize

!
    this%bdyfile=trim(bdyfile)
    this%ndims_bc=8
    this%halo_integrate=4
    
    this%foursides(1)='_bottom'
    this%foursides(2)='_top'
    this%foursides(3)='_right'
    this%foursides(4)='_left'

    allocate(this%dimsize_bc(this%ndims_bc))
    allocate(this%dimname_bc(this%ndims_bc))
    
    this%dimname_bc(1)='lon'
    this%dimname_bc(2)='lat'
    this%dimname_bc(3)='lonp'
    this%dimname_bc(4)='latm'
    this%dimname_bc(5)='halo'
    this%dimname_bc(6)='halop'
    this%dimname_bc(7)='lev'
    this%dimname_bc(8)='levp'

    call fv3io%open(trim(bdyfile),'w',0)
    do k=1,this%ndims_bc
       call fv3io%get_dim(trim(this%dimname_bc(k)),dimsize)
       this%dimsize_bc(k)=dimsize
    enddo
    call fv3io%close()
    write(6,*) 'Boundary file =',trim(this%bdyfile)
    write(6,*) 'list of boundary dimensions='
    do k=1,this%ndims_bc
       write(6,'(I5,A10,I5)') k,this%dimname_bc(k),this%dimsize_bc(k)
    enddo

  end subroutine init

  subroutine close(this)
    implicit none
    class(io_fv3lam_bdy) :: this

    this%bdyfile=""
    this%ndims_bc=0
    this%halo_integrate=0
    if(allocated(this%dimsize_bc)) deallocate(this%dimsize_bc)
    if(allocated(this%dimname_bc)) deallocate(this%dimname_bc)

    if(allocated(this%i_bottom)) deallocate(this%i_bottom)
    if(allocated(this%j_bottom)) deallocate(this%j_bottom)
    if(allocated(this%i_top)) deallocate(this%i_top)
    if(allocated(this%j_top)) deallocate(this%j_top)
    if(allocated(this%i_right)) deallocate(this%i_right)
    if(allocated(this%j_right)) deallocate(this%j_right)
    if(allocated(this%i_left)) deallocate(this%i_left)
    if(allocated(this%j_left)) deallocate(this%j_left)

    if(allocated(this%i_w_bottom)) deallocate(this%i_w_bottom)
    if(allocated(this%j_w_bottom)) deallocate(this%j_w_bottom)
    if(allocated(this%i_w_top)) deallocate(this%i_w_top)
    if(allocated(this%j_w_top)) deallocate(this%j_w_top)
    if(allocated(this%i_w_right)) deallocate(this%i_w_right)
    if(allocated(this%j_w_right)) deallocate(this%j_w_right)
    if(allocated(this%i_w_left)) deallocate(this%i_w_left)
    if(allocated(this%j_w_left)) deallocate(this%j_w_left)

    if(allocated(this%i_s_bottom)) deallocate(this%i_s_bottom)
    if(allocated(this%j_s_bottom)) deallocate(this%j_s_bottom)
    if(allocated(this%i_s_top)) deallocate(this%i_s_top)
    if(allocated(this%j_s_top)) deallocate(this%j_s_top)
    if(allocated(this%i_s_right)) deallocate(this%i_s_right)
    if(allocated(this%j_s_right)) deallocate(this%j_s_right)
    if(allocated(this%i_s_left)) deallocate(this%i_s_left)
    if(allocated(this%j_s_left)) deallocate(this%j_s_left)

    if(allocated(this%bdy_bottom)) deallocate(this%bdy_bottom)
    if(allocated(this%bdy_top)) deallocate(this%bdy_top)
    if(allocated(this%bdy_right)) deallocate(this%bdy_right)
    if(allocated(this%bdy_left)) deallocate(this%bdy_left)


  end subroutine close

  subroutine read_bdy_ij(this)
!
    use module_ncio, only: ncio
    implicit none
    class(io_fv3lam_bdy) :: this
    type(ncio) :: fv3io
!
    integer :: lon,lat,lonp,latm,halo,halop
!   
  

    lon=this%dimsize_bc(1)
    lat=this%dimsize_bc(2)
    lonp=this%dimsize_bc(3)
    latm=this%dimsize_bc(4)
    halo=this%dimsize_bc(5)
    halop=this%dimsize_bc(6)
!
! i/j distribution of the boundary 
    this%i_top_bottom_size=lon
    this%j_top_bottom_size=halo
    this%i_right_left_size=halo
    this%j_right_left_size=lat
    allocate(this%i_bottom(this%i_top_bottom_size))
    allocate(this%j_bottom(this%j_top_bottom_size))
    allocate(this%i_top(this%i_top_bottom_size))
    allocate(this%j_top(this%j_top_bottom_size))
    allocate(this%i_right(this%i_right_left_size))
    allocate(this%j_right(this%j_right_left_size))
    allocate(this%i_left(this%i_right_left_size))
    allocate(this%j_left(this%j_right_left_size))
! for u/v w
    this%iw_top_bottom_size=lonp
    this%jw_top_bottom_size=halo
    this%iw_right_left_size=halop
    this%jw_right_left_size=lat
    allocate(this%i_w_bottom(this%iw_top_bottom_size))
    allocate(this%j_w_bottom(this%jw_top_bottom_size))
    allocate(this%i_w_top(this%iw_top_bottom_size))
    allocate(this%j_w_top(this%jw_top_bottom_size))
    allocate(this%i_w_right(this%iw_right_left_size))
    allocate(this%j_w_right(this%jw_right_left_size))
    allocate(this%i_w_left(this%iw_right_left_size))
    allocate(this%j_w_left(this%jw_right_left_size))
! for u/v s
    this%is_top_bottom_size=lon
    this%js_top_bottom_size=halop
    this%is_right_left_size=halo
    this%js_right_left_size=latm
    allocate(this%i_s_bottom(this%is_top_bottom_size))
    allocate(this%j_s_bottom(this%js_top_bottom_size))
    allocate(this%i_s_top(this%is_top_bottom_size))
    allocate(this%j_s_top(this%js_top_bottom_size))
    allocate(this%i_s_right(this%is_right_left_size))
    allocate(this%j_s_right(this%js_right_left_size))
    allocate(this%i_s_left(this%is_right_left_size))
    allocate(this%j_s_left(this%js_right_left_size))

    write(6,'(a,4I10)') 'scaler=',this%i_top_bottom_size,this%j_top_bottom_size,&
                         this%i_right_left_size,this%j_right_left_size
    write(6,'(a,4I10)') 'u/v w=',this%iw_top_bottom_size,this%jw_top_bottom_size,&
                         this%iw_right_left_size,this%jw_right_left_size
    write(6,'(a,4I10)') 'u/v s=',this%is_top_bottom_size,this%js_top_bottom_size,&
                         this%is_right_left_size,this%js_right_left_size
!
    call fv3io%open(trim(this%bdyfile),'r',0)
!
    call fv3io%get_var("i_bottom",this%i_top_bottom_size,this%i_bottom)
    call fv3io%get_var("j_bottom",this%j_top_bottom_size,this%j_bottom)
    call fv3io%get_var("i_top",this%i_top_bottom_size,this%i_top)
    call fv3io%get_var("j_top",this%j_top_bottom_size,this%j_top)
    call fv3io%get_var("i_right",this%i_right_left_size,this%i_right)
    call fv3io%get_var("j_right",this%j_right_left_size,this%j_right)
    call fv3io%get_var("i_left",this%i_right_left_size,this%i_left)
    call fv3io%get_var("j_left",this%j_right_left_size,this%j_left)
! u/v w
    call fv3io%get_var("i_w_bottom",this%iw_top_bottom_size,this%i_w_bottom)
    call fv3io%get_var("j_w_bottom",this%jw_top_bottom_size,this%j_w_bottom)
    call fv3io%get_var("i_w_top",this%iw_top_bottom_size,this%i_w_top)
    call fv3io%get_var("j_w_top",this%jw_top_bottom_size,this%j_w_top)
    call fv3io%get_var("i_w_right",this%iw_right_left_size,this%i_w_right)
    call fv3io%get_var("j_w_right",this%jw_right_left_size,this%j_w_right)
    call fv3io%get_var("i_w_left",this%iw_right_left_size,this%i_w_left)
    call fv3io%get_var("j_w_left",this%jw_right_left_size,this%j_w_left)
! u/v s
    call fv3io%get_var("i_s_bottom",this%is_top_bottom_size,this%i_s_bottom)
    call fv3io%get_var("j_s_bottom",this%js_top_bottom_size,this%j_s_bottom)
    call fv3io%get_var("i_s_top",this%is_top_bottom_size,this%i_s_top)
    call fv3io%get_var("j_s_top",this%js_top_bottom_size,this%j_s_top)
    call fv3io%get_var("i_s_right",this%is_right_left_size,this%i_s_right)
    call fv3io%get_var("j_s_right",this%js_right_left_size,this%j_s_right)
    call fv3io%get_var("i_s_left",this%is_right_left_size,this%i_s_left)
    call fv3io%get_var("j_s_left",this%js_right_left_size,this%j_s_left)
    call fv3io%close
!
!    write(6,*) this%i_bottom
!    write(6,*) this%j_bottom

  end subroutine read_bdy_ij

  subroutine read_bdy(this,varname)
!
    use module_ncio, only: ncio
    implicit none
    character(len=*), intent(in) :: varname
    class(io_fv3lam_bdy) :: this
!
    type(ncio) :: fv3io
!
    integer :: nxi,nyj,nzk
    integer :: i,j,k
    real, allocatable :: r2d4b(:,:)
    character(len=80) :: bdyname
!
! figure out vertical dimension
!
    if(trim(varname)=="zh") then
       nzk=this%dimsize_bc(8)
    elseif(trim(varname)=="ps") then
       nzk=1
    else
       nzk=this%dimsize_bc(7)
    endif
!
    call fv3io%open(trim(this%bdyfile),'r',0)
!
! bottom and top
    if(trim(varname(1:1))=='u' .or. trim(varname(1:1))=='v') then
       if(trim(varname(3:3)) == "w") then
          nxi=this%iw_top_bottom_size
          nyj=this%jw_top_bottom_size
       else
          nxi=this%is_top_bottom_size
          nyj=this%js_top_bottom_size
       endif
    else
       nxi=this%i_top_bottom_size
       nyj=this%j_top_bottom_size
    endif
    if(allocated(this%bdy_bottom)) deallocate(this%bdy_bottom)
    allocate(this%bdy_bottom(nxi,nyj,nzk))
    if(allocated(this%bdy_top)) deallocate(this%bdy_top)
    allocate(this%bdy_top(nxi,nyj,nzk))
    
    if(trim(varname)=='delp' .or. trim(varname)=='delz') then
       this%bdy_bottom=0.0
       this%bdy_top=0.0
    else
       if(nzk==1) then
          allocate(r2d4b(nxi,nyj))
          call fv3io%get_var(trim(varname)//"_bottom",nxi,nyj,r2d4b)
          this%bdy_bottom(:,:,1)=r2d4b
          call fv3io%get_var(trim(varname)//"_top",nxi,nyj,r2d4b)
          this%bdy_top(:,:,1)=r2d4b
          deallocate(r2d4b)
       else
          call fv3io%get_var(trim(varname)//"_bottom",nxi,nyj,nzk,this%bdy_bottom)
          call fv3io%get_var(trim(varname)//"_top",nxi,nyj,nzk,this%bdy_top)
       endif
    endif
! right and left
    if(trim(varname(1:1))=='u' .or. trim(varname(1:1))=='v') then
       if(trim(varname(3:3)) == "w") then
          nxi=this%iw_right_left_size
          nyj=this%jw_right_left_size
       else
          nxi=this%is_right_left_size
          nyj=this%js_right_left_size
       endif
    else
       nxi=this%i_right_left_size
       nyj=this%j_right_left_size
    endif
    if(allocated(this%bdy_right)) deallocate(this%bdy_right)
    allocate(this%bdy_right(nxi,nyj,nzk)) 
    if(allocated(this%bdy_left)) deallocate(this%bdy_left)
    allocate(this%bdy_left(nxi,nyj,nzk))

    if(trim(varname)=='delp' .or. trim(varname)=='delz') then
       this%bdy_right=0.0
       this%bdy_left=0.0
    else
       if(nzk==1) then
          allocate(r2d4b(nxi,nyj))
          call fv3io%get_var(trim(varname)//"_right",nxi,nyj,r2d4b)
          this%bdy_right(:,:,1)=r2d4b
          call fv3io%get_var(trim(varname)//"_left",nxi,nyj,r2d4b)
          this%bdy_left(:,:,1)=r2d4b
          deallocate(r2d4b)
       else
          call fv3io%get_var(trim(varname)//"_right",nxi,nyj,nzk,this%bdy_right)
          call fv3io%get_var(trim(varname)//"_left",nxi,nyj,nzk,this%bdy_left)
       endif
    endif
!
    call fv3io%close
    write(6,*) 'read boundary =',trim(varname), " bottom, top, right, left "
    do k=1,nzk
       write(6,'(I10,10e13.5)') k,maxval(this%bdy_bottom(:,:,k)),minval(this%bdy_bottom(:,:,k)), &
              maxval(this%bdy_top(:,:,k)),minval(this%bdy_top(:,:,k)), &
              maxval(this%bdy_right(:,:,k)),minval(this%bdy_right(:,:,k)), &
              maxval(this%bdy_left(:,:,k)),minval(this%bdy_left(:,:,k))
    enddo

  end subroutine read_bdy

  subroutine update_bdy(this,varname)
!
    use module_ncio, only: ncio
    implicit none
    character(len=*), intent(in) :: varname
    class(io_fv3lam_bdy) :: this
!
    type(ncio) :: fv3io
!
    integer :: nxi,nyj,nzk
    integer :: i,j,k
    real, allocatable :: r2d4b(:,:)
    real, allocatable :: r3d4b(:,:,:)
    character(len=80) :: bdyname
!
! figure out vertical dimension
!
    if(trim(varname)=="zh") then
       nzk=this%dimsize_bc(8)
    elseif(trim(varname)=="ps") then
       nzk=1
    else
       nzk=this%dimsize_bc(7)
    endif
!
    call fv3io%open(trim(this%bdyfile_new),'w',0)
!
! bottom and top
    if(trim(varname(1:1))=='u' .or. trim(varname(1:1))=='v') then
       if(trim(varname(3:3)) == "w") then
          nxi=this%iw_top_bottom_size
          nyj=this%jw_top_bottom_size
       else
          nxi=this%is_top_bottom_size
          nyj=this%js_top_bottom_size
       endif
    else
       nxi=this%i_top_bottom_size
       nyj=this%j_top_bottom_size
    endif
    write(6,*) "update boundary =",trim(varname),nxi,nyj,nzk
    write(6,*) "         bottom =",maxval(this%bdy_bottom),minval(this%bdy_bottom)
    write(6,*) "         top    =",maxval(this%bdy_top),minval(this%bdy_top)
    if(nzk==1) then
       allocate(r2d4b(nxi,nyj))
       r2d4b=this%bdy_bottom(:,:,1)
       call fv3io%replace_var(trim(varname)//"_bottom",nxi,nyj,r2d4b)
       r2d4b=this%bdy_top(:,:,1)
       call fv3io%replace_var(trim(varname)//"_top",nxi,nyj,r2d4b)
       deallocate(r2d4b)
    else
       allocate(r3d4b(nxi,nyj,nzk-1))
       r3d4b(:,:,1:nzk-1)=this%bdy_bottom(:,:,2:nzk)
       call fv3io%replace_var(trim(varname)//"_bottom",nxi,nyj,nzk-1,r3d4b)
       r3d4b(:,:,1:nzk-1)=this%bdy_top(:,:,2:nzk)
       call fv3io%replace_var(trim(varname)//"_top",nxi,nyj,nzk-1,r3d4b)
       deallocate(r3d4b)
    endif
! right and left
    if(trim(varname(1:1))=='u' .or. trim(varname(1:1))=='v') then
       if(trim(varname(3:3)) == "w") then
          nxi=this%iw_right_left_size
          nyj=this%jw_right_left_size
       else
          nxi=this%is_right_left_size
          nyj=this%js_right_left_size
       endif
    else
       nxi=this%i_right_left_size
       nyj=this%j_right_left_size
    endif
    write(6,*) "update boundary =",trim(varname),nxi,nyj
    write(6,*) "          right =",maxval(this%bdy_right),minval(this%bdy_right)
    write(6,*) "          left  =",maxval(this%bdy_left),minval(this%bdy_left)
    if(nzk==1) then
       allocate(r2d4b(nxi,nyj))
       r2d4b=this%bdy_right(:,:,1)
       call fv3io%replace_var(trim(varname)//"_right",nxi,nyj,r2d4b)
       r2d4b=this%bdy_left(:,:,1)
       call fv3io%replace_var(trim(varname)//"_left",nxi,nyj,r2d4b)
       deallocate(r2d4b)
    else
       allocate(r3d4b(nxi,nyj,nzk-1))
       r3d4b(:,:,1:nzk-1)=this%bdy_right(:,:,2:nzk)
       call fv3io%replace_var(trim(varname)//"_right",nxi,nyj,nzk-1,r3d4b)
       r3d4b(:,:,1:nzk-1)=this%bdy_left(:,:,2:nzk)
       call fv3io%replace_var(trim(varname)//"_left",nxi,nyj,nzk-1,r3d4b)
       deallocate(r3d4b)
    endif
!
    call fv3io%close

  end subroutine update_bdy

  subroutine create_new_bdy(this,bdyfile_new)
!
!-----------------------------------------------------------------------
!***  Create a new BC file and prepare its dimensions and variables.
!***  The number of layers is one less than in the original BC file
!***  since the top dummy layer is removed.  All fields will be on
!***  the forecast model layers and not input model layers.
!-----------------------------------------------------------------------
!
    use module_ncio, only: ncio
    use netcdf
    implicit none

    character(len=*),intent(in) :: bdyfile_new
    class(io_fv3lam_bdy) :: this
!
    type(ncio) :: fv3io
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: ncid_bc,ncid_bc_new

      integer,allocatable :: dimsize_bc(:)
      integer,allocatable :: dimid_bc(:)

      integer :: var_id,var_id_new
      character(len=50) :: name_att,name_var
      integer :: ndims,nctype,num_vars_bc,natts
      integer,dimension(1:3) :: dimids=(/0,0,0/)                        &
                               ,dimids_north=(/0,0,0/)                  &
                               ,dimids_south=(/0,0,0/)                  &
                               ,dimids_east =(/0,0,0/)                  &
                               ,dimids_west =(/0,0,0/)
      integer :: n,na,var_id_t
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!

      this%bdyfile_new=trim(bdyfile_new)

!-----------------------------------------------------------------------
!***  Create the new BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_create(bdyfile_new,nf90_clobber,ncid_bc_new))
!
!-----------------------------------------------------------------------
!***  The variables' dimensions.
!-----------------------------------------------------------------------
!
      allocate(dimsize_bc(this%ndims_bc))
      allocate(dimid_bc(this%ndims_bc))
      dimsize_bc=this%dimsize_bc
      dimid_bc=0
      dimsize_bc(7)=this%dimsize_bc(7)-1
      dimsize_bc(8)=this%dimsize_bc(8)-1
!
!-----------------------------------------------------------------------
!***  Define the dimensions in the new BC file.
!-----------------------------------------------------------------------
!
      do n=1,this%ndims_bc
        call check(nf90_def_dim(ncid_bc_new                             &
                               ,this%dimname_bc(n)                           &
                               ,dimsize_bc(n)                           &
                               ,dimid_bc(n)))
      enddo

      deallocate(dimsize_bc)
      deallocate(dimid_bc)
!
!-----------------------------------------------------------------------
!***  How many variables are in a normal BC file?  Get that number
!***  then loop through the variables in a normal BC file and define
!***  those in the new post-GSI BC file.  The exception is that zh
!***  will be skipped since it will be replaced by delz.
!-----------------------------------------------------------------------
!
      call check(nf90_open(this%bdyfile,nf90_nowrite,ncid_bc))
      call check(nf90_inquire(ncid_bc                                   &
                             ,nvariables=num_vars_bc))                     !<--Total # of vbls in a normal BC file.
!
!-----------------------------------------------------------------------
!***  The new file's variables must be defined while that file is
!***  in define mode.  Define each variable in the new file using
!***  those in the original file but exclude zh since instead we
!***  we will want delz.
!-----------------------------------------------------------------------
!
      var_id_new=0

      do n=1,num_vars_bc
        var_id=n
        call check(nf90_inquire_variable(ncid_bc,var_id,name_var,nctype &  !<--Name and type of this variable
                  ,ndims,dimids,natts))                                    !<--# of dimensions and attributes in this variable
!
        if(name_var(1:2)/='zh')then                                        !<--We do not need zh in the new BC file
          var_id_new=var_id_new+1
          call check(nf90_def_var(ncid_bc_new                           &
                                 ,name_var                              &  !<--The variable's name
                                 ,nctype                                &  !<--The variable's type
                                 ,dimids(1:ndims)                       &  !<--The IDs of the variable's dimensions
                                 ,var_id_new))                             !<--The variable's ID
!
!-----------------------------------------------------------------------
!***  Copy each variable's attributes to the new file's variable.
!-----------------------------------------------------------------------
!
          if(natts>0)then
            do na=1,natts
              call check(nf90_inq_attname(ncid_bc,var_id,na,name_att))     !<--Get the attribute's name and ID from restart.
              call check(nf90_copy_att(ncid_bc                          &
                                      ,var_id                           &
                                      ,name_att                         &
                                      ,ncid_bc_new                      &
                                      ,var_id_new))
            enddo
          endif
        endif
!
      enddo

!
!-----------------------------------------------------------------------
!***  We want to add delp and delz as variables in the new BC file. 
!***  In the original regional FV3 the boundary values for delp are 
!***  computed during the vertical remapping of the input data.  In 
!***  the DA process the remapping and wind rotation are bypassed so 
!***  the delp and delz values must be present in the BC file to fill
!***  the BC arrays in the model.  Use the same dimensions as t.
!-----------------------------------------------------------------------
!
!-----------------------
!***  Delp on all sides
!-----------------------
!
      call check(nf90_inq_varid(ncid_bc,'t_bottom',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_north))
!<-- Use T's dimension IDs
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_bottom'                             &  !<--Define delp on the domain bottom (north)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_north(1:ndims)                     &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp bottombndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_top',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_south))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_top'                                &  !<--Define delp on the domain top (south)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_south(1:ndims)                     &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp top bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_left',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_east))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_left'                               &  !<--Define delp on the domain's left (east)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_east(1:ndims)                      &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp left bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_right',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_west))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_right'                              &  !<--Define delp on the domain's right (west)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_west(1:ndims)                      &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp right bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))

!
!-----------------------
!***  Delz on all sides
!-----------------------
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_bottom'                             &  !<--Define delz on the domain bottom (north)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_north(1:ndims)                     &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz bottombndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_top'                                &  !<--Define delz on the domain top (south)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_south(1:ndims)                     &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz top bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_left'                               &  !<--Define delz on the domain's left (east)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_east(1:ndims)                      &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz left bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_right'                              &  !<--Define delz on the domain's right (west)
                             ,NF90_FLOAT                                &  !<--The variable's type
                             ,dimids_west(1:ndims)                      &  !<--The IDs of the variable's dimensions
                             ,var_id))                                     !<--The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz right bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
!-----------------------------------------------------------------------
!***  We are finished with definitions so put the new BC file
!***  into data mode.
!-----------------------------------------------------------------------
!
      call check(nf90_enddef(ncid_bc_new))
      call check(nf90_close(ncid_bc))
!
!-----------------------------------------------------------------------
!
! fill i and j list
!
      
!
      call fv3io%open(trim(this%bdyfile_new),'w',0)
!
      call fv3io%replace_var("i_bottom",this%i_top_bottom_size,this%i_bottom)
      call fv3io%replace_var("j_bottom",this%j_top_bottom_size,this%j_bottom)
      call fv3io%replace_var("i_top",this%i_top_bottom_size,this%i_top)
      call fv3io%replace_var("j_top",this%j_top_bottom_size,this%j_top)
      call fv3io%replace_var("i_right",this%i_right_left_size,this%i_right)
      call fv3io%replace_var("j_right",this%j_right_left_size,this%j_right)
      call fv3io%replace_var("i_left",this%i_right_left_size,this%i_left)
      call fv3io%replace_var("j_left",this%j_right_left_size,this%j_left)
! u/v w
      call fv3io%replace_var("i_w_bottom",this%iw_top_bottom_size,this%i_w_bottom)
      call fv3io%replace_var("j_w_bottom",this%jw_top_bottom_size,this%j_w_bottom)
      call fv3io%replace_var("i_w_top",this%iw_top_bottom_size,this%i_w_top)
      call fv3io%replace_var("j_w_top",this%jw_top_bottom_size,this%j_w_top)
      call fv3io%replace_var("i_w_right",this%iw_right_left_size,this%i_w_right)
      call fv3io%replace_var("j_w_right",this%jw_right_left_size,this%j_w_right)
      call fv3io%replace_var("i_w_left",this%iw_right_left_size,this%i_w_left)
      call fv3io%replace_var("j_w_left",this%jw_right_left_size,this%j_w_left)
! u/v s
      call fv3io%replace_var("i_s_bottom",this%is_top_bottom_size,this%i_s_bottom)
      call fv3io%replace_var("j_s_bottom",this%js_top_bottom_size,this%j_s_bottom)
      call fv3io%replace_var("i_s_top",this%is_top_bottom_size,this%i_s_top)
      call fv3io%replace_var("j_s_top",this%js_top_bottom_size,this%j_s_top)
      call fv3io%replace_var("i_s_right",this%is_right_left_size,this%i_s_right)
      call fv3io%replace_var("j_s_right",this%js_right_left_size,this%j_s_right)
      call fv3io%replace_var("i_s_left",this%is_right_left_size,this%i_s_left)
      call fv3io%replace_var("j_s_left",this%js_right_left_size,this%j_s_left)
      call fv3io%close

  end subroutine create_new_bdy

  subroutine check(status)
     use netcdf, only: nf90_noerr,nf90_strerror
     integer, intent ( in) :: status

     if(status /= nf90_noerr) then
        print *,'ncdf error ', trim(nf90_strerror(status))
        stop
     end if
  end subroutine check

end module module_io_fv3lam_bdy
