module module_update_bc
!$$$   module documentation block
!                .      .    .                                       .
! module:   update boundary condition based on the analysis results
! from GSD for RR
!   prgmmr: Hu              org: gsd                date: 2022-03-18
!
! abstract: module for updating surface, soil, moisture from GSD for RR
!
! program history log:
!   2022-03-18  Hu
!
! subroutines included:
!
! Variable Definitions:

  implicit none

! set default to private
  private
! set subroutines to public
  public :: update_bc_4side
  public :: update_bc_4side_gsiwind
  public :: update_bc_4side_wind
  public :: update_bc_fv3uv2earth
! set passed variables to public

contains

subroutine update_bc_4side(fv3bdy,fv3bk,varname)

  use kinds, only: r_kind,i_kind
  use module_io_fv3lam_bdy , only : io_fv3lam_bdy
  use module_io_fv3lam_bk , only : io_fv3lam_bk

  implicit none

  character(len=*), intent(in) :: varname
  type(io_fv3lam_bdy) :: fv3bdy
  type(io_fv3lam_bk) :: fv3bk

!
  integer :: nlon,nlat,nz
  integer :: nxi,nyj
  integer :: i,j,k
  integer :: ii,jj,kk
  integer :: is,ie,js,je
  real,allocatable :: r2d4b(:,:)
  real :: delta
  integer :: fill_halo
!
  nlon=fv3bk%nlon
  nlat=fv3bk%nlat
  nz=fv3bk%nlvl
!
  if(trim(varname)=='zh') then
     nz=nz+1
  elseif(trim(varname)=='ps') then
     nz=1
  endif

  fill_halo=0       ! add delta of the outside restart grid
  if(trim(varname)=='w') then
    fill_halo=1     ! use the outside restart grid as value
  endif
  
!!!!!!!!!!!!!!!!
! bottom
! start grid of the background in boundary area
     nxi=fv3bdy%i_top_bottom_size
     nyj=fv3bdy%j_top_bottom_size
     allocate(r2d4b(nlon,nlat))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bk%field3d(:,:,k)
        call get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_bottom,&
             fv3bdy%j_bottom,fv3bdy%bdy_bottom(:,:,kk),r2d4b,k)
     enddo ! k
     deallocate(r2d4b)

!!!!!!!!!!!!!!!!
! top
! start grid of the background in boundary area
     nxi=fv3bdy%i_top_bottom_size
     nyj=fv3bdy%j_top_bottom_size
     allocate(r2d4b(nlon,nlat))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bk%field3d(:,:,k)
        call get_bc_top(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_top,&
             fv3bdy%j_top,fv3bdy%bdy_top(:,:,kk),r2d4b,k)
     enddo ! k
     deallocate(r2d4b)

!!!!!!!!!!!!!!!!
! left
! start grid of the background in boundary area
     nxi=fv3bdy%i_right_left_size
     nyj=fv3bdy%j_right_left_size

     allocate(r2d4b(nlon,nlat))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bk%field3d(:,:,k)
        call get_bc_left(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_left, &
             fv3bdy%j_left,fv3bdy%bdy_left(:,:,kk),r2d4b,k)
     enddo ! k
     deallocate(r2d4b)
!!!!!!!!!!!!!!!!
! right
! start grid of the background in boundary area
     nxi=fv3bdy%i_right_left_size
     nyj=fv3bdy%j_right_left_size

     allocate(r2d4b(nlon,nlat))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bk%field3d(:,:,k)
        call get_bc_right(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_right, &
             fv3bdy%j_right,fv3bdy%bdy_right(:,:,kk),r2d4b,k)
     enddo ! k
     deallocate(r2d4b)
!!!!!!!!!!
!
end subroutine update_bc_4side

subroutine update_bc_fv3uv2earth(fv3bdy,nlon,nlat,nz,grid_reverse_flag,ud,vd)

  use kinds, only: r_kind,i_kind
  use module_io_fv3lam_bdy ,only : io_fv3lam_bdy
!  use module_io_fv3lam_bk , only : io_fv3lam_bk
  use mod_fv3lam_wind,      only : fv3uv2earth,earthuv2fv3
  use mod_fv3lam_wind,      only : reverse_grid_r_uv,reverse_grid_r_uv_earth

  implicit none

  type(io_fv3lam_bdy) :: fv3bdy
!
  logical,intent(in) :: grid_reverse_flag
  integer,intent(in) :: nlon,nlat,nz
  real(r_kind),intent(inout) :: ud(nlon,nlat+1,nz)
  real(r_kind),intent(inout) :: vd(nlon+1,nlat,nz)
  
  real(r_kind),allocatable :: ua(:,:)
  real(r_kind),allocatable :: va(:,:)
  integer :: k
!
!
  allocate(ua(nlon,nlat),va(nlon,nlat))

  do k=1,nz
     if(.not.grid_reverse_flag) then
        call reverse_grid_r_uv(ud(:,:,k),nlon,nlat+1,1)
        call reverse_grid_r_uv(vd(:,:,k),nlon+1,nlat,1)
     endif

     call fv3uv2earth(ud(:,:,k),vd(:,:,k),nlon,nlat,ua,va)

     if(.not.grid_reverse_flag) then
        call reverse_grid_r_uv_earth(ua,nlon,nlat,1)
        call reverse_grid_r_uv_earth(va,nlon,nlat,1)
     endif

     call agrid_to_dgrid_u(nlon,nlat,ua,ud(:,:,k))
     call agrid_to_dgrid_v(nlon,nlat,va,vd(:,:,k))
  enddo ! k
! 
  deallocate(ua,va)
  write(6,*) 'ud on earth'
  do k=1,nz
      write(6,*) k,maxval(ud(:,:,k)),minval(ud(:,:,k))
  enddo
  write(6,*) 'vd on earth'
  do k=1,nz
      write(6,*) k,maxval(vd(:,:,k)),minval(vd(:,:,k))
  enddo

end subroutine update_bc_fv3uv2earth

subroutine update_bc_4side_wind(fv3bdy,nlon_in,nlat_in,nz,uv,varname)

  use kinds, only: r_kind,i_kind
  use module_io_fv3lam_bdy , only : io_fv3lam_bdy

  implicit none

  character(len=*), intent(in) :: varname
  type(io_fv3lam_bdy) :: fv3bdy
!
  integer,intent(in) :: nlon_in,nlat_in,nz
  real(r_kind),intent(in) :: uv(nlon_in,nlat_in,nz) 
  integer :: nlon,nlat

  integer :: nxi,nyj
  integer :: i,j,k
  integer :: ii,jj,kk
  integer :: is,ie,js,je
  real,allocatable :: r2d4b(:,:)
  real,allocatable :: r2dwind(:,:)
  real :: delta
  integer :: fill_halo
!
  fill_halo=0     ! use the outside restart grid as value
  
  if(trim(varname) == 'u_w' .or. trim(varname) == 'v_w') then
!
!  process each level
!
     nlon=nlon_in
     nlat=nlat_in
     if(trim(varname) == 'u_w') then  ! make from u_s
        nlon=nlon_in+1
        nlat=nlat_in-1
     endif
     
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'v_w') then
           r2dwind=uv(:,:,k)
        elseif(trim(varname) == 'u_w') then
           allocate(r2d4b(nlon-1,nlat+1))
           r2d4b=uv(:,:,k)
           call dgrid_to_cgrid_u(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
!
        write(6,*) trim(varname),k,maxval(r2dwind),minval(r2dwind)
!!!!!!!!!!!!!!!!!!!!
! bottom u/v w
!
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        call get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_bottom,&
             fv3bdy%j_w_bottom,fv3bdy%bdy_bottom(:,:,kk),r2dwind,k)

!!!!!!!!!!!!!!!!!
! top u/v w
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        call get_bc_top(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_top,&
             fv3bdy%j_w_top,fv3bdy%bdy_top(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!
! left u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        call get_bc_left(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_left,&
             fv3bdy%j_w_left,fv3bdy%bdy_left(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!
! right u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        call get_bc_right(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_right,&
             fv3bdy%j_w_right,fv3bdy%bdy_right(:,:,kk),r2dwind,k)
     enddo

     deallocate(r2dwind)
  
  elseif(trim(varname) == 'u_s' .or. trim(varname) == 'v_s') then
!
!  process each level
!    
     nlon=nlon_in
     nlat=nlat_in
     if(trim(varname) == 'v_s') then  ! make from v_w
        nlon=nlon_in-1
        nlat=nlat_in+1
     endif
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'u_s') then
           r2dwind=uv(:,:,k)
        elseif(trim(varname) == 'v_s') then
           allocate(r2d4b(nlon+1,nlat-1))
           r2d4b=uv(:,:,k)
           call dgrid_to_cgrid_v(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
        write(6,*) trim(varname),k,maxval(r2dwind),minval(r2dwind)
!!!!!!!!!!!!!!!!!
! bottom u/v s
!
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        call get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_bottom,&
             fv3bdy%j_s_bottom,fv3bdy%bdy_bottom(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!!!!!!!
! top u/v s
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        call get_bc_top(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_top,&
             fv3bdy%j_s_top,fv3bdy%bdy_top(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!!!!
! left u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        call get_bc_left(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_left,&
             fv3bdy%j_s_left,fv3bdy%bdy_left(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!
! right u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        call get_bc_right(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_right,&
             fv3bdy%j_s_right,fv3bdy%bdy_right(:,:,kk),r2dwind,k)
     enddo

     deallocate(r2dwind)

  endif

end subroutine update_bc_4side_wind

subroutine update_bc_4side_gsiwind(fv3bdy,fv3bk,varname)

  use kinds, only: r_kind,i_kind
  use module_io_fv3lam_bdy , only : io_fv3lam_bdy
  use module_io_fv3lam_bk , only : io_fv3lam_bk

  implicit none

  character(len=*), intent(in) :: varname
  type(io_fv3lam_bdy) :: fv3bdy
  type(io_fv3lam_bk) :: fv3bk

!
  integer :: nlon,nlat,nz
  integer :: nxi,nyj
  integer :: i,j,k
  integer :: ii,jj,kk
  integer :: is,ie,js,je
  real,allocatable :: r2d4b(:,:)
  real,allocatable :: r2dwind(:,:)
  real :: delta
  integer :: fill_halo
!
  nlon=fv3bk%nlon
  nlat=fv3bk%nlat
  nz=fv3bk%nlvl
!
  if(trim(varname)=='u_s') then
     nlat=nlat+1
  else if(trim(varname)=='u_w') then
     nlon=nlon+1
  elseif(trim(varname)=='v_s') then
     nlat=nlat+1
  elseif(trim(varname)=='v_w') then
     nlon=nlon+1
  endif

  fill_halo=1     ! use the outside restart grid as value
  
  if(trim(varname) == 'u_w' .or. trim(varname) == 'v_w') then
!
!  process each level
!
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'v_w') then
           r2dwind=fv3bk%field3d(:,:,k)
        elseif(trim(varname) == 'u_w') then
           allocate(r2d4b(nlon-1,nlat+1))
           r2d4b=fv3bk%field3d(:,:,k)
           call dgrid_to_cgrid_u(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
!!!!!!!!!!!!!!!!!!!!
! bottom u/v w
!
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        call get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_bottom,&
             fv3bdy%j_w_bottom,fv3bdy%bdy_bottom(:,:,kk),r2dwind,k)

!!!!!!!!!!!!!!!!!
! top u/v w
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        call get_bc_top(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_top,&
             fv3bdy%j_w_top,fv3bdy%bdy_top(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!
! left u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        call get_bc_left(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_left,&
             fv3bdy%j_w_left,fv3bdy%bdy_left(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!
! right u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        call get_bc_right(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_w_right,&
             fv3bdy%j_w_right,fv3bdy%bdy_right(:,:,kk),r2dwind,k)
     enddo

     deallocate(r2dwind)

  elseif(trim(varname) == 'u_s' .or. trim(varname) == 'v_s') then
!
!  process each level
!    
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'u_s') then
           r2dwind=fv3bk%field3d(:,:,k)
        elseif(trim(varname) == 'v_s') then
           allocate(r2d4b(nlon+1,nlat-1))
           r2d4b=fv3bk%field3d(:,:,k)
           call dgrid_to_cgrid_v(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
!!!!!!!!!!!!!!!!!
! bottom u/v s
!
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        call get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_bottom,&
             fv3bdy%j_s_bottom,fv3bdy%bdy_bottom(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!!!!!!!
! top u/v s
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        call get_bc_top(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_top,&
             fv3bdy%j_s_top,fv3bdy%bdy_top(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!!!!!
! left u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        call get_bc_left(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_left,&
             fv3bdy%j_s_left,fv3bdy%bdy_left(:,:,kk),r2dwind,k)
!!!!!!!!!!!!!!!
! right u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        call get_bc_right(nxi,nyj,nlon,nlat,fill_halo,fv3bdy%i_s_right,&
             fv3bdy%j_s_right,fv3bdy%bdy_right(:,:,kk),r2dwind,k)
     enddo

     deallocate(r2dwind)

  endif

end subroutine update_bc_4side_gsiwind

subroutine agrid_to_cgrid_u(nlon,nlat,ua,uc)
!
!-----------------------------------------------------------------------
! from A grid to D grid u
!-----------------------------------------------------------------------
!
  use kinds, only: r_kind,i_kind
  implicit none

      integer,intent(in)  :: nlon,nlat
      real(r_kind), intent(in)    :: ua(nlon,nlat)
      real(r_kind), intent(inout) :: uc(nlon+1,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

     do j=1,nlat
        uc(1,j)=ua(1,j)
        do i=2,nlon
           uc(i,j)=0.5*(ua(i,j)+ua(i-1,j))
        end do
        uc(nlon+1,j)=ua(nlon,j)
     end do
!
end subroutine agrid_to_cgrid_u

subroutine agrid_to_cgrid_v(nlon,nlat,va,vc)
!
!-----------------------------------------------------------------------
! from A grid to D grid v
!-----------------------------------------------------------------------
!
  use kinds, only: r_kind,i_kind
  implicit none

      integer,intent(in)  :: nlon,nlat
      real(r_kind), intent(in)    :: va(nlon,nlat)
      real(r_kind), intent(inout) :: vc(nlon,nlat+1)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

     j=1
     do i=1,nlon
        vc(i,j)= va(i,j)
     end do

     do j=2,nlat
        do i=1,nlon
           vc(i,j)=0.5*(va(i,j)+va(i,j-1))
        end do
     end do

     j=nlat
     do i=1,nlon
        vc(i,j+1)= va(i,j)
     end do
!
end subroutine agrid_to_cgrid_v

subroutine agrid_to_dgrid_u(nlon,nlat,ua,ud)
!
!-----------------------------------------------------------------------
! from A grid to D grid u
!-----------------------------------------------------------------------
!
  use kinds, only: r_kind,i_kind
  implicit none

      integer,intent(in)  :: nlon,nlat
      real(r_kind), intent(in)    :: ua(nlon,nlat)
      real(r_kind), intent(inout) :: ud(nlon,nlat+1)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

     j=1
     do i=1,nlon
        ud(i,j)= ua(i,j)
     end do

     do j=2,nlat
        do i=1,nlon
           ud(i,j)=0.5*(ua(i,j)+ua(i,j-1))
        end do
     end do

     j=nlat
     do i=1,nlon
        ud(i,j+1)= ua(i,j)
     end do
!
end subroutine agrid_to_dgrid_u

subroutine agrid_to_dgrid_v(nlon,nlat,va,vd)
!
!-----------------------------------------------------------------------
! from A grid to D grid v
!-----------------------------------------------------------------------
!
  use kinds, only: r_kind,i_kind
  implicit none

      integer,intent(in)  :: nlon,nlat
      real(r_kind), intent(in)    :: va(nlon,nlat)
      real(r_kind), intent(inout) :: vd(nlon+1,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

     do j=1,nlat
        vd(1,j)=va(1,j)
        do i=2,nlon
           vd(i,j)=0.5*(va(i,j)+va(i-1,j))
        end do
        vd(nlon+1,j)=va(nlon,j)
     end do
!
end subroutine agrid_to_dgrid_v

subroutine dgrid_to_cgrid_u(nlon,nlat,ud,uc)
!
!-----------------------------------------------------------------------
!***  The GSI updates the D-grid winds but the BC file also needs the
!***  C-grid winds.   Compute the C-grid winds in the boundary rows
!***  by interpolating from the D-grid winds and write them into the 
!***  BC file.
!-----------------------------------------------------------------------
!
  use kinds, only: r_kind,i_kind
  implicit none

      integer,intent(in)  :: nlon,nlat
      real, intent(in)    :: ud(nlon-1,nlat+1)
      real, intent(inout) :: uc(nlon,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

      do j=1,nlat
          do i=2,nlon-1
            uc(i,j)=0.25*(ud(i-1,j)+ud(i-1,j+1)      &
                         +ud(i  ,j)+ud(i  ,j+1))
          enddo
          uc(1,j)=2.*uc(2,j)-uc(3,j)
          uc(nlon,j)=2.*uc(nlon-1,j)-uc(nlon-2,j)
      enddo

!
end subroutine dgrid_to_cgrid_u

subroutine dgrid_to_cgrid_v(nlon,nlat,vd,vc)
!
!-----------------------------------------------------------------------
!***  The GSI updates the D-grid winds but the BC file also needs the
!***  C-grid winds.   Compute the C-grid winds in the boundary rows
!***  by interpolating from the D-grid winds and write them into the 
!***  BC file.
!-----------------------------------------------------------------------
!
  implicit none

      integer,intent(in)  :: nlon,nlat
      real, intent(in)    :: vd(nlon+1,nlat-1)
      real, intent(inout) :: vc(nlon,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

      do j=2,nlat-1
          do i=1,nlon
            vc(i,j)=0.25*(vd(i  ,j-1)+vd(i  ,j)    &
                         +vd(i+1,j-1)+vd(i+1,j))
          enddo
      enddo
      do i=1,nlon
        vc(i,1)=2.*vc(i,2)-vc(i,3)
        vc(i,nlat)=2.*vc(i,nlat-1)-vc(i,nlat-2)
      enddo

!
end subroutine dgrid_to_cgrid_v

subroutine get_bc_bottom(nxi,nyj,nlon,nlat,fill_halo,&
                         i_bottom,j_bottom,bdy_bottom,field2d,k)
!
!  process bottom boundary 
!
  use kinds, only: r_kind,i_kind
  implicit none
!
  integer,intent(in) :: nlon,nlat
  integer,intent(in) :: nxi,nyj
  integer,intent(in) :: i_bottom(nxi),j_bottom(nyj)
  integer,intent(in) :: fill_halo
  real,intent(in)    :: field2d(nlon,nlat)
  real,intent(inout) :: bdy_bottom(nxi,nyj)
  integer,intent(in) :: k
!
  real,allocatable :: r2d4b(:,:)
  integer :: is,ie,js,je

  real :: delta
  integer :: i,j
  integer :: ii,jj
!

  je=nyj
  do jj=1,nyj
     do ii=1,nxi
           i=i_bottom(ii)
           j=j_bottom(jj)
           if( i==1 .and. j==1) then
              is=ii
              js=jj
           endif
           if( i==nlon .and. j==1) then
              ie=ii
           endif
      enddo
   enddo

   if(k==1) write(6,'(a50,10I5)') 'update bottom is,ie,js,je,nxi,nyj', &
                          is,ie,js,je,nxi,nyj
!
  allocate(r2d4b(nxi,nyj))
  r2d4b=bdy_bottom

        do jj=1,nyj
           do ii=1,nxi
              i=i_bottom(ii)
              j=j_bottom(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=field2d(i,j)
              endif
           enddo
        enddo
!
! halo update
        if(fill_halo==1) then
           do i=is,ie
              do j=1,js-1
                 r2d4b(i,j)=r2d4b(i,js)
              enddo
           enddo
           do j=1,nyj
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(is,j)
              enddo
              do i=ie+1,nxi
                 r2d4b(i,j)=r2d4b(ie,j)
              enddo
           enddo
        elseif(fill_halo==0) then
! get bottom halo update
           do i=is,ie
              delta=r2d4b(i,js)-bdy_bottom(i,js)
              do j=1,js-1
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo

           do j=1,nyj
! get left halo update
              delta=r2d4b(is,j)-bdy_bottom(is,j)
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
! get right halo update
              delta=r2d4b(ie,j)-bdy_bottom(ie,j)
              do i=ie+1,nxi
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo
        else
           ! no halo update, keep halo as is
        endif

   bdy_bottom=r2d4b
   deallocate(r2d4b)
!
end subroutine get_bc_bottom

subroutine get_bc_top(nxi,nyj,nlon,nlat,fill_halo,&
                         i_top,j_top,bdy_top,field2d,k)
!
!  process top boundary 
!
  use kinds, only: r_kind,i_kind
  implicit none
!
  integer,intent(in) :: nlon,nlat
  integer,intent(in) :: nxi,nyj
  integer,intent(in) :: i_top(nxi),j_top(nyj)
  integer,intent(in) :: fill_halo
  real,intent(in)    :: field2d(nlon,nlat)
  real,intent(inout) :: bdy_top(nxi,nyj)
  integer,intent(in) :: k
!
  real,allocatable :: r2d4b(:,:)
  integer :: is,ie,js,je

  real :: delta
  integer :: i,j
  integer :: ii,jj
!


     je=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=i_top(ii)
           j=j_top(jj)
           if( i==1 .and. j==nlat) then
              is=ii
              js=jj
           endif
           if( i==nlon .and. j==nlat) then
              ie=ii
           endif
        enddo
     enddo
     if(k==1) write(6,'(a50,10I5)') 'update top is,ie,js,je,nxi,nyj:',is,ie,js,je,nxi,nyj
     allocate(r2d4b(nxi,nyj))

        r2d4b=bdy_top
        do jj=1,nyj
           do ii=1,nxi
              i=i_top(ii)
              j=j_top(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=field2d(i,j)
              endif
           enddo
        enddo
! halo update
        if(fill_halo==1) then
           do i=is,ie
              do j=js+1,je
                 r2d4b(i,j)=r2d4b(i,js)
              enddo
           enddo
           do j=1,nyj
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(is,j)
              enddo
              do i=ie+1,nxi
                 r2d4b(i,j)=r2d4b(ie,j)
              enddo
           enddo
        elseif(fill_halo==0) then
! get top halo update
           do i=is,ie
              delta=r2d4b(i,js)-bdy_top(i,js)
              do j=js+1,je
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo
!
           do j=1,nyj
! get left halo update
              delta=r2d4b(is,j)-bdy_top(is,j)
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
! get right halo update
              delta=r2d4b(ie,j)-bdy_top(ie,j)
              do i=ie+1,nxi
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo
        else
! 
        endif
!
        bdy_top=r2d4b
     deallocate(r2d4b)

end subroutine get_bc_top

subroutine get_bc_left(nxi,nyj,nlon,nlat,fill_halo,&
                         i_left,j_left,bdy_left,field2d,k)
!
!  process top boundary 
!
  use kinds, only: r_kind,i_kind
  implicit none
!
  integer,intent(in) :: nlon,nlat
  integer,intent(in) :: nxi,nyj
  integer,intent(in) :: i_left(nxi),j_left(nyj)
  integer,intent(in) :: fill_halo
  real,intent(in)    :: field2d(nlon,nlat)
  real,intent(inout) :: bdy_left(nxi,nyj)
  integer,intent(in) :: k
!
  real,allocatable :: r2d4b(:,:)
  integer :: is,ie,js,je

  real :: delta
  integer :: i,j
  integer :: ii,jj


     js=1
     je=nyj
     ie=nxi
     do ii=1,nxi
        i=i_left(ii)
        if( i==1) is=ii
     enddo

     if(k==1) write(6,'(a50,10I5)') 'update left is,ie,js,je,nxi,nyj:',is,ie,js,je,nxi,nyj
     allocate(r2d4b(nxi,nyj))

        r2d4b=bdy_left
        do jj=1,nyj
           do ii=1,nxi
              i=i_left(ii)
              j=j_left(jj)
              if( (i >=1 .and. i<=nlon) .and.    &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=field2d(i,j)
              endif
           enddo
        enddo
!
        if(fill_halo==1) then
           do j=1,nyj
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(is,j)
              enddo
           enddo
        else if(fill_halo==0) then
           do j=1,nyj
! get left halo update
              delta=r2d4b(is,j)-bdy_left(is,j)
              do i=1,is-1
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo
        else
        endif
!
        bdy_left=r2d4b

     deallocate(r2d4b)

end subroutine get_bc_left

subroutine get_bc_right(nxi,nyj,nlon,nlat,fill_halo,&
                         i_right,j_right,bdy_right,field2d,k)
!
!  process top boundary 
!
  use kinds, only: r_kind,i_kind
  implicit none
!
  integer,intent(in) :: nlon,nlat
  integer,intent(in) :: nxi,nyj
  integer,intent(in) :: i_right(nxi),j_right(nyj)
  integer,intent(in) :: fill_halo
  real,intent(in)    :: field2d(nlon,nlat)
  real,intent(inout) :: bdy_right(nxi,nyj)
  integer,intent(in) :: k
!
  real,allocatable :: r2d4b(:,:)
  integer :: is,ie,js,je

  real :: delta
  integer :: i,j
  integer :: ii,jj

     js=1
     je=nyj
     ie=nxi
     do ii=1,nxi
        i=i_right(ii)
        if(i==nlon) is=ii
     enddo
     if(k==1) write(6,'(a50,10I5)') 'update right is,ie,js,je,nxi,nyj:',is,ie,js,je,nxi,nyj

     allocate(r2d4b(nxi,nyj))
        r2d4b=bdy_right
        do jj=1,nyj
           do ii=1,nxi
              i=i_right(ii)
              j=j_right(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=field2d(i,j)
              endif
           enddo
        enddo
!
        if(fill_halo==1) then
           do j=1,nyj
              do i=is+1,ie
                 r2d4b(i,j)=r2d4b(is,j)
              enddo
           enddo
        else if(fill_halo==0) then
           do j=1,nyj
! get right halo update
              delta=r2d4b(is,j)-bdy_right(is,j)
              do i=is+1,ie
                 r2d4b(i,j)=r2d4b(i,j)+delta
              enddo
           enddo
        else
        endif
!
        bdy_right=r2d4b

     deallocate(r2d4b)

end subroutine get_bc_right

end module module_update_bc
