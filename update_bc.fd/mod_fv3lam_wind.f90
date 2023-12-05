module mod_fv3lam_wind
!$$$ module documentation block
!           .      .    .                                       .
! module:   mod_fv3_lola
!   prgmmr: parrish
!
! abstract:  This module contains routines to interpolate from a single
!             fv3 D grid tile to a rotated lat-lon analysis grid which completely
!             covers the fv3 tile.  Points beyond the fv3 tile are
!             filled with nearest fv3 edge values, but have no actual
!             impact on the analysis.
!
! program history log:
!   2017-02-24  parrish--initial documentation (patterned after
!   mod_fv3_to_a.f90)
!   2017-10-10  wu w - setup interpolation and trnsform coeff in generate_anl_grid
!                      add routines earthuv2fv3, fv3uv2earth, fv3_h_to_ll
!                        fv3_ll_to_h
!   2019-11-01  wu   - add checks in generate_anl_grid to present the mean
!                      longitude correctly to fix problem near lon=0
!
! subroutines included:
!   sub generate_anl_grid
!   sub earthuv2fv3
!   sub fv3uv2earth
!   sub fv3_h_to_ll
!   sub fv3_ll_to_h
!   sub rotate2deg 
!   sub unrotate2deg 
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

!       DIAGRAM:  D-Grid layout:
!
!   1                  nx
!   .                   .   (U,H)
!
! 1                     nx +1
! .                       .    (V)

!   U   U   U   U   U   U        + ny +1 (for U)
! V H V H V H V H V H V H V      + ny    (for V,H)
!   U   U   U   U   U   U                            xh(i) = i            dx=1
! V H V H V H V H V H V H V                          xu(i) = i
!   U   U   U   U   U   U                            xv(i) = i-0.5
! V H V H V H V H V H V H V
!   U   U   U   U   U   U                            yh(j) = j            dy=1
! V H V H V H V H V H V H V                          yu(j) = j-0.5
!   U   U   U   U   U   U                            yv(j) = j
! V H V H V H V H V H V H V
!   U   U   U   U   U   U
! V H V H V H V H V H V H V      + 1     (for V,H)
!   U   U   U   U   U   U        + 1     (for U)

! U(nx ,ny +1),V(nx +1,ny ),H(nx ,ny )

  use kinds, only: r_kind,i_kind
  implicit none
!
  private
  public :: initial_wind_convert,fv3uv2earth,earthuv2fv3
  public :: reverse_grid_r_uv,reverse_grid_r_uv_earth
  public :: reverse_grid_r
  public :: cangu,sangu,cangv,sangv,nx,ny

  integer(i_kind) nx,ny
  real(r_kind) ,allocatable,dimension(:,:):: cangu,sangu,cangv,sangv
  

contains

subroutine initial_wind_convert(nx,ny,grid_lon,grid_lat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    generate_anl_grid
!   prgmmr: parrish
!
! abstract:  define rotated lat-lon analysis grid which is centered on fv3 tile 
!             and oriented to completely cover the tile.
!
! program history log:
!   2022-03-28  Hu   based on gsi
!
!   input argument list:
!    nx, ny               - number of cells = nx*ny 
!    grid_lon ,grid_lat   - longitudes and latitudes of fv3 grid cell corners
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds, only: r_kind,i_kind
  use constants, only: one,half
  implicit none

  integer(i_kind), intent(in) :: nx,ny                 ! fv3 tile x- and y-dimensions
  real(r_kind)   , intent(in) :: grid_lon(nx+1,ny+1)   ! fv3 cell corner longitudes
  real(r_kind)   , intent(in) :: grid_lat(nx+1,ny+1)   ! fv3 cell corner latitudes

  integer(i_kind) i,j
  real(r_kind) x(nx+1,ny+1),y(nx+1,ny+1),z(nx+1,ny+1), xr,yr,zr,xu,yu,zu,rlat,rlon
  real(r_kind) xv,yv,zv,vval
  real(r_kind) uval,ewval,nsval
  real(r_kind) diff,sq180
  real(r_kind) pi,deg2rad,rad2deg
!

  pi      = acos(-one)
  deg2rad = pi/180.0_r_kind
  rad2deg = one/deg2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! find coefficients for wind conversion btw FV3 & earth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(cangu)) deallocate(cangu)
  if(allocated(sangu)) deallocate(sangu)
  if(allocated(cangv)) deallocate(cangv)
  if(allocated(sangv)) deallocate(sangv)
  allocate  (   cangu(nx,ny+1) )
  allocate  (   sangu(nx,ny+1) )
  allocate  (   cangv(nx+1,ny) )
  allocate  (   sangv(nx+1,ny) )

!   1.  compute x,y,z at cell cornor from grid_lon, grid_lat

  do j=1,ny+1
     do i=1,nx+1
        x(i,j)=cos(grid_lat(i,j)*deg2rad)*cos(grid_lon(i,j)*deg2rad)
        y(i,j)=cos(grid_lat(i,j)*deg2rad)*sin(grid_lon(i,j)*deg2rad)
        z(i,j)=sin(grid_lat(i,j)*deg2rad)
     enddo
  enddo

!  2   find angles to E-W and N-S for U edges
  sq180=180._r_kind**2 
  do j=1,ny+1
     do i=1,nx
!      center lat/lon of the edge 
        rlat=half*(grid_lat(i,j)+grid_lat(i+1,j))
        diff=(grid_lon(i,j)-grid_lon(i+1,j))**2
        if(diff < sq180)then
           rlon=half*(grid_lon(i,j)+grid_lon(i+1,j))
        else
           rlon=half*(grid_lon(i,j)+grid_lon(i+1,j)-360._r_kind)
        endif
!    vector to center of the edge
        xr=cos(rlat*deg2rad)*cos(rlon*deg2rad)
        yr=cos(rlat*deg2rad)*sin(rlon*deg2rad)
        zr=sin(rlat*deg2rad)
!     vector of the edge
        xu= x(i+1,j)-x(i,j)
        yu= y(i+1,j)-y(i,j)
        zu= z(i+1,j)-z(i,j)
!    find angle with cross product
        uval=sqrt((xu**2+yu**2+zu**2))
        ewval=sqrt((xr**2+yr**2))
        nsval=sqrt((xr*zr)**2+(zr*yr)**2+(xr*xr+yr*yr)**2)
        cangu(i,j)=(-yr*xu+xr*yu)/ewval/uval
        sangu(i,j)=(-xr*zr*xu-zr*yr*yu+(xr*xr+yr*yr)*zu) / nsval/uval
     enddo
  enddo
 
!  3   find angles to E-W and N-S for V edges
  do j=1,ny
     do i=1,nx+1
        rlat=half*(grid_lat(i,j)+grid_lat(i,j+1))
        diff=(grid_lon(i,j)-grid_lon(i,j+1))**2
        if(diff < sq180)then
           rlon=half*(grid_lon(i,j)+grid_lon(i,j+1))
        else
           rlon=half*(grid_lon(i,j)+grid_lon(i,j+1)-360._r_kind)
        endif
        xr=cos(rlat*deg2rad)*cos(rlon*deg2rad)
        yr=cos(rlat*deg2rad)*sin(rlon*deg2rad)
        zr=sin(rlat*deg2rad)
        xv= x(i,j+1)-x(i,j)
        yv= y(i,j+1)-y(i,j)
        zv= z(i,j+1)-z(i,j)
        vval=sqrt((xv**2+yv**2+zv**2))
        ewval=sqrt((xr**2+yr**2))
        nsval=sqrt((xr*zr)**2+(zr*yr)**2+(xr*xr+yr*yr)**2)
        cangv(i,j)=(-yr*xv+xr*yv)/ewval/vval
        sangv(i,j)=(-xr*zr*xv-zr*yr*yv+(xr*xr+yr*yr)*zv) / nsval/vval
     enddo
  enddo

end subroutine initial_wind_convert

subroutine earthuv2fv3(u,v,nx,ny,u_out,v_out)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    earthuv2fv3
!   prgmmr: wu                      2017-06-15
!
! abstract: project earth UV to fv3 UV and interpolate to edge of the cell
!
! program history log:
!   
!
!   input argument list:
!    u,v -  earth wind components at center of the cell
!    nx,ny - dimensions
!
!   output argument list:
!    u_out,v_out - output fv3 winds on the cell boundaries
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block  

  use kinds, only: r_kind,i_kind
  use constants, only: half
  implicit none

  integer(i_kind), intent(in   ) :: nx,ny                 ! fv3 tile x- and y-dimensions
  real(r_kind),intent(in   ) :: u(nx,ny),v(nx,ny)
  real(r_kind),intent(  out) :: u_out(nx,ny+1),v_out(nx+1,ny)
  integer(i_kind) i,j


!!!!!!! earth u/v to covariant u/v
  j=1
  do i=1,nx
     u_out(i,j)= u(i,j)*cangu(i,j)+v(i,j)*sangu(i,j)
  end do

  do j=2,ny
     do i=1,nx
        u_out(i,j)=half   *( (u(i,j)+u(i,j-1))*cangu(i,j)+(v(i,j)+v(i,j-1))*sangu(i,j) )
     end do
  end do
  j=ny
  do i=1,nx
     u_out(i,j+1)= u(i,j)*cangu(i,j+1)+v(i,j)*sangu(i,j+1)
  end do

  do j=1,ny
     v_out(1,j)=u(1,j)*cangv(1,j)+v(1,j)*sangv(1,j)
     do i=2,nx
        v_out(i,j)=half   *( (u(i,j)+u(i-1,j))*cangv(i,j)+(v(i,j)+v(i-1,j))*sangv(i,j) )
     end do
     v_out(nx+1,j)=u(nx,j)*cangv(nx+1,j)+v(nx,j)*sangv(nx+1,j)
  end do
end subroutine earthuv2fv3

subroutine fv3uv2earth(u,v,nx,ny,u_out,v_out)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fv3uv2earth
!   prgmmr: wu                      2017-06-15
!
! abstract: project fv3 UV to earth UV and interpolate to the center of the cells
!
! program history log:
!   
!
!   input argument list:
!    u,v - fv3 winds on the cell boundaries
!    nx,ny - dimensions
!
!   output argument list:
!    u_out,v_out - output earth wind components at center of the cell
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block  

  use kinds, only: r_kind,i_kind
  use constants, only: half
  implicit none

  integer(i_kind), intent(in   ) :: nx,ny                 ! fv3 tile x- and y-dimensions
  real(r_kind),intent(in   ) :: u(nx,ny+1),v(nx+1,ny)
  real(r_kind),intent(  out) :: u_out(nx,ny),v_out(nx,ny)
  integer(i_kind) i,j

  do j=1,ny
     do i=1,nx
        u_out(i,j)=half *( (u(i,j)*sangv(i,j)-v(i,j)*sangu(i,j))/(cangu(i,j)*sangv(i,j)-sangu(i,j)*cangv(i,j)) &
                       +(u(i,j+1)*sangv(i+1,j)-v(i+1,j)*sangu(i,j+1))/(cangu(i,j+1)*sangv(i+1,j)-sangu(i,j+1)*cangv(i+1,j)))
        v_out(i,j)=half *( (u(i,j)*cangv(i,j)-v(i,j)*cangu(i,j))/(sangu(i,j)*cangv(i,j)-cangu(i,j)*sangv(i,j)) &
                       +(u(i,j+1)*cangv(i+1,j)-v(i+1,j)*cangu(i,j+1))/(sangu(i,j+1)*cangv(i+1,j)-cangu(i,j+1)*sangv(i+1,j)))
     end do
  end do
  return
end subroutine fv3uv2earth

subroutine reverse_grid_r_uv(grid,nx,ny,nz)
!
!  reverse the first two dimension of the array grid
!
    use kinds, only: r_kind,i_kind

    implicit none
    integer(i_kind), intent(in     ) :: nx,ny,nz
    real(r_kind),    intent(inout  ) :: grid(nx,ny,nz)
    real(r_kind)                     :: tmp_grid(nx,ny)
    integer(i_kind)                  :: i,j,k
!
    do k=1,nz
       tmp_grid(:,:)=grid(:,:,k)
       do j=1,ny
          do i=1,nx
             grid(i,j,k)=-tmp_grid(nx+1-i,ny+1-j)
          enddo
       enddo
    enddo

end subroutine reverse_grid_r_uv

subroutine reverse_grid_r_uv_earth(grid,nx,ny,nz)
!
!  reverse the first two dimension of the array grid
!
    use kinds, only: r_kind,i_kind

    implicit none
    integer(i_kind), intent(in     ) :: nx,ny,nz
    real(r_kind),    intent(inout  ) :: grid(nx,ny,nz)
    real(r_kind)                     :: tmp_grid(nx,ny)
    integer(i_kind)                  :: i,j,k
!
    do k=1,nz
       tmp_grid(:,:)=grid(:,:,k)
       do j=1,ny
          do i=1,nx
             grid(i,j,k)=tmp_grid(nx+1-i,ny+1-j)
          enddo
       enddo
    enddo

end subroutine reverse_grid_r_uv_earth


subroutine reverse_grid_r(grid,nx,ny,nz)
!
!  reverse the first two dimension of the array grid
!
    use kinds, only: r_kind,i_kind

    implicit none
    integer(i_kind),  intent(in     ) :: nx,ny,nz
    real(r_kind),     intent(inout  ) :: grid(nx,ny,nz)
    real(r_kind)                      :: tmp_grid(nx,ny)
    integer(i_kind)                   :: i,j,k
!
    do k=1,nz
       tmp_grid(:,:)=grid(:,:,k)
       do j=1,ny
          do i=1,nx
             grid(i,j,k)=tmp_grid(nx+1-i,ny+1-j)
          enddo
       enddo
    enddo

end subroutine reverse_grid_r

end module mod_fv3lam_wind
