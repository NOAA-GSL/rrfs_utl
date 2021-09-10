program use_raphrrr_sfc
!
  use kinds, only: i_kind,r_kind,r_single,i_byte
  use module_ncio, only : ncio
  use module_map_utils, only : map_util
  use module_surface, only : use_surface 
  use mpi
!
  implicit none
!
  type(ncio)     :: raphrrr,rrfs
  type(map_util) :: map
  type(use_surface) :: sfc
!
!  namelist files
!
  character*80 :: rapfile
  character*80 :: hrrrfile
  character*80 :: hrrr_akfile
  character*80 :: rrfsfile
  namelist/setup/ rapfile,hrrrfile,hrrr_akfile,rrfsfile
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! define grid
  integer :: nx_rrfs,ny_rrfs
  integer :: nx_rap,ny_rap
  integer :: nx_hrrr,ny_hrrr
  integer :: nx_hrrrak,ny_hrrrak
  integer :: nz_rrfs,nz_raphrrr
!
! define map
  real(r_single),allocatable,target :: rlon2d_rrfs(:,:),rlat2d_rrfs(:,:)
  real(r_single),allocatable,target :: rlon2d_raphrrr(:,:),rlat2d_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_rrfs(:,:)
  real(r_single),allocatable,target :: tmp2d4b(:,:)
  real(r_kind),  allocatable,target :: tmp2d8b(:,:)

  integer :: i,j,k
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
  if(mype==0) then
!
     rapfile='missing'
     hrrrfile='missing'
     hrrr_akfile='missing'
     rrfsfile='missing'
     open(15, file='use_raphrrr_sfc.namelist')
        read(15,setup)
     close(15)
     write(*,setup)

! read in rrfs latlon
     call rrfs%open('fv3_grid_spec',"r",200)
     call rrfs%get_dim("grid_xt",nx_rrfs)
     call rrfs%get_dim("grid_yt",ny_rrfs)
     write(*,*) 'nx_rrfs,ny_rrfs=',nx_rrfs,ny_rrfs

     allocate(rlon2d_rrfs(nx_rrfs,ny_rrfs))
     allocate(rlat2d_rrfs(nx_rrfs,ny_rrfs))
     call rrfs%get_var("grid_lont",nx_rrfs,ny_rrfs,rlon2d_rrfs)
     call rrfs%get_var("grid_latt",nx_rrfs,ny_rrfs,rlat2d_rrfs)
     call rrfs%close()
!
! read in rrfs land mask
     call rrfs%open(trim(rrfsfile),"r",200)
     allocate(landmask_rrfs(nx_rrfs,ny_rrfs))
     allocate(tmp2d8b(nx_rrfs,ny_rrfs))
     call rrfs%get_var("slmsk",nx_rrfs,ny_rrfs,tmp2d8b)
     do j=1,ny_rrfs
       do i=1,nx_rrfs
          landmask_rrfs(i,j)=int(tmp2d8b(i,j))
          if(landmask_rrfs(i,j) >=2 ) landmask_rrfs(i,j)=1
       enddo
     enddo
     deallocate(tmp2d8b)

     call rrfs%get_dim("zaxis_1",nz_rrfs)
     call rrfs%close()
     write(*,*) maxval(landmask_rrfs), minval(landmask_rrfs)
!
! read in rap latlon
     call raphrrr%open(trim(rapfile),"r",200)
     call raphrrr%get_dim("west_east",nx_rap)
     call raphrrr%get_dim("south_north",ny_rap)
     call raphrrr%get_dim("soil_layers_stag",nz_raphrrr)
     write(*,*) 'nx_rap,ny_rap=',nx_rap,ny_rap,nz_raphrrr
     if(nz_raphrrr /= nz_rrfs) then
        write(*,*) "Error in vertical level =", nz_raphrrr,nz_rrfs
        stop 123
     endif
!
     allocate(rlon2d_raphrrr(nx_rap,ny_rap))
     allocate(rlat2d_raphrrr(nx_rap,ny_rap))
     call raphrrr%get_var("XLONG",nx_rap,ny_rap,rlon2d_raphrrr)
     call raphrrr%get_var("XLAT",nx_rap,ny_rap,rlat2d_raphrrr)
     call map%init_general_transform(nx_rap,ny_rap,rlat2d_raphrrr,rlon2d_raphrrr)
     deallocate(rlon2d_raphrrr)
     deallocate(rlat2d_raphrrr)

     allocate(landmask_raphrrr(nx_rap,ny_rap))
     allocate(tmp2d4b(nx_rap,ny_rap))
     call raphrrr%get_var("LANDMASK",nx_rap,ny_rap,tmp2d4b)
     do j=1,ny_rap
       do i=1,nx_rap
          landmask_raphrrr(i,j)=int(tmp2d4b(i,j))
       enddo
     enddo
     call raphrrr%get_var("SNOW",nx_rap,ny_rap,tmp2d4b)
     do j=1,ny_rap
       do i=1,nx_rap
          if(tmp2d4b(i,j) > 0.01 ) landmask_raphrrr(i,j)=2  ! snow coverage
       enddo
     enddo
     deallocate(tmp2d4b)
     call raphrrr%close()
!
! initial sfc and map index
     call sfc%init(nx_rrfs,ny_rrfs,nz_rrfs,nx_rap,ny_rap,4)
     call sfc%build_mapindex(map,rlon2d_rrfs,rlat2d_rrfs,landmask_rrfs,landmask_raphrrr)
     call sfc%set_varname()
     call sfc%use_sfc(rapfile,rrfsfile)
     call sfc%close()

! read in hrrr latlon
     call raphrrr%open(trim(hrrrfile),"r",200)
     call raphrrr%get_dim("west_east",nx_hrrr)
     call raphrrr%get_dim("south_north",ny_hrrr)
     write(*,*) 'nx_hrrr,ny_hrrr=',nx_hrrr,ny_hrrr
     call raphrrr%close()

! read in hrrr ak latlon
     call raphrrr%open(trim(hrrr_akfile),"r",200)
     call raphrrr%get_dim("west_east",nx_hrrrak)
     call raphrrr%get_dim("south_north",ny_hrrrak)
     write(*,*) 'nx_hrrrak,ny_hrrrak=',nx_hrrrak,ny_hrrrak
     call raphrrr%close()

     write(6,*) "=== USE_RAPHRRR_SFC REPROCCESS SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!
end program use_raphrrr_sfc
