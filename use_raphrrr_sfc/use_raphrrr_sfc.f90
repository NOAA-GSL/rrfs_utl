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
  character*80 :: rrfsfile_read
  logical :: do_lake_surgery
  namelist/setup/ rapfile,hrrrfile,hrrr_akfile,rrfsfile,do_lake_surgery
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
  integer :: nz_rrfs_lake,nz_raphrrr_lake
  integer :: nz_rrfs_snow,nz_raphrrr_snow
!
! define map
  real(r_single),allocatable,target :: rlon2d_rrfs(:,:),rlat2d_rrfs(:,:)
  real(r_single),allocatable,target :: rlon2d_raphrrr(:,:),rlat2d_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_rrfs(:,:)
  integer(i_byte),allocatable :: lakemask_raphrrr(:,:)
  integer(i_byte),allocatable :: lakemask_rrfs(:,:)
  real(r_single),allocatable,target :: tmp2d4b(:,:)

  integer :: i,j,k,n
  character*80 :: raphrrrfile
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
     do_lake_surgery=.false.
     open(15, file='use_raphrrr_sfc.namelist')
        read(15,setup)
     close(15)
     write(*,setup)

     rrfsfile_read=trim(rrfsfile)//"_read"
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
     allocate(lakemask_rrfs(nx_rrfs,ny_rrfs))
     allocate(tmp2d4b(nx_rrfs,ny_rrfs))
     call rrfs%get_var("slmsk",nx_rrfs,ny_rrfs,tmp2d4b)
     do j=1,ny_rrfs
       do i=1,nx_rrfs
          landmask_rrfs(i,j)=int(tmp2d4b(i,j))
          if(landmask_rrfs(i,j) >=2 ) landmask_rrfs(i,j)=1
       enddo
     enddo
     if(do_lake_surgery) then
       call rrfs%get_var("clm_lake_initialized",nx_rrfs,ny_rrfs,tmp2d4b)
       do j=1,ny_rrfs
         do i=1,nx_rrfs
            lakemask_rrfs(i,j)=int(tmp2d4b(i,j))
         enddo
       enddo
     else
       lakemask_rrfs(:,:)=0
     endif
     deallocate(tmp2d4b)

     call rrfs%get_dim("zaxis_1",nz_rrfs)
     if(do_lake_surgery) then
        call rrfs%get_dim("levlake_clm_lake",nz_rrfs_lake)
        call rrfs%get_dim("levsnowsoil1_clm_lake",nz_rrfs_snow)
     else
        nz_rrfs_lake=-99  
        nz_rrfs_snow=-99
     endif
     call rrfs%close()
     write(*,*) "landmask=",maxval(landmask_rrfs), minval(landmask_rrfs)
     write(*,*) "lakemask=",maxval(lakemask_rrfs), minval(lakemask_rrfs)
     write(*,*) "nz_rrfs=",nz_rrfs
     write(*,*) "nz_rrfs_lake=",nz_rrfs_lake
     write(*,*) "nz_rrfs_snow=",nz_rrfs_snow
!
! use RAP
! read in rap latlon
     do n=1,3

        raphrrrfile='missing'
        if(n==1 .and. trim(rapfile) /= 'missing' ) then
           raphrrrfile=rapfile
        elseif(n==2 .and. trim(hrrrfile) /= 'missing' ) then
           raphrrrfile=hrrrfile
        elseif(n==3 .and. trim(hrrr_akfile) /= 'missing' ) then
           raphrrrfile=hrrr_akfile
        else
           cycle
        endif
        write(*,*) "===>"
        write(*,*) "===> tansfer surface fields from ",trim(raphrrrfile)
        write(*,*) "===>"
        call raphrrr%open(trim(raphrrrfile),"r",200)
        call raphrrr%get_dim("west_east",nx_rap)
        call raphrrr%get_dim("south_north",ny_rap)
        call raphrrr%get_dim("soil_layers_stag",nz_raphrrr)
        if(do_lake_surgery) then
           call raphrrr%get_dim("soil_levels_or_lake_levels_stag",nz_raphrrr_lake)
           call raphrrr%get_dim("snow_and_soil_levels_stag",nz_raphrrr_snow)
        else
           nz_raphrrr_lake=-99
           nz_raphrrr_snow=-99
        endif
        write(*,*) 'nx_rap,ny_rap=',nx_rap,ny_rap,nz_raphrrr
        write(*,*) 'nz_raphrrr_lake,nz_raphrrr_snow=',nz_raphrrr_lake,nz_raphrrr_snow
        if(nz_raphrrr /= nz_rrfs .or. nz_raphrrr_lake /= nz_rrfs_lake .or. &
           nz_raphrrr_snow /= nz_rrfs_snow) then
           write(*,*) "Error in vertical level =", nz_raphrrr,nz_rrfs
           write(*,*) "Error in vertical level lake=", nz_raphrrr_lake,nz_rrfs_lake
           write(*,*) "Error in vertical level snow=", nz_raphrrr_snow,nz_rrfs_snow
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
        allocate(lakemask_raphrrr(nx_rap,ny_rap))
        allocate(tmp2d4b(nx_rap,ny_rap))
        call raphrrr%get_var("LANDMASK",nx_rap,ny_rap,tmp2d4b)
        do j=1,ny_rap
          do i=1,nx_rap
             landmask_raphrrr(i,j)=int(tmp2d4b(i,j))
          enddo
        enddo
        call raphrrr%get_var("LAKEMASK",nx_rap,ny_rap,tmp2d4b)
        do j=1,ny_rap
          do i=1,nx_rap
             lakemask_raphrrr(i,j)=int(tmp2d4b(i,j))
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
        call sfc%init(nx_rrfs,ny_rrfs,nz_rrfs,nz_rrfs_lake,nz_rrfs_snow,nx_rap,ny_rap,4)
        call sfc%build_mapindex(map,rlon2d_rrfs,rlat2d_rrfs,landmask_rrfs,landmask_raphrrr)
        if(do_lake_surgery) &
           call sfc%build_lakeindex(map,rlon2d_rrfs,rlat2d_rrfs,lakemask_rrfs,lakemask_raphrrr)
        call sfc%set_varname()
        call sfc%use_sfc(raphrrrfile,rrfsfile,rrfsfile_read)
        if(do_lake_surgery) &
           call sfc%use_lake(raphrrrfile,rrfsfile,rrfsfile_read)
        if(n==2) call sfc%remove_snow(raphrrrfile,rrfsfile,rlat2d_rrfs)
        call sfc%close()
! release memory
        deallocate(landmask_raphrrr)
        deallocate(lakemask_raphrrr)
        call map%destory_general_transform()

     enddo ! n
!
     deallocate(landmask_rrfs)
     deallocate(lakemask_rrfs)
!
! use HRRR
! read in hrrr latlon
!     call raphrrr%open(trim(hrrrfile),"r",200)
!     call raphrrr%get_dim("west_east",nx_hrrr)
!     call raphrrr%get_dim("south_north",ny_hrrr)
!!     call raphrrr%get_dim("soil_layers_stag",nz_raphrrr)
!     write(*,*) 'nx_hrrr,ny_hrrr=',nx_hrrr,ny_hrrr,nz_raphrrr
!     if(nz_raphrrr /= nz_rrfs) then
!        write(*,*) "Error in vertical level =", nz_raphrrr,nz_rrfs
!        stop 123
!     endif
!!     call raphrrr%close()

! read in hrrr ak latlon
!     call raphrrr%open(trim(hrrr_akfile),"r",200)
!     call raphrrr%get_dim("west_east",nx_hrrrak)
!     call raphrrr%get_dim("south_north",ny_hrrrak)
!     call raphrrr%get_dim("soil_layers_stag",nz_raphrrr)
!     write(*,*) 'nx_hrrrak,ny_hrrrak=',nx_hrrrak,ny_hrrrak,nz_raphrrr
!     if(nz_raphrrr /= nz_rrfs) then
!!        write(*,*) "Error in vertical level =", nz_raphrrr,nz_rrfs
!        stop 123
!     endif
!     call raphrrr%close()
!
     write(6,*) "=== USE_RAPHRRR_SFC REPROCCESS SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!
end program use_raphrrr_sfc
