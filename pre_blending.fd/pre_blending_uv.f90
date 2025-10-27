PROGRAM pre_blending
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2025-05-20
!
! ABSTRACT: 
!     This appllication calculates ensemble mean and recenter
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!
! REMARKS:
!
! ATTRIBUTES:
!
!$$$
!
!_____________________________________________________________________

  use mpi
  use netcdf 
  use module_ncfile_stat, only : ncfile_stat
  use module_mpi_arrange, only : mpi_io_arrange
  use chgres_winds_mod
  use remap_scalar_mod
  use remap_dwinds_mod

  use general_sub2grid_simple_mod, only: general_sub2grid_create_info
  use general_sub2grid_simple_mod, only: general_sub2grid_destroy_info
  use general_sub2grid_simple_mod, only: general_grid2sub
  use general_sub2grid_simple_mod, only: general_sub2grid
  use general_sub2grid_simple_mod, only: sub2grid_info


!  use netcdf, only: nf90_open,nf90_close,nf90_noerr
!  use netcdf, only: nf90_put_var
!  use netcdf, only: nf90_nowrite,nf90_write
!  use netcdf, only: nf90_inq_varid
!  use netcdf, only: nf90_strerror

  implicit none
! 
  type(ncfile_stat) :: ncfs_all
  type(mpi_io_arrange) :: mpiioarg
  type(sub2grid_info) :: s
!
! MPI variables
  integer :: npe, mype, ierror
! sub communicator
  integer :: new_comm
  integer :: color, key
!
! namelist
!
  integer, parameter    :: filename_len=100
  integer             :: numvar(2)
  character(len=200)  :: varlist(2)

  character (len=filename_len)   :: filecold(2)
  character (len=filename_len)   :: akbk,akbk_cold

  logical :: ifexist
  integer :: ncioid
!
! MPI distribution array
  integer :: mype_fileid
  character(len=20) :: mype_varname,local_varname
  integer :: mype_vartype
  integer :: mype_nx,mype_ny
  integer :: mype_lbegin,mype_lend
  integer,allocatable,dimension(:) :: lvl2core

  integer :: lon2, lat2, nsig
  integer :: mype_istart,mype_jstart

!
! array for wind
 real(kind=8),allocatable, dimension(:,:,:)     :: u_s,v_s     !
 real(kind=8),allocatable, dimension(:,:,:)     :: u_w,v_w     !
 real(kind=8),allocatable, dimension(:,:,:)     :: ud_local    !(3950, 2701, 66)
 real(kind=8),allocatable, dimension(:,:,:)     :: vd_local    !(3951, 2700, 66)
 real(kind=8),allocatable, dimension(:,:)       :: gridx,gridy !(nlon,nlat)

 real(kind=8),allocatable, dimension(:)         ::  ak0      !(67)
 real(kind=8),allocatable, dimension(:)         ::  bk0      !(67)
 real(kind=8),allocatable, dimension(:,:)       ::  psc      !(3950, 2700)
 real(kind=8),allocatable, dimension(:,:)       ::  Atm_ps   !(3950, 2700)
 real(kind=8),allocatable, dimension(:,:,:)     ::  ud       !(3950, 2701, 66)
 real(kind=8),allocatable, dimension(:,:,:)     ::  vd       !(3951, 2700, 66)
 real(kind=8),allocatable, dimension(:)         ::  Atm_ak   !(66)
 real(kind=8),allocatable, dimension(:)         ::  Atm_bk   !(66)
 real(kind=8),allocatable, dimension(:,:,:)     ::  Atm_u    !(3950, 2701, 65)
 real(kind=8),allocatable, dimension(:,:,:)     ::  Atm_v    !(3951, 2700, 65)

 real(kind=8), allocatable :: ps_local(:,:), zh_local(:,:,:)
 real(kind=8), allocatable :: omga_local(:,:,:), delp_local(:,:,:)
 real(kind=8), allocatable :: t_local(:,:,:), qa_local(:,:,:,:)
 real(kind=8), allocatable :: Atm_phis_local(:,:)

 real(kind=8), allocatable :: Atm_delp(:,:,:)
 real(kind=8), allocatable :: Atm_pt(:,:,:), Atm_q(:,:,:,:)

 real,allocatable, dimension(:,:)     ::  d2r4
 real,allocatable, dimension(:,:,:)   ::  d3r4
 real,allocatable, dimension(:,:,:)   ::  sub_vars
 !
 real,allocatable, dimension(:,:,:)   ::  d3r4_us,d3r4_vs,d3r4_uw
 integer :: us_1st,vs_1st,uw_1st,vw_1st

!
  integer :: n,i,j,k,iret,ilev,ierr
  integer :: grid_cdfid,cdfid,oldMode
! dimensions
  TYPE(MPI_Status) :: istatus
  integer :: nlev, nlat, nlon, nlatp, nlonp, nlevp
  integer :: dimid_lev, dimid_lat, dimid_lon, dimid_latp, dimid_lonp, dimid_nlev
  integer :: local_nx,local_ny
! Variable IDs and I/O arrays
  integer :: varid, dimids(4), start(4), count(4), chunksizes(4)
  integer :: nlevid
  integer :: ntotalcore,num_fields
  integer,allocatable :: kbegin(:),kend(:)
  logical :: create_new_file
  character(len=20),allocatable :: varname(:)

!
!**********************************************************************
!**********************************************************************
!
!            END OF DECLARATIONS....start of program

! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)
!
  filecold(1)='out.atm.tile7.nc'
!
  call check(nf90_open(trim(filecold(1)), IOR(NF90_NOWRITE, NF90_MPIIO), cdfid, &
                       comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))
    call check(nf90_inq_dimid(cdfid, "lev", dimid_lev))
    call check(nf90_inquire_dimension(cdfid, dimid_lev, len=nlev))    ! 66
    call check(nf90_inq_dimid(cdfid, "lat", dimid_lat))
    call check(nf90_inquire_dimension(cdfid, dimid_lat, len=nlat))    ! 2700
    call check(nf90_inq_dimid(cdfid, "lon", dimid_lon))
    call check(nf90_inquire_dimension(cdfid, dimid_lon, len=nlon))    ! 3950
    call check(nf90_inq_dimid(cdfid, "latp", dimid_latp))
    call check(nf90_inquire_dimension(cdfid, dimid_latp, len=nlatp))  ! 2701
    call check(nf90_inq_dimid(cdfid, "lonp", dimid_lonp))
    call check(nf90_inquire_dimension(cdfid, dimid_lonp, len=nlonp))  ! 3951
    call check(nf90_inq_dimid(cdfid, "levp", dimid_nlev))
    call check(nf90_inquire_dimension(cdfid, dimid_nlev, len=nlevp))  ! 67
  call check(nf90_close(cdfid))

  if (mype == 0) then
     write(*,*) nlat,nlon,nlatp,nlonp,nlev,nlevp
     akbk_cold="gfs_ctrl.nc"
     call check(nf90_open(trim(akbk_cold), NF90_NOWRITE, grid_cdfid))
     allocate(ak0(nlevp), bk0(nlevp))
     call check(nf90_inq_varid(grid_cdfid, "vcoord", varid))
     start(1:2) = [1, 1]
     count(1:2) = [nlevp,1]
     call check(nf90_get_var(grid_cdfid, varid, ak0, start=start(1:2), count=count(1:2)))
     start(1:2) = [1, 2]
     call check(nf90_get_var(grid_cdfid, varid, bk0, start=start(1:2), count=count(1:2)))
     call check(nf90_close(grid_cdfid))

     akbk="fv_core.res.nc"
     call check(nf90_open(trim(akbk), NF90_NOWRITE, grid_cdfid))
     allocate(Atm_ak(nlev), Atm_bk(nlev))
     call check(nf90_inq_varid(grid_cdfid, "ak", varid))
     call check(nf90_get_var(grid_cdfid, varid, Atm_ak))
     call check(nf90_inq_varid(grid_cdfid, "bk", varid))
     call check(nf90_get_var(grid_cdfid, varid, Atm_bk))
     call check(nf90_close(grid_cdfid))

     Atm_ak(:)=ak0(2:nlevp)
     Atm_bk(:)=bk0(2:nlevp)
     ak0(1)=1.000000000000000E-009_8
     bk0(1)=1.000000000000000E-009_8
  else
     allocate(ak0(nlevp), bk0(nlevp))
     allocate(Atm_ak(nlev), Atm_bk(nlev))
  endif
  call MPI_Bcast(ak0, nlevp, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(bk0, nlevp, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Atm_ak, nlev, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Atm_bk, nlev, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
!
!-------------------------------------------------------------------
! now working on scalars 
!-------------------------------------------------------------------
!
!
! get atm_ps, this will be used for wind
  if(mype==0) then
      open(12,file='Atm_ps.bin',form='unformatted',status='old')
         read(12) mype_nx,mype_ny
         allocate(Atm_ps(mype_nx,mype_ny))
         read(12) Atm_ps
     close(12)
     write(*,*) mype, "ps=",maxval(Atm_ps), minval(Atm_ps)
  endif
!  
!-------------------------------------------------------------------
! now working on wind 
!-------------------------------------------------------------------
!
  filecold(1)='out.atm.tile7.nc'
!
!
  allocate(gridx(nlon,nlat))
  allocate(gridy(nlon,nlat))
  allocate(psc(nlon,nlat))

  allocate(d3r4(nlon,nlat,1))
  if(mype==0) then
     call check(nf90_open(trim(filecold(1)), nf90_nowrite, cdfid))

     start = [1, 1, 1, 1]
     count = [nlon, nlat, 1, 1]

     call check(nf90_inq_varid(cdfid, "ps", varid))
     call check(nf90_get_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
     psc=d3r4(:,:,1)

     call check(nf90_inq_varid(cdfid, "geolon", varid))
     call check(nf90_get_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
     gridx=d3r4(:,:,1)

     call check(nf90_inq_varid(cdfid, "geolat", varid))
     call check(nf90_get_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
     gridy=d3r4(:,:,1)
     write(*,*) maxval(psc),minval(psc)
     write(*,*) maxval(gridx),minval(gridx)
     write(*,*) maxval(gridy),minval(gridy)

     call check(nf90_close(cdfid))
  endif
  deallocate(d3r4)
  write(*,*)"mype1=",mype, "read ps geolon geolat for U and V"

  call mpi_barrier(MPI_COMM_WORLD,ierror)
  call MPI_Bcast(psc, nlon*nlat, MPI_DOUBLE , 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gridx, nlon*nlat, MPI_DOUBLE , 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gridy, nlon*nlat, MPI_DOUBLE , 0, MPI_COMM_WORLD, ierr)
  if(mype==3) then
     write(*,*) 'psc=',maxval(psc),minval(psc)
     write(*,*) 'gridx=',maxval(gridx),minval(gridx)
     write(*,*) 'gridy=',maxval(gridy),minval(gridy)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierror)

  numvar(1)=4
  varlist(1)='u_s v_s u_w v_w'
  if(mype==0) then
!
! find dimension of each field
!
     call ncfs_all%init(1,filecold, numvar, varlist)
     call ncfs_all%fill_dims()
!
!  distibute variables to each core
!
     call mpiioarg%init(npe)
     call mpiioarg%arrange(ncfs_all)

     call ncfs_all%close()
  endif

  call MPI_Scatter(mpiioarg%fileid, 1, mpi_integer, mype_fileid, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%varname, 20, mpi_character, mype_varname, 20, mpi_character, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%vartype, 1, mpi_integer, mype_vartype, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%nx, 1, mpi_integer, mype_nx, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%ny, 1, mpi_integer, mype_ny, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%lvlbegin, 1, mpi_integer, mype_lbegin, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%lvlend, 1, mpi_integer, mype_lend, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  !
  allocate(lvl2core(100))
  lvl2core=-1
  us_1st=-1
  vs_1st=-1
  uw_1st=-1
  vw_1st=-1
  if(mype==0) then
     do n=1,npe
        do k=mpiioarg%lvlbegin(n),mpiioarg%lvlend(n)
           if(lvl2core(k)==-1) lvl2core(k)=n
        enddo
        if(us_1st==-1 .and. trim(mpiioarg%varname(n))=='u_s') us_1st=n-1
        if(vs_1st==-1 .and. trim(mpiioarg%varname(n))=='v_s') vs_1st=n-1
        if(uw_1st==-1 .and. trim(mpiioarg%varname(n))=='u_w') uw_1st=n-1
        if(vw_1st==-1 .and. trim(mpiioarg%varname(n))=='v_w') vw_1st=n-1
     enddo
  endif
  call MPI_Bcast(lvl2core, 100, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(us_1st, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(vs_1st, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(uw_1st, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Bcast(vw_1st, 1, mpi_integer, 0, MPI_COMM_WORLD,ierror)

  if(mype==0) call mpiioarg%close()
  if(mype==1) then
     write(6,*) "The first core of u_s,v_s,u_w,v_w:",us_1st,vs_1st,uw_1st,vw_1st
     write(6,*) "The core id matchies level:",lvl2core
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)

! Create sub-communicator to handle each file
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= 1) then
     color = mype_fileid
  else
     color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(mpi_comm_world,color,key,new_comm,ierror)
  if ( ierror /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',ierror
     call mpi_abort(mpi_comm_world,101,ierror)
  endif
  if(npe <= nlev ) then
     write(6,'(a,i5)')'***ERROR*** need more than nlev cores to run nlev = ',nlev
     call mpi_abort(mpi_comm_world,102,1)
  endif
!
! read 2D field from each file using sub communicator
!
  allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))
  if (MPI_COMM_NULL /= new_comm) then

     iret=nf90_open(trim(filecold(1)),nf90_nowrite,ncioid,comm=new_comm,info=MPI_INFO_NULL)
     if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(filecold(1)), ' Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           call mpi_abort(new_comm,103,1)
     endif

     call mype_read(ncioid,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)
     do k=mype_lbegin,mype_lend
        write(6,'(a,a10,i10,2f20.6)') 'reading wind=',trim(adjustl(mype_varname)),k,maxval(d3r4(:,:,:)),minval(d3r4(:,:,:))
     enddo

     call check(nf90_close(ncioid))
  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
!  send u_s,v_s,u_w,v_w to the core with v_w. So those 4 fields in the same 
!
  if(mype==0) write(6,*) "start collect u_s, v_s, u_w, v_w in the same level"
  if(trim(mype_varname)=='v_w') then
     allocate(d3r4_us(nlon,nlatp,mype_lbegin:mype_lend))
     allocate(d3r4_vs(nlon,nlatp,mype_lbegin:mype_lend))
     allocate(d3r4_uw(nlonp,nlat,mype_lbegin:mype_lend))
  endif
! send u_s to v_w cores
  if(trim(mype_varname)=='u_s') then
!     write(6,*) 'send u_s', mype, vw_1st+mype,nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Send(d3r4, nlon*nlatp*(mype_lend-mype_lbegin+1), MPI_real,  vw_1st+mype, 0, mpi_comm_world, ierror)
  endif
  if(trim(mype_varname)=='v_w') then
!     write(6,*) 'receive  u_s', mype, mype-vw_1st, nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Recv(d3r4_us, nlon*nlatp*(mype_lend-mype_lbegin+1), MPI_real, mype-vw_1st, 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)

! send v_s to v_w cores
  if(trim(mype_varname)=='v_s') then
!     write(6,*) 'send v_s', mype, vw_1st+(mype-vs_1st),nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Send(d3r4, nlon*nlatp*(mype_lend-mype_lbegin+1), MPI_real,  vw_1st+(mype-vs_1st), 0, mpi_comm_world, ierror)
  endif
  if(trim(mype_varname)=='v_w') then
!     write(6,*) 'receive  v_s', mype, vs_1st+(mype-vw_1st), nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Recv(d3r4_vs, nlon*nlatp*(mype_lend-mype_lbegin+1), MPI_real, vs_1st+(mype-vw_1st), 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)

! send u_s to v_w cores
  if(trim(mype_varname)=='u_w') then
!     write(6,*) 'send u_w', mype, vw_1st+(mype-uw_1st),nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Send(d3r4, nlonp*nlat*(mype_lend-mype_lbegin+1), MPI_real,  vw_1st+(mype-uw_1st), 0, mpi_comm_world, ierror)
  endif
  if(trim(mype_varname)=='v_w') then
!     write(6,*) 'receive  u_w', mype, uw_1st+(mype-vw_1st), nlon*nlatp*(mype_lend-mype_lbegin+1)
     call MPI_Recv(d3r4_uw, nlonp*nlat*(mype_lend-mype_lbegin+1), MPI_real, uw_1st+(mype-vw_1st), 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!  if(trim(mype_varname)=='v_w') then
!     do k=mype_lbegin,mype_lend
!        write(6,'(2I10,8f20.6)') mype,k,maxval(d3r4(:,:,k)),minval(d3r4(:,:,k)),maxval(d3r4_us(:,:,k)),minval(d3r4_us(:,:,k)),&
!                          maxval(d3r4_vs(:,:,k)),minval(d3r4_vs(:,:,k)),maxval(d3r4_uw(:,:,k)),minval(d3r4_uw(:,:,k))
!     enddo
!  endif

  allocate(u_s(nlon,nlatp,1))
  allocate(v_s(nlon,nlatp,1))
  allocate(u_w(nlonp,nlat,1))
  allocate(v_w(nlonp,nlat,1))

! send u_s to each level u_s (1-66 level)
  allocate(d2r4(nlon,nlatp))
  if(trim(mype_varname)=='v_w') then
     do k=mype_lbegin,mype_lend
        call MPI_Send(d3r4_us(:,:,k), nlon*nlatp, MPI_real,  k-1, 0, mpi_comm_world, ierror)
     enddo
  endif
  if(mype < nlev) then
     call MPI_Recv(d2r4, nlon*nlatp, MPI_real, lvl2core(mype+1)+vw_1st-1, 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
     u_s(:,:,1)=d2r4
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  deallocate(d2r4)
! send v_s to each level v_s (1-66 level)
  allocate(d2r4(nlon,nlatp))
  if(trim(mype_varname)=='v_w') then
     do k=mype_lbegin,mype_lend
        call MPI_Send(d3r4_vs(:,:,k), nlon*nlatp, MPI_real,  k-1, 0, mpi_comm_world, ierror)
     enddo
  endif
  if(mype < nlev) then
     call MPI_Recv(d2r4, nlon*nlatp, MPI_real, lvl2core(mype+1)+vw_1st-1, 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
     v_s(:,:,1)=d2r4
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  deallocate(d2r4)
! send u_w to each level u_w
  allocate(d2r4(nlonp,nlat))
  if(trim(mype_varname)=='v_w') then
     do k=mype_lbegin,mype_lend
        call MPI_Send(d3r4_uw(:,:,k), nlonp*nlat, MPI_real,  k-1, 0, mpi_comm_world, ierror)
     enddo
  endif
  if(mype < nlev) then
     call MPI_Recv(d2r4, nlonp*nlat, MPI_real, lvl2core(mype+1)+vw_1st-1, 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
     u_w(:,:,1)=d2r4
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  deallocate(d2r4)
! send v_w to each level v_w
  allocate(d2r4(nlonp,nlat))
  if(trim(mype_varname)=='v_w') then
     do k=mype_lbegin,mype_lend
        call MPI_Send(d3r4(:,:,k), nlonp*nlat, MPI_real,  k-1, 0, mpi_comm_world, ierror)
     enddo
  endif
  if(mype < nlev) then
     call MPI_Recv(d2r4, nlonp*nlat, MPI_real, lvl2core(mype+1)+vw_1st-1, 0, mpi_comm_world,MPI_STATUS_IGNORE,ierror)
     v_w(:,:,1)=d2r4
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  deallocate(d2r4)
  deallocate(d3r4)
  if(trim(mype_varname)=='v_w') then
     deallocate(d3r4_us)
     deallocate(d3r4_vs)
     deallocate(d3r4_uw)
  endif

  if(mype < nlev) then
     if(mype==0) write(6,'(10a15)') 'level/core id', 'u_s max', 'u_s min', 'v_s max', 'v_s min', 'u_w max', 'u_w min', 'v_w max', 'v_w min'
     write(6,'(I15,8f15.6)') mype,maxval(u_s),minval(u_s),maxval(v_s),minval(v_s),maxval(u_w),minval(u_w),maxval(v_w),minval(v_w)
  endif

  allocate(ud_local(nlon,nlatp,1))
  allocate(vd_local(nlonp,nlat,1))
!
! Create sub-communicator to handle each file
  key=mype+1
  if(mype < nlev) then
     color = 1
  else
     color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(mpi_comm_world,color,key,new_comm,ierror)
  if ( ierror /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',ierror
     call mpi_abort(mpi_comm_world,104,ierror)
  endif

  ud_local = 0.0
  vd_local = 0.0
  if (MPI_COMM_NULL /= new_comm) then
     write(*,*) " chgres winds ",mype
     call chgres_winds_main(gridx, gridy, u_s, v_s, u_w, v_w, ud_local, vd_local)
  endif

  deallocate(u_s)
  deallocate(v_s)
  deallocate(u_w)
  deallocate(v_w)
!
! collect all ud vd to root
!
  if(mype==0) then
     allocate(ud(nlon,nlatp,nlev))
     allocate(vd(nlonp,nlat,nlev))
  endif
!
  if (MPI_COMM_NULL /= new_comm) then
     write(*,*) " gather winds ",mype,maxval(ud_local),maxval(vd_local)
     call MPI_Gather(ud_local, nlon*nlatp, MPI_DOUBLE, ud, nlon*nlatp, MPI_DOUBLE, 0, new_comm, ierr)
     call MPI_Gather(vd_local, nlonp*nlat, MPI_DOUBLE, vd, nlonp*nlat, MPI_DOUBLE, 0, new_comm, ierr)
  end if
!
  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  deallocate(ud_local)
  deallocate(vd_local)
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  !
  if(mype==0) then
    do k=1,nlev
      write(*,*) "ud=",k,maxval(ud(:,:,k)),minval(ud(:,:,k))
    enddo
    do k=1,nlev
      write(*,*) "vd=",k,maxval(vd(:,:,k)),minval(vd(:,:,k))
    enddo
!
!           allocate(Atm_ps(nlon,nlat))
    allocate(Atm_u(nlon,nlatp,nlev-1))
    allocate(Atm_v(nlonp,nlat,nlev-1))
!           Atm_ps=psc
!
! wind vertical remap
!
    write(*,*) " vertical remap the wind"
    call remap_dwinds_main(nlev, nlev-1, ak0, bk0, Atm_ak, Atm_bk, psc, ud, vd, &
                      1, nlon, 1, nlat, Atm_u, Atm_v, Atm_ps)
    deallocate(ud)
    deallocate(vd)
!
    do k=1,nlev-1
       write(*,*) "Atm_u=",k,maxval(Atm_u(:,:,k)),minval(Atm_u(:,:,k))
    enddo
    do k=1,nlev-1
       write(*,*) "Atm_v=",k,maxval(Atm_v(:,:,k)),minval(Atm_v(:,:,k))
    enddo

!          call check(nf90_open(trim(filecold(1)),IOR(NF90_WRITE, NF90_MPIIO), cdfid))
!          call check(nf90_inq_dimid(cdfid, "lat", dimid_lat))
!          call check(nf90_inq_dimid(cdfid, "lon", dimid_lon))
!          call check(nf90_inq_dimid(cdfid, "latp", dimid_latp))
!          call check(nf90_inq_dimid(cdfid, "lonp", dimid_lonp))

!          call check(nf90_create("cold2warm_uv.nc",IOR(nf90_netcdf4, nf90_clobber), ncid=cdfid))
          call check(nf90_open("cold2warm_all.nc",nf90_write,cdfid))
          call check(nf90_inq_dimid(cdfid, "lat", dimid_lat))
          call check(nf90_inq_dimid(cdfid, "lon", dimid_lon))
          call check(nf90_inq_dimid(cdfid, "nlev",nlevid))
          call check(nf90_redef(cdfid))
!
! define nlev, and u and v
!          call check( nf90_def_dim(cdfid, "lat",  nlat,   dimid_lat))
!          call check( nf90_def_dim(cdfid, "lon",  nlon,   dimid_lon))
          call check( nf90_def_dim(cdfid, "latp", nlatp,  dimid_latp))
          call check( nf90_def_dim(cdfid, "lonp", nlonp,  dimid_lonp))
!          call check( nf90_def_dim(cdfid, "nlev", nlev-1, nlevid))
        ! u_cold2fv3 (nlev, latp, lon)
          dimids(1:3) = [dimid_lon, dimid_latp, nlevid]
          chunksizes(1:4) = [nlon, nlatp, 1, 1]
          call check(nf90_def_var(cdfid, "u_cold2fv3", NF90_FLOAT, dimids(1:3), varid))
        !  call check(nf90_def_var_chunking(cdfid, varid, NF90_CHUNKED, chunksizes))

        ! v_cold2fv3 (nlev, lat, lonp)
          dimids(1:3) = [dimid_lonp, dimid_lat, nlevid]
          chunksizes(1:4) = [nlonp, nlat, 1, 1]
          call check(nf90_def_var(cdfid, "v_cold2fv3", NF90_FLOAT, dimids(1:3), varid))
        !  call check(nf90_def_var_chunking(cdfid, varid, NF90_CHUNKED, chunksizes))

          call check(nf90_enddef(cdfid))
!
           ! Write u_cold2fv3
          call check(nf90_inq_varid(cdfid, "u_cold2fv3", varid))
          start = [1, 1, 1, 0]
          count = [nlon, nlatp, nlev-1, 0]
          allocate(d3r4(nlon, nlatp, nlev-1))
          d3r4=Atm_u
          call check(nf90_put_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
          deallocate(d3r4)
!
           ! Write v_cold2fv3
          call check(nf90_inq_varid(cdfid, "v_cold2fv3", varid))
          start = [1, 1, 1, 0]
          count = [nlonp, nlat, nlev-1, 0]
          allocate(d3r4(nlonp, nlat, nlev-1))
          d3r4=Atm_v
          call check(nf90_put_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
          deallocate(d3r4)

          iret=nf90_close(cdfid)

          deallocate(Atm_u)
          deallocate(Atm_v)
  endif
!
! release memory
!
  deallocate(psc)
  deallocate(gridx)
  deallocate(gridy)
  if(mype==0) deallocate(Atm_ps)
  if(mype==0) write(*,*) "done with wind remap"

  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
  if(mype==0)  write(6,*) "=== RRFS PRE_BLENDING SUCCESS ==="
  call MPI_FINALIZE(ierror)
!    
contains
    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, "NetCDF error on rank ", mype, ": ", trim(nf90_strerror(status))
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
    end subroutine check

END PROGRAM pre_blending

SUBROUTINE mype_read(ncioid,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)

  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_get_var,nf90_inq_varid
!
  integer,intent(in) :: ncioid
!
! MPI distribution array
  character(len=20),intent(in) :: mype_varname
  integer,intent(in) :: mype_vartype
  integer,intent(in) :: mype_nx,mype_ny
  integer,intent(in) :: mype_lbegin,mype_lend
  real(4),intent(inout) :: d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend)
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:)

  integer :: startloc(3)
  integer :: countloc(3)
  integer :: var_id
  integer :: ilev
!
!
     if(mype_vartype==5) then
        allocate(tmpd3r4(mype_nx,mype_ny,1))
     elseif(mype_vartype==6) then
        allocate(tmpd3r8(mype_nx,mype_ny,1))
     else
        write(6,*) 'Warning, unknown datatype'
     endif
!
! now read in each fields from fv3 file
!
     do ilev=mype_lbegin,mype_lend
        startloc=(/1,1,ilev/)
        countloc=(/mype_nx,mype_ny,1/)

        iret=nf90_inq_varid(ncioid,trim(adjustl(mype_varname)),var_id)
        if(mype_vartype==5) then
           iret=nf90_get_var(ncioid,var_id,tmpd3r4,start=startloc,count=countloc)
           d3r4(:,:,ilev)=tmpd3r4(:,:,1)
        elseif(mype_vartype==6) then
           iret=nf90_get_var(ncioid,var_id,tmpd3r8,start=startloc,count=countloc)
           d3r4(:,:,ilev)=tmpd3r8(:,:,1)
        endif
     enddo  ! ilev
! release memory

     if(mype_vartype==5) then
        deallocate(tmpd3r4)
     elseif(mype_vartype==6) then
        deallocate(tmpd3r8)
     else
        write(6,*) 'Warning, unknown datatype'
     endif

END SUBROUTINE mype_read

SUBROUTINE mype_write(ncioid,nsig,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)

  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_put_var,nf90_inq_varid
!
  integer,intent(in) :: ncioid
!
! MPI distribution array
  character(len=20),intent(in) :: mype_varname
  integer,intent(in) :: mype_vartype
  integer,intent(in) :: mype_nx,mype_ny
  integer,intent(in) :: mype_lbegin,mype_lend
  real(4),intent(inout) :: d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend)
  integer,intent(in) :: nsig
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:)

  integer :: startloc(3)
  integer :: countloc(3)
  integer :: var_id
  integer :: ilev
  character(len=20) :: local_varname

  logical :: if_write
  integer :: i,j,nlevel
!
!
     if_write=.true.
!
     if(trim(adjustl(mype_varname))=="delp") local_varname="delp_cold2fv3"
     if(trim(adjustl(mype_varname))=="t") local_varname="t_cold2fv3"
     if(trim(adjustl(mype_varname))=="sphum") local_varname="sphum_cold2fv3"

     ilev=mype_lbegin
     nlevel=mype_lend-mype_lbegin+1

     if(mype_lend==nsig) then
       if( nlevel == 1)  then
          nlevel=0
          if_write=.false.
       else
          nlevel=nlevel-1
       endif
     endif
         
     if(if_write) then

        startloc=(/1,1,ilev/)
        countloc=(/mype_nx,mype_ny,nlevel/)

!        write(6,'(a,a20,I5,2f15.6,7I5)') 'writing =',trim(adjustl(mype_varname)), &
!                   mype_vartype,maxval(d3r4(:,:,ilev:ilev+nlevel-1)),minval(d3r4(:,:,ilev:ilev+nlevel-1)),&
!                   nlevel,mype_lbegin,mype_lend,ilev,ilev+nlevel-1,mype_nx,mype_ny
        iret=nf90_inq_varid(ncioid,trim(adjustl(local_varname)),var_id)
        if(mype_vartype==5) then
           allocate(tmpd3r4(mype_nx,mype_ny,nleve1))
           tmpd3r4(:,:,1:nleve1)=d3r4(:,:,ilev:ilev+nlevel-1)
           iret=nf90_put_var(ncioid,var_id,values=tmpd3r4,start=startloc,count=countloc)
           deallocate(tmpd3r4)
        elseif(mype_vartype==6) then
           allocate(tmpd3r8(mype_nx,mype_ny,nleve1))
           tmpd3r8(:,:,1:nlevel)=d3r4(:,:,ilev:ilev+nlevel-1)
           iret=nf90_put_var(ncioid,var_id,values=tmpd3r8,start=startloc,count=countloc)
           deallocate(tmpd3r8)
        else
           write(6,*) 'Warning, unknown datatype'
        endif
        write(*,*) "finished==",trim(adjustl(local_varname)),ilev

     endif

END SUBROUTINE mype_write

SUBROUTINE mype_write2(ncioid,nsig,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)

  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_put_var,nf90_inq_varid
!
  integer,intent(in) :: ncioid
!
! MPI distribution array
  character(len=20),intent(in) :: mype_varname
  integer,intent(in) :: mype_vartype
  integer,intent(in) :: mype_nx,mype_ny
  integer,intent(in) :: mype_lbegin,mype_lend
  real(4),intent(in) :: d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend)
  integer,intent(in) :: nsig
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:,:)

  integer :: startloc(4)
  integer :: countloc(4)
  integer :: var_id
  integer :: ilev
  character(len=20) :: local_varname
!
!
     if(mype_vartype==5) then
        allocate(tmpd3r4(mype_nx,mype_ny,1,1))
     elseif(mype_vartype==6) then
        allocate(tmpd3r8(mype_nx,mype_ny,1,1))
     else
        write(6,*) 'Warning, unknown datatype'
     endif
!
!
!
     do ilev=mype_lbegin,mype_lend

        startloc=(/1,1,ilev,1/)
        countloc=(/mype_nx,mype_ny,1,1/)

        if(trim(adjustl(mype_varname))=="delp") local_varname="delp_cold2fv3"
        if(trim(adjustl(mype_varname))=="t") local_varname="t_cold2fv3"
        if(trim(adjustl(mype_varname))=="sphum") local_varname="sphum_cold2fv3"
        iret=nf90_inq_varid(ncioid,trim(adjustl(local_varname)),var_id)
        if(ilev < nsig) then
!           write(6,'(a,a20,I5,2f15.6)') 'writing =',trim(adjustl(mype_varname)), &
!                ilev,maxval(d3r4(:,:,ilev)),minval(d3r4(:,:,ilev))
           if(mype_vartype==5) then
              tmpd3r4(:,:,1,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r4,start=startloc,count=countloc)
           elseif(mype_vartype==6) then
              tmpd3r8(:,:,1,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r8,start=startloc,count=countloc)
           endif
        endif
     enddo  ! ilev

! release memory

     if(mype_vartype==5) then
        deallocate(tmpd3r4)
     elseif(mype_vartype==6) then
        deallocate(tmpd3r8)
     else
        write(6,*) 'Warning, unknown datatype'
     endif

END SUBROUTINE mype_write2


subroutine reorg(lon,lat,fin,fout)
  implicit none
  integer, intent(in) :: lon,lat
  real,intent(in) :: fin(lat,lon)
  real(8),intent(inout) :: fout(lon,lat)
  integer :: i,j

  do i=1,lon
     do j=1,lat
        fout(i,j)=fin(j,i)
     enddo
  enddo

end subroutine reorg

subroutine reorg_ad(lon,lat,fin,fout)
  implicit none
  integer, intent(in) :: lon,lat
  real,intent(inout) :: fin(lat,lon)
  real(8),intent(in) :: fout(lon,lat)
  integer :: i,j

  do i=1,lon
     do j=1,lat
        fin(j,i)=fout(i,j)
     enddo
  enddo

end subroutine reorg_ad

