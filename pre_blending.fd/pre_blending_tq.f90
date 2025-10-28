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
  logical :: lexist
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
! check file existing status
!
  filecold(1)='out.atm.tile7.nc'
  inquire(file=trim(filecold(1)),exist=lexist)
  if(lexist) then
    if(mype==0) write(6,*) "file out.atm.tile7.nc exists"
  else
     write(6,*) "Error, out.atm.tile7.nc does not exist",mype
     call mpi_abort(mpi_comm_world,101,ierror)
  endif

  inquire(file="gfs_ctrl.nc",exist=lexist)
  if(lexist) then
    if(mype==0) write(6,*) "file gfs_ctrl.nc exists"
  else
     write(6,*) "Error, gfs_ctrl.nc does not exist",mype
     call mpi_abort(mpi_comm_world,101,ierror)
  endif

  inquire(file="fv_core.res.nc",exist=lexist)
  if(lexist) then
    if(mype==0) write(6,*) "file fv_core.res.nc exists"
  else
     write(6,*) "Error, fv_core.res.nc does not exist",mype
     call mpi_abort(mpi_comm_world,101,ierror)
  endif

  inquire(file="C3463_oro_data.tile7.halo0.nc",exist=lexist)
  if(lexist) then
    if(mype==0) write(6,*) "file C3463_oro_data.tile7.halo0.nc exists"
  else
     write(6,*) "Error, C3463_oro_data.tile7.halo0.nc does not exist",mype
     call mpi_abort(mpi_comm_world,101,ierror)
  endif
!
! read dimensions
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
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
!-------------------------------------------------------------------
! now working on scalars 
!-------------------------------------------------------------------
!
!
  numvar(1)=6
  numvar(2)=1
  varlist(1)='ps o3mr delp t sphum zh'
  varlist(2)='orog_filt'
  filecold(1)='out.atm.tile7.nc'
  filecold(2)='C3463_oro_data.tile7.halo0.nc'
  if(mype==0) then
!
! find dimension of each field
!
     call ncfs_all%init(2,filecold, numvar, varlist)
     call ncfs_all%fill_dims()
!
!  distibute variables to each core
!
     call mpiioarg%init(npe)
     call mpiioarg%arrange(ncfs_all)
     ntotalcore=mpiioarg%ntotalcore
     num_fields=ncfs_all%num_totalvl
     nsig=nlev

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
  call MPI_Bcast(ntotalcore, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(num_fields, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(nsig, 1, mpi_integer, 0, mpi_comm_world, ierror)
  allocate(kbegin(ntotalcore))
  allocate(kend(ntotalcore))
  allocate(varname(ntotalcore))
  if(mype==0) then
     kbegin=mpiioarg%lvlbegin
     kend=mpiioarg%lvlend
     varname=mpiioarg%varname
  endif
  call MPI_Bcast(kbegin, ntotalcore, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(kend, ntotalcore, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(varname, ntotalcore*20, mpi_character, 0, mpi_comm_world, ierror)

  if(mype==0) call mpiioarg%close()
  call mpi_barrier(MPI_COMM_WORLD,ierror)

! Create sub-communicator to handle each file
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= 2) then
     color = mype_fileid
  else
     color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(mpi_comm_world,color,key,new_comm,ierror)
  if ( ierror /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',ierror
     call mpi_abort(mpi_comm_world,101,ierror)
  endif
!
! read 2D field from each file using sub communicator
!
  allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))
  if (MPI_COMM_NULL /= new_comm) then

     iret=nf90_open(trim(filecold(mype_fileid)),nf90_nowrite,ncioid,comm=new_comm,info=MPI_INFO_NULL)
     if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(filecold(1)), ' Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
     endif

     call mype_read(ncioid,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)
     write(6,'(a10,2f20.6)') trim(adjustl(mype_varname)),maxval(d3r4(:,:,:)),minval(d3r4(:,:,:))

     call check(nf90_close(ncioid))
  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  if(mype==0) write(6,*)"all scalers are read in====="
  call mpi_barrier(MPI_COMM_WORLD,ierror)

  call general_sub2grid_create_info(s,mype,ntotalcore,mype_nx,mype_ny,nsig,num_fields,kbegin,kend)

  mype_istart=s%istart(mype+1)
  mype_jstart=s%jstart(mype+1)

  allocate(sub_vars(s%lat2,s%lon2,s%num_fields))
  call general_grid2sub(s,d3r4,sub_vars)
  deallocate(d3r4)

  lon2=s%lon2
  lat2=s%lat2
  nsig=s%nsig

  allocate(ps_local(lon2,lat2))
  allocate(delp_local(lon2,lat2,nsig))
  allocate(zh_local(lon2,lat2,nsig+1))
  allocate(omga_local(lon2,lat2,nsig))
  allocate(t_local(lon2,lat2,nsig))
  allocate(qa_local(lon2,lat2,nsig,1))
  allocate(Atm_phis_local(lon2,lat2))

  write(*,*) mype,lon2,lat2,nsig,mype_istart,mype_jstart
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  i=0
  do n=1,ntotalcore
     do ilev=kbegin(n),kend(n)
        i=i+1
        if(trim(varname(n))=="ps") call reorg(lon2,lat2,sub_vars(:,:,i),ps_local(:,:))
        if(trim(varname(n))=="orog_filt") then
            call reorg(lon2,lat2,sub_vars(:,:,i),Atm_phis_local(:,:))
            Atm_phis_local(:,:)=Atm_phis_local(:,:)*9.80665
        endif

        k=ilev
        if(trim(varname(n))=="o3mr") call reorg(lon2,lat2,sub_vars(:,:,i),omga_local(:,:,k))
        if(trim(varname(n))=="delp") call reorg(lon2,lat2,sub_vars(:,:,i),delp_local(:,:,k))
        if(trim(varname(n))=="t") call reorg(lon2,lat2,sub_vars(:,:,i),t_local(:,:,k))
        if(trim(varname(n))=="sphum") call reorg(lon2,lat2,sub_vars(:,:,i),qa_local(:,:,k,1))
        if(trim(varname(n))=="zh") call reorg(lon2,lat2,sub_vars(:,:,i),zh_local(:,:,k))
     enddo
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierror)
  deallocate(sub_vars)
  if(mype==0) then
    write(*,*) "ps=",maxval(ps_local),minval(ps_local)
    write(*,*) "orog_filt=",maxval(Atm_phis_local),minval(Atm_phis_local)
    do k=1,nsig
      write(*,*) "w=",k,maxval(omga_local(:,:,k)),minval(omga_local(:,:,k))
    enddo
    do k=1,nsig
      write(*,*) "delp=",k,maxval(delp_local(:,:,k)),minval(delp_local(:,:,k))
    enddo
    do k=1,nsig
      write(*,*) "t=",k,maxval(t_local(:,:,k)),minval(t_local(:,:,k))
    enddo
    do k=1,nsig
      write(*,*) "sphum=",k,maxval(qa_local(:,:,k,1)),minval(qa_local(:,:,k,1))
    enddo
    do k=1,nsig+1
      write(*,*) "zh=",k,maxval(zh_local(:,:,k)),minval(zh_local(:,:,k))
    enddo
  endif

  allocate(Atm_ps(lon2,lat2))
  allocate(Atm_delp(lon2,lat2,nsig-1))
  allocate(Atm_pt(lon2,lat2,nsig-1))
  allocate(Atm_q(lon2,lat2,nsig-1,1))

  call remap_scalar_main(nlev, nlev-1, 1, ak0, bk0, Atm_ak, Atm_bk, ps_local, qa_local, &
                           zh_local, omga_local, t_local, 1, lon2, 1, lat2, &
                           Atm_pt, Atm_q, Atm_delp, Atm_phis_local, Atm_ps)

  ps_local=Atm_ps                    
  do k=1,nsig-1
     qa_local(:,:,k,1)=Atm_q(:,:,k,1)
     t_local(:,:,k)=Atm_pt(:,:,k)
     delp_local(:,:,k)=Atm_delp(:,:,k)
  enddo

  if(mype==0) then
    do k=1,nsig-1
      write(*,*) "pt=",k,maxval(Atm_pt(:,:,k)),minval(Atm_pt(:,:,k))
    enddo
    do k=1,nsig-1
      write(*,*) "q=",k,maxval(Atm_q(:,:,k,1)),minval(Atm_q(:,:,k,1))
    enddo
    do k=1,nsig-1
      write(*,*) "delp=",k,maxval(Atm_delp(:,:,k)),minval(Atm_delp(:,:,k))
    enddo
  endif
  deallocate(Atm_ps)
  deallocate(Atm_delp)
  deallocate(Atm_pt)
  deallocate(Atm_q)
!
!
!
  allocate(sub_vars(lat2,lon2,s%num_fields))
  sub_vars=0.0

  i=0
  do n=1,ntotalcore
     do ilev=kbegin(n),kend(n)
        i=i+1
        k=ilev
        if(k<nsig) then
           if(trim(varname(n))=="delp") call reorg_ad(lon2,lat2,sub_vars(:,:,i),delp_local(:,:,k))
           if(trim(varname(n))=="t") call reorg_ad(lon2,lat2,sub_vars(:,:,i),t_local(:,:,k))
           if(trim(varname(n))=="sphum") call reorg_ad(lon2,lat2,sub_vars(:,:,i),qa_local(:,:,k,1))
           if(trim(varname(n))=="ps") call reorg_ad(lon2,lat2,sub_vars(:,:,i),ps_local(:,:))
        endif
     enddo
  enddo
!
!  release memory
!
  deallocate(Atm_phis_local)
  deallocate(qa_local)
  deallocate(t_local)
  deallocate(omga_local)
  deallocate(zh_local)
  deallocate(delp_local)
  deallocate(ps_local)
!
! distribute from sub to full grid
!
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))
  call general_sub2grid(s,sub_vars,d3r4)
  deallocate(sub_vars)

  call general_sub2grid_destroy_info(s)
!  write(6,'(a10,2i10,2f15.7)') trim(adjustl(mype_varname)),mype_lbegin,mype_lend,maxval(d3r4(:,:,:)),minval(d3r4(:,:,:))

  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
!  write to the file
!
! get atm_ps, this will be used for wind
  allocate(Atm_ps(mype_nx,mype_ny))
  if(trim(adjustl(mype_varname))=="ps") then
     Atm_ps(:,:)=d3r4(:,:,mype_lbegin)
     write(*,*) mype, "ps=",maxval(Atm_ps), minval(Atm_ps)
     open(12,file='Atm_ps.bin',form='unformatted',status='unknown')
         write(12) mype_nx,mype_ny
         write(12) Atm_ps
     close(12)
  endif
!
! Create sub-communicator to handle each file
!
  create_new_file=.false.
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= 1 ) then
     if(trim(adjustl(mype_varname))=="delp" .or. &
        trim(adjustl(mype_varname))=="t" .or. &
        trim(adjustl(mype_varname))=="sphum") then
        color = mype_fileid
     else
        color = MPI_UNDEFINED
     endif
  else
     color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(mpi_comm_world,color,key,new_comm,ierror)
  if ( ierror /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',ierror
     call mpi_abort(mpi_comm_world,101,ierror)
  endif
!
! read 2D field from each file using sub communicator
!
  if (MPI_COMM_NULL /= new_comm) then

     if(create_new_file) then
        call check(nf90_create("cold2warm_all.nc",IOR(nf90_netcdf4, nf90_clobber), &
                               comm=new_comm, info=MPI_INFO_NULL, ncid=cdfid))
        if(iret/=nf90_noerr) then
            write(6,*)' problem creating cold2warm_all.nc ', ', Status =',iret
            write(6,*)  nf90_strerror(iret)
            call flush(6)
            stop(444)
        endif

        call check(nf90_set_fill(cdfid, NF90_NOFILL, oldMode))

        call check(nf90_redef(cdfid))
          call check( nf90_def_dim(cdfid, "lat",  nlat,   dimid_lat))
          call check( nf90_def_dim(cdfid, "lon",  nlon,   dimid_lon))
          call check( nf90_def_dim(cdfid, "nlev", nlev-1, nlevid))
        ! t_cold2fv3 (nlev, lat, lon)
          dimids(1:3) = [dimid_lon, dimid_lat, nlevid]
          chunksizes(1:4) = [nlon, nlat, 1, 1]
          call check(nf90_def_var(cdfid, "t_cold2fv3", NF90_FLOAT, dimids(1:3), varid))
          call check(nf90_def_var_chunking(cdfid, varid, NF90_CHUNKED, chunksizes(1:3)))
          call check(nf90_var_par_access(cdfid, varid, NF90_COLLECTIVE))

        ! delp_cold2fv3 (nlev, lat, lon)
          call check(nf90_def_var(cdfid, "delp_cold2fv3", NF90_FLOAT, dimids(1:3), varid))
          call check(nf90_def_var_chunking(cdfid, varid, NF90_CHUNKED, chunksizes(1:3)))
          call check(nf90_var_par_access(cdfid, varid, NF90_COLLECTIVE))

        ! sphum_cold2fv3 (nlev, lat, lon)
          call check(nf90_def_var(cdfid, "sphum_cold2fv3", NF90_FLOAT, dimids(1:3), varid))
          call check(nf90_def_var_chunking(cdfid, varid, NF90_CHUNKED, chunksizes(1:3)))
          call check(nf90_var_par_access(cdfid, varid, NF90_COLLECTIVE))

        call check(nf90_enddef(cdfid))
        call mype_write(cdfid,nsig,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)
        call check(nf90_close(cdfid))
     else
        iret=nf90_open("cold2warm_all.nc",nf90_write,cdfid,comm=new_comm,info=MPI_INFO_NULL)
        if(iret/=nf90_noerr) then
            write(6,*)' problem opening cold2warm_all.nc', ' Status =',iret
            write(6,*)  nf90_strerror(iret)
            call flush(6)
            stop(444)
        endif
        call mype_write2(cdfid,nsig,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)
        call check(nf90_close(cdfid))
     endif

  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  deallocate(d3r4)
!  
  if(mype==0) write(*,*) "done with T Q and delp remap"
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
  if(mype==0)  write(6,*) "=== RRFS PRE_BLENDING T,Q, and delp SUCCESS ==="
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

