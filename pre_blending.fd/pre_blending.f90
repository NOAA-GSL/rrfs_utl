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


!  use netcdf, only: nf90_open,nf90_close,nf90_noerr
!  use netcdf, only: nf90_put_var
!  use netcdf, only: nf90_nowrite,nf90_write
!  use netcdf, only: nf90_inq_varid
!  use netcdf, only: nf90_strerror

  implicit none
! 
  type(ncfile_stat) :: ncfs_all
  type(mpi_io_arrange) :: mpiioarg
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
  integer             :: numvar(1)
  character(len=200)  :: varlist(1)

  character (len=filename_len)   :: filecold(1)
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

 real,allocatable, dimension(:,:)     ::  d2r4
 real,allocatable, dimension(:,:,:)   ::  d3r4

!
  integer :: i,j,k,iret,ilev,ierr
  integer :: grid_cdfid,cdfid
! dimensions
  integer :: nlev, nlat, nlon, nlatp, nlonp, nlevp
  integer :: dimid_lev, dimid_lat, dimid_lon, dimid_latp, dimid_lonp, dimid_nlev
  integer :: local_nx,local_ny
! Variable IDs and I/O arrays
  integer :: varid, dimids(4), start(4), count(4)
  integer :: nlevid

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
!
  filecold(1)='out.atm.tile7.nc'
  call check(nf90_open(trim(filecold(1)), IOR(NF90_WRITE, NF90_MPIIO), cdfid, &
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
!     allocate(vcoord(nlevp,2))
     call check(nf90_inq_varid(grid_cdfid, "vcoord", varid))
!     call check(nf90_get_var(grid_cdfid, varid, vcoord))
!     ak0=vcoord(:,1)
!     bk0=vcoord(:,2)
!     deallocate(vcoord)
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
  else
     allocate(ak0(nlevp), bk0(nlevp))
     allocate(Atm_ak(nlev), Atm_bk(nlev))
  endif
  call MPI_Bcast(ak0, nlevp, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(bk0, nlevp, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Atm_ak, nlev, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Atm_bk, nlev, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
!  if(mype==10) then
!     write(*,*) ak0
!     write(*,*) bk0
!     write(*,*) Atm_ak
!     write(*,*) Atm_bk
!  endif
!
!  get namelist
!
  numvar(1)=1
  varlist(1)='u_s' ! v_s u_w v_w'

!
  allocate(gridx(nlon,nlat))
  allocate(gridy(nlon,nlat))
  allocate(psc(nlon,nlat))

  call check(nf90_open(trim(filecold(1)), nf90_nowrite, cdfid, &
                       comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))
  allocate(d2r4(nlon,nlat))
  if(mype==0) then
     call check(nf90_inq_varid(cdfid, "ps", varid))
     start = [1, 1, 0, 0]
     count = [nlon,nlat, 0, 0]
     call check(nf90_get_var(cdfid, varid, d2r4, start=start(1:2), count=count(1:2)))
     psc=d2r4
  endif
  if(mype==1) then
     call check(nf90_inq_varid(cdfid, "geolon", varid))
     start = [1, 1, 0, 0]
     count = [nlon, nlat, 0, 0]
     call check(nf90_get_var(cdfid, varid, d2r4, start=start(1:2), count=count(1:2)))
     gridx=d2r4
  endif
  if(mype==2) then
     call check(nf90_inq_varid(cdfid, "geolat", varid))
     start = [1, 1, 0, 0]
     count = [nlon, nlat, 0, 0]
     call check(nf90_get_var(cdfid, varid, d2r4, start=start(1:2), count=count(1:2)))
     gridy=d2r4
  endif
  deallocate(d2r4)
  call check(nf90_close(cdfid))

  call MPI_Bcast(psc, nlon*nlat, MPI_DOUBLE , 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gridx, nlon*nlat, MPI_DOUBLE , 1, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gridy, nlon*nlat, MPI_DOUBLE , 2, MPI_COMM_WORLD, ierr)
  if(mype==10) then
     write(*,*) 'psc=',maxval(psc),minval(psc)
     write(*,*) 'gridx=',maxval(gridx),minval(gridx)
     write(*,*) 'gridy=',maxval(gridy),minval(gridy)
  endif

  numvar(1)=1
  varlist(1)='u_s' ! v_s u_w v_w'
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

  if(mype==0) call mpiioarg%close()
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
!
! read 2D field from each file using sub communicator
!
  if (MPI_COMM_NULL /= new_comm) then

     iret=nf90_open(trim(filecold(1)),nf90_nowrite,ncioid,comm=new_comm,info=MPI_INFO_NULL)
     if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(filecold(1)), ' Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
     endif

     allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))

     if(trim(mype_varname)=='u_s')then
        allocate(u_s(nlon,nlatp,mype_lbegin:mype_lend))
        allocate(v_s(nlon,nlatp,mype_lbegin:mype_lend))
        allocate(u_w(nlonp,nlat,mype_lbegin:mype_lend))
        allocate(v_w(nlonp,nlat,mype_lbegin:mype_lend))
        allocate(ud_local(nlon,nlatp,mype_lbegin:mype_lend))
        allocate(vd_local(nlonp,nlat,mype_lbegin:mype_lend))

        local_nx=nlon
        local_ny=nlatp
        local_varname="u_s"
        call mype_read(ncioid,local_nx,local_ny,mype_lbegin,mype_lend,mype_vartype,local_varname,d3r4)
        u_s=d3r4
     !   write(6,'(a10,2f12.6)') trim(adjustl(local_varname)),maxval(u_s(:,:,:)),minval(u_s(:,:,:))
        local_varname="v_s"
        call mype_read(ncioid,local_nx,local_ny,mype_lbegin,mype_lend,mype_vartype,local_varname,d3r4)
        v_s=d3r4
     !   write(6,'(a10,2f12.6)') trim(adjustl(local_varname)),maxval(v_s(:,:,:)),minval(v_s(:,:,:))

        local_nx=nlonp
        local_ny=nlat
        local_varname="u_w"
        call mype_read(ncioid,local_nx,local_ny,mype_lbegin,mype_lend,mype_vartype,local_varname,d3r4)
        u_w=d3r4
     !   write(6,'(a10,2f12.6)') trim(adjustl(local_varname)),maxval(u_w(:,:,:)),minval(u_w(:,:,:))
        local_varname="v_w"
        call mype_read(ncioid,local_nx,local_ny,mype_lbegin,mype_lend,mype_vartype,local_varname,d3r4)
        v_w=d3r4
     !   write(6,'(a10,2f12.6)') trim(adjustl(local_varname)),maxval(v_w(:,:,:)),minval(v_w(:,:,:))

        ud_local = 0.0
        vd_local = 0.0
        call chgres_winds_main(gridx, gridy, u_s, v_s, u_w, v_w, ud_local, vd_local)
!        write(6,'(a10,I10,2f12.6)') 'ud_local ',mype,maxval(ud_local(:,:,:)),minval(ud_local(:,:,:))
!        write(6,'(a10,I10,2f12.6)') 'vd_local ',mype,maxval(vd_local(:,:,:)),minval(vd_local(:,:,:))

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
        call MPI_Gather(ud_local, nlon*nlatp, MPI_DOUBLE, ud, nlon*nlatp, MPI_DOUBLE, 0, new_comm, ierr)
        call MPI_Gather(vd_local, nlonp*nlat, MPI_DOUBLE, vd, nlonp*nlat, MPI_DOUBLE, 0, new_comm, ierr)
!
        deallocate(ud_local)
        deallocate(vd_local)
        if(mype==0) then
           do k=1,nlev
              write(*,*) "ud=",k,maxval(ud(:,:,k)),minval(ud(:,:,k))
           enddo
           do k=1,nlev
              write(*,*) "ud=",k,maxval(vd(:,:,k)),minval(vd(:,:,k))
           enddo
!
           allocate(Atm_ps(nlon,nlat))
           allocate(Atm_u(nlon,nlatp,nlev-1))
           allocate(Atm_v(nlonp,nlat,nlev-1))
           Atm_ps=psc
!
! wind vertical remap
!
           write(*,*) " vertical remap the wind"
           call remap_dwinds_main(nlev, nlevp, ak0, bk0, Atm_ak, Atm_bk, psc, ud, vd, &
                           1, nlon, 1, nlat, Atm_u, Atm_v, Atm_ps)
!
           deallocate(ud)
           deallocate(vd)
           deallocate(Atm_ps)
        endif
     else
        call mype_read(ncioid,mype_nx,mype_ny,mype_lbegin,mype_lend,mype_vartype,mype_varname,d3r4)
     endif

     deallocate(d3r4)
     iret=nf90_close(ncioid)
!
     if(mype==0) then
          do k=1,nlev
              write(*,*) "Atm_u=",k,maxval(Atm_u(:,:,k)),minval(Atm_u(:,:,k))
          enddo
          do k=1,nlev
              write(*,*) "Atm_v=",k,maxval(Atm_v(:,:,k)),minval(Atm_v(:,:,k))
          enddo

          call check(nf90_open(trim(filecold(1)),IOR(NF90_WRITE, NF90_MPIIO), cdfid))

          call check(nf90_inq_dimid(cdfid, "lat", dimid_lat))
          call check(nf90_inq_dimid(cdfid, "lon", dimid_lon))
          call check(nf90_inq_dimid(cdfid, "latp", dimid_latp))
          call check(nf90_inq_dimid(cdfid, "lonp", dimid_lonp))
          call check(nf90_redef(cdfid))
!
! define nlev, and u and v
          call check( nf90_def_dim(cdfid, "nlev", nlev-1, nlevid))
        ! u_cold2fv3 (nlev, latp, lon)
          dimids(1:3) = [dimid_lon, dimid_latp, nlevid]
          call check(nf90_def_var(cdfid, "u_cold2fv3", NF90_FLOAT, dimids(1:3), varid))

        ! v_cold2fv3 (nlev, lat, lonp)
          dimids(1:3) = [dimid_lonp, dimid_lat, nlevid]
          call check(nf90_def_var(cdfid, "v_cold2fv3", NF90_FLOAT, dimids(1:3), varid))

          call check(nf90_enddef(cdfid))
!
write(*,*) 'cehck 5'
           ! Write u_cold2fv3
          call check(nf90_inq_varid(cdfid, "u_cold2fv3", varid))
          start = [1, 1, 1, 0]
          count = [nlon, nlatp, nlev-1, 0]
          allocate(d3r4(nlon, nlatp, nlev-1))
          d3r4=Atm_u
write(*,*) 'cehck 6'
          call check(nf90_put_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
          deallocate(d3r4)
!
write(*,*) 'cehck 7'
           ! Write v_cold2fv3
          call check(nf90_inq_varid(cdfid, "v_cold2fv3", varid))
          start = [1, 1, 1, 0]
          count = [nlonp, nlat, nlev-1, 0]
          allocate(d3r4(nlonp, nlat, nlev-1))
          d3r4=Atm_v
          call check(nf90_put_var(cdfid, varid, d3r4, start=start(1:3), count=count(1:3)))
          deallocate(d3r4)

          deallocate(Atm_u)
          deallocate(Atm_v)

          iret=nf90_close(cdfid)
     endif
!
! release memory
!
     deallocate(psc)
     deallocate(gridx)
     deallocate(gridy)
  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif


  call mpi_barrier(MPI_COMM_WORLD,ierror)
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
     !   write(6,'(a10,I5,2f12.6)') trim(adjustl(mype_varname)),ilev,maxval(d3r4(:,:,ilev)),minval(d3r4(:,:,ilev))
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
