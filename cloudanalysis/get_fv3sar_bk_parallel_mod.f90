module get_fv3sar_bk_parall_mod
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! program:  get_fv3sar_bk_mod driver to get fv3sar background
! calculation
!
!   PRGMMR: Ming Hu          ORG: GSD/AMB        DATE: 2020-8-10
!
! ABSTRACT: 
!
! PROGRAM HISTORY LOG:
!    2020-08-10  Hu  Add NCO document block
!
! USAGE:
!   INPUT FILES: 
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 
!   MACHINE:  Linux cluster (JET)
!
!$$$
!
!_____________________________________________________________________
!

  use kinds,   only: r_single,i_kind, r_kind
  use constants, only: init_constants,init_constants_derived
  use constants, only: rd,h1000,rd_over_cp,grav
  use gsi_rfv3io_tten_mod, only: gsi_rfv3io_get_grid_specs
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_read,gsi_fv3ncdf2d_read
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_write
  use gsi_rfv3io_tten_mod, only: nlon_regional,nlat_regional,nsig_regional
  use gsi_rfv3io_tten_mod, only: eta1_ll
  use namelist_mod, only: fv3_io_layout_y
  use rapidrefresh_cldsurf_mod, only: l_cld_uncertainty

  use general_sub2grid_simple_mod, only: sub2grid_info

  implicit none
  private

  type(sub2grid_info) :: s

  public :: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  public :: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg,ges_qnr,ges_qni,ges_qnc,ges_qcf
  public :: unc_ql,unc_qi,unc_qr,unc_qs,unc_qg
  public :: aeta1_ll,pt_ll
  public :: pt_regional
  public :: xlon,xlat,xland,soiltbk
  public :: fv3_io_layout_end
!
  public :: read_fv3sar_bk_full
  public :: release_mem_fv3sar
  public :: update_fv3sar
  public :: update_fv3sar_unc
  public :: read_fv3sar_init
!
  public :: lon2, lat2, nsig
  public :: mype_istart,mype_jstart
!
! MPI variables
  integer :: mypeLocal
! MPI distribution array
  integer :: mype_fileid
  character(len=20) :: mype_varname
  integer :: mype_vartype
  integer :: mype_nx,mype_ny
  integer :: mype_lbegin,mype_lend
  integer :: lon2, lat2, nsig
  integer :: mype_istart,mype_jstart
  integer :: num_fields
  integer :: ntotalcore
  integer,allocatable :: kbegin(:),kend(:)
  character(len=20),allocatable :: varname(:)
!
! background
!
  integer, parameter :: maxcores=10
  real(r_single)  pt_regional
  
  integer(i_kind) :: iunit
  data iunit / 15 /

  real(r_single),allocatable:: t_bk(:,:,:)
  real(r_single),allocatable:: h_bk(:,:,:)
  real(r_single),allocatable:: p_bk(:,:,:)
  real(r_single),allocatable:: ps_bk(:,:)
  real(r_single),allocatable:: zh(:,:)
  real(r_single),allocatable:: q_bk(:,:,:)

  real(r_single),allocatable:: pblh(:,:)         ! PBL height (grid coordinate)

  real(r_kind) :: pt_ll
  real(r_kind),allocatable :: aeta1_ll(:)   !
!
  integer :: fv3sar_bg_opt
  integer(i_kind),allocatable:: fv3_io_layout_begin(:)
  integer(i_kind),allocatable:: fv3_io_layout_end(:)
  integer(i_kind),allocatable:: io_layout_tmp(:)
!
! hydrometeors
  real(r_single),allocatable :: ges_ql(:,:,:)  ! cloud water
  real(r_single),allocatable :: ges_qi(:,:,:)  ! could ice
  real(r_single),allocatable :: ges_qr(:,:,:)  ! rain
  real(r_single),allocatable :: ges_qs(:,:,:)  ! snow
  real(r_single),allocatable :: ges_qg(:,:,:)  ! graupel
  real(r_single),allocatable :: ges_qnr(:,:,:) ! rain number concentration
  real(r_single),allocatable :: ges_qni(:,:,:) ! cloud ice number concentration
  real(r_single),allocatable :: ges_qnc(:,:,:) ! cloud water number concentration
  real(r_single),allocatable :: ges_qcf(:,:,:) ! cloud fraction
!
! uncertainties
  real(r_single),allocatable :: unc_ql(:,:,:)  ! cloud water
  real(r_single),allocatable :: unc_qi(:,:,:)  ! could ice
  real(r_single),allocatable :: unc_qr(:,:,:)  ! rain
  real(r_single),allocatable :: unc_qs(:,:,:)  ! snow
  real(r_single),allocatable :: unc_qg(:,:,:)  ! graupel

! fix files
  real(r_single),allocatable :: xlon(:,:)
  real(r_single),allocatable :: xlat(:,:)
  real(r_single),allocatable :: xland(:,:)
  real(r_single),allocatable :: soiltbk(:,:)
!
! backgrond files
  character(len=80) :: dynvars   !='fv3_dynvars'
  character(len=80) :: tracers   !='fv3_tracer'
  character(len=80) :: tracers_unc   !='fv3_tracer_unc'
  character(len=80) :: sfcvars   !='fv3_sfcdata'
  character(len=80) :: phyvars   !='fv3_phydata'
  character(len=80) :: gridspec  !='fv3_grid_spec'
  character(len=80) :: ak_bk     !='fv3_akbk'
!
  integer, parameter    :: max_num_file = 10
  integer, parameter    :: filename_len=200
  character(len=filename_len)  :: filelist(max_num_file)
  integer                      :: numvar(max_num_file)
  character(len=filename_len)  :: varlist(max_num_file)
  integer :: totalnumfiles
!
  integer :: rfv3io_mype

contains
!

  subroutine read_fv3sar_init(fv3sar_bg_opt_in,mype,npe)
!
  use mpi
  use module_ncio, only : ncio
  use module_ncfile_stat, only : ncfile_stat
  use module_mpi_arrange, only : mpi_io_arrange
! 
  implicit none
  integer,intent(in) :: fv3sar_bg_opt_in
  integer,intent(in) :: mype
  integer,intent(in) :: npe

  type(ncfile_stat) :: ncfs_all
  type(mpi_io_arrange) :: mpiioarg
  type(ncio) :: geo

  integer(i_kind) :: ix,iy,nlat
  character(len=80) :: thisfv3file
  integer :: ierr,ierror
  real(r_kind),allocatable :: eta1_ll0(:)
  integer(i_kind),allocatable:: fv3_io_layout_begin0(:)
  integer(i_kind),allocatable:: fv3_io_layout_end0(:)
!
  fv3sar_bg_opt=fv3sar_bg_opt_in    ! 0=restart, 1=input
  gridspec='fv3_grid_spec'
  ak_bk='fv3_akbk'
  dynvars='fv3_dynvars'
  tracers='fv3_tracer'
  tracers_unc='fv3_tracer_unc'
  sfcvars='fv3_sfcdata'
  phyvars='fv3_phydata'
!
  if(mype==0) then
    call gsi_rfv3io_get_grid_specs(trim(gridspec),trim(ak_bk),ierr)
  endif
  call MPI_Bcast(nlon_regional, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(nlat_regional, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(nsig_regional, 1, mpi_integer, 0, mpi_comm_world, ierror)
  if(.not.allocated(eta1_ll)) allocate(eta1_ll(nsig_regional+1))
  call MPI_Bcast(eta1_ll, nsig_regional+1, mpi_double, 0, mpi_comm_world, ierror)
!
  allocate(fv3_io_layout_begin(fv3_io_layout_y))
  allocate(fv3_io_layout_end(fv3_io_layout_y))
  if(mype==0) then
    fv3_io_layout_begin=0
    fv3_io_layout_end=0

    if(fv3_io_layout_y > 1) then
       iy=0
       do ix=1,fv3_io_layout_y
           write(thisfv3file,'(a,a1,I4.4)') trim(gridspec),".",ix-1
           call geo%open(trim(thisfv3file),'r',200)
           call geo%get_dim("grid_yt",nlat)
           call geo%close
           fv3_io_layout_begin(ix)=iy+1
           fv3_io_layout_end(ix)=iy+nlat
           iy=fv3_io_layout_end(ix)
       enddo  ! find dimension
    else
         fv3_io_layout_begin=1
         fv3_io_layout_end=nlat_regional
    endif
    write(6,'(a,20I5)') '  end index for each subdomain',fv3_io_layout_end
    
  endif
  call MPI_Bcast(fv3_io_layout_begin, fv3_io_layout_y, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(fv3_io_layout_end, fv3_io_layout_y, mpi_integer, 0, mpi_comm_world, ierror)

  filelist(1)=trim(gridspec)
  filelist(2)=trim(sfcvars)
  filelist(3)=trim(dynvars)
  filelist(4)=trim(tracers)
  filelist(5)=trim(phyvars)
  
  if(mype==0) then
    numvar(1)=2
    varlist(1)='grid_lont grid_latt'
    if(fv3sar_bg_opt==0) then
       numvar(2)=3
       varlist(2)='t2m slmsk tsfc'
       numvar(3)=4
       varlist(3)='DZ T delp phis'
       numvar(4)=9
       varlist(4)='sphum liq_wat ice_wat rainwat snowwat graupel water_nc ice_nc rain_nc'
       numvar(5)=1
       varlist(5)='mynn_3d_cldfra_bl'
       totalnumfiles=5
    else
       numvar(2)=3
       varlist(2)='t2m slmsk tisfc'
       numvar(3)=5
       varlist(3)='zh t delp phis ps'
       numvar(4)=6
       varlist(4)='sphum liq_wat ice_wat rainwat snowwat graupel'
       totalnumfiles=4
    endif

!
! find dimension of each field
!
     call ncfs_all%init(totalnumfiles,filelist(1:totalnumfiles), &
                        numvar(1:totalnumfiles),varlist(1:totalnumfiles))
     call ncfs_all%fill_dims()
!
!  distibute variables to each core
!
     call mpiioarg%init(npe)
     call mpiioarg%arrange(ncfs_all)
     num_fields=ncfs_all%num_totalvl
     if(fv3sar_bg_opt==0) then
       nsig=mpiioarg%nlvlmax
     else
       nsig=mpiioarg%nlvlmax-1
     endif
     ntotalcore=mpiioarg%ntotalcore

     call ncfs_all%close()
  endif

  write(6,*) 'background type =',fv3sar_bg_opt
  call MPI_Bcast(num_fields, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(nsig, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(ntotalcore, 1, mpi_integer, 0, mpi_comm_world, ierror)
  call MPI_Bcast(totalnumfiles, 1, mpi_integer, 0, mpi_comm_world, ierror)

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
  
  call MPI_Scatter(mpiioarg%fileid, 1, mpi_integer, mype_fileid, 1, &
                   mpi_integer,0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%varname, 20, mpi_character, mype_varname, 20, &
                   mpi_character, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%vartype, 1, mpi_integer, mype_vartype, 1, &
                   mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%nx, 1, mpi_integer, mype_nx, 1, mpi_integer, 0, &
                   MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%ny, 1, mpi_integer, mype_ny, 1, mpi_integer, 0, &
                   MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%lvlbegin, 1, mpi_integer, mype_lbegin, 1, &
                   mpi_integer, 0, MPI_COMM_WORLD,ierror)
  call MPI_Scatter(mpiioarg%lvlend, 1, mpi_integer, mype_lend, 1, mpi_integer, &
                   0, MPI_COMM_WORLD,ierror)

  if(mype==0) call mpiioarg%close()

  end subroutine read_fv3sar_init

  subroutine read_fv3sar_bk_full(mype)
!
  use mpi
  use netcdf, only: nf90_open,nf90_close,nf90_noerr
  use netcdf, only: nf90_get_var
  use netcdf, only: nf90_nowrite
  use netcdf, only: nf90_Inquire_Dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_strerror

  use general_sub2grid_simple_mod, only: general_sub2grid_create_info
  use general_sub2grid_simple_mod, only: general_sub2grid_destroy_info
  use general_sub2grid_simple_mod, only: general_grid2sub

  use constants, only: rd,h1000,rd_over_cp,grav

  implicit none
  integer,intent(in) :: mype
!
! sub communicator
  integer :: new_comm
  integer :: color, key
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:)
  real(4),allocatable :: d3r4(:,:,:)
  real(4),allocatable :: sub_vars(:,:,:)

  integer :: startloc(3)
  integer :: countloc(3)
  integer :: ncioid,var_id
!
!
  character(len=80) :: filename
  integer :: n,i,j,k,iret,ilev,iens
  integer :: ierror
!
! ===============================================================================
!
! Create sub-communicator to handle each file
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= totalnumfiles) then
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
     if(mype_vartype==5) then
        allocate(tmpd3r4(mype_nx,mype_ny,1))
     elseif(mype_vartype==6) then
        allocate(tmpd3r8(mype_nx,mype_ny,1))
     else
        write(6,*) 'Warning, unknown datatype'
     endif

     filename = trim(filelist(mype_fileid))
     iret=nf90_open(trim(filename),nf90_nowrite,ncioid,comm=new_comm,info=MPI_INFO_NULL)
     if(iret/=nf90_noerr) then
         write(6,*)' problem opening ', trim(filename),'fileid=',mype_fileid,', Status =',iret
         write(6,*)  nf90_strerror(iret)
         call flush(6)
         stop(333)
     endif

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
        write(6,'(a,a20,I5,2f20.7)') 'reading =',trim(adjustl(mype_varname)), &
                   ilev,maxval(d3r4(:,:,ilev)),minval(d3r4(:,:,ilev))
     enddo  ! ilev
     iret=nf90_close(ncioid)

     if(mype_vartype==5) then
        deallocate(tmpd3r4)
     elseif(mype_vartype==6) then
        deallocate(tmpd3r8)
     else
        write(6,*) 'Warning, unknown datatype'
     endif
  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierror)

  write(6,'(a5,2x,10a10)') "core","varname","lvlbegin","lvlend"
  do k=1,ntotalcore
       write(6,'(I5,2x,a10,10I10)') k,trim(varname(k)),kbegin(k),kend(k)
  enddo
  write(6,*)"======================================================================="

  call general_sub2grid_create_info(s,mype,ntotalcore,nlat_regional,nlon_regional,nsig,num_fields,kbegin,kend)

  mype_istart=s%istart(mype+1)
  mype_jstart=s%jstart(mype+1)

  allocate(sub_vars(s%lat2,s%lon2,s%num_fields))
  call general_grid2sub(s,d3r4,sub_vars)
  deallocate(d3r4)

  lon2=s%lon2
  lat2=s%lat2
  nsig=s%nsig
! hydrometeors
  allocate(ges_ql(lon2,lat2,nsig))   ! cloud water
  allocate(ges_qi(lon2,lat2,nsig))   ! cloud ice
  allocate(ges_qr(lon2,lat2,nsig))   ! rain
  allocate(ges_qs(lon2,lat2,nsig))   ! snow
  allocate(ges_qg(lon2,lat2,nsig))   ! graupel
  allocate(ges_qnr(lon2,lat2,nsig))   ! rain number concentration
  allocate(ges_qni(lon2,lat2,nsig))   ! cloud ice number concentration
  allocate(ges_qnc(lon2,lat2,nsig))   ! cloud water number concentration
  allocate(ges_qcf(lon2,lat2,nsig))   ! cloud fraction
  if(l_cld_uncertainty) then
!   hydrometeor uncertainties
    allocate(unc_ql(lon2,lat2,nsig))   ! cloud water
    allocate(unc_qi(lon2,lat2,nsig))   ! cloud ice
    allocate(unc_qr(lon2,lat2,nsig))   ! rain
    allocate(unc_qs(lon2,lat2,nsig))   ! snow
    allocate(unc_qg(lon2,lat2,nsig))   ! graupel
  endif
!
! fix files
  allocate(xlon(lon2,lat2))
  allocate(xlat(lon2,lat2))
  allocate(xland(lon2,lat2))
  allocate(soiltbk(lon2,lat2))
!
! dynamic 
  allocate(t_bk(lon2,lat2,nsig))
  allocate(h_bk(lon2,lat2,nsig))
  allocate(p_bk(lon2,lat2,nsig))
  allocate(q_bk(lon2,lat2,nsig))
  allocate(ps_bk(lon2,lat2))
  allocate(zh(lon2,lat2))
  allocate(pblh(lon2,lat2))
!
  i=0
  allocate(tmpd3r4(lon2,lat2,nsig+1))
  do n=1,ntotalcore
     do ilev=kbegin(n),kend(n)
        i=i+1
        if(trim(varname(n))=="grid_lont") call reorg(lon2,lat2,sub_vars(:,:,i),xlon(:,:))
        if(trim(varname(n))=="grid_latt") call reorg(lon2,lat2,sub_vars(:,:,i),xlat(:,:))
        if(trim(varname(n))=="slmsk")     call reorg(lon2,lat2,sub_vars(:,:,i),xland(:,:))
        if(trim(varname(n))=="tsfc" .or. trim(varname(n))=="tisfc") &
                                          call reorg(lon2,lat2,sub_vars(:,:,i),soiltbk(:,:))
        k=nsig-ilev+1
        if(trim(varname(n))=="liq_wat")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_ql(:,:,k))
        if(trim(varname(n))=="ice_wat")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qi(:,:,k))
        if(trim(varname(n))=="rainwat")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qr(:,:,k))
        if(trim(varname(n))=="snowwat")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qs(:,:,k))
        if(trim(varname(n))=="graupel")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qg(:,:,k))
        if(trim(varname(n))=="water_nc")  call reorg(lon2,lat2,sub_vars(:,:,i),ges_qnc(:,:,k))
        if(trim(varname(n))=="ice_nc")    call reorg(lon2,lat2,sub_vars(:,:,i),ges_qni(:,:,k))
        if(trim(varname(n))=="rain_nc")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qnr(:,:,k))
        if(trim(varname(n))=="mynn_3d_cldfra_bl")   call reorg(lon2,lat2,sub_vars(:,:,i),ges_qcf(:,:,k))

        if(trim(varname(n))=="T" .or. trim(varname(n))=="t") call reorg(lon2,lat2,sub_vars(:,:,i),t_bk(:,:,k))
        if(trim(varname(n))=="DZ")      call reorg(lon2,lat2,sub_vars(:,:,i),h_bk(:,:,k))
        if(trim(varname(n))=="delp")    call reorg(lon2,lat2,sub_vars(:,:,i),p_bk(:,:,k))
        if(trim(varname(n))=="sphum")   call reorg(lon2,lat2,sub_vars(:,:,i),q_bk(:,:,k))

        if(trim(varname(n))=="phis")    call reorg(lon2,lat2,sub_vars(:,:,i),pblh(:,:))
        if(trim(varname(n))=="ps")      call reorg(lon2,lat2,sub_vars(:,:,i),ps_bk(:,:))

        k=nsig+1-ilev+1
        if(trim(varname(n))=="zh")      call reorg(lon2,lat2,sub_vars(:,:,i),tmpd3r4(:,:,k))
     enddo
  enddo

  deallocate(sub_vars)
!
  write(6,*) 'xlon=',maxval(xlon(:,:)),minval(xlon(:,:))
  write(6,*) 'xlat=',maxval(xlat(:,:)),minval(xlat(:,:))
  write(6,*) 'xlan=',maxval(xland(:,:)),minval(xland(:,:))
  write(6,*) 'soiltbk=',maxval(soiltbk(:,:)),minval(soiltbk(:,:))
  do k=1,s%nsig
     write(6,*) 'ql=',k,maxval(ges_ql(:,:,k)),minval(ges_ql(:,:,k))
  enddo
  do k=1,s%nsig
     write(6,*) 'qi=',k,maxval(ges_qi(:,:,k)),minval(ges_qi(:,:,k))
  enddo
!
  zh=pblh/grav
  write(6,*) 'zh=',maxval(zh),minval(zh)
!  
  if(fv3sar_bg_opt==1) then
     do k=1,nsig
        h_bk(:,:,k)=(tmpd3r4(:,:,k)+tmpd3r4(:,:,k+1))*0.5_r_kind
     enddo
     do k=1,nsig
        h_bk(:,:,k)=h_bk(:,:,k)-tmpd3r4(:,:,1)
     enddo
  else
     tmpd3r4(:,:,1)=0
     do k=2,nsig+1
        tmpd3r4(:,:,k)=tmpd3r4(:,:,k-1)-h_bk(:,:,k-1)
     enddo
     do k=1,nsig
        h_bk(:,:,k)=(tmpd3r4(:,:,k+1)+tmpd3r4(:,:,k))*0.5_r_kind
     enddo
  endif
  do k=1,nsig
     write(6,*) 'Z==',k,maxval(h_bk(:,:,k)),minval(h_bk(:,:,k))
  enddo
!
!   get pressure
  tmpd3r4=0.0
  tmpd3r4(:,:,s%nsig+1)=eta1_ll(s%nsig+1)*1000.0_r_kind
  do k=nsig,1,-1
     tmpd3r4(:,:,k)=tmpd3r4(:,:,k+1)+p_bk(:,:,k)
  enddo

  if(fv3sar_bg_opt==0) then
     ps_bk(:,:)=tmpd3r4(:,:,1)*0.01_r_kind
     write(6,*) 'Ps(restart)=',maxval(ps_bk),minval(ps_bk)
  else
     ps_bk=ps_bk*0.01_r_kind
     write(6,*) 'Ps(input)=',maxval(ps_bk),minval(ps_bk)
  endif
  do k=1,nsig
        p_bk(:,:,k)=(tmpd3r4(:,:,k)+tmpd3r4(:,:,k+1))*0.5_r_kind*0.01_r_kind
  enddo
  do k=1,nsig
     write(6,*) 'P==',k,maxval(p_bk(:,:,k)),minval(p_bk(:,:,k))
  enddo
 
  deallocate(tmpd3r4)
!
!    get temperature
  do k=1,nsig
     t_bk(:,:,k)=t_bk(:,:,k)*(h1000/p_bk(:,:,k))**rd_over_cp
     write(6,*) 'T==',k,maxval(t_bk(:,:,k)),minval(t_bk(:,:,k))
  enddo

!    get moisture 
  do k=1,nsig
     q_bk(:,:,k)=q_bk(:,:,k)/(1.0_r_kind-q_bk(:,:,k)) ! conver to mixing ratio
     write(6,*) 'q==',k,maxval(q_bk(:,:,k)),minval(q_bk(:,:,k))
  enddo
!
! 2.5 calculate PBL height
! 
   pblh=0.0
   call calc_pbl_height(lat2,lon2,s%nsig,q_bk,t_bk,p_bk,pblh)
   write(6,*) 'max.min, pblh=',maxval(pblh), minval(pblh)
!
!  call general_sub2grid_destroy_info(s)
!
end subroutine read_fv3sar_bk_full

subroutine reorg(lon,lat,fin,fout)
  implicit none
  integer, intent(in) :: lon,lat
  real(r_single),intent(in) :: fin(lat,lon)
  real(r_single),intent(inout) :: fout(lon,lat)
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
  real(r_single),intent(inout) :: fin(lat,lon)
  real(r_single),intent(in) :: fout(lon,lat)
  integer :: i,j

  do i=1,lon
     do j=1,lat
        fin(j,i)=fout(i,j)
     enddo
  enddo

end subroutine reorg_ad

subroutine release_mem_fv3sar
  implicit none
!
!  release memory
!
    ! write(6,*) 'core', mype ,',finished, now release memory bk'
     if(allocated(t_bk)) deallocate(t_bk)
     if(allocated(h_bk)) deallocate(h_bk)
     if(allocated(p_bk)) deallocate(p_bk)
     if(allocated(ps_bk))deallocate(ps_bk)
     if(allocated(zh))   deallocate(zh)
     if(allocated(q_bk)) deallocate(q_bk)
     if(allocated(pblh)) deallocate(pblh)
!
!  release memory
!
     write(6,*) 'core', 1 ,',finished, now release memory or hydrometeors'
     if(allocated(ges_ql))  deallocate(ges_ql)
     if(allocated(ges_qi))  deallocate(ges_qi)
     if(allocated(ges_qr))  deallocate(ges_qr)
     if(allocated(ges_qs))  deallocate(ges_qs)
     if(allocated(ges_qg))  deallocate(ges_qg)
     if(allocated(ges_qnr)) deallocate(ges_qnr)
     if(allocated(ges_qni)) deallocate(ges_qni)
     if(allocated(ges_qnc)) deallocate(ges_qnc)
     if(l_cld_uncertainty) then
       write(6,*) 'core', 1 ,', release memory for hydrometeor uncertainties'
       if(allocated(unc_ql))  deallocate(unc_ql)
       if(allocated(unc_qi))  deallocate(unc_qi)
       if(allocated(unc_qr))  deallocate(unc_qr)
       if(allocated(unc_qs))  deallocate(unc_qs)
       if(allocated(unc_qg))  deallocate(unc_qg)
     endif
!
!  release memory
!
     write(6,*) 'core', 1 ,',finished, now release memory for fix variables'
     if(allocated(xlon))  deallocate(xlon)
     if(allocated(xlat))  deallocate(xlat)
     if(allocated(xland))  deallocate(xland)
     if(allocated(soiltbk))  deallocate(soiltbk)

     if(allocated(fv3_io_layout_end))  deallocate(fv3_io_layout_end)
     if(allocated(fv3_io_layout_begin))  deallocate(fv3_io_layout_begin)
  
end subroutine release_mem_fv3sar

subroutine update_fv3sar(mype)

  use mpi
  use netcdf, only: nf90_open,nf90_close,nf90_noerr
  use netcdf, only: nf90_write,nf90_put_var
  use netcdf, only: nf90_Inquire_Dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_strerror

  use general_sub2grid_simple_mod, only: general_sub2grid_destroy_info
  use general_sub2grid_simple_mod, only: general_sub2grid
  implicit none
!
  integer,intent(in) :: mype
!
! sub communicator
  integer :: new_comm
  integer :: color, key
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:)
  real(4),allocatable :: d3r4(:,:,:)
  real(r_single),allocatable :: sub_vars(:,:,:)

  integer :: startloc(3)
  integer :: countloc(3)
  integer :: ncioid,var_id
!
  character(len=80) :: filename
  integer :: n,i,j,k,iret,ilev,iens
  integer :: ierror
!
  integer :: lon2,lat2,nsig
!
!
!
  lon2=s%lon2
  lat2=s%lat2
  nsig=s%nsig
!
  allocate(sub_vars(lat2,lon2,s%num_fields))
  sub_vars=0.0

  i=0
  do n=1,ntotalcore
     do ilev=kbegin(n),kend(n)
        i=i+1
        k=nsig-ilev+1
        if(trim(varname(n))=="liq_wat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_ql(:,:,k))
        if(trim(varname(n))=="ice_wat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qi(:,:,k))
        if(trim(varname(n))=="rainwat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qr(:,:,k))
        if(trim(varname(n))=="snowwat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qs(:,:,k))
        if(trim(varname(n))=="graupel")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qg(:,:,k))
        if(trim(varname(n))=="water_nc")  call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qnc(:,:,k))
        if(trim(varname(n))=="ice_nc")    call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qni(:,:,k))
        if(trim(varname(n))=="rain_nc")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qnr(:,:,k))
        if(trim(varname(n))=="T" .or. trim(varname(n))=="t") call reorg_ad(lon2,lat2,sub_vars(:,:,i),t_bk(:,:,k))
        if(trim(varname(n))=="sphum")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),q_bk(:,:,k))

     enddo
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierror)
  allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))
  call general_sub2grid(s,sub_vars,d3r4)
  deallocate(sub_vars)

! Create sub-communicator to handle each file
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= totalnumfiles) then
     if(trim(adjustl(mype_varname))=="liq_wat" .or. &
        trim(adjustl(mype_varname))=="ice_wat" .or. &
        trim(adjustl(mype_varname))=="rainwat" .or. &
        trim(adjustl(mype_varname))=="snowwat" .or. &
        trim(adjustl(mype_varname))=="graupel" .or. &
        trim(adjustl(mype_varname))=="water_nc" .or. &
        trim(adjustl(mype_varname))=="ice_nc" .or. &
        trim(adjustl(mype_varname))=="rain_nc" .or. &
        trim(adjustl(mype_varname))=="T" .or. &
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

        if(mype_vartype==5) then
           allocate(tmpd3r4(mype_nx,mype_ny,1))
        elseif(mype_vartype==6) then
           allocate(tmpd3r8(mype_nx,mype_ny,1))
        else
           write(6,*) 'Warning, unknown datatype'
        endif

        filename = trim(filelist(mype_fileid))
        iret=nf90_open(trim(filename),nf90_write,ncioid,comm=new_comm,info=MPI_INFO_NULL)
        if(iret/=nf90_noerr) then
            write(6,*)' problem opening ', trim(filename),'fileid=',mype_fileid,', Status =',iret
            write(6,*)  nf90_strerror(iret)
            call flush(6)
            stop(444)
        endif

        do ilev=mype_lbegin,mype_lend
           write(6,'(a,a20,I5,2f15.6)') 'writing =',trim(adjustl(mype_varname)), &
                   ilev,maxval(d3r4(:,:,ilev)),minval(d3r4(:,:,ilev))

           startloc=(/1,1,ilev/)
           countloc=(/mype_nx,mype_ny,1/)

           iret=nf90_inq_varid(ncioid,trim(adjustl(mype_varname)),var_id)
           if(mype_vartype==5) then
              tmpd3r4(:,:,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r4,start=startloc,count=countloc)
           elseif(mype_vartype==6) then
              tmpd3r8(:,:,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r8,start=startloc,count=countloc)
           endif
        enddo  ! ilev
        iret=nf90_close(ncioid)

        if(mype_vartype==5) then
           deallocate(tmpd3r4)
        elseif(mype_vartype==6) then
           deallocate(tmpd3r8)
        else
           write(6,*) 'Warning, unknown datatype'
        endif

  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  deallocate(d3r4)

!  call general_sub2grid_destroy_info(s)

end subroutine update_fv3sar

subroutine update_fv3sar_unc(mype)

  use mpi
  use netcdf, only: nf90_open,nf90_close,nf90_noerr
  use netcdf, only: nf90_write,nf90_put_var
  use netcdf, only: nf90_Inquire_Dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_strerror

  use general_sub2grid_simple_mod, only: general_sub2grid_destroy_info
  use general_sub2grid_simple_mod, only: general_sub2grid
  implicit none
!
  integer,intent(in) :: mype
!
! sub communicator
  integer :: new_comm
  integer :: color, key
!
! array
  real(4),allocatable :: tmpd3r4(:,:,:)
  real(8),allocatable :: tmpd3r8(:,:,:)
  real(4),allocatable :: d3r4(:,:,:)
  real(r_single),allocatable :: sub_vars(:,:,:)

  integer :: startloc(3)
  integer :: countloc(3)
  integer :: ncioid,var_id
!
  character(len=80) :: filename
  integer :: n,i,j,k,iret,ilev,iens
  integer :: ierror
!
  integer :: lon2,lat2,nsig
!
!
  write(6,*)'get_fv3sar_bk_parallel_mod: begin update_fv3sar_unc'
!
  lon2=s%lon2
  lat2=s%lat2
  nsig=s%nsig
!
  allocate(sub_vars(lat2,lon2,s%num_fields))
  sub_vars=0.0

  i=0
  do n=1,ntotalcore
     do ilev=kbegin(n),kend(n)
        i=i+1
        k=nsig-ilev+1
        if(trim(varname(n))=="liq_wat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_ql(:,:,k))
        if(trim(varname(n))=="ice_wat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qi(:,:,k))
        if(trim(varname(n))=="rainwat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qr(:,:,k))
        if(trim(varname(n))=="snowwat")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qs(:,:,k))
        if(trim(varname(n))=="graupel")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qg(:,:,k))
        if(trim(varname(n))=="water_nc")  call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qnc(:,:,k))
        if(trim(varname(n))=="ice_nc")    call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qni(:,:,k))
        if(trim(varname(n))=="rain_nc")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),ges_qnr(:,:,k))
        if(trim(varname(n))=="T" .or. trim(varname(n))=="t") call reorg_ad(lon2,lat2,sub_vars(:,:,i),t_bk(:,:,k))
        if(trim(varname(n))=="sphum")   call reorg_ad(lon2,lat2,sub_vars(:,:,i),q_bk(:,:,k))

     enddo
  enddo

  call mpi_barrier(MPI_COMM_WORLD,ierror)
  allocate(d3r4(mype_nx,mype_ny,mype_lbegin:mype_lend))
  call general_sub2grid(s,sub_vars,d3r4)
  deallocate(sub_vars)

! Create sub-communicator to handle each file
  key=mype+1
  if(mype_fileid > 0 .and. mype_fileid <= totalnumfiles) then
     if(trim(adjustl(mype_varname))=="liq_wat" .or. &
        trim(adjustl(mype_varname))=="ice_wat" .or. &
        trim(adjustl(mype_varname))=="rainwat" .or. &
        trim(adjustl(mype_varname))=="snowwat" .or. &
        trim(adjustl(mype_varname))=="graupel" .or. &
        trim(adjustl(mype_varname))=="water_nc" .or. &
        trim(adjustl(mype_varname))=="ice_nc" .or. &
        trim(adjustl(mype_varname))=="rain_nc" .or. &
        trim(adjustl(mype_varname))=="T" .or. &
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

        if(mype_vartype==5) then
           allocate(tmpd3r4(mype_nx,mype_ny,1))
        elseif(mype_vartype==6) then
           allocate(tmpd3r8(mype_nx,mype_ny,1))
        else
           write(6,*) 'Warning, unknown datatype'
        endif

        filename = trim(filelist(mype_fileid))
        write(6,*)'get_fv3sar_bk_parallel_mod: filename: ', filename
        if(filename==trim(tracers)) filename=trim(tracers_unc)
        write(6,*)'get_fv3sar_bk_parallel_mod: filename: ', filename
        iret=nf90_open(trim(filename),nf90_write,ncioid,comm=new_comm,info=MPI_INFO_NULL)
        if(iret/=nf90_noerr) then
            write(6,*)' problem opening ', trim(filename),'fileid=',mype_fileid,', Status =',iret
            write(6,*)  nf90_strerror(iret)
            call flush(6)
            stop(444)
        endif

        do ilev=mype_lbegin,mype_lend
           write(6,'(a,a20,I5,2f15.6)') 'writing =',trim(adjustl(mype_varname)), &
                   ilev,maxval(d3r4(:,:,ilev)),minval(d3r4(:,:,ilev))

           startloc=(/1,1,ilev/)
           countloc=(/mype_nx,mype_ny,1/)

           iret=nf90_inq_varid(ncioid,trim(adjustl(mype_varname)),var_id)
           if(mype_vartype==5) then
              tmpd3r4(:,:,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r4,start=startloc,count=countloc)
           elseif(mype_vartype==6) then
              tmpd3r8(:,:,1)=d3r4(:,:,ilev)
              iret=nf90_put_var(ncioid,var_id,tmpd3r8,start=startloc,count=countloc)
           endif
        enddo  ! ilev
        iret=nf90_close(ncioid)

        if(mype_vartype==5) then
           deallocate(tmpd3r4)
        elseif(mype_vartype==6) then
           deallocate(tmpd3r8)
        else
           write(6,*) 'Warning, unknown datatype'
        endif

  endif

  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_free(new_comm,iret)
  endif

  deallocate(d3r4)

  call general_sub2grid_destroy_info(s)

end subroutine update_fv3sar_unc

end module get_fv3sar_bk_parall_mod
