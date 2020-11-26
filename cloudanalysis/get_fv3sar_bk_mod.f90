module get_fv3sar_bk_mod
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
  use gsi_rfv3io_tten_mod, only: bg_fv3regfilenameg,fv3sar_bg_opt
  use gsi_rfv3io_tten_mod, only: rfv3io_mype
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_read,gsi_fv3ncdf2d_read
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_write
  use gsi_rfv3io_tten_mod, only: nlon_regional,nlat_regional,nsig_regional
  use gsi_rfv3io_tten_mod, only: eta1_ll
  use mpi_mod, only: npe, mype,mpi_comm_world
  use mpi_mod, only: mpi_finish

  implicit none
  private

  public :: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  public :: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg,ges_qnr,ges_qni,ges_qnc
  public :: aeta1_ll,pt_ll
  public :: pt_regional
  public :: xlon,xlat,xland,soiltbk
!
  public :: read_fv3sar_bk
  public :: release_mem_fv3sar_bk
  public :: read_fv3sar_hydr
  public :: release_mem_fv3sar_hydr
  public :: read_fv3sar_fix
  public :: release_mem_fv3sar_fix
  public :: update_fv3sar
!
! MPI variables
  integer :: mypeLocal
!
! background
!
  integer, parameter :: maxcores=1
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

! hydrometeors
  real(r_single),allocatable :: ges_ql(:,:,:)  ! cloud water
  real(r_single),allocatable :: ges_qi(:,:,:)  ! could ice
  real(r_single),allocatable :: ges_qr(:,:,:)  ! rain
  real(r_single),allocatable :: ges_qs(:,:,:)  ! snow
  real(r_single),allocatable :: ges_qg(:,:,:)  ! graupel
  real(r_single),allocatable :: ges_qnr(:,:,:) ! rain number concentration
  real(r_single),allocatable :: ges_qni(:,:,:) ! cloud ice number concentration
  real(r_single),allocatable :: ges_qnc(:,:,:) ! cloud water number concentration

! fix files
  real(r_single),allocatable :: xlon(:,:)
  real(r_single),allocatable :: xlat(:,:)
  real(r_single),allocatable :: xland(:,:)
  real(r_single),allocatable :: soiltbk(:,:)
!
! backgrond files
  character(len=:),allocatable :: dynvars   !='fv3_dynvars'
  character(len=:),allocatable :: tracers   !='fv3_tracer'
  character(len=:),allocatable :: sfcvars   !='fv3_sfcdata'
!

contains
!
  subroutine read_fv3sar_bk
!
  implicit none
  integer, parameter :: maxcores=1
  real(r_single),allocatable :: field1(:)
!
  integer(i_kind) :: mype_2d,mype_t,mype_q,mype_p,mype_ql
  integer(i_kind) :: i,j,k
  integer(i_kind) :: ierror

!
! ===============================================================================
!
! MPI setup
  write(6,*) 'total cores for this run is ',npe
  if(npe < maxcores) then
     write(6,*) 'ERROR, this run must use ',maxcores,' or more cores !!!'
     call mpi_finish
     stop 1234
  endif

  mypeLocal=mype+1
  if(mypeLocal <= maxcores) then
!
! 2.0  open and read background dimesion
!          
     rfv3io_mype=mype
     fv3sar_bg_opt=0    ! 0=restart, 1=input
     call bg_fv3regfilenameg%init
     mype_t=mype
     mype_q=mype
     mype_p=mype
     mype_ql=mype
     mype_2d=mype

     dynvars= bg_fv3regfilenameg%dynvars
     tracers= bg_fv3regfilenameg%tracers
     sfcvars= bg_fv3regfilenameg%sfcdata

     write(6,*) 'read in background ==========>'

! 2.1 read in background fields
!
!     t_bk        - 3D background potential temperature (K)
!     ps_bk       - 2D background surface pressure (hPa)
!     q           - 3D moisture (water vapor mixing ratio kg/kg)
!     zh          - terrain
!     pbk         - 3D background pressure  (hPa)
!     hbk         - 3D height above the ground (not the sea level)
!
     call gsi_rfv3io_get_grid_specs(bg_fv3regfilenameg,ierror)

     allocate(ps_bk(nlon_regional,nlat_regional))
     allocate(zh(nlon_regional,nlat_regional))
     call gsi_fv3ncdf2d_read(dynvars,'phis','PHIS',ps_bk,mype_2d)
     zh=ps_bk/grav
     write(6,*) 'zh=',maxval(zh),minval(zh)
     call gsi_fv3ncdf2d_read(sfcvars,'t2m','T2M',ps_bk,mype_2d)
     write(6,*) 't2m=',maxval(ps_bk),minval(ps_bk)
!
!    t_bk and q_bk are used as temperal arrary to calculate H and P
     allocate(t_bk(nlon_regional,nlat_regional,nsig_regional))
     allocate(q_bk(nlon_regional,nlat_regional,nsig_regional+1))

!    get height
     allocate(h_bk(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(dynvars,'DZ','zh',t_bk,mype_t)
     if(t_bk(1,1,nsig_regional) > 10000.0_r_single) then
        do k=1,nsig_regional
           h_bk(:,:,k)=t_bk(:,:,k) - zh(:,:)
        enddo
        fv3sar_bg_opt=1
     else
        q_bk(:,:,1)=0
        do k=2,nsig_regional+1
           do j=1,nlat_regional
              do i=1,nlon_regional
                 q_bk(i,j,k)=q_bk(i,j,k-1)-t_bk(i,j,k-1)
              enddo
           enddo
        enddo 
        do k=1,nsig_regional
           do j=1,nlat_regional
              do i=1,nlon_regional
                 h_bk(i,j,k)=(q_bk(i,j,k+1)+q_bk(i,j,k))/2.0_r_kind
              enddo
           enddo
        enddo 
     endif
     do k=1,nsig_regional
        write(6,*) 'Z==',k,maxval(h_bk(:,:,k)),minval(h_bk(:,:,k))
     enddo 

!    get pressure
     allocate(p_bk(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(dynvars,'delp','DELP',t_bk,mype_t)
     q_bk(:,:,nsig_regional+1)=eta1_ll(nsig_regional+1)*1000.0_r_kind
     do k=nsig_regional,1,-1
        q_bk(:,:,k)=q_bk(:,:,k+1)+t_bk(:,:,k)
     enddo 
     if(fv3sar_bg_opt==0) then
        ps_bk(:,:)=q_bk(:,:,1)*0.01_r_kind
        write(6,*) 'Ps(restart)=',maxval(ps_bk),minval(ps_bk)
     else
        call gsi_fv3ncdf2d_read(dynvars,'ps','PS',ps_bk,mype_2d)
        ps_bk=ps_bk*0.01_r_kind
        write(6,*) 'Ps(input)=',maxval(ps_bk),minval(ps_bk)
     endif
     do k=1,nsig_regional
        p_bk(:,:,k)=(q_bk(:,:,k)+q_bk(:,:,k+1))*0.5_r_kind*0.01_r_kind
     enddo 
     do k=1,nsig_regional
        write(6,*) 'P==',k,maxval(p_bk(:,:,k)),minval(p_bk(:,:,k))
     enddo 

!    get temperature
     call gsi_fv3ncdf_read(dynvars,'T','t',t_bk,mype_t)
     do k=1,nsig_regional
        t_bk(:,:,k)=t_bk(:,:,k)*(h1000/p_bk(:,:,k))**rd_over_cp
        write(6,*) 'T==',k,maxval(t_bk(:,:,k)),minval(t_bk(:,:,k))
     enddo 

!    get moisture 
     deallocate(q_bk)
     allocate(q_bk(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'sphum','SPHUM',q_bk,mype_t)
     do k=1,nsig_regional
        q_bk(:,:,k)=q_bk(:,:,k)/(1.0_r_kind-q_bk(:,:,k)) ! conver to mixing ratio
        write(6,*) 'q==',k,maxval(q_bk(:,:,k)),minval(q_bk(:,:,k))
     enddo 
! 
! 2.5 calculate PBL height
! 
     allocate(pblh(nlon_regional,nlat_regional))
     call calc_pbl_height(nlon_regional,nlat_regional,nsig_regional,q_bk,t_bk,p_bk,pblh)
     write(6,*) 'max.min, pblh=',maxval(pblh), minval(pblh)

!     deallocate(q_bk)
!     deallocate(ps_bk)

!
! 5.10 update
!

!     call gsi_fv3ncdf_write(tracers,'radar_tten',ges_tten,mype_t)
!
  endif ! mypeLocal <= maxcores

  call MPI_Barrier(mpi_comm_world, ierror)

end subroutine read_fv3sar_bk

subroutine release_mem_fv3sar_bk
  implicit none
!
!  release memory
!
     write(6,*) 'core', mype ,',finished, now release memory bk'
     if(allocated(t_bk)) deallocate(t_bk)
     if(allocated(h_bk)) deallocate(h_bk)
     if(allocated(p_bk)) deallocate(p_bk)
     if(allocated(ps_bk))deallocate(ps_bk)
     if(allocated(zh))   deallocate(zh)
     if(allocated(q_bk)) deallocate(q_bk)
     if(allocated(pblh)) deallocate(pblh)

end subroutine release_mem_fv3sar_bk

subroutine read_fv3sar_hydr
!
  implicit none
  integer, parameter :: maxcores=1
  integer(i_kind) :: mype_2d,mype_t,mype_q,mype_p,mype_ql
  integer(i_kind) :: i,j,k
  integer(i_kind) :: ierror

!
! ===============================================================================
!
! MPI setup
  write(6,*) 'total cores for this run is ',npe
  if(npe < maxcores) then
     write(6,*) 'ERROR, this run must use ',maxcores,' or more cores !!!'
     call mpi_finish
     stop 1234
  endif

  mypeLocal=mype+1
  if(mypeLocal <= maxcores) then
!
! 2.0  open and read background dimesion
!          
     rfv3io_mype=mype

     write(6,*) 'read in hybrometeors==========>'
! 2.1 read in background fields
!
!    get cloud water
     allocate(ges_ql(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'LIQ_WAT','liq_wat',ges_ql,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'ql==',k,maxval(ges_ql(:,:,k)),minval(ges_ql(:,:,k))
     enddo 
!
!    get cloud ice
     allocate(ges_qi(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'ICE_WAT','ice_wat',ges_qi,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qi==',k,maxval(ges_qi(:,:,k)),minval(ges_qi(:,:,k))
     enddo 
!
!    get rain water
     allocate(ges_qr(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'RAINWAT','rainwat',ges_qr,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qr==',k,maxval(ges_qr(:,:,k)),minval(ges_qr(:,:,k))
     enddo 
!
!    get snow water
     allocate(ges_qs(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'SNOWWAT','snowwat',ges_qs,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qs==',k,maxval(ges_qs(:,:,k)),minval(ges_qs(:,:,k))
     enddo 
!
!    get grapel
     allocate(ges_qg(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'GRAUPEL','graupel',ges_qg,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qg==',k,maxval(ges_qg(:,:,k)),minval(ges_qg(:,:,k))
     enddo 
!
!    get cloud ice number concentration
     allocate(ges_qni(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'ICE_NC','ice_nc',ges_qni,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qni==',k,maxval(ges_qni(:,:,k)),minval(ges_qni(:,:,k))
     enddo 
!
!    get rain water number concentration
     allocate(ges_qnr(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'RAIN_NC','rain_nc',ges_qnr,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qnr==',k,maxval(ges_qnr(:,:,k)),minval(ges_qnr(:,:,k))
     enddo 
!
!    get cloud water number concentration
     allocate(ges_qnc(nlon_regional,nlat_regional,nsig_regional))
     call gsi_fv3ncdf_read(tracers,'WATER_NC','water_nc',ges_qnc,rfv3io_mype)
     do k=1,nsig_regional
        write(6,*) 'qnc==',k,maxval(ges_qnc(:,:,k)),minval(ges_qnc(:,:,k))
     enddo 

  endif ! mypeLocal <= maxcores

  call MPI_Barrier(mpi_comm_world, ierror)

end subroutine read_fv3sar_hydr

subroutine release_mem_fv3sar_hydr
!
!  release memory
!
     write(6,*) 'core', mype ,',finished, now release memory or hydrometeors'
     if(allocated(ges_ql))  deallocate(ges_ql)
     if(allocated(ges_qi))  deallocate(ges_qi)
     if(allocated(ges_qr))  deallocate(ges_qr)
     if(allocated(ges_qs))  deallocate(ges_qs)
     if(allocated(ges_qg))  deallocate(ges_qg)
     if(allocated(ges_qnr)) deallocate(ges_qnr)
     if(allocated(ges_qni)) deallocate(ges_qni)
     if(allocated(ges_qnc)) deallocate(ges_qnc)

end subroutine release_mem_fv3sar_hydr

subroutine read_fv3sar_fix
!
  implicit none
  integer, parameter :: maxcores=1
  integer(i_kind) :: i,j,k
  integer(i_kind) :: ierror
  integer(i_kind) :: rfv2io_mype
  integer :: NCID
  CHARACTER*180   geofile

!
! ===============================================================================
!
! MPI setup
  write(6,*) 'total cores for this run is ',npe
  if(npe < maxcores) then
     write(6,*) 'ERROR, this run must use ',maxcores,' or more cores !!!'
     call mpi_finish
     stop 1234
  endif

  mypeLocal=mype+1
  if(mypeLocal <= maxcores) then
!
! 2.0  open and read background dimesion
!          
     rfv2io_mype=mype
     write(6,*) 'read in latlon and surface ==========>'

! 2.1 read in fix variables from background fields
!
     allocate (xlon(nlon_regional,nlat_regional))
     allocate (xlat(nlon_regional,nlat_regional))
     write(geofile,'(a,a)') './', 'fv3_grid_spec'
     call OPEN_geo(geofile, NCID)
     call GET_geo_sngl_fv3sar(NCID,nlon_regional,nlat_regional,xlat,xlon)
     call CLOSE_geo(NCID)
     write(*,*) 'FV3SAR grid'
     write(*,*) 'nlonfv3,nlatfv3=',nlon_regional,nlat_regional
     write(*,*) 'max, min lon=', maxval(xlon),minval(xlon)
     write(*,*) 'max, min lat=', maxval(xlat),minval(xlat)

     allocate (xland(nlon_regional,nlat_regional))
     call gsi_fv3ncdf2d_read(sfcvars,'SLMSK','slmsk',xland,rfv2io_mype)
     write(*,*) 'land mask min and max=',minval(xland),maxval(xland)

     allocate (soiltbk(nlon_regional,nlat_regional))
     call gsi_fv3ncdf2d_read(sfcvars,'TISFC','tisfc',soiltbk,rfv2io_mype)
     write(*,*) 'soil temperature min and max=',minval(soiltbk),maxval(soiltbk)
  endif
!

end subroutine read_fv3sar_fix

subroutine release_mem_fv3sar_fix
  implicit none
!
!  release memory
!
     write(6,*) 'core', mype ,',finished, now release memory for fix variables'
     if(allocated(xlon))  deallocate(xlon)
     if(allocated(xlat))  deallocate(xlat)
     if(allocated(xland))  deallocate(xland)
     if(allocated(soiltbk))  deallocate(soiltbk)
  
end subroutine release_mem_fv3sar_fix

subroutine update_fv3sar

  implicit none
  integer :: rfv3io_mype
     rfv3io_mype=mype
     call gsi_fv3ncdf_write(tracers,'liq_wat',ges_ql,rfv3io_mype)
     call gsi_fv3ncdf_write(tracers,'ice_wat',ges_qi,rfv3io_mype)

end subroutine update_fv3sar

end module get_fv3sar_bk_mod
