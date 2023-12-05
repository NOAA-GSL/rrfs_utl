program gsdcloudanalysis_ref2tten
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! program:  gsdcloudanalysis_ref2tten   driver for radar tten
! calculation
!
!   PRGMMR: Ming Hu          ORG: GSD/AMB        DATE: 2014-03-28
!
! ABSTRACT: 
!
! PROGRAM HISTORY LOG:
!    2014-03-28  Hu  Add NCO document block
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
!   MACHINE:  Linux cluster 
!
!$$$
!
!_____________________________________________________________________
!

  use mpi
  use kinds,   only: r_single,i_kind, r_kind
  use constants, only: init_constants,init_constants_derived
  use constants, only: rd,h1000,rd_over_cp,grav,half
  use gsi_rfv3io_tten_mod, only: gsi_rfv3io_get_grid_specs
  use gsi_rfv3io_tten_mod, only: bg_fv3regfilenameg,fv3sar_bg_opt
  use gsi_rfv3io_tten_mod, only: rfv3io_mype
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_read,gsi_fv3ncdf2d_read
  use gsi_rfv3io_tten_mod, only: gsi_fv3ncdf_write,gsi_fv3ncdf_append
  use gsi_rfv3io_tten_mod, only: nlon_regional,nlat_regional,nsig_regional
  use gsi_rfv3io_tten_mod, only: eta1_ll

  implicit none
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! dimension
!
  integer(i_kind)    :: rad_missing
  integer(i_kind)    :: lightning_missing
!
! background
!
  integer :: ntime
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
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
!  radar observation : 3D reflectvity in RR grid
!
  character (len=80) :: radarfile
  character (len=80) :: lightningfile
  real(r_single),allocatable :: ref_mos_3d(:,:,:)
  real(r_single),allocatable :: ref_mosaic31(:,:,:)
  INTEGER(i_kind)          :: Nmsclvl_radar,nlon_radar,nlat_radar
  real(r_single)  ::  krad_bot          ! radar bottom level
  INTEGER(i_kind)          :: iunit_radar
  INTEGER(i_kind)          :: iunit_lightning
  INTEGER(i_kind)          :: iunit_satcast
!
!  satcast observation : 2D in RR grid
!
  character (len=80) :: scfile
  real(r_single),allocatable :: satcast_cr(:,:)
  INTEGER(i_kind)            :: nlon_sc, nlat_sc
!
!
!
  real(r_single),allocatable :: ges_tten(:,:,:)
  real(r_kind) :: dfi_lhtp
  logical         :: l_tten_for_convection_only
  logical         :: l_convection_suppress
  real(r_single) :: dfi_radar_latent_heat_time_period
  REAL(r_kind) :: convection_refl_threshold     ! units dBZ
  integer      :: fv3_io_layout_y
  integer      :: timelevel
!
  namelist/setup/ dfi_radar_latent_heat_time_period,convection_refl_threshold, &
                  l_tten_for_convection_only,l_convection_suppress,  &
                  fv3_io_layout_y,timelevel
!
!
  real(r_single),allocatable :: field1(:)
!
  INTEGER(i_kind) :: miss_obs_int
  REAL(r_kind)    :: miss_obs_real
  PARAMETER ( miss_obs_int = -99999  )
  PARAMETER ( miss_obs_real = -99999.0_r_single )

!  real(r_single), allocatable :: sat_ctp(:,:)
!
!  misc.
!
  integer(i_kind),allocatable:: fv3_io_layout_end(:),io_layout_tmp(:)
  integer(i_kind)            :: ybegin,yend
  INTEGER(i_kind) :: i,j,k,n
  integer(i_kind)            :: header1,header2,header3
  integer(i_kind)            :: numlight
  integer(i_kind)            :: nlat_lightning
  integer(i_kind)            :: nlon_lightning
  real(r_single),allocatable   :: lightning_in(:,:)
  real(r_single),allocatable   :: lightning(:,:)
  integer(i_kind) :: mype_2d,mype_t,mype_q,mype_p,mype_ql

  character(len=80) :: dynvars   !='fv3_dynvars'
  character(len=80) :: tracers   !='fv3_tracer'
  character(len=80) :: sfcvars   !='fv3_sfcdata'
  character(len=80) :: phyvars   !='fv3_phydata'
  character(len=80) :: gridspec  !='fv3_grid_spec'

  logical :: ifexist

!
! ===============================================================================
!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  timelevel=1
  fv3_io_layout_y=1
  dfi_radar_latent_heat_time_period=20.0_r_single
  convection_refl_threshold=28.0_r_kind
  l_tten_for_convection_only=.true.
  l_convection_suppress=.false.
!
!  read namelist
!
  inquire(file='namelist.ref2tten', EXIST=ifexist )
  if(ifexist) then
    open(10,file='namelist.ref2tten',status='old')
       read(10,setup)
    close(10)
    if(mype==0) then
      write(*,*) 'Namelist setup are:'
      write(*,setup)
    endif
  else
    if(mype==0) then
      write(*,*) 'No namelist file exist, use default values'
      write(*,*) "dfi_radar_latent_heat_time_period=",  &
                                    dfi_radar_latent_heat_time_period
      write(*,*) "convection_refl_threshold=",convection_refl_threshold
      write(*,*) "l_tten_for_convection_only=",l_tten_for_convection_only
      write(*,*) "fv3_io_layout_y=",fv3_io_layout_y
      write(*,*) "timelevel=",timelevel
    endif
  endif
!  mype=0
  if(mype==0) write(6,*) 'total cores for this run is ',npe
  if(npe /= fv3_io_layout_y) then
     write(6,*) 'ERROR, this run must use ', fv3_io_layout_y,' cores !!!'
     call MPI_FINALIZE(ierror)
     stop 1234
  endif

  mypeLocal=mype+1
  if(mypeLocal <= fv3_io_layout_y) then
!  standard output file for each core
     write(radarfile,'(a,I2.2)') 'stdout_refltotten.',mype
     open(6, file=trim(radarfile),form='formatted',status='unknown')
     write(6,*) '===> deal with subdomain = ', mype
!
! 2.0  open and read background dimesion
!          
     rfv3io_mype=mype
     fv3sar_bg_opt=0    ! 0=restart, 1=input
     if(fv3_io_layout_y==1) then
        call bg_fv3regfilenameg%init
     else
        write(gridspec,'(a,I4.4)') 'fv3_grid_spec.',mype
        write(dynvars,'(a,I4.4)') 'fv3_dynvars.',mype
        write(tracers,'(a,I4.4)') 'fv3_tracer.',mype
        write(sfcvars,'(a,I4.4)') 'fv3_sfcdata.',mype
        write(phyvars,'(a,I4.4)') 'fv3_phydata.',mype
        call bg_fv3regfilenameg%init(grid_spec_input=trim(gridspec), &
                 dynvars_input=trim(dynvars),tracers_input=trim(tracers),&
                 sfcdata_input=trim(sfcvars),phydata_input=trim(phyvars))
     endif
     mype_t=mype
     mype_q=mype
     mype_p=mype
     mype_ql=mype
     mype_2d=mype

     dynvars= trim(bg_fv3regfilenameg%dynvars)
     tracers= trim(bg_fv3regfilenameg%tracers)
     sfcvars= trim(bg_fv3regfilenameg%sfcdata)
     phyvars= trim(bg_fv3regfilenameg%phydata)
     write(6,*) trim(dynvars)
     write(6,*) trim(tracers)
     write(6,*) trim(sfcvars)
     write(6,*) trim(phyvars)

     call init_constants(.true.)
     call init_constants_derived
     krad_bot=7.0_r_single
     Nmsclvl_radar = -999
     iunit_radar=25
     iunit_lightning=26
     iunit_satcast=27
     rad_missing=1
     lightning_missing=1
!
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

     allocate(fv3_io_layout_end(fv3_io_layout_y))
     fv3_io_layout_end=nlat_regional
     if(fv3_io_layout_y > 1) then
       allocate(io_layout_tmp(fv3_io_layout_y))
       io_layout_tmp=0
       fv3_io_layout_end=0
       fv3_io_layout_end(mypeLocal)=nlat_regional 
       call mpi_reduce(fv3_io_layout_end,io_layout_tmp, fv3_io_layout_y,&
                      mpi_integer, mpi_max, 0,MPI_COMM_WORLD, ierror)
       if(mype==0) then
         fv3_io_layout_end(1)=io_layout_tmp(1)
         do i=2,fv3_io_layout_y
           fv3_io_layout_end(i)=fv3_io_layout_end(i-1) + io_layout_tmp(i)
         enddo
       endif
       deallocate(io_layout_tmp)
       call mpi_bcast(fv3_io_layout_end,fv3_io_layout_y,mpi_integer,0,MPI_COMM_WORLD, ierror)
       if(mype==0) write(*,*) 'end of each subdomain ',mype,fv3_io_layout_end
     endif

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

     deallocate(q_bk)
     deallocate(ps_bk)

!     allocate(sat_ctp(nlon_regional,nlat_regional))
!     sat_ctp=miss_obs_real
     do ntime=1,timelevel

        n=ntime

        if(fv3_io_layout_y==1) then
           write(radarfile,'(a,I2.2)') 'RefInGSI3D.dat_',n
        else
           write(radarfile,'(a,I4.4,a1,I2.2)') 'RefInGSI3D.dat.',mype,'_',n
        endif
        write(6,*)
        write(6,*) 'processing ',trim(radarfile)

        inquire(file=trim(radarfile), exist=ifexist)
        if (ifexist) then
           open(iunit_radar,file=trim(radarfile),form='unformatted',status='old')
              read(iunit_radar) Nmsclvl_radar,nlon_radar,nlat_radar
              allocate(ref_mosaic31(nlon_regional,nlat_regional,Nmsclvl_radar))
              read(iunit_radar) ref_mosaic31
           close(iunit_radar)
           write(6,*) 'Nmsclvl_radar,nlon_radar,nlat_radar',  &
                    Nmsclvl_radar,nlon_radar,nlat_radar
           do k=1,Nmsclvl_radar
              write(6,*) 'ref_mosaic31=',k,maxval(ref_mosaic31(:,:,k)), &
                                        minval(ref_mosaic31(:,:,k))
           enddo
           rad_missing=0
        endif
        if (.not.ifexist) then
           write(6,*) 'WARNING: RADAR FILE MISSING:', radarfile
           rad_missing=1
        endif
!
! Read in lightning data
!
        write(lightningfile,'(a,I2.2)') 'LightningInGSI.dat_',n
        write(6,*)
        write(6,*) 'processing ',trim(lightningfile)

        inquire(file=trim(lightningfile), exist=ifexist)
        if (ifexist) then
           open(iunit_lightning,file=trim(lightningfile),form='unformatted',status='old')
              read(iunit_lightning) header1,nlon_lightning,nlat_lightning,numlight,header2,header3
              allocate(lightning_in(header1,numlight))
              lightning_in=-9999.0_r_single
              read(iunit_lightning) lightning_in
           close(iunit_lightning)
           write(6,*) 'finished read ',trim(lightningfile), numlight
           allocate(lightning(nlon_regional,nlat_regional))
           lightning=-9999.0_r_single
           ybegin=1
           yend=nlat_regional
           if(fv3_io_layout_y > 1) then
              ybegin=1
              if(mypeLocal>1) ybegin=fv3_io_layout_end(mypeLocal-1)+1
              yend=fv3_io_layout_end(mypeLocal)
           endif
           call read_Lightning2cld(nlon_regional,nlat_regional,header1,numlight,ybegin,yend,lightning_in,lightning)
           deallocate(lightning_in)
           lightning_missing=0
        endif
        if(.not.ifexist) then
           write(6,*) 'WARNING: LIGHTNING FILE MISSING:', lightningfile
           lightning_missing=1
        endif
!
!  2.6 vertical interpolation of radar reflectivity
!
        allocate(ref_mos_3d(nlon_regional,nlat_regional,nsig_regional))
        ref_mos_3d=miss_obs_real
!    EJ: only call this part if radar data are not missing
        if ( rad_missing == 0 ) then
           call vinterp_radar_ref(nlon_regional,nlat_regional,nsig_regional,Nmsclvl_radar, &
                          ref_mos_3d,ref_mosaic31,h_bk,zh)
           deallocate( ref_mosaic31 )
        endif
        do k=1,nsig_regional
           write(6,*) 'vinterp ref_mos_3d=',k,maxval(ref_mos_3d(:,:,k)), &
                                              minval(ref_mos_3d(:,:,k))
        enddo

        call build_missing_REFcone(nlon_regional,nlat_regional,nsig_regional, &
                                krad_bot,ref_mos_3d,h_bk,pblh)
        do k=1,nsig_regional
           write(6,*) 'refcon ref_mos_3d=',k,maxval(ref_mos_3d(:,:,k)), &
                                             minval(ref_mos_3d(:,:,k))
        enddo

!
! Convert lightning flash rate to reflectivities 
!
!
        if(lightning_missing == 0) then
           write(6,*) 'calling convert_lghtn2ref'
!    
           call convert_lghtn2ref(mype,nlon_regional,nlat_regional,nsig_regional,ref_mos_3d,lightning,h_bk)
           deallocate( lightning )

           write(6,*) 'done adding lightning data'
           do k=1,nsig_regional
              write(6,*) 'lightning',k,' ref_mos_3d=',k,maxval(ref_mos_3d(:,:,k)), &
                                                     minval(ref_mos_3d(:,:,k))
           enddo
        endif
!
!
!     write(6,*) 'calling convert_stcst2ref'
!
!  convert satcast cooling rates to reflectivities
!
!     call convert_stcst2ref(nlon_regional,nlat_regional,nsig_regional, ref_mos_3d,satcast_cr,h_bk)
!     deallocate( satcast_cr )
!     write(6,*) 'done adding satcast data'
 
!     do k=1,nsig_regional
!       write(6,*) 'satcast',k,' ref_mos_3d=',k,maxval(ref_mos_3d(:,:,k)), &
!                                               minval(ref_mos_3d(:,:,k))
!     enddo

!
! 4.10 radar temperature tendency for DFI
!
        allocate(ges_tten(nlon_regional,nlat_regional,nsig_regional))
        ges_tten=-20.0_r_kind
        ges_tten(:,:,nsig_regional)=-10.0_r_kind

        dfi_lhtp=dfi_radar_latent_heat_time_period
        call radar_ref2tten(nlon_regional,nlat_regional,nsig_regional,ref_mos_3d,& 
                     p_bk,t_bk,ges_tten,dfi_lhtp,krad_bot,pblh,  &
                     l_tten_for_convection_only,convection_refl_threshold, &
                     fv3_io_layout_y,fv3_io_layout_end,mype)
! for debug only     if(ntime==2) ges_tten=ref_mos_3d
        deallocate(ref_mos_3d)
        do k=1,nsig_regional
           write(6,*) 'ges_tten=',k,maxval(ges_tten(:,:,k)), &
                                 minval(ges_tten(:,:,k))
        enddo
        if(.not.l_convection_suppress) ges_tten(:,:,nsig_regional)=ges_tten(:,:,nsig_regional-1)

!
! 5.10 update
!

        if(ntime==1) call gsi_fv3ncdf_append(phyvars,'radar_tten',ges_tten,mype_t)
        if(ntime==2) call gsi_fv3ncdf_append(phyvars,'radar_tten_2',ges_tten,mype_t)
        if(ntime==3) call gsi_fv3ncdf_append(phyvars,'radar_tten_3',ges_tten,mype_t)
        if(ntime==4) call gsi_fv3ncdf_append(phyvars,'radar_tten_4',ges_tten,mype_t)

! release memory for this time level
        deallocate(ges_tten)

     enddo  ! ntime
!
!  release memory
!
     write(6,*) 'core', mype ,',finished, now release memory'
!  deallocate(sat_ctp)
     deallocate(h_bk)
     deallocate(t_bk)
     deallocate(p_bk)
     deallocate(pblh)
     deallocate(zh)

  endif ! mypeLocal <= maxcores

  call MPI_Barrier(mpi_comm_world, ierror)
  close(6)
  if(mype==0)  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="
  call MPI_FINALIZE(ierror)

end program gsdcloudanalysis_ref2tten

SUBROUTINE wrf_debug( level , str )
!  USE module_wrf_error
  IMPLICIT NONE
  CHARACTER*(*) str
  INTEGER , INTENT (IN) :: level
  INTEGER               :: debug_level
  CHARACTER (LEN=256) :: time_str
  CHARACTER (LEN=256) :: grid_str
  CHARACTER (LEN=512) :: out_str
!  CALL get_wrf_debug_level( debug_level )
  IF ( level .LE. debug_level ) THEN
    ! old behavior
!      CALL wrf_message( str )
  ENDIF
  write(*,*) 'wrf_debug called !'
  RETURN
END SUBROUTINE wrf_debug

