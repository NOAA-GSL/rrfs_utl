subroutine update_snowice_fv3lam(snowiceRR,xland,nlon,nlat,id,fv3_io_layout_y)
!$$$  documentation block
!                .      .    .                                       .
!   update_snowice_fv3lam: update snow ice based on ims observations
!   prgmmr: Ming Hu                 date: 2022-02-15
!
! program history log:
!
! 2009-07-27: make consistent of all land surface parameters and using fraction seac ice
! if ( fractional_seaice == 0 ) then
!    xice_threshold = 0.5
! else if ( fractional_seaice == 1 ) then
!    xice_threshold = 0.02
! endif
!
!
! For sea-ice:
!                     ivgtyp(i,j)=15     int IVGTYP / =15 for MODIS
!                     landmask(i,j)=2.   float LANDMASK
!                     isltyp(i,j)=16.    int ISLTYP
! 
! For water:
!                     ivgtyp(i,j)=17  / =17 and 21 for MODIS
!                     landmask(i,j)=0.
!                     isltyp(i,j)=14.
! 
 
!   input argument list:
!       snowRR: snow seaice  cover
!       nlon:  x dimension
!       nlat:  Y dimension
!       id:   subdomain id
!       fv3_io_layout_y: total subdomain number
!
! attributes:
!   language: f90
!
!$$$

  use kinds, only: r_single,i_kind
  use module_ncio, only: ncio
  implicit none

  type(ncio) :: fv3grid
!
  integer,intent(in) :: nlon, nlat
  real,intent(in)    :: snowiceRR(nlon,nlat)
  real,intent(in)    :: xland(nlon,nlat)
  integer,intent(in) :: id,fv3_io_layout_y

! Declare local parameters

  character(len=120) :: flnm1
  character (len=5) :: lutype
  integer(i_kind) :: l, n
  
  real(8),allocatable :: tmp8b3d(:,:,:),tmp8b2d(:,:)
  real,allocatable :: tmp4b3d(:,:,:),tmp4b2d(:,:)

  character (len=31) :: rmse_var
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
  integer(i_kind) nlon_regional,nlat_regional,nsig_regional
  integer(i_kind) nsig_soil_regional
  real(r_single),allocatable::precip(:,:)
  real(r_single),allocatable::surftemp(:,:)
  real(r_single),allocatable::tskin(:,:)
  real(r_single),allocatable::tsnow(:,:)
  real(r_single),allocatable::soilmoisture(:,:,:)
  real(r_single),allocatable::soiltemp(:,:,:)

  real(r_single),allocatable::snow(:,:)
  real(r_single),allocatable::snowh(:,:)
  real(r_single),allocatable::snowc(:,:)
  real(r_single),allocatable::seaice(:,:)

  real(r_single),allocatable::snowRRbk(:,:)
  real(r_single),allocatable::snowhRRbk(:,:)
  real(r_single),allocatable::snowcRRbk(:,:)
  real(r_single),allocatable::tskinRRbk(:,:)
  real(r_single),allocatable::tsnowRRbk(:,:)
  real(r_single),allocatable::soiltempRRbk(:,:,:)
  real(r_single),allocatable::surftempRRbk(:,:)
!  surface parameters
  real(r_single),allocatable::landmask(:,:)
!
!
  integer(i_kind)   :: halo
  real(r_single)    :: xice_threshold
  integer(i_kind)   :: fractional_seaice
!
  integer(i_kind)   :: num_seaice2water, num_water2seaice
  integer(i_kind)   :: numtrimsnow, numbuildsnow, numusetrim
  integer   :: MSLANDID
!
  real snowtr,snowhtr,snowctr,snowav,snowhav,snowcav,tskinav,tsnowav,surftempav,soilt1av,soilt2av,soilt3av
  real snowsum,snowhsum,snowcsum,tskinsum,tsnowsum,surftempsum,soilt1sum,soilt2sum,soilt3sum
  real rhosn, snowtrimsum, snowbuiltsum
  integer nsoil,ii,jj,itr,jtr,ist,iend,jst,jend,numnb, numbuildmin
  integer :: nxlocal,nylocal
!
  integer :: i,j,k
  integer :: iii,jjj,num,numh,i4,j4
  real    :: newvalue, newvalueh
!
! =============================================================
!
  halo=2
  fractional_seaice=1 ! should always be =1 in RRFS
  if ( fractional_seaice == 0 ) then
    xice_threshold = 0.5
    write(*,*) ' do not use fraction sea ice'
  else if ( fractional_seaice == 1 ) then
    xice_threshold = 0.02 ! default RRFS min_seaice=1.0e-11, default WRF = 0.02
    write(*,*) ' use fraction sea ice'
  endif
!  
  open(12,file="coupler.res",form='formatted')
     read(12,*)
     read(12,*)
     read(12,*) iyear,imonth,iday,ihour,iminute,isecond
     write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond
  close(12)
!
!-- calculate snow and graupal precipiation over surface
!  (note: remove liquid precip for snow trimming to work correctly)

  if(fv3_io_layout_y > 1 ) then
    write(flnm1,'(a,I4.4)') 'fv_tracer.res.tile1.nc.',id-1
  else
    flnm1='fv_tracer.res.tile1.nc'
  endif 
  call fv3grid%open(trim(flnm1),"r",200)
  call fv3grid%get_dim("xaxis_1",nlon_regional)
  call fv3grid%get_dim("yaxis_1",nlat_regional)
  call fv3grid%get_dim("zaxis_1",nsig_regional)
  if(nlon_regional/=nlon .or. nlat_regional/=nlat) then
     write(6,*) 'Wrong dimension=',nlon_regional,nlat_regional,nlon,nlat
     call fv3grid%close() 
     stop 1234
  endif

  allocate(tmp4b3d(nlon_regional,nlat_regional,nsig_regional))
  allocate(precip(nlon_regional,nlat_regional))
  precip=0.0
  call fv3grid%get_var("snowwat",nlon,nlat,nsig_regional,tmp4b3d)
  precip(:,:)=tmp4b3d(:,:,nsig_regional)
  call fv3grid%get_var("graupel",nlon,nlat,nsig_regional,tmp4b3d)
  precip(:,:)=precip(:,:)+tmp4b3d(:,:,nsig_regional)
  deallocate(tmp4b3d)
  call fv3grid%close() 
  write(6,*) 'precipiation on ground=',maxval(precip),minval(precip)

  allocate(surftemp(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  surftemp=0.0
  if(1==1) then   ! use 1st level atmosphere temperature
     if(fv3_io_layout_y > 1 ) then
        write(flnm1,'(a,I4.4)') 'fv_core.res.tile1.nc.',id-1
     else
        flnm1='fv_core.res.tile1.nc'
     endif
     call fv3grid%open(trim(flnm1),"r",200)
     call fv3grid%get_dim("xaxis_1",nlon_regional)
     call fv3grid%get_dim("yaxis_2",nlat_regional)
     call fv3grid%get_dim("zaxis_1",nsig_regional)
     if(nlon_regional/=nlon .or. nlat_regional/=nlat) then
        write(6,*) 'Wrong dimension=',nlon_regional,nlat_regional,nlon,nlat
        call fv3grid%close() 
        stop 1234
     endif

     allocate(tmp4b3d(nlon_regional,nlat_regional,nsig_regional))
     call fv3grid%get_var("T",nlon,nlat,nsig_regional,tmp4b3d)
     surftemp(1:nlon,1:nlat)=tmp4b3d(1:nlon,1:nlat,nsig_regional)
     deallocate(tmp4b3d)
     call fv3grid%close()
     write(6,*) 'surface temperature =',maxval(surftemp),minval(surftemp)
  endif

  if(fv3_io_layout_y > 1 ) then
     write(flnm1,'(a,I4.4)') 'sfc_data.nc.',id-1
  else
     flnm1='sfc_data.nc'
  endif
  call fv3grid%open(trim(flnm1),"r",200)
  call fv3grid%get_dim("xaxis_1",nlon_regional)
  call fv3grid%get_dim("yaxis_1",nlat_regional)
  call fv3grid%get_dim("zaxis_1",nsig_soil_regional)
  if(nlon_regional/=nlon .or. nlat_regional/=nlat) then
     write(6,*) 'Wrong dimension=',nlon_regional,nlat_regional,nlon,nlat
     call fv3grid%close() 
     stop 1234
  endif

  allocate(tmp8b3d(nlon_regional,nlat_regional,nsig_soil_regional))
  allocate(tmp8b2d(nlon_regional,nlat_regional))

  write(6,*)' nlon,lat,sig_regional=',nlon_regional,nlat_regional,nsig_regional

  allocate(snow(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(snowh(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(snowc(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(seaice(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(tskin(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(tsnow(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(landmask(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  snow=0.0
  snowh=0.0
  snowc=0.0
  seaice=0.0
  tskin=0.0
  tsnow=0.0
  landmask=0
  
!  rmse_var='TSK'
  call fv3grid%get_var("tsfcl",nlon,nlat,tmp8b2d)
  tskin(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)
  write(6,*)' max,min skin temp (K)=',maxval(tskin),minval(tskin)
!
!  rmse_var='SOILT1'
  call fv3grid%get_var("tsnow_land",nlon,nlat,tmp8b2d)
  tsnow(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)
  write(6,*)' max,min snow temp (K)=',maxval(tsnow),minval(tsnow)
!
!  rmse_var='SNOW' [mm]
  call fv3grid%get_var("weasdl",nlon,nlat,tmp8b2d)
  snow(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)
  write(6,*)' max,min SNOW=',maxval(snow),minval(snow)
!
!  rmse_var='SNOWH'
!  call fv3grid%get_var("snwdph",nlon,nlat,tmp8b2d)
!??? not used in this code???  snowdp=tmp8b2d(:,:)
!  write(6,*)' max,min SNOWH=',maxval(snowh),minval(snowh)
!  rmse_var='SNOWH' [mm]
  call fv3grid%get_var("snodl",nlon,nlat,tmp8b2d)
  snowh(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)*1.e-3 ! convert to [m]
  write(6,*)' max,min SNOWH=',maxval(snowh),minval(snowh)
!
!  rmse_var='SNOWC' [fraction]
  call fv3grid%get_var("sncovr",nlon,nlat,tmp8b2d)
  snowc(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)
  write(6,*)' max,min SNOWC=',maxval(snowc),minval(snowc)
!
!  rmse_var='SEAICE'
  call fv3grid%get_var("fice",nlon,nlat,tmp8b2d)
  seaice(1:nlon_regional,1:nlat_regional)=tmp8b2d(:,:)
  write(6,*)' max,min SEAICE=',maxval(seaice),minval(seaice)
!
  allocate(soiltemp(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo,nsig_soil_regional))
  allocate(soiltempRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo,nsig_soil_regional))
!  rmse_var='TSLB'
  call fv3grid%get_var("tslb",nlon,nlat,nsig_soil_regional,tmp8b3d)
  do k=1,nsig_soil_regional
     soiltemp(1:nlon_regional,1:nlat_regional,k)=tmp8b3d(:,:,k)
     write(6,*)' max,min TSLB=',k, maxval(soiltemp(:,:,k)),minval(soiltemp(:,:,k))
  enddo
!
  deallocate(tmp8b3d)
!
!  rmse_var='LANDMASK' 0 - water, 1 - land, 2 - ice
  call fv3grid%get_var("slmsk",nlon,nlat,tmp8b2d)
  landmask(1:nlon_regional,1:nlat_regional)=tmp8b2d
  write(6,*)' max,min LANDMASK=',maxval(landmask),minval(landmask)
!
  deallocate(tmp8b2d)
  call fv3grid%close()
!
! fill halo first 
!
  do j=1,halo
     do i=1,nlon_regional
        tskin(i,1-j)=tskin(i,1)
        tsnow(i,1-j)=tsnow(i,1)
        snow(i,1-j)=snow(i,1)
        snowh(i,1-j)=snowh(i,1)
        snowc(i,1-j)=snowc(i,1)
        seaice(i,1-j)=seaice(i,1)
        soiltemp(i,1-j,:)=soiltemp(i,1,:)
        landmask(i,1-j)=landmask(i,1)
        surftemp(i,1-j)=surftemp(i,1)
        tskin(i,nlat_regional+j)=tskin(i,nlat_regional)
        tsnow(i,nlat_regional+j)=tsnow(i,nlat_regional)
        snow(i,nlat_regional+j)=snow(i,nlat_regional)
        snowh(i,nlat_regional+j)=snowh(i,nlat_regional)
        snowc(i,nlat_regional+j)=snowc(i,nlat_regional)
        seaice(i,nlat_regional+j)=seaice(i,nlat_regional)
        soiltemp(i,nlat_regional+j,:)=soiltemp(i,nlat_regional,:)
        landmask(i,nlat_regional+j)=landmask(i,nlat_regional)
        surftemp(i,nlat_regional+j)=surftemp(i,nlat_regional)
     enddo
  enddo

  if(fv3_io_layout_y > 1 ) then
! fill in north halo

     if(id <= fv3_io_layout_y-2) then


        if(1==1) then   ! use 1st level atmosphere temperature
           write(flnm1,'(a,I4.4)') 'fv_core.res.tile1.nc.',id
           call fv3grid%open(trim(flnm1),"r",200)
           call fv3grid%get_dim("xaxis_1",nxlocal)
           call fv3grid%get_dim("yaxis_2",nylocal)
           if(nlon_regional/=nxlocal) then
              write(6,*) 'Wrong dimension=',nlon_regional,nxlocal
              call fv3grid%close()
              stop 1234
           endif

           allocate(tmp4b3d(nxlocal,nylocal,nsig_regional))
           call fv3grid%get_var("T",nxlocal,nylocal,nsig_regional,tmp4b3d)
           surftemp(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp4b3d(:,1:halo,nsig_regional)
           deallocate(tmp4b3d)
           call fv3grid%close()
        endif

        write(flnm1,'(a,I4.4)') 'sfc_data.nc.',id
        call fv3grid%open(trim(flnm1),"r",200)
        call fv3grid%get_dim("xaxis_1",nxlocal)
        call fv3grid%get_dim("yaxis_1",nylocal)
        if(nlon_regional/=nxlocal) then
           write(6,*) 'Wrong dimension=',nlon_regional,nxlocal
           call fv3grid%close()
           stop 1234
        endif

        allocate(tmp8b3d(nxlocal,nylocal,nsig_soil_regional))
        allocate(tmp8b2d(nxlocal,nylocal))

        call fv3grid%get_var("tsfcl",nxlocal,nylocal,tmp8b2d)
        tskin(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("tsnow_land",nxlocal,nylocal,tmp8b2d)
        tsnow(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("weasdl",nxlocal,nylocal,tmp8b2d)
        snow(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("snodl",nxlocal,nylocal,tmp8b2d)
        snowh(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("sncovr",nxlocal,nylocal,tmp8b2d)
        snowc(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("fice",nxlocal,nylocal,tmp8b2d)
        seaice(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)
!
        call fv3grid%get_var("tslb",nxlocal,nylocal,nsig_soil_regional,tmp8b3d)
        do k=1,nsig_soil_regional
           soiltemp(1:nlon_regional,nlat_regional+1:nlat_regional+halo,k)=tmp8b3d(:,1:halo,k)
        enddo
!
        deallocate(tmp8b3d)
!
        call fv3grid%get_var("slmsk",nxlocal,nylocal,tmp8b2d)
        landmask(1:nlon_regional,nlat_regional+1:nlat_regional+halo)=tmp8b2d(:,1:halo)

        deallocate(tmp8b2d)
        call fv3grid%close()
     endif

! fill in south halo
     if(id >= 2) then

        if(1==1) then   ! use 1st level atmosphere temperature
           write(flnm1,'(a,I4.4)') 'fv_core.res.tile1.nc.',id-2
           call fv3grid%open(trim(flnm1),"r",200)
           call fv3grid%get_dim("xaxis_1",nxlocal)
           call fv3grid%get_dim("yaxis_2",nylocal)
           if(nlon_regional/=nxlocal) then
              write(6,*) 'Wrong dimension=',nlon_regional,nxlocal
              call fv3grid%close()
              stop 1234
           endif

           allocate(tmp4b3d(nxlocal,nylocal,nsig_regional))
           call fv3grid%get_var("T",nxlocal,nylocal,nsig_regional,tmp4b3d)
           surftemp(1:nlon_regional,1-halo:0)=tmp4b3d(:,nylocal-halo+1:nylocal,nsig_regional)
           deallocate(tmp4b3d)
           call fv3grid%close()
        endif

        write(flnm1,'(a,I4.4)') 'sfc_data.nc.',id-2
        call fv3grid%open(trim(flnm1),"r",200)
        call fv3grid%get_dim("xaxis_1",nxlocal)
        call fv3grid%get_dim("yaxis_1",nylocal)
        if(nlon_regional/=nxlocal) then
           write(6,*) 'Wrong dimension=',nlon_regional,nxlocal
           call fv3grid%close()
           stop 1234
        endif

        allocate(tmp8b3d(nxlocal,nylocal,nsig_soil_regional))
        allocate(tmp8b2d(nxlocal,nylocal))

        call fv3grid%get_var("tsfcl",nxlocal,nylocal,tmp8b2d)
        tskin(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("tsnow_land",nxlocal,nylocal,tmp8b2d)
        tsnow(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("weasdl",nxlocal,nylocal,tmp8b2d)
        snow(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("snodl",nxlocal,nylocal,tmp8b2d)
        snowh(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("sncovr",nxlocal,nylocal,tmp8b2d)
        snowc(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("fice",nxlocal,nylocal,tmp8b2d)
        seaice(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)
!
        call fv3grid%get_var("tslb",nxlocal,nylocal,nsig_soil_regional,tmp8b3d)
        do k=1,nsig_soil_regional
           soiltemp(1:nlon_regional,1-halo:0,k)=tmp8b3d(:,nylocal-halo+1:nylocal,k)
        enddo
!
        deallocate(tmp8b3d)
!
        call fv3grid%get_var("slmsk",nxlocal,nylocal,tmp8b2d)
        landmask(1:nlon_regional,1-halo:0)=tmp8b2d(:,nylocal-halo+1:nylocal)

        deallocate(tmp8b2d)
        call fv3grid%close()
     endif

  endif

  do i=1,halo
     do j=1-halo,nlat_regional+halo
        tskin(1-i,j)=tskin(1,j)
        tsnow(1-i,j)=tsnow(1,j)
        snow(1-i,j)=snow(1,j)
        snowh(1-i,j)=snowh(1,j)
        snowc(1-i,j)=snowc(1,j)
        seaice(1-i,j)=seaice(1,j)
        soiltemp(1-i,j,:)=soiltemp(1,j,:)
        landmask(1-i,j)=landmask(1,j)
        surftemp(1-i,j)=surftemp(1,j)
        tskin(nlon_regional+i,j)=tskin(nlon_regional,j)
        tsnow(nlon_regional+i,j)=tsnow(nlon_regional,j)
        snow(nlon_regional+i,j)=snow(nlon_regional,j)
        snowh(nlon_regional+i,j)=snowh(nlon_regional,j)
        snowc(nlon_regional+i,j)=snowc(nlon_regional,j)
        seaice(nlon_regional+i,j)=seaice(nlon_regional,j)
        soiltemp(nlon_regional+i,j,:)=soiltemp(nlon_regional,j,:)
        landmask(nlon_regional+i,j)=landmask(nlon_regional,j)
        surftemp(nlon_regional+i,j)=surftemp(nlon_regional,j)
     enddo
  enddo
!
!
  write(6,*) '================================================='
!
! save the RR background snow in snowRRbk
!

  allocate(snowRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(snowhRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(snowcRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(tskinRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(tsnowRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))
  allocate(surftempRRbk(1-halo:nlon_regional+halo,1-halo:nlat_regional+halo))

  snowRRbk=snow
  snowhRRbk=snowh
  snowcRRbk=snowc
  tskinRRbk=tskin
  tsnowRRbk=tsnow
  soiltempRRbk=soiltemp
  surftempRRbk=surftemp
!
!  trim snow
!

  snowtrimsum=0.
  snowbuiltsum=0.

  numtrimsnow=0
  numbuildsnow=0
  numusetrim=0
  numbuildmin=0
  itr = 0
  jtr = 0

  DO J=1,nlat
  DO I=1,nlon
    !-- check on snow temperature
    if(snowc(i,j) > 0.) then
      tsnow(i,j)=min(273.15,0.5*(tskin(i,j)+soiltemp(i,j,1)))
    else
      tsnow(i,j)=tskin(i,j) 
    endif
          !print *, 'Temperature inside snow is corrected ', i,j, tsnow(i,j),snowc(i,j)

    if(landmask(i,j) == 1 ) then  ! on land
      if(snowiceRR(i,j) < 1.0e-12 .and. snow(i,j) > 0.0 ) then 
      !-- No snow in IMS but there is snow in RRFS
          write(6,*) 'trim snow',i,j,snow(i,j),snowh(i,j),precip(i,j),surftemp(i,j),snowiceRR(i,j) 
        if(snow(i,j) > 1.)then 
        !-- use trimmed snow only when it is > 1 mm
          numtrimsnow=numtrimsnow+1
          itr=i
          jtr=j

          ! save values of snow to be trimmed
          snowtr=snow(i,j)
          snowhtr=snowh(i,j)
          snowctr=snowc(i,j)
          snowtrimsum=snowtrimsum+snow(i,j)
        endif

          ! trim snow
          snow(i,j) = 0.0_4
          snowh(i,j) = 0.0_4
          snowc(i,j) = 0.0_4
      endif

      !tgs snow building
      if(snowiceRR(i,j) > 1.0e-12 .and. snow(i,j) == 0.0 ) then   !  underforecasted snow
      ! Snow in IMS, but no snow in RRFS
        if(surftemp(i,j) < 278.0 ) then   
        !-- build only if T < 278K
           write(6,*) 'build snow at i,j',i,j,'precip,surftemp,snowiceRR',precip(i,j),surftemp(i,j),snowiceRR(i,j)

           snowsum = 0.
           snowhsum = 0.
           snowcsum = 0.
           tskinsum = 0.
           tsnowsum = 0.
           soilt1sum = 0.
           soilt2sum = 0.
           soilt3sum = 0.
           surftempsum = 0.

           numnb=0
           ist=i-2
           iend=i+2
           jst=j-2
           jend=j+2
           do ii=ist,iend
           do jj=jst,jend
             if(landmask(ii,jj) == 1) then  ! land
               if(ii.eq.itr.and.jj.eq.jtr) then 
               !-- snow trimmed at the neighbor point
                 print *,'trimmed snow at ii,jj', ii,jj,itr,jtr
                 numnb=100
               endif
               if( numnb== 100) exit

               if(snowRRbk(ii,jj) > 1.) then ! swe > 1 mm
                 numnb = numnb + 1
                 snowsum = snowsum + snowRRbk(ii,jj)
                 snowhsum = snowhsum + snowhRRbk(ii,jj)
                 snowcsum = snowcsum + snowcRRbk(ii,jj) 
                 tskinsum = tskinsum + tskinRRbk(ii,jj)
                 tsnowsum = tsnowsum + tsnowRRbk(ii,jj)
                 soilt1sum = soilt1sum + soiltempRRbk(ii,jj,1)
                 soilt2sum = soilt2sum + soiltempRRbk(ii,jj,2)
                 soilt3sum = soilt3sum + soiltempRRbk(ii,jj,3)
                 !write(*,*) id,ii,jj
                 !write(*,*) surftempRRbk(ii,jj)
                 surftempsum = surftempsum + surftempRRbk(ii,jj)
               endif
             endif
           enddo
             if( numnb == 100) exit
           enddo

           !-- compute averages for all neighbor land points
           if( (numnb.ge.1) .and. (numnb .ne. 100)) then
           !-- at least one neighbor with snow and no points with trimmed snow
             snowav=snowsum/numnb
             snowhav=snowhsum/numnb
             snowcav=snowcsum/numnb
             tskinav=tskinsum/numnb
             tsnowav=tsnowsum/numnb
             soilt1av=soilt1sum/numnb
             soilt2av=soilt2sum/numnb
             soilt3av=soilt3sum/numnb
             surftempav=surftempsum/numnb
             write(6,*) 'snow neighbors found, numnb =',numnb, &
               'snowsum,snowav,snowhav,snowcav,tskinav,tsnowav,soilt1av,soilt2av,soilt3av,surftempav', &
                snowsum,snowav,snowhav,snowcav,tskinav,tsnowav,soilt1av,soilt2av,soilt3av,surftempav
           endif

           numbuildsnow=numbuildsnow+1
           if(numnb == 100) then ! use point with trimmed snow
             if(snowhtr > 1.e-12) then
               numusetrim=numusetrim+1
               write(6,*) 'trimmed snow at itr,jtr',itr,jtr,'is used to build snow at point i,j',i,j
               write(6,*) 'snowtr, snowhtr, snowctr, tskin(itr,jtr), tsnow(itr,jtr)', &
                           snowtr, snowhtr, snowctr,tskin(itr,jtr),tsnow(itr,jtr)
               rhosn=max(58.8,min(500.,snowtr/snowhtr))
               snow(i,j) = max(1.,snowtr) ! not less than 1 mm SWE
               snowh(i,j) = snow(i,j)/rhosn
               snowc(i,j) = min(1.,snow(i,j)/32.)
               tskin(i,j) = tskin(itr,jtr)
               tsnow(i,j) = min(tsnow(itr,jtr),272.)
               soiltemp(i,j,1) = min(soiltemp(itr,jtr,1),272.)
               soiltemp(i,j,2) = min(soiltemp(itr,jtr,2),272.5)
               soiltemp(i,j,3) = min(soiltemp(itr,jtr,3),273.)
             else
             !tgs 22apr15 - this warning is OK in the cold-start run
             !-- This warning in the cycled run indicates a problem.
               write(6,*) 'WARNING in snow build from the neighbor-point trimmed snow '
               write(6,*) 'Set snow to min value,j,snowhtr',i,j,snowhtr
               numbuildmin=numbuildmin+1
               snow(i,j) = 1.0
               snowh(i,j) = 1.0/250. ! rhosn=250.,snowh[m]=snow[mm]/rhosn
               snowc(i,j) = min(1.,snow(i,j)/32.) ! snowc=1 if snow=32mm 
               tskin(i,j) = min(tskin(i,j),272.)
               tsnow(i,j) = min(tsnow(i,j),272.)
               soiltemp(i,j,1) = min(soiltemp(i,j,1),272.)
               soiltemp(i,j,2) = min(soiltemp(i,j,2),272.5)
               soiltemp(i,j,3) = min(soiltemp(i,j,3),273.)
             endif
           else

             if(numnb.ge.1) then ! use neighbor's average
               if(snowhav > 1.e-12 .and. snowav > 1.e-12) then
                 write(6,*)'build snow based on neighbor points ',numnb
                 rhosn=max(58.8,min(500.,snowav/snowhav))
                 snow(i,j) = max(1.,snowav)
                 snowh(i,j) = snow(i,j)/rhosn
                 snowc(i,j) = min(1.,snow(i,j)/32.)
                 tskin(i,j) = min(min(tskinav,tskin(i,j)),272.)
                 tsnow(i,j) = min(min(tsnowav,tsnow(i,j)),272.)
                 soiltemp(i,j,1) = min(min(soilt1av,soiltemp(i,j,1)),272.)
                 soiltemp(i,j,2) = min(min(soilt2av,soiltemp(i,j,2)),272.5)
                 soiltemp(i,j,3) = min(min(soilt3av,soiltemp(i,j,3)),273.)
               else
               !tgs 22apr15 - this warning is OK in the cold-start run
               !-- This warning in the cycled run indicates a problem.
                 write(6,*) ' WARNING in snow build from the neighbors average '
                 write(6,*) 'Set snow to min value - i,j,snowhav,rhosn',i,j,snowhav,rhosn
                 numbuildmin=numbuildmin+1
                 snow(i,j) = 1.0
                 snowh(i,j) = 1.0/250. ! rhosn=250.,snowh[m]=snow[mm]/rhosn
                 snowc(i,j) = min(1.,snow(i,j)/32.) ! snowc=1 if snow=32mm 
                 tskin(i,j) = min(tskin(i,j),272.)
                 tsnow(i,j) = min(tsnow(i,j),272.)
                 soiltemp(i,j,1) = min(soiltemp(i,j,1),272.)
                 soiltemp(i,j,2) = min(soiltemp(i,j,2),272.5)
                 soiltemp(i,j,3) = min(soiltemp(i,j,3),273.)
               endif
             else ! no neighbors with snow
               write(6,*) 'set snow to min value - 1mm of SWE'
               numbuildmin=numbuildmin+1
               snow(i,j) = 1.0  
               snowh(i,j) = 1.0/250. ! rhosn=250.,snowh[mm]=snow[mm]*1.e3/rhosn
               snowc(i,j) = min(1.,snow(i,j)/32.) ! snowc=1 if snow=32mm 
               tskin(i,j) = min(tskin(i,j),272.)
               tsnow(i,j) = min(tsnow(i,j),272.)
               soiltemp(i,j,1) = min(soiltemp(i,j,1),272.)
               soiltemp(i,j,2) = min(soiltemp(i,j,2),272.5)
               soiltemp(i,j,3) = min(soiltemp(i,j,3),273.)
             endif
           endif  !  if(numnb == 100) then
           snowbuiltsum=snowbuiltsum+snow(i,j)
           write(6,*) 'BUILD - snow,snowh,snowc,tskin,tsnow,soiltemp1,soiltemp2,soiltemp3', &
              i,j,snow(i,j),snowh(i,j),snowc(i,j),tskin(i,j),tsnow(i,j),soiltemp(i,j,1),soiltemp(i,j,2),soiltemp(i,j,3)       
        endif
      endif
    endif

! limit snow depth not to exceed 5.e4 mm (50 m)
    if((snowh(i,j) >= 0.0_4 .and. snowh(i,j) <=50.) .and. (snow(i,j)  <=20000. .and. snow(i,j)  >=0.0_4) ) then
    elseif(snowh(i,j) < 0.0_4 .or. snow(i,j)  < 0.0_4) then
      snowh(i,j)=0.0_4
      snow(i,j) = 0.0_4
    elseif(snowh(i,j) > 50. .or. snow(i,j)  > 20000.) then
      write(6,*) 'Huge snow value i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
      newvalue=0.0_4
      newvalueh=0.0_4
      num=0
      numh=0
      do jjj=j-1,j+1
        do iii=i-1,i+1
          !write(6,*) iii,jjj,snowh(iii,jjj),snow(iii,jjj)
          if(iii .ne. i .and. jjj .ne. j) then
            newvalue=newvalue+snow(iii,jjj)
            newvalueh=newvalueh+snowh(iii,jjj)
            num=num+1
          endif
        enddo
      enddo
      if(num > 0 .and. newvalue < 100000.0 .and. newvalueh < 2.e5) then
        snow(i,j)=newvalue/num
        snowh(i,j)=newvalueh/num
      else
        snow(i,j)=snow(iii,jjj)
        snowh(i,j)=snowh(iii,jjj)
      endif

      write(6,*)'Corrected snow value i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
    else
      write(6,*) '===>Error<===: strange point i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
      snowh(i,j) = 0.0_4
      snow(i,j)  = 0.0_4
      snowc(i,j) = 0.0_4
    endif
! check consistency of snow variables after snow trim
    if((snow(i,j) <0.0_4 .and. abs(snow(i,j) ) <1.0e-10 .and. snowh(i,j) > 0.0_4)                   &
      .or. (snowh(i,j) <0.0_4 .and. abs(snowh(i,j) ) <1.0e-10 .and.snow(i,j) > 0.0_4)) then
      write(6,*) 'Inconsistency of snow and snowh AFTER snow trim at i,j,snow,snowh', i,j,snow(i,j),snowh(i,j)
      snow(i,j)  = 0.0_4
      snowh(i,j) = 0.0_4
      snowc(i,j) = 0.0_4
      write(6,*) 'Corrected snow and snowh at i,j,snow,snowh',i,j,snow(i,j),snowh(i,j)
    endif 
  ENDDO
  ENDDO

  write(6,*) 'SUMMARY on snow trim/build:'
  write(6,*) 'grid point with trimmed snow: ', numtrimsnow
  write(6,*) 'grid point with built snow: ', numbuildsnow
  write(6,*) 'grid point with built snow from the "trimmed" neighbor: ', numusetrim
  write(6,*) 'grid point with built min (=1 mm) snow: ', numbuildmin
  write(6,*) 'total SWE trimmed:',snowtrimsum,'[mm]'
  write(6,*) 'total SWE built:',snowbuiltsum, '[mm]'
!
!
!  get rid of snow on water
!
  DO J=1,nlat
  DO I=1,nlon
    if( landmask(i,j) == 0 .and. snow(i,j) > 0.0_4) then    ! snow on open water
      snow(i,j) = 0.0_4
      snowh(i,j) = 0.0_4
      snowc(i,j) = 0.0_4
    endif
  ENDDO
  ENDDO

!
  deallocate(snowRRbk)
  deallocate(snowhRRbk)
  deallocate(snowcRRbk)
  deallocate(tskinRRbk)
  deallocate(tsnowRRbk)
  deallocate(surftempRRbk)

!
!           update mass core netcdf file with snow,snowh,snowc
!
  write(6,*) ' ================== '
  write(6,*) ' trim snow and replace ice '
  write(6,*) ' ================== '

  if(fv3_io_layout_y > 1 ) then
     write(flnm1,'(a,I4.4)') 'sfc_data.nc.',id-1
  else
     flnm1='sfc_data.nc.0000'
  endif
  call fv3grid%open(trim(flnm1),"w",200)
  call fv3grid%get_dim("xaxis_1",nlon_regional)
  call fv3grid%get_dim("yaxis_1",nlat_regional)
  call fv3grid%get_dim("zaxis_1",nsig_soil_regional)
  if(nlon_regional/=nlon .or. nlat_regional/=nlat) then
     write(6,*) 'Wrong dimension=',nlon_regional,nlat_regional,nlon,nlat
     call fv3grid%close()
     stop 1234
  endif

  allocate(tmp8b3d(nlon_regional,nlat_regional,nsig_soil_regional))
  allocate(tmp8b2d(nlon_regional,nlat_regional))

!  rmse_var='TSLB'
  do k=1,nsig_soil_regional
     tmp8b3d(:,:,k)=soiltemp(1:nlon_regional,1:nlat_regional,k)
     write(6,*)' max,min TSLB=',k, maxval(tmp8b3d(:,:,k)),minval(tmp8b3d(:,:,k))
  enddo
  call fv3grid%replace_var("tslb",nlon,nlat,nsig_soil_regional,tmp8b3d)
  deallocate(tmp8b3d)
! SNCOVR
  tmp8b2d=snowc(1:nlon_regional,1:nlat_regional)
  do j=1,nlat_regional
     do i =1,nlon_regional
        if(tmp8b2d(i,j)<1.0e-10) tmp8b2d(i,j)=0.0_8
     enddo
  enddo
  write(6,*)' max,min snowc=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("sncovr",nlon,nlat,tmp8b2d)
! 'SNODL'
  tmp8b2d=snowh(1:nlon_regional,1:nlat_regional)*1.e3 ! convert to [mm]
  do j=1,nlat_regional
     do i =1,nlon_regional
        if(tmp8b2d(i,j)<1.0e-10) tmp8b2d(i,j)=0.0_8
     enddo
  enddo
  write(6,*)' max,min snowh=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("snodl",nlon,nlat,tmp8b2d)
! 'SNWDPH'
  tmp8b2d=snowh(1:nlon_regional,1:nlat_regional)*1.e3 ! convert to [mm]
  do j=1,nlat_regional
     do i =1,nlon_regional
        if(tmp8b2d(i,j)<1.0e-10) tmp8b2d(i,j)=0.0_8
     enddo
  enddo
  write(6,*)' max,min snowh=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("snwdph",nlon,nlat,tmp8b2d)
! 'WEASDL'
  tmp8b2d=snow(1:nlon_regional,1:nlat_regional)
  do j=1,nlat_regional
     do i =1,nlon_regional
        if(tmp8b2d(i,j)<1.0e-10) tmp8b2d(i,j)=0.0_8
     enddo
  enddo
  write(6,*)' max,min snow=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("weasdl",nlon,nlat,tmp8b2d)
! TSFC
  tmp8b2d=tskin(1:nlon_regional,1:nlat_regional)
  write(6,*)' max,min tsk=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsfc",nlon,nlat,tmp8b2d)
! TSFCL
  tmp8b2d=tskin(1:nlon_regional,1:nlat_regional)
  write(6,*)' max,min tsk=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsfcl",nlon,nlat,tmp8b2d)
! TSNOW_LAND
  tmp8b2d=tsnow(1:nlon_regional,1:nlat_regional)
  write(6,*)' max,min soilt1=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsnow_land",nlon,nlat,tmp8b2d)
!
  deallocate(tmp8b2d)
  call fv3grid%close()
  
  deallocate(snow)
  deallocate(snowh)
  deallocate(snowc)
  deallocate(seaice)
  deallocate(tskin)
  deallocate(tsnow)

  deallocate(landmask)

  deallocate(soiltemp)
  deallocate(soiltempRRbk)
  deallocate(precip)
  deallocate(surftemp)

end subroutine update_snowice_fv3lam
