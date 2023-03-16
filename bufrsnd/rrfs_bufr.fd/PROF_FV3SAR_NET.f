      SUBROUTINE PROF_FV3SAR_NET(filedyn,filephys,datestr,ITAG,INCR,
     &                           NSTAT)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  PROF        PROFILE SOUNDINGS
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-04-22
C
C ABSTRACT:  THIS ROUTINE GENERATES THE RAW PROFILE SOUNDING
C            OUTPUT FILES FROM THE FORECAST RESTRT FILE AND
C            AUXILIARY FILES
C
C PROGRAM HISTORY LOG:
C   99-04-22  T BLACK - ORIGINATOR
C   02-07-01  G MANIKIN - FIXED PROBLEM WITH DHCNVC AND DHRAIN
C                          COMPUTATIONS - SEE COMMENTS BELOW
C   03-04-01  M PYLE - BEGAN CONVERTING FOR WRF
C   19-10-23  M PYLE - BEGAN CONVERTING FOR FV3SAR NETCDF output
C   20-02-04  M Pyle - Parallelized version reading just the station data to 
C                      improve performance for compressed netCDF output from model.
C   22-07-13  M Pyle - pass in NSTAT limit rather than having as a hard-
C                      coded parameter statement.
C
C USAGE:  CALL PROF FROM PROGRAM POST0
C
C   INPUT ARGUMENT LIST:
C     NHB    - THE UNIT NUMBER FOR READING THE NHB FILE
C     LRSTRT - THE UNIT NUMBER FOR READING THE RESTRT FILE
C     ITAG   - THE FORECAST HOUR WE ARE DEALING WITH
C     LCLAS1 - THE UNIT NUMBER FOR WRITING THE PROFILE DATA
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C
C-----------------------------------------------------------------------
c      use vrbls3d
c      use vrbls2d
c      use soil
c      use masks
C
       use netcdf
      include 'mpif.h'

!      INCLUDE "parmeta"
      INCLUDE "parmsoil"
C-----------------------------------------------------------------------
                             P A R A M E T E R
     & (LCL1ML=15,LCL1SL=52)
!       NWORDM=(LCL1ML+1)*LM+2*LCL1SL
!     &, LRECPR=4*(8+9+LCL1ML*LM+LCL1SL))

C-----------------------------------------------------------------------
C
C    PARMS FOR HOURLY PROFILER OUTPUT
C      NSTAT - MAX NUMBER OF STATIONS
C      NWORDM - DIMENSION OF OUTPUT ARRAY, MUST BE LARGE ENOUGH
C          TO HOLD ALL VARIABLES
C          (MAX NO MULTI-LAYER VARIABLES*LM + NO OF SINGLE LAYER VARS)
C      LCL1ML - NUMBER OF MULTI-LAYER VARIABLES OUTPUT FOR CLASS 1
C      LCL1SL - NUMBER OF SINGLE LAYER VARIABLES OUTPUT FOR CLASS 1
C
C------------------------------------------------------------------------
                             P A R A M E T E R
     & (ITB=76,JTB=134)
                             P A R A M E T E R
     & (A2=17.2693882,A3=273.16,A4=35.86,PQ0=379.90516
     &, DTR=1.74532925E-2
     &, G=9.81,GI=1./G,RD=287.04,CP=1004.6,CAPA=RD/CP,RHCRIT=0.9999
     &, con_rd=287.05,con_rv=461.5,con_fvirt=con_rv/con_rd-1.
     &, spval=1.e+10)


C------------------------------------------------------------------------
      real, allocatable:: STNLAT(:),STNLON(:)

C
       REAL, ALLOCATABLE::
     & RES(:),FIS(:),THS(:),HBOT(:)
     &,CFRACL(:),CFRACM(:),CFRACH(:)
     &,SOILTB(:),SFCEXC(:),SMSTAV(:),SMSTOT(:)
     &,Z0(:),CZEN(:),CZMEAN(:)
     &,ACPREC(:),CUPREC(:),ACSNOW(:),ACSNOM(:)
     &,SSROFF(:),BGROFF(:),SFCSHX(:),SFCLHX(:)
     &,SUBSHX(:),SNOPCX(:),ASWIN(:),ASWOUT(:)
     &,ASWTOA(:),ALWIN(:),ALWOUT(:),ALWTOA(:)
     &,TSHLTR(:),QSHLTR(:)
     &,TH10(:),Q10(:),U10(:),V10(:)
     &,TLMIN(:),TLMAX(:)
     &,SMC(:,:),CMC(:),STC(:,:),SH2O(:,:)
     &,VEGFRC(:),POTFLX(:),PSLP(:)
     &,EGRID2(:),SM(:),SICE(:)
     &,FACTR(:),PSFC(:)
     &,STATPR(:),STACPR(:),STAEVP(:)
     &,STAPOT(:),STASHX(:),STASUB(:),STAPCX(:)
     &,STASWI(:),STASWO(:),STALWI(:),STALWO(:)
     &,STALWT(:),STASWT(:),STASNM(:),STASRF(:)
     &,STABRF(:),STASNO(:),SNO(:),SR(:)
     &,ACPREC0(:),CUPREC0(:),SFCLHX0(:),POTFLX0(:)
     &,SFCSHX0(:),SUBSHX0(:),SNOPCX0(:),ASWIN0(:)
     &,ASWOUT0(:),ALWIN0(:),ALWOUT0(:),ALWTOA0(:)
     &,ASWTOA0(:),ACSNOW0(:),ACSNOM0(:),SSROFF0(:)
     &,BGROFF0(:)

C

      REAL, ALLOCATABLE:: T(:,:),Q(:,:),U(:,:),V(:,:),Q2(:,:)
     &,           CWM(:,:),TRAIN(:,:),TCUCN(:,:)
     &,           RSWTT(:,:),RLWTT(:,:),RTOP(:,:)
     &,           OMGA(:,:),t_hold(:,:)
     &,           PINT(:,:),UL(:),DPRES(:,:),DZ(:,:),CLDFRA(:,:)
     &,           DUM3D(:,:),DUM3D2(:,:),DUM3D3(:,:), DUM3D4(:,:)
C  

      REAL, ALLOCATABLE:: DHCNVC(:,:),DHRAIN(:,:),STADHC(:),STADHR(:),
     &                 TCUCN0(:,:),TRAIN0(:,:),VARTMP(:),VARTMPSTA(:)

      REAL,ALLOCATABLE:: DUM(:,:,:),DUMMY(:,:),DUMMY2(:,:)

        


      REAL :: DUMSC
      REAL, ALLOCATABLE:: PRODAT(:),FPACK(:)

      INTEGER, ALLOCATABLE:: LMH(:,:),IW(:,:),ICNT(:),IDSPL(:)


       REAL, ALLOCATABLE :: PMID(:,:),W(:,:), WH(:,:),
     &                      PDS(:)

	real, allocatable:: CROT(:),SROT(:),ak(:)
C------------------------------------------------------------------------
       INTEGER, ALLOCATABLE :: IDSTN(:),IHINDX(:),JHINDX(:)
       INTEGER, ALLOCATABLE :: IVINDX(:),JVINDX(:)
       INTEGER:: GDS(200),IDAT(3)
C
C------------------------------------------------------------------------
                             L O G I C A L
     & RUN,RESTRT,FRST
C------------------------------------------------------------------------
                             C H A R A C T E R
     & RSTFIL*90,RESTHR*4,LABEL*32,CISTAT*8
     &,FNAME*90,ENVAR*90,BLANK*4

	CHARACTER(LEN=8), ALLOCATABLE :: CIDSTN_SAVE(:),
     &                                   CIDSTN(:)

C	new stuff
      character(len=31) :: VarName,varin
	character(len=180) :: fileName,filedyn,filephys
	character(len=180) :: fileNamehold
      integer ::  DataHandle
      character(len=19):: startdate,datestr,datestrold
      character(len=2):: fhroldchar

	real:: rinc(5)
	integer:: IDATE(8),JDATE(8),im,jm,lm,K
        integer :: varid, ncid_dyn, ncid_phys

C------------------------------------------------------------------------
      DATA BLANK/'    '/
C------------------------------------------------------------------------
C***
C***  READ IN THE INFORMATION FILE ABOUT THE SOUNDINGS
C***

c	write(6,*) 'filename= ', filename
c	write(6,*) 'startedate= ', startdate

!	datestr=startdate

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)

      ALLOCATE(CIDSTN(NSTAT))
      ALLOCATE(STNLAT(NSTAT))
      ALLOCATE(STNLON(NSTAT))
      ALLOCATE(IDSTN(NSTAT))
      ALLOCATE(IHINDX(NSTAT))
      ALLOCATE(JHINDX(NSTAT))
      ALLOCATE(IVINDX(NSTAT))
      ALLOCATE(JVINDX(NSTAT))

	if (MYPE .eq. 0) then

      REWIND 19
      READ(19)NUMSTA,IDSTN,STNLAT,STNLON
     1,       IHINDX,JHINDX,IVINDX,JVINDX,CIDSTN
	write(0,*) 'read from unit 19'
!mp
	allocate(CIDSTN_SAVE(NUMSTA))
	DO N=1,NUMSTA
	CIDSTN_SAVE(N)=CIDSTN(N)
	ENDDO
!mp
        endif

       call mpi_bcast(NUMSTA,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

       call mpi_bcast(IDSTN,NSTAT,MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
       call mpi_bcast(STNLAT,NSTAT,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call mpi_bcast(STNLON,NSTAT,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call mpi_bcast(IHINDX,NSTAT,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call mpi_bcast(JHINDX,NSTAT,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call mpi_bcast(IVINDX,NSTAT,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call mpi_bcast(JVINDX,NSTAT,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
!	write(6,*) 'STNLAT(1), STNLON(1): ', STNLAT(1), STNLON(1)
!	write(6,*) 'IHINDX(1),JHINDX(1): ', IHINDX(1),JHINDX(1)
!	write(6,*) 'IVINDX(1),JVINDX(1): ', IVINDX(1),JVINDX(1)

	if (mype .eq. 0) then
      WRITE(6,20)NUMSTA
         endif
   20 FORMAT('INIT:  NUMBER OF PROFILE STATIONS ',I5)

	if (ITAG .eq. 0 .and. mype .eq. 0) then
      WRITE(6,30)(IDSTN(N),STNLAT(N)/DTR,STNLON(N)/DTR
     1,               IHINDX(N),JHINDX(N),IVINDX(N),JVINDX(N)
     2,               CIDSTN(N),N=1,NUMSTA)
	endif
   30 FORMAT(2X,I6,2F8.2,4I8,4X,A8)

c	if (ITAG .eq. 0) then
	  FRST=.TRUE.
c	else 
c	  FRST=.FALSE.
c	endif


!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

       if ( frst ) then
         frst = .false.
	write(6,*) 'filedyn early in PROF= ', trim(filedyn)
	write(6,*) 'filephys early in PROF= ', trim(filephys)

         call check( nf90_open(filedyn, NF90_NOWRITE, ncid_dyn) )
         call check( nf90_open(filephys, NF90_NOWRITE, ncid_phys) )

       endif

C Getting start time

        varname='time'
        call check( nf90_inq_varid(ncid_dyn,trim(varname),varid))
!        write(0,*) ncid_dyn,trim(varname),varid
        call check( nf90_get_var(ncid_dyn,varid,ihr))

        call check( nf_get_att_text(ncid_dyn,varid,'units',varin))
!        write(0,*) 'varin is: ', varin, ' and end'

        read(varin,101)idate(1),idate(2),idate(3),idate(4),idate(5)
 101    format(T13,i4,1x,i2,1x,i2,1x,i2,1x,i2)


        iyear=IDATE(1)
        imn=IDATE(2)
        iday=IDATE(3)
        ihrst=IDATE(4)

 15   format(i4,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr =',iyear,imn,iday,ihrst

      IFHR=ITAG

          call check( nf90_inq_dimid(ncid_dyn,'grid_xt',varid))
          call check( nf90_inquire_dimension(ncid_dyn,varid,len=im))

          call check( nf90_inq_dimid(ncid_dyn,'grid_yt',varid) )
          call check( nf90_inquire_dimension(ncid_dyn,varid,len=jm))

          call check(nf90_inq_dimid(ncid_dyn,'pfull',varid))
          call check(nf90_inquire_dimension(ncid_dyn,varid,len=lm))

          if (allocated(ak)) deallocate(ak); allocate(ak(LM))
!
!
! PHY file
!
!

!        float acond(time, grid_yt, grid_xt) ;
!                acond:long_name = "Aerodynamic conductance" ;
!        float albdo_ave(time, grid_yt, grid_xt) ;
!                albdo_ave:long_name = "surface albedo" ;
!        float alnsf(time, grid_yt, grid_xt) ;
!                alnsf:long_name = "mean nir albedo with strong cosz dependency" ;
!        float alnwf(time, grid_yt, grid_xt) ;
!                alnwf:long_name = "mean nir albedo with weak cosz dependency" ;
!        float alvsf(time, grid_yt, grid_xt) ;
!                alvsf:long_name = "mean vis albedo with strong cosz dependency" ;
!        float alvwf(time, grid_yt, grid_xt) ;
!               alvwf:long_name = "mean vis albedo with weak cosz dependency" ;
!        float c0(time, grid_yt, grid_xt) ;
!                c0:long_name = "nsst coefficient1 to calculate d(tz)/d(ts)" ;
!        float cd(time, grid_yt, grid_xt) ;
!                cd:long_name = "nsst coefficient2 to calculate d(tz)/d(ts)" ;
!        float cduvb_ave(time, grid_yt, grid_xt) ;
!                cduvb_ave:long_name = "Clear sky UV-B Downward Solar Flux" ;
!        float cnwat(time, grid_yt, grid_xt) ;
!                cnwat:long_name = "canopy water (cnwat in gfs data)" ;
!        float cpofp(time, grid_yt, grid_xt) ;
!                cpofp:long_name = "Percent frozen precipitation" ;
!        float cprat_ave(time, grid_yt, grid_xt) ;
!                cprat_ave:long_name = "averaged surface convective precipitation rate" ;
!        float cpratb_ave(time, grid_yt, grid_xt) ;
!                cpratb_ave:long_name = "averaged bucket surface convective precipitation rate" ;
!        float csdlf(time, grid_yt, grid_xt) ;
!                csdlf:long_name = "Clear Sky Downward Long Wave Flux" ;
!        float csdsf(time, grid_yt, grid_xt) ;
!                csdsf:long_name = "Clear Sky Downward Short Wave Flux" ;
!        float csulf(time, grid_yt, grid_xt) ;
!                csulf:long_name = "Clear Sky Upward Long Wave Flux" ;
!        float csulftoa(time, grid_yt, grid_xt) ;
!                csulftoa:long_name = "Clear Sky Upward Long Wave Flux at toa" ;
!        float csusf(time, grid_yt, grid_xt) ;
!                csusf:long_name = "Clear Sky Upward Short Wave Flux" ;
!        float csusftoa(time, grid_yt, grid_xt) ;
!                csusftoa:long_name = "Clear Sky Upward Short Wave Flux at toa" ;
!        float cwork_aveclm(time, grid_yt, grid_xt) ;
!                cwork_aveclm:long_name = "cloud work function (valid only with sas)" ;
!        float dconv(time, grid_yt, grid_xt) ;
!                dconv:long_name = "nsst thickness of free convection layer" ;
!        float dlwrf(time, grid_yt, grid_xt) ;
!                dlwrf:long_name = "instantaneous surface downward longwave flux" ;
!       float dlwrf_ave(time, grid_yt, grid_xt) ;
!                dlwrf_ave:long_name = "surface downward longwave flux" ;
!       float dswrf(time, grid_yt, grid_xt) ;
!                dswrf:long_name = "instantaneous surface downward shortwave flux" ;
!        float dswrf_ave(time, grid_yt, grid_xt) ;
!                dswrf_ave:long_name = "averaged surface downward shortwave flux" ;
!        float dswrf_avetoa(time, grid_yt, grid_xt) ;
!                dswrf_avetoa:long_name = "top of atmos downward shortwave flux" ;
!        float dtcool(time, grid_yt, grid_xt) ;
!                dtcool:long_name = "nsst sub-layer cooling amount" ;
!        float duvb_ave(time, grid_yt, grid_xt) ;
!                duvb_ave:long_name = "UV-B Downward Solar Flux" ;
!        float evbs_ave(time, grid_yt, grid_xt) ;
!                evbs_ave:long_name = "Direct Evaporation from Bare Soil - GFS lsm" ;
!        float evcw_ave(time, grid_yt, grid_xt) ;
!                evcw_ave:long_name = "Canopy water evaporation - GFS lsm" ;
!        float f10m(time, grid_yt, grid_xt) ;
!                f10m:long_name = "10-meter wind speed divided by lowest model wind speed" ;
!        float facsf(time, grid_yt, grid_xt) ;
!                facsf:long_name = "fractional coverage with strong cosz dependency" ;
!        float facwf(time, grid_yt, grid_xt) ;
!                facwf:long_name = "fractional coverage with weak cosz dependency" ;
!        float ffhh(time, grid_yt, grid_xt) ;
!                ffhh:long_name = "fh parameter from PBL scheme" ;
!        float ffmm(time, grid_yt, grid_xt) ;
!                ffmm:long_name = "fm parameter from PBL scheme" ;
!        float fldcp(time, grid_yt, grid_xt) ;
!                fldcp:long_name = "Field Capacity (volumetric)" ;
!       float fricv(time, grid_yt, grid_xt) ;
!                fricv:long_name = "uustar surface frictional wind" ;
!        float gflux(time, grid_yt, grid_xt) ;
!                gflux:long_name = "instantaneous surface ground heat flux" ;
!        float gflux_ave(time, grid_yt, grid_xt) ;
!                gflux_ave:long_name = "surface ground heat flux" ;
!        float hgt_hyblev1(time, grid_yt, grid_xt) ;
!                hgt_hyblev1:long_name = "layer 1 height" ;
!        float hpbl(time, grid_yt, grid_xt) ;
!                hpbl:long_name = "surface planetary boundary layer height" ;
!        float icec(time, grid_yt, grid_xt) ;
!                icec:long_name = "surface ice concentration (ice=1; no ice=0)" ;
!        float icetk(time, grid_yt, grid_xt) ;
!                icetk:long_name = "sea ice thickness (icetk in gfs_data)" ;
!        float land(time, grid_yt, grid_xt) ;
!                land:long_name = "sea-land-ice mask (0-sea, 1-land, 2-ice)" ;
!        float lhtfl(time, grid_yt, grid_xt) ;
!                lhtfl:long_name = "instantaneous surface latent heat net flux" ;
!        float lhtfl_ave(time, grid_yt, grid_xt) ;
!                lhtfl_ave:long_name = "surface latent heat flux" ;
!        float nbdsf_ave(time, grid_yt, grid_xt) ;
!                nbdsf_ave:long_name = "Near IR Beam Downward Solar Flux" ;
!        float nddsf_ave(time, grid_yt, grid_xt) ;
!                nddsf_ave:long_name = "Near IR Diffuse Downward Solar Flux" ;
!        float orog(time, grid_yt, grid_xt) ;
!                orog:long_name = "surface geopotential height" ;
!        float pevpr(time, grid_yt, grid_xt) ;
!                pevpr:long_name = "instantaneous surface potential evaporation" ;
!        float pevpr_ave(time, grid_yt, grid_xt) ;
!                pevpr_ave:long_name = "averaged potential evaporation rate" ;
!        float prate_ave(time, grid_yt, grid_xt) ;
!                prate_ave:long_name = "surface precipitation rate" ;
!        float prateb_ave(time, grid_yt, grid_xt) ;
!                prateb_ave:long_name = "bucket surface precipitation rate" ;
!        float pres_avehcb(time, grid_yt, grid_xt) ;
!                pres_avehcb:long_name = "pressure high cloud bottom level" ;
!        float pres_avehct(time, grid_yt, grid_xt) ;
!                pres_avehct:long_name = "pressure high cloud top level" ;
!        float pres_avelcb(time, grid_yt, grid_xt) ;
!                pres_avelcb:long_name = "pressure low cloud bottom level" ;
!        float pres_avelct(time, grid_yt, grid_xt) ;
!                pres_avelct:long_name = "pressure low cloud top level" ;
!        float pres_avemcb(time, grid_yt, grid_xt) ;
!                pres_avemcb:long_name = "pressure middle cloud bottom level" ;
!        float pres_avemct(time, grid_yt, grid_xt) ;
!                pres_avemct:long_name = "pressure middle cloud top level" ;
!        float prescnvclb(time, grid_yt, grid_xt) ;
!                prescnvclb:long_name = "pressure at convective cloud bottom level" ;
!        float prescnvclt(time, grid_yt, grid_xt) ;
!                prescnvclt:long_name = "pressure at convective cloud top level" ;
!        float pressfc(time, grid_yt, grid_xt) ;
!                pressfc:long_name = "surface pressure" ;
!        float pwatclm(time, grid_yt, grid_xt) ;
!                pwatclm:long_name = "atmos column precipitable water" ;
!        float qrain(time, grid_yt, grid_xt) ;
!                qrain:long_name = "nsst sensible heat flux due to rainfall" ;
!        float refdmax(time, grid_yt, grid_xt) ;
!                refdmax:long_name = "max hourly 1-km agl reflectivity" ;
!        float refdmax263k(time, grid_yt, grid_xt) ;
!                refdmax263k:long_name = "max hourly -10C reflectivity" ;
!        float refl_10cm(time, pfull, grid_yt, grid_xt) ;
!                refl_10cm:long_name = "Radar reflectivity" ;
!        float rh02max(time, grid_yt, grid_xt) ;
!                rh02max:long_name = "max hourly 2m RH" ;
!        float rh02min(time, grid_yt, grid_xt) ;
!                rh02min:long_name = "min hourly 2m RH" ;
!        float sbsno_ave(time, grid_yt, grid_xt) ;
!                sbsno_ave:long_name = "Sublimation (evaporation from snow) - GFS lsm" ;
!        float sfcr(time, grid_yt, grid_xt) ;
!                sfcr:long_name = "surface roughness" ;
!        float sfexc(time, grid_yt, grid_xt) ;
!                sfexc:long_name = "Exchange Coefficient" ;
!        float shdmax(time, grid_yt, grid_xt) ;
!                shdmax:long_name = "maximum fractional coverage of green vegetation" ;
!        float shdmin(time, grid_yt, grid_xt) ;
!                shdmin:long_name = "minimum fractional coverage of green vegetation" ;
!        float shtfl(time, grid_yt, grid_xt) ;
!                shtfl:long_name = "instantaneous surface sensible heat net flux" ;
!        float shtfl_ave(time, grid_yt, grid_xt) ;
!                shtfl_ave:long_name = "surface sensible heat flux" ;
!        float sltyp(time, grid_yt, grid_xt) ;
!                sltyp:long_name = "surface slope type" ;
!        float snoalb(time, grid_yt, grid_xt) ;
!                snoalb:long_name = "maximum snow albedo in fraction" ;
!        float snod(time, grid_yt, grid_xt) ;
!                snod:long_name = "surface snow depth" ;
!        float snohf(time, grid_yt, grid_xt) ;
!                snohf:long_name = "Snow Phase Change Heat Flux - GFS lsm" ;
!        float snowc_ave(time, grid_yt, grid_xt) ;
!                snowc_ave:long_name = "snow cover - GFS lsm" ;
!        float soill1(time, grid_yt, grid_xt) ;
!                soill1:long_name = "liquid soil mositure at layer-1" ;
!        float soill2(time, grid_yt, grid_xt) ;
!                soill2:long_name = "liquid soil mositure at layer-2" ;
!        float soill3(time, grid_yt, grid_xt) ;
!                soill3:long_name = "liquid soil mositure at layer-3" ;
!        float soill4(time, grid_yt, grid_xt) ;
!                soill4:long_name = "liquid soil mositure at layer-4" ;
!        float soilm(time, grid_yt, grid_xt) ;
!                soilm:long_name = "total column soil oisture content" ;
!        float soilt1(time, grid_yt, grid_xt) ;
!                soilt1:long_name = "soil temperature 0-10cm" ;
!        float soilt2(time, grid_yt, grid_xt) ;
!                soilt2:long_name = "soil temperature 10-40cm" ;
!        float soilt3(time, grid_yt, grid_xt) ;
!                soilt3:long_name = "soil temperature 40-100cm" ;
!        float soilt4(time, grid_yt, grid_xt) ;
!                soilt4:long_name = "soil temperature 100-200cm" ;
!        float soilw1(time, grid_yt, grid_xt) ;
!                soilw1:long_name = "volumetric soil moisture 0-10cm" ;
!        float soilw2(time, grid_yt, grid_xt) ;
!                soilw2:long_name = "volumetric soil moisture 10-40cm" ;
!        float soilw3(time, grid_yt, grid_xt) ;
!                soilw3:long_name = "volumetric soil moisture 40-100cm" ;
!        float soilw4(time, grid_yt, grid_xt) ;
!                soilw4:long_name = "volumetric soil moisture 100-200cm" ;
!        float sotyp(time, grid_yt, grid_xt) ;
!                sotyp:long_name = "soil type in integer 1-9" ;
!        float spd10max(time, grid_yt, grid_xt) ;
!                spd10max:long_name = "hourly maximum wind speed" ;
!        float spfh2m(time, grid_yt, grid_xt) ;
!                spfh2m:long_name = "2m specific humidity" ;
!        float spfh_hyblev1(time, grid_yt, grid_xt) ;
!                spfh_hyblev1:long_name = "layer 1 specific humidity" ;
!        float spfhmax_max2m(time, grid_yt, grid_xt) ;
!                spfhmax_max2m:long_name = "maximum specific humidity" ;
!        float spfhmin_min2m(time, grid_yt, grid_xt) ;
!                spfhmin_min2m:long_name = "minimum specific humidity" ;
!        float ssrun_acc(time, grid_yt, grid_xt) ;
!                ssrun_acc:long_name = "surface storm water runoff - GFS lsm" ;
!        float sunsd_acc(time, grid_yt, grid_xt) ;
!                sunsd_acc:long_name = "Sunshine Duration" ;
!        float t02max(time, grid_yt, grid_xt) ;
!                t02max:long_name = "max hourly 2m Temperature" ;
!        float t02min(time, grid_yt, grid_xt) ;
!                t02min:long_name = "min hourly 2m Temperature" ;
!        float tcdc_avebndcl(time, grid_yt, grid_xt) ;
!                tcdc_avebndcl:long_name = "boundary layer cloud layer total cloud cover" ;
!        float tcdc_aveclm(time, grid_yt, grid_xt) ;
!                tcdc_aveclm:long_name = "atmos column total cloud cover" ;
!        float tcdc_avehcl(time, grid_yt, grid_xt) ;
!                tcdc_avehcl:long_name = "high cloud level total cloud cover" ;
!        float tcdc_avelcl(time, grid_yt, grid_xt) ;
!                tcdc_avelcl:long_name = "low cloud level total cloud cover" ;
!        float tcdc_avemcl(time, grid_yt, grid_xt) ;
!                tcdc_avemcl:long_name = "mid cloud level total cloud cover" ;
!        float tcdccnvcl(time, grid_yt, grid_xt) ;
!                tcdccnvcl:long_name = "convective cloud layer total cloud cover" ;
!        float tg3(time, grid_yt, grid_xt) ;
!                tg3:long_name = "deep soil temperature" ;
!        float tisfc(time, grid_yt, grid_xt) ;
!                tisfc:long_name = "surface temperature over ice fraction" ;
!        float tmax_max2m(time, grid_yt, grid_xt) ;
!                tmax_max2m:long_name = "max temperature at 2m height" ;
!        float tmin_min2m(time, grid_yt, grid_xt) ;
!                tmin_min2m:long_name = "min temperature at 2m height" ;
!        float tmp2m(time, grid_yt, grid_xt) ;
!                tmp2m:long_name = "2m temperature" ;
!        float tmp_avehct(time, grid_yt, grid_xt) ;
!                tmp_avehct:long_name = "temperature high cloud top level" ;
!        float tmp_avelct(time, grid_yt, grid_xt) ;
!                tmp_avelct:long_name = "temperature low cloud top level" ;
!        float tmp_avemct(time, grid_yt, grid_xt) ;
!                tmp_avemct:long_name = "temperature middle cloud top level" ;
!        float tmp_hyblev1(time, grid_yt, grid_xt) ;
!                tmp_hyblev1:long_name = "layer 1 temperature" ;
!        float tmpsfc(time, grid_yt, grid_xt) ;
!                tmpsfc:long_name = "surface temperature" ;
!        float tprcp(time, grid_yt, grid_xt) ;
!                tprcp:long_name = "total precipitation" ;
!        float trans_ave(time, grid_yt, grid_xt) ;
!                trans_ave:long_name = "transpiration - GFS lsm" ;
!        float tref(time, grid_yt, grid_xt) ;
!                tref:long_name = "nsst reference or foundation temperature" ;
!        float u-gwd_ave(time, grid_yt, grid_xt) ;
!                u-gwd_ave:long_name = "surface zonal gravity wave stress" ;
!        float u10max(time, grid_yt, grid_xt) ;
!                u10max:long_name = "hourly maximum (magnitude) u-wind" ;
!        float uflx_ave(time, grid_yt, grid_xt) ;
!                uflx_ave:long_name = "surface zonal momentum flux" ;
!        float ugrd10m(time, grid_yt, grid_xt) ;
!                ugrd10m:long_name = "10 meter u wind" ;
!        float ugrd_hyblev1(time, grid_yt, grid_xt) ;
!                ugrd_hyblev1:long_name = "layer 1 zonal wind" ;
!        float ulwrf(time, grid_yt, grid_xt) ;
!                ulwrf:long_name = "instantaneous surface upward longwave flux" ;
!        float ulwrf_ave(time, grid_yt, grid_xt) ;
!                ulwrf_ave:long_name = "surface upward longwave flux" ;
!        float ulwrf_avetoa(time, grid_yt, grid_xt) ;
!                ulwrf_avetoa:long_name = "top of atmos upward longwave flux" ;
!        float uswrf(time, grid_yt, grid_xt) ;
!                uswrf:long_name = "instantaneous surface upward shortwave flux" ;
!        float uswrf_ave(time, grid_yt, grid_xt) ;
!                uswrf_ave:long_name = "averaged surface upward shortwave flux" ;
!        float uswrf_avetoa(time, grid_yt, grid_xt) ;
!                uswrf_avetoa:long_name = "top of atmos upward shortwave flux" ;
!        float v-gwd_ave(time, grid_yt, grid_xt) ;
!                v-gwd_ave:long_name = "surface meridional gravity wave stress" ;
!        float v10max(time, grid_yt, grid_xt) ;
!                v10max:long_name = "hourly maximum (magnitude) v-wind" ;
!        float vbdsf_ave(time, grid_yt, grid_xt) ;
!                vbdsf_ave:long_name = "Visible Beam Downward Solar Flux" ;
!        float vddsf_ave(time, grid_yt, grid_xt) ;
!                vddsf_ave:long_name = "Visible Diffuse Downward Solar Flux" ;
!        float veg(time, grid_yt, grid_xt) ;
!                veg:long_name = "vegetation fraction" ;
!        float vflx_ave(time, grid_yt, grid_xt) ;
!                vflx_ave:long_name = "surface meridional momentum flux" ;
!        float vgrd10m(time, grid_yt, grid_xt) ;
!                vgrd10m:long_name = "10 meter v wind" ;
!        float vgrd_hyblev1(time, grid_yt, grid_xt) ;
!                vgrd_hyblev1:long_name = "layer 1 meridional wind" ;
!        float vtype(time, grid_yt, grid_xt) ;
!                vtype:long_name = "vegetation type in integer 1-13" ;
!        float w0(time, grid_yt, grid_xt) ;
!                w0:long_name = "nsst coefficient3 to calculate d(tz)/d(ts)" ;
!        float watr_acc(time, grid_yt, grid_xt) ;
!                watr_acc:long_name = "total water runoff" ;
!        float wd(time, grid_yt, grid_xt) ;
!                wd:long_name = "nsst coefficient4 to calculate d(tz)/d(ts)" ;
!        float weasd(time, grid_yt, grid_xt) ;
!                weasd:long_name = "surface snow water equivalent" ;
!        float wilt(time, grid_yt, grid_xt) ;
!                wilt:long_name = "wiltimg point (volumetric)" ;
!        float xs(time, grid_yt, grid_xt) ;
!                xs:long_name = "nsst salinity content in diurnal thermocline layer" ;
!        float xt(time, grid_yt, grid_xt) ;
!                xt:long_name = "nsst heat content in diurnal thermocline layer" ;
!        float xtts(time, grid_yt, grid_xt) ;
!                xtts:long_name = "nsst d(xt)/d(ts)" ;
!        float xu(time, grid_yt, grid_xt) ;
!                xu:long_name = "nsst u-current content in diurnal thermocline layer" ;
!        float xv(time, grid_yt, grid_xt) ;
!                xv:long_name = "nsst v-current content in diurnal thermocline layer" ;
!        float xz(time, grid_yt, grid_xt) ;
!                xz:long_name = "nsst diurnal thermocline layer thickness" ;
!        float xzts(time, grid_yt, grid_xt) ;
!                xzts:long_name = "nsst d(xt)/d(ts)" ;
!        float zc(time, grid_yt, grid_xt) ;
!                zc:long_name = "nsst sub-layer cooling thickness" ;
!


!
! DYN
!

! 3D
!
!        float cld_amt(time, pfull, grid_yt, grid_xt) ;
!                cld_amt:long_name = "cloud amount" ;
!        float clwmr(time, pfull, grid_yt, grid_xt) ;
!                clwmr:long_name = "cloud water mixing ratio" ;
!        float delz(time, pfull, grid_yt, grid_xt) ;
!                delz:long_name = "height thickness" ;
!        float dpres(time, pfull, grid_yt, grid_xt) ;
!                dpres:long_name = "pressure thickness" ;
!        float dzdt(time, pfull, grid_yt, grid_xt) ;
!                dzdt:long_name = "vertical wind" ;
!        float grle(time, pfull, grid_yt, grid_xt) ;
!                grle:long_name = "graupel mixing ratio" ;
!        float icmr(time, pfull, grid_yt, grid_xt) ;
!                icmr:long_name = "cloud ice mixing ratio" ;
!        float o3mr(time, pfull, grid_yt, grid_xt) ;
!                o3mr:long_name = "ozone mixing ratio" ;
!        float rwmr(time, pfull, grid_yt, grid_xt) ;
!                rwmr:long_name = "rain mixing ratio" ;
!        float snmr(time, pfull, grid_yt, grid_xt) ;
!                snmr:long_name = "snow mixing ratio" ;
!        float spfh(time, pfull, grid_yt, grid_xt) ;
!                spfh:long_name = "specific humidity" ;
!        float tmp(time, pfull, grid_yt, grid_xt) ;
!                tmp:long_name = "temperature" ;
!        float ugrd(time, pfull, grid_yt, grid_xt) ;
!                ugrd:long_name = "zonal wind" ;
!        float vgrd(time, pfull, grid_yt, grid_xt) ;
!                vgrd:long_name = "meridional wind" ;

! 2D 
!
!        float hgtsfc(time, grid_yt, grid_xt) ;
!                hgtsfc:long_name = "surface geopotential height" ;
!        float pressfc(time, grid_yt, grid_xt) ;
!                pressfc:long_name = "surface pressure" ;
!        float dnvvelmax(time, grid_yt, grid_xt) ;
!                dnvvelmax:long_name = "Max hourly downdraft velocity" ;
!        float srh01(time, grid_yt, grid_xt) ;
!                srh01:long_name = "0-1km srh" ;
!        float srh03(time, grid_yt, grid_xt) ;
!                srh03:long_name = "0-3km srh" ;
!        float uhmax03(time, grid_yt, grid_xt) ;
!                uhmax03:long_name = "Max hourly 0-3km updraft helicity" ;
!        float uhmax25(time, grid_yt, grid_xt) ;
!                uhmax25:long_name = "Max hourly 2-5km updraft helicity" ;
!        float uhmin03(time, grid_yt, grid_xt) ;
!                uhmin03:long_name = "Max hourly 0-3km updraft helicity" ;
!        float uhmin25(time, grid_yt, grid_xt) ;
!                uhmin25:long_name = "Max hourly 2-5km updraft helicity" ;
!        float upvvelmax(time, grid_yt, grid_xt) ;
!                upvvelmax:long_name = "Max hourly updraft velocity" ;
!        float ustm(time, grid_yt, grid_xt) ;
!                ustm:long_name = "u comp of storm motion" ;
!        float vstm(time, grid_yt, grid_xt) ;
!                vstm:long_name = "v comp of storm motion" ;
!        float maxvort01(time, grid_yt, grid_xt) ;
!                maxvort01:long_name = "Max hourly 0-1km vert. vorticity" ;
!        float maxvort02(time, grid_yt, grid_xt) ;
!                maxvort02:long_name = "Max hourly 0-2km vert. vorticity" ;
!        float maxvorthy1(time, grid_yt, grid_xt) ;
!                maxvorthy1:long_name = "Max hourly lev 1 vert vort." ;
!
!
!
!


        write(6,*) 'allocate with IM, JM, LM: ', IM, JM, LM


      allocate(ICNT(0:NPES-1),IDSPL(0:NPES-1))

       do I=0,npes-1
      call para_range(1,NUMSTA,npes,I,
     &  NSTART,NEND)
	ICNT(I)=(NEND-NSTART)+1
        IDSPL(I)= NSTART-1 ! ?
	if (mype .eq. 0) then
	write(0,*) 'mype, ICNT, IDSPL: ', I,ICNT(I),IDSPL(I)
        endif
       enddo

      call para_range(1,NUMSTA,npes,mype,
     &  NSTART,NEND)

      write(0,*) 'mype, NSTART, NEND: ', mype, NSTART, NEND,ICNT(MYPE)


!!!!!
       ALLOCATE(RES(NUMSTA),FIS(NUMSTA),THS(NUMSTA),HBOT(NUMSTA))
       ALLOCATE(CFRACL(NUMSTA),CFRACM(NUMSTA),CFRACH(NUMSTA))
       ALLOCATE(SOILTB(NUMSTA),SFCEXC(NUMSTA))
       ALLOCATE(SMSTAV(NUMSTA),SMSTOT(NUMSTA))
       ALLOCATE(Z0(NUMSTA),CZEN(NUMSTA),CZMEAN(NUMSTA))
       ALLOCATE(ACPREC(NUMSTA))
       ALLOCATE(CUPREC(NUMSTA),ACSNOW(NUMSTA),ACSNOM(NUMSTA))
       ALLOCATE(SSROFF(NUMSTA),BGROFF(NUMSTA),SFCSHX(NUMSTA))
       ALLOCATE(SFCLHX(NUMSTA),SUBSHX(NUMSTA),SNOPCX(NUMSTA))
       ALLOCATE(ASWIN(NUMSTA),ASWOUT(NUMSTA),ASWTOA(NUMSTA))
       ALLOCATE(ALWIN(NUMSTA),ALWOUT(NUMSTA),ALWTOA(NUMSTA))
       ALLOCATE(TSHLTR(NUMSTA),QSHLTR(NUMSTA))
       ALLOCATE(TH10(NUMSTA),Q10(NUMSTA),U10(NUMSTA),V10(NUMSTA))
       ALLOCATE(TLMIN(NUMSTA),TLMAX(NUMSTA),SMC(NUMSTA,NSOIL))
       ALLOCATE(CMC(NUMSTA),STC(NUMSTA,NSOIL),SH2O(NUMSTA,NSOIL))
       ALLOCATE(VEGFRC(NUMSTA),POTFLX(NUMSTA),PSLP(NUMSTA))
       ALLOCATE(EGRID2(NUMSTA),SM(NUMSTA),SICE(NUMSTA))
       ALLOCATE(FACTR(NUMSTA))
       ALLOCATE(STATPR(NUMSTA),STACPR(NUMSTA),STAEVP(NUMSTA))
       ALLOCATE(STAPOT(NUMSTA),STASHX(NUMSTA),STASUB(NUMSTA))
       ALLOCATE(STAPCX(NUMSTA),STASWI(NUMSTA),STASWO(NUMSTA))
       ALLOCATE(STALWI(NUMSTA),STALWO(NUMSTA),STALWT(NUMSTA))
       ALLOCATE(STASWT(NUMSTA),STASNM(NUMSTA),STASRF(NUMSTA))
       ALLOCATE(STABRF(NUMSTA),STASNO(NUMSTA),ACPREC0(NUMSTA))
       ALLOCATE(CUPREC0(NUMSTA),SFCLHX0(NUMSTA),POTFLX0(NUMSTA))
       ALLOCATE(SFCSHX0(NUMSTA),SUBSHX0(NUMSTA),SNOPCX0(NUMSTA))
       ALLOCATE(ASWIN0(NUMSTA),ASWOUT0(NUMSTA),ALWIN0(NUMSTA))
       ALLOCATE(ALWOUT0(NUMSTA),ALWTOA0(NUMSTA),ASWTOA0(NUMSTA))
       ALLOCATE(ACSNOW0(NUMSTA),ACSNOM0(NUMSTA),SSROFF0(NUMSTA))
       ALLOCATE(BGROFF0(NUMSTA),PSFC(NUMSTA),SNO(NUMSTA),SR(NUMSTA))


	ALLOCATE(T(NUMSTA,LM))
	ALLOCATE(Q(NUMSTA,LM))
	ALLOCATE(U(NUMSTA,LM))
	ALLOCATE(V(NUMSTA,LM))
	ALLOCATE(Q2(NUMSTA,LM))
	ALLOCATE(CWM(NUMSTA,LM))
	ALLOCATE(TRAIN(NUMSTA,LM))
	ALLOCATE(TCUCN(NUMSTA,LM))
	ALLOCATE(RSWTT(NUMSTA,LM))
	ALLOCATE(RLWTT(NUMSTA,LM))
	ALLOCATE(CLDFRA(NUMSTA,LM))
	ALLOCATE(RTOP(NUMSTA,LM))
	ALLOCATE(OMGA(NUMSTA,LM))
	ALLOCATE(t_hold(NUMSTA,LM))
	ALLOCATE(PINT(NUMSTA,LM+1))
        ALLOCATE(W(NUMSTA,LM))
        ALLOCATE(WH(NUMSTA,LM))
        ALLOCATE(IW(NUMSTA,LM))
        ALLOCATE(DPRES(NUMSTA,LM))
        ALLOCATE(DZ(NUMSTA,LM))

        ALLOCATE(STADHC(LM))
        ALLOCATE(STADHR(LM))
        ALLOCATE(VARTMP(LM))
        ALLOCATE(VARTMPSTA(NUMSTA))
        ALLOCATE(DHRAIN(LM,NUMSTA))
        ALLOCATE(DHCNVC(LM,NUMSTA))
        ALLOCATE(TCUCN0(LM,NUMSTA))
        ALLOCATE(TRAIN0(LM,NUMSTA))


! former parameter statements
        NWORDM=(LCL1ML+1)*LM+2*LCL1SL
        LRECPR=4*(8+9+LCL1ML*LM+LCL1SL)
! former parameter statements

        if (allocated(FPACK)) deallocate(FPACK); allocate(FPACK(NWORDM))
        if (allocated(PRODAT)) deallocate(PRODAT);
     &                  allocate(PRODAT(NWORDM))
        if (ALLOCATED(DUM)) deallocate(DUM);
     &                          allocate(DUM(IM,JM,4))
        if (ALLOCATED(DUMMY)) deallocate(DUMMY);
     &                          allocate(DUMMY(IM,JM))
        if (ALLOCATED(DUMMY2)) deallocate(DUMMY2);
     &                          allocate(DUMMY2(IM,JM))

        if (ALLOCATED(DUM3D)) deallocate(DUM3D);
     &                          allocate(DUM3D(NUMSTA,LM))
        if (ALLOCATED(DUM3D2)) deallocate(DUM3D2);
     &                          allocate(DUM3D2(NUMSTA,LM))
        if (ALLOCATED(DUM3D3)) deallocate(DUM3D3);
     &                          allocate(DUM3D3(NUMSTA,LM))
        if (ALLOCATED(DUM3D4)) deallocate(DUM3D4);
     &                          allocate(DUM3D4(NUMSTA,LM))

        if (ALLOCATED(LMH)) deallocate(LMH);
     &                          allocate(LMH(IM,JM))


!!!!!


!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"lon1", wbd))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"lat1", sbd))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"dlon", dxval))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"dlat", dyval))
!       call check(nf90_get_att(ncid_dyn,NF90_GLOBAL,"cen_lat",cenlat))
!       TPH0D=cenlat
!       cenlat=nint(1000.0*cenlat)
!       call check(nf90_get_att(ncid_dyn,NF90_GLOBAL,"cen_lon",cenlon))
!       TLM0D=cenlon
!       cenlon=nint(1000.0*cenlon)
!
!
!
!        write(6,*) 'dxval= ', dxval
!        write(6,*) 'dyval= ', dyval
!        write(6,*) 'cenlat= ', cenlat
!        write(6,*) 'cenlon= ', cenlon

!        truelat1=0.
!        truelat2=0.
!        write(6,*) 'truelat1= ', truelat1
!        write(6,*) 'truelat2= ', truelat2
        maptype=5
        write(6,*) 'maptype is ', maptype
!need to get DT
       call check (nf90_get_att(ncid_phys, NF90_GLOBAL,"dtp", DT))
        print*,'DT= ',DT

! get 3-D variables
      print*,'im,jm,lm= ',im,jm,lm
c
	write(6,*) 'DateStr: ', DateStr

	! start reading 3d netcdf output
 
       varname='ugrd'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)

       U(N,:)=VARTMP(:)

	if (N .eq. 20) then
           write(0,*) ' U(20,20): ', U(N,20)
        endif

       ENDDO
       varname='vgrd'

       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)

       V(N,:)=VARTMP(:)

	if (N .eq. 20) then
           write(0,*) 'V(20,20): ', V(N,20)
        endif

       ENDDO

!       enddo


! W defined over LM, not LM+1??

       varname='dzdt'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       W(N,:)=VARTMP(:)
       WH(N,:)=VARTMP(:)
	if (N .eq. 20) then
           write(0,*) 'W(20,20): ', W(N,20)
        endif

       ENDDO

!      DO L = 1,LM
!      DO N=NSTART,NEND
!	I=IHINDX(N)	
!	J=JHINDX(N)
!            WH(N,L) = (W(N,L)+W(N,L+1))*0.5
!      END DO
!      END DO


       varname='tmp'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       DUM3D(N,:)=VARTMP(:)
	if (N .eq. 20) then
           write(0,*) 'T(20,20): ', DUM3D(N,20)
        endif
       ENDDO
       

       varname='spfh'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       DUM3D2(N,:)=VARTMP(:)
	if (N .eq. 10) then
           do L=1,LM
           write(0,*) 'L,Q(10,L): ', L, DUM3D2(N,L)
           enddo
        endif
       ENDDO


       varname='dpres'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       DUM3D3(N,:)=VARTMP(:)
	if (N .eq. 20) then
           write(0,*) 'DP(20,20): ', DUM3D3(N,20)
        endif
       ENDDO


       varname='delz'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       DUM3D4(N,:)=VARTMP(:)
	if (N .eq. 20) then
           write(0,*) 'DZ(20,20): ', DUM3D4(N,20)
        endif
       ENDDO


        varname='pressfc'
       DO N=NSTART,NEND
       call read_netcdf_2d(ncid_dyn,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
       PSFC(N)=DUMSC
       ENDDO

       write(0,*) 'past PSFC def with min/max: ', 
     &           minval(psfc),maxval(psfc)

	DO L=1,LM
	DO N=NSTART,NEND
	 Q(N,L)=DUM3D2(N,L)  ! need to convert in any way?
         T_hold( N , L ) = dum3d ( N, L ) 
         DPRES(N,L) = DUM3D3(N,L)
         DZ(N,L) = abs(DUM3D4(N,L))
	ENDDO
	ENDDO

	if (mype .eq. 0) then
       write(0,*) 'past big def block'
	write(0,*) 'Q(1,1),T_hold(1,1),dpres(1,1),dz(1,1): ', 
     &              Q(1,1),T_hold(1,1),dpres(1,1),dz(1,1)
	endif


! think interface vs. midlayer pressure - this below might not be correct.

	if (allocated(PMID)) deallocate(PMID)
	allocate(PMID(NUMSTA,LM))

!	DO N=NSTART,NEND
!         I=IHINDX(N)
!         J=JHINDX(N)

! test going to post style PMID calc.  
!was         PMID(N,LM)=PSFC(N)-0.5*DPRES(N,LM)
!	do L=LM-1,1,-1
!         PMID(N,L)=PMID(N,L+1)-(DPRES(N,L))
!	if (PMID(N,L) .eq. -9999.) then
!	  write(0,*) 'bad at N, L: ', N, L
!          write(0,*) 'PMID(N,L+1),DPRES(N,L): ', PMID(N,L+1),DPRES(N,L)
!was	endif

!	if (N .eq. 1) then
!	 write(0,*) 'L, DPRES(N,L), PMID(N,L): ', L, DPRES(N,L), 
!     &                                          PMID(N,L)
!        endif

!was        enddo

!        ENDDO



        if (allocated(PDS)) deallocate(PDS)
        allocate(PDS(NUMSTA))

! Did I mistranslate something wrong here?  the PMID/DPRES business seems weird.

      do L=1,LM
	DO N=NSTART,NEND
           if(DPRES(N,L)/=spval .and. T_hold(N,L)/=spval .and. 
     &     Q(N,L)/=spval .and. DZ(N,L)/=spval)then
           PMID(N,L)=con_rd*DPRES(N,L)* 
     &          T_hold(N,L)*(Q(N,L)*con_fvirt+1.0)/G/DZ(N,L)

!	if (N .eq. NSTART .or. N .eq. NEND ) then
!        write(0,*) 'N,L,DP,T,Q,DZ,PMID: ',N,L,DPRES(N,L),T_hold(N,L),
!     &                                  Q(N,L),DZ(N,L),PMID(N,L)
!	endif

           else
            PMID(N,L)=spval

           end if
!! dong add missing value
           if (w(N,L) < spval) then
            omga(N,L)=(-1.)*w(N,L)*DPRES(N,L)/DZ(N,L)
           else
            omga(N,L) = spval
           end if
        ENDDO 
       ENDDO


       varname='clwmr'
       DO N=NSTART,NEND
       call read_netcdf_3d(ncid_dyn,ifhr
     & ,varname,IHINDX(N),JHINDX(N),LM,VARTMP)
       CWM(N,:)=VARTMP(:)
       ENDDO


        DO L = 1, LM
        DO N=NSTART,NEND
          Q2(N,L)=0.
        ENDDO
        ENDDO


      VarName='soilt1'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        STC(N,1) = DUMSC
        ENDDO


      VarName='soilt2'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        STC(N,2) = DUMSC
        ENDDO

      VarName='soilt3'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        STC(N,3) = DUMSC
        ENDDO


      VarName='soilt4'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        STC(N,4) = DUMSC
        ENDDO


      VarName='soilw1'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SMC(N,1) = DUMSC
        ENDDO

      VarName='soilw2'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SMC(N,2) = DUMSC
        ENDDO

      VarName='soilw3'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SMC(N,3) = DUMSC
        ENDDO

      VarName='soilw4'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SMC(N,4) = DUMSC
        ENDDO


      VarName='soill1'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SH2O(N,1) = DUMSC
        ENDDO

      VarName='soill2'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SH2O(N,2) = DUMSC
        ENDDO

      VarName='soill3'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SH2O(N,3) = DUMSC
        ENDDO

      VarName='soill4'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        SH2O(N,4) = DUMSC
        ENDDO


      VarName='spfh2m'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
! should QSHLTER be specific humidity?
       QSHLTR(N)=DUMSC
       ENDDO

      VarName='tmp2m'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        TSHLTR(N)=DUMSC
        ENDDO

      VarName='ugrd10m'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        U10(N)=DUMSC
         ENDDO

      VarName='vgrd10m'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        V10(N)=DUMSC
         ENDDO

      DO N=NSTART,NEND
        TH10(N)=-9999.
        Q10(N)=-9999.
       END DO

      VarName='cnwat'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        CMC(N)=DUMSC ! canopy water
         ENDDO

c
c reading SMSTAV
!?      VarName='SMSTAV'
      DO N=NSTART,NEND
        SMSTAV(N)=-9999.
      ENDDO


      VarName='veg'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        VEGFRC(N)=DUMSC ! canopy water
        ENDDO

      VarName='snod'
        DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        ACSNOW(N)=DUMSC
        ACSNOM(N)=-9999.
       ENDDO

	DO L=1,LM
	DO N=NSTART,NEND
         I=IHINDX(N)
         J=JHINDX(N)
         T(N,L)=t_hold(N,L)
 	 OMGA(N,L)=-WH(N,L)*PMID(N,L)*G/
     &              (RD*T(N,L)*(1+.608*Q(N,L)))

!!!!! CONSTRAIN Q TO A PARTICULAR RH VALUE, FOLLOWING CALRH OF WRFPOST
           QC= (PQ0/PMID(N,L)) *EXP(A2*(T(N,L)-A3)/(T(N,L)-A4))
           RH=Q(N,L)/QC

           IF (RH .gt. RHCRIT) THEN
           IF (RH .gt. 1.02) THEN
           write(6,*) 'reducing RH from: ', RH, ' at N,L: ', N,L
	   ENDIF
           Q(N,L)=0.999*RHCRIT*QC
           ENDIF
!!!!! END RH CONSTRAIN

	ENDDO
	ENDDO


      VarName='hgtsfc'
      DO N=NSTART,NEND
       call read_netcdf_2d(ncid_dyn,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        RES(N)=1.0
        FIS(N)=(DUMSC)*G
      ENDDO

C***  PDS IS SURFACE PRESSURE.
!      VarName='pressfc'
      DO N=NSTART,NEND
!       call read_netcdf_2d(ncid_dyn,ifhr,IHINDX(N),JHINDX(N)
!     & ,varname,DUMSC)
        PDS(N)=PSFC(N)
        PINT(N,LM+1)=PSFC(N)
      ENDDO

       call check(nf90_get_att(ncid_dyn,NF90_GLOBAL,"ak",ak))

       PT=ak(1)

	write(6,*) 'returned P_TOP into PT as : ', PT

        DO N=NSTART,NEND
!         I=IHINDX(N)
!         J=JHINDX(N)
!         PINT (N,LM+1)=art(N)+PT
!         PINT (N,1)=PT

	THS(N)=TSHLTR(N)*(100000./PINT(N,LM+1))**CAPA

!!! constrain surface RH

           QC=(PQ0/PINT(N,LM+1))*EXP(A2*(TSHLTR(N)-A3)/(TSHLTR(N)-A4))
           RH=QSHLTR(N)/QC
           IF (RH .gt. RHCRIT) THEN
           write(6,*) 'reducing surface RH from: ', RH, ' at N: ', N
           QSHLTR(N)=0.999*RHCRIT*QC
           ENDIF


        ENDDO


CC
CC RAINC is "ACCUMULATED TOTAL CUMULUS PRECIPITATION"
CC RAINNC is "ACCUMULATED TOTAL GRID SCALE PRECIPITATION"

        write(6,*) 'getting prate_ave'
      VarName='prate_ave'
      DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        ACPREC(N)=DUMSC*3600.*float(ITAG)*0.001
        CUPREC(N)=0.
       ENDDO

	write(0,*) 'maxval(ACPREC): ', maxval(ACPREC)


! XLAND 1 land 2 sea
        varname='land'
       DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
	if (DUMSC .lt. 0.5) then
            SM(N)=1.
            SICE(N)=0.
        elseif (DUMSC .ge. 0.5 .and. DUMSC .lt. 1.5) then
            SM(N)=0.
            SICE(N)=0.
        elseif (DUMSC .gt. 1.5) then
            SM(N)=0.
            SICE(N)=1.
        endif

       ENDDO


!      VarName='HFX'
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)
!      VarName='QFX'
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)
!
!       VarName='LH'
!       call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &   IM,1,JM,1,IM,JS,JE,1)
!      VarName='SNOWC'
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)

!        call ext_ncd_ioclose(DataHandle)

       call check(nf90_close(ncid_dyn))
       call check(nf90_close(ncid_phys))


!!!!!!!!!!!!!!!!! END INSERT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

C------------------------------------------------------------------------
C***
C***  READ QUANTITIES NEEDED FROM THE NHB FILE
C***


!!	ICE available in wrfinput file...zero out for now

      DO N=NSTART,NEND
	SICE(N)=0.0
      ENDDO
C

	write(6,*) 'set LMH to : ', LM
	write(6,*) 'IM,jm: ', Im,jm
	do J=1,JM
	 do I=1,IM
	   LMH(I,J)=LM
	 enddo
	enddo

C

C
C	GET ROTATION ANGLES FOR WIND
C

	write(6,*) 'numsta= ', numsta

	ALLOCATE(CROT(NUMSTA),SROT(NUMSTA))

	CROT=0.	
	SROT=0.

      NTSPH=INT(3600./DT+0.50)



C
C------------------------------------------------------------------------
C

	DO L=1,LM
        DO N=NSTART,NEND
	  TRAIN(N,L)=-9999.
	  TCUCN(N,L)=-9999.
        ENDDO
        ENDDO

      DO N=NSTART,NEND
	Z0(N)=-9999.
        HBOT(N)=-9999.
        CFRACL(N)=-9999.
        CFRACM(N)=-9999.
	CZEN(N)=-9999.
      ENDDO


      DO N=NSTART,NEND
        SNO(N)=-9999. ! many "sno" type variables...which do we need here?
	PSLP(N)=-9999.
	CFRACH(N)=-9999.
	SMSTOT(N)=-9999.
	SFCEXC(N)=-9999.
	CZMEAN(N)=-9999.
        SR(N)=-9999.
      ENDDO

      DO N=NSTART,NEND
        SSROFF(N)=-9999.
        BGROFF(N)=-9999.
       SFCSHX(N)=-9999.
       SFCLHX(N)=-9999.
       SUBSHX(N)=-9999.
       SNOPCX(N)=-9999.
      ENDDO

      DO N=NSTART,NEND
        ASWIN(N)=-9999.
        ASWOUT(N)=-9999.
        ASWTOA(N)=-9999.
        ALWIN(N)=-9999.
        ALWOUT(N)=-9999.
        ALWTOA(N)=-9999.
      ENDDO

      DO N=NSTART,NEND
	POTFLX(N)=-9999.
	TLMIN(N)=-9999.
	TLMAX(N)=-9999.
      ENDDO
C------------------------------------------------------------------------
C***
C***  READ RADIATIVE TEMPERATURE TENDENCIES
C***
      DO L=1,LM
        DO N=NSTART,NEND
         RSWTT(N,L)=-9999.
         RLWTT(N,L)=-9999.
        ENDDO
      ENDDO
C
c     CLOSE(LRSTRT)
C------------------------------------------------------------------------
C***
C***  THE FORECAST HOUR
C***  
c     IFHR=NTSD/NTSPH

	IFHR=ITAG
	write(6,*) 'IFHR: ', IFHR
C------------------------------------------------------------------------
      IF(ITAG.GT.0)THEN
	write(6,*) 'working on preceding file'
C***
C***  GENERATE THE NAME OF THE PRECEDING RESTRT FILE
C***
        ITAG0=ITAG-INCR

	IFHR=ITAG0
	write(6,*) 'older hour IFHR: ', IFHR


!! this needs to change significantly
	write(0,*) 'ITAG0: ', ITAG0

        RINC(5)=0.
        write(fhroldchar,301) ITAG0
        write(0,*) 'fhroldchar: ', fhroldchar
 301  format(i2.2)

        len=index(filedyn,' ')-1
	write(6,*) 'LEN= ', LEN
        filedyn(len-4:len-3)=fhroldchar
	write(6,*) 'filedyn older in PROF: ', trim(filedyn), '_END'

        len=index(filephys,' ')-1
	write(6,*) 'LEN= ', LEN
        filephys(len-4:len-3)=fhroldchar
	write(6,*) 'filephys older in PROF: ', trim(filephys), '_END'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INSERT READ FROM ABOVE ONCE WORK OUT KINKS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call check( nf90_open(filedyn, NF90_NOWRITE, ncid_dyn) )
         call check( nf90_open(filephys, NF90_NOWRITE, ncid_phys) )



C Getting start time

        varname='time'
        iret = nf90_inq_varid(ncid_dyn,trim(varname),varid)
        write(0,*) ncid_dyn,trim(varname),varid
        iret = nf90_get_var(ncid_dyn,varid,ihr)

        iret = nf_get_att_text(ncid_dyn,varid,'units',varin)
        write(0,*) 'varin is: ', varin, ' and end'

        read(varin,101)idate(1),idate(2),idate(3),idate(4),idate(5)


        iyear=IDATE(1)
        imn=IDATE(2)
        iday=IDATE(3)
        ihrst=IDATE(4)


!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"lon1", wbd))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"lat1", sbd))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"dlon", dxval))
!       call check (nf90_get_att(ncid_dyn, NF90_GLOBAL,"dlat", dyval))
!       call check(nf90_get_att(ncid_dyn,NF90_GLOBAL,"cen_lat",cenlat))
!       cenlat=nint(1000.0*cenlat)
!       call check(nf90_get_att(ncid_dyn,NF90_GLOBAL,"cen_lon",cenlon))
!       cenlon=nint(1000.0*cenlon)
!
!

!        write(6,*) 'dxval= ', dxval
!        write(6,*) 'dyval= ', dyval
!        write(6,*) 'cenlat= ', cenlat
!        write(6,*) 'cenlon= ', cenlon

        truelat1=0.
        truelat2=0.
        write(6,*) 'truelat1= ', truelat1
        write(6,*) 'truelat2= ', truelat2
        maptype=5
        write(6,*) 'maptype is ', maptype
!need to get DT
       call check (nf90_get_att(ncid_phys, NF90_GLOBAL,"dtp", DT))
        print*,'DT= ',DT

! get 3-D variables
      print*,'im,jm,lm= ',im,jm,lm
c

        varname='snod'
       DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        ACSNOW0(N)=DUMSC
        ACSNOM0(N)=-9999.
       ENDDO




CC
CC RAINC is "ACCUMULATED TOTAL CUMULUS PRECIPITATION"
CC RAINNC is "ACCUMULATED TOTAL GRID SCALE PRECIPITATION"

        write(6,*) 'getting tprcp'
      VarName='prate_ave'
       DO N=NSTART,NEND
       call read_netcdf_2d(ncid_phys,ifhr,IHINDX(N),JHINDX(N)
     & ,varname,DUMSC)
        ACPREC0(N)=DUMSC*3600.*float(ITAG0)*0.001
        CUPREC0(N)=0.
       ENDDO

!	write(0,*) 'maxval(ACPREC0): ', maxval(ACPREC0)


	write(6,*) 'done reading old file'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END INSERT READ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C
        DO L=1,LM
        DO N=NSTART,NEND
c           TRAIN0(N,L)=DUM(IHINDX(N),JHINDX(N),1)
c           TCUCN0(N,L)=DUM(IHINDX(N),JHINDX(N),2)
	TRAIN0(L,N)=-9999.
	TCUCN0(L,N)=-9999.
          ENDDO
C
c       ENDDO      
        ENDDO
C

!!!
!!!	THESE RUNOFF ASSIGNMENTS COULD BE *WRONG* !!!!!!
!!!

        DO N=NSTART,NEND
          SSROFF0(N)=DUM(IHINDX(N),JHINDX(N),3)
          BGROFF0(N)=DUM(IHINDX(N),JHINDX(N),4)
        ENDDO
C
        DO N=NSTART,NEND
         SFCSHX0(N)=-9999.
         SFCLHX0(N)=-9999.
         SUBSHX0(N)=-9999.
         SNOPCX0(N)=-9999.
        ENDDO
C
        DO N=NSTART,NEND
        ASWIN0(N)=-9999.
        ASWOUT0(N)=-9999.
        ASWTOA0(N)=-9999.
        ALWIN0(N)=-9999.
        ALWOUT0(N)=-9999.
        ALWTOA0(N)=-9999.
          POTFLX0(N)=-9999.
        ENDDO
C
C
      ENDIF
C
c     CLOSE(LRSTRT)

	write(0,*) 'down to here (a)'
C
C
Cmp 	IDEALLY, WON'T NEED MANY MODS BELOW THIS POINT
C
C

C------------------------------------------------------------------------
C***
C***  ALL THE DATA IS NOW IN.
C***  CALCULATE CLOUD FRACTION AND CLOUD WATER/ICE ID NUMBER.
C***
C------------------------------------------------------------------------
      UTIM=1.
      US=1.
      CCLIMIT=1.E-3
      CLIMIT =1.E-20
C-----------------------------------------------------------------------
!$OMP parallel do 
      DO N=NSTART,NEND
        IW(N,1)=-9999
        CLDFRA(N,1)=-9999.
      ENDDO

!	write(6,*) 'here b'
C
C------------------QW, QI AND QINT--------------------------------------
C

!!!
!!! skip section for now
!!!

      DO 220 L=2,LM
C
!$OMP parallel do private(cwmkl,fiq,hh,iwkl,lml,pp,qc,qi,qint,qkl,qw,
!$OMP*                    rqkl,tkl,tmt0,tmt15,u00kl)
      DO 210 N=NSTART,NEND
	IW(N,L)=-9999
        CLDFRA(N,L)=-9999.

!      LML=LM-LMH(IHINDX(N),JHINDX(N))
!	write(6,*) 'LML, IHINDX,JHINDX,LMH: ', IHINDX(N), 
!     &            JHINDX(N),LMH(IHINDX(N),JHINDX(N))
!      TKL=T(N,L)
!      QKL=Q(N,L)
!      CWMKL=CWM(N,L)
!      TMT0=(TKL-273.16)*HH
!      TMT15=AMIN1(TMT0,-15.)*HH
!      AI=0.008855
!      BI=1.
C
!      IF(TMT0.LT.-20.)THEN
!        AI=0.007225
!        BI=0.9674
!      ENDIF
C
!      PP=PMID(N,L)
!      QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
!      QI=QW*(BI+AI*AMIN1(TMT0,0.))
!      QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
!      IF(TMT0.LE.-40.)QINT=QI
C
C-------------------ICE-WATER ID NUMBER IW------------------------------
C

!	no defs for UL
!
!	so no U00KL,FIQ,IW,CWM.  Skip whole section
!
!	write(6,*) 'here c'
!	write(6,*) 'L+LML: ', L+LML
!	write(6,*) 'UL(L+LML): ', UL(L+LML)

!      U00KL=U00(N)+UL(L+LML)*(0.95-U00(N))*UTIM
!      IF(TMT0.LT.-15.)THEN
!        FIQ=QKL-U00KL*QI
!        IF(FIQ.GT.0..OR.CWMKL.GT.CLIMIT)THEN
!          IW(N,L)=1
!        ELSE
!          IW(N,L)=0
!        ENDIF
!      ENDIF
C
!      IF(TMT0.GE.0.)THEN
!        IW(N,L)=0
!      ENDIF
C
!      IF(TMT0.LT.0..AND.TMT0.GE.-15.)THEN
!        IW(N,L)=0
!        IF(IW(N,L-1).EQ.1.AND.CWMKL.GT.CLIMIT)IW(N,L)=1
!      ENDIF
C
!      IWKL=IW(N,L)
C----------------THE SATUATION SPECIFIC HUMIDITY------------------------
!      FIW=FLOAT(IWKL)
!      QC=(1.-FIW)*QINT+FIW*QI
C----------------THE RELATIVE HUMIDITY----------------------------------
!      IF(QC.LE.0.)THEN
!         RQKL=0.
!       ELSE
!         RQKL=QKL/QC
!      ENDIF
C----------------CLOUD COVER RATIO CCR----------------------------------
!      IF(RQKL.GE.0.9999)THEN
!        CCR(N,L)=AMIN1(US,RQKL)
!      ELSE
!        ARG=-1000.*CWMKL/(US-RQKL)
!        ARG=AMAX1(ARG,-25.)
!        CCR(N,L)= RQKL*(1.-EXP(ARG))
!      ENDIF
C----------------------------------------------------------------------
!	write(6,*) 'here d'
  210                 CONTINUE
  220                 CONTINUE
  221                 continue
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C***
C***  BEGIN THE PROFILE POSTING CODE.
C***
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C***
C***  USE ZERO IN ACCUMULATION ARRAYS AT APPROPRIATE TIMES
C***
       IF(ITAG .eq. 0) THEN

!	write(6,*) 'here (2)'
C
C
C what would appropriate if test be here?
C
C
        DO N=NSTART,NEND
C
C*** ZERO ACCUMLATION ARRAYS.
C
          STATPR(N)=0.
          STACPR(N)=0.
          STAEVP(N)=0.
          STAPOT(N)=0.
          STASHX(N)=0.
          STASUB(N)=0.
          STAPCX(N)=0.
          STASWI(N)=0.
          STASWO(N)=0.
          STALWI(N)=0.
          STALWO(N)=0.
          STALWT(N)=0.
          STASWT(N)=0.
          STASNM(N)=0.
          STASNO(N)=0.
          STASRF(N)=0.
          STABRF(N)=0.
          DO L=1,LM
            DHCNVC(L,N)=0.
            DHRAIN(L,N)=0.
          ENDDO
        ENDDO
C
        GO TO 300
       ENDIF
C---------------------------------------------------------------------
C***
C***  WE MUST CHECK TO SEE IF WE ARE 1 HOUR AFTER ANY OF THE 
C***  ACCUMULATION BUCKETS HAVE BEEN EMPTIED.  IF WE ARE AT SUCH A 
C***  TIME THEN WE NEED TO SET TO ZERO THE VARIABLES USED TO HOLD
C***  THE PRECEDING HOUR'S VALUES.
C***
C---------------------------------------------------------------------
C
C
C 	At this point, accumulation buckets are a foreign concept in
C	the WRF model.  
C
C
c     TIME=(NTSD-1)*DT
c     RESET0=TIME-(NTSD/NPREC)*NPREC*DT
c     RESET1=(NPHS-1)*DT+3600.

!	TIME=IFCST
!	write(6,*) 'here (3)'

	RESET0=25.  ! designed to prevent resets.  Reconsider later


C
c      IF(MOD(NTSD,NPREC).GE.NPHS.AND.RESET0.LE.RESET1)THEN
c        DO N=NSTART,NEND
c          STATPR(N)=0.
c          STACPR(N)=0.
c          STASNM(N)=0.
c          STASNO(N)=0.
c          STASRF(N)=0.
c          STABRF(N)=0.
c        ENDDO
c      ELSE
!	write(6,*) 'set STATPR'
        DO N=NSTART,NEND
          STATPR(N)=ACPREC0(N)*1.E3
!	if (ACPREC0(N) .gt. 0) then
!	write(6,*) 'N,ACPREC0(N),STATPR(N): ', N,
!     &			ACPREC0(N),STATPR(N)
!	endif
          STACPR(N)=CUPREC0(N)*1.E3
          STASNM(N)=ACSNOM0(N)*1.E3
          STASNO(N)=ACSNOW0(N)*1.E3
          STASRF(N)=SSROFF0(N)*1.E3
          STABRF(N)=BGROFF0(N)*1.E3
        ENDDO

!	write(6,*) 'past set'
c      ENDIF          
C
c     RESET0=TIME-(NTSD/NRDSW)*NRDSW*DT
c     IF(MOD(NTSD,NRDSW).GE.NPHS.AND.RESET0.LE.RESET1)THEN
c       DO N=NSTART,NEND
c         STASWI(N)=0.
c         STASWO(N)=0.
c         STASWT(N)=0.
c       ENDDO
c     ELSE
        DO N=NSTART,NEND
          STASWI(N)=ASWIN0(N)
          STASWO(N)=ASWOUT0(N)
          STASWT(N)=ASWTOA0(N)
        ENDDO
c     ENDIF
C
c     RESET0=TIME-(NTSD/NRDLW)*NRDLW*DT
c     IF(MOD(NTSD,NRDLW).GE.NPHS.AND.RESET0.LE.RESET1)THEN
c       DO N=NSTART,NEND
c         STALWI(N)=0.
c         STALWO(N)=0.
c         STALWT(N)=0.
c       ENDDO
c     ELSE
        DO N=NSTART,NEND
          STALWI(N)=ALWIN0(N)
          STALWO(N)=ALWOUT0(N)
          STALWT(N)=-ALWTOA0(N)
        ENDDO
c     ENDIF
C
c     RESET0=TIME-(NTSD/NSRFC)*NSRFC*DT
c     IF(MOD(NTSD,NSRFC).GE.NPHS.AND.RESET0.LE.RESET1)THEN
c       DO N=NSTART,NEND
c         STAEVP(N)=0.
c         STAPOT(N)=0.
c         STASHX(N)=0.
c         STASUB(N)=0.
c         STAPCX(N)=0.
c       ENDDO
c     ELSE
        DO N=NSTART,NEND
          STAEVP(N)=SFCLHX0(N)
          STAPOT(N)=POTFLX0(N)
          STASHX(N)=SFCSHX0(N)
          STASUB(N)=SUBSHX0(N)
          STAPCX(N)=SNOPCX0(N)
        ENDDO
c     ENDIF
C
c     RESET0=TIME-(NTSD/NHEAT)*NHEAT*DT
c     IF(MOD(NTSD,NHEAT).GE.NCNVC.AND.RESET0.LE.RESET1)THEN
c       DO N=NSTART,NEND
c         DO L=1,LM
c           DHCNVC(L,N)=0.
c           DHRAIN(L,N)=0.
c         ENDDO
c       ENDDO
c     ELSE
       DO N=NSTART,NEND
         DO L=1,LM
            DHCNVC(L,N)=TCUCN0(L,N)
            DHRAIN(L,N)=TRAIN0(L,N)
         ENDDO
       ENDDO
c     ENDIF
 
C------------------------------------------------------------------
  300 CONTINUE
C------------------------------------------------------------------
C
C***  FOR ROTATION OF WINDS FROM E-GRID TO GEODETIC ORIENTATION
C***  WE NEED THE TWO QUANTITIES BELOW.
C
!      SINPH0=SIN(TPH0D*DTR)
!      COSPH0=COS(TPH0D*DTR)
C
C***  INITIAL CALCULATIONS/PREPARATIONS.  WE LOAD SEVERAL
C***  ARRAYS WITH PROFILE VARIABLES.
C
!$OMP parallel do
      DO N=NSTART,NEND
        IF(CZMEAN(N).GT.0.)THEN
          FACTR(N)=CZEN(N)/CZMEAN(N)
        ELSE
          FACTR(N)=0.
        ENDIF
      ENDDO
C
C***  ADJUST SHORTAVE TENDENCIES TO ACCOUNT FOR CHANGE OF SOLAR POSITION
C***  BETWEEN CALLS TO RADIATION
C
!$OMP parallel do
      DO L=1,LM
        DO N=NSTART,NEND
          RSWTT(N,L)=RSWTT(N,L)*FACTR(N)
        ENDDO
      ENDDO
C
C***  COMPUTE RTOP
C
!$OMP parallel do
      DO L=1,LM
        DO N=NSTART,NEND
          APEL=PMID(N,L)
          RTOP(N,L)=RD*T(N,L)*(1.+0.608*Q(N,L))/APEL
        ENDDO
      ENDDO
C
C***  EGRID2 IS THE SURFACE TEMPERATURE.
C
!$OMP parallel do 
      DO N=NSTART,NEND
        EGRID2(N)= THS(N)*(PDS(N)*1.E-5)**CAPA
        IF(ACPREC(N).LT.0.)ACPREC(N)=0.
        IF(CUPREC(N).LT.0.)CUPREC(N)=0.
      ENDDO
C
C***  SET CYCLE, DATE, AND FORECAST TIME.
C
c      IHR  =NTSD/NTSPH+0.5

         IDATE(2)=imn
         IDATE(3)=iday
         IDATE(1)=iyear
         IDATE(4)=ihrst
         IDATE(5)=ihrst

	IDAT(3)=IDATE(1)
	IDAT(1)=IDATE(2)
	IDAT(2)=IDATE(3)
	
      IYR  =IDAT(3)
      IMNTH=IDAT(1)
      IDAY =IDAT(2)
c      IFCST=(NTSD-1)*DT

!      IFCST=NTSPH*ITAG
      IFCST=3600*ITAG
	IHR=ITAG
	write(6,*) 'IFCST: ', IFCST
C
      WRITE(6,*)' POST PROFILE FOR ',
     1                       IYR,IMNTH,IDAY,IHR
	write(6,*) 'IHRST= ', IHRST
C
C***  SET RTSPH,RTSCU,RTSRA TO 1. OVER THE NUMBER OF TIMES THE
C***  VARIOUS PHYSICS ROUTINES HAVE BEEN
C***  CALLED SINCE LAST OUTPUT OF PROFILER DATA.  NECESSARY FOR
C***  CORRECT AVERAGING OF VARIABLES.
C
C	write(6,*) 'APHTIM, ACUTIM, ARATIM were: ', 
C     &                APHTIM, ACUTIM, ARATIM
      
	APHTIM=0.
	ACUTIM=0.
	ARATIM=0.

      IF(APHTIM.GT.0.)THEN
        RTSPH=1./APHTIM
      ELSE
        RTSPH=1.
      ENDIF
C
      IF(ACUTIM.GT.0.)THEN
        RTSCU=1./ACUTIM
      ELSE
        RTSCU=1.
      ENDIF
C
      IF(ARATIM.GT.0.)THEN
        RTSRA=1./ARATIM
      ELSE
        RTSRA=1.
      ENDIF
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C***
C***  OUTPUT PROFILE DATA.  THE FOLLOWING LOOP IS OVER ALL PROFILE SITES.
C***
C--------------------------------------------------------------------------

!! gather all onto root task?

	write(0,*) 'call collect_all - mype, NSTART: ',mype,nstart

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call collect_all(RES,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(FIS,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(THS,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(HBOT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(CFRACL,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(CFRACM,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(CFRACH,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(SOILTB,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SFCEXC,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(SMSTAV,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SMSTOT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(Z0,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
! skip CZEN, CZMEAN as not written to profdat

      call collect_all(ACPREC,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(CUPREC,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ACSNOW,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ACSNOM,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(SNO,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SR,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(SSROFF,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(BGROFF,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SFCSHX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(SFCLHX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SUBSHX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SNOPCX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(ASWIN,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ASWOUT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ASWTOA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(ALWIN,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ALWOUT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(ALWTOA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(PDS,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(TSHLTR,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(QSHLTR,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(TH10,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(Q10,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(U10,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(V10,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(TLMIN,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(TLMAX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(CMC,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(VEGFRC,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(POTFLX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(PSLP,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(EGRID2,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SM,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(SICE,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)

      call collect_all(STATPR,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STACPR,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STAEVP,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STAPOT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASHX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASUB,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STAPCX,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASWI,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASWO,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASWT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STALWI,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STALWO,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STALWT,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASNM,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASRF,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STABRF,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      call collect_all(STASNO,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)


      do L=1,NSOIL

      VARTMPSTA(:)=SMC(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      SMC(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=STC(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      STC(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=SH2O(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      SH2O(:,L)=VARTMPSTA(:)

      enddo

!	if (MYPE .eq. 0) then
!	do I=1,NUMSTA
!         write(0,*) 'I, STC(I,1:4): ', I, STC(I,1),STC(I,2),
!     *                                    STC(I,3),STC(I,4)
!        enddo
!        endif

       do L=1,LM

       
      VARTMPSTA(:)=T(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      T(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=Q(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      Q(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=U(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      U(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=V(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      V(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=Q2(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      Q2(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=CWM(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      CWM(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=TRAIN(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      TRAIN(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=TCUCN(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      TCUCN(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=RSWTT(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      RSWTT(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=RLWTT(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      RLWTT(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=CLDFRA(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      CLDFRA(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=RTOP(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      RTOP(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=PMID(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      PMID(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=OMGA(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      OMGA(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=t_hold(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      t_hold(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=W(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      W(:,L)=VARTMPSTA(:)

      VARTMPSTA(:)=WH(:,L)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      WH(:,L)=VARTMPSTA(:)


      VARTMPSTA(:)=DHRAIN(L,:)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      DHRAIN(L,:)=VARTMPSTA(:)

      VARTMPSTA(:)=DHCNVC(L,:)
      call collect_all(VARTMPSTA,NSTART,NEND,icnt,idspl,
     &                       NUMSTA,mype,npes,ierr)
      DHCNVC(L,:)=VARTMPSTA(:)

      enddo ! on LM

! only write on root task

	write(0,*) 'task: ', mype, ' to edge of big root task write block'

	if (MYPE .eq. 0) then

	write(0,*) 'task: ', mype, ' IS INSIDE big root task write block'
	LCLAS1=79


	write(6,*) 'open output file with RECL: ', LRECPR
      OPEN(UNIT=LCLAS1,ACCESS='DIRECT',RECL=LRECPR,IOSTAT=IER)
C--------------------------------------------------------------------------
	write(6,*) 'RECORD LENGTH = ', LRECPR

      DO 1000 N=1,NUMSTA
C
C***  ZERO OUTPUT ARRAY.
C
      DO K=1,NWORDM
        PRODAT(K)=0.
        FPACK(K) =0.
      ENDDO
C
C***  CONSTRUCT HEADER FOR CURRENT PROFILE SITE.  THE HEADER CONTAINS
C***  THE FOLLOWING INFORMATION:  PACKED CYCLE-DATE, FORECAST TIME,
C***  INTEGER STATION ID, STATION LATITUDE, STATION LONGITUDE, STATION
C***  ELEVATION, NUMBER OF VERTICAL LEVELS IN PROFILE, NUMBER OF MULTI-
C***  LEVEL PARAMETERS, NUMBER OF SINGLE LEVEL PARAMETERS, TOTAL LENGTH
C***  (IN WORDS) OF MULTI- AND SINGLE LEVEL DATA, PROFILE CLASS FLAG,
C***  AND A DUMMY WORD FOR FUTURE USE.
C
      IH=IHINDX(N)
      JH=JHINDX(N)
      LMHK     = LMH(IH,JH)
      NWORD2   = 2*LMHK
      NWORD3   = 3*LMHK
      NWORD4   = 4*LMHK
      NWORD5   = 5*LMHK
      NWORD6   = 6*LMHK
      NWORD7   = 7*LMHK
      NWORD8   = 8*LMHK
      NWORD9   = 9*LMHK
      NWORD10  = 10*LMHK
      NWORD11  = 11*LMHK
      NWORD12  = 12*LMHK
      NWORD13  = 13*LMHK
      ISTAT    = IDSTN(N)
      CISTAT   = CIDSTN_SAVE(N)
      FPACK(1) = STNLAT(N)/DTR
      FPACK(2) = -STNLON(N)/DTR
      IF(FPACK(2).LT.-180.)FPACK(2)=FPACK(2)+360.
      IF(FPACK(2).GT. 180.)FPACK(2)=FPACK(2)-360.

!	if (mod(N,20) .eq. 0) then
!	write(0,*) 'N, IDSTN,CIDSTN_SAVE, STNLAT, STNLON: ', 
!     &             N, IDSTN(N), CIDSTN_SAVE(N), FPACK(1:2)
!	endif

      FPACK(3) = FIS(N)*GI
      FPACK(4) = FLOAT(LMHK)
      FPACK(5) = LCL1ML -2 !20080708: B Zhou don't store 13 and 14th sounding
      FPACK(6) = LCL1SL
      FPACK(7) = 9+FPACK(5)*FPACK(4)+FPACK(6)
	if (N .eq. 1) then
         write(0,*) 'FPACK(1:7): ', FPACK(1:7)
        endif
      FPACK(8) = 999.
      FPACK(9) = 999.

C
C***  WIND ROTATION SINES AND COSINES
C
C      DLM    = STNLON(N)+TLM0D*DTR
C      XX     = COSPH0*COS(STNLAT(N))*COS(DLM)
C     1        +SINPH0*SIN(STNLAT(N))
C      YY     = -COS(STNLAT(N))*SIN(DLM)
C      TLON   = ATAN(YY/XX)
C      ALPHA  = ASIN(SINPH0*SIN(TLON)/COS(STNLAT(N)))
C      SINALP = SIN(ALPHA)
C      COSALP = COS(ALPHA)

C      SINALP = SROT(N)
C      COSALP = CROT(N)

C
C------------------------------------------------------------------
C***  EXTRACT PRESSURE AND TEMPERATURE PROFILES.
C***  EXTRACT/ROTATE U AND V WIND COMPONENT PROFILES.
C***  EXTRACT SPECIFIC HUMIDITY AND TEMPERATURE TENDENCY.
C***  EXTRACT CLOUD WATER, HEATING DUE TO CONVECTION, LARGE
C***  SCALE RAIN, SHORT WAVE RADIATION, LONG WAVE RADIATION,
C***  AND CLOUD FRACTION.
C------------------------------------------------------------------
C
      DO LV=1,LMHK
        LVL=LMHK-LV+1
        PRODAT(LVL)      = PMID(N,LV)
	if (N .eq. NSTART .or. N .eq. 36) then
	write(6,*) 'PRODAT definition, PMID: ', N,LV,LVL,PMID(N,LV)
	endif
	

        PRODAT(LMHK+LVL) = T(N,LV)

C***  ROTATE WINDS
C
        UT     = U(N,LV)
        VT     = V(N,LV)
!        PRODAT(NWORD2+LVL) = UT*COSALP+VT*SINALP
!        PRODAT(NWORD3+LVL) = VT*COSALP-UT*SINALP

        PRODAT(NWORD2+LVL) = UT
        PRODAT(NWORD3+LVL) = VT

C
        PRODAT(NWORD4+LVL) = Q(N,LV)
C
        IF(RTOP(N,LV).GT.1.E-12) THEN
           PRODAT(NWORD5+LVL) = OMGA(N,LV)
	ENDIF

        IF(IW(N,LV).GT.0.5)THEN
          PRODAT(NWORD6+LVL) = -CWM(N,LV)
        ELSE
          PRODAT(NWORD6+LVL) = CWM(N,LV)
        ENDIF

C
        PRODAT(NWORD7+LVL) = TCUCN(N,LV)
        PRODAT(NWORD8+LVL) = TRAIN(N,LV)
        PRODAT(NWORD9+LVL) = RSWTT(N,LV)
        PRODAT(NWORD10+LVL)= RLWTT(N,LV)
        PRODAT(NWORD11+LVL)= CLDFRA(N,LV)*100.

C
        IF(LV.EQ.1)THEN
          PRODAT(NWORD12+LVL)=Q2(N,LV)
        ELSE
          PRODAT(NWORD12+LVL)=(Q2(N,LV)+Q2(N,LV-1))*0.5
        ENDIF
      ENDDO

C
C***  MODIFY ACCUMLATIONS SO AS TO REPRESENT ACCUMULATED
C***  CHANGE SINCE LAST PROFILE OUTPUT TIME.
C
CGSM  MODIFIED CODE TO ACCOUNT FOR DHCNVC AND DHRAIN BEING
C       COMPUTED FROM TOP DOWN WHILE PRODAT IS FILLED FROM
C       BOTTOM UP 
C

      DO LL=1,LMHK
        LVL=LMHK-LL+1
        STADHC(LL) = PRODAT(NWORD7+LL) - DHCNVC(LVL,N)
        STADHR(LL) = PRODAT(NWORD8+LL) - DHRAIN(LVL,N)
C
        DHCNVC(LVL,N) = PRODAT(NWORD7+LL)
        DHRAIN(LVL,N) = PRODAT(NWORD8+LL)
C
Ctmp        IF(MOD(NTSD,NHEAT).LT.NCNVC)THEN
          DHCNVC(LVL,N) = 0.
          DHRAIN(LVL,N) = 0.
Ctmp        ENDIF
      ENDDO
C
C***  EXTRACT SINGLE LEVEL DATA.   EGRID2 IS SURFACE TEMPERATURE.
C
      PRODAT(NWORD13+1)  = PSLP  (N)
      PRODAT(NWORD13+2)  = PDS   (N)
      PRODAT(NWORD13+3)  = EGRID2(N)
      PRODAT(NWORD13+4)  = TLMIN (N)
      PRODAT(NWORD13+5)  = TLMAX (N)
      PRODAT(NWORD13+6)  = SMSTAV(N)*100.
      PRODAT(NWORD13+7)  = ACPREC(N)*1000.
      PRODAT(NWORD13+8)  = CUPREC(N)*1000.
      PRODAT(NWORD13+27) = Z0    (N)
C
      STAPRX=PRODAT(NWORD13+7)-STATPR(N)
      STACRX=PRODAT(NWORD13+8)-STACPR(N)

	if (STAPRX .gt. 5.) then
	write(6,*) '1hr precip: ',  N,STAPRX
	endif
C
C***  ROTATE WINDS
C
      UT     = U10(N)
      VT     = V10(N)
!      PRODAT(NWORD13+28) = UT*COSALP+VT*SINALP
!      PRODAT(NWORD13+29) = VT*COSALP-UT*SINALP

      PRODAT(NWORD13+28) = UT
      PRODAT(NWORD13+29) = VT
C
      PRODAT(NWORD13+30) = TH10  (N)
      PRODAT(NWORD13+31) = Q10   (N)
      PRODAT(NWORD13+32) = TSHLTR(N)
      PRODAT(NWORD13+33) = QSHLTR(N)
      PRODAT(NWORD13+34) = SFCEXC(N)
      PRODAT(NWORD13+35) = VEGFRC(N)
      PRODAT(NWORD13+36) = CMC   (N)*1000.
      PRODAT(NWORD13+37) = SMC   (N,1)
      PRODAT(NWORD13+38) = SMC   (N,2)
      PRODAT(NWORD13+39) = SMC   (N,3)
      PRODAT(NWORD13+40) = SMC   (N,4)
      PRODAT(NWORD13+41) = STC   (N,1)
      PRODAT(NWORD13+42) = STC   (N,2)
      PRODAT(NWORD13+43) = STC   (N,3)
      PRODAT(NWORD13+44) = STC   (N,4)
      PRODAT(NWORD13+45) = SM    (N) + SICE(N)
      PRODAT(NWORD13+46) = CFRACL(N)*100.
      PRODAT(NWORD13+47) = CFRACM(N)*100.
      PRODAT(NWORD13+48) = CFRACH(N)*100.
      PRODAT(NWORD13+49) = SR    (N)*100.
      PRODAT(NWORD13+50) = NINT(HBOT(N))
C
      PRODAT(NWORD13+9)   = SFCLHX(N)
      PRODAT(NWORD13+10)  = POTFLX(N)
      PRODAT(NWORD13+11)  = SFCSHX(N)
      PRODAT(NWORD13+12)  = SUBSHX(N)
      PRODAT(NWORD13+13)  = SNOPCX(N)
      PRODAT(NWORD13+14)  = ASWIN (N)
      PRODAT(NWORD13+15)  = ASWOUT(N)
      PRODAT(NWORD13+16)  = ALWIN (N)
      PRODAT(NWORD13+17)  = ALWOUT(N)
      PRODAT(NWORD13+18)  =-ALWTOA(N)
      PRODAT(NWORD13+19)  = ASWTOA(N)
      PRODAT(NWORD13+20)  = ACSNOW(N)*1000.
      PRODAT(NWORD13+21)  = SMSTOT(N)*1000.
      PRODAT(NWORD13+22)  = SNO   (N)*1000.
      PRODAT(NWORD13+23)  = ACSNOM(N)*1000.
      PRODAT(NWORD13+24)  = SSROFF(N)*1000.
      PRODAT(NWORD13+25)  = BGROFF(N)*1000.
      PRODAT(NWORD13+26)  = SOILTB(N)
C
C***  ACCUMULATED CHANGE SINCE LAST PROFILE OUTPUT TIME.
C
      PSFCEVP  = PRODAT(NWORD13+9 ) - STAEVP(N)
      PPOTEVP  = PRODAT(NWORD13+10) - STAPOT(N)
      PSFCSHX  = PRODAT(NWORD13+11) - STASHX(N)
      PSFCSUB  = PRODAT(NWORD13+12) - STASUB(N)
      PSNOPCX  = PRODAT(NWORD13+13) - STAPCX(N)
      PRSWIN   = PRODAT(NWORD13+14) - STASWI(N)
      PRSWOUT  = PRODAT(NWORD13+15) - STASWO(N)
      PRLWIN   = PRODAT(NWORD13+16) - STALWI(N)
      PRLWOUT  = PRODAT(NWORD13+17) - STALWO(N)
      PRLWTOA  = PRODAT(NWORD13+18) - STALWT(N)
      PRSWTOA  = PRODAT(NWORD13+19) - STASWT(N)
      PACSNOW  = PRODAT(NWORD13+20) - STASNO(N)
      PACSNOM  = PRODAT(NWORD13+23) - STASNM(N)
      PSSROFF  = PRODAT(NWORD13+24) - STASRF(N)
      PBGROFF  = PRODAT(NWORD13+25) - STABRF(N)
C***
C***  TRANSFER STATION PROFILE DATA TO "PACKED" OUTPUT ARRAY.
C***
      NN   = 0
      NLEN = FPACK(7)
C	write(6,*) 'NWORD13+41,NWORD13+32 ', NWORD13+41,NWORD13+32
C	write(6,*) 'SOIL TEMP ', PRODAT(NWORD13+41)
C        write(6,*) 'SHELT TEMP ', PRODAT(NWORD13+32) 
C
      DO NL = 10,NLEN
        NN = NL-9
        FPACK(NL) = PRODAT(NN)
      ENDDO
C
C***  REPLACE ACCUMULATED QUANTITIES WITH ACCUMULATION
C***  SINCE LAST PROFILE OUTPUT TIME.
C
      DO LL = 1,LMHK
!        FPACK(9+NWORD7+LL) = STADHC(LL)*RTSCU
!        FPACK(9+NWORD8+LL) = STADHR(LL)*RTSRA
        FPACK(9+NWORD7+LL) = -9999.
        FPACK(9+NWORD8+LL) = -9999.
      ENDDO
C
      FPACK(9+NWORD13+7)  = STAPRX
!	write(6,*) 'precip written to FPACK element: ', 9+NWORD13+7
      FPACK(9+NWORD13+8)  = STACRX
      FPACK(9+NWORD13+9)  = PSFCEVP * RTSPH
      FPACK(9+NWORD13+10) = PPOTEVP * RTSPH
      FPACK(9+NWORD13+11) = PSFCSHX * RTSPH
      FPACK(9+NWORD13+12) = PSFCSUB * RTSPH
      FPACK(9+NWORD13+13) = PSNOPCX * RTSPH
      FPACK(9+NWORD13+14) = PRSWIN  * RTSPH
      FPACK(9+NWORD13+15) = PRSWOUT * RTSPH
      FPACK(9+NWORD13+16) = PRLWIN  * RTSPH
      FPACK(9+NWORD13+17) = PRLWOUT * RTSPH
      FPACK(9+NWORD13+18) = PRLWTOA * RTSPH
      FPACK(9+NWORD13+19) = PRSWTOA * RTSPH
      FPACK(9+NWORD13+20) = PACSNOW
      FPACK(9+NWORD13+23) = PACSNOM
      FPACK(9+NWORD13+24) = PSSROFF
      FPACK(9+NWORD13+25) = PBGROFF
C
!      IF(RESTRT)THEN
      IF(ITAG .eq. 0)THEN
        DO LL = 1,LMHK
          FPACK(9+NWORD7+LL) = 0.
          FPACK(9+NWORD8+LL) = 0.
        ENDDO
C
        FPACK(9+NWORD13+7)  = 0.
        FPACK(9+NWORD13+8)  = 0.
        FPACK(9+NWORD13+9)  = 0.
        FPACK(9+NWORD13+10) = 0.
        FPACK(9+NWORD13+11) = 0.
        FPACK(9+NWORD13+12) = 0.
        FPACK(9+NWORD13+13) = 0.
        FPACK(9+NWORD13+14) = 0.
        FPACK(9+NWORD13+15) = 0.
        FPACK(9+NWORD13+16) = 0.
        FPACK(9+NWORD13+17) = 0.
        FPACK(9+NWORD13+18) = 0.
        FPACK(9+NWORD13+19) = 0.
        FPACK(9+NWORD13+20) = 0.
        FPACK(9+NWORD13+23) = 0.
        FPACK(9+NWORD13+24) = 0.
        FPACK(9+NWORD13+25) = 0.
      ENDIF
C---------------------------------------------------------------------
C***
C***  WRITE PROFILE DATA
C***
      
!       NREC=(IFHR/INCR)*NUMSTA+N
       NREC=N

      WRITE(LCLAS1,REC=NREC)IHRST,IDAT,IFCST,ISTAT,CISTAT
     1,                    (FPACK(NL),NL=1,NLEN)

	if (N .le. 1) then
	write(0,*) 'IHRST,IDAT,IFCST,ISTAT,CISTAT: ', 
     &          IHRST,IDAT,IFCST,ISTAT,CISTAT

	do I=1,NLEN
	write(0,*) 'N, I, FPACK(I): ', N, I, FPACK(I)
        enddo

	endif



C---------------------------------------------------------------------
 1000 CONTINUE
      CLOSE(LCLAS1)
      endif ! mype=0

	write(0,*) 'this task: ', mype, ' is to the end deallocs'

	DEALLOCATE(T,Q,U,V,Q2,CWM,TRAIN,TCUCN)
	DEALLOCATE(RSWTT,RLWTT,CLDFRA,RTOP,OMGA)
	DEALLOCATE(t_hold,PINT)

        DEALLOCATE(DHCNVC,DHRAIN,STADHC,STADHR,TCUCN0,TRAIN0)
        DEALLOCATE(DUM,DUMMY,DUMMY2,DUM3D,DUM3D2,DUM3D3)
        DEALLOCATE(PRODAT,FPACK,LMH,VARTMP)
        DEALLOCATE(CIDSTN,STNLAT,STNLON)
        DEALLOCATE(IDSTN,IHINDX,JHINDX,IVINDX,JVINDX)

        DEALLOCATE(
     & RES,FIS,THS,HBOT
     &,CFRACL,CFRACM,CFRACH
     &,SOILTB,SFCEXC,SMSTAV,SMSTOT
     &,Z0,CZEN,CZMEAN
     &,ACPREC,CUPREC,ACSNOW,ACSNOM
     &,SSROFF,BGROFF,SFCSHX,SFCLHX
     &,SUBSHX,SNOPCX,ASWIN,ASWOUT
     &,ASWTOA,ALWIN,ALWOUT,ALWTOA
     &,TSHLTR,QSHLTR
     &,TH10,Q10,U10,V10
     &,TLMIN,TLMAX
     &,SMC,CMC,STC,SH2O
     &,VEGFRC,POTFLX,PSLP,PSFC
     &,EGRID2,SM,SICE
     &,FACTR
     &,STATPR,STACPR,STAEVP
     &,STAPOT,STASHX,STASUB,STAPCX
     &,STASWI,STASWO,STALWI,STALWO
     &,STALWT,STASWT,STASNM,STASRF
     &,STABRF,STASNO,SNO,SR
     &,ACPREC0,CUPREC0,SFCLHX0,POTFLX0
     &,SFCSHX0,SUBSHX0,SNOPCX0,ASWIN0
     &,ASWOUT0,ALWIN0,ALWOUT0,ALWTOA0
     &,ASWTOA0,ACSNOW0,ACSNOM0,SSROFF0
     &,BGROFF0,ICNT,IDSPL)


      call mpi_finalize(ierr)
C
C***  END OF PROFILE SITE LOOP
C
C***  END PROFILE POSTING CODE.
C---------------------------------------------------------------------
      END

      subroutine check(status)
      use netcdf
      integer, intent ( in) :: status

       if(status /= nf90_noerr) then
         print*, 'trouble... status: ', status
         stop "Stopped"
       end if
      end subroutine check



! --------------------------------


      subroutine read_netcdf_3d(ncid,ifhr,
     & VarName,itarg,jtarg,LM,buf) 

      use netcdf
      implicit none
      character(len=31),intent(in) :: VarName
      integer,intent(in) :: ncid,ifhr,itarg,jtarg,LM
      real,intent(out)   :: buf(LM)
      integer            :: iret,i,j,jj,varid
      real,parameter     :: spval_netcdf=-1.e+10

        iret = nf90_inq_varid(ncid,trim(varname),varid)

        iret = nf90_get_var(ncid,varid,buf,
     &       start=(/itarg,jtarg,1,1/),
     &       count=(/1,1,LM,1/))
        if (iret /= 0) then
          print*,VarName," not found ", IRET
          buf=-9999.
        end if

      end subroutine read_netcdf_3d

! --------------------------------

      subroutine read_netcdf_2d(ncid,ifhr,itarg,jtarg,
     &  VarName,bufvar)

      use netcdf
      implicit none
      character(len=31),intent(in) :: VarName
      integer,intent(in) :: ncid,ifhr,itarg,jtarg
      real,intent(out)   :: bufvar(1,1)
      integer            :: iret,i,j,jj,varid

        iret = nf90_inq_varid(ncid,trim(varname),varid)
!        write(0,*) ncid,trim(varname),varid

        iret = nf90_get_var(ncid,varid,bufvar,
     &       start=(/itarg,jtarg,1/),
     &       count=(/1,1,1/))

        if (iret /= 0) then
          print*,VarName, " not found -Assigned missing values"
              bufvar = -9999.
        endif


      end subroutine read_netcdf_2d

! --------------------------------

      subroutine collect_all(var,buf_start,buf_end,icnt,idspl,
     &                       ntot,mype,npes,ierr)
     
      include 'mpif.h'
      real var(ntot),varout(ntot)
      integer :: icnt(0:NPES-1),idspl(0:NPES-1),ierr
      integer :: buf_start, buf_end, ntot, mype, npes

!	write(0,*) 'size(var): ', size(var)
!        write(0,*) 'buf_start in collect_all: ', buf_start
!	write(0,*) 'mype, icnt(MYPE),idspl(MYPE): ', 
!     &             mype, icnt(MYPE),idspl(MYPE)
!	write(0,*) 'what is MPI_COMM_WORLD: ', MPI_COMM_WORLD

      call mpi_gatherv(var(buf_start),icnt(mype),MPI_REAL,
     &            varout,icnt,idspl,MPI_REAL,0,MPI_COMM_WORLD,ierr)

!	write(0,*) 'mype, ierr from mpi_gatherv: ', mype, ierr

	if (mype == 0) then
          var=varout
        endif

      end subroutine collect_all

! --------------------------------
