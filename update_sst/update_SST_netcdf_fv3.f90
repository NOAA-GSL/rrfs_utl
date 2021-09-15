subroutine update_SST_netcdf_fv3 (sstRR, glat, glon, nlon, nlat, xland, vegtyp, ilake,iice,iyear,imonth,iday,ihour)
!$$$  documentation block
!                .      .    .                                       .
!   update_SST_netcdf_mass: read SST from wrf mass netcdf old background file
!           and update SST in wrf mass background file
!   prgmmr: Ming Hu                 date: 2008-02-25
!   updates: Tanya Smirnova         date: 2010-10-11
!
! program history log:
!
!
!   input argument list:
!       sstRR: sst
!       nlon:  x dimension
!       nlat:  Y dimension
!
! attributes:
!   language: f90
!
!$$$

  use kinds, only: r_single,i_kind, r_kind
  use mpi
  use constants, only: init_constants,init_constants_derived
  use constants, only: rd,h1000,rd_over_cp,grav,half
  use gsi_rfv3io_sst_mod, only: gsi_rfv3io_get_grid_specs
  use gsi_rfv3io_sst_mod, only: bg_fv3regfilenameg,fv3sar_bg_opt
  use gsi_rfv3io_sst_mod, only: rfv3io_mype
  use gsi_rfv3io_sst_mod, only: gsi_fv3ncdf_read,gsi_fv3ncdf2d_read
  use gsi_rfv3io_sst_mod, only: gsi_fv3ncdf_write,gsi_fv3ncdf_append
  use gsi_rfv3io_sst_mod, only: gsi_fv3ncdf_append2d
  use gsi_rfv3io_sst_mod, only: nlon_regional,nlat_regional,nsig_regional
  use gsi_rfv3io_sst_mod, only: eta1_ll

  implicit none

!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror

  INCLUDE 'netcdf.inc'

  integer, parameter :: WRF_INTEGER = 106
!
  integer :: nlon, nlat, ilake, iice
  real  :: sstRR(nlon,nlat)
  real  :: glat(nlon,nlat)
  real  :: glon(nlon,nlat)
  real  :: xland(nlon,nlat) ! tgs is it slmsk in FV3? xland(water)=0
  real  :: vegtyp(nlon,nlat)

! Declare local parameters

  character(len=120) :: flnm1
  
  integer(i_kind) :: i,j,k
  
  character (len=80), dimension(3)  ::  dimnames
  
  integer(i_kind) :: l, n
  
  integer(i_kind) :: ierr, ier, Status, Status_next_time
  character(len=:),allocatable :: sfcvars   !='sfc_data.nc'

!  character (len=31) :: rmse_var
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
!  integer(i_kind) nlon_regional,nlat_regional,nsig_regional
  real(r_single),allocatable::field2(:,:)
  real(r_single),allocatable::field3(:,:,:)
  real(r_single),allocatable::surftemp(:,:)
  real(r_single),allocatable::temp2m(:,:)
  real(r_single),allocatable::sst(:,:)
  real(r_single),allocatable::lu_index(:,:)
  !tgs FV3 lake variables are lake_depth and lake_fraction
  real(r_single),allocatable::lake_frac(:,:)
  real(r_single),allocatable::lake_depth(:,:)
  real(r_single),allocatable::tsea(:,:)
  real(r_single),allocatable::tsfc(:,:)

  real(r_single)    :: time, time1, time2
  real(r_single)    :: a, b

! Lakes from RUC model
! -- Great Salt Lake lake surface temps
        REAL salt_lake_lst (13)
        data salt_lake_lst &
        /1.,3.,6.,13.,17.,20.,26.,25.,20.,14.,9.,3.,1./
! -- Salton Sea - California
        REAL salton_lst (13)
        data salton_lst &
        / 12.8, 12.8, 17.2, &
         21.1, 24.4, 26.7,  &
         30.0, 31.7, 29.4,  &
         25.0, 20.6, 15.0,  &
         12.8/
! -- Lake Champlain - Vermont
        REAL champ_lst (13)
        data champ_lst&
       /  1.3,  0.6,  1.0,  &
          3.0,  7.5, 15.5,  &
         20.5, 21.8, 18.2,  &
         13.0,  8.2,  4.5,  &
          1.3/

        real xc1,yc1, xc2,yc2
        integer isup,jsup, iwin,jwin, isalton,jsalton

      integer julm(13)
      data julm/0,31,59,90,120,151,181,212,243,273,304,334,365/

        INTEGER  mon1,mon2,day1,day2,juld
        real rday,wght1,wght2
!20aug18 - lake model
        integer :: NCID,in_SF_LAKE_PHYSICS
!
        integer :: nlon_regional_smc, nlat_regional_smc, nsig_regional_smc

  
  nlon_regional=nlon
  nlat_regional=nlat
!
! ===============================================================================
!

  flnm1='sfc_data.nc'

!tgs - in FV3 lkm=1 is turning on the lake model, we'll do a change later when
!      lake model is on
  in_SF_LAKE_PHYSICS=0
!
!-------------  get grid info

!  open and read background dimesion

  fv3sar_bg_opt=1    ! 0=restart, 1=input
  call bg_fv3regfilenameg%init
  mype=0

  sfcvars= bg_fv3regfilenameg%sfcdata

!  read in background fields
  call gsi_rfv3io_get_grid_specs(bg_fv3regfilenameg,ierror)  

!
  write(6,*) '================================================='
  sfcvars="./sfc_data.nc"
  write(*,*)"nlon_regional, nlat_regional", nlon_regional, nlat_regional

  allocate(surftemp(nlon_regional,nlat_regional))
  allocate(tsea(nlon_regional,nlat_regional))
  allocate(tsfc(nlon_regional,nlat_regional))
  allocate(temp2m(nlon_regional,nlat_regional))
  allocate(sst(nlon_regional,nlat_regional))
  allocate(lu_index(nlon_regional,nlat_regional))
!tgs - so far we do not have lake info in NA RRFS. The variables in FV3 are
!      lake_depth and lake_frac
  allocate(lake_frac(nlon_regional,nlat_regional)) 
  allocate(lake_depth(nlon_regional,nlat_regional))

  allocate(field2(nlon_regional,nlat_regional))
!
  write(6,*) '================================================='
  call gsi_fv3ncdf2d_read(sfcvars,'tsea','TSEA',field2,mype)
  sst=field2(:,:)
  write(*,*)'Done with reading TSEA from file ',sfcvars
  write(6,*)' max,min bck skin temp tsea(K)=',maxval(sst),minval(sst)
       write(6,*)'background skin temp tsea(170,170)', sst(170,170)
       write(6,*)'new  sstRR(170,170)', sstRR(170,170)
       write(6,*)'Winnipeg skin temp tsea(516,412)', sst(516,412)
       write(6,*)'Winnipeg sstRR(516,412)', sstRR(516,412)
!
  write(6,*) '================================================='
  call gsi_fv3ncdf2d_read(sfcvars,'tsfc','TSFC',field2,mype)
  surftemp=field2(:,:)
  write(*,*)'Done with reading TSFC from file ',sfcvars
  write(6,*)' max,min bck skin temp tsfc(K)=',maxval(surftemp),minval(surftemp)
       write(6,*)'background skin temp tsfc(170,170)', surftemp(170,170)
       write(6,*)'Winnipeg skin temp tsfc(516,412)', surftemp(516,412)
!
  write(6,*) '================================================='
  call gsi_fv3ncdf2d_read(sfcvars,'t2m','T2M',field2,mype)
  temp2m=field2(:,:)
  write(*,*)'Done with reading T2M from file ',sfcvars
  write(6,*)' max,min bck 2m temp (K)=',maxval(temp2m),minval(temp2m)
       write(6,*)'background t2 temp(292,258)', temp2m(292,258)
       write(6,*)'new  sstRR(170,170)', sstRR(292,258)
!
  write(6,*) '================================================='
!  rmse_var='LAKE_DEPTH'
  !lake_depth=0
  !write(6,*)' max,min bck lake_depth (K)=',maxval(lake_depth),minval(lake_depth)
!
  write(6,*) '================================================='
!  rmse_var='VTYPE'
  call gsi_fv3ncdf2d_read(sfcvars,'vtype','VTYPE',field2,mype)
  lu_index=field2(:,:)
  write(*,*)'Done with reading VTYPE (LU_INDEX) from file ',sfcvars
  write(6,*)' max,min VTYPE (LU_INDEXi)=',maxval(lu_index),minval(lu_index)
!
  write(6,*) '================================================='
!  rmse_var='LAKE_FRAC'
  lake_frac=0
  !write(6,*)' max,min LAKEMASK=',maxval(lake_frac),minval(lake_frac)
!
! Compute weight for the current date
       juld = julm(imonth) + iday
       if(juld.le.15) juld=juld+365

       mon2 = imonth
       if(iday.gt.15) mon2 = mon2 + 1
       if(mon2.eq.1) mon2=13
       mon1=mon2-1
! **** Assume data valid at 15th of month
       day2=julm(mon2)+15
       day1=julm(mon1)+15
       rday=juld
       wght1=(day2-rday)/float(day2-day1)
       wght2=(rday-day1)/float(day2-day1)
       write(6,*)'Date weights =',wght1,wght2

!
!  update skin temperature over water
!
if(1==1) then  ! turn off , use GFS SST
! find i,j for a point in northern Lake Superior
  DO J=1,nlat
  DO I=1,nlon
   if((glat(i,j)>48.4 .and. glat(i,j)<49.6) .and. (glon(i,j)<-87.9 .and. glon(i,j)>-88.1)) then
     isup=i
     jsup=j
     print *,' Lake Superior --> i,j,glat(i,j),glon(i,j)',i,j,glat(i,j),glon(i,j), &
     'vegtyp(i,j)=',vegtyp(i,j),'lu_index(i,j)',lu_index(i,j),xland(i,j)
     goto 99
   endif
  ENDDO
  ENDDO

99  continue

! find i,j for a point in northern Lake Winnipeg
  DO J=1,nlat
  DO I=1,nlon
   if((glat(i,j)>53.3 .and. glat(i,j)<53.7) .and. (glon(i,j)<-98.3 .and.  glon(i,j)>-98.7)) then
     iwin=i
     jwin=j
     print *,' Lake Winnipeg --> i,j,glat(i,j),glon(i,j)',i,j,glat(i,j),glon(i,j), &
     'vegtyp(i,j)=',vegtyp(i,j),'lu_index(i,j)',lu_index(i,j),xland(i,j)
     goto 999
   endif
  ENDDO
  ENDDO

999  continue

  write(*,*) 'in_SF_LAKE_PHYSICS:', in_SF_LAKE_PHYSICS
       
  DO J=1,nlat
  DO I=1,nlon
    if( xland(i,j) < 0.00001 ) then    ! water, xland = 0
    ! only unfrozen water points (sea or lakes)
       if(in_SF_LAKE_PHYSICS == 0 .and. lake_depth(i,j) > 0.) then
       ! --- CLM lake model is off, use climatology for several lakes
       ! --- Great Salt Lake, Utah Lake -- Utah
            if (glat(i,j).gt.39.5 .and. glat(i,j).lt.42. .and.  &
               glon(i,j).gt.-114..and. glon(i,j).lt.-111.) then
               write(6,*)'Global data Salt Lake temp',i,j,sstRR(i,j)
               sstRR(i,j) = 273.15 + wght1*salt_lake_lst(mon1)  &
                       +wght2*salt_lake_lst(mon2)
                write(6,*)'Climatology Salt Lake temp',i,j,sstRR(i,j)  &
                          ,glat(i,j),glon(i,j)
            end if

            ! --- Salton Sea -- California
            if (glat(i,j).gt.33. .and. glat(i,j).lt.33.7 .and.  &
                glon(i,j).gt.-116.3 .and. glon(i,j).lt.-115.3) then
            write(6,*)'Global data Salton Sea temp',i,j,sstRR(i,j)
            sstRR(i,j) = 273.15 + wght1*salton_lst(mon1)  &  
                       +wght2*salton_lst(mon2)
            write(6,*)'Climatology Salton Sea temp',i,j,sstRR(i,j)  &
                ,glat(i,j),glon(i,j)
              isalton=i
              jsalton=j
            end if

            ! --- Lake Champlain -- Vermont
            if (glat(i,j).gt.44. .and. glat(i,j).lt.45.2 .and.  &
               glon(i,j).gt.-74. .and. glon(i,j).lt.-73.) then
            write(6,*)'Global data Lake Champlain temp',i,j,sstRR(i,j)
            sstRR(i,j) = 273.15 + wght1*champ_lst(mon1)  &
                       +wght2*champ_lst(mon2)
            write(6,*)'Climatology Lake Champlain temp',i,j,sstRR(i,j)  &
                ,glat(i,j),glon(i,j)
            end if
            ! --- For Lake Nipigon, use point for n. Lake Superior
            !   -- Lake Nipigon is deep!
            if (glat(i,j).gt.49. .and. glat(i,j).lt.51. .and. &
               glon(i,j).gt.-90. .and. glon(i,j).lt.-87.) then
               write(6,*)'Global data Lake Nipigon temp',i,j,sstRR(i,j)
                sstRR(i,j) = sstRR(isup,jsup)
                write(6,*)'Lake Nipigon temp',i,j,sstRR(i,j) &
                 ,glat(i,j),glon(i,j)
            end if

     if(1 == 2) then
! --- For Lake of the Woods and other
!      Minnesota lakes, use point for n. Lake Winnipeg
!    -- These lakes are NOT DEEP!
            if (glat(i,j).gt.46. .and. glat(i,j).lt.50. .and.  &
               glon(i,j).gt.-96. .and. glon(i,j).lt.-93.) then
                write(6,*)'Global data Minnesota lake temp',i,j,sstRR(i,j)
                sstRR(i,j) = sstRR(iwin,jwin)
                write(6,*)'Minnesota lake temp',i,j,sstRR(i,j) &
                 ,glat(i,j),glon(i,j)
            end if

! --- For Canadian lakes, including Winnipeg, Manitoba, Winnipegosis,
!      use point for n. Lake Winnipeg
!    -- These lakes are NOT DEEP!
            if (glat(i,j).gt.50. .and. glat(i,j).lt.68. .and.  &
               glon(i,j).gt.-148. .and. glon(i,j).lt.-48.) then
                write(6,*)'Global data Canadian lake temp',i,j,sstRR(i,j)
                sstRR(i,j) = sstRR(iwin,jwin)
                write(6,*)'Canadian lake temp',i,j,sstRR(i,j) &
                 ,glat(i,j),glon(i,j)
            end if
! --- For lakes in Washington, Oregon, Nevada
!      use point for n. Lake Winnipeg (?????)
!    -- These lakes are NOT DEEP!
            if (glat(i,j) > 33.8 .and. glat(i,j) < 50. .and.  &
               glon(i,j) < -114. ) then
                write(6,*)'Global data US west lake temp',i,j,sstRR(i,j)
                sstRR(i,j) = sstRR(iwin,jwin)
!                sstRR(i,j) = sstRR(isalton,jsalton)
                write(6,*)'US west lake temp',i,j,sstRR(i,j) &
                 ,glat(i,j),glon(i,j)
            end if

     endif ! 1 == 2
      endif   ! lakes

!-- update skin temp for water points
! in_SF_LAKE_PHYSICS == 0 --> lakemask == 0, SST is updated at all water points
! in_SF_LAKE_PHYSICS == 1 --> lakemask == 1 where CLM lake model is applied, SST is
! not updated at these points.
           if( lake_frac(i,j) == 0.) then
               surftemp(i,j) = sstRR(i,j)
           endif

           if(sstRR(i,j) > 400. .or. sstRR(i,j) < 100. ) then
             print *,'Bad SST at point i,j',i,j,sstRR(i,j)
           endif

    elseif (xland(i,j) > 1.99) then ! ice, xland = 2
    ! frozen water - xland = 2
          if(in_SF_LAKE_PHYSICS == 0) then
          ! CLM lake model is turned off
            if( lake_depth(i,j) > 0. .and. sstRR(i,j) > 273.) then
            ! -- frozen lakes with CLM lake model turned off.
              print *,'Ice lake cannnot have SST > 274K'
              !print *,'i,j,sstRR,lu_index,vegtyp =' ,i,j,sstRR(i,j),lu_index(i,j),vegtyp(i,j)

              !set skin temp of frozen lakes to 2-m temperature
              sstRR(i,j)= min(273.15,temp2m(i,j))

              !update skin temp for frozen lakes
              surftemp(i,j)=sstRR(i,j)
            endif
          endif ! in_SF_LAKE_PHYSICS == 0

    endif  ! water or ice

    ! update SST 
    sst(i,j)=sstRR(i,j)
  ENDDO
  ENDDO
  write(*,*) 'Skin temperature updated with current SST'
       write(6,*)' updated skin temp(170,170)', surftemp(170,170)
       write(6,*)' max,min skin temp =',maxval(surftemp),minval(surftemp)
endif  ! 1==1, when 1==2 SST update is turned off
!
!
!           update fv3 netcdf file with new SST
!
  write(6,*) ' ============================= '
  write(6,*) ' update SST in background file '
  write(6,*) ' ============================= '
  flnm1='sfc_data.nc'
  sfcvars='./sfc_data.nc'
     
!-------------  get date info

  write(6,*) ' Update SST in background at time:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour

  write(6,*) '================================================='
  write(6,*)' max,min skin temp =',maxval(surftemp),minval(surftemp)
  call gsi_fv3ncdf_append2d(sfcvars,'tsea',sst,mype)
  call gsi_fv3ncdf_append2d(sfcvars,'tsfc',surftemp,mype)
  deallocate(field2)

  call MPI_Barrier(mpi_comm_world, ierror)
  close(6)
  if(mype==0)  write(*,*) "=== FV3 UPDATESST SUCCESS ==="

end subroutine update_SST_netcdf_fv3


