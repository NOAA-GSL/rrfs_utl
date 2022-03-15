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
!                     vegcat(i,j)=24     not used
!                     ivgtyp(i,j)=24     int IVGTYP / =15 for MODIS
!                     lu_index(i,j)=24   float LU_INDEX / =15 for MODIS
!                     landmask(i,j)=1.   float LANDMASK
!                     xland_rr(i,j)=1.      float XLAND
!                     isltyp(i,j)=16.    int ISLTYP
! 
! For water:
!                     vegcat(i,j)=16   / =17 and 21 (inland) for MODIS
!                     ivgtyp(i,j)=16  / =17 and 21 for MODIS
!                     lu_index(i,j)=16  / =17 and 21 for MODIS
!                     landmask(i,j)=0.
!                     xland_rr(i,j)=2.
!                     isltyp(i,j)=14.
! 

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
  real(r_single),allocatable::xland_rr(:,:)
  real(r_single),allocatable::lu_index(:,:)
  integer,allocatable:: ivgtyp(:,:)
  integer,allocatable:: isltyp(:,:)
!
!
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
!
  integer :: i,j,k
  integer :: iii,jjj,num,numh,i4,j4
  real    :: newvalue, newvalueh
!
! =============================================================
!
  fractional_seaice=1
  if ( fractional_seaice == 0 ) then
    xice_threshold = 0.5
    write(*,*) ' do not use fraction sea ice'
  else if ( fractional_seaice == 1 ) then
    xice_threshold = 0.02
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
! calculate snow and graupal precipiation over surface
!tgs 10 mar 2013 - remove liquid precip for snow trimming to work correctly

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

  allocate(surftemp(nlon_regional,nlat_regional))
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
     surftemp(:,:)=tmp4b3d(:,:,nsig_regional)
     deallocate(tmp4b3d)
     call fv3grid%close()
     write(6,*) 'surface temperature =',maxval(surftemp),minval(surftemp)
  endif

  if(fv3_io_layout_y > 1 ) then
     write(flnm1,'(a,I4.4)') 'sfc_data.nc.',id-1
  else
     flnm1='sfc_data.nc.0000'
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

  allocate(snow(nlon_regional,nlat_regional))
  allocate(snowh(nlon_regional,nlat_regional))
  allocate(snowc(nlon_regional,nlat_regional))
  allocate(seaice(nlon_regional,nlat_regional))
  allocate(tskin(nlon_regional,nlat_regional))
  allocate(tsnow(nlon_regional,nlat_regional))

  allocate(landmask(nlon_regional,nlat_regional))
  allocate(xland_rr(nlon_regional,nlat_regional))
  allocate(lu_index(nlon_regional,nlat_regional))
  allocate(ivgtyp(nlon_regional,nlat_regional))
  allocate(isltyp(nlon_regional,nlat_regional))
  
!  rmse_var='TSK'
  call fv3grid%get_var("tsfcl",nlon,nlat,tmp8b2d)
  tskin=tmp8b2d(:,:)
  write(6,*)' max,min skin temp (K)=',maxval(tskin),minval(tskin)
!
!  rmse_var='SOILT1'
  call fv3grid%get_var("tsnow_land",nlon,nlat,tmp8b2d)
  tsnow=tmp8b2d(:,:)
  write(6,*)' max,min snow temp (K)=',maxval(tsnow),minval(tsnow)
!
!  rmse_var='SNOW'
  call fv3grid%get_var("weasdl",nlon,nlat,tmp8b2d)
  snow=tmp8b2d(:,:)
  write(6,*)' max,min SNOW=',maxval(snow),minval(snow)
!
!  rmse_var='SNOWH'
!  call fv3grid%get_var("snwdph",nlon,nlat,tmp8b2d)
!??? not used in this code???  snowdp=tmp8b2d(:,:)
!  write(6,*)' max,min SNOWH=',maxval(snowh),minval(snowh)
!  rmse_var='SNOWH'
  call fv3grid%get_var("snodl",nlon,nlat,tmp8b2d)
  snowh=tmp8b2d(:,:)
  write(6,*)' max,min SNOWH=',maxval(snowh),minval(snowh)
!
!  rmse_var='SNOWC'
  call fv3grid%get_var("sncovr",nlon,nlat,tmp8b2d)
  snowc=tmp8b2d(:,:)
  write(6,*)' max,min SNOWC=',maxval(snowc),minval(snowc)
!
!  rmse_var='SEAICE'
  call fv3grid%get_var("fice",nlon,nlat,tmp8b2d)
  seaice=tmp8b2d(:,:)
  write(6,*)' max,min SEAICE=',maxval(seaice),minval(seaice)
!
  allocate(soiltemp(nlon_regional,nlat_regional,nsig_soil_regional))
  allocate(soiltempRRbk(nlon_regional,nlat_regional,nsig_soil_regional))
!  rmse_var='TSLB'
  call fv3grid%get_var("tslb",nlon,nlat,nsig_soil_regional,tmp8b3d)
  do k=1,nsig_soil_regional
     soiltemp(:,:,k)=tmp8b3d(:,:,k)
     write(6,*)' max,min TSLB=',k, maxval(soiltemp(:,:,k)),minval(soiltemp(:,:,k))
  enddo
!
  deallocate(tmp8b3d)
!
!  rmse_var='XLAND'
  call fv3grid%get_var("slmsk",nlon,nlat,tmp8b2d)
  xland_rr=tmp8b2d
  write(6,*)' max,min XLAND=',maxval(xland_rr),minval(xland_rr)
!
!  rmse_var='LANDMASK'
  call fv3grid%get_var("slmsk",nlon,nlat,tmp8b2d)
  landmask=tmp8b2d
  write(6,*)' max,min LANDMASK=',maxval(landmask),minval(landmask)
!
!  rmse_var='IVGTYP'
  call fv3grid%get_var("vtype",nlon,nlat,tmp8b2d)
  ivgtyp=int(tmp8b2d)
  write(6,*)' max,min IVGTYP=',maxval(ivgtyp),minval(ivgtyp)
!
!  rmse_var='ISLTYP'
  call fv3grid%get_var("stype",nlon,nlat,tmp8b2d)
  isltyp=int(tmp8b2d)
  write(6,*)' max,min ISLTYP=',maxval(isltyp),minval(isltyp)

  deallocate(tmp8b2d)
  call fv3grid%close()

  write(6,*) '================================================='
!
! save the RR background snow in snowRRbk
!

  allocate(snowRRbk(nlon_regional,nlat_regional))
  allocate(snowhRRbk(nlon_regional,nlat_regional))
  allocate(snowcRRbk(nlon_regional,nlat_regional))
  allocate(tskinRRbk(nlon_regional,nlat_regional))
  allocate(tsnowRRbk(nlon_regional,nlat_regional))
  allocate(surftempRRbk(nlon_regional,nlat_regional))

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
  DO J=1,nlat
  DO I=1,nlon
! xland is the RR land/water mask from the geo* file - no effect from sea ice, =1 for land, 0 - water.
    if(landmask(i,j) == 1 ) then  ! on land
    !if(int(xland(i,j)+0.01) == 1 .and. int(seaice(i,j)+0.01) == 0  ) then  ! on land
      if(snowiceRR(i,j) < 1.0e-12 .and. snow(i,j) > 0.0 ) then   ! over forecast snow ?
!tgs may be increase 274K to 276K? Sometimes 100-200mm of snow trimmed with 274.
!tgs 13 April 2012 - change 276K to 280K
!tgs 10 March 2013 - not enough snow trimming in KS - turn temp threshold back
!                    to 276 K. 
!      if(precip(i,j) < 1.0e-12 .and. surftemp(i,j) > 280.0 ) then   ! make sure 
!      if(precip(i,j) < 1.0e-12 .and. surftemp(i,j) > 276.0 ) then   ! make sure 
        if(precip(i,j) < 1.0e-12) then   ! make sure 
          write(6,*) 'trim snow',i,j,snow(i,j),precip(i,j),surftemp(i,j),snowiceRR(i,j) 
          numtrimsnow=numtrimsnow+1
          itr=i
          jtr=j

! save values of snow to be trimmed
          snowtr=snow(i,j)
          snowhtr=snowh(i,j)
          snowctr=snowc(i,j)
          snowtrimsum=snowtrimsum+snow(i,j)
! trim snow
          snow(i,j) = 0.0
          snowh(i,j) = 0.0
          snowc(i,j) = 0.0
        endif
      endif

!tgs snow building
      if(snowiceRR(i,j) > 1.0e-12 .and. snow(i,j) == 0.0 ) then   !  underforecasted snow
        if(surftemp(i,j) < 278.0 ) then   
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
           ist=max(1,i-2)
           iend=min(nlon,i+2)
           jst=max(1,j-2)
           jend=min(nlat,j+2)
           do ii=ist,iend
           do jj=jst,jend
             if(landmask(ii,jj) == 1) then  ! land
             !if(int(xland(ii,jj)+0.01) == 1) then  ! land
               if(ii.eq.itr.and.jj.eq.jtr) then 
! snow trimmed at the neighbor point
                 numnb=100
               endif
               if( numnb== 100) exit

               if(snowRRbk(ii,jj) > 1.) then
                 numnb = numnb + 1
                 snowsum = snowsum + snowRRbk(ii,jj)
                 snowhsum = snowhsum + snowhRRbk(ii,jj)
                 snowcsum = snowcsum + snowcRRbk(ii,jj) 
                 tskinsum = tskinsum + tskinRRbk(ii,jj)
                 tsnowsum = tsnowsum + tsnowRRbk(ii,jj)
                 soilt1sum = soilt1sum + soiltempRRbk(ii,jj,1)
                 soilt2sum = soilt2sum + soiltempRRbk(ii,jj,2)
                 soilt3sum = soilt3sum + soiltempRRbk(ii,jj,3)
                 surftempsum = surftempsum + surftempRRbk(ii,jj)
               endif
             endif
           enddo
             if( numnb == 100) exit
           enddo

! compute averages for all neighbor land points
           if( (numnb.ge.1) .and. (numnb .ne. 100)) then
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
             numusetrim=numusetrim+1
             write(6,*) 'trimmed snow at itr,jtr',itr,jtr,'is used to build snow at point i,j',i,j
             write(6,*) 'snowtr, snowhtr, snowctr, tskin(itr,jtr), tsnow(itr,jtr)', &
                  snowtr, snowhtr, snowctr,tskin(itr,jtr),tsnow(itr,jtr)
             if(snowhtr > 1.e-12) then
!                rhosn=max(76.9,min(500.,snowtr/snowhtr))
!tgs 26jun18 - consistency with the changed limits of snow density in RUC LSM.
! bug fix 12mar2019                rhosn=max(58.8,min(500.,snowav/snowhav))
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
!tgs 22apr15 - this warning is OK if the GFS background snow is getting trimmed (cold-start).
! This warning in the cycled RAP and HRRR indicates a problem.
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

             if(numnb.ge.1) then
               if(snowhav > 1.e-12 .and. snowav > 1.e-12) then
                 write(6,*)'build snow based on neighbor points ',numnb
!                  rhosn=max(76.9,min(500.,snowav/snowhav))
!tgs 26jun18 - consistency with the changed limits of snow density in RUC LSM.
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
!tgs 22apr15 - this warning is OK if the GFS background snow is getting trimmed (cold-start).
! This warning in the cycled RAP and HRRR indicates a problem.
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
             else
               write(6,*) 'set snow to min value'
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
           endif  !  if(numnb == 100) then
           snowbuiltsum=snowbuiltsum+snow(i,j)
           write(6,*) 'BUILD - snow,snowh,snowc,tskin,tsnow,soiltemp1,soiltemp2,soiltemp3', &
              i,j,snow(i,j),snowh(i,j),snowc(i,j),tskin(i,j),tsnow(i,j),soiltemp(i,j,1),soiltemp(i,j,2),soiltemp(i,j,3)       
        endif
      endif
    endif

! limit snow depth not to exceed 50 m
    if((snowh(i,j) >= 0. .and. snowh(i,j) <=50.0) .and. (snow(i,j)  <=20000. .and. snow(i,j)  >=0.) ) then
    elseif(snowh(i,j) < 0. .or. snow(i,j)  < 0.) then
      snowh(i,j)=0.
      snow(i,j) = 0.
    elseif(snowh(i,j) > 50. .or. snow(i,j)  > 20000.) then
      write(6,*) 'Huge snow value i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
      newvalue=0.0
      newvalueh=0.0
      num=0
      numh=0
      do jjj=j-1,j+1
        do iii=i-1,i+1
          !write(6,*) iii,jjj,snowh(iii,jjj),snow(iii,jjj)
          if(iii .ne. i .and. jjj .ne. j) then
            i4=min(max(iii,1),nlon)
            j4=min(max(jjj,1),nlat)
            newvalue=newvalue+snow(i4,j4)
            newvalueh=newvalueh+snowh(i4,j4)
            num=num+1
          endif
        enddo
      enddo
      if(num > 0 .and. newvalue < 100000.0 .and. newvalueh < 200.0) then
        snow(i,j)=newvalue/num
        snowh(i,j)=newvalueh/num
      else
        i4=min(max(i-1,1),nlon)
        j4=min(max(j-1,1),nlat)
        snow(i,j)=snow(i4,j4)
        snowh(i,j)=snowh(i4,j4)
      endif

      write(6,*)'Corrected snow value i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
    else
      write(6,*) '===>Error<===: strange point i,j,snowh(i,j),snow(i,j)',i,j,snowh(i,j),snow(i,j)
      snowh(i,j) = 0.0
      snow(i,j)  = 0.0
      snowc(i,j) = 0.0
    endif
! check consistency of snow variables after snow trim
    if((snow(i,j) <= 0..and.snowh(i,j) > 0.) .or. (snowh(i,j) <=0..and.snow(i,j) > 0.)) then
      write(6,*) 'Inconsistency of snow and snowh AFTER snow trim at i,j,snow,snowh', i,j,snow(i,j),snowh(i,j)
      snow(i,j)  = 0.
      snowh(i,j) = 0.
      snowc(i,j) = 0.
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
!  replace seaice and xland
!
if(1==2) then  ! turn off , use GFS sea ice
  num_seaice2water=0
  num_water2seaice=0
  DO J=1,nlat
  DO I=1,nlon
    if( int(xland(i,j)+0.01) == 0 ) then    ! water
      if(seaice(i,j) >= xice_threshold .and. snowiceRR(i,j) < xice_threshold ) then  ! turn old seaice into water
! For water:
!for MODIS
!            ivgtyp(i,j)=luse(i,j)
!            lu_index(i,j)=luse(i,j)
        ivgtyp(i,j)=17
        lu_index(i,j)=17
        landmask(i,j)=0.
        xland_rr(i,j)=2.
        isltyp(i,j)=14
        seaice(i,j)=snowiceRR(i,j)
        num_seaice2water = num_seaice2water + 1
      elseif(seaice(i,j) < xice_threshold .and. snowiceRR(i,j) >= xice_threshold .and. surftemp(i,j) < 280.) then  ! turn old water into seaice
! for sea ice
!for MODIS
        ivgtyp(i,j)=15
        lu_index(i,j)=15
        landmask(i,j)=1.   
        xland_rr(i,j)=1.     
        isltyp(i,j)=16 
        seaice(i,j)=snowiceRR(i,j)
        num_water2seaice=num_water2seaice+1
       else

      endif
    else
!land - nothing to do here
!switch to MODIS for land
!             ivgtyp(i,j)=luse(i,j)
!             lu_index(i,j)=luse(i,j)
! make sure landmask and xland are consistent for land
!             landmask(i,j)=1.
!             xland_rr(i,j)=1.
    endif
  ENDDO
  ENDDO
  write(6,*) 'SUMMARY on seaice:'
  write(6,*) 'grid point from old seaice into water: ', num_seaice2water
  write(6,*) 'grid point from old water  into seaice: ', num_water2seaice 
endif ! 1==2

!!! Security check for consistency of all land surface parameters on water/ice:
  DO J=1,nlat
  DO I=1,nlon
    if(landmask(i,j) == 0 ) then    ! water
    !if( int(xland(i,j)+0.01) == 0 ) then    ! water
      if(seaice(i,j) < xice_threshold) then
!water
! for MODIS water category is 17
            ivgtyp(i,j)=0 !17
            !lu_index(i,j)=17
            landmask(i,j)=0.
            !xland_rr(i,j)=2.
            isltyp(i,j)=0 !14 ! STASGO water
      else
!ice
!for MODIS ice category is 15
             ivgtyp(i,j)=15
             !lu_index(i,j)=15
             landmask(i,j)=2
             !xland_rr(i,j)=1.
             isltyp(i,j)=16 ! STASGO ice
      endif
    endif  ! water
  ENDDO
  ENDDO
!
!  get rid of snow on water
!
  DO J=1,nlat
  DO I=1,nlon
    if( landmask(i,j) == 0 ) then    ! water
    !if( int(xland(i,j)+0.01) == 0 ) then    ! water
      !do k=1,nsoil
      !   soilmoisture(i,j,k)=1.
      !enddo
      if( seaice(i,j) < 0.001 .and. snow(i,j) > 0.0 ) then  ! snow on water
        snow(i,j) = 0.0
        snowh(i,j) = 0.0
        snowc(i,j) = 0.0
      endif
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
     tmp8b3d(:,:,k)=soiltemp(:,:,k)
     write(6,*)' max,min TSLB=',k, maxval(tmp8b3d(:,:,k)),minval(tmp8b3d(:,:,k))
  enddo
  call fv3grid%replace_var("tslb",nlon,nlat,nsig_soil_regional,tmp8b3d)
  deallocate(tmp8b3d)
! SNCOVR
  tmp8b2d=snowc
  write(6,*)' max,min snowc=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("sncovr",nlon,nlat,tmp8b2d)
! 'SNODL'
  tmp8b2d=snowh
  write(6,*)' max,min snowh=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("snodl",nlon,nlat,tmp8b2d)
! 'SNWDPH'
  tmp8b2d=snowh
  write(6,*)' max,min snowh=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("snwdph",nlon,nlat,tmp8b2d)
! 'WEASDL'
  tmp8b2d=snow
  write(6,*)' max,min snow=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("weasdl",nlon,nlat,tmp8b2d)
! TSFC
  tmp8b2d=tskin
  write(6,*)' max,min tsk=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsfc",nlon,nlat,tmp8b2d)
! TSFCL
  tmp8b2d=tskin
  write(6,*)' max,min tsk=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsfcl",nlon,nlat,tmp8b2d)
! TSNOW_LAND
  tmp8b2d=tsnow
  write(6,*)' max,min soilt1=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("tsnow_land",nlon,nlat,tmp8b2d)
! slmsk
  tmp8b2d=landmask
  write(6,*)' max,min landmask=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("slmsk",nlon,nlat,tmp8b2d)
! STYPE
  tmp8b2d=isltyp
  write(6,*)' max,min isltyp=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("stype",nlon,nlat,tmp8b2d)
! VTYPE
  tmp8b2d=ivgtyp
  write(6,*)' max,min ivgtyp=',maxval(tmp8b2d),minval(tmp8b2d)
  call fv3grid%replace_var("vtype",nlon,nlat,tmp8b2d)
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
  deallocate(xland_rr)
  deallocate(lu_index)
  deallocate(ivgtyp)
  deallocate(isltyp)

  deallocate(soiltemp)
  deallocate(soiltempRRbk)
  deallocate(precip)
  deallocate(surftemp)

end subroutine update_snowice_fv3lam
