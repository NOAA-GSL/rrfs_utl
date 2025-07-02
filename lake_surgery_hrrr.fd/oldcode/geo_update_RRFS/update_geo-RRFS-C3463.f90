
subroutine update_geo_rrfs (nlon, nlat)

  use kinds, only: r_single,i_kind
  use module_ncio, only : ncio
  implicit none


  INCLUDE 'netcdf.inc'

  type(ncio)     :: rrfs

  integer, parameter :: WRF_INTEGER = 106
!
  integer :: nlon, nlat
  real  :: sstRR(nlon,nlat)

! Declare local parameters

  character(len=120) :: flnm1,flnm2
  character(len=19)  :: DateStr1
  integer(i_kind)            :: dh1
 
  integer(i_kind) :: i,j,k, ifound
  integer(i_kind) :: ndim1
  integer(i_kind) :: WrfType
  integer(i_kind), dimension(4)  :: start_index, end_index
  character (len= 4) :: staggering=' N/A'
  character (len= 3) :: ordering
  character (len=31) :: name,name1,name2,name3,name4,name5
  
  character (len=80), dimension(3)  ::  dimnames
  character (len=80) :: SysDepInfo
  
  integer(i_kind) :: l, n, nk, mk, ll, ipb,ipe,jpb,jpe, iter

  integer imax,jmax, nx, ny, nz
  real(4),allocatable :: field2d4b(:,:)
  real(4),allocatable :: field    (:,:)
  !  data imax /1799/, jmax /1059/
    data imax /3950/, jmax /2700/
  !data imax /1820/, jmax /1092/
  
  integer(i_kind) :: ierr, ier, Status, Status_next_time

! rmse stuff
  
  character (len=31) :: rmse_var
  integer(i_kind) iyear,imonth,iday,ihour,iminute,isecond
  integer(i_kind) nlon_regional,nlat_regional,nsig_regional
  real(r_single),allocatable::field2(:,:)
  real(r_single),allocatable::field3(:,:,:)
  real(r_single),allocatable::field3f(:,:,:)
  real(r_single),allocatable::landmask(:,:)
  real(r_single),allocatable::lake_frac(:,:)
  real(r_single),allocatable::hgt(:,:)
  real(r_single),allocatable::lake_d(:,:)
  real(r_single),allocatable::lu_index(:,:)
  real(r_single),allocatable::landusefraction(:,:,:)
  real(r_single),allocatable::landusefraction_21(:,:,:)
  real(r_single),allocatable::vegfrac(:,:,:)
  real(r_single),allocatable::lai(:,:,:)
  real(r_single),allocatable::alb(:,:,:)
  real(r_single),allocatable::snoalb(:,:)
  real(r_single),allocatable::soiltyptop(:,:,:)
  real(r_single),allocatable::soiltypbot(:,:,:)
!
  real(r_single),allocatable::landmask_flnm2(:,:)
  real(r_single),allocatable::hgt_flnm2(:,:)
  real(r_single),allocatable::stcat(:,:)
  real(r_single),allocatable::sbcat(:,:)
  real(r_single),allocatable::soiltem(:,:)
  real(r_single),allocatable::slope(:,:)
  real(r_single),allocatable::lu_index_flnm2(:,:)
  real(r_single),allocatable::landusefraction_flnm2(:,:,:)
  real(r_single),allocatable::vegfrac_flnm2(:,:,:)
  real(r_single),allocatable::lai_flnm2(:,:,:)
  real(r_single),allocatable::alb_flnm2(:,:,:)
  real(r_single),allocatable::soiltyptop_flnm2(:,:,:)
  real(r_single),allocatable::soiltypbot_flnm2(:,:,:)

  integer(i_kind) :: ichange(2000,2000)
  integer(i_kind) :: inear  (2000,2000)

  integer(i_kind) wrf_real

  real(r_single)    :: time, time1, time2, hav
  real(r_single)    :: a, b, iswater
  integer(i_kind), dimension(4)  :: start_index1,  end_index1
  integer(i_kind) :: nlakes, ilake, icen, jcen, nchange

  integer(i_kind) :: ilakecen(8), jlakecen(8)

! *** i/j dimensions for RRFS C3463 domain - 3950 x 2700
   data ilakecen /2614,2674,2776,2803,2845,2908,2741,1945/ 
   data jlakecen /1192,1057,1123,1168,1045,1144,1167, 822/
!  1/2/3/4/5/6/7 = Superior, Michigan,Huron, Georgian Bay, Eric, Ontario, N. Channel
!         8 = Bonneville

! - i/j for RRFS CONUS C3359 domain - 1820 x 1092
!  data ilakecen /1159,1193,1308,1343,1354,1433,1281, 462/ 
! data jlakecen / 905, 784, 811, 844, 713, 794, 855, 681/

  nlakes = 8
  wrf_real=104

   flnm1 = 'orofile.nc'
! read in rrfs latlon

!   call rrfs%open('gfsice.sfc_data.nc',"r",200)
    call rrfs%open('orofile.nc',"r",200)

    call rrfs%get_dim("lon",nx)

    call rrfs%get_dim("lat",ny)

    nz=1
!    call rrfs%get_dim("zaxis_1",nz)

    print *,  'nx_rrfs,ny_rrfs=',nx,ny,nz

    allocate(field2d4b(nx,ny))

    call rrfs%get_var('orog_raw',nx,ny,field2d4b)
    write(*,*) maxval(field2d4b),minval(field2d4b)


  rmse_var='LANDUSEF'

    call rrfs%get_var('slmsk',nx,ny,field2d4b)
    write(*,*) maxval(field2d4b),minval(field2d4b)
  landmask=field2d4b

    call rrfs%get_var('lake_frac',nx,ny,field2d4b)
    write(*,*) maxval(field2d4b),minval(field2d4b)
  lake_frac = field2d4b

    call rrfs%get_var('lake_depth',nx,ny,field2d4b)
    write(*,*) maxval(field2d4b),minval(field2d4b)
  lake_d=field2d4b



    deallocate(field2d4b)
    call rrfs%close()

      print *,' We got this far'
!

! loop around Great Lakes with Rhoomba to remove lake bathymetry for these large lakes
!  ilakecen(1) = 580


!   do ilake = 1,5
    do ilake = 1,7
     icen = ilakecen(ilake)
     jcen = jlakecen(ilake)
     inear (icen  ,jcen) = 1
!    if (lake_d(icen,jcen).lt.10.000001   &
!     .and. lake_d(icen,jcen).gt. 9.99999) write (6,*) 'Lake #',ilake,' center not in lake'
!    if (lake_d(icen,jcen).eq.10.) write (6,*) 'Lake #',ilake,' center not in lake'
     if (landmask(icen,jcen).eq.1.) write (6,*) 'Lake #',ilake,' center not in lake'
   nchange = 0

    do i=icen-1,icen+1
    do j=jcen-1,jcen+1
        write(6,*)'Lake',ilake,'  lake_d(i,j)',i,j,lake_d(i,j)
    enddo
    enddo

   do iter=0,120
!    do j=max(1,jcen-itea)r,min(jmax,jcen+iter)
!    do j=jcen-iter,min(jmax,jcen+iter)
!    do i=icen-iter,min(imax,icen+iter)
     do j=jcen-iter,jcen+iter
     do i=icen-iter,icen+iter
       if (ichange(i,j) >ilake+1) then
          write (6,*) "Too close, i/j = ", i,j,ichange(i,j),lake_d(i,j)
          go to 222
! --- skip to next grid point

       end if
!       if ((lake_d(i,j).gt.10.00001 .or. lake_d(i,j).lt. 9.99999) .and. inear(i,j).eq.1) then
!      if (lake_d(i,j).ne.10. .and. inear(i,j).eq.1) then
!      if (lake_d(i,j).ne.10. .and. inear(i,j).eq.1 .and. lu_index(i,j).eq.21) then
       if (landmask(i,j).ne.1. .and. inear(i,j).eq.1) then
         nchange = nchange + 1
!        write(6,*) 'Change lake_d,landmask',i,j,lake_d(i,j),landmask(i,j)
         inear (i+1,j ) = 1
         inear (i+1,j+1 ) = 1
         inear (i-1,j ) = 1
         inear (i-1,j-1 ) = 1
         inear (i  ,j+1) = 1
         inear (i+1,j+1) = 1
         inear (i-1,j+1) = 1
         inear (i  ,j-1) = 1
         inear (i-1,j-1) = 1
         inear (i+1,j-1) = 1
         
! inear(i,j)=1  means that this point is a candidate for another lake point
!  It is *adjacent* to another lake point.  It means that these points are *eligible*
!     to be evaluated in the next do-loop interation.
!     In this way, the Rhoomba search continues.
!  NOTE: in this version, the search continues over the full 100 dx in i/j directions
!     even if no more lake-water points are discovered.  i.e., it is not dependent (yet)
!     on *contiguous* lake points.
!  Possible improvement - stop searching when full lake edge (coast) has been discovered.
!    This can be done when inear=0 for all points in a perimeter "iter" points from the initial point.

!  These are the new values for the Great Lakes points being removed from the lu_index=21 points
         lake_d(i,j) = 10.
! ***** Key line below
         landmask(i,j) = 1.
         lake_d  (i,j) = 0.
         lake_frac (i,j) = 0.
         ichange(i,j) = ilake
!        hgt(i,j) = 1.
       else 
!        write(6,*) 'Shore encountered',i,j
         ichange(i,j) = ilake+1
       end if
222    continue
     end do     
     end do

   end do
   write(6,*) 'Lake #',ilake, '  No. of points changed = ',nchange
   nchange = 0  
   end do  ! ilake

!   go to 876

!  Search for a box for the Bonneville Salt Flats
    ilake = 8
     icen = ilakecen(ilake)
     jcen = jlakecen(ilake)
     print *,' '
     print *,' Remove Bonneville Salt Flats'
   do j=jcen-7,jcen+9
   do i=icen-8,icen+7
       if (landmask(i,j).ne.1.                      ) then
         nchange = nchange + 1
         write(6,*) 'Change lake_d,landmask',i,j,lake_d(i,j),landmask(i,j)
         landmask(i,j) = 1.
         lake_frac(i,j) = 0.
         lake_d   (i,j) = 0.
       end if
   end do
   end do

   write(6,*) 'Bonneville-  No. of points changed = ',nchange

876   continue

!  test writing out revised orofile
    print *, ' open rev.orofile.nc'
!   call rrfs%open('rev-orofile.nc',"w",200)
    call rrfs%open('rev-orofile.nc',"w",400)
!   call rrfs%write_var('slmsk',nx,ny,landmask)
!   call rrfs%replace_var('slmsk',nx,ny,landmask)
    call rrfs%replace_var('lake_depth',nx,ny,lake_d  )
    call rrfs%replace_var('lake_frac',nx,ny,lake_frac)
    call rrfs%close()

   stop 999


if(1==1) then
  write(6,*) '================================================='
  rmse_var='ALBEDO12M'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
!  deallocate(field3f)
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  allocate(alb(nlon_regional,nlat_regional,nsig_regional))
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )

  alb=field3f

    do i=74,80
     do j=269,275
!       write(6,*)'before update alb(i,j)', i,j,alb(i,j,1:nsig_regional)
     enddo
    enddo
endif !(1==2)
!
  call ext_ncd_ioclose(dh1, Status)

! read in fields from the second geo* file that
! doeas not have a lake in the
! middle of the Big Island in Hawaii
! bad points: 76,275 and 77,275
! call ext_ncd_open_for_read( trim(flnm2), 0, 0, "", dh1, Status)
  if ( Status /= 0 )then
     write(6,*)'CONVERT_NETCDF_MASS:  problem with flnm2 = ',&
          trim(flnm2),', Status = ', Status
     stop 74 
  endif

!-------------  get date info

  call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*)' Skin temp data from background file at time:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond
!
!-------------  get grid info

  rmse_var='LANDUSEF'

  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )                !DEDE

  write(6,*)' dh1  = ',dh1         !DEDE
  write(6,*)'rmse_var = ',trim(rmse_var)
  write(6,*)'ndim1 = ',ndim1
  write(6,*)'ordering = ',trim(ordering)
  write(6,*)'staggering = ',trim(staggering)
  write(6,*)'start_index = ',start_index
  write(6,*)'end_index = ',end_index
  write(6,*)'WrfType = ',WrfType
  write(6,*)'ierr  = ',ierr   !DEDE

  nlon_regional=end_index(1)
  nlat_regional=end_index(2)
  nsig_regional=end_index(3)
  write(6,*)' nlon,lat,sig_regional=',nlon_regional,nlat_regional,nsig_regional
  allocate(landmask_flnm2(nlon_regional,nlat_regional))
  allocate(lu_index_flnm2(nlon_regional,nlat_regional))
  allocate(landusefraction_flnm2(nlon_regional,nlat_regional,nsig_regional))
  allocate(hgt_flnm2(nlon_regional,nlat_regional))

if(1==2) then
  write(6,*) '================================================='
  rmse_var='LU_INDEX'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  write(6,*)' max,min LU_INDEX=',maxval(field2),minval(field2)
  lu_index_flnm2=field2
   do i=76,78
     do j=275,275
        write(6,*)'LU_INDEX_flnm2(i,j)',i,j,LU_INDEX_flnm2(i,j)
     enddo
    enddo
endif ! 1==2
!
  write(6,*) '================================================='
  rmse_var='LANDMASK'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  end_index(3) =1
  write(6,*)' end_index=',end_index
  call ext_ncd_read_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
  landmask_flnm2=field2
   do i=76,77
     do j=275,275
!        write(6,*)'LANDMASK_flnm2(i,j)',i,j,landmask_flnm2(i,j)
     enddo
    enddo
!
  call ext_ncd_ioclose(dh1, Status)
!
if(1==1) then

  DO J=1,nlat
  DO I=1,nlon
      iswater = 21.   ! MODIS 30" with lakes
! reduce the area of the Great Salt Lake based on flnm2 file
    if(landmask_flnm2(i,j)==1. .and. landmask(i,j)==0. ) then
 print *,'Correct Salt Lake i,j',i,j
       ifound=0
        DO LL=1,16
          JPE = MIN (nlat-1, J+LL)
          JPB = MAX (1 , J-LL)
          IPE = MIN (nlon-3, I+LL)
          IPB = MAX (1 , I-LL)
!
          DO NK=IPB,IPE
           DO MK=JPB,JPE
             IF ( landmask(nk,mk)==1. .and. ifound ==0 ) THEN
               alb(i,j,:)=alb (nk,mk,:)
               lai(i,j,:)=lai(nk,mk,:)
               vegfrac(i,j,:)=vegfrac(nk,mk,:)
               landusefraction(i,j,:)=landusefraction(nk,mk,:)
               soiltyptop(i,j,:)=soiltyptop(nk,mk,:)
               soiltypbot(i,j,:)=soiltypbot(nk,mk,:)
               lu_index(i,j)=lu_index(nk,mk)
               snoalb(i,j)=snoalb(nk,mk)
               stcat(i,j)=stcat(nk,mk)
               sbcat(i,j)=sbcat(nk,mk)
               soiltem(i,j)=soiltem(nk,mk)
               slope(i,j)=slope(nk,mk)
               landmask(i,j)= 1.
               lake_d(i,j)= 10.

              ifound=1
            ENDIF
           ENDDO  ! MK
          ENDDO  ! NK
        ENDDO  ! LL
     if(ifound==0) print *,'Land point is not found',i,j
    endif

  ENDDO
  ENDDO
endif ! 1==2

if(1==2) then
! replace topo with the smoothed topography
! 1 - South America
  DO J=1,185
  DO I=775,954
     hgt(i,j)=hgt_flnm2(i,j)
  ENDDO
  ENDDO
! 2 - Greenland and Alps
  DO J=550,834
  DO I=500,954
     hgt(i,j)=hgt_flnm2(i,j)
  ENDDO
  ENDDO
! 3 - Hawaii
  DO J=250,320
  DO I=1,100
     hgt(i,j)=hgt_flnm2(i,j)
  ENDDO
  ENDDO
endif ! 1==2

!
!
!           update mass core netcdf file with new LU_INDEX
!
  write(6,*) ' ============================= '
  write(6,*) ' update GREENFRAC in geo file '
!  write(6,*) ' update LU_INDEX and LANDUSEF in geo file '
  write(6,*) ' ============================= '
  flnm1='geo_em.d01-mod2.nc'
  call ext_ncd_open_for_update( trim(flnm1), 0, 0, "", dh1, Status)
  if ( Status /= 0 )then
     write(6,*)'UPDATE_NETCDF_MASS:  problem with flnm1 = ',&
          trim(flnm1),', Status = ', Status
     stop 75
  endif
     
!-------------  get date info

  call ext_ncd_get_next_time(dh1, DateStr1, Status_next_time)
  read(DateStr1,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
  write(6,*) ' Update LU_INDEX in background at time:'
  write(6,*)' iy,m,d,h,m,s=',iyear,imonth,iday,ihour,iminute,isecond
!
!-------------  get grid info
  rmse_var='SOILCBOT'
  call ext_ncd_get_var_info (dh1,rmse_var,ndim1,ordering,staggering, &
                               start_index,end_index1, WrfType, ierr    )
  if( (nlon_regional .ne. end_index1(1)) .or.    &
      (nlat_regional .ne. end_index1(2)) ) then
      write(6,*) ' Dimensions do not match!!!'
      write(6,*)' nlon,lat=',nlon_regional,nlat_regional
      stop 123
  endif

! ********************* write out modified fields

if(1==1) then
  write(6,*) '================================================='
  field2=lu_index
!  write(6,*)' max,min LU_INDEX  =',maxval(field2),minval(field2)
  rmse_var='LU_INDEX'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )

  write(6,*) '================================================='
  field2=landmask
!  write(6,*)' max,min landmask  =',maxval(field2),minval(field2)
  rmse_var='LANDMASK'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=lake_d
  rmse_var='LAKE_DEPTH'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=hgt
  rmse_var='HGT_M'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=snoalb
  rmse_var='SNOALB'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=stcat
  rmse_var='SCT_DOM'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=sbcat
  rmse_var='SCB_DOM'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=soiltem
  rmse_var='SOILTEMP'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  field2=slope
  rmse_var='SLOPECAT'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='

if(1==2) then
  field2=hgt
  write(6,*)' max,min HGT_M  =',maxval(field2),minval(field2)
  rmse_var='HGT_M'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  end_index1(3) =1
  write(6,*)' end_index1=',end_index1
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field2,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
endif ! 1==2

  write(6,*) '================================================='
  rmse_var='SOILCTOP'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  deallocate(field3f)
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  field3f=soiltyptop
  do k=1,nsig_regional
!    write(6,*)' max,min SOILCTOP=',k, maxval(field3f(:,:,k)),minval(field3f(:,:,k))
  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  write(6,*) '================================================='
  rmse_var='SOILCBOT'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  deallocate(field3f)
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  field3f=soiltypbot
  do k=1,nsig_regional
!    write(6,*)' max,min SOILCBOT=',k, maxval(field3f(:,:,k)),minval(field3f(:,:,k))
  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )

  write(6,*) '================================================='
  rmse_var='ALBEDO12M'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  deallocate(field3f)
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  field3f=alb
  do k=1,nsig_regional
!    write(6,*)' max,min ALBEDO12M=',k, maxval(field3f(:,:,k)),minval(field3f(:,:,k))
  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  deallocate(field3f)
endif ! 1==2
  write(6,*) '================================================='
  rmse_var='GREENFRAC'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  field3f=vegfrac
  do k=1,nsig_regional
!    write(6,*)' max,min GREENFRAC=',k, maxval(field3f(:,:,k)),minval(field3f(:,:,k))
  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
  deallocate(field3f)
  write(6,*) '================================================='
if(1==1) then
  rmse_var='LAI12M'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index1, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index1=',end_index1
!  deallocate(field3f)
  nsig_regional=end_index1(3)
  allocate(field3f(nlon_regional,nlat_regional,nsig_regional))
  field3f=lai
  do k=1,nsig_regional
!    write(6,*)' max,min LAI12M=',k, maxval(field3f(:,:,k)),minval(field3f(:,:,k))
  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3f,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index1,               & !dom
       start_index,end_index1,               & !mem
       start_index,end_index1,               & !pat
       ierr                                 )
!
  write(6,*) '================================================='
  rmse_var='LANDUSEF'
  call ext_ncd_get_var_info (dh1,trim(rmse_var),ndim1,ordering,staggering, &
       start_index,end_index, WrfType, ierr    )
  write(6,*)' rmse_var=',trim(rmse_var)
  write(6,*)' ordering=',ordering
  write(6,*)' WrfType,WRF_REAL=',WrfType,WRF_REAL
  write(6,*)' ndim1=',ndim1
  write(6,*)' staggering=',staggering
  write(6,*)' start_index=',start_index
  write(6,*)' end_index=',end_index
  deallocate(field3)
  nsig_regional=end_index(3)
  allocate(field3(nlon_regional,nlat_regional,nsig_regional))
  do k=1,nsig_regional
  field3(:,:,k)=landusefraction(:,:,k)
  enddo
!  do k=1,nsig_regional
!    write(6,*)' max,min LANDUSEF=',k, maxval(field3(:,:,k)),minval(field3(:,:,k))
!  enddo
  call ext_ncd_write_field(dh1,DateStr1,TRIM(rmse_var),              &
       field3,WRF_REAL,0,0,0,ordering,           &
       staggering, dimnames ,               &
       start_index,end_index,               & !dom
       start_index,end_index,               & !mem
       start_index,end_index,               & !pat
       ierr                                 )
endif ! 1==2

  call ext_ncd_ioclose(dh1, Status)

end subroutine update_geo_rrfs

SUBROUTINE wrf_debug( level , str )
! USE module_wrf_error
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

