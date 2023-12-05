program process_NSSL_mosaic
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2007-12-17
!
! ABSTRACT: 
!     This routine read in NSSL reflectiivty mosaic fiels and 
!     interpolate them into GSI mass grid
!
!     tversion=8  : NSSL 8 tiles netcdf
!     tversion=81 : NCEP 8 tiles binary
!     tversion=4  : NSSL 4 tiles binary
!     tversion=1  : NSSL 1 tile grib2
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  mosaic_files
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use mpi
  use kinds, only: r_kind,i_kind
  use module_read_NSSL_refmosaic, only: read_nsslref
  use module_ncio, only : ncio

  implicit none
!
  INCLUDE 'netcdf.inc'
!
  type(ncio) :: rrfs
  type(read_nsslref) :: readref 

!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
  character*256 output_file
!
!  grid
  integer(i_kind) :: nlon,nlat
  real,allocatable:: xlon(:,:)    !
  real,allocatable:: ylat(:,:)    !
  REAL, allocatable :: ref3d(:,:,:)   ! 3D reflectivity
  REAL, allocatable :: ref0(:,:,:)   ! 3D reflectivity
  REAL(r_kind), allocatable :: ref3d_column(:,:)   ! 3D reflectivity in column
  CHARACTER*80   geofile
!
!  namelist files
!
  integer      ::  tversion
  integer      ::  fv3_io_layout_y
  character*10 :: analysis_time
  CHARACTER*180   dataPath
  namelist/setup/ tversion,analysis_time,dataPath,fv3_io_layout_y
  integer(i_kind)  ::  idate
!
!  namelist and other variables for netcdf output
!
!  output_netcdf             logical controlling whether netcdf output file should be created
!  max_height                maximum height (m MSL) for data to be retained
!  use_clear_air_type        logical controlling whether to output clear-air (non-precipitation) reflectivity obs
!  precip_dbz_thresh         threshold (dBZ) for minimum reflectivity that is considered precipitation
!  clear_air_dbz_thresh      threshold (dBZ) for maximum reflectivity that is considered clear air
!  clear_air_dbz_value       value (dBZ) assigned to clear-air reflectivity obs
!  precip_dbz_horiz_skip     horizontal thinning factor for reflectivity data in precipitation
!  precip_dbz_vert_skip      vertical thinning factor for reflectivity data in precipitation
!  clear_air_dbz_horiz_skip  horizontal thinning factor for clear air reflectivity data
!  clear_air_dbz_vert_skip   vertical thinning factor for clear air reflectivity data
!
  logical :: output_netcdf = .false.
  real :: max_height = 20000.0
  logical :: use_clear_air_type = .false.
  real :: precip_dbz_thresh = 15.0
  real :: clear_air_dbz_thresh = 0.0
  real :: clear_air_dbz_value = 0.0
  integer :: precip_dbz_horiz_skip = 0
  integer :: precip_dbz_vert_skip = 0
  integer :: clear_air_dbz_horiz_skip = 0
  integer :: clear_air_dbz_vert_skip = 0
  namelist/setup_netcdf/ output_netcdf, max_height,                       &
                         use_clear_air_type, precip_dbz_thresh,           &
                         clear_air_dbz_thresh, clear_air_dbz_value,       &
                         precip_dbz_horiz_skip, precip_dbz_vert_skip,     &
                         clear_air_dbz_horiz_skip, clear_air_dbz_vert_skip
  logical, allocatable :: precip_ob(:,:,:)
  logical, allocatable :: clear_air_ob(:,:,:)
  integer, parameter :: maxMosaiclvl=33
  real :: height_real(maxMosaiclvl)
  integer :: levelheight(maxMosaiclvl)
  data levelheight /500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500,         &
                    2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, &
                    7500, 8000, 8500, 9000, 10000, 11000, 12000, 13000, 14000,  &
                    15000, 16000, 17000, 18000, 19000/
  integer num_precip_obs, num_clear_air_obs, num_obs
  logical :: fileexist
!
!  ** misc
!      
  character*80 outfile
  character*256 outfile_netcdf
  integer i,ii,j,jj,k,kk
  integer :: id
  INTEGER(i_kind)  ::  maxlvl
  INTEGER(i_kind)  ::  numlvl,numref
  integer :: maxcores

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror) 
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) write(*,*) mype, 'deal with mosaic'

  fv3_io_layout_y=1
  datapath="./"
  open(15, file='namelist.mosaic')
    read(15,setup)
  close(15)

  inquire(file='namelist.mosaic_netcdf', exist=fileexist)
  if(fileexist) then
    if(mype==0) write(*,*) 'reading namelist.mosaic_netcdf'
    open(15, file='namelist.mosaic_netcdf')
    read(15, setup_netcdf)
    close(15)
  endif

  if(mype==0) then
    write(6,*)
    write(6,*) 'tversion = ', tversion
    write(6,*) 'analysis_time = ', analysis_time
    write(6,*) 'dataPath = ', dataPath
    write(6,*) 'fv3_io_layout_y = ', fv3_io_layout_y
    write(6,*) 'output_netcdf = ', output_netcdf
    write(6,*) 'max_height = ', max_height
    write(6,*) 'use_clear_air_type = ', use_clear_air_type
    write(6,*) 'precip_dbz_thresh = ', precip_dbz_thresh
    write(6,*) 'clear_air_dbz_thresh = ', clear_air_dbz_thresh
    write(6,*) 'clear_air_dbz_value = ', clear_air_dbz_value
    write(6,*) 'precip_dbz_horiz_skip = ', precip_dbz_horiz_skip
    write(6,*) 'precip_dbz_vert_skip = ', precip_dbz_vert_skip
    write(6,*) 'clear_air_dbz_horiz_skip = ', clear_air_dbz_horiz_skip
    write(6,*) 'clear_air_dbz_vert_skip = ', clear_air_dbz_vert_skip
    write(6,*)
  endif

!
!  safty check for cores used in this run
!
  read(analysis_time,'(I10)') idate
  if(mype==0) write(6,*) 'cycle time is :', idate

  if( tversion == 8 .or. tversion == 14) then
     maxcores=8
  elseif( tversion == 81 ) then
     maxcores=8
  elseif( tversion == 4 ) then
     maxcores=4
  elseif( tversion == 1 ) then
     maxcores=33
  else
     write(*,*) 'unknow tversion !'
     stop 1234
  endif

  if(mype==0) write(6,*) 'total cores for this run is ',npe
  if(npe < maxcores) then
     write(6,*) 'ERROR, this run must use ',maxcores,' or more cores !!!'
     call MPI_FINALIZE(ierror)
     stop 1234
  endif
!
! read NSSL mosaic 
!
  mypeLocal=mype+1
  call readref%init(tversion,mypeLocal,datapath)
!
! deal with certain tile
!
  call readref%readtile(mypeLocal)
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
  maxlvl=readref%maxlvl
!
!  
  do id=0,fv3_io_layout_y-1
!
! get domain dimension
!
     if(fv3_io_layout_y==1) then
        write(geofile,'(a,a)') './', 'fv3sar_grid_spec.nc'
     else
        write(geofile,'(a,a,I4.4)') './', 'fv3sar_grid_spec.nc.',id
     endif
     call rrfs%open(trim(geofile),"r",0)
     call rrfs%get_dim("grid_xt",nlon)
     call rrfs%get_dim("grid_yt",nlat)
     allocate(xlon(nlon,nlat))
     allocate(ylat(nlon,nlat))
     call rrfs%get_var("grid_lont",nlon,nlat,xlon)
     call rrfs%get_var("grid_latt",nlon,nlat,ylat)
     if(mype==0) then
       write(*,*) 'FV3LAM grid from ',trim(geofile)
       write(*,*) 'nx_rrfs,ny_rrfs=',nlon,nlat
       write(*,*) 'max, min lon=', maxval(xlon),minval(xlon)
       write(*,*) 'max, min lat=', maxval(ylat),minval(ylat)
     endif
     call rrfs%close()

     allocate(ref3d(nlon,nlat,maxlvl))
     ref3d=-999.0
     call mosaic2grid(readref,nlon,nlat,maxlvl,xlon,ylat,ref3d)
     call mpi_barrier(MPI_COMM_WORLD,ierror)
!
!  collect data from all processes to root (0)
!
     if(mype==0) then
        allocate( ref0(nlon,nlat,maxlvl) )
     endif
     call MPI_REDUCE(ref3d, ref0, nlon*nlat*maxlvl, MPI_REAL, MPI_MAX, 0, &
                     MPI_COMM_WORLD, ierror)
     deallocate(ref3d)
!
     if(mype==0) then
        if(fv3_io_layout_y==1) then
           write(outfile,'(a,a)') './', 'RefInGSI3D.dat'
           write(outfile_netcdf,'(a,a)') './', 'Gridded_ref.nc'
        else
           write(outfile,'(a,a,I4.4)') './', 'RefInGSI3D.dat.',id
           write(outfile_netcdf,'(a,a,I4.4,a)') './', 'Gridded_ref.',id,'.nc'
        endif
        OPEN(10,file=trim(outfile),form='unformatted')
           write(10) maxlvl,nlon,nlat
           write(10) ref0
        close(10)
        DO k=1,maxlvl
           write(*,*) k,maxval(ref0(:,:,k)),minval(ref0(:,:,k))
        ENDDO

! turn this part off to speed up the process for RRFS.
       if(1==2) then
!
        allocate(ref3d_column(maxlvl+2,nlon*nlat))
        ref3d_column=-999.0
        numref=0
        DO j=1,nlat
        DO i=1,nlon
          numlvl=0
          DO k=1,maxlvl
            if(abs(ref0(i,j,k)) < 888.0 ) numlvl=numlvl+1
          ENDDO
          if(numlvl > 0 ) then
            numref=numref+1
            ref3d_column(1,numref)=float(i)
            ref3d_column(2,numref)=float(j)
            DO k=1,maxlvl
               ref3d_column(2+k,numref)=ref0(i,j,k)
            ENDDO
          endif
        ENDDO
        ENDDO

        write(*,*) 'Dump out results', numref, 'out of', nlon*nlat
        OPEN(10,file='./'//'RefInGSI.dat',form='unformatted')
         write(10) maxlvl,nlon,nlat,numref,1,2
         write(10) ((ref3d_column(k,i),k=1,maxlvl+2),i=1,numref)
        close(10)
  
        write(*,*) 'Start write_bufr_nsslref'
        call write_bufr_nsslref(maxlvl,nlon,nlat,numref,ref3d_column,idate)
        deallocate(ref3d_column)
       endif

       if ( output_netcdf .and. (maxlvl.eq.maxMosaiclvl) ) then

         allocate( precip_ob(nlon,nlat,maxlvl) )
         allocate( clear_air_ob(nlon,nlat,maxlvl) )

         ! Don't produce any netcdf radar observations along the lateral boundaries
         ref0(1,:,:) = -999.0
         ref0(nlon,:,:) = -999.0
         ref0(:,1,:) = -999.0
         ref0(:,nlat,:) = -999.0

         ! Identify precip and clear-air reflectivity observations
         precip_ob(:,:,:) = .false.
         clear_air_ob(:,:,:) = .false.
         num_precip_obs = 0
         num_clear_air_obs = 0
         do j=2,nlat-1
           do i=2,nlon-1
             do k=1,maxlvl
               if ( (levelheight(k) .le. max_height) .and. (ref0(i,j,k) .ge. precip_dbz_thresh) ) then
                 precip_ob(i,j,k) = .true.
                 num_precip_obs = num_precip_obs + 1
               else if ( use_clear_air_type .and. (levelheight(k) .le. max_height) .and. &
                         (ref0(i,j,k) .gt. -900.0) .and. (ref0(i,j,k) .le. clear_air_dbz_thresh) ) then
                 clear_air_ob(i,j,k) = .true.
                 ref0(i,j,k) = clear_air_dbz_value
                 num_clear_air_obs = num_clear_air_obs + 1
               endif
             enddo
           enddo
         enddo
         write(*,*) 'number of precip obs found, before thinning, = ', num_precip_obs
         write(*,*) 'number of clear air obs found, before thinning, = ', num_clear_air_obs

         ! Thin precip reflectivity observations
         if (precip_dbz_vert_skip .gt. 0) then
           do k=1,maxlvl
             if (mod(k-1, precip_dbz_vert_skip+1) .ne. 0) then
               write(*,*) 'Thinning:  removing precip obs at level ', k
               precip_ob(:,:,k) = .false.
             endif
           enddo
         endif
         if (precip_dbz_horiz_skip .gt. 0) then
           write(*,*) 'Horizontal thinning of precip obs'
           do j=2,nlat-1
             do i=2,nlon-1
               do k=1,maxlvl
                 if (precip_ob(i,j,k)) then
                   do jj=max(2, j-precip_dbz_horiz_skip), min(nlat-1, j+precip_dbz_horiz_skip)
                     do ii=max(2, i-precip_dbz_horiz_skip), min(nlon-1, i+precip_dbz_horiz_skip)
                       precip_ob(ii,jj,k) = .false.
                     enddo
                   enddo
                   precip_ob(i,j,k) = .true.
                 endif
               enddo
             enddo
           enddo
         endif

         ! Thin clear-air reflectivity observations
         if (use_clear_air_type .and. (clear_air_dbz_vert_skip .gt. 0) ) then
           do k=1,maxlvl
             if (mod(k-1, clear_air_dbz_vert_skip+1) .ne. 0) then
               write(*,*) 'Thinning:  removing clear air obs at level ', k
               clear_air_ob(:,:,k) = .false.
             endif
           enddo
         endif
         if (use_clear_air_type .and. (clear_air_dbz_vert_skip .lt. 0) ) then
           write(*,*) 'Thinning:  removing clear air obs at all but two levels'
           clear_air_ob(:, :, 1:12) = .false.
           clear_air_ob(:, :, 14:21) = .false.
           clear_air_ob(:, :, 23:maxlvl) = .false.
         endif
         if (use_clear_air_type .and. (clear_air_dbz_horiz_skip .gt. 0) ) then
           do j=2,nlat-1
             do i=2,nlon-1
               do k=1,maxlvl
                 if (clear_air_ob(i,j,k)) then
                   do jj=max(2, j-clear_air_dbz_horiz_skip), min(nlat-1, j+clear_air_dbz_horiz_skip)
                     do ii=max(2, i-clear_air_dbz_horiz_skip), min(nlon-1, i+clear_air_dbz_horiz_skip)
                       clear_air_ob(ii,jj,k) = .false.
                     enddo
                   enddo
                   clear_air_ob(i,j,k) = .true.
                 endif
               enddo
             enddo
           enddo
         endif

         ! Count number of valid obs
         num_obs = 0
         do j=2,nlat-1
           do i=2,nlon-1
             do k=1,maxlvl
               if ( precip_ob(i,j,k) .or. clear_air_ob(i,j,k) ) then
                 num_obs = num_obs + 1
               endif
             enddo
           enddo
         enddo
         write(*,*) 'num_obs = ', num_obs

         ! Write obs to netcdf file
         do k=1,maxlvl
           height_real(k) = levelheight(k)
           do j=2,nlat-1
             do i=2,nlon-1
               if ( .not. precip_ob(i,j,k) .and. .not. clear_air_ob(i,j,k) ) then
                 ref0(i,j,k) = -999.0
               endif
             enddo
           enddo
         enddo
         call write_netcdf_nsslref( outfile_netcdf,maxlvl,nlon,nlat,ref0,idate,xlon,ylat,height_real )

         write(*,*) 'Finish netcdf output'

         deallocate(precip_ob)
         deallocate(clear_air_ob)

       else if (output_netcdf) then

         write(*,*) 'unknown vertical levels'
         write(*,*) 'maxlvl = ', maxlvl
         write(*,*) 'maxMosaiclvl = ', maxMosaiclvl
         write(*,*) 'no netcdf output'

       endif ! output_netcdf

       deallocate(ref0)

     endif ! mype==0

     deallocate(xlon)
     deallocate(ylat)

  enddo ! id

  call readref%close()

  if(mype==0)  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

  call MPI_FINALIZE(ierror)
!
end program process_NSSL_mosaic

subroutine mosaic2grid(readref,nlon,nlat,maxlvl,xlon,ylat,ref3d)
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-01-20
!
! ABSTRACT: 
!     This routine interpolate NSSL reflectiivty mosaic fields to grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  mosaic_files
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use kinds, only: r_kind,i_kind
  use module_read_NSSL_refmosaic, only: read_nsslref

  implicit none
!
  type(read_nsslref),intent(in) :: readref 
  integer,intent(in) :: nlon,nlat,maxlvl
  real,intent(in)    :: xlon(nlon,nlat)    !
  real,intent(in)    :: ylat(nlon,nlat)    !
  real,intent(inout) :: ref3d(nlon,nlat,maxlvl)   ! 3D reflectivity
!
!  For reflectiivty mosaic
!
  REAL   :: lonMin,latMin,lonMax,latMax
  REAL*8 :: dlon,dlat
!
!  ** misc
!      
  integer i,j,k,kk
  integer :: tversion

  REAL ::  rlat,rlon
  INTEGER  :: ip,jp,ipp1,jpp1
  REAL ::  rip,rjp
  REAL ::  dip,djp
  REAL ::  w1,w2,w3,w4
  REAL ::  ref1,ref2,ref3,ref4,refl_ltng
  real, allocatable :: mscValue(:,:)

!**********************************************************************
!
!
  tversion=readref%tversion
  dlon=readref%dlon
  dlat=readref%dlat
  latMax=readref%latMax
  lonMax=readref%lonMax
  lonMin=readref%lonMin
  latMin=readref%latMin
!
  if(readref%if_fileexist) then
!
     allocate(mscValue(readref%mscNlon,readref%mscNlat))
     do k=1, readref%mscNlev
          if(tversion == 1) then
             kk=readref%ilevel
             mscValue(:,:) = readref%mscValue3d(:,:,1)
          else
             kk=k
             mscValue(:,:) = readref%mscValue3d(:,:,k)
          endif
!
          DO j=1,nlat
          DO i=1,nlon
             rlat=ylat(i,j)
             rlon=xlon(i,j)

             if(tversion == 14 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(rlat-latMin)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             elseif(tversion == 8 .or. tversion == 81) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(latMax-rlat)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp 
             elseif(tversion == 4 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(rlat-latMin)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             elseif(tversion == 1 ) then
               if(rlon<0.0) rlon=360.0+rlon
               rip=(rlon-lonMin)/dlon+1
               rjp=(latMax-rlat)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             else
               write(*,*) ' Unknown Mosaic format !!'
               stop 123
             endif
             if( ip >= 1 .and. ip <= readref%mscNlon ) then
             if( jp >= 1 .and. jp <= readref%mscNlat ) then
! inside mosaic domain
               ipp1=min(ip+1,readref%mscNlon)
               jpp1=min(jp+1,readref%mscNlat)
               w1=(1.0-dip)*(1.0-djp)
               w2=dip*(1.0-djp)
               w3=dip*djp
               w4=(1.0-dip)*djp
               ref1=mscValue(ip,jp)
               ref2=mscValue(ipp1,jp)
               ref3=mscValue(ipp1,jpp1)
               ref4=mscValue(ip,jpp1)
               if(ref1 > readref%rthresh_ref .and. ref2 > readref%rthresh_ref .and.  &
                  ref3 > readref%rthresh_ref .and. ref4 > readref%rthresh_ref ) then
                  ref3d(i,j,kk)=(ref1*w1+ref2*w2+ref3*w3+ref4*w4)/float(readref%var_scale)
               elseif(ref1 > readref%rthresh_miss .and. ref2 > readref%rthresh_miss .and.  &
                  ref3 > readref%rthresh_miss .and. ref4 > readref%rthresh_miss ) then
                  ref3d(i,j,kk)=-99.0   ! clear
               else
                  ref3d(i,j,kk)=-999.0  ! no observation
               endif
             endif
             endif
          ENDDO
          ENDDO
      ENDDO  ! mscNlev

      deallocate(mscValue)
  else
     ref3d=-999.0
!     write(*,*) trim(readref%mosaicfile), '   does not exist!!!'
  endif

end subroutine mosaic2grid
