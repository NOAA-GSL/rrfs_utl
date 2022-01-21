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
!  ** misc
!      
  character*80 outfile
  integer i,j,k,kk
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
  open(15, file='mosaic.namelist')
    read(15,setup)
  close(15)
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
     deallocate(xlon)
     deallocate(ylat)
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
        else
           write(outfile,'(a,a,I4.4)') './', 'RefInGSI3D.dat.',id
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
        DO j=2,nlat-1
        DO i=2,nlon-1
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
       deallocate(ref0)
     endif
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
