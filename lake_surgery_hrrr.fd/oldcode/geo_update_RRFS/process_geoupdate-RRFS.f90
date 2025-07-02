PROGRAM process_SST
!
!   PRGMMR: Ming Hu  nd Tanya Smirnova   ORG: GSD        DATE: 2010-09-25
!
! ABSTRACT: 
!     This routine reads in SST
!
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  imssnow
!
!   OUTPUT FILES:  RRimssnow
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

!  use mpi

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror

! RR grid
  integer :: nlon,nlat
! RR in LC
!  parameter (nlon=648,nlat=647)
! RR in RLL
!  parameter (nlon=758,nlat=567)
! parameter (nlon=3640, nlat=2520)
  parameter (nlon=3950, nlat=2700)
  real,allocatable  :: xlon(:,:)    !
  real,allocatable  :: ylat(:,:)    !
  real,allocatable  :: xland(:,:)   !
  real,allocatable  :: vegtyp(:,:)   !
  real,allocatable  :: landusefr(:,:,:)   !
  character(len=120) :: flnm1,flnm2


  real, allocatable :: sstRR(:,:)    ! sst in RR 
  real, allocatable :: sstGlobal(:,:)  ! sst from global dataset
  integer, allocatable :: imaskSST(:,:)

!
  character*80 input_file
  integer :: istatus
  integer :: i,j,iwater,ilake,iice
!
  INTEGER :: iyear, imonth, iday, ihr

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
!  call MPI_INIT(ierror)
!  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
!  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
  if(mype==0) then

!  flnm1 = 'C3445_oro_data.tile7.halo0.nc'
   flnm1 = 'orofile.nc'

!  call GET_DIM_ATT_geo(flnm1,nlon,nlat)
! call GET_DIM_ATT_geo('./RRFS_noGreatLakes.nc',nlon,nlat)
!  call GET_DIM_ATT_geo('./geo_em.d01.nam_terrainsmooth.nc',nlon,nlat)
  write(*,*) 'grid dimension =',nlon,nlat
  allocate(xlon(nlon,nlat))
  allocate(ylat(nlon,nlat)) 
  allocate(xland(nlon,nlat))
  allocate(vegtyp(nlon,nlat))
  allocate(landusefr(nlon,nlat,21))

  call update_geo_rrfs(nlon, nlat)
!
  endif ! mype==0

!  call MPI_FINALIZE(ierror)
!
END PROGRAM process_SST
