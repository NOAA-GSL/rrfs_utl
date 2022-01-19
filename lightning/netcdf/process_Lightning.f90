program process_Lightning
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2022-01-18
!
! ABSTRACT: 
!     This routine read in lightning data and 
!     map them into FV3LAM ESG grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  NLDN lightning data
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
  use module_ncio, only : ncio
  use pesg, only : gtoxm_ak_dd,xmtog_ak_dd
  use module_esggrid_util, only: edp,esggrid_util

  implicit none
  INCLUDE 'netcdf.inc'
!
  type(ncio) :: rrfs
  type(esggrid_util) :: esggrid
!
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror

  real     :: rad2deg = 180.0/3.1415926
!
  character*256 output_file
!
!  grid and map
  CHARACTER*180   geofile

! esg grid
  real(edp)              :: dlat,dlon
  real(edp)              :: xi,yj

!  grid
  integer(i_kind) :: nlonfv3,nlatfv3
  real,allocatable:: xlonfv3(:,:)    !
  real,allocatable:: ylatfv3(:,:)    !

!
!  For lightning data
!
  INTEGER ::    numStrike
  CHARACTER*180   lightsngle
  real,allocatable:: llon(:)    !
  real,allocatable:: llat(:)    !
  integer,allocatable:: ltime(:)   !
  integer,allocatable:: lstrike(:) !
  character*21,allocatable:: ctime(:)   !
  real :: rtmp
  integer,allocatable:: lquality(:) !

  REAL, allocatable :: lightning(:,:)   ! lightning  strakes
  REAL(r_kind), allocatable :: lightning_out(:,:)   ! lightning  strakes

  integer :: numNLDN_all, numNLDN_used
  integer :: numAlaska_all, numAlask_used
!
!! Declare namelists 
!
! SETUP (general control namelist) :
!
  character*10 :: analysis_time
  integer      :: NLDN_filenum
  logical      :: IfAlaska
  character(len=20) :: grid_type
  namelist/setup/analysis_time, NLDN_filenum,IfAlaska,grid_type
!
!  ** misc
      
  CHARACTER*180   workpath

  real :: LAT_LL_P,LON_LL_P
  real(r_kind) :: XC,YC

  integer i,j,igrid,jgrid,nt
  integer :: ii,jj

  integer :: NCID, istatus

  integer :: numlightning,idate,iii


!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
!  call MPI_INIT(ierror) 
!  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
!  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  numNLDN_all=0
  numNLDN_used=0
  numAlaska_all=0
  numAlask_used=0
  iii=0
  read(5,setup)

! define esg grid
  
  call esggrid%init(grid_type)

!
! set geogrid fle name
!
  write(geofile,'(a,a)') './', 'fv3sar_grid_spec.nc'
  call rrfs%open(trim(geofile),"r",200)
  call rrfs%get_dim("grid_xt",nlonfv3)
  call rrfs%get_dim("grid_yt",nlatfv3)
  write(*,*) 'nx_rrfs,ny_rrfs=',nlonfv3,nlatfv3

!  allocate(xlonfv3(nlonfv3,nlatfv3))
!  allocate(ylatfv3(nlonfv3,nlatfv3))
!  call rrfs%get_var("grid_lont",nlonfv3,nlatfv3,xlonfv3)
!  call rrfs%get_var("grid_latt",nlonfv3,nlatfv3,ylatfv3)
!  write(*,*) 'FV3SAR grid'
!  write(*,*) 'nlonfv3,nlatfv3=', nlonfv3,nlatfv3
!  write(*,*) 'max, min lon=', maxval(xlonfv3),minval(xlonfv3)
!  write(*,*) 'max, min lat=', maxval(ylatfv3),minval(ylatfv3)
  call rrfs%close()

  allocate(lightning(nlonfv3,nlatfv3))
  lightning=0

!
!  process NLDN data
!
  DO nt=1,NLDN_filenum
    write(lightsngle,'(a,I1)') 'NLDN_lightning_',nt
    write(*,*) trim(lightsngle)
    call ifexist_file(trim(lightsngle),istatus)
    if (ISTATUS .NE. NF_NOERR) CYCLE

    call GET_DIM_ATT_NLDN(lightsngle,numStrike)
    write(*,*) 'number of strikes=', nt, numStrike

    allocate(llon(numStrike))
    allocate(llat(numStrike))
    allocate(ltime(numStrike))
    allocate(lStrike(numStrike))

    call GET_lightning_NLDN(lightsngle,numStrike,llon,llat,ltime,lStrike)
!
!  check quality
!
    allocate(lquality(numStrike))
    lquality = 0    ! 0 good data,  > 0 bad data
    call Check_NLDN(numStrike,llon,llat,ltime,lstrike,lquality)

    do i=1,numStrike

      if(lquality(i) == 0 ) then
        dlon=llon(i)
        dlat=llat(i)
        call esggrid%lltoxy(dlon,dlat,xc,yc)

        igrid = int(XC+0.5)
        jgrid = int(YC+0.5)
        if( (igrid > 0 .and. igrid< nlonfv3).and.  &
            (jgrid > 0 .and. jgrid< nlatfv3)) then 
            lightning(igrid,jgrid) = lightning(igrid,jgrid) + 1
            numNLDN_used=numNLDN_used+1
        endif
      endif

    enddo

    deallocate(llon)
    deallocate(llat)
    deallocate(ltime)
    deallocate(lStrike)
    deallocate(lquality)
  enddo ! nt
!
!  statistic
!
  write(*,*) ' The total number of NLDN data is:', numNLDN_all
  write(*,*) ' The number of NLDN data used is:', numNLDN_used

!  for FV3 LAM

     allocate(lightning_out(3,nlonfv3*nlatfv3))
     numlightning=0
     do j=1,nlatfv3
     do i=1,nlonfv3
       if(lightning(i,j) > 0 ) then
         numlightning=numlightning+1
         lightning_out(1,numlightning)=float(i)
         lightning_out(2,numlightning)=float(j)
         lightning_out(3,numlightning)=lightning(i,j)
         if(lightning_out(3,numlightning) > 1000.0 ) then
            lightning_out(3,numlightning)=1000.0
            write(6,*) 'high lightning strokes=',lightning(i,j),i,j
         endif
       endif
     enddo
     enddo

     write(*,*) 'Write out results for FV3 LAM:',numlightning
     OPEN(10,file='LightningInFV3LAM.dat',form='unformatted')
      write(10) 3,nlonfv3,nlatfv3,numlightning,1,2
      write(10) ((real(lightning_out(i,j)),i=1,3),j=1,numlightning)
      write(10) lightning
     close(10)

!     deallocate(xlonfv3)
!     deallocate(ylatfv3)
!  
   read(analysis_time,'(I10)') idate
   write(6,*) 'cycle time is :', idate

   write(6,*) ' write lightning in BUFR'
   call write_bufr_lightning(1,nlonfv3,nlatfv3,numlightning,lightning_out,idate)
     deallocate(lightning_out)

!  call MPI_FINALIZE(ierror)
!
end program process_Lightning 
