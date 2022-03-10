program  process_metar_cloud
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2009-09-04
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products and 
!     interpolate them into GSI mass grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  
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
  use kinds, only: r_kind,i_kind,r_single
  use module_ncio, only : ncio

  use cld_parm_array_mod, only: region_dy,region_dx
  use cld_parm_array_mod, only: metar_impact_radius
  use cld_parm_array_mod, only: l_metar_impact_radius_change, &
                            metar_impact_radius_max,metar_impact_radius_min,&
                            metar_impact_radius_max_height,metar_impact_radius_min_height
  use cld_parm_array_mod, only: init_cld_parm

  use cld_parm_array_mod, only : obstype, sis, nchanl,nreal,ilat,ilon,ndata
  use cld_parm_array_mod, only : cdata_regular

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! MPI variables
  integer :: npe, mype,ierror
!
  type(ncio) :: rrfs
!
!  grid
  integer(i_kind) :: nlon,nlat
!
  real(r_kind),allocatable,dimension(:,:):: cdata_fv3
  integer,allocatable,dimension(:,:):: index_regular
!

  CHARACTER*180   geofile
!
!  For NASA LaRC 
!
  CHARACTER*180   workPath
!     ****VARIABLES FOR THIS NETCDF FILE****
!
! namelist
!
  integer :: analysis_time
  integer :: analysis_minute
  character(len=100) :: prepbufrfile
  real(r_kind)       :: twindin
  character(len=20) :: grid_type
  namelist/setup/ grid_type,analysis_time,analysis_minute,prepbufrfile,twindin,&
                  metar_impact_radius,l_metar_impact_radius_change, &
                  metar_impact_radius_max,metar_impact_radius_min, &
                  metar_impact_radius_max_height,metar_impact_radius_min_height
!
!  ** misc
      
  integer i,j,ii,jj,ij

  integer :: NCID
  logical :: ifexist
  integer :: lunout

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) then
     lunout=15
     call init_cld_parm
!
!  get namelist
!
     analysis_time=2018051718
     prepbufrfile='prepbufr'
     analysis_minute=0
     twindin=0.5
     grid_type="none"
 
     inquire(file='namelist.metarcld', EXIST=ifexist )
     if(ifexist) then
       open(10,file='namelist.metarcld',status='old')
          read(10,setup)
       close(10)
       write(*,*) 'Namelist setup are:'
       write(*,setup)
     else
       write(*,*) 'No namelist file exist, use default values'
       write(*,*) analysis_time
     endif
!
! set geogrid fle name
!
     write(geofile,'(a,a)') './', 'fv3sar_grid_spec.nc'
     call rrfs%open(trim(geofile),"r",200)
     call rrfs%get_dim("grid_xt",nlon)
     call rrfs%get_dim("grid_yt",nlat)
     write(*,*) 'nx_rrfs,ny_rrfs=',nlon,nlat
!
     call read_prepbufr_metarcld(prepbufrfile,analysis_time,analysis_minute,&
                                 twindin,nlon,nlat,grid_type)

     write(*,*) 'number of cloud on fv3 grid=',ndata
     write(*,*) obstype,nreal,nchanl,ilat,ilon,sis
     open(lunout,file='fv3_metarcloud.bin',form='unformatted')
        write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndata
        write(lunout) cdata_regular
     close(lunout)
!
     deallocate(cdata_regular)
!
     write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

  endif ! if mype==0 

  call MPI_FINALIZE(ierror)

end program process_metar_cloud
