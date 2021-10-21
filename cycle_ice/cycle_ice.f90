PROGRAM cycle_ice
!
! read in GFS ice and repalce ice in RRFS  
!
!   Ming Hu, 2021-10-20
!
!
  use mpi
  use module_ncio, only : ncio
!
  implicit none
!
  type(ncio)     :: rrfs
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
  integer :: nx,ny
!
!  from wrf netcdf 
  real(8), allocatable,target :: field2d8b(:,:)
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
! NCEP LSF has to use all cores allocated to run this application 
! but this if check can make sure only one core run through the real code.
if(mype==0) then
!
!========
!  
! read in rrfs latlon
    call rrfs%open('old.sfc_data.nc',"r",200)
    call rrfs%get_dim("xaxis_1",nx)
    call rrfs%get_dim("yaxis_1",ny)
    write(*,*) 'nx_rrfs,ny_rrfs=',nx,ny

    allocate(field2d8b(nx,ny))
    call rrfs%get_var("hice",nx,ny,field2d8b)
    call rrfs%close()
!
    call rrfs%open('sfc_data.nc',"w",200)
    call rrfs%replace_var("hice",nx,ny,field2d8b)
    call rrfs%close()
!
!  fice
    call rrfs%open('old.sfc_data.nc',"r",200)
    call rrfs%get_var("fice",nx,ny,field2d8b)
    call rrfs%close()
    call rrfs%open('sfc_data.nc',"w",200)
    call rrfs%replace_var("fice",nx,ny,field2d8b)
    call rrfs%close()
    deallocate(field2d8b)
!
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

endif ! mype==0

call MPI_FINALIZE(ierror)
!
end program 
