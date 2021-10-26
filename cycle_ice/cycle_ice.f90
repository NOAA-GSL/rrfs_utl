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
  integer :: nx,ny,nz,n,i,j,k
!
!  from wrf netcdf 
  real(8), allocatable,target :: field2d8b(:,:)
  real(8), allocatable,target :: field3d8b(:,:,:)
  real(4), allocatable :: field2d(:,:)
  character(len=30) :: varname
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
    call rrfs%get_dim("zaxis_1",nz)
    write(*,*) 'nx_rrfs,ny_rrfs=',nx,ny,nz
    call rrfs%close()

    allocate(field2d8b(nx,ny))

    do n=1,6
       if(n==1) varname='hice' 
       if(n==2) varname='fice' 
       if(n==3) varname='slmsk' 
       if(n==4) varname='qwv_surf_ice' 
       if(n==5) varname='clw_surf_ice' 
       if(n==6) varname='sfalb_ice' 
       call rrfs%open('old.sfc_data.nc',"r",200)
       call rrfs%get_var(trim(varname),nx,ny,field2d8b)
       call rrfs%close()
!
       call rrfs%open('sfc_data.nc',"w",200)
       call rrfs%replace_var(trim(varname),nx,ny,field2d8b)
       call rrfs%close()
    enddo
!
!
    allocate(field2d(nx,ny))
!
    call rrfs%open('sfc_data.nc',"r",200)
    call rrfs%get_var('fice',nx,ny,field2d8b)
    field2d=field2d8b
    call rrfs%get_var('tisfc',nx,ny,field2d8b)
    call rrfs%close()
    do j=1,ny
      do i =1,nx
        if(field2d(i,j) > 1.0e-10) field2d8b(i,j)=min(271.2, field2d8b(i,j))
      enddo
    enddo
!
    call rrfs%open('sfc_data.nc',"w",200)
    call rrfs%replace_var('tisfc',nx,ny,field2d8b)
    call rrfs%close()
    deallocate(field2d8b)
!
!
    allocate(field3d8b(nx,ny,nz))
    call rrfs%open('sfc_data.nc',"r",200)
    call rrfs%get_var('tiice',nx,ny,nz,field3d8b)
    call rrfs%close()
!
    do j=1,ny
      do i =1,nx
        if(field2d(i,j) > 1.0e-10) then
          do k=1,nz
            field3d8b(i,j,k)=min(271.2, field3d8b(i,j,k))
          enddo
        endif
      enddo
    enddo
!
    call rrfs%open('sfc_data.nc',"w",200)
    call rrfs%replace_var('tiice',nx,ny,nz,field3d8b)
    call rrfs%close()
!
    deallocate(field2d)
    deallocate(field3d8b)
!
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

endif ! mype==0

call MPI_FINALIZE(ierror)
!
end program 
