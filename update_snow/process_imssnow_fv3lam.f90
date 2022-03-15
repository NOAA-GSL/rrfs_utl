PROGRAM process_imssnow
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-02-11
!
! ABSTRACT: 
!     This appllication uses NESDIS SNOW/ICE data from a grib file to
!          update FV3LAM ice and snow fields 
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  imssnow
!   OUTPUT FILES: updated surface file
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

  use mpi
  use module_imssnow, only : type_imssnow
  use module_ncio, only: ncio

  implicit none
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
  integer,allocatable :: mpi_layout_begin(:),mpi_layout_end(:)

  type(ncio) :: fv3grid
  type(type_imssnow) :: imssnow
!
! grid
!
  integer :: nlon,nlat
  real,allocatable :: xlon(:,:)    !
  real,allocatable :: ylat(:,:)    !
  integer,allocatable :: xland(:,:)   !
  real,allocatable :: xlandIMS(:,:)   !
  real,allocatable :: snowice(:,:)    ! snow/ice in RR 
  real(8),allocatable :: tmp8b(:,:)
!
!
  integer :: num_args
  character(len=40), dimension(:), allocatable :: args
  integer            :: fv3_io_layout_y
  integer,allocatable :: fv3_layout_begin(:),fv3_layout_end(:)
!
  character(len=80) :: snowfile
  integer :: ireturn
  integer :: id,ix,iy
  character(len=4) :: clayout_y
  character(len=40) :: thisfv3file
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)
!
  write(*,*) 'number of cores=',npe

!  Get file names from command line arguements
   num_args = command_argument_count()
   if(num_args==1) then
     allocate(args(num_args))
     do ix = 1, num_args
       call get_command_argument(ix,args(ix))
     end do
     clayout_y=trim(args(1))
     read(clayout_y,'(I2)') fv3_io_layout_y
   else
     fv3_io_layout_y=1
   endif
   write(*,*) 'subdomain number=',fv3_io_layout_y

   if(mype < fv3_io_layout_y) then
      write(thisfv3file,'(a,I3.3)') 'stdout_update_snowice.',mype
      open(6, file=trim(thisfv3file),form='formatted',status='unknown')
      write(6,*) '===> start process for core = ', mype
   endif

! distribute subdomain to cores
   allocate(mpi_layout_begin(fv3_io_layout_y))
   allocate(mpi_layout_end(fv3_io_layout_y))
   mypeLocal=mype+1
   if(npe==1) then
      mpi_layout_begin(mypeLocal)=1
      mpi_layout_end(mypeLocal)=fv3_io_layout_y
   else
      mpi_layout_begin=0
      mpi_layout_end=0
      do ix=1,fv3_io_layout_y
         iy=mod(ix-1,npe)+1
         mpi_layout_end(iy)=mpi_layout_end(iy)+1
      enddo
      iy=0
      do ix=1,npe
         if(mpi_layout_end(ix) > 0) then
           mpi_layout_begin(ix)=iy+1
           mpi_layout_end(ix)=iy+mpi_layout_end(ix)
           iy=mpi_layout_end(ix)
         else
           mpi_layout_begin(ix)=0
           mpi_layout_end(ix)=0
         endif
      enddo
   endif
   write(6,'(a)') 'begin and end domain index for each core:'
   write(6,'(a20,20I5)') 'mpi_layout_begin=',mpi_layout_begin
   write(6,'(a20,20I5)') 'mpi_layout_end=',mpi_layout_end

  if(mypeLocal <= fv3_io_layout_y) then

     snowfile='imssnow2'
     call imssnow%init(snowfile,ireturn)
     if(ireturn==0) then
        call imssnow%read_imssnow()

        do id=mpi_layout_begin(mypeLocal),mpi_layout_end(mypeLocal)

           write(6,*) 'process subdomain ',id,' with core ',mype
           if(fv3_io_layout_y > 1) then
             write(thisfv3file,'(a,I4.4)') 'fv3_grid_spec.',id-1
           else
             thisfv3file='fv3_grid_spec'
           endif
           call fv3grid%open(trim(thisfv3file),'r',200)
           call fv3grid%get_dim("grid_xt",nlon)
           call fv3grid%get_dim("grid_yt",nlat)
           write(6,*) 'grid dimension =',nlon,nlat
           allocate(xlon(nlon,nlat))
           allocate(ylat(nlon,nlat))
           call fv3grid%get_var("grid_lont",nlon,nlat,xlon)
           call fv3grid%get_var("grid_latt",nlon,nlat,ylat)
           call fv3grid%close
! read xland
           if(fv3_io_layout_y > 1) then
             write(thisfv3file,'(a,I4.4)') 'sfc_data.nc.',id-1
           else
             thisfv3file='sfc_data.nc'
           endif
           allocate(xland(nlon,nlat))
           allocate(tmp8b(nlon,nlat))
           call fv3grid%open(trim(thisfv3file),'r',200)
           call fv3grid%get_var("slmsk",nlon,nlat,tmp8b)
           call fv3grid%close
           xland=int(tmp8b)
           deallocate(tmp8b)
!
! map to grid
!
           allocate(snowice(nlon,nlat))
           allocate(xlandIMS(nlon,nlat))
           snowice=0  ! pecentage
           xlandIMS=0
           call imssnow%map2grid(xland,nlon,nlat,xlon,ylat,snowice,xlandIMS)
           deallocate(xlon)
           deallocate(ylat)
           deallocate(xland)
           write(6,*) 'get snow on grid=',maxval(snowice),minval(snowice)
!
!
!  trim snow cover field based on NESDIS snow cover data
!
           call update_snow_fv3lam(snowice, xland, nlon,nlat,id,fv3_io_layout_y)
!
           deallocate(snowice)
           deallocate(xlandIMS)
        enddo ! domain

        call imssnow%close()
     else
        write(6,*) 'file is not exist: ',trim(snowfile)
        write(6,*) 'stop update snow/ice'
     endif

  endif !mypeLocal <= fv3_io_layout_y

  deallocate(mpi_layout_begin)
  deallocate(mpi_layout_end)

!
  call MPI_FINALIZE(ierror)
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="
!

END PROGRAM process_imssnow

