PROGRAM check_process_imssnow
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

  use module_ncio, only: ncio

  implicit none
! MPI variables
  type(ncio) :: fv3grid
!
! grid
!
  integer :: nlon,nlat
  real,allocatable :: xlon(:,:)    !
  real,allocatable :: ylat(:,:)    !
  real(8),allocatable :: slmsk(:,:)
  real(8),allocatable :: sncovr(:,:)
  real(8),allocatable :: snodl(:,:)
  real(8),allocatable :: snwdph(:,:)
  real(8),allocatable :: weasdl(:,:)
  real(8),allocatable :: tsfc(:,:)
  real(8),allocatable :: tsfcl(:,:)
  real(8),allocatable :: tsnow_land(:,:)
  real(8),allocatable :: tslb(:,:,:)
!
!
  integer :: num_args
  character(len=40), dimension(:), allocatable :: args
  integer            :: fv3_io_layout_y
!
  character(len=80) :: snowfile
  integer :: ireturn
  integer :: id,ix,iy
  integer :: i,j
  character(len=4) :: clayout_y
  character(len=40) :: thisfv3file
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program
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

! distribute subdomain to cores
   !do id=1, fv3_io_layout_y
   do id=5, 7
      write(6,*) 'process subdomain ',id
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
      do j=1,nlat
      do i=1,nlon
         if(abs(ylat(i,j)-42.05)<0.05 .and. abs(xlon(i,j)-286.75) <0.05) then
            write(*,*) i,j,xlon(i,j),ylat(i,j)
         endif
      enddo
      enddo
! read xland
      if(fv3_io_layout_y > 1) then
         write(thisfv3file,'(a,I4.4)') 'sfc_data.nc.',id-1
      else
         thisfv3file='sfc_data.nc'
      endif
      allocate(slmsk(nlon,nlat))
      allocate(sncovr(nlon,nlat))
      allocate(snodl(nlon,nlat))
      allocate(snwdph(nlon,nlat))
      allocate(weasdl(nlon,nlat))
      allocate(tsfc(nlon,nlat))
      allocate(tsfcl(nlon,nlat))
      allocate(tsnow_land(nlon,nlat))
      allocate(tslb(nlon,nlat,9))
      call fv3grid%open(trim(thisfv3file),'r',200)
      call fv3grid%get_var("slmsk",nlon,nlat,slmsk)
      call fv3grid%get_var("sncovr",nlon,nlat,sncovr)
      call fv3grid%get_var("snodl",nlon,nlat,snodl)
      call fv3grid%get_var("snwdph",nlon,nlat,snwdph)
      call fv3grid%get_var("weasdl",nlon,nlat,weasdl)
      call fv3grid%get_var("tsfc",nlon,nlat,tsfc)
      call fv3grid%get_var("tsfcl",nlon,nlat,tsfcl)
      call fv3grid%get_var("tsnow_land",nlon,nlat,tsnow_land)
      call fv3grid%get_var("tslb",nlon,nlat,9,tslb)
      do j=1,nlat
      do i=1,nlon
!         if(abs(ylat(i,j)-41.0215186174087)<0.1 .and. abs(xlon(i,j)-284.509251283993) <0.1) then
         if((abs(i-2742)<3 .and. abs(j-1) < 3 .and. id==6) .or. &
            (abs(i-2742)<3 .and. abs(j-nlat)< 3 .and. id==5))  then
            write(*,*) 
            write(*,'(3I5,2f12.5)') id,i,j,xlon(i,j),ylat(i,j)
            write(*,*) 'slmsk,sncovr,snodl=',slmsk(i,j),sncovr(i,j)
            write(*,*) '      snodl,snwdph=',snodl(i,j),snwdph(i,j)
            write(*,*) '       weasdl,tsfc=',weasdl(i,j),tsfc(i,j)
            write(*,*) '  tsfcl,tsnow_land=',tsfcl(i,j),tsnow_land(i,j)
            write(*,*) '              tslb=',tslb(i,j,1:9)
         endif
      enddo
      enddo
      call fv3grid%close
      
      deallocate(slmsk)
      deallocate(sncovr)
      deallocate(snodl)
      deallocate(snwdph)
      deallocate(weasdl)
      deallocate(tsfc)
      deallocate(tsfcl)
      deallocate(tsnow_land)
      deallocate(tslb)
      deallocate(xlon)
      deallocate(ylat)
   enddo
!

END PROGRAM check_process_imssnow

