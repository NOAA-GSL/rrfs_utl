PROGRAM gen_cs
!
! generate cross section
!
!   Ming Hu, 2021-11-27
!
!
  use mpi
  use module_ncio, only : ncio
!
  implicit none
!
  type(ncio)     :: model
!
  real,parameter :: g=9.80665
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
!
  integer :: istart,iend,jstart,jend,jcs
  real :: cstop,dz
  character (len=250) :: bkfile,anafile,csfile
  namelist/setup/ istart,iend,jstart,jend,cstop,dz,jcs,&
                  bkfile,anafile,csfile
!
  logical :: ifexist
!
  integer :: nx,ny,nz,csnx,csnz
  integer :: n,i,j,k
!
!  
  real(4), allocatable :: field3d4b(:,:,:)
  real(4), allocatable :: tmp3d4b(:,:,:)
!
  real(4), allocatable :: ana_bk_cs2d(:,:)
  real(4), allocatable :: cs2d(:,:)
!
  real(4), allocatable :: height(:,:)
  real(4), allocatable :: csheight(:)
  real(4), allocatable :: weight_cs2d(:,:)
  character(len=30) :: varname
!
  real(4) :: wght, csh
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
   istart=150
   iend=250
   jstart=540
   jend=640
   jcs=(jend-jstart)/2
   cstop=2000.0
   dz=20.0
   bkfile="wrf_inout_bk"
   anafile="wrf_inout"
   csfile="csfile"
!
!  read namelist
!
   inquire(file='namelist_cs', EXIST=ifexist )
   if(ifexist) then
      open(10,file='namelist_cs',status='old')
         read(10,setup)
      close(10)
      write(*,*) 'Namelist setup are:'
      write(*,setup)
   else
      write(*,*) 'No namelist file exist, use default values'
   endif
!
!========
!  
! read in model dimension
   call model%open(trim(bkfile),"r",200)
   call model%get_dim("west_east",nx)
   call model%get_dim("south_north",ny)
   call model%get_dim("bottom_top",nz)
   write(*,*) 'nx_model,ny_model=',nx,ny,nz
   csnx=iend-istart+1
   csnz=int(cstop/dz)+1

   allocate(height(csnx,nz))
   allocate(csheight(csnz))
   allocate(field3d4b(nx,ny,nz+1))
   allocate(tmp3d4b(nx,ny,nz+1))

   call model%get_var("PH",nx,ny,nz+1,field3d4b)
   call model%get_var("PHB",nx,ny,nz+1,tmp3d4b)
   field3d4b=(field3d4b+tmp3d4b)/g
   call model%close()

   do k=1,nz
      do i=istart,iend
        height(i-istart+1,k)=(field3d4b(i,jcs,k+1)+field3d4b(i,jcs,k))*0.5
      enddo
   enddo

   write(*,*) 'height=',height(1,:)
   
   csheight(1)=10.0
   do k=2,csnz
      csheight(k)=csheight(k-1)+dz
   enddo
   write(*,*) 'csheight=',csheight
! calculate position of each cs point
   allocate(weight_cs2d(csnx,csnz))
   do j=1,csnz
      csh=csheight(j)
      do i=1,csnx
        wght=-1.0
        do k=1,nz-1
          if(csh >= height(i,k) .and. csh < height(i,k+1)) then
             wght=k+(csh-height(i,k))/(height(i,k+1)-height(i,k))
          endif
        enddo
        if(csh < height(i,1) .and. (csh+dz*0.5) > height(i,1)) wght=1.0
        weight_cs2d(i,j)=wght 
      enddo
   enddo
!
   deallocate(field3d4b,tmp3d4b)
   allocate(field3d4b(nx,ny,nz))
   allocate(tmp3d4b(nx,ny,nz))
  
   allocate(ana_bk_cs2d(csnx,nz))
   allocate(cs2d(csnx,csnz))

   open(11,file=trim(csfile),form='unformatted')
   write(11) csnx,csnz
   do n=1,2
      if(n==1) varname='T' 
      if(n==2) varname='QVAPOR' 
      call model%open(trim(bkfile),"r",200)
      call model%get_var(trim(varname),nx,ny,nz,tmp3d4b)
      call model%close()
!
      call model%open(trim(anafile),"w",200)
      call model%get_var(trim(varname),nx,ny,nz,field3d4b)
      call model%close()
!
      do k=1,nz
         do i=istart,iend
           ana_bk_cs2d(i-istart+1,k)=field3d4b(i,jcs,k)-tmp3d4b(i,jcs,k)
         enddo
      enddo
!
      if(n==2) ana_bk_cs2d=ana_bk_cs2d*1000.0
!      
      cs2d=-999.0
      do j=1,csnz
        do i=1,csnx
           wght=weight_cs2d(i,j)
           if(wght > 0.1 .and. wght < nz) then
              k=int(wght)
              wght=wght-k
              cs2d(i,j)=ana_bk_cs2d(i,k+1)*wght + ana_bk_cs2d(i,k)*(1.0-wght)
           end if
        enddo
     enddo 
     write(11) trim(varname)
     write(11) cs2d
   enddo
!
   close(11)
   deallocate(ana_bk_cs2d)
   deallocate(cs2d)
   deallocate(field3d4b)
   deallocate(tmp3d4b)
!
  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

endif ! mype==0

call MPI_FINALIZE(ierror)
!
end program 
