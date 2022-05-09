module module_mpi_arrange
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-03-08
!
! ABSTRACT: 
!     This module figure out the piece of fv3lam file each core to read
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES: 
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!
!$$$
!
!_____________________________________________________________________

  implicit none

  integer,parameter :: max_varname_length=20
!
! Rset default to private
!

  private

  public :: mpi_io_arrange

  type :: mpi_io_arrange

      integer,allocatable :: fileid(:)
      character(len=max_varname_length),allocatable :: varname(:)
      integer,allocatable :: vartype(:)
      integer,allocatable :: lvlbegin(:),lvlend(:)
      integer,allocatable :: nx(:),ny(:)
      integer :: ntotalcore

    contains
      procedure :: init
      procedure :: arrange
      procedure :: close
  end type mpi_io_arrange
!
! constants
!
contains

  subroutine init(this,ntotalcore)
!                .      .    .                                       .
! subprogram: 
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    implicit none

    integer, intent(in)          :: ntotalcore
    class(mpi_io_arrange) :: this
!
    this%ntotalcore=ntotalcore

    allocate(this%fileid(ntotalcore))
    allocate(this%varname(ntotalcore))
    allocate(this%vartype(ntotalcore))
    allocate(this%lvlbegin(ntotalcore))
    allocate(this%lvlend(ntotalcore))
    allocate(this%nx(ntotalcore))
    allocate(this%ny(ntotalcore))

  end subroutine init

  subroutine close(this)
    implicit none
    class(mpi_io_arrange) :: this

    if(allocated(this%fileid))   deallocate(this%fileid)
    if(allocated(this%varname))  deallocate(this%varname)
    if(allocated(this%vartype))  deallocate(this%vartype)
    if(allocated(this%lvlbegin)) deallocate(this%lvlbegin)
    if(allocated(this%lvlend))   deallocate(this%lvlend)
    if(allocated(this%nx))       deallocate(this%nx)
    if(allocated(this%ny))       deallocate(this%ny)
    this%ntotalcore=0

  end subroutine close

  subroutine arrange(this,ncfs_dyn,ncfs_tracer,ncfs_sfc)
!                .      .    .                                       .
! subprogram: 
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    use module_ncfile_stat, only : ncfile_stat
    implicit none

    class(mpi_io_arrange) :: this
    type(ncfile_stat),intent(in) :: ncfs_dyn
    type(ncfile_stat),intent(in) :: ncfs_tracer
    type(ncfile_stat),intent(in) :: ncfs_sfc
!
! MPI distribution array
    integer,allocatable :: fileid(:)
    character(len=max_varname_length),allocatable :: varname(:)
    integer,allocatable :: vartype(:)
    integer,allocatable :: lvlbegin(:),lvlend(:)
    integer,allocatable :: nx(:),ny(:)
    integer :: ntotalcore
    integer :: nlvl2d,nlvl3d
    integer :: nlvlcore
    integer :: nlvlmax
    integer,allocatable :: nlvl3d_list(:)
    integer :: mynlvl3d
!
    integer :: i,j,k,nn,n3d,nz
!
! start  allocate memory
!
    ntotalcore=this%ntotalcore
    allocate(fileid(ntotalcore))
    allocate(varname(ntotalcore))
    allocate(vartype(ntotalcore))
    allocate(lvlbegin(ntotalcore))
    allocate(lvlend(ntotalcore))
    allocate(nx(ntotalcore))
    allocate(ny(ntotalcore))
!
! find number of 3d and 2 d fields
!
    nlvl2d=0
    nlvl3d=0
    nlvlmax=0
    if(ncfs_tracer%num_totalvl > 0) then
        do k=1,ncfs_tracer%numvar
           if(ncfs_tracer%num_dim(k) == 3) nlvl3d=nlvl3d+1
           if(ncfs_tracer%num_dim(k) == 2) nlvl2d=nlvl2d+1
           if(nlvlmax < ncfs_tracer%dim_3(k)) nlvlmax=ncfs_tracer%dim_3(k)
        enddo
    endif
    if(ncfs_dyn%num_totalvl > 0) then
        do k=1,ncfs_dyn%numvar
           if(ncfs_dyn%num_dim(k) == 3) nlvl3d=nlvl3d+1
           if(ncfs_dyn%num_dim(k) == 2) nlvl2d=nlvl2d+1
           if(nlvlmax < ncfs_dyn%dim_3(k)) nlvlmax=ncfs_dyn%dim_3(k)
        enddo
    endif
    if(ncfs_sfc%num_totalvl > 0) then
        do k=1,ncfs_sfc%numvar
           if(ncfs_sfc%num_dim(k) == 3) nlvl3d=nlvl3d+1
           if(ncfs_sfc%num_dim(k) == 2) nlvl2d=nlvl2d+1
           if(nlvlmax < ncfs_sfc%dim_3(k)) nlvlmax=ncfs_sfc%dim_3(k)
        enddo
    endif
    write(6,*) 'total 2d level=',nlvl2d,' 3d level=',nlvl3d
    write(6,*) 'max level =',nlvlmax
!
!  check if we have enought cores
!
    if( ntotalcore-nlvl2d < nlvl3d ) then
        write(6,*) 'not enough cores for the paralleli IO',ntotalcore-nlvl2d,nlvl3d
        stop 123
    endif
!
! decide how many cores can be used for each 3d field
!
    allocate(nlvl3d_list(nlvl3d))
    nlvlcore=(ntotalcore-nlvl2d)/nlvl3d
    nlvl3d_list=nlvlcore
    nlvlcore=(ntotalcore-nlvl2d) - nlvl3d*nlvlcore 
    if(nlvlcore > 0 ) then
        do k=1,nlvlcore
          nlvl3d_list(k)=nlvl3d_list(k)+1
        enddo
    endif
    write(6,*) nlvl3d_list
!
!  decide which levels of a variable to read for each core
!
    nn=0
    n3d=0
    if(ncfs_dyn%numvar > 0) then
        do i=1,ncfs_dyn%numvar
           if(ncfs_dyn%num_dim(i)==3) then
              nz=ncfs_dyn%dim_3(i)
              n3d=n3d+1
              mynlvl3d=min(nlvl3d_list(n3d),nz)
              nlvlcore=nz/mynlvl3d
              j=nz-nlvlcore*mynlvl3d
              do k=1,mynlvl3d
                 nn=nn+1
                 fileid(nn)=1
                 varname(nn)=trim(ncfs_dyn%list_varname(i))
                 vartype(nn)=ncfs_dyn%vartype(i)
                 nx(nn)=ncfs_dyn%dim_1(i)
                 ny(nn)=ncfs_dyn%dim_2(i)
                 if(k==1) then
                    lvlbegin(nn)=1
                 else
                    lvlbegin(nn)=lvlend(nn-1)+1
                 endif
                 lvlend(nn)=lvlbegin(nn)+nlvlcore-1
                 if(k<=j) lvlend(nn)=lvlend(nn)+1
                 if(k==mynlvl3d) lvlend(nn)=min(lvlend(nn),nz)
              enddo
           elseif(ncfs_dyn%num_dim(i)==2) then
              nn=nn+1
              varname(nn)=trim(ncfs_dyn%list_varname(i))
              vartype(nn)=ncfs_dyn%vartype(i)
              nx(nn)=ncfs_dyn%dim_1(i)
              ny(nn)=ncfs_dyn%dim_2(i)
              lvlbegin(nn)=1
              lvlend(nn)=1
              fileid(nn)=1
           endif
        enddo
    endif

    if(ncfs_tracer%numvar > 0) then
        do i=1,ncfs_tracer%numvar
           if(ncfs_tracer%num_dim(i)==3) then
              nz=ncfs_tracer%dim_3(i)
              n3d=n3d+1
              mynlvl3d=min(nlvl3d_list(n3d),nz)
              nlvlcore=nz/mynlvl3d
              j=nz-nlvlcore*mynlvl3d
              do k=1,mynlvl3d
                 nn=nn+1
                 fileid(nn)=2
                 varname(nn)=trim(ncfs_tracer%list_varname(i))
                 vartype(nn)=ncfs_tracer%vartype(i)
                 nx(nn)=ncfs_tracer%dim_1(i)
                 ny(nn)=ncfs_tracer%dim_2(i)
                 if(k==1) then
                    lvlbegin(nn)=1
                 else
                    lvlbegin(nn)=lvlend(nn-1)+1
                 endif
                 lvlend(nn)=lvlbegin(nn)+nlvlcore-1
                 if(k<=j) lvlend(nn)=lvlend(nn)+1
                 if(k==mynlvl3d) lvlend(nn)=min(lvlend(nn),nz)
              enddo
           elseif(ncfs_tracer%num_dim(i)==2) then
              nn=nn+1
              varname(nn)=trim(ncfs_tracer%list_varname(i))
              vartype(nn)=ncfs_tracer%vartype(i)
              nx(nn)=ncfs_tracer%dim_1(i)
              ny(nn)=ncfs_tracer%dim_2(i)
              lvlbegin(nn)=1
              lvlend(nn)=1
              fileid(nn)=2
           endif
        enddo
    endif

    if(ncfs_sfc%numvar > 0) then
        do i=1,ncfs_sfc%numvar
           if(ncfs_sfc%num_dim(i)==3) then
              nz=ncfs_sfc%dim_3(i)
              n3d=n3d+1
              mynlvl3d=min(nlvl3d_list(n3d),nz)
              nlvlcore=nz/mynlvl3d
              j=nz-nlvlcore*mynlvl3d
              do k=1,mynlvl3d
                 nn=nn+1
                 fileid(nn)=3
                 varname(nn)=trim(ncfs_sfc%list_varname(i))
                 vartype(nn)=ncfs_sfc%vartype(i)
                 nx(nn)=ncfs_sfc%dim_1(i)
                 ny(nn)=ncfs_sfc%dim_2(i)
                 if(k==1) then
                    lvlbegin(nn)=1
                 else
                    lvlbegin(nn)=lvlend(nn-1)+1
                 endif
                 lvlend(nn)=lvlbegin(nn)+nlvlcore-1
                 if(k<=j) lvlend(nn)=lvlend(nn)+1
                 if(k==mynlvl3d) lvlend(nn)=min(lvlend(nn),nz)
              enddo
           elseif(ncfs_sfc%num_dim(i)==2) then
              nn=nn+1
              varname(nn)=trim(ncfs_sfc%list_varname(i))
              vartype(nn)=ncfs_sfc%vartype(i)
              nx(nn)=ncfs_sfc%dim_1(i)
              ny(nn)=ncfs_sfc%dim_2(i)
              lvlbegin(nn)=1
              lvlend(nn)=1
              fileid(nn)=3
           endif
        enddo
    endif

    do k=1,ntotalcore
       write(6,'(2I5,2x,a10,10I10)') k,fileid(k),varname(k),vartype(k),nx(k),ny(k),lvlbegin(k),lvlend(k) 
    enddo
!
!  save results
!
    this%fileid=fileid
    this%varname=varname
    this%vartype=vartype
    this%lvlbegin=lvlbegin
    this%lvlend=lvlend
    this%nx=nx
    this%ny=ny
!
!  release local arrary
!
    deallocate(fileid)
    deallocate(varname)
    deallocate(vartype)
    deallocate(lvlbegin)
    deallocate(lvlend)
    deallocate(nx)
    deallocate(ny)
    deallocate(nlvl3d_list)

  end subroutine arrange

end module module_mpi_arrange
