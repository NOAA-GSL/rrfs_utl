module module_ncfile_stat
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-03-08
!
! ABSTRACT: 
!     This module read fv3lam fields and figure out dimension of each fields
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

  public :: ncfile_stat

  type :: ncfile_stat

      character(len=80) :: filename
      integer :: numvar
      character(len=max_varname_length),allocatable :: list_varname(:)
      integer,allocatable :: num_dim(:)
      integer,allocatable :: dim_1(:)
      integer,allocatable :: dim_2(:)
      integer,allocatable :: dim_3(:)
      integer,allocatable :: vartype(:)
      integer :: num_totalvl

    contains
      procedure :: init
      procedure :: fill_dims
      procedure :: close
  end type ncfile_stat
!
! constants
!
contains

  subroutine init(this,filein,numvar,varlist)
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

    character(len=*), intent(in) :: filein
    integer, intent(in)          :: numvar
    character(len=*), intent(in) :: varlist
    class(ncfile_stat) :: this
!
    integer :: n
!
    this%filename=trim(filein)
    this%numvar=numvar
    if( numvar>0 ) then
       allocate(this%list_varname(this%numvar))
       allocate(this%dim_1(this%numvar))
       allocate(this%dim_2(this%numvar))
       allocate(this%dim_3(this%numvar))
       allocate(this%num_dim(this%numvar))
       allocate(this%vartype(this%numvar))
       this%dim_1=1
       this%dim_2=1
       this%dim_3=1
       this%num_dim=1
       this%vartype=0
       read(varlist,*) this%list_varname
    endif
    this%num_totalvl=0

    write(6,*) 'process file: ',trim(this%filename)
!    write(6,*) 'variable will be process in this file:'
!    write(6,*) this% list_varname
!
  end subroutine init

  subroutine close(this)
    implicit none
    class(ncfile_stat) :: this

    if(this%numvar>0) then
       deallocate(this%list_varname)
       deallocate(this%dim_1)
       deallocate(this%dim_2)
       deallocate(this%dim_3)
       deallocate(this%num_dim)
       deallocate(this%vartype)
    endif
    this%numvar=0
    this%num_totalvl=0

  end subroutine close

  subroutine fill_dims(this)

    use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
    use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
    use netcdf, only: nf90_inquire_variable,nf90_inq_varid

    implicit none
    class(ncfile_stat) :: this

    integer :: ncid
    integer :: ndimensions,nvariables,nattributes,unlimiteddimid
    integer,allocatable :: dim_len(:)
    integer,allocatable :: dim_id(:)
    integer :: k,len,ndim
    character(len=max_varname_length) :: name
    integer :: varid,iret,xtype
    logical :: ifexist
!
!
    if(this%numvar<=0) return
     
    inquire(file=trim(this%filename),exist=ifexist )
    if(.not.ifexist) then
      write(6,*) 'file does not exist ',trim(this%filename)
      this%num_totalvl=0
      return 
    endif
!    
    iret=nf90_open(this%filename,nf90_nowrite,ncid)

    iret=nf90_inquire(ncid,ndimensions,nvariables,nattributes,unlimiteddimid)
    if(allocated(dim_len)) deallocate(dim_len)
    allocate(dim_len(ndimensions))

    do k=1,ndimensions
       iret=nf90_inquire_dimension(ncid,k,name,len)
       dim_len(k)=len
    enddo

    do k=1,this%numvar
       name=this%list_varname(k)
       iret = nf90_inq_varid(ncid, trim(name), VarId)
       iret=nf90_inquire_variable(ncid,VarId,ndims=ndim)

       if(allocated(dim_id    )) deallocate(dim_id    )
       allocate(dim_id(ndim))
       iret = nf90_inquire_variable(ncid, VarId, dimids = dim_id(1:ndim))
       if(ndim==4) then
          if(dim_len(dim_id(4))==1) then
              this%dim_1(k)=dim_len(dim_id(1))
              this%dim_2(k)=dim_len(dim_id(2))
              this%dim_3(k)=dim_len(dim_id(3))
              this%num_dim(k)=3
          else
              write(6,*) '4th dimension must be 1',trim(name),dim_len(dim_id(4))
              stop 123
          endif
       elseif(ndim==3) then
          if(dim_len(dim_id(3))==1) then
              this%dim_1(k)=dim_len(dim_id(1))
              this%dim_2(k)=dim_len(dim_id(2))
              this%dim_3(k)=1
              this%num_dim(k)=2
          else
              this%dim_1(k)=dim_len(dim_id(1))
              this%dim_2(k)=dim_len(dim_id(2))
              this%dim_3(k)=dim_len(dim_id(3))
              this%num_dim(k)=3
          endif
       elseif(ndim==2) then
          this%dim_1(k)=dim_len(dim_id(1))
          this%dim_2(k)=dim_len(dim_id(2))
          this%dim_3(k)=1
          this%num_dim(k)=2
       else
          write(6,*) 'something wrong with dimension size, ',trim(name),ndim
          stop 123
       endif
          
       iret = nf90_inquire_variable(ncid, VarId, xtype = xtype)
       this%vartype(k)=xtype
    enddo

    iret=nf90_close(ncid)

    this%num_totalvl=0 
    do k=1,this%numvar
       this%num_totalvl=this%num_totalvl+this%dim_3(k)
       write(6,'(a20,10I10)') this%list_varname(k),this%num_dim(k),this%dim_1(k), &
                     this%dim_2(k),this%dim_3(k),this%vartype(k)
    enddo
    write(6,*) 'Total levels=',this%num_totalvl
!  if(xtype==NF90_FLOAT) then

  end subroutine fill_dims

end module module_ncfile_stat
