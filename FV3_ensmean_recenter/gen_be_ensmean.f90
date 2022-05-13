program gen_be_ensmean
!
!---------------------------------------------------------------------- 
!  Purpose: Calculate ensemble mean file from input FV3LAM NETCDF input
!  ensemble members.
!
!  2021-04 Yongming Wang and X. Wang - Initial codes from WRFDA were changed for FV3LAM
!                              - Enable parallel ensemble IO
!                                poc: xuguang.wang@ou.edu
!
!  2022-05 Ming Hu -  Add function to process multiple variables at one run
!----------------------------------------------------------------------

   use netcdf 
   implicit none

   integer, parameter    :: max_num_file = 5
   integer, parameter    :: filename_len=200
!
! namelist
  integer :: fv3_io_layout_y
  character(len=filename_len)  :: filebase
  character(len=filename_len)  :: filetail(max_num_file)
  integer                      :: numvar(max_num_file)
  character(len=filename_len)  :: varlist(max_num_file)
  integer :: ens_size
  logical    :: l_write_mean             ! if write ensmeble mean
  logical    :: l_recenter               ! if recenter

  namelist/setup/ ens_size,fv3_io_layout_y,l_write_mean,l_recenter, &
                  filebase,filetail,&
                  numvar,varlist

   character (len=filename_len)   :: directory                 ! General filename stub.
   character (len=filename_len)   :: filename                  ! General filename stub.

   integer :: totalnumvar
   character (len=20), allocatable  ::  tailist_all(:)
   character (len=20), allocatable  ::  varlist_all(:)
   character (len=20)    :: varname                       ! Variable to search for.

   ! === variables for mpi 
   integer                :: iret, mype, npe, mype1, orig_group, new_group, new_comm
   integer,allocatable,dimension(:) :: new_group_members
   integer :: new_rank,new_size
   integer,allocatable,dimension(:,:) :: dis_group
   integer               :: groupsize
   integer               :: numgroup
   integer               :: num_cores
   integer               :: num_iteration
   integer               :: color
   integer               :: num_files
!
!
   integer               :: i,j,k
   logical               :: ifexist

   !
  
  include 'mpif.h'

! Initialize mpi, mype is process number, npe is total number of processes.
  call mpi_init(iret)
  call mpi_comm_rank(mpi_comm_world,mype,iret)
  call mpi_comm_size(mpi_comm_world,npe,iret)

  mype1 = mype + 1

  if (mype == 0) print*, "Calculate ensemble mean for FV3LAM"

  ! Get user input from command line

  directory='./'

!
!  get namelist
!
  numvar=0
  varlist=''
  filetail=''
  ens_size=1
  fv3_io_layout_y=1
  filebase='fv3sar_tile1'
  l_write_mean=.true.
  l_recenter=.false.


  inquire(file='namelist.ens', EXIST=ifexist )
  if(ifexist) then
    open(10,file='namelist.ens',status='old')
       read(10,setup)
    close(10)
  else
     write(*,*) 'No namelist file exist, use default values'
  endif

  if(mype==0) then
    write(*,*) 'Namelist setup are:'
    write(*,setup)
  endif

  !filename = trim(adjustl(directory)) // trim(adjustl(filebase))
  filename = trim(adjustl(filebase))
!
! find how many variables to process
!
  totalnumvar=0
  do k=1,max_num_file
     totalnumvar=totalnumvar+numvar(k)
  enddo
  if(totalnumvar <=0) then
     write(6,'(a)')'***ERROR***  varaible number is 0'
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  endif
  allocate(tailist_all(totalnumvar))
  allocate(varlist_all(totalnumvar))

  i=1
  num_files=0
  num_iteration=0
  do k=1,max_num_file
     if(numvar(k) > 0) then
        tailist_all(i:i+numvar(k)-1)=trim(filetail(k))
        read(varlist(k),*) varlist_all(i:i+numvar(k)-1)
        i=i+numvar(k)
        num_files=num_files+1
        if(num_iteration < numvar(k)) num_iteration = numvar(k)
     endif
  enddo 
  if(mype==0) then
     write(6,*) 'total variable number=',totalnumvar,' file number=',num_files
     write(6,'(50a10)') (trim(tailist_all(i)),i=1,totalnumvar)
     write(6,'(50a10)') (trim(varlist_all(i)),i=1,totalnumvar)
  endif

!
!  find group size and number
!
  groupsize=ens_size+1
  numgroup=npe/groupsize
  if(mype==0) write(6,*) 'number of groups=',numgroup,' group size=',groupsize
  if ( groupsize < 1 ) then
     write(6,'(a,i4)')'***ERROR***  npe too small.  npe = ',npe
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  end if
  if( numgroup /= num_files) then
     write(6,'(a,i4)')'***ERROR***  numgroup does not equal to num_files = '
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  end if

! find how many cores will be used
  num_cores=groupsize*numgroup

! find how many iterations needed to cover all the varaibles
!  num_iteration=totalnumvar/numgroup+1

! how to distribute variables to group and iteration
  allocate(dis_group(numgroup,num_iteration))
  dis_group=0
!  j=0
!  do k=1,num_iteration
!     do i=1,numgroup
!        j=j+1
!        if(j <=totalnumvar) then
!          dis_group(i,k)=j
!        endif
!     enddo
!  enddo
   k=1
   i=1
   dis_group(1,1)=1
   do j=2,totalnumvar
     if(trim(tailist_all(j)) == trim(tailist_all(j-1))) then
       k=k+1
       dis_group(i,k)=j
     else
       i=i+1
       k=1
       dis_group(i,k)=j
     endif
   enddo

! Create sub-communicator to handle number of cases (nanals)
  allocate(new_group_members(num_cores))
  do k=1,num_cores
     new_group_members(k)=k-1
  end do

  if(mype < num_cores) then
     color = mype/groupsize+1
  else
     color = MPI_UNDEFINED
  endif

  call MPI_Comm_split(mpi_comm_world,color,new_group_members,new_comm,iret)
!  call mpi_comm_group(mpi_comm_world,orig_group,iret)
!  call mpi_group_incl(orig_group,groupsize,new_group_members,new_group,iret)
!  call mpi_comm_create(mpi_comm_world,new_group,new_comm,iret)
  if ( iret /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',iret
     call mpi_abort(mpi_comm_world,101,iret)
  endif

  ! Process input files (one file per task)
  if (MPI_COMM_NULL /= new_comm) then
     call MPI_Comm_size(new_comm, new_size, iret)
     call MPI_Comm_rank(new_comm, new_rank, iret)
     do k=1,num_iteration
        j=dis_group(color,k)
        if(j>=1 .and. j<=totalnumvar) then
           call ncio_ensmean_recenter(ens_size,new_rank,new_comm,l_write_mean,l_recenter,&
                                      varlist_all(j),filename,tailist_all(j))
           call mpi_barrier(new_comm,iret)
        endif
     enddo
  else
    write(6,'(a,i5)') 'No files to process for mpi task = ',mype
  end if
  call mpi_barrier(mpi_comm_world,iret)

  deallocate(new_group_members)

!  call MPI_Group_free(orig_group,iret)
  if (MPI_COMM_NULL /= new_comm) then
!     call MPI_Group_free(new_group,iret)
     call MPI_Comm_free(new_comm,iret)
  endif

  if ( mype == 0 ) write(6,'(a)')'ensmean_done'

  deallocate(dis_group)
  deallocate(tailist_all)
  deallocate(varlist_all)

  call mpi_finalize(iret)

end program gen_be_ensmean

