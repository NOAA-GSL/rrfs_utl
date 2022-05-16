subroutine ncio_ensmean_recenter(ensize,mype,new_comm,l_write_mean,l_recenter,varname,filename,filetail)
!
!---------------------------------------------------------------------- 
!  Purpose: Calculate ensemble mean file from input FV3LAM NETCDF input
!  ensemble members.
!
!  2021-04 Yongming Wang and X. Wang - Initial codes from WRFDA were changed for FV3LAM
!                              - Enable parallel ensemble IO
!                                poc: xuguang.wang@ou.edu
!  2022-05 Ming Hu -  Use a subroutine to read ensemble, calculate mean
!                      and recenter for one variable.
!
!----------------------------------------------------------------------

   use netcdf 
   implicit none

   integer, parameter    :: max_num_dims = 4          ! Maximum number of dimensions.

   integer,intent(in)    :: ensize                   ! size of ensemble
   integer,intent(in)    :: mype                     ! rank 
   integer,intent(in)    :: new_comm                 ! group communicator
   character (len=*),intent(in)   :: filename        ! General filename stub.
   character (len=*),intent(in)   :: filetail        ! file type
   character (len=*),intent(inout):: varname         ! Variable to search for.
   logical,intent(in)    :: l_write_mean             ! if write ensmeble mean
   logical,intent(in)    :: l_recenter               ! if recenter
!
!
!
   character (len=3)     :: ce                        ! Member index -> character.

   integer               :: nsize                     ! size of array
   integer               :: nsize2d                     ! size of array
   integer               :: k,i                       ! Loop counters.
   integer               :: length                    ! Filename length.
   integer               :: rcode                     ! NETCDF return code.
   integer               :: cdfid                     ! NETCDF file IDs.
   integer               :: cdfid_mean                ! NETCDF file IDs.
   integer               :: id_var                    ! NETCDF variable ID.
   integer               :: ivtype                    ! 4=integer, 5=real, 6=d.p.
   integer               :: ndims                     ! Number of field dimensions.
   integer               :: natts                     ! Number of field attributes.
   real                  :: rnanals                   ! 1 / ensemble size.

   integer               :: dimids(1:10)              ! Array of dimension IDs.
   integer               :: dims(1:max_num_dims)      ! Array of dimensions.
   integer               :: istart(1:max_num_dims)    ! Array of dimension starts.
   integer               :: iend(1:max_num_dims)      ! Array of dimension ends.

   real (kind=4), allocatable     :: data_r(:,:,:)       ! real array.
   real (kind=4), allocatable     :: data_r_diff(:,:,:)  ! real array mean.
   real (kind=4), allocatable     :: data_r_mean(:,:)  ! real array mean.
   real (kind=8), allocatable     :: data_d(:,:,:)       ! double 
   real (kind=8), allocatable     :: data_d_diff(:,:,:)  ! double mean      
   real (kind=8), allocatable     :: data_d_mean(:,:)  ! double mean      

   !
   character (len=200)   :: input_file                  ! General filename stub.
   integer :: iret
   logical :: l_positive
  
   
   include 'mpif.h'

   rnanals=1.0_8/ensize

!  Open file:
   if ( mype == 0 ) then
      if(l_recenter) then
         input_file = 'control_'//trim(filetail)
      else
         input_file = trim(filename)//'_'//trim(filetail)
      endif
   else
      write(UNIT=ce,FMT='(i3.3)') mype
      input_file = trim(filename)//'_mem'//trim(ce)//'_'//trim(filetail)
   endif
   if( mype <=1) print *, 'APM read ',mype,trim(input_file),' ', trim(varname)
   rcode = nf90_open( trim(input_file), NF90_NOWRITE, cdfid )

   rcode = nf90_inq_varid ( cdfid, trim(varname), id_var )
!  Check variable is in file:
   if ( rcode /= 0 ) then
      write(6,FMT='(A,A,I10)') &
            varname, ' variable is not in input file ',mype
   end if

!  Get dimension information for this field:
   dimids = 0
   dims = 0
   rcode = nf90_inquire_variable( cdfid, id_var, varname, ivtype, ndims, dimids, natts )
   do i = 1, ndims
      rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
   end do
   istart = 1
   iend = 1
   do i = 1, ndims-1
      iend(i) = dims(i)
   end do
!
!  Allocate and initialize data:
   if ( ivtype == 5 ) then
      allocate( data_r(iend(1),iend(2),iend(3)))
      allocate( data_r_diff(iend(1),iend(2),iend(3)))
      data_r = 0.0
      data_r_diff = 0.0
   elseif ( ivtype == 6 ) then
      allocate( data_d(iend(1),iend(2),iend(3)))
      allocate( data_d_diff(iend(1),iend(2),iend(3)))
      data_d = 0.0
      data_d_diff = 0.0
   else
      write(6,'(A,A)') varname, ' variable is not real/double type'
   end if

!  Calculate accumulating mean and variance:
   if(ivtype == 5) then
      call ncvgt( cdfid, id_var, istart, iend, data_r, rcode)
      if( varname == "ref_f3d" ) data_r = max(data_r,0.0)
      if( mype==0 ) then  ! save control 
         data_r_diff=data_r
         data_r=0.0
      endif
   elseif(ivtype == 6) then
       call ncvgt( cdfid, id_var, istart, iend, data_d, rcode)
       if( varname == "ref_f3d" ) data_d = max(data_d,0.0)
      if( mype==0 ) then  ! save control 
         data_d_diff=data_d
         data_d=0.0
      endif
   endif
   nsize=iend(1)*iend(2)*iend(3)
   nsize2d=iend(1)*iend(2)

   if(ivtype == 5) then
     if(mype==0) then
        allocate( data_r_mean(iend(1),iend(2)))
     endif
     do k=1,iend(3)
        call mpi_reduce(data_r(:,:,k),data_r_mean,nsize2d,mpi_real,mpi_sum,0,new_comm,iret)
        if(mype==0) then
           data_r(:,:,k)=data_r_mean
        endif
     enddo
     if(mype==0) then
        deallocate(data_r_mean)
     endif
   elseif(ivtype == 6) then
     if(mype==0) then
        allocate( data_d_mean(iend(1),iend(2)))
     endif
     do k=1,iend(3)
        call mpi_reduce(data_d(:,:,k),data_d_mean,nsize2d,mpi_double,mpi_sum,0,new_comm,iret)
        if(mype==0) then
           data_d(:,:,k)=data_d_mean
        endif
     enddo
     if(mype==0) then
        deallocate(data_d_mean)
     endif
   endif

   rcode = nf90_close( cdfid )

   if(mype==0) then
      if(ivtype == 5) then
        data_r = data_r * rnanals
        write(6,*) trim(varname),' ens mean=',maxval(data_r),minval(data_r)
      elseif(ivtype == 6) then
        data_d = data_d * rnanals
        write(6,*) trim(varname),' ens mean=',maxval(data_d),minval(data_d)
      endif

      if(l_write_mean) then
         input_file = trim(filename)//'_'//trim(filetail)
         print*,'mean file=',trim(input_file),' ', trim(varname)
         rcode = nf90_open(trim(input_file), NF90_WRITE, cdfid_mean )
         if ( rcode /= 0) then
            write(6,FMT='(A,A,1x,A,I10)') &
             ' Error opening netcdf file ', trim(input_file),trim(varname),rcode
         end if
         if(ivtype == 5) then
           call ncvpt( cdfid_mean, id_var, istart, iend, data_r, rcode)
         elseif(ivtype == 6) then
           call ncvpt( cdfid_mean, id_var, istart, iend, data_d, rcode)
         endif

         rcode = nf90_close( cdfid_mean )
      endif
   end if

   if(l_recenter) then
      l_positive=.false.
      if( trim(varname)=="ref_f3d" .or. trim(varname)=="sphum" .or. trim(varname)=="liq_wat" .or.   &
          trim(varname)=="rainwat" .or. trim(varname)=="ice_wat" .or. trim(varname)=="snowwat" .or. &
          trim(varname)=="graupel" )then
         l_positive=.true.
      endif
      if(mype==0) then
         if ( ivtype == 5 ) then
            if(l_positive) data_r=max(data_r, 0.0)
            data_r_diff=data_r_diff-data_r
            write(6,*) trim(varname),' difference =',maxval(data_r_diff),minval(data_r_diff)
         else
            if(l_positive) data_d=max(data_d, 0.0)
            data_d_diff=data_d_diff-data_d
            write(6,*) trim(varname),' difference =',maxval(data_d_diff),minval(data_d_diff)
         endif
      endif

      if(ivtype == 5) then
         call mpi_bcast(data_r_diff,nsize,mpi_real,0,new_comm,iret)
      elseif(ivtype == 6) then
         call mpi_bcast(data_d_diff,nsize,mpi_double,0,new_comm,iret)
      endif
!
      if(mype > 0) then
!
!  get the values after recenter
         if ( ivtype == 5 ) then
            data_r = data_r + data_r_diff
            if(l_positive) data_r=max(data_r, 0.0)
         elseif ( ivtype == 6 ) then
            data_d = data_d + data_d_diff
            if(l_positive) data_d=max(data_d, 0.0)
         endif
!
! update each member
!
!  Open file:
         write(UNIT=ce,FMT='(i3.3)') mype
         input_file ='rec_'//trim(filename)//'_mem'//trim(ce)//'_'//trim(filetail)
         if( mype <=1) print *, 'APM write ',trim(input_file), ' ',trim(varname)
         rcode = nf90_open( trim(input_file), NF90_WRITE, cdfid )
         if ( rcode /= 0 ) then
            write(UNIT=6,FMT='(A,A,I10)') trim(input_file), ' problem to open for writing ',mype
         end if

!  Get variable ID:
         rcode = nf90_inq_varid ( cdfid, varname, id_var )
!  Check variable is in file:
         if ( rcode /= 0 ) then
            write(UNIT=6,FMT='(A,A,I10)') trim(varname), ' variable is not in input file ',mype
         end if

!  Get dimension information for this field:
         dimids = 0
         dims = 0
         rcode = nf90_inquire_variable( cdfid, id_var, varname, ivtype, ndims, dimids, natts )
         do i = 1, ndims
             rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
         end do
         istart = 1
         iend = 1
         do i = 1, ndims-1
            iend(i) = dims(i)
         end do

         if ( ivtype == 5 ) then
            call ncvpt( cdfid, id_var, istart, iend, data_r, rcode)
         elseif ( ivtype == 6 ) then
            call ncvpt( cdfid, id_var, istart, iend, data_d, rcode)
         endif

         rcode = nf90_close( cdfid )

      end if
   end if

   if(ivtype == 5) then
      deallocate( data_r )
      deallocate( data_r_diff )
   elseif(ivtype == 6) then
      deallocate( data_d_diff )
      deallocate( data_d )
   endif

end subroutine ncio_ensmean_recenter

