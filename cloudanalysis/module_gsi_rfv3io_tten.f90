module gsi_rfv3io_tten_mod
!$$$   module documentation block
!             .      .    .                                       .
! module:     gsi_rfv3io_mod
!   prgmmr:
!
! abstract: IO routines for regional FV3
!
  use kinds, only: r_kind,i_kind,r_single
  use netcdf, only: nf90_write,nf90_inq_varid
  use netcdf, only: nf90_put_var,nf90_get_var
  use netcdf, only: nf90_open,nf90_close,nf90_noerr
  use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
  use netcdf, only: nf90_inquire_variable

  implicit none
  integer(i_kind) rfv3io_mype
  integer(i_kind) nx,ny,nz 
  integer(i_kind) lat2,lon2,nsig
  integer(i_kind) nlon_regional,nlat_regional,nsig_regional
 
  real,parameter :: grav=9.8

!    directory names (hardwired for now)
  type type_fv3regfilenameg
      character(len=:),allocatable :: grid_spec     !='fv3_grid_spec'
      character(len=:),allocatable :: ak_bk         !='fv3_akbk'
      character(len=:),allocatable :: dynvars       !='fv3_dynvars'
      character(len=:),allocatable :: tracers       !='fv3_tracer'
      character(len=:),allocatable :: tracers_unc   !='fv3_tracer_unc'
      character(len=:),allocatable :: sfcdata       !='fv3_sfcdata'
      character(len=:),allocatable :: phydata       !='fv3_phydata'
      character(len=:),allocatable :: couplerres    !='coupler.res'
      contains
      procedure , pass(this):: init=>fv3regfilename_init
  end type type_fv3regfilenameg

  integer(i_kind):: fv3sar_bg_opt=0
  type(type_fv3regfilenameg):: bg_fv3regfilenameg


! set default to private
  private
! set subroutines to public
  public :: nlon_regional,nlat_regional,nsig_regional
  public :: gsi_rfv3io_get_grid_specs
  public :: type_fv3regfilenameg
  public :: bg_fv3regfilenameg
  public :: fv3sar_bg_opt
  public :: rfv3io_mype


  public :: gsi_fv3ncdf_read
  public :: gsi_fv3ncdf2d_read
  public :: gsi_fv3ncdf_write

  public :: aeta1_ll,aeta2_ll
  public :: eta1_ll,eta2_ll

  real(r_kind),allocatable :: aeta1_ll(:),aeta2_ll(:)
  real(r_kind),allocatable :: eta1_ll(:),eta2_ll(:)

contains

  subroutine fv3regfilename_init(this,grid_spec_input,ak_bk_input,dynvars_input,&
                      tracers_input,tracers_unc_input,sfcdata_input,phydata_input,couplerres_input)
  implicit None
  class(type_fv3regfilenameg),intent(inout):: this

  character(*),optional :: grid_spec_input,ak_bk_input,dynvars_input, &
                      tracers_input,tracers_unc_input,sfcdata_input,phydata_input,couplerres_input

  if(present(grid_spec_input))then
    this%grid_spec=grid_spec_input
  else
    this%grid_spec='fv3_grid_spec'
  endif

  if(present(ak_bk_input))then
    this%ak_bk=ak_bk_input
  else
    this%ak_bk='fv3_akbk'
  endif

  if(present(dynvars_input))then
    this%dynvars=dynvars_input
  else
    this%dynvars='fv3_dynvars'
  endif

  if(present(tracers_input))then
    this%tracers=tracers_input
  else
    this%tracers='fv3_tracer'
  endif

  if(present(tracers_unc_input))then
    this%tracers_unc=tracers_unc_input
  else
    this%tracers_unc='fv3_tracer_unc'
  endif

  if(present(sfcdata_input))then
    this%sfcdata=sfcdata_input
  else
    this%sfcdata='fv3_sfcdata'
  endif

  if(present(phydata_input))then
    this%phydata=phydata_input
  else
    this%phydata='fv3_phydata'
  endif

  if(present(couplerres_input))then
    this%couplerres=couplerres_input
  else
    this%couplerres='coupler.res'
  endif

  end subroutine fv3regfilename_init

subroutine gsi_rfv3io_get_grid_specs(fv3filenamegin,ierr)


  implicit none
  integer(i_kind),intent(  out) :: ierr
  type (type_fv3regfilenameg),intent(in) :: fv3filenamegin

  integer(i_kind) gfile_grid_spec

  character(len=:),allocatable    :: grid_spec
  character(len=:),allocatable    :: ak_bk
  integer(i_kind) i,k,ndimensions,iret,nvariables,nattributes,unlimiteddimid
  integer(i_kind) len,gfile_loc
  character(len=128) :: name

  real(r_kind),allocatable :: ak(:),bk(:),abk_fv3(:)

  grid_spec=fv3filenamegin%grid_spec
  ak_bk=fv3filenamegin%ak_bk
!!!!!!!!!!    grid_spec  !!!!!!!!!!!!!!!
    ierr=0
    iret=nf90_open(trim(grid_spec),nf90_nowrite,gfile_grid_spec)
    if(iret/=nf90_noerr) then
       write(6,*)' gsi_rfv3io_get_grid_specs: problem opening', &
                   trim(grid_spec),', Status = ',iret
       ierr=1
       return
    endif

    iret=nf90_inquire(gfile_grid_spec,ndimensions,nvariables,nattributes,unlimiteddimid)
    gfile_loc=gfile_grid_spec
    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       if(trim(name)=='grid_xt') nx=len
       if(trim(name)=='grid_yt') ny=len
    enddo
    nlon_regional=nx
    nlat_regional=ny
    lat2=ny
    lon2=nx
    if(rfv3io_mype==0)write(6,*),'nx,ny=',nx,ny

!!!    get nx,ny,grid_lon,grid_lont,grid_lat,grid_latt,nz,ak,bk

    iret=nf90_open(trim(ak_bk),nf90_nowrite,gfile_loc)
    if(iret/=nf90_noerr) then
       write(6,*)'gsi_rfv3io_get_grid_specs: problem opening ', &
                 trim(ak_bk),',Status = ',iret
       ierr=1
       return
    endif
    iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
    do k=1,ndimensions
       iret=nf90_inquire_dimension(gfile_loc,k,name,len)
       if(trim(name)=='xaxis_1') nz=len
    enddo
    if(rfv3io_mype==0)write(6,'(" nz=",i5)') nz

    nsig=nz-1
    nsig_regional=nz-1

!!!    get ak,bk

    allocate(aeta1_ll(nsig),aeta2_ll(nsig))
    allocate(eta1_ll(nsig+1),eta2_ll(nsig+1))
    allocate(ak(nz),bk(nz),abk_fv3(nz))


    do k=ndimensions+1,nvariables
       iret=nf90_inquire_variable(gfile_loc,k,name,len)
       if(trim(name)=='ak'.or.trim(name)=='AK') then
          iret=nf90_get_var(gfile_loc,k,abk_fv3)
          do i=1,nz
             ak(i)=abk_fv3(nz+1-i)
          enddo
       endif
       if(trim(name)=='bk'.or.trim(name)=='BK') then
          iret=nf90_get_var(gfile_loc,k,abk_fv3)
          do i=1,nz
             bk(i)=abk_fv3(nz+1-i)
          enddo
       endif
    enddo
    iret=nf90_close(gfile_loc)

!!!!! change unit of ak 
    do i=1,nsig+1
       eta1_ll(i)=ak(i)*0.001_r_kind
       eta2_ll(i)=bk(i)
    enddo
    do i=1,nsig
       aeta1_ll(i)=0.5_r_kind*(ak(i)+ak(i+1))*0.001_r_kind
       aeta2_ll(i)=0.5_r_kind*(bk(i)+bk(i+1))
    enddo
    if(rfv3io_mype==0)then
       do i=1,nz
          write(6,'(" ak,bk(",i3,") = ",2f17.6)') i,ak(i),bk(i)
       enddo
    endif


    deallocate (ak,bk,abk_fv3)

    return
end subroutine gsi_rfv3io_get_grid_specs

subroutine gsi_fv3ncdf2d_read(filenamein,varname,varname2,work_sub,mype_io)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsi_fv3ncdf2d_read       
!   prgmmr: wu w             org: np22                date: 2017-10-17
!
! abstract: read in 2d fields from fv3_sfcdata file in mype_2d 
!                Scatter the field to each PE 
! program history log:
!   input argument list:
!     it    - time index for 2d fields
!
!   output argument list:
!     ges_z - surface elevation
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block

    implicit none

    character(*)     ,intent(in   ) :: varname,varname2,filenamein
    real(r_single)   ,intent(out  ) :: work_sub(lon2,lat2) 
    integer(i_kind)  ,intent(in   ) :: mype_io

    character(len=128) :: name
    real(r_kind),allocatable,dimension(:,:)   :: f2d
    integer(i_kind),allocatable,dimension(:)  :: dim_id,dim

    integer(i_kind) iret,gfile_loc,i,k,len,ndim
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    if(rfv3io_mype == mype_io) then

       work_sub=-99999.0

       iret=nf90_open(filenamein,nf90_nowrite,gfile_loc)
       if(iret/=nf90_noerr) then
          write(6,*)' problem opening3 ',trim(filenamein),', Status = ',iret
          return
       endif

       iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
       allocate(dim(ndimensions))
       do k=1,ndimensions
          iret=nf90_inquire_dimension(gfile_loc,k,name,len)
          dim(k)=len
       enddo
       !write(*,*) ndimensions,nvariables,nattributes,unlimiteddimid

       do i=1,nvariables
          iret=nf90_inquire_variable(gfile_loc,i,name,len)
          if( trim(name)==varname.or.trim(name)==varname2 ) then
             iret=nf90_inquire_variable(gfile_loc,i,ndims=ndim)
             if(allocated(dim_id    )) deallocate(dim_id    )
             allocate(dim_id(ndim))
             iret=nf90_inquire_variable(gfile_loc,i,dimids=dim_id)
             !write(*,*) trim(name),ndim,dim_id

             if(allocated(f2d       )) deallocate(f2d       )
             if(ndim==2) then
                allocate(f2d(dim(dim_id(1)),dim(dim_id(2))))
             elseif(ndim==3 .and. dim(dim_id(3))==1) then
                allocate(f2d(dim(dim_id(1)),dim(dim_id(2))))
             else
                write(6,*) 'unknow dimension szie=',ndim
             endif

             iret=nf90_get_var(gfile_loc,i,f2d)

             if(lon2==dim(dim_id(1)) .and. lat2==dim(dim_id(2))) then
                work_sub(:,:)=f2d(:,:)
             else
                write(6,*) 'horizontal mismatch=',lon2,lat2
                write(6,*) 'dimension reads in=',dim(dim_id(1)),dim(dim_id(2))
                stop 123
             endif
             deallocate (f2d)

             deallocate (dim_id,dim)
             exit
          endif
       enddo ! i
       iret=nf90_close(gfile_loc)

    endif  ! mype

    return
end subroutine gsi_fv3ncdf2d_read

subroutine gsi_fv3ncdf_read(filenamein,varname,varname2,work_sub,mype_io)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsi_fv3ncdf_read       
!   prgmmr: wu               org: np22                date: 2017-10-10
!
! abstract: read in a field from a netcdf FV3 file in mype_io
!          then scatter the field to each PE 
! program history log:
!
!   input argument list:
!     filename    - file name to read from       
!     varname     - variable name to read in
!     varname2    - variable name to read in
!     mype_io     - pe to read in the field
!
!   output argument list:
!     work_sub    - output sub domain field
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$  end documentation block


    implicit none
    character(*)      ,intent(in   ) :: varname,varname2,filenamein
    real(r_single)    ,intent(out  ) :: work_sub(lon2,lat2,nsig) 
    integer(i_kind)   ,intent(in   ) :: mype_io

    character(len=128)                          :: name
    real(r_kind),allocatable,dimension(:,:,:)   :: f3d
    integer(i_kind),allocatable,dimension(:)    :: dim_id,dim

    integer(i_kind) n,ns,k,len,ndim
    integer(i_kind) gfile_loc,iret
    integer(i_kind) kk
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    if(rfv3io_mype==mype_io ) then

       work_sub=-99999.0

       iret=nf90_open(trim(filenamein),nf90_nowrite,gfile_loc)
       if(iret/=nf90_noerr) then
          write(6,*)' gsi_fv3ncdf_read: problem opening ',trim(filenamein),gfile_loc,', Status = ',iret
          write(6,*)' gsi_fv3ncdf_read:problem opening5 with varnam ',trim(varname)
          return
       endif

       iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
       allocate(dim(ndimensions))

       do k=1,ndimensions
          iret=nf90_inquire_dimension(gfile_loc,k,name,len)
          dim(k)=len
       enddo

       do k=1,nvariables
          iret=nf90_inquire_variable(gfile_loc,k,name,len)
          !write(*,*) k,trim(name),varname,varname2
          if(trim(name)==varname .or. trim(name)==varname2) then
             iret=nf90_inquire_variable(gfile_loc,k,ndims=ndim)
             if(allocated(dim_id    )) deallocate(dim_id    )
             allocate(dim_id(ndim))
             iret=nf90_inquire_variable(gfile_loc,k,dimids=dim_id)
             !write(*,*) trim(name),ndim,dim_id

             if(ndim==3) then
                if(allocated(f3d        )) deallocate(f3d        )
                allocate(f3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
             elseif(ndim==4 .and. dim(dim_id(4))==1) then
                if(allocated(f3d        )) deallocate(f3d        )
                allocate(f3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
             else
                write(6,*) 'unknow dimension szie=',ndim
             endif

             iret=nf90_get_var(gfile_loc,k,f3d)
             !do kk=1,dim(dim_id(3))
             !   write(6,*) kk,maxval(f3d(:,:,kk)),minval(f3d(:,:,kk))
             !enddo

             if(lon2==dim(dim_id(1)) .and. lat2==dim(dim_id(2))) then
                if(nsig == dim(dim_id(3))) then
                   do kk=1,nsig
                     work_sub(:,:,kk)=f3d(:,:,nsig-kk+1)
                   enddo
                elseif(nsig+1== dim(dim_id(3))) then
                   do kk=1,nsig
                     work_sub(:,:,kk)=f3d(:,:,nsig+1-kk+1)
                   enddo
                elseif(nsig+2==dim(dim_id(3)) .and. trim(name)=='zh') then
                   do kk=1,nsig
                     work_sub(:,:,kk)=( f3d(:,:,nsig+2-kk+1)+ &
                                        f3d(:,:,nsig+2-kk) )*0.5_r_kind
                   enddo
                else
                   write(6,*) 'vertical mismatch=',nsig,nsig+1
                   write(6,*) 'dimension reads in=',dim(dim_id(3))
                   stop 123
                endif
             else
                write(6,*) 'horizontal mismatch=',lon2,lat2
                write(6,*) 'dimension reads in=',dim(dim_id(1)),dim(dim_id(2))
                stop 123
             endif

             deallocate (f3d,dim,dim_id)
             exit
          endif
       enddo     !   k
       iret=nf90_close(gfile_loc)

    endif !mype

    return
end subroutine gsi_fv3ncdf_read

subroutine gsi_fv3ncdf_write(filename,varname,var,mype_io)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    gsi_nemsio_write
!   pgrmmr: wu
!
! abstract:
!
! program history log:
!
!   input argument list:
!    varu,varv
!    add_saved
!    mype     - mpi task id
!    mype_io
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    implicit none

    real(r_single)   ,intent(in   ) :: var(lon2,lat2,nsig)
    integer(i_kind),intent(in   ) :: mype_io
    character(*)   ,intent(in   ) :: varname,filename

    integer(i_kind) :: VarId,gfile_loc
    integer(i_kind),allocatable,dimension(:)    :: dim_id,dim
    character(len=128)                          :: name
    integer(i_kind) iret,k,kk,len,ndim
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    real(r_kind),allocatable,dimension(:,:,:)   :: f3d

    if(rfv3io_mype==mype_io) then

       call check( nf90_open(trim(filename),nf90_write,gfile_loc) )

       iret=nf90_inquire(gfile_loc,ndimensions,nvariables,nattributes,unlimiteddimid)
       allocate(dim(ndimensions))

       do k=1,ndimensions
          iret=nf90_inquire_dimension(gfile_loc,k,name,len)
          dim(k)=len
       enddo

       call check( nf90_inq_varid(gfile_loc,trim(varname),VarId) )

       iret=nf90_inquire_variable(gfile_loc,VarId,ndims=ndim)
       if(allocated(dim_id    )) deallocate(dim_id    )
       allocate(dim_id(ndim))
       iret=nf90_inquire_variable(gfile_loc,VarId,dimids=dim_id)

       if(ndim==3 .or. (ndim==4 .and. dim(dim_id(4))==1) ) then
          if(allocated(f3d        )) deallocate(f3d        )
          allocate(f3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

          if(lon2==dim(dim_id(1)) .and. lat2==dim(dim_id(2))) then
             if(nsig == dim(dim_id(3))) then
                do kk=1,nsig
                   f3d(:,:,nsig-kk+1)=var(:,:,kk)
                enddo
             elseif(nsig+1== dim(dim_id(3))) then
                do kk=1,nsig
                  f3d(:,:,nsig+1-kk+1)=var(:,:,kk)
                enddo
                f3d(:,:,1)=var(:,:,nsig)
             else
                write(6,*) 'vertical mismatch=',nsig,nsig+1
                write(6,*) 'dimension reads in=',dim(dim_id(3))
                stop 123
             endif
          else
             write(6,*) 'horizontal mismatch=',lon2,lat2
             write(6,*) 'dimension reads in=',dim(dim_id(1)),dim(dim_id(2))
             stop 123
          endif
       else
          write(6,*) 'unknow dimension szie=',ndim
       endif

       write(6,*) 'write out ',trim(varname),' to ',trim(filename)
       call check( nf90_put_var(gfile_loc,VarId,f3d) )
       call check( nf90_close(gfile_loc) )
    end if !mype_io

end subroutine gsi_fv3ncdf_write

subroutine check(status)
    use kinds, only: i_kind
    use netcdf, only: nf90_noerr,nf90_strerror
    integer(i_kind), intent ( in) :: status

    if(status /= nf90_noerr) then
       print *,'ncdf error ', trim(nf90_strerror(status))
       stop  
    end if
end subroutine check

end module gsi_rfv3io_tten_mod
