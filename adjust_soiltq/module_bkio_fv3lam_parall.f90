module module_bkio_fv3lam_parall
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-03-08
!
! ABSTRACT: 
!     This module read and write fv3lam fields for soil ajustment
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

  use kinds, only: r_kind,r_single
  implicit none
!
! Rset default to private
!
  
  private

  public :: bkio_fv3lam

  type :: bkio_fv3lam

      character(len=80)   :: sfcfile
      character(len=80)   :: gridfile
      character(len=80)   :: tracerfile
      character(len=80)   :: dynfile

      integer             :: iyear,imonth,iday,ihour,iminute
      integer             :: nlon
      integer             :: nlat
      integer             :: nsoil
      real(r_kind)        :: p_top
      integer             :: fv3_io_layout_y
      integer,allocatable :: fv3_layout_begin(:)
      integer,allocatable :: fv3_layout_end(:)

! ges_xlon: longitude in degree
! ges_xlat: latitudde in degree
! sno     : snow depth (snodl)
! sncovr  : snow coverage (sncovr)
! landmask: land mask (0,1,2)
! ges_p1  : lowest model level Pressure (cp caculcted from delp using 2mb top)
! ges_t1  : lowest model level T (K: t)
! ges_q1  : lowest model level Q (kg/kg: specific humidity: from sphum)
! sumqc   : maximum qi and qc in a column  (kg/kg: mixing ratio?: from liq_wat,ice_wat)
! qsatg   : lowest model level Q saturation (kg/kg, specific humidity: calculated)
! ges_tsk    :  surface temperature (K: from tsfcl)
! ges_soilt1 :  top level soil temperature (K: from tsnow_land)
! ges_qvg :  surface water vapor mixing ratio (kg/kg: from qwv_surf_land)
! ges_qcg :  surface cloud water mixing ratio (kg/kg: from clw_surf_land)
! ges_tslb   :  soil temperature (K: from tslb)
! ges_smois  :  soil moisture (mixing ratio: from smois)


      real(r_kind),allocatable,dimension(:,:)  :: ges_xlon
      real(r_kind),allocatable,dimension(:,:)  :: ges_xlat
      real(r_kind),allocatable,dimension(:,:)  :: coast_prox
      real(r_kind),allocatable,dimension(:,:)  :: sno
      real(r_kind),allocatable,dimension(:,:)  :: sncovr
      !real(r_kind),allocatable,dimension(:,:)  :: landmask
      real(r_kind),dimension(:,:),pointer  :: landmask

      real(r_kind),dimension(:,:),pointer  :: ges_p1
      real(r_kind),dimension(:,:),pointer  :: ges_t1
      real(r_kind),dimension(:,:),pointer  :: ges_q1
      real(r_kind),allocatable,dimension(:,:)  :: sumqc
      real(r_kind),allocatable,dimension(:,:)  :: qsatg

      real(r_kind),dimension(:,:),pointer  :: ges_tsk   
      real(r_kind),dimension(:,:),pointer  :: tsk_comp   
      real(r_kind),dimension(:,:),pointer  :: ges_soilt1
      real(r_kind),dimension(:,:),pointer  :: ges_qvg
      real(r_kind),dimension(:,:),pointer  :: ges_qcg
      real(r_kind),dimension(:,:,:),pointer:: ges_tslb
      real(r_kind),dimension(:,:,:),pointer:: ges_smois

      integer :: is_t,is_q
      real(r_kind),dimension(:,:  ),pointer :: tinc      =>NULL()
      real(r_kind),dimension(:,:  ),pointer :: qinc      =>NULL()

    contains
      procedure :: init
      procedure :: setup_grid
      procedure :: read_ges
      procedure :: update_soil
      procedure :: close
  end type bkio_fv3lam
!
! constants
!
contains

  subroutine init(this,fv3_io_layout_y,iyear,imonth,iday,ihour,iminute)
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

    integer,intent(in) :: fv3_io_layout_y
    integer,intent(in) :: iyear,imonth,iday,ihour,iminute

    class(bkio_fv3lam) :: this

!
    this%fv3_io_layout_y=fv3_io_layout_y
    this%iyear=iyear
    this%imonth=imonth
    this%iday=iday
    this%ihour=ihour
    this%iminute=iminute
!
    this%sfcfile='sfc_data.nc'
    this%gridfile='fv3_grid_spec'
    this%tracerfile='fv_tracer.res.tile1.nc'
    this%dynfile='fv_core.res.tile1.nc'
!
    this%is_t=1
    this%is_q=1

  end subroutine init

  subroutine close(this)
    implicit none
    class(bkio_fv3lam) :: this

    this%fv3_io_layout_y=0
    if(allocated(this%ges_xlon)) deallocate(this%ges_xlon)
    if(allocated(this%ges_xlat)) deallocate(this%ges_xlat)
    if(allocated(this%fv3_layout_begin)) deallocate(this%fv3_layout_begin)
    if(allocated(this%fv3_layout_end)) deallocate(this%fv3_layout_end)

    if(allocated(this%coast_prox)) deallocate(this%coast_prox)
    if(allocated(this%sno)) deallocate(this%sno)
    if(allocated(this%sncovr)) deallocate(this%sncovr)
    if(associated(this%landmask)) deallocate(this%landmask)

    if(associated(this%ges_p1)) deallocate(this%ges_p1)
    if(associated(this%ges_t1)) deallocate(this%ges_t1)
    if(associated(this%ges_q1)) deallocate(this%ges_q1)
    if(allocated(this%sumqc)) deallocate(this%sumqc)
    if(allocated(this%qsatg)) deallocate(this%qsatg)

    if(associated(this%ges_tsk)) deallocate(this%ges_tsk)
    if(associated(this%tsk_comp)) deallocate(this%tsk_comp)
    if(associated(this%ges_soilt1)) deallocate(this%ges_soilt1)
    if(associated(this%ges_qvg)) deallocate(this%ges_qvg)
    if(associated(this%ges_qcg)) deallocate(this%ges_qcg)
    if(associated(this%ges_tslb)) deallocate(this%ges_tslb)
    if(associated(this%ges_smois)) deallocate(this%ges_smois)

    if(associated(this%tinc)) deallocate(this%tinc)
    if(associated(this%qinc)) deallocate(this%qinc)

  end subroutine close

  subroutine update_soil(this)
    use module_ncio, only: ncio
    implicit none
    class(bkio_fv3lam) :: this
    type(ncio) :: fv3io

    character(len=80) :: thisfv3file
    integer :: id,fv3_io_layout_y
    integer :: nlon,nlat,nlat_local,nz
    integer :: i,j,k
    real(r_kind),allocatable :: r2d8b(:,:)
    real(r_kind),allocatable :: r3d8b(:,:,:)

    nlon=this%nlon
    nlat=this%nlat
    nz=this%nsoil
    fv3_io_layout_y=this%fv3_io_layout_y

    write(6,*) 'update soil fields========>'
    write(6,*) 'this%tsk=',maxval(this%ges_tsk),minval(this%ges_tsk)
    write(6,*) 'this%soilt1=',maxval(this%ges_soilt1),minval(this%ges_soilt1)
    write(6,*) 'this%qvg=',maxval(this%ges_qvg),minval(this%ges_qvg)
    do k=1,nz
       write(6,*) 'this%ges_tslb=',maxval(this%ges_tslb(:,:,k)),minval(this%ges_tslb(:,:,k))
    enddo
    do k=1,nz
       write(6,*) 'this%ges_smois=',maxval(this%ges_smois(:,:,k)),minval(this%ges_smois(:,:,k))
    enddo

!    this%ges_smois(:,:,9)=this%coast_prox
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%sfcfile),".",id-1
       else
          thisfv3file=trim(this%sfcfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

       call fv3io%open(trim(thisfv3file),'w',200)
! 
       allocate(r2d8b(nlon,nlat_local))
       r2d8b(:,:)=this%ges_tsk(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))
       call fv3io%replace_var("tsfcl",nlon,nlat_local,r2d8b)

       r2d8b(:,:)=this%tsk_comp(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))
       call fv3io%replace_var("tsfc",nlon,nlat_local,r2d8b)

       r2d8b(:,:)=this%ges_soilt1(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))
       call fv3io%replace_var("tsnow_land",nlon,nlat_local,r2d8b)

       r2d8b(:,:)=this%ges_qvg(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))
       call fv3io%replace_var("qwv_surf_land",nlon,nlat_local,r2d8b)

       r2d8b(:,:)=this%ges_qcg(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))
       call fv3io%replace_var("clw_surf_land",nlon,nlat_local,r2d8b)

       deallocate(r2d8b)

       allocate(r3d8b(nlon,nlat_local,nz))
       r3d8b(:,:,:)=this%ges_tslb(:,this%fv3_layout_begin(id):this%fv3_layout_end(id),:)
       call fv3io%replace_var("tslb",nlon,nlat_local,nz,r3d8b)
       r3d8b(:,:,:)=this%ges_smois(:,this%fv3_layout_begin(id):this%fv3_layout_end(id),:)
       call fv3io%replace_var("smois",nlon,nlat_local,nz,r3d8b)

       deallocate(r3d8b)
       call fv3io%close
    enddo

  end subroutine update_soil

  subroutine setup_grid(this,mype)
    use module_ncio, only: ncio
    implicit none
    class(bkio_fv3lam) :: this
    integer,intent(in) :: mype
    type(ncio) :: fv3io
!
    character(len=80) :: thisfv3file
    integer :: id,iy
    integer :: nlon,nlat
    real,allocatable :: r2d4b(:,:)
!
!
    allocate(this%fv3_layout_begin(this%fv3_io_layout_y))
    allocate(this%fv3_layout_end(this%fv3_io_layout_y))
    this%fv3_layout_begin=0
    this%fv3_layout_end=0
!
    iy=0
    do id=1,this%fv3_io_layout_y
       if(this%fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%gridfile),".",id-1
       else
          thisfv3file=trim(this%gridfile)
       endif
       call fv3io%open(trim(thisfv3file),'r',0)
       call fv3io%get_dim("grid_xt",nlon)
       call fv3io%get_dim("grid_yt",nlat)
       if(mype==0) write(6,*) 'grid dimension =',id,nlon,nlat
       call fv3io%close
       this%fv3_layout_begin(id)=iy+1
       iy=iy+nlat
       this%fv3_layout_end(id)=iy
    enddo
    this%nlon=nlon
    this%nlat=this%fv3_layout_end(this%fv3_io_layout_y)
    if(mype==0) then
       write(6,'(a,2I10)') " nlon,nlat=",this%nlon,this%nlat
       write(6,'(a20,20I6)') "fv3_layout_begin=",this%fv3_layout_begin
       write(6,'(a20,20I6)') "fv3_layout_end=",this%fv3_layout_end
    endif
!
    if (mype==0) then
       allocate(this%ges_xlon(this%nlon,this%nlat))
       allocate(this%ges_xlat(this%nlon,this%nlat))
       do id=1,this%fv3_io_layout_y
          if(this%fv3_io_layout_y > 1) then
             write(thisfv3file,'(a,a,I4.4)') trim(this%gridfile),".",id-1
          else
             thisfv3file=trim(this%gridfile)
          endif
          nlon=this%nlon
          nlat=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1
          allocate(r2d4b(nlon,nlat))
          call fv3io%open(trim(thisfv3file),'r',0)
          call fv3io%get_var("grid_lont",nlon,nlat,r2d4b)
          this%ges_xlon(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d4b(:,:)
          call fv3io%get_var("grid_latt",nlon,nlat,r2d4b)
          this%ges_xlat(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d4b(:,:)
          call fv3io%close
          deallocate(r2d4b)
       enddo
       write(6,*) "lon=",maxval(this%ges_xlon),minval(this%ges_xlon)
       write(6,*) "lat=",maxval(this%ges_xlat),minval(this%ges_xlat)
    endif

  end subroutine setup_grid

  subroutine read_ges(this,mype)
!
    use module_ncio, only: ncio
    use mpi

    use netcdf, only: nf90_open,nf90_close,nf90_noerr
    use netcdf, only: nf90_get_var,nf90_put_var
    use netcdf, only: nf90_nowrite,nf90_write
    use netcdf, only: nf90_Inquire_Dimension
    use netcdf, only: nf90_inq_varid
    use netcdf, only: nf90_strerror

    implicit none
    class(bkio_fv3lam) :: this
    integer,intent(in) :: mype
    type(ncio) :: fv3io
!
    integer :: startloc(3)
    integer :: countloc(3)
    integer :: ncioid,var_id
    integer :: iret
!
    character(len=80) :: thisfv3file
    integer :: id,fv3_io_layout_y
    integer :: nlon,nlat,nlat_local,nz
    integer :: i,j,k
    real(r_kind),allocatable :: r2d8b(:,:)
    real(r_kind),allocatable :: r3d8b(:,:,:)
    real,allocatable :: r3d4b(:,:,:)
    real,allocatable :: r2d4b(:,:)
    real,allocatable :: r3d4b2(:,:,:)
    real,allocatable :: r1d4b(:)
!
!
!
    fv3_io_layout_y=this%fv3_io_layout_y   
    nlon=this%nlon
    nlat=this%fv3_layout_end(fv3_io_layout_y)
!
! read model top pressure
if(mype==2) then
    call fv3io%open('./fv_core.res.nc','r',200)
    call fv3io%get_dim("xaxis_1",nz)
    allocate(r1d4b(nz))
    call fv3io%get_var("ak",nz,r1d4b)
    this%p_top=r1d4b(1)
    call fv3io%close
    deallocate(r1d4b)
    write(6,*) 'model top pressure (pa)=',this%p_top
endif
!
    allocate(this%ges_p1(nlon,nlat))
    allocate(this%ges_t1(nlon,nlat))
    allocate(this%ges_q1(nlon,nlat))
    allocate(this%sumqc(nlon,nlat))
    allocate(this%qsatg(nlon,nlat))
    this%ges_p1=0.0
    this%ges_t1=0.0
    this%ges_q1=0.0
    this%sumqc=0.0
    this%qsatg=0.0
    allocate(this%tinc(nlon,nlat))
    allocate(this%qinc(nlon,nlat))
    this%tinc=0.0
    this%qinc=0.0

if(mype==1) then
! read q, qc and qi from tracer    
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%tracerfile),".",id-1
       else
          thisfv3file=trim(this%tracerfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

! use specific humidity 
       iret=nf90_open(trim(thisfv3file),nf90_nowrite,ncioid,info=MPI_INFO_NULL)
       if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(thisfv3file),', Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
       endif
       iret=nf90_inq_varid(ncioid,"zaxis_1",var_id)
       iret=nf90_Inquire_Dimension(ncioid, var_id, len = nz)

       allocate(r3d4b(nlon,nlat_local,1))
       startloc=(/1,1,nz/)
       countloc=(/nlon,nlat_local,1/)

       iret=nf90_inq_varid(ncioid,"sphum",var_id)
       iret=nf90_get_var(ncioid,var_id,r3d4b,start=startloc,count=countloc)
       this%ges_q1(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r3d4b(:,:,1)

       deallocate(r3d4b)
       iret=nf90_close(ncioid)
    enddo
    this%sumqc=0.0
    write(6,*) 'sumqc=',maxval(this%sumqc),minval(this%sumqc)
    write(6,*) 'ges_q1=',maxval(this%ges_q1),minval(this%ges_q1)
endif
!
! read t and pressure from dyn
!
if(mype==2) then
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%dynfile),".",id-1
       else
          thisfv3file=trim(this%dynfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

       call fv3io%open(trim(thisfv3file),'r',0)
       call fv3io%get_dim("zaxis_1",nz)
! Pressure
       allocate(r3d4b(nlon,nlat_local,nz+1))
       allocate(r2d8b(nlon,nlat_local))
       call fv3io%get_var("delp",nlon,nlat_local,nz,r3d4b(:,:,1:nz))
       r3d4b(:,:,nz+1)=this%p_top
       do k=nz,1,-1
          r3d4b(:,:,k)=r3d4b(:,:,k)+r3d4b(:,:,k+1)
       enddo
       r2d8b(:,:)=(r3d4b(:,:,1)+r3d4b(:,:,2))*0.5_8*0.001_8  ! cb
       this%ges_p1(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)
       deallocate(r3d4b)
       deallocate(r2d8b)
       call fv3io%close
! T
       iret=nf90_open(trim(thisfv3file),nf90_nowrite,ncioid,info=MPI_INFO_NULL)
       if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(thisfv3file),', Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
       endif

       iret=nf90_inq_varid(ncioid,"zaxis_1",var_id)
       iret=nf90_Inquire_Dimension(ncioid, var_id, len = nz)

       allocate(r3d4b(nlon,nlat_local,1))
       startloc=(/1,1,nz/)
       countloc=(/nlon,nlat_local,1/)

       iret=nf90_inq_varid(ncioid,"T",var_id)
       iret=nf90_get_var(ncioid,var_id,r3d4b,start=startloc,count=countloc)
       this%ges_t1(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r3d4b(:,:,1)

       deallocate(r3d4b)
       iret=nf90_close(ncioid)
    enddo
    write(6,*) 'ges_p1=',maxval(this%ges_p1),minval(this%ges_p1)
    write(6,*) 'ges_t1=',maxval(this%ges_t1),minval(this%ges_t1)
!
    call genqsat_2m(this%qsatg,this%ges_t1,this%ges_p1,nlon,nlat,1,.true.)
    write(6,*) 'qsat=',maxval(this%qsatg),minval(this%qsatg)
endif
!
! read q from tracer    
if(mype==3) then
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,a,I4.4)') "bk_",trim(this%tracerfile),".",id-1
       else
          thisfv3file="bk_"//trim(this%tracerfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

! use specific humidity 
       iret=nf90_open(trim(thisfv3file),nf90_nowrite,ncioid,info=MPI_INFO_NULL)
       if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(thisfv3file),', Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
       endif

       iret=nf90_inq_varid(ncioid,"zaxis_1",var_id)
       iret=nf90_Inquire_Dimension(ncioid, var_id, len = nz)

       allocate(r3d4b(nlon,nlat_local,1))
       startloc=(/1,1,nz/)
       countloc=(/nlon,nlat_local,1/)

       iret=nf90_inq_varid(ncioid,"sphum",var_id)
       iret=nf90_get_var(ncioid,var_id,r3d4b,start=startloc,count=countloc)
       this%qinc(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r3d4b(:,:,1)

       deallocate(r3d4b)
       iret=nf90_close(ncioid)
    enddo
    write(6,*) 'q bk=',maxval(this%qinc),minval(this%qinc)
endif
!
if(mype==4) then
    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,a,I4.4)') "bk_",trim(this%dynfile),".",id-1
       else
          thisfv3file="bk_"//trim(this%dynfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1
! T
       iret=nf90_open(trim(thisfv3file),nf90_nowrite,ncioid,info=MPI_INFO_NULL)
       if(iret/=nf90_noerr) then
           write(6,*)' problem opening ', trim(thisfv3file),', Status =',iret
           write(6,*)  nf90_strerror(iret)
           call flush(6)
           stop 333
       endif
       iret=nf90_inq_varid(ncioid,"zaxis_1",var_id)
       iret=nf90_Inquire_Dimension(ncioid, var_id, len = nz)

       allocate(r3d4b(nlon,nlat_local,1))
       startloc=(/1,1,nz/)
       countloc=(/nlon,nlat_local,1/)

       iret=nf90_inq_varid(ncioid,"T",var_id)
       iret=nf90_get_var(ncioid,var_id,r3d4b,start=startloc,count=countloc)
       this%tinc(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r3d4b(:,:,1)

       deallocate(r3d4b)
       iret=nf90_close(ncioid)
    enddo
    write(6,*) 't bk=',maxval(this%tinc),minval(this%tinc)
endif
!
!  read surface
!

if(mype==0) then
    allocate(this%coast_prox(nlon,nlat))
    allocate(this%sno(nlon,nlat))
    allocate(this%sncovr(nlon,nlat))
    allocate(this%landmask(nlon,nlat))

    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%sfcfile),".",id-1
       else
          thisfv3file=trim(this%sfcfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

       call fv3io%open(trim(thisfv3file),'r',200)
       call fv3io%get_dim("zaxis_1",nz)
! slmsk: 0 - water, 1 - land, 2 - ice
       allocate(r2d8b(nlon,nlat_local))
       call fv3io%get_var("slmsk",nlon,nlat_local,r2d8b)
       this%landmask(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%close
       deallocate(r2d8b)
    enddo
    this%nsoil=nz
    write(6,*) 'slmsk=',maxval(this%landmask),minval(this%landmask)

!    do j=1,nlat
!       do i=1,nlon
!          if( abs(this%landmask(i,j)) < 1.0e-5) then
!             this%qinc(i,j)=0.0_r_kind
!             this%tinc(i,j)=0.0_r_kind
!          endif
!       enddo
!    enddo

    call gsl_gen_coast_prox(nlon,nlat,this%landmask,this%coast_prox)
    write(6,*) 'coast_prox=',maxval(this%coast_prox),minval(this%coast_prox)

    allocate(this%ges_tsk(nlon,nlat))
    allocate(this%tsk_comp(nlon,nlat))
    allocate(this%ges_soilt1(nlon,nlat))
    allocate(this%ges_qvg(nlon,nlat))
    allocate(this%ges_qcg(nlon,nlat))

    allocate(this%ges_tslb(nlon,nlat,nz))
    allocate(this%ges_smois(nlon,nlat,nz))

    do id=1,fv3_io_layout_y
       if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,I4.4)') trim(this%sfcfile),".",id-1
       else
          thisfv3file=trim(this%sfcfile)
       endif
       nlat_local=this%fv3_layout_end(id)-this%fv3_layout_begin(id)+1

       call fv3io%open(trim(thisfv3file),'r',200)
! 
       allocate(r2d8b(nlon,nlat_local))
       call fv3io%get_var("tsfcl",nlon,nlat_local,r2d8b)
       !-- skin temperature on land
       this%ges_tsk(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%get_var("tsfc",nlon,nlat_local,r2d8b)
       !-- skin temperature composite
       this%tsk_comp(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%get_var("tsnow_land",nlon,nlat_local,r2d8b)
       !-- snow temperautre on land
       this%ges_soilt1(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%get_var("qwv_surf_land",nlon,nlat_local,r2d8b)
       !-- snow temperautre on land
       this%ges_qvg(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%get_var("clw_surf_land",nlon,nlat_local,r2d8b)
       !-- snow temperautre on land
       this%ges_qcg(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       call fv3io%get_var("snodl",nlon,nlat_local,r2d8b) 
       !--  snodl is snow depth on land, units [mm], convert to [m]
       this%sno(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)*1.e-3

       call fv3io%get_var("sncovr",nlon,nlat_local,r2d8b) 
       !-- snow cover: 0-1
       this%sncovr(:,this%fv3_layout_begin(id):this%fv3_layout_end(id))=r2d8b(:,:)

       deallocate(r2d8b)

       allocate(r3d8b(nlon,nlat_local,nz))
       call fv3io%get_var("tslb",nlon,nlat_local,nz,r3d8b)
       this%ges_tslb(:,this%fv3_layout_begin(id):this%fv3_layout_end(id),:)=r3d8b(:,:,:)
       call fv3io%get_var("smois",nlon,nlat_local,nz,r3d8b)
       this%ges_smois(:,this%fv3_layout_begin(id):this%fv3_layout_end(id),:)=r3d8b(:,:,:)

       deallocate(r3d8b)
       call fv3io%close
    enddo
    write(6,*) 'this%tsk=',maxval(this%ges_tsk),minval(this%ges_tsk)
    write(6,*) 'this%soilt1=',maxval(this%ges_soilt1),minval(this%ges_soilt1)
    write(6,*) 'this%qvg=',maxval(this%ges_qvg),minval(this%ges_qvg)
    do k=1,nz
       write(6,*) 'this%ges_tslb=',maxval(this%ges_tslb(:,:,k)),minval(this%ges_tslb(:,:,k))
    enddo
    do k=1,nz
       write(6,*) 'this%ges_smois=',maxval(this%ges_smois(:,:,k)),minval(this%ges_smois(:,:,k))
    enddo
endif

    call mpi_barrier(mpi_comm_world,iret)
if(mype==1) then
    call MPI_Send(this%sumqc, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
    call MPI_Send(this%ges_q1, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
endif
if(mype==2) then
    call MPI_Send(this%ges_p1, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
    call MPI_Send(this%ges_t1, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
    call MPI_Send(this%qsatg, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
endif

if(mype==3) then
    call MPI_Send(this%qinc, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
endif
if(mype==4) then
    call MPI_Send(this%tinc, nlon*nlat, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,iret)
endif

if(mype==0) then
   call MPI_Recv(this%tinc, nlon*nlat, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%qinc, nlon*nlat, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%sumqc, nlon*nlat, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%ges_q1, nlon*nlat, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%ges_p1, nlon*nlat, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%ges_t1, nlon*nlat, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)
   call MPI_Recv(this%qsatg, nlon*nlat, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,iret)

   this%qinc=this%ges_q1-this%qinc
   this%tinc=this%ges_t1-this%tinc
   write(6,*) 'mype=0 qinc=',maxval(this%qinc),minval(this%qinc)
   write(6,*) 'mype=0 tinc=',maxval(this%tinc),minval(this%tinc)
   write(6,*) 'mype=0 sumqc=',maxval(this%sumqc),minval(this%sumqc)
   write(6,*) 'mype=0 ges_q1=',maxval(this%ges_q1),minval(this%ges_q1)
   write(6,*) 'mype=0 ges_p1=',maxval(this%ges_p1),minval(this%ges_p1)
   write(6,*) 'mype=0 ges_t1=',maxval(this%ges_t1),minval(this%ges_t1)
   write(6,*) 'mype=0 qsatg=',maxval(this%qsatg),minval(this%qsatg)
endif

  end subroutine read_ges

subroutine gsl_gen_coast_prox(nlon,nlat,isli,coast_prox)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gsd_gen_coast_prox calculate coast proximity based on isli
!   prgmmr: Hu          org: GSD                date: 2015-01-14
!
! abstract:  This routine does the following things:
!              1) calculate coast proximity based on isli 
! 
! 
! program history log:
!   2014-01-22  Hu - original code
!
!   input argument list:
!
!   output argument list:
!
!   comments:
!
! attributes:
!$$$
  use kinds, only: r_kind,i_kind

  implicit none

! Declare passed variables

  integer, intent(in) :: nlon,nlat
  real(r_kind),intent(in) :: isli(nlon,nlat)
  real(r_kind),intent(inout) :: coast_prox(nlon,nlat)

  integer(i_kind) i,j,ico
  integer(i_kind) ia,ib,ja,jb,ic,jc,nco,nip

!*******************************************************************************
!
! water, land, seaice index
!  isli = 0 water, =1 land, =2 sea ice (on water)
  coast_prox=0.0_r_kind
  ico = 3
  do j=1,nlat
     ja = max(1   ,j-ico)
     jb = min(nlat,j+ico+1)
     do i=1,nlon
       if (abs(isli(i,j)-1.0_r_kind) <0.001_r_kind .or. &
           abs(isli(i,j)-2.0_r_kind) <0.001_r_kind ) then
           ia = max(1   ,i-ico)
           ib = min(nlon,i+ico+1)
           nco = 0
           nip = 0
           do jc=ja,jb
           do ic=ia,ib
              if (abs(isli(ic,jc)-1.0_r_kind) <0.001_r_kind .or. &
                  abs(isli(ic,jc)-2.0_r_kind) <0.001_r_kind ) nco = nco+1
                nip = nip+1
            end do
            end do
            coast_prox(i,j) = float(nco)/float (nip)
         end if
     end do
  end do
!
  return
end subroutine gsl_gen_coast_prox


end module module_bkio_fv3lam_parall
