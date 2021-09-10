!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
module module_surface

  use kinds, only: r_kind,r_single,i_kind,i_byte,i_byte2
  implicit none

  public :: use_surface
!
!-------------------------------------------------------------------------

! set default to private
  private
  type :: use_surface
      integer :: nvar
      integer :: nvar3d
      character(len=20),allocatable :: var_rrfs(:),var_rap(:)
      integer :: halo
      integer :: nlat,nlon,nlev
      integer :: nlat_target,nlon_target
      integer(i_byte2), allocatable :: index_x(:,:)
      integer(i_byte2), allocatable :: index_y(:,:)
      integer(i_byte2), allocatable :: index_x_nomatch(:,:)
      integer(i_byte2), allocatable :: index_y_nomatch(:,:)
    contains
      procedure :: init
      procedure :: build_mapindex
      procedure :: set_varname
      procedure :: use_sfc
      procedure :: close
  end type use_surface
!
! constants
!
contains
   
  subroutine set_varname(this)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!   output argument list:
!
!
! 3D fields
! 1. soil moisture SMOIS (9 levels) to smois (9 levels)
! 2. Soil temperature TSLB (9 levels) to sh2o (9 levels). 
!    2.1:    smfr=smois-sh2o
! 3. Liquid soil water SH2O (9 levels) to tslb(9 levels) and tiice (9 levels)
!
! 2D fields
! 4. TSK to tsfc
! 5. SNOW to weasdl and sheleg
! 6. SNOWH to snodl
! 7, SOILT1 to tsnow_land and tsnow_ice (when seaice > 0.025)
! 8. CANWAT to canopy
! 9. SNOWC to sncovr and sncovr_ice (when seaice > 0.025)
! 10. SEAICE to fice
! 11. QVG to qwv_surf_ice and qwv_surf_land
! 12. QCG to  clw_surf_ice and clw_surf_land
! 13. slmsk = 2 when seaice > 0.025
!
    implicit none

    class(use_surface) :: this
!
    this%nvar=12
    this%nvar3d=3
!
    if(allocated(this%var_rrfs)) deallocate(this%var_rrfs)
    allocate(this%var_rrfs(this%nvar))
    if(allocated(this%var_rap)) deallocate(this%var_rap)
    allocate(this%var_rap(this%nvar))
!
! Soil moisture: 3D float
    this%var_rap(1)="SMOIS"
    this%var_rrfs(1)="smois"  !  double  smc?
! Liquid soil water: 3D float
    this%var_rap(2)="SH2O"
    this%var_rrfs(2)="sh2o"    ! double slc?
!    this%var_rrfs(2)="smfr"    ! double =smois-sh2o

! Soil temperature: 3D float
    this%var_rap(3)="TSLB"
    this%var_rrfs(3)="tslb"   ! double stc?
!    this%var_rrfs(3)="tiice"   ! double

! SURFACE SKIN TEMPERATURE: 2D float
    this%var_rap(4)="TSK"
    this%var_rrfs(4)="tsfc"   ! double tsfcl?

! SNOW WATER EQUIVALENT: 2D float
    this%var_rap(5)="SNOW"
    this%var_rrfs(5)="weasdl"  ! double 
!    this%var_rrfs(5)="sheleg"  ! double 

! PHYSICAL SNOW DEPTH: 2D float
    this%var_rap(6)="SNOWH"
    this%var_rrfs(6)="snodl"   !double

! TEMPERATURE INSIDE SNOW: 2D float
! tsnow_land over land and tsnow_ice over sea ice
    this%var_rap(7)="SOILT1"
    this%var_rrfs(7)="tsnow_land" ! double
!    this%var_rrfs(7)="tsnow_ice" ! double

! CANOPY WATER: 2D float
    this%var_rap(8)="CANWAT"
    this%var_rrfs(8)="canopy"

! FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER): 2D float
! sncovr over land and sncovr_ice over sea ice
    this%var_rap(9)="SNOWC"
    this%var_rrfs(9)="sncovr"   ! double
!    this%var_rrfs(9)="sncovr_ice"   ! double

! SEA ICE FLAG: 2D float
    this%var_rap(10)="SEAICE"
    this%var_rrfs(10)="fice"  !  double

! WATER VAPOR MIXING RATIO AT THE SURFACE
    this%var_rap(11)="QVG"
    this%var_rrfs(11)="qwv_surf_land"  !  double
!    this%var_rrfs(11)="qwv_surf_ice"  !  double

! CLOUD WATER MIXING RATIO AT THE GROUND SURFACE
    this%var_rap(12)="QCG"
    this%var_rrfs(12)="clw_surf_land"  !  double
!    this%var_rrfs(12)="clw_surf_ice"  !  double
!
!
  end subroutine set_varname

  subroutine use_sfc(this,rapfile,rrfsfile)
!                .      .    .                                       .
! subprogram:   build_mapindex
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
    use module_ncio, only : ncio
  
    implicit none

    character*80,intent(in) :: rapfile
    character*80,intent(in) :: rrfsfile

    class(use_surface) :: this
    type(ncio)     :: raphrrr,rrfs
!
    real(r_single),allocatable,target :: fld2d4b(:,:)
    real(r_single),allocatable,target :: fld3d4b(:,:,:)
!
    real(r_single),allocatable,target :: tmp2d4b(:,:)
    real(r_single),allocatable,target :: tmp3d4b(:,:,:)
    real(r_kind),  allocatable,target :: tmp2d8b(:,:)
    real(r_kind),  allocatable,target :: tmp3d8b(:,:,:)
    integer(i_byte),allocatable :: landmask_rrfs(:,:)
!
    integer  :: nx_rap,ny_rap
    integer  :: nx_rrfs,ny_rrfs,nz_rrfs
    integer  :: i,j,ix,jx
    character(len=20) :: thisvar_rrfs,thisvar_rap
    integer  :: k,kx 
!

    nx_rrfs=this%nlon
    ny_rrfs=this%nlat
    nz_rrfs=this%nlev
    nx_rap=this%nlon_target
    ny_rap=this%nlat_target
!
    call raphrrr%open(trim(rapfile),"r",200)
    allocate(fld2d4b(nx_rap,ny_rap))
    call raphrrr%get_var("SEAICE",nx_rap,ny_rap,fld2d4b)
!
    call rrfs%open(trim(rrfsfile),"r",200)
    allocate(landmask_rrfs(nx_rrfs,ny_rrfs))
    allocate(tmp2d8b(nx_rrfs,ny_rrfs))
    call rrfs%get_var("slmsk",nx_rrfs,ny_rrfs,tmp2d8b)
    do j=1,ny_rrfs
       do i=1,nx_rrfs
          ix=this%index_x(i,j)
          jx=this%index_y(i,j)
          if( (ix >= 1 .and. ix <= nx_rap) .and. &
              (jx >= 1 .and. jx <= ny_rap) ) then
             if(fld2d4b(ix,jx) > 0.025_r_single) then
                tmp2d8b(i,j)=2.0_r_kind
             endif
          endif
        enddo
    enddo
    call rrfs%close()

    call rrfs%open(trim(rrfsfile),"w",200)
    call rrfs%replace_var("slmsk",nx_rrfs,ny_rrfs,tmp2d8b)

    do j=1,ny_rrfs
      do i=1,nx_rrfs
         landmask_rrfs(i,j)=int(tmp2d8b(i,j))
      enddo
    enddo

    deallocate(tmp2d8b)
    call rrfs%close()
!
!
    do k=1,this%nvar
       thisvar_rrfs=this%var_rrfs(k)
       thisvar_rap=this%var_rap(k)
       write(*,*) "=============================================="
       write(*,'(4a)') "Working to replace rrfs variable ",trim(thisvar_rrfs), &
                  " from RAP/HRRR ",trim(thisvar_rap)
       if(k <= this%nvar3d) then
          allocate(tmp3d4b(nx_rap,ny_rap,nz_rrfs))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,nz_rrfs,tmp3d4b)

          allocate(tmp3d8b(nx_rrfs,ny_rrfs,nz_rrfs))
          call rrfs%open(trim(rrfsfile),"r",200)
          call rrfs%get_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,nz_rrfs,tmp3d8b)
          call rrfs%close()

          do j=1,ny_rrfs
            do i=1,nx_rrfs
               ix=this%index_x(i,j)
               jx=this%index_y(i,j)
               if( (ix >= 1 .and. ix <= nx_rap) .and. &
                   (jx >= 1 .and. jx <= ny_rap) ) then
                  do kx=1,nz_rrfs
                     tmp3d8b(i,j,kx)=tmp3d4b(ix,jx,kx)
                  enddo
               endif
            enddo
          enddo
          deallocate(tmp3d4b)

          call rrfs%open(trim(rrfsfile),"w",200)
          call rrfs%replace_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,nz_rrfs,tmp3d8b)
          if(trim(thisvar_rrfs) == "tslb") then
             call rrfs%replace_var("tiice",nx_rrfs,ny_rrfs,nz_rrfs,tmp3d8b)
          elseif(trim(thisvar_rrfs) == "smois") then
             allocate(fld3d4b(nx_rrfs,ny_rrfs,nz_rrfs))
             fld3d4b=tmp3d8b  ! save for calculating (smois-sh2o)
          elseif(trim(thisvar_rrfs) == "sh2o") then
             tmp3d8b=fld3d4b - tmp3d8b
             deallocate(fld3d4b)
             call rrfs%replace_var("smfr",nx_rrfs,ny_rrfs,nz_rrfs,tmp3d8b)
          endif
          deallocate(tmp3d8b)
          call rrfs%close()
       else
          allocate(tmp2d4b(nx_rap,ny_rap))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,tmp2d4b)

          allocate(tmp2d8b(nx_rrfs,ny_rrfs))
          call rrfs%open(trim(rrfsfile),"r",200)
          call rrfs%get_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,tmp2d8b)
          call rrfs%close()
        
          if(maxval(tmp2d8b) - minval(tmp2d8b) < 1.0e-10) then
            do j=1,ny_rrfs
              do i=1,nx_rrfs
                ix=this%index_x_nomatch(i,j)
                jx=this%index_y_nomatch(i,j)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                  tmp2d8b(i,j)=tmp2d4b(ix,jx)
                endif
              enddo
            enddo
          else
            do j=1,ny_rrfs
              do i=1,nx_rrfs
                ix=this%index_x(i,j)
                jx=this%index_y(i,j)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                  tmp2d8b(i,j)=tmp2d4b(ix,jx)
                endif
              enddo
            enddo
          endif
          deallocate(tmp2d4b)

          call rrfs%open(trim(rrfsfile),"w",200)
          call rrfs%replace_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,tmp2d8b)
          if(trim(thisvar_rrfs) == "weasdl") then
            call rrfs%replace_var("sheleg",nx_rrfs,ny_rrfs,tmp2d8b)
          elseif(trim(thisvar_rrfs) == "qwv_surf_land") then
            call rrfs%replace_var("qwv_surf_ice",nx_rrfs,ny_rrfs,tmp2d8b)
          elseif(trim(thisvar_rrfs) == "clw_surf_land") then
            call rrfs%replace_var("clw_surf_ice",nx_rrfs,ny_rrfs,tmp2d8b)
          endif
          deallocate(tmp2d8b)
          call rrfs%close()
!
! for SOILT1 and SNOWC
          if(trim(thisvar_rap) == "SOILT1" .or. trim(thisvar_rap) == "SNOWC") then
             if(trim(thisvar_rap) == "SOILT1") thisvar_rrfs="tsnow_ice"
             if(trim(thisvar_rap) == "SNOWC")  thisvar_rrfs="sncovr_ice"

             allocate(tmp2d4b(nx_rap,ny_rap))
             call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,tmp2d4b)

             allocate(tmp2d8b(nx_rrfs,ny_rrfs))
             call rrfs%open(trim(rrfsfile),"r",200)
             call rrfs%get_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,tmp2d8b)
             call rrfs%close()

             if(maxval(tmp2d8b) - minval(tmp2d8b) < 1.0e-10) then
               do j=1,ny_rrfs
                 do i=1,nx_rrfs
                   ix=this%index_x_nomatch(i,j)
                   jx=this%index_y_nomatch(i,j)
                   if( (ix >= 1 .and. ix <= nx_rap) .and. &
                       (jx >= 1 .and. jx <= ny_rap) ) then
                      if(fld2d4b(ix,jx) > 0.025) then
                        tmp2d8b(i,j)=tmp2d4b(ix,jx)
                      endif
                   endif
                 enddo
               enddo
             else
               do j=1,ny_rrfs
                 do i=1,nx_rrfs
                    ix=this%index_x(i,j)
                    jx=this%index_y(i,j)
                    if( (ix >= 1 .and. ix <= nx_rap) .and. &
                        (jx >= 1 .and. jx <= ny_rap) ) then
                       if(fld2d4b(ix,jx) > 0.025) then
                          tmp2d8b(i,j)=tmp2d4b(ix,jx)
                       endif
                    endif
                 enddo
               enddo
             endif
             deallocate(tmp2d4b)

             call rrfs%open(trim(rrfsfile),"w",200)
             call rrfs%replace_var(trim(thisvar_rrfs),nx_rrfs,ny_rrfs,tmp2d8b)
             deallocate(tmp2d8b)
             call rrfs%close()
          endif
!
       endif
    enddo
!
    deallocate(fld2d4b)
    deallocate(landmask_rrfs)
    call raphrrr%close()

  end subroutine use_sfc

  subroutine build_mapindex(this,map,rlon,rlat,landmask_this,landmask_target)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    map: grid convert class
!    landmask_this: land mask for this grid
!    landmask_target: land mask for target grid
!
!   output argument list:
!
    use module_map_utils, only : map_util
    implicit none

    class(use_surface) :: this
    type(map_util),  intent(in) :: map
    integer(i_byte), intent(in) :: landmask_this(this%nlon,this%nlat)
    integer(i_byte), intent(in) :: landmask_target(this%nlon_target,this%nlat_target)
    real(r_single),  intent(in) :: rlon(this%nlon,this%nlat)
    real(r_single),  intent(in) :: rlat(this%nlon,this%nlat)
!
    real(r_single) :: xc,yc
    integer  :: ixc,jyc
    integer  :: i,j,ix,jx
    integer  :: ii,jj,ixx,jxx
    real(r_single) :: dist,thisdist
!
    do j=1,this%nlat
    do i=1,this%nlon
       call map%tll2xy(rlon(i,j),rlat(i,j),xc,yc)
       ixc=int(xc+0.5)
       jyc=int(yc+0.5)
       if( (ixc > 0 .and. ixc < this%nlon_target + 1) .and. &
           (jyc > 0 .and. jyc < this%nlat_target + 1) ) then
          ix=min(max(ixc,1),this%nlon_target)
          jx=min(max(jyc,1),this%nlat_target)
          this%index_x_nomatch(i,j)=ix
          this%index_y_nomatch(i,j)=jx 
          if( landmask_target(ix,jx) == 2 ) then
! snow/ice in RAP HRRR
             this%index_x(i,j)=ix
             this%index_y(i,j)=jx 
          elseif( landmask_this(i,j)==landmask_target(ix,jx) ) then
! match water and land
             this%index_x(i,j)=ix
             this%index_y(i,j)=jx 
          else
! dewl with not matched water and land
             dist=99999.0
             do jj=max(1,jx-this%halo),min(jx+this%halo,this%nlat_target)
               do ii=max(1,ix-this%halo),min(ix+this%halo,this%nlon_target)
                 if(landmask_this(i,j)==landmask_target(ii,jj)) then
                    thisdist=(jj-jx)*(jj-jx)+(ii-ix)*(ii-ix)
                    thisdist=sqrt(thisdist)
                    if(thisdist <= dist) then
                       dist=thisdist
                       ixx=ii
                       jxx=jj
                    endif
                 endif
               enddo
             enddo
             if(dist < 88888.0) then
               this%index_x(i,j)=ixx
               this%index_y(i,j)=jxx
             else
               !write(*,*) i,j,this%index_x(i,j),this%index_y(i,j)
             endif
          endif
       endif
    enddo
    enddo
!     
!
  end subroutine build_mapindex

  subroutine init(this,nlon,nlat,nlev,nlon_target,nlat_target,halo)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    nlon,nlat - grid dimension
!
!   output argument list:
!
    implicit none

    class(use_surface) :: this
    integer,intent(in) :: halo
    integer,intent(in) :: nlon,nlat,nlev
    integer,intent(in) :: nlon_target,nlat_target

    this%halo=halo
    this%nlon=nlon
    this%nlat=nlat
    this%nlev=nlev
    this%nlon_target=nlon_target
    this%nlat_target=nlat_target
    allocate(this%index_x(this%nlon,this%nlat))
    allocate(this%index_y(this%nlon,this%nlat))
    this%index_x=-999
    this%index_y=-999
    allocate(this%index_x_nomatch(this%nlon,this%nlat))
    allocate(this%index_y_nomatch(this%nlon,this%nlat))
    this%index_x=-999
    this%index_y=-999

  end subroutine init

  subroutine close(this)
!                .      .    .                                       .
! subprogram:    close
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

    class(use_surface) :: this

    this%nvar=0
    this%nvar3d=0
    this%halo=0
    this%nlon=0
    this%nlat=0
    this%nlev=0
    this%nlon_target=0
    this%nlat_target=0

    deallocate(this%index_x)
    deallocate(this%index_y)
    deallocate(this%index_x_nomatch)
    deallocate(this%index_y_nomatch)

    deallocate(this%var_rrfs)
    deallocate(this%var_rap)

  end subroutine close

end module module_surface
