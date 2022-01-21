module module_read_NSSL_refmosaic
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-01-20
!
! ABSTRACT: 
!     This routine read in NSSL reflectivity mosaic fields  
!     from different resources:
!
!     tversion=8  : NSSL 8 tiles netcdf
!     tversion=81 : NCEP 8 tiles binary
!     tversion=4  : NSSL 4 tiles binary
!     tversion=14 : NSSL 4 tiles netcdf
!     tversion=1  : NSSL 1 tile grib2 for 33 levels
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  mosaic_files
!
!   OUTPUT FILES:
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
!
  use kinds, only: r_kind,i_kind

  implicit none

  public :: read_nsslref
!
!-------------------------------------------------------------------------

! set default to private
  private
  type :: read_nsslref
     integer           ::   mscNlon   ! number of longitude of mosaic data
     integer           ::   mscNlat   ! number of latitude of mosaic data
     integer           ::   mscNlev   ! number of vertical levels of mosaic data
     real, allocatable :: mscValue3d(:,:,:)    ! reflectivity

     real              :: lonMin,latMin,lonMax,latMax
     real*8            :: dlon,dlat
     integer           :: maxlvl

     integer           :: ilevel
     integer,allocatable :: levelheight(:)

     real     ::  rthresh_ref,rthresh_miss 
     integer  ::  tversion

     character(len=256) ::   mosaicfile
     logical            :: if_fileexist=.false.
     integer            :: var_scale
    contains
      procedure :: init
      procedure :: readtile
      procedure :: close
  end type read_nsslref
!
! constants
!
contains

  subroutine init(this,tversion,mypelocal,datapath)
!                .      .    .                                       .
! subprogram:    init
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
!
    implicit none

    class(read_nsslref) :: this
    integer,intent(in)  :: tversion
    integer,intent(in)  :: mypelocal
    character(len=180),intent(in) :: dataPath

    integer :: k,n,nlevel
    character*256   filenameall(200)
!
!**********************************************************************
!
    this%tversion=tversion

    if( tversion == 8 .or. tversion == 14) then
       this%maxlvl = 31
       this%rthresh_ref=-500.0
       this%rthresh_miss=-5000.0
    elseif( tversion == 81 ) then
       this%maxlvl = 31
       this%rthresh_ref=-90.0
       this%rthresh_miss=-900.0
    elseif( tversion == 4 ) then
       this%maxlvl = 33
       this%rthresh_ref=-500.0
       this%rthresh_miss=-5000.0
    elseif( tversion == 1 ) then
       this%maxlvl = 33
       this%rthresh_ref=-90.0
       this%rthresh_miss=-900.0
    else
       write(*,*) 'unknow tversion !'
       stop 1234
    endif

    allocate(this%levelheight(this%maxlvl))
    this%levelheight(1)=500
    do k=2,11 
      this%levelheight(k)=this%levelheight(k-1)+250
    enddo
    do k=12,23 
      this%levelheight(k)=this%levelheight(k-1)+500
    enddo
    do k=24,this%maxlvl 
      this%levelheight(k)=this%levelheight(k-1)+1000
    enddo
!    write(*,*) 'set mosaic level=',this%levelheight

    this%mosaicfile='none'
    if ( tversion == 1 ) then
       open(10,file='filelist_mrms',form='formatted',err=300)
       do n=1,200
          read(10,'(a)',err=200,end=400) filenameall(n)
       enddo
300    write(6,*) 'read_grib2 open filelist_mrms failed ',mypelocal
       stop(555)
200    write(6,*) 'read_grib2 read msmr file failed ',n,mypelocal
       stop(555)
400    nlevel=n-1
       close(10)

       if(nlevel .gt. this%maxlvl) then
          write(*,*) 'vertical level is too large:',nlevel,this%maxlvl
          stop 666
       endif
       if(mypeLocal <= this%maxlvl) then
          this%mosaicfile=trim(filenameall(mypeLocal))
          write(*,*) 'process level:',mypeLocal,trim(this%mosaicfile)
       endif

    elseif(tversion==81) then
       if(mypeLocal <= 8) then
          write(this%mosaicfile,'(a,a,I1)') trim(dataPath), 'mosaic_t',mypeLocal
          write(*,*) 'process tile:',trim(this%mosaicfile)
       endif
    else
       if(mypeLocal <= tversion) then
          write(this%mosaicfile,'(a,a,I1)') trim(dataPath), 'mosaic_t',mypeLocal
          write(*,*) 'process tile:',trim(this%mosaicfile)
       endif
    endif

  end subroutine init

  subroutine close(this)
!                .      .    .                                       .
! subprogram:    init
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
!
    implicit none

    class(read_nsslref) :: this

    this%tversion=0
    this%mosaicfile='none'
    this%maxlvl=0
    if(allocated(this%levelheight)) deallocate(this%levelheight)
    if(allocated(this%mscValue3d)) deallocate(this%mscValue3d)

  end subroutine close
!

  subroutine readtile(this,mypelocal)
!                .      .    .                                       .
! subprogram:    init
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
!
    implicit none
!
    include 'netcdf.inc'

    class(read_nsslref) :: this
    integer, intent(in) :: mypelocal
!
    integer  :: tversion
    logical  :: fileexist
    integer  :: var_scale

    integer  :: NCID
    integer  :: ntot, ntot2d, mt
    integer  :: nx,ny,nz
    integer  :: yr, mo, da, hr, mn, sc
    real*8   :: rdx,rdy
    real*4   :: rdx4,rdy4
    real     :: rlatmax,rlonmin

    integer ::   mscNlon   ! number of longitude of mosaic data
    integer ::   mscNlat   ! number of latitude of mosaic data
    integer ::   mscNlev   ! number of vertical levels of mosaic data
    real, allocatable :: msclon(:)        ! longitude of mosaic data
    real, allocatable :: msclat(:)        ! latitude of mosaic data
    real, allocatable :: msclev(:)        ! level of mosaic data
    real, allocatable :: mscValue(:,:)    ! reflectivity
    integer  :: height
    integer  :: worklevel
    integer  :: status
    
    real   :: lonMin,latMin,lonMax,latMax
    real*8 :: dlon,dlat

    integer  :: i,j,k
    integer*2, dimension(:),   allocatable  :: var

    character*256   mosaicfile 
!
!   deal with certain tile
!
    mosaicfile=trim(this%mosaicfile)
    tversion=this%tversion
    fileexist=.false.
   
    if( tversion == 8 .or. tversion == 14) then
       call ifexist_file(mosaicfile,status)
       fileexist=status .EQ. NF_NOERR
    elseif( tversion == 81) then
       inquire(file=trim(mosaicfile),exist=fileexist)
       if(mypeLocal > 8) fileexist=.false.
    elseif(tversion == 4) then
       open(99,file=trim(mosaicfile),form='unformatted',access='direct',&
             recl=6*4,status='old',err=225)
       rewind(99)
       read(99,rec=1,err=225) yr, mo, da, hr, mn, sc
       fileexist=.true.
225    continue
       close(99)
    elseif(tversion == 1) then
       inquire(file=trim(mosaicfile),exist=fileexist)
       if(mypeLocal > this%maxlvl) fileexist=.false.
    endif

    if(.not.fileexist) then
       this%if_fileexist=.false.
       return
    else
       this%if_fileexist=.true.
    endif

    if( tversion == 14 ) then
       call GET_DIM_ATT_Mosaic(mosaicfile,mscNlon,mscNlat,mscNlev, &
                 lonMin,latMin,lonMax,latMax,dlon,dlat)
       var_scale=10
    elseif( tversion == 8 ) then
         call GET_DIM_ATT_Mosaic8(mosaicfile,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat)
         var_scale=10
    elseif( tversion == 81 ) then
         call read_ncep_binary_head(mypelocal-1,mosaicfile,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,rdx4,rdy4)
         nx=mscNlon
         ny=mscNlat
         nz=mscNlev
         dlon=rdx4
         dlat=rdy4
         var_scale=1
    elseif( tversion == 4 ) then
         call read_head_Mosaic4(mosaicfile,nx,ny,nz,rlonmin,rlatmax,&
                   rdx,rdy,var_scale)
         mscNlon=nx
         mscNlat=ny
         mscNlev=nz
         dlon=rdx
         dlat=rdy
         lonMin=rlonmin
         lonMax=lonMin+dlon*(mscNlon-1)
         latMax=rlatmax
         latMin=latMax-dlat*(mscNlat-1)
    elseif( tversion == 1 ) then
         call read_grib2_head(mosaicfile,nx,ny,nz,rlonmin,rlatmax,&
                   rdx,rdy)
         var_scale=1
         mscNlon=nx
         mscNlat=ny
         mscNlev=nz
         dlon=rdx
         dlat=rdy
         lonMin=rlonmin
         lonMax=lonMin+dlon*(mscNlon-1)
         latMax=rlatmax
         latMin=latMax-dlat*(mscNlat-1)
    else
         write(*,*) ' unknown tile version !!!'
         stop 123
    endif

    if(mypelocal==1) then
       write(*,*) 'mscNlon,mscNlat,mscNlev'
       write(*,*) mscNlon,mscNlat,mscNlev
       write(*,*) 'dlon,dlat,lonMin,lonMax,latMax,latMin'
       write(*,*) dlon,dlat,lonMin,lonMax,latMax,latMin
    endif

    this%var_scale=var_scale
    this%mscNlon=mscNlon
    this%mscNlat=mscNlat
    this%mscNlev=mscNlev

    allocate(msclon(this%mscNlon))
    allocate(msclat(this%mscNlat))
    allocate(msclev(this%mscNlev))

    do i=1,mscNlon
       msclon(i)=lonMin+(i-1)*dlon
    enddo
    do i=1,mscNlat
       msclat(i)=latMin+(i-1)*dlat
    enddo
    deallocate(msclon)
    deallocate(msclat)
    deallocate(msclev)
!
!  ingest mosaic file and interpolation
!  
    if(tversion == 1) mscNlev=1
    allocate(this%mscValue3d(mscNlon,mscNlat,mscNlev))
    allocate(mscValue(mscNlon,mscNlat))

    if( tversion == 8 .or. tversion == 14) then
       call OPEN_Mosaic(mosaicfile, NCID)

       if(tversion == 14 ) then
          call Check_DIM_ATT_Mosaic(NCID,mscNlon,mscNlat,mscNlev,  &
               lonMin,latMin,lonMax,latMax,dlon,dlat)
       elseif(tversion == 8 ) then
          call Check_DIM_ATT_Mosaic8(NCID,mscNlon,mscNlat,mscNlev,  &
               lonMin,latMin,lonMax,latMax,dlon,dlat)
       endif
       write(*,*) mscNlon,mscNlat,mscNlev
       write(*,*) 'Area of tile=',lonMin,latMin,lonMax,latMax,dlon,dlat
!
! we used to do it level by level. This needs a big memory to process
!
       do k=1, mscNlev
          call  GET_Mosaic_sngl_Mosaic(NCID,mscNlon,mscNlat,k,mscValue)
          do j=1,ny
             do i=1,nx
                this%mscValue3d(i,j,k)=mscValue(i,j)
             enddo
          enddo
       enddo  ! mscNlev

       call CLOSE_Mosaic(NCID)
    elseif(tversion == 81) then
       call read_ncep_binary_value(mypelocal-1,mosaicfile,mscNlon,mscNlat,mscNlev, &
                   this%mscValue3d)
    elseif(tversion == 4) then
       ntot = nx*ny*nz
       allocate(var(ntot))
       call read_data_Mosaic4(mosaicfile,ntot,var)
       do k=1,mscNlev
          ntot2d=nx*ny*(k-1)
          do j=1,ny
             do i=1,nx
                this%mscValue3d(i,j,k) = var(ntot2d+(j-1)*nx+i)
             enddo
          enddo
       enddo
       deallocate(var)
    elseif(tversion == 1) then
       ntot = nx*ny*nz
       call read_grib2_sngle(mosaicfile,ntot,height,mscValue)
       if(this%levelheight(mypeLocal) .eq. height) then
            worklevel=mypeLocal
       else
            worklevel=0
            do k=1,this%maxlvl
               if (this%levelheight(k) .eq. height) worklevel=k
            enddo
            if(worklevel==0) then
               write(6,*) 'Error, cannot find working level', &
                         mypeLocal,this%levelheight(mypeLocal), height
               stop 12345
            else
               write(6,*) 'Find new level for core ',mypeLocal,' new level is',worklevel   
            endif
       endif
       this%ilevel=worklevel
       this%mscValue3d(:,:,1)=mscValue(:,:)
       write(*,*) 'level max min height',mypeLocal,maxval(this%mscValue3d),minval(this%mscValue3d),height
    else
       write(*,*) 'unknown type'
       stop 1234
    endif

    deallocate(mscValue)

    this%lonMin=lonMin
    this%latMin=latMin
    this%lonMax=lonMax
    this%latMax=latMax
    this%dlon=dlon
    this%dlat=dlat

  end subroutine readtile

end module module_read_NSSL_refmosaic
