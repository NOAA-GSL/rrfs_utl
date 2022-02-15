module module_imssnow
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-02-10
!
! ABSTRACT: 
!     This module includes tools to read in NESDIS NESDIS SNOW/ICE data from a grib file and  
!     map them into RR mass grid
!
! 
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  imssnow
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

  use grib2_read_snow_mod, only : read_grib2_head_dim,read_grib2_head_time
  use grib2_read_snow_mod, only : read_grib2_sngle

  implicit none
!
  public :: type_imssnow
!
!-------------------------------------------------------------------------

! set default to private
  private
  type :: type_imssnow
      integer :: nxobs
      integer :: nyobs
      integer :: iyear,imonth,iday,ihr,imm
      integer :: kgds(200)
      character(len=100) :: snowfile

      real, allocatable :: snowicec(:,:)
      integer, allocatable :: maskims(:,:)

    contains
      procedure :: init
      procedure :: read_imssnow
      procedure :: map2grid
      procedure :: close
  end type type_imssnow
!
! constants
!
contains


  subroutine init(this,snowfile,iflag)
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

    class(type_imssnow)         :: this
    character(len=*),intent(in) :: snowfile
    integer,intent(inout)       :: iflag
!
!  for grib2
!
    real :: rlatmin,rlonmin
    real*8  :: rdx,rdy
    integer :: nxobs,nyobs
    logical :: fileexist
    integer :: kgds(200)
    integer :: iyear, imonth, iday, ihr, imm

    iflag=0 
    inquire(file=trim(snowfile),exist=fileexist)
    if(fileexist) then
      call read_grib2_head_dim(snowfile,nxobs,nyobs,rlonmin,rlatmin,rdx,rdy,kgds)
      write(6,*) 'grib2 file grid =',nxobs,nyobs,rlonmin,rlatmin,rdx,rdy
      write(6,*) 'kgds=',kgds(1:11)
!
      call read_grib2_head_time(snowfile,iyear,imonth,iday,ihr,imm)
      write(6,'(a,5I5)') 'date: iyear,imonth,iday,ihr,imm=',&
                      iyear,imonth,iday,ihr,imm

    else
      iflag=1
      return
    endif

    this%kgds=kgds
    this%snowfile=trim(snowfile)
    this%nxobs=nxobs
    this%nyobs=nyobs

    allocate(this%snowicec(this%nxobs,this%nyobs))
    allocate(this%maskims(this%nxobs,this%nyobs))
    this%snowicec=0.0
    this%maskims=0

    this%iyear=iyear
    this%imonth=imonth
    this%iday=iday
    this%ihr=ihr
    this%imm=imm

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

    class(type_imssnow) :: this

    this%snowfile=""
    this%nxobs=0
    this%nyobs=0

    deallocate(this%snowicec)
    deallocate(this%maskims)

  end subroutine close

  subroutine read_imssnow(this)
!                .      .    .                                       .
! subprogram:    read_imssnow
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

    class(type_imssnow)         :: this
!
!  for grib2 ims snow
!
    integer :: ntot
    integer :: i,j
!
! sea ice and snow, and mask
    real, allocatable :: ICEC(:,:)
    logical*1, allocatable :: lmask(:,:)

!**********************************************************************
!
!            END OF DECLARATIONS....start of program

    allocate(ICEC(this%nxobs,this%nyobs))
    allocate(lmask(this%nxobs,this%nyobs))
    ICEC=0
    this%snowicec=0

! read field
    ntot=this%nxobs*this%nyobs
    call read_grib2_sngle(this%snowfile,ntot,ICEC,this%SNOWICEC,lmask)
    write(6,*) 'read in ice =',maxval(icec),minval(icec)
    write(6,*) 'read in snow =',maxval(this%SNOWICEC),minval(this%SNOWICEC)
!
    do J = 1, this%nxobs
      do I = 1, this%nyobs
         IF (icec(I,J) > 0. .OR. this%snowicec(I,J) > 0.) THEN
            this%snowicec(I,J) = 1.
         ELSE
            this%snowicec(I,J) = 0.
         ENDIF

!tgs Logical lmask --> .true. for land, .false. for water
         IF (lmask(I,J)) THEN
            this%maskims(i,j) = 1    ! land
         ELSE
            this%maskims(i,j) = 0    ! water
         ENDIF
      enddo
    enddo

    deallocate(icec)
    deallocate(lmask)

  end subroutine read_imssnow

  subroutine map2grid(this,xlandRR,nlonRR,nlatRR,xlonRR,ylatRR,snowiceRR,xlandIMS)
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2009-04-15
!
! ABSTRACT:
!     This routine map NESDIS SNOW/ICE data to RR grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  imssnow
!
!   OUTPUT FILES:  RRimssnow
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

    use gdswzd_mod, only : gdswzd
    implicit none

    class(type_imssnow)         :: this
!  grid
    integer, intent(in) :: nlonRR,nlatRR
    integer, intent(in):: xlandRR(nlonRR,nlatRR)
    real, intent(in):: xlonRR(nlonRR,nlatRR)    !
    real, intent(in):: ylatRR(nlonRR,nlatRR)    !
!
    real, intent(inout) :: snowiceRR(nlonRR,nlatRR)    !
    real, intent(inout) :: xlandIMS(nlonRR,nlatRR)
!
    INTEGER JPDS(200),JGDS(200),KPDS_SRC(200),KGDS_SRC(200)
!
    REAL :: DTR
    REAL :: XPNMC8,YPNMC8,ENNMC8,ALNMC8,ORIENT8
    REAL :: XPNMCAF,YPNMCAF,ENNMCAF,ALNMCAF,ORIENTAF

    real :: YYLAT(1),XLONG(1), RM, RAD, XX(1),YY(1),X, Y
    integer :: IS,IP1,JS,JP1, IPOINT, JPOINT
!
    integer :: iland,KOUNT
    real    :: XRATIO,YRATIO,AREA11,AREA21,AREA12,AREA22,AREA
    integer :: i,j,k,ifound,LL,JPE,JPB,IPE,IPB,NK,MK
!
!
    integer                   :: nret,nout
    real                      :: dum
    real, parameter           :: undefined_value = -1.0

!
    kgds_src=this%kgds
    write(6,*) 'NESDIS grid info: kgds_src', kgds_src(1:11)
    if(this%nxobs==6144 .and. this%nyobs==6144) then

    else
       write(6,*) 'wrong dimension=',6144
    endif

    nout=0
    DO j=1,nlatRR
    DO i=1,nlonRR
       YYLAT(1)=ylatRR(i,j)
       XLONG(1)=xlonRR(i,j)
       call gdswzd(kgds_src,-1,1,undefined_value,XX,YY, &
                   XLONG,YYLAT,nret)
       X=XX(1)
       Y=YY(1)
       if (nret /= 1) then
         nout=nout+1
!        snowiceRR(i,j) = 0.
         cycle
       else
! check if X and Y are defined
	 if(X == undefined_value .or. Y == undefined_value) then
!           snowiceRR(i,j) = 0.0 
            cycle
         endif
	 if(X < 1.00001 .or. Y < 1.00001) then
            cycle
         endif
         IS  = NINT(X)
         IP1 = IS + 1
         JS  = NINT(Y)
         JP1 = JS + 1

         ipoint = nint(X)
         jpoint = nint(Y)

       end if

!  AND ONLY SEA POINTS ARE INTERPOLATED TO SEA POINTS (FOR ICE)
!  (NESDIS/IMS LAND MASK: SEA=0,LAND=1, WHILE THE RR MASK in XLAND IS: SEA=2,
!  LAND=1).
!

       ILAND = 1
       IF( xlandRR(i,j) == 0 ) ILAND=0

!tgs - buggy      IF( int(maskims(i,j)) == ILAND ) THEN
       IF( this%maskims(is,js) == ILAND ) THEN
          snowiceRR(i,j) = this%snowicec(IPOINT,JPOINT)
       ELSE
!
!  NEAREST NEIGHBOR NOT SAME SFC TYPE, SO USE ALL 4 SURROUNDING POINTS
!
          KOUNT = 0
!
          XRATIO = X - REAL(IS)
          YRATIO = Y - REAL(JS)
!
          AREA11 = (1.0E0 - XRATIO) * (1.0E0 - YRATIO)
          AREA21 = XRATIO * (1.0E0 - YRATIO)
          AREA12 = (1.0E0 - XRATIO) * YRATIO
          AREA22 = XRATIO * YRATIO
!
!
          IF( this%maskims(IS, JS) .EQ. ILAND) THEN
             KOUNT  = KOUNT + 1
             AREA   = AREA11
             IPOINT = IS
             JPOINT = JS
          END IF
!
          IF( this%maskims(IS, JP1) .EQ. ILAND ) THEN
             KOUNT = KOUNT +1
             IF (KOUNT .EQ. 1) THEN
                IPOINT = IS
                JPOINT = JP1
             ELSEIF (AREA12 .GT. AREA) THEN
                AREA   = AREA12
                IPOINT = IS
                JPOINT = JP1
             END IF
          END IF
!
          IF( this%maskims(IP1, JS) .EQ. ILAND ) THEN
             KOUNT = KOUNT + 1
             IF (KOUNT .EQ. 1) THEN
                AREA   = AREA21
                IPOINT = IP1
                JPOINT = JS
             ELSEIF (AREA21 .GT. AREA) THEN
                AREA   = AREA21
                IPOINT = IP1
                JPOINT = JS
             END IF
          END IF
!
!
          IF( this%maskims(IP1, JP1) .EQ. ILAND ) THEN
             KOUNT = KOUNT + 1
             IF (KOUNT .EQ. 1) THEN
                AREA   = AREA22
                IPOINT = IP1
                JPOINT = JP1
             ELSEIF (AREA22 .GT. AREA) THEN
                AREA   = AREA22
                IPOINT = IP1
                JPOINT = JP1
             END IF
          END IF
!
!     DETERMINE SNO/ICE USING NEAREST NEIGHBOR WITH SAME SFC TYPE 
!
          IF(KOUNT .GT. 0) THEN
             snowiceRR(i,j) = this%snowicec(IPOINT,JPOINT)
          ELSE
!
!         NO IMMEDIATELY SURROUNDING POINTS IN THE 6144 X 6144 FIELD OF
!         SNOW/ICE HAVE THE SAME LAND-SEA TYPE AS THE model POINT.  THE

!         model POINT MAY BE SMALL ISLAND OR LAKE OR SMALL BAY OR PENNIN.
!         (INVARIABLY A SMALL LAKE IN ETA GRID)
!         SO EXPAND SEARCH RADIUS AND TAKE FIRST SFC TYPE MATCH
!
             IPOINT = NINT(X)
             JPOINT = NINT(Y)
!
!  Define the frame (no. of grid points) over which to search for
!    a matching land/water type from IMS data for the model gridpoint.
             ifound=0
             DO LL=1,16
                JPE = MIN (6144, JPOINT+LL)
                JPB = MAX (1 , JPOINT-LL)
                IPE = MIN (6144, IPOINT+LL)
                IPB = MAX (1 , IPOINT-LL)
!
                DO NK=IPB,IPE
                DO MK=JPB,JPE
                   IF ( this%maskims(nk,mk) == ILAND .and. ifound ==0 ) THEN
                      snowiceRR(i,j) = this%snowicec(NK,MK)
                      ifound=1
                   ENDIF
                ENDDO  ! MK
                ENDDO  ! NK
             ENDDO  ! LL
!
!  NO LAND/SEA MASK MATCHES FOUND, SO 
!     A) NORTH OF 55N, WE ASSIGN SNOW/ICE IRRESPECTIVE OF SFC TYPE
!     B) SOUTH OF 55N, WE KEEP A PRIORI ZERO DEFAULT
!   (THE "B" OPTION BEST FOR WARMER LATS OF U.S., WHERE THIS CONDITION 
!   IS VIRTUALLY ALWAYS A SMALL ETA LAKE WITH NO COUNTERPART WATER 
!   NEARBY IN THE NESDIS/IMS GRID, E.G., SALT LAKE, WHERE WE MUST
!   AVOID GETTING SEA-ICE OWING TO SURROUNDING SNOW COVER)
!
             IF ( ifound==0 ) THEN
                snowiceRR(i,j) = this%snowicec(IPOINT,JPOINT)
             ENDIF
          ENDIF   !  KOUNT .GT. 0
!
       ENDIF
!
!tgs - save NESDIS 4-km land/water mask
       xlandIMS(i,j)=float(this%maskims(IPOINT,JPOINT))
!

    ENDDO  ! nlon
    ENDDO  ! nlat
                                                                                                                       
    print*,"- WARNING: NUM of MODEL POINT OUTSIDE NESDIS GRID.",nout

end subroutine map2grid

end module module_imssnow
