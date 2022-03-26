PROGRAM check_process_imssnow
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-02-11
!
! ABSTRACT: 
!     This appllication uses NESDIS SNOW/ICE data from a grib file to
!          update FV3LAM ice and snow fields 
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

  use mpi
  use kinds, only: r_kind
  use module_io_fv3lam_bdy , only : io_fv3lam_bdy
  use module_update_bc     , only : update_bc_4side,update_bc_4side_gsiwind
  use module_update_bc     , only : update_bc_fv3uv2earth,update_bc_4side_wind
  use module_io_fv3lam_bk  , only : io_fv3lam_bk
  use mod_fv3lam_wind,       only : initial_wind_convert,reverse_grid_r

  implicit none
! MPI variables
  type(io_fv3lam_bdy) :: fv3bdy
  type(io_fv3lam_bk)  :: fv3bk
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! namelist
  integer :: fv3_io_layout_y
  integer :: bdy_update_type
  integer :: grid_type_fv3_regional
  namelist/setup/ fv3_io_layout_y,bdy_update_type,grid_type_fv3_regional
  logical :: ifexist
!
  integer :: numbdy
  character(len=20),allocatable :: bdyvar(:)
  integer :: nn
  logical :: grid_reverse_flag

  real(r_kind),allocatable:: ud(:,:,:)

!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program

! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
!  get namelist
!
  fv3_io_layout_y=1
  bdy_update_type=1  ! 1 update original bdy
                     ! 2 generate GSI bdy
  grid_type_fv3_regional=2

  inquire(file='namelist.updatebc', EXIST=ifexist )
  if(ifexist) then
    open(10,file='namelist.updatebc',status='old')
       read(10,setup)
    close(10)
    write(*,*) 'Namelist setup are:'
    write(*,setup)
  else
     write(*,*) 'No namelist file exist, use default values'
  endif

! check the grid type
  if( grid_type_fv3_regional == 1 ) then
     if(mype==0) write(6,*) 'FV3 regional input grid is  E-W N-S grid'
     grid_reverse_flag=.true.    ! grid is revered comparing to usual map view
  else if(grid_type_fv3_regional == 2) then
     if(mype==0) write(6,*) 'FV3 regional input grid is  W-E S-N grid'
     grid_reverse_flag=.false.   ! grid orientated just like we see on map view    
  else
     write(6,*) 'Error: FV3 regional input grid is unknown grid'
     stop 123
  endif

  if(bdy_update_type==1) then
     numbdy=14
     allocate(bdyvar(numbdy))
     bdyvar(1)="ps"
     bdyvar(2)="t"
     bdyvar(3)="w"
     bdyvar(4)="zh"
     bdyvar(5)="sphum"
     bdyvar(6)="liq_wat"
     bdyvar(7)="o3mr"
     bdyvar(8)="ice_wat"
     bdyvar(9)="rainwat"
     bdyvar(10)="snowwat"
     bdyvar(11)="graupel"
     bdyvar(12)="ice_aero"
     bdyvar(13)="liq_aero"
     bdyvar(14)="uv"

  elseif (bdy_update_type==2) then
     numbdy=16
     allocate(bdyvar(numbdy))
     bdyvar(1)="ps"
     bdyvar(2)="t"
     bdyvar(3)="w"
     bdyvar(4)="sphum"
     bdyvar(5)="liq_wat"
     bdyvar(6)="o3mr"
     bdyvar(7)="ice_wat"
     bdyvar(8)="rainwat"
     bdyvar(9)="snowwat"
     bdyvar(10)="graupel"
     bdyvar(11)="ice_aero"
     bdyvar(12)="liq_aero"
     bdyvar(13)="u"
     bdyvar(14)="v"
     bdyvar(15)="delp"
     bdyvar(16)="delz"
  else
     write(*,*) 'ERROR: bdy_update_type need to be 1 or 2'
     call MPI_FINALIZE(ierror)
     stop
  endif
  
  if(mype==0) then
     call fv3bk%init(fv3_io_layout_y)
     call fv3bk%setup_grid()  

     call fv3bdy%init('gfs_bndy.tile7.000.nc')
     call fv3bdy%read_bdy_ij()

     if(bdy_update_type==2) call fv3bdy%create_new_bdy('gfs_bndy.tile7.000_gsi.nc')

     do nn=1,numbdy

         if(trim(bdyvar(nn)) == "u" .or. trim(bdyvar(nn)) == "v" .or. &
            trim(bdyvar(nn)) == "uv") then
            if(bdy_update_type==2) then
               if(trim(bdyvar(nn)) == "u" ) then
                  write(6,*) 'update bdy ',trim(bdyvar(nn)),' from u'
                  call fv3bk%read_field(trim(bdyvar(nn)))  

                  call fv3bdy%read_bdy("u_s")
                  call update_bc_4side_gsiwind(fv3bdy,fv3bk,"u_s")
                  call fv3bdy%update_bdy_gsi("u_s")

                  call fv3bdy%read_bdy("u_w")
                  call update_bc_4side_gsiwind(fv3bdy,fv3bk,"u_w")
                  call fv3bdy%update_bdy_gsi("u_w")

               else if(trim(bdyvar(nn)) == "v") then
                  write(6,*) 'update bdy ',trim(bdyvar(nn)),' from v'
                  call fv3bk%read_field(trim(bdyvar(nn)))  

                  call fv3bdy%read_bdy("v_w")
                  call update_bc_4side_gsiwind(fv3bdy,fv3bk,"v_w")
                  call fv3bdy%update_bdy_gsi("v_w")

                  call fv3bdy%read_bdy("v_s")
                  call update_bc_4side_gsiwind(fv3bdy,fv3bk,"v_s")
                  call fv3bdy%update_bdy_gsi("v_s")
               endif
            elseif(trim(bdyvar(nn)) == "uv") then
               call fv3bk%setup_grid_latlon()
               if(grid_type_fv3_regional == 2) then
                  call reverse_grid_r(fv3bk%grid_lon,fv3bk%nlon+1,fv3bk%nlat+1,1)
                  call reverse_grid_r(fv3bk%grid_lat,fv3bk%nlon+1,fv3bk%nlat+1,1)
               endif
               call initial_wind_convert(fv3bk%nlon,fv3bk%nlat,fv3bk%grid_lon,fv3bk%grid_lat)
               if(allocated(ud)) deallocate(ud)
               allocate(ud(fv3bk%nlon,fv3bk%nlat+1,fv3bk%nlvl))
               call fv3bk%read_field("u")
               ud=fv3bk%field3d
               call fv3bk%read_field("v")
               call update_bc_fv3uv2earth(fv3bdy,fv3bk%nlon,fv3bk%nlat,fv3bk%nlvl,grid_reverse_flag,ud,fv3bk%field3d)

               call fv3bdy%read_bdy("u_s")
               call update_bc_4side_wind(fv3bdy,fv3bk%nlon,fv3bk%nlat+1,fv3bk%nlvl,ud,"u_s")
               call fv3bdy%update_bdy("u_s")

               call fv3bdy%read_bdy("u_w")
               call update_bc_4side_wind(fv3bdy,fv3bk%nlon,fv3bk%nlat+1,fv3bk%nlvl,ud,"u_w")
               call fv3bdy%update_bdy("u_w")

               call fv3bdy%read_bdy("v_w")
               call update_bc_4side_wind(fv3bdy,fv3bk%nlon+1,fv3bk%nlat,fv3bk%nlvl,fv3bk%field3d,"v_w")
               call fv3bdy%update_bdy("v_w")

               call fv3bdy%read_bdy("v_s")
               call update_bc_4side_wind(fv3bdy,fv3bk%nlon+1,fv3bk%nlat,fv3bk%nlvl,fv3bk%field3d,"v_s")
               call fv3bdy%update_bdy("v_s")

               if(allocated(ud)) deallocate(ud)
            else
               write(6,*) 'wrong wind variable name',trim(bdyvar(nn))
            endif
         else
            write(6,*) 'update bdy ',trim(bdyvar(nn))
            call fv3bk%read_field(trim(bdyvar(nn)))
            call fv3bdy%read_bdy(trim(bdyvar(nn)))
            call update_bc_4side(fv3bdy,fv3bk,trim(bdyvar(nn)))
            if(bdy_update_type==2) then
               call fv3bdy%update_bdy_gsi(trim(bdyvar(nn)))
            else
               call fv3bdy%update_bdy(trim(bdyvar(nn)))
            endif
          endif

     enddo ! nn

     call fv3bdy%close()
     call fv3bk%close()

     write(6,*) "=== RRFS UPDATE BC SUCCESS ==="

  endif ! mype==0

  deallocate(bdyvar)

  call MPI_FINALIZE(ierror)
!

END PROGRAM check_process_imssnow

