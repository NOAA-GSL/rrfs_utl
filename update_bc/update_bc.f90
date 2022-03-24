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
  use module_update_bc     , only : update_bc_4side
  use module_io_fv3lam_bk  , only : io_fv3lam_bk

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
  namelist/setup/ fv3_io_layout_y
  logical :: ifexist
!
  integer, parameter :: numbdy=6
  character(len=20) :: bdyvar(numbdy)
  character(len=20) :: bkvar(numbdy)
  integer :: nn
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program

! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  bdyvar(1)="ps"
  bdyvar(2)="zh"
  bdyvar(3)="t"
  bdyvar(4)="sphum"
  bdyvar(5)="u"
  bdyvar(6)="v"
  bkvar(1)="delp"
  bkvar(2)="DZ"
  bkvar(3)="T"
  bkvar(4)="sphum"
  bkvar(5)="u"
  bkvar(6)="v"
  
!
! NCEP LSF has to use all cores allocated to run this application 
! but this if check can make sure only one core run through the real code.
  if(mype==0) then
!
!  get namelist
!
     fv3_io_layout_y=1

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

     call fv3bk%init(fv3_io_layout_y)
     call fv3bk%setup_grid()  

     call fv3bdy%init('gfs_bndy.tile7.000.nc')
     call fv3bdy%read_bdy_ij()

     do nn=1,numbdy

         if(trim(bdyvar(nn)) == "u") then
            write(6,*) 'update bdy ',trim(bdyvar(nn)),' from u'
            call fv3bk%read_field(trim(bdyvar(nn)))  

            call fv3bdy%read_bdy("u_s")
            call update_bc_4side(fv3bdy,fv3bk,"u_s")
            call fv3bdy%update_bdy("u_s")

            call fv3bdy%read_bdy("u_w")
            call update_bc_4side(fv3bdy,fv3bk,"u_w")
            call fv3bdy%update_bdy("u_w")

         else if(trim(bdyvar(nn)) == "v") then
            write(6,*) 'update bdy ',trim(bdyvar(nn)),' from v'
            call fv3bk%read_field(trim(bdyvar(nn)))  

            call fv3bdy%read_bdy("v_w")
            call update_bc_4side(fv3bdy,fv3bk,"v_w")
            call fv3bdy%update_bdy("v_w")

            call fv3bdy%read_bdy("v_s")
            call update_bc_4side(fv3bdy,fv3bk,"v_s")
            call fv3bdy%update_bdy("v_s")

         else
            write(6,*) 'update bdy ',trim(bdyvar(nn)),' from ',trim(bkvar(nn))
            call fv3bk%read_field(trim(bkvar(nn)))
            call fv3bdy%read_bdy(trim(bdyvar(nn)))
            call update_bc_4side(fv3bdy,fv3bk,trim(bdyvar(nn)))
            call fv3bdy%update_bdy(trim(bdyvar(nn)))
          endif

     enddo ! nn

     call fv3bdy%close()
     call fv3bk%close()

    write(6,*) "=== RRFS UPDATE BC SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!

END PROGRAM check_process_imssnow

