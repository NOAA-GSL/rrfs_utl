module mpi_mod
!$$$   module documentation block
!                .      .    .                                       .
! module:  rapid refresh module
! prgmmr:  Ming Hu             org: GSD/AMB           date: 2008-06-04
!
! abstract: 
!      This module contains code to get a namelis
!
! program history log:
!   2020-08-10 Hu           initial build
! 
! Subroutines Included:
!   sub load_namelist  - load namelist variables 
!
! Variable Definitions:
!
! attributes:
!   language: f90
!   machine:  linux cluster (wjet)
!
!$$$ end documentation block

  use kinds, only: i_kind
  use mpi

  implicit none

! set default to private
  private
 
! MPI variables
  integer :: npe, mype
!
! set subroutines to public
  public :: mpi_setup, mpi_finish
  public :: npe, mype
  public :: mpi_comm_world
  public :: mpi_integer, mpi_sum
  public :: MPI_Barrier

contains

  subroutine mpi_setup
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  init_rapidrefresh_cldsurf
! prgmmr:  Ming Hu             org: GSD/AMB           date: 2008-06-04
!
! abstract:  set defaults for RR related variables
!
! program history log:
!   2008-06-03  Hu        initial build for cloud analysis
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  liunx cluster (Wjet)
!
!$$$
!
    implicit none
    integer:: ierror

!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  end subroutine mpi_setup

  subroutine mpi_finish
    implicit none
    integer:: ierror
    call MPI_FINALIZE(ierror)
  end subroutine mpi_finish

end module mpi_mod
