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
  use module_bkio_fv3lam, only : bkio_fv3lam
  use constants, only : init_constants,init_constants_derived
  use gsl_update_mod, only: gsl_update_soil_tq

  implicit none
! MPI variables
  type(bkio_fv3lam) :: fv3bk
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! namelist
  integer :: iyear,imonth,iday,ihour,iminute
  integer :: fv3_io_layout_y
  integer :: is_t,is_q
  namelist/setup/ fv3_io_layout_y,iyear,imonth,iday,ihour,iminute,&
                  is_t,is_q
  logical :: ifexist

  real(r_kind),allocatable,dimension(:,:)  :: deltaT
  real(r_kind),allocatable,dimension(:,:)  :: delta
  real(r_kind),allocatable,dimension(:,:,:)  :: ges_smois
  real(r_kind),allocatable,dimension(:,:,:)  :: ges_tslb
  real(r_kind),allocatable,dimension(:,:)  :: ges_tsk
  real(r_kind),allocatable,dimension(:,:)  :: ges_soilt1
  integer :: k
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program

! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
! NCEP LSF has to use all cores allocated to run this application 
! but this if check can make sure only one core run through the real code.
  if(mype==0) then
!
!  get namelist
!
     fv3_io_layout_y=1
     iyear=2022
     imonth=1
     iday=1
     ihour=0
     iminute=0
     is_t=1
     is_q=1

     inquire(file='namelist.soiltq', EXIST=ifexist )
     if(ifexist) then
       open(10,file='namelist.soiltq',status='old')
          read(10,setup)
       close(10)
       write(*,*) 'Namelist setup are:'
       write(*,setup)
     else
       write(*,*) 'No namelist file exist, use default values'
       write(*,*) iyear,imonth,iday,ihour,iminute
     endif

     call init_constants(.true.)
     call init_constants_derived
!
     call fv3bk%init(fv3_io_layout_y,iyear,imonth,iday,ihour,iminute)
     call fv3bk%setup_grid()
     call fv3bk%read_ges()
!
     allocate(ges_smois(fv3bk%nlon,fv3bk%nlat,fv3bk%nsoil))
     allocate(ges_tslb(fv3bk%nlon,fv3bk%nlat,fv3bk%nsoil))
     allocate(ges_tsk(fv3bk%nlon,fv3bk%nlat))
     allocate(ges_soilt1(fv3bk%nlon,fv3bk%nlat))
     ges_smois=fv3bk%ges_smois
     ges_tslb=fv3bk%ges_tslb
     ges_tsk=fv3bk%ges_tsk
     ges_soilt1=fv3bk%ges_soilt1
!
     call gsl_update_soil_tq(fv3bk)
!
     allocate(deltaT(fv3bk%nlon,fv3bk%nlat))
     allocate(delta(fv3bk%nlon,fv3bk%nlat))
     deltaT=fv3bk%ges_t1-fv3bk%ges_tsk
     do k=1,fv3bk%nsoil
        delta(:,:)=fv3bk%ges_smois(:,:,k)-ges_smois(:,:,k)
        call unfill_fv3_grid2t_ldmk(delta,fv3bk%nlon,fv3bk%nlat, &
                                    fv3bk%landmask,fv3bk%sno,deltaT,2)
        fv3bk%ges_smois(:,:,k)=ges_smois(:,:,k)+delta(:,:)

        delta(:,:)=fv3bk%ges_tslb(:,:,k)-ges_tslb(:,:,k)
        call unfill_fv3_grid2t_ldmk(delta,fv3bk%nlon,fv3bk%nlat, &
                                    fv3bk%landmask,fv3bk%sno,deltaT,3)
        fv3bk%ges_tslb(:,:,k)=ges_tslb(:,:,k)+delta(:,:)
     enddo

     delta(:,:)=fv3bk%ges_tsk(:,:)-ges_tsk(:,:)
     call unfill_fv3_grid2t_ldmk(delta,fv3bk%nlon,fv3bk%nlat, &
                                 fv3bk%landmask,fv3bk%sno,deltaT,1)
     fv3bk%ges_tsk(:,:)=ges_tsk(:,:)+delta(:,:)

     delta(:,:)=fv3bk%ges_soilt1(:,:)-ges_soilt1(:,:)
     call unfill_fv3_grid2t_ldmk(delta,fv3bk%nlon,fv3bk%nlat, &
                                 fv3bk%landmask,fv3bk%sno,deltaT,4)
     fv3bk%ges_soilt1(:,:)=ges_soilt1(:,:)+delta(:,:)

     deallocate(ges_soilt1)
     deallocate(ges_tsk)
     deallocate(ges_tslb)
     deallocate(ges_smois)
     deallocate(delta)
     deallocate(deltaT)

     call fv3bk%update_soil()
     call fv3bk%close
!
    write(6,*) "=== RRFS ADJUST SOIL T/Q SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!

END PROGRAM check_process_imssnow

