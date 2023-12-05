SUBROUTINE read_Lightning2cld(obsfile,nlon,nlat,istart,jstart,lightning,&
                              istat_lightning)
!
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_Lightning2cld     read in lightning flash rate  
!
!   PRGMMR: Ming Hu          ORG: GSL/AVID        DATE: 2022-02-02
!
! ABSTRACT: 
!  This subroutine read in lightning flash rate
!
! PROGRAM HISTORY LOG:
!    2022-02-02  Hu  Add NCO document block
!
!
!   input argument list:
!     mype        - processor ID
!     lunin       - unit in which data are read in
!     nlon        - no. of lons on subdomain (buffer points on ends)
!     nlat        - no. of lats on subdomain (buffer points on ends)
!     ybegin      - begin Y index for the domain
!     yend        - end Y index for the domain

!   output argument list:
!     lightning   - lightning flash rate in analysis grid
!
! USAGE:
!   INPUT FILES: 
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 
!   MACHINE:  Linux cluster (WJET)
!
!$$$
!
!_____________________________________________________________________
!

  use kinds, only: r_kind,i_kind, r_single
  implicit none

  character(len=*), intent(in) :: obsfile
  integer(i_kind),intent(in) :: nlon,nlat
  integer(i_kind),intent(in) :: istart,jstart

  real(r_single), intent(inout):: lightning(nlon,nlat)
  integer(i_kind),intent(inout):: istat_lightning
!
!  local
!
  real(r_single),allocatable :: lightning_in(:,:)

  integer(i_kind) :: lunin
  integer(i_kind) :: ilat1s,ilon1s
  integer(i_kind) :: i,ii,jj
  integer(i_kind) :: header1,nlon_lightning,nlat_lightning,numlight,header2,header3
  logical :: fileexist
  integer(i_kind) :: num

!
!
!
  lightning=-9999.0_r_single
  lunin=12
  istat_lightning=0
  inquire(file=trim(obsfile),exist=fileexist)
  if(fileexist) then
     write(6,*) 'read in lightning from=',trim(obsfile)
     open(lunin, file=trim(obsfile),form='unformatted',status='old')

     read(lunin) header1,nlon_lightning,nlat_lightning,numlight,&
                        header2,header3
     write(6,*) header1,nlon_lightning,nlat_lightning,numlight,header2,header3
     allocate(lightning_in(header1,numlight))
     lightning_in=-9999.0_r_single
     read(lunin) lightning_in
     close(lunin)
  else
     write(6,*) 'problem open lightning file=',trim(obsfile)
     return
  endif

  ilon1s=header2
  ilat1s=header3

  write(6,*) "process lightning obs for domain ", nlon,nlat,istart,jstart
  if(header1/= 4) then
     write(*,*) 'ERROR: not match expected item size ',header1
     stop 1234
  endif

  do i=1,numlight,max(1,numlight/5)
    write(6,'(10f10.2)') lightning_in(1:header1,i)
  enddo

  num=0
  DO i=1,numlight
    ii=int(lightning_in(ilon1s,i)-jstart+2+0.001_r_single)
    jj=int(lightning_in(ilat1s,i)-istart+2+0.001_r_single)

    if( (ii >= 1 .and. ii <= nlon ) .and. &
        (jj >= 1 .and. jj <= nlat ) ) then
      lightning(ii,jj)=lightning_in(4,i)
      num=num+1
    else
!      write(6,*) 'lightning obs outside analysis domain ',i,ii,jj
    endif
  ENDDO

  deallocate(lightning_in)

  if(num>0) istat_lightning=1
  write(6,*) 'read in ligthning number=',num,istat_lightning

END SUBROUTINE read_Lightning2cld
