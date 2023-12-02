SUBROUTINE read_Lightning2cld(nlon,nlat,numitem,numlight,ybegin,yend,light_in,lightning)
!
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_NESDIS     read in lightning flash rate  
!
!   PRGMMR: Ming Hu          ORG: GSD/AMB        DATE: 2008-11-30
!
! ABSTRACT: 
!  This subroutine read in lightning flash rate
!
! PROGRAM HISTORY LOG:
!    2009-01-20  Hu  Add NCO document block
!
!
!   input argument list:
!     nlon        - no. of lons on subdomain (buffer points on ends)
!     nlat        - no. of lats on subdomain (buffer points on ends)
!     numitem     - number of item in each observation
!     numlight    - number of observation
!     ybegin      - begin Y index for the domain
!     yend        - end Y index for the domain
!
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

  INTEGER(i_kind),intent(in) :: nlon,nlat
  INTEGER(i_kind),intent(in) :: numitem
  INTEGER(i_kind),intent(in) :: numlight 
  INTEGER(i_kind),intent(in) :: ybegin,yend
  real(r_single),intent(in)  :: light_in(numitem,numlight)
  real(r_single), intent(out):: lightning(nlon,nlat)
!
!  local
!
  integer(i_kind):: ilat1s,ilon1s
  INTEGER(i_kind) :: i,j,ii,jj
!
  ilon1s=1
  ilat1s=2

  write(6,*) "process lightning obs for domain ", nlon,nlat,ybegin,yend
  if(numitem /= 4) then
     write(*,*) 'ERROR: not match expected item size ',numitem
     stop 1234
  endif

  DO i=1,numlight
    ii=int(light_in(ilon1s,i)+0.001_r_single)
    jj=int(light_in(ilat1s,i)+0.001_r_single)

    if( (ii >= 1 .and. ii <= nlon ) .and. &
        (jj >= ybegin .and. jj <= yend ) ) then
      jj=jj-ybegin+1
      lightning(ii,jj)=light_in(4,i)
    else
!      write(6,*) 'lightning obs outside analysis domain ',i,ii,jj
    endif
  ENDDO

END SUBROUTINE read_Lightning2cld
