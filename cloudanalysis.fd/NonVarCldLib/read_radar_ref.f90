SUBROUTINE read_radar_ref(mype,lunin,istart,jstart,   &
                         nlon,nlat,Nmsclvl,numref,ref_mosaic31)
!
!
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_NESDIS     read in radar reflectivity    
!
!   PRGMMR: Ming Hu          ORG: GSD/AMB        DATE: 2006-11-30
!
! ABSTRACT: 
!  This subroutine read in radar reflectivity
!
! PROGRAM HISTORY LOG:
!    2009-01-20  Hu  Add NCO document block
!
!
!   input argument list:
!     mype        - processor ID
!     lunin       - unit in which data are read in
!     jstart      - start lon of the whole array on each pe
!     istart      - start lat of the whole array on each pe
!     nlon        - no. of lons on subdomain (buffer points on ends)
!     nlat        - no. of lats on subdomain (buffer points on ends)
!     numref      - number of observation
!
!   output argument list:
!     Nmsclvl     - vertical level of radar observation ref_mosaic31
!     ref_mosaic31- radar reflectivity horizontally in analysis grid and 
!                       vertically in mosaic grid (height)
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
  use kinds, only: r_kind,i_kind
  implicit none

  INTEGER(i_kind),intent(in) :: mype
  INTEGER(i_kind),intent(in) :: nlon,nlat
  integer(i_kind),intent(in) :: lunin
  integer(i_kind),intent(in) :: istart
  integer(i_kind),intent(in) :: jstart
  INTEGER(i_kind),intent(in) :: numref

  INTEGER(i_kind),intent(out):: Nmsclvl
  real(r_kind),   intent(out):: ref_mosaic31(nlon,nlat,31)
!
!  local 
!
  real(r_kind),allocatable :: ref_in(:,:)

  character(10) :: obstype
  integer(i_kind):: nreal,nchanl,ilat1s,ilon1s
  character(20) :: isis

  INTEGER(i_kind) :: i, ii,jj, k
  INTEGER(i_kind) :: ib,jb

!
  ib=jstart   ! begin i point of this domain
  jb=istart   ! begin j point of this domain

  read(lunin) obstype,isis,nreal,nchanl

  ilon1s=1
  ilat1s=2
  Nmsclvl = nreal - 2
  IF( Nmsclvl .ne. 21 .and. Nmsclvl .ne.31) then
     write(6,*) ' read_radar_ref: ',      &
                'vertical dimesion inconsistent when read in reflectivty mosaic'
     write(6,*) 'read in:',Nmsclvl
     write(6,*) 'need:', 21, 'or', 31
     stop 118
  ENDIF
  allocate( ref_in(nreal,numref) )
  ref_mosaic31=-9999.0_r_kind

  read(lunin)  ref_in
  DO i=1,numref
    ii=int(ref_in(ilon1s,i)+0.001_r_kind) - ib + 2
    jj=int(ref_in(ilat1s,i)+0.001_r_kind) - jb + 2
    if( ( ii >= 1 .and. ii <= nlon ) .and.   &
        ( jj >= 1 .and. jj <= nlat ) ) then
       DO k=1,Nmsclvl
         ref_mosaic31(ii,jj,k)=ref_in(2+k,i)
       ENDDO
    else
       write(6,*) 'read_radar_ref: Error ii or jj:',mype,ii,jj,i,ib,jb
    endif
  ENDDO
  deallocate(ref_in)

END SUBROUTINE read_radar_ref

SUBROUTINE read_radar_ref_bin(mype,lunin,istart,jstart,nlon,nlat,Nmsclvl,ref_mosaic31)
!
!
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_NESDIS     read in radar reflectivity    
!
!   PRGMMR: Ming Hu          ORG: GSL/AVID        DATE: 2022-02-02
!
! ABSTRACT: 
!  This subroutine read in radar reflectivity
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
!     Nmsclvl     - vertical level of radar observation ref_mosaic31
!
!   output argument list:
!     ref_mosaic31- radar reflectivity horizontally in analysis grid and 
!                       vertically in mosaic grid (height)
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
  use kinds, only: r_single,i_kind
  implicit none

  INTEGER(i_kind),intent(in) :: mype
  INTEGER(i_kind),intent(in) :: nlon,nlat
  INTEGER(i_kind),intent(in) :: istart,jstart
  integer(i_kind),intent(in) :: lunin
  INTEGER(i_kind),intent(in) :: Nmsclvl
  real(r_single), intent(inout):: ref_mosaic31(nlon,nlat,Nmsclvl)
  real(r_single), allocatable  :: ref_tmp(:,:,:)
  
!
!  local 
!
  INTEGER(i_kind) :: Nmsclvl_radar,nlon_radar,nlat_radar
  integer :: i,j,k
!
  read(lunin) Nmsclvl_radar,nlon_radar,nlat_radar
  if( Nmsclvl_radar==Nmsclvl)  then
     write(6,*) 'read reflectivity dimension=',Nmsclvl_radar,nlon_radar,nlat_radar
     write(6,*) 'read reflectivity dimension=',istart,jstart,nlon,nlat
     allocate(ref_tmp(nlon_radar,nlat_radar,Nmsclvl_radar))
     read(lunin) ref_tmp

     do j=1,nlat
       do i=1,nlon
          ref_mosaic31(i,j,:)=ref_tmp(min(max(jstart+i-2,1),nlon_radar),min(max(istart+j-2,1),nlat_radar),:)
       enddo
     enddo

     deallocate(ref_tmp)
     do k=1,Nmsclvl_radar
        write(6,*) 'ref_mosaic31=',k,maxval(ref_mosaic31(:,:,k)), &
                                     minval(ref_mosaic31(:,:,k))
     enddo
  else
     write(6,*) 'mismatch dimension in radar reflectivity obs:'
     write(6,*) 'obs=',nlon_radar,nlat_radar,Nmsclvl_radar
     write(6,*) 'input=',nlon,nlat,Nmsclvl
  endif

end SUBROUTINE read_radar_ref_bin
