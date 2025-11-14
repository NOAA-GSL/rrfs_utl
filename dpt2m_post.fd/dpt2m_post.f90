!       read 6 hour forecast from z500 field from file 'FILE.grb2'
program dpt2m_post

use wgrib2api

implicit none
real, allocatable :: qsfc(:,:),qlvl(:,:),zlvl(:,:),ter(:,:),q2m_ori(:,:)
real, allocatable :: q2m(:,:),t2m(:,:),psfc(:,:),dpt2m(:,:),dpt2m_ori(:,:)
character(len=200) :: filenm,metadata
integer :: nx,ny,i,j,iret
real, allocatable :: land(:,:)
real :: z2m,h_ratio,weight,alpha,qv,tem

character(len=200) :: arg1,arg2,invfile

write (*,*)"Usage: dpt2m_post input_grib2 output_grib2"

call get_command_argument(1, arg1)
call get_command_argument(2, arg2)

filenm=arg1
invfile=trim(arg1)//'.inv'
iret = grb2_mk_inv(filenm,invfile)
if (iret.ne.0) stop 1

iret = grb2_inq(filenm,invfile,':SPFH:surface:',data2=qsfc, nx=nx, ny=ny)
if (iret.ne.1) stop 2                                         ! error if not 1 match
allocate(q2m(nx,ny))
allocate(dpt2m(nx,ny))

! Read in grib2
iret = grb2_inq(filenm,invfile,':SPFH:1 hybrid level:',data2=qlvl)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':SPFH:2 m above ground:',data2=q2m_ori)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':PRES:surface:',data2=psfc)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':TMP:2 m above ground:',data2=t2m)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':DPT:2 m above ground:',desc=metadata,data2=dpt2m_ori)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':HGT:1 hybrid level:',data2=zlvl)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':HGT:surface:',data2=ter)
if (iret.ne.1) stop 2                                         ! error if not 1 match

iret = grb2_inq(filenm,invfile,':LAND:surface:',data2=land)
if (iret.ne.1) stop 2                                         ! error if not 1 match


zlvl=zlvl-ter

! define 2m height
z2m=2.0

! shape parameter for exponential interploation
alpha=16.0

! initialize 2m dew point
dpt2m=dpt2m_ori

do i=1,nx
  do j=1,ny
    if (t2m(i,j) .le. 9999.9 .and. land(i,j) .gt. 0.7) then
! exponential interpolation between lowest level and surface
!      h_ratio=z2m/zlvl(i,j)
!      weight = 1 - exp(-alpha * h_ratio)
!      weight = max(0.0, min(weight, 1.0))
!      q2m(i,j) = (1-weight) * qsfc(i,j) + weight * qlvl(i,j)
! q2m limiter
      q2m(i,j) = q2m_ori(i,j)
      q2m(i,j) = min(q2m(i,j),1.05*(qlvl(i,j)))

! cap q2m between 0 and 0.04 kg/kg
      q2m(i,j) = max(0.0, min(q2m(i,j), 0.04))

      ! convert q2m to dpt2m
      qv  = max(1.0e-8,(q2m(i,j)/(1.-q2m(i,j))))
      tem = max(psfc(i,j) * qv/( 0.622+0.378 *qv), 1.0e-8)
      if (abs(log(tem/611.2)) .lt. 1e-20) then ! avoid divide by zero
        dpt2m(i,j) = 273.15
      else
        dpt2m(i,j) = 243.5/( ( 17.67 /       &
    &             log(tem/611.2) ) - 1) + 273.15 
      end if
      dpt2m(i,j) = min(dpt2m(i,j),t2m(i,j))
    end if
  enddo
enddo
    
iret = grb2_wrt(trim(arg2),filenm,1,data2=dpt2m,meta=metadata)
if (iret.ne.0) stop 4

stop
end
