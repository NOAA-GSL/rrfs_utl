module gsl_update_mod
!$$$   module documentation block
!                .      .    .                                       .
! module:    gsd_update_mod module for updating surface, soil, moisture from GSD for RR
!   prgmmr: Hu              org: gsd                date: 2012-01-12
!
! abstract: module for updating surface, soil, moisture from GSD for RR
!
! program history log:
!   2012-01-12  Hu
!   2015-01-12  Hu  fix the bug in coast proximity calculation in subdomain
!   2015-01-14  Hu  do T soil nudging over snow
!   2015-01-15  Hu  move the land/sea mask check to fine grid update step
!
! subroutines included:
!   sub gsd_update_soil_tq  - change surface and soil based on analysis increment
!
! Variable Definitions:

  implicit none

! set default to private
  private
! set subroutines to public
  public :: gsl_update_soil_tq
! set passed variables to public

contains

subroutine gsl_update_soil_tq(fv3bk)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    update_surface    change surface and soil based on analysis increment
!   prgmmr: Hu          org: GSD                date: 2011-08-18
!
! abstract:  This routine does the following things:
!              1) add lowest level t increment to T2 
! 
! 
! program history log:
!   1990-10-06  parrish - original code
!   2013-10-19  todling - metguess now holds background
!
!   input argument list:
!    tinc : first level temperature analysis increment
!    qinc : first level moisture analysis increment
!
!   output argument list:
!
!   comments:
!
! attributes:
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: zero,one,fv,one_tenth,deg2rad,pi
  use constants, only: partialSnowThreshold,t0c,qmin
  use module_bkio_fv3lam_parall , only : bkio_fv3lam

  implicit none

  type(bkio_fv3lam) :: fv3bk

! Declare passed variables
  integer(i_kind) ::  is_t,is_q
  integer ::  it   ! guess time level

! Declare local variables
  integer(i_kind)  :: gmt,nday,iyear,imonth,iday
  real(r_kind)     :: declin
  real(r_kind)     :: hrang,xxlat

  integer(i_kind) i,j,k,ier,istatus
  integer(i_kind) :: nlon,nlat,nsig_soil
  real(r_kind) :: ainc,tinct
  real(r_kind) :: coast_fac,temp,temp_fac,dts_min,tincf
  real(r_kind) :: snowthreshold
  real(r_kind) :: rhgues,sumqc
  real(r_kind) :: qs_surf
! 
  real(r_kind),allocatable,dimension(:,:) :: csza
  real(r_kind),allocatable,dimension(:,:) :: qsat_surf
  real(r_kind),dimension(:,:  ),pointer :: tinc      =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: qinc      =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: qsatg     =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_tsk   =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_soilt1=>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_qvg   =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_qcg   =>NULL()
  real(r_kind),dimension(:,:,:),pointer :: ges_tslb  =>NULL()
  real(r_kind),dimension(:,:,:),pointer :: ges_smois =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_q     =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_tsen  =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: ges_psurf =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: tsk_comp  =>NULL()
  real(r_kind),dimension(:,:  ),pointer :: landmask  =>NULL()

  
!*******************************************************************************
!
  nlon=fv3bk%nlon
  nlat=fv3bk%nlat
  nsig_soil=fv3bk%nsoil
  is_t=fv3bk%is_t
  is_q=fv3bk%is_q
  it=1

  snowthreshold=1.0e-10_r_kind
 
  tinc=>fv3bk%tinc 
  qinc=>fv3bk%qinc 

!                                                                              
!   calculation solar declination
! 
  iyear=fv3bk%iyear   
  imonth=fv3bk%imonth
  iday=fv3bk%iday
  call w3fs13(iyear,imonth,iday,nday)
  declin=deg2rad*23.45_r_kind*sin(2.0_r_kind*pi*(284+nday)/365.0_r_kind)
!
!  csza = fraction of solar constant (cos of zenith angle)
   allocate(csza(nlon,nlat))
   gmt = fv3bk%ihour 
   do j=1,nlat
     do i=1,nlon   
       hrang=15._r_kind*gmt*deg2rad + fv3bk%ges_xlon(i,j)*deg2rad - pi
       xxlat=fv3bk%ges_xlat(i,j)*deg2rad
       csza(i,j)=sin(xxlat)*sin(declin)                &
                +cos(xxlat)*cos(declin)*cos(hrang)
     end do
   end do

  if( is_t > 0) then
!     --------------------------------------------
! --- Increment top level of soil temp and snow temp
!       ONLY AT LAND POINTS according to
!       sfc temp increment.  
!     --------------------------------------------

! -- modifications to reintroduce soil temperature nudging
!   -- Stan and Tanya - 15 July 2004
!  -- allow cooling of soil only (not warming)
!      - allow only up to 1.0 K (half of negative ainc for temp)
!      - don't allow if snow water is > 6 mm

!     do it=1,nfldsig
        ges_tsk=>fv3bk%ges_tsk
        ges_soilt1=>fv3bk%ges_soilt1
        ges_qvg=>fv3bk%ges_qvg
        ges_qcg=>fv3bk%ges_qcg
        ges_smois=>fv3bk%ges_smois
        ges_tslb=>fv3bk%ges_tslb
        ges_q=>fv3bk%ges_q1
        ges_tsen=>fv3bk%ges_t1
        tsk_comp=>fv3bk%tsk_comp
        ges_psurf=>fv3bk%ges_p1

        landmask=>fv3bk%landmask

        do j=1,nlat
           do i=1,nlon
              ainc=tinc(i,j)
              coast_fac = max(0._r_kind,(fv3bk%coast_prox(i,j)-0.5_r_kind))/0.5_r_kind
              temp = ges_tsen(i,j)

! *** Increase soil adjustment by factor of 2.0 if temps are 
!       25 deg C, and 2.5 if temp is 32.5 C .
!             -- Stan B. (John, Tanya) - 31 July 2006

              temp_fac  = 1._r_kind+min (1.5_r_kind,max(0._r_kind,(temp-283._r_kind)/15._r_kind))
              dts_min = -2._r_kind
! -- Allow soil temp cooling to be up to 2.5 * 0.6 = 1.5 X 2 C ( = 3K)
              dts_min = dts_min*temp_fac*0.6_r_kind

! mhu, Jan 15,2015: move the land/sea masck check to fine grid update step
              tincf = ainc*temp_fac*coast_fac
! mhu and Tanya: Jan 14, 2015: do T soil nudging over snow
              if(nsig_soil == 9) then
! - top level soil temp
                 ges_tslb(i,j,1) = ges_tslb(i,j,1) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.6_r_kind)) 
! - 0-1 cm level -  soil temp
                 ges_tslb(i,j,2) = ges_tslb(i,j,2) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.55_r_kind))
! - 1-4 cm level -  soil temp
                 ges_tslb(i,j,3) = ges_tslb(i,j,3) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.4_r_kind))
! - 4-10 cm level -  soil temp
                 ges_tslb(i,j,4) = ges_tslb(i,j,4) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.3_r_kind))
! - 10-30 cm level -  soil temp
                 ges_tslb(i,j,5) = ges_tslb(i,j,5) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.2_r_kind))
              else
! - top level soil temp
                 ges_tslb(i,j,1) = ges_tslb(i,j,1) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.6_r_kind))
! - 0-5 cm level -  soil temp
                 ges_tslb(i,j,2) = ges_tslb(i,j,2) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.4_r_kind))
! - 5-20 cm level -  soil temp
                 ges_tslb(i,j,3) = ges_tslb(i,j,3) +   &
                                 min(3._r_kind,max(dts_min,tincf*0.2_r_kind))
              endif
              if (fv3bk%sncovr(i,j) < partialSnowThreshold) THEN
! When grid cell is partially covered with snow or snow-free - always update TSK and SOILT1
                 ges_tsk(i,j)    = ges_tsk(i,j)    + min(3._r_kind,max(dts_min,tincf*0.6_r_kind))
                 ges_soilt1(i,j) = ges_soilt1(i,j) + min(3._r_kind,max(dts_min,tincf*0.6_r_kind))
                 tsk_comp(i,j)   = tsk_comp(i,j)   + min(3._r_kind,max(dts_min,tincf*0.6_r_kind))
              else  
! grid cell is fully covered with snow
                 if(tincf < zero) then
! always adjust TSK and SOILT1 when tincf < 0 - cooling
                    ges_tsk(i,j)    = ges_tsk(i,j)    + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind))
                    ges_soilt1(i,j) = ges_soilt1(i,j) + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind))
                    tsk_comp(i,j)   = tsk_comp(i,j)   + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind))
                 else
! if ticnf > 0 - warming, then adjust snow TSK and SOILT1 only if TSK < t0c (273 K).
! If TSK > t0c(273 K) most likely due to melting process, then leave TSK and SOILT1 unchanged.
                    if(ges_tsk(i,j) < t0c ) then
                       ges_tsk(i,j)    = min(t0c,ges_tsk(i,j)    + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind)))
                       ges_soilt1(i,j) = min(t0c,ges_soilt1(i,j) + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind)))
                       tsk_comp(i,j)   = min(t0c,tsk_comp(i,j)   + min(3._r_kind,max(-2._r_kind,tincf*0.6_r_kind)))
                    endif ! tsk < 273 K
                 endif ! tincf < 0.

              endif ! sno(i,j,it) < 32
           end do
        end do
!     end do ! it

! Compute surface saturated specific humidity qsat_surf for updated skin
! temperature
       allocate(qsat_surf(nlon,nlat))
       call genqsat_2m(qsat_surf,ges_tsk,ges_psurf,nlon,nlat,1,.true.)
       write(6,*) 'qsat_surf=',maxval(qsat_surf),minval(qsat_surf)

!---------------------------------------------------------
!  Nudge soil moisture
!     Tanya S. and Stan B. - 21 July 2004 - first version
!---------------------------------------------------------

! Compute saturation specific humidity.

     if( is_q > 0) then
        ges_q=>fv3bk%ges_q1
        ges_smois=>fv3bk%ges_smois

        do j=1,nlat
           do i=1,nlon
              tinct=tinc(i,j)
              ainc=qinc(i,j)/fv3bk%qsatg(i,j)  ! analysis increment in RH

! -- update surface qvg and qcg
              qs_surf = qsat_surf(i,j)/(1. - qsat_surf(i,j)) ! mix.ratio
              if(landmask(i,j)== 1)then
              !-- land
                ges_qvg(i,j) = min(qs_surf,max(qmin,ges_qvg(i,j)+(qinc(i,j)/(1-qinc(i,j)))))
              else
              !-- water or ice
                ges_qvg(i,j) = qs_surf
              endif

              if(ges_qvg(i,j) < qs_surf) then
              !-- undersaturated
                ges_qcg(i,j) = 0.
              endif

! -- use overall limits based on k level
              ainc = max(-0.3_r_kind,min(0.3_r_kind,ainc))

! -- When background is already dry and RH increment
!      is negative (toward drier still), limit ainc further.
              rhgues=ges_q(i,j)/fv3bk%qsatg(i,j)
              if (rhgues < 0.2_r_kind .and. ainc < 0.0_r_kind ) then
                  ainc=ainc*rhgues/0.2_r_kind
              end if
              if (rhgues < 0.4_r_kind .and. ainc < 0.0_r_kind ) then
                  ainc=ainc*rhgues/0.4_r_kind
              end if

              ainc = max(-0.15_r_kind,min(0.15_r_kind,ainc))

! - Only do nudging over land, if daytime (defined as
!          cos of sun zenith angle > 0.1), and if
!          sfc temp increment is negative (meaning that
!          background sfc temp was too warm)

! -- some adjustments below to soil moisture adjustment,
!      which seems to have resulted in too much moistening
!      overall.  Stan B. - 24 Oct 04 - 04z

! mhu, Jan 15,2015: move the land/sea masck check to fine grid update step
!              if (isli(i,j,it) == 1 .and. csza(i,j) > 0.3_r_kind) then
              if (csza(i,j) > 0.3_r_kind) then
                 sumqc=0 !fv3bk%sumqc(i,j)
                 if( sumqc < 1.0e-6_r_kind) then
                 if( fv3bk%sno(i,j) < snowthreshold ) then  ! don't do the moisture adjustment if there is snow     

                    if (tinct < -0.15_r_kind) then

! - top level soil moisture
! -- mod - 3/15/13
!      increase moistening from factor of 0.2 to 0.3
                       ges_smois(i,j,1) = min (max(ges_smois(i,j,1),ges_smois(i,j,2)), &
                                   ges_smois(i,j,1) + min(0.03_r_kind,max(0._r_kind,(ainc*0.3_r_kind))))
                       ges_smois(i,j,2) = min (max(ges_smois(i,j,2),ges_smois(i,j,3)), &
                                   ges_smois(i,j,2) + min(0.03_r_kind,max(0._r_kind,(ainc*0.3_r_kind))))
                       if(nsig_soil == 9) then
                          ges_smois(i,j,3) = min (max(ges_smois(i,j,3),ges_smois(i,j,4)), &
                                   ges_smois(i,j,3) + min(0.03_r_kind,max(0._r_kind,(ainc*0.3_r_kind))))
                          ges_smois(i,j,4) = min (max(ges_smois(i,j,4),ges_smois(i,j,5)), &  
                                   ges_smois(i,j,4) + min(0.03_r_kind,max(0._r_kind,(ainc*0.2_r_kind))))
                          ges_smois(i,j,5) = min (max(ges_smois(i,j,5),ges_smois(i,j,6)), &
                                   ges_smois(i,j,5) + min(0.03_r_kind,max(0._r_kind,(ainc*0.1_r_kind))))
                       endif
! -- above logic
!     7/26/04 - 
!       previously - min was sm1_p (level 2)
!       now   - min is (max of level 1 and level 2)
!       Implication - If level 1 was already more moist than
!       level 2, don't force level 1 SM back down to level 2.
! -- mod - 5/1/05
!      Decrease moistening from factor of 0.2 to 0.1
! -- mod - 3/15/13
!      increase moistening from factor of 0.1 to 0.3
                    endif

                    if (tinct >  0.15_r_kind) then
! - top level soil moisture
! -- addition 5/1/05
!     Now also dry soil if tinc is positive (warming)
!      and the RH_inc is negative.
                        ges_smois(i,j,1) = max(0.0_r_kind,ges_smois(i,j,1) + & 
                                                  max(-0.03_r_kind,min(0._r_kind,(ainc*0.3_r_kind))))
                        ges_smois(i,j,2) = max(0.0_r_kind,ges_smois(i,j,2) + & 
                                                  max(-0.03_r_kind,min(0._r_kind,(ainc*0.3_r_kind))))
                        if(nsig_soil == 9) then
                           ges_smois(i,j,3) = max(0.0_r_kind,ges_smois(i,j,3) + & 
                                                  max(-0.03_r_kind,min(0._r_kind,(ainc*0.3_r_kind))))
                           ges_smois(i,j,4) = max(0.0_r_kind,ges_smois(i,j,4) + & 
                                                  max(-0.03_r_kind,min(0._r_kind,(ainc*0.2_r_kind))))
                           ges_smois(i,j,5) = max(0.0_r_kind,ges_smois(i,j,5) + &
                                                  max(-0.03_r_kind,min(0._r_kind,(ainc*0.1_r_kind))))


                        endif
                    END IF
                 endif  !  sno(i,j,it) < snowthreshold
                 endif  !  sumqc < 1.0e-6_r_kind
              endif
           end do
        end do
!     end do ! it
     endif ! is_q > 0
  endif

  deallocate(csza)
  return
end subroutine gsl_update_soil_tq

end module gsl_update_mod
