module namelist_mod
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

  use kinds, only: r_kind, i_kind

  use rapidrefresh_cldsurf_mod, only: init_rapidrefresh_cldsurf, &
                            dfi_radar_latent_heat_time_period,metar_impact_radius,&
                            metar_impact_radius_lowcloud,l_gsd_terrain_match_surftobs,&
                            l_metar_impact_radius_change, &
                            metar_impact_radius_max,metar_impact_radius_min,&
                            metar_impact_radius_max_height,metar_impact_radius_min_height,&
                            l_sfcobserror_ramp_t, l_sfcobserror_ramp_q, &
                            l_pbl_pseudo_surfobst,l_pbl_pseudo_surfobsq,l_pbl_pseudo_surfobsuv,&
                            pblh_ration,pps_press_incr,l_gsd_limit_ocean_q, &
                            l_pw_hgt_adjust, l_limit_pw_innov, max_innov_pct, &
                            l_cleansnow_warmts,l_conserve_thetaV,r_cleansnow_warmts_threshold,&
                            i_conserve_thetav_iternum,l_gsd_soiltq_nudge,l_cld_bld,cld_bld_hgt, &
                            build_cloud_frac_p, clear_cloud_frac_p,       &
                            l_hydrometeor_bkio,nesdis_npts_rad, &
                            iclean_hydro_withRef,iclean_hydro_withRef_allcol, &
                            i_use_2mq4b,i_use_2mt4b,i_gsdcldanal_type,i_gsdsfc_uselist,&
                            i_lightpcp,i_sfct_gross,l_use_hydroretrieval_all,l_numconc,l_closeobs,&
                            i_coastline,i_gsdqc,qv_max_inc,ioption,l_precip_clear_only,l_fog_off,&
                            cld_bld_coverage,cld_clr_coverage,&
                            i_cloud_q_innovation,i_ens_mean,DTsTmax,&
                            i_T_Q_adjust,l_saturate_bkCloud,l_rtma3d,i_precip_vertical_check,&
                            l_qnr_from_qr,n0_rain,l_cld_uncertainty


  implicit none

! set default to private
  private
    namelist/rapidrefresh_cldsurf/dfi_radar_latent_heat_time_period, &
                                metar_impact_radius,metar_impact_radius_lowcloud,&
                                l_metar_impact_radius_change,metar_impact_radius_max,&
                                metar_impact_radius_min,metar_impact_radius_max_height,&
                                metar_impact_radius_min_height,l_gsd_terrain_match_surftobs,&
                                l_sfcobserror_ramp_t,l_sfcobserror_ramp_q, &
                                l_pbl_pseudo_surfobst,l_pbl_pseudo_surfobsq,l_pbl_pseudo_surfobsuv,&
                                pblh_ration,pps_press_incr,l_gsd_limit_ocean_q,&
                                l_pw_hgt_adjust, l_limit_pw_innov,max_innov_pct, &
                                l_cleansnow_warmts,l_conserve_thetaV,r_cleansnow_warmts_threshold,&
                                i_conserve_thetav_iternum,l_gsd_soiltq_nudge,l_cld_bld,cld_bld_hgt, &
                                build_cloud_frac_p, clear_cloud_frac_p,   &
                                nesdis_npts_rad, &
                                iclean_hydro_withRef,iclean_hydro_withRef_allcol,&
                                i_use_2mq4b,i_use_2mt4b,i_gsdcldanal_type,i_gsdsfc_uselist,&
                                i_lightpcp,i_sfct_gross,l_use_hydroretrieval_all,l_numconc,l_closeobs,&
                                i_coastline,i_gsdqc,qv_max_inc,ioption,l_precip_clear_only,l_fog_off,&
                                cld_bld_coverage,cld_clr_coverage,&
                                i_cloud_q_innovation,i_ens_mean,DTsTmax, &
                                i_T_Q_adjust,l_saturate_bkCloud,l_rtma3d,i_precip_vertical_check,&
                                l_qnr_from_qr,n0_rain,l_cld_uncertainty

    integer :: ios
    integer :: iyear,imonth,iday,ihour,iminute,isecond
    integer :: fv3_io_layout_y
    integer :: fv3sar_bg_opt
    namelist/setup/iyear,imonth,iday,ihour,iminute,isecond,fv3_io_layout_y,fv3sar_bg_opt
! set subroutines to public
  public :: load_namelist
  public :: iyear,imonth,iday,ihour,iminute,isecond
  public :: fv3_io_layout_y
  public :: fv3sar_bg_opt

contains

  subroutine load_namelist(mype)
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
    use kinds, only: i_kind 
!    use mpi_mod, only : mype
!
    implicit none
   
    integer, intent(in) :: mype
    integer:: ios

    iyear=2020
    imonth=07
    iday=29
    ihour=10
    iminute=0
    isecond=0
    fv3_io_layout_y=1
    fv3sar_bg_opt=0 

    open(11,file='gsiparm.anl')

    read(11,setup,iostat=ios)
    read(11,rapidrefresh_cldsurf,iostat=ios)
    if(ios/=0) then
      write(*,*) ' Error read rapidrefresh_cldsurf!!'
      stop 123
    endif
    close(11)

    if(mype==0) then
       write(6,setup)
       write(6,rapidrefresh_cldsurf)
    endif
    return
  end subroutine load_namelist

end module namelist_mod
