subroutine unfill_fv3_grid2t_ldmk(b,nx,ny,landmask,           &
                                  sncovr,deltaT,i_snowt_check)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    unfill_mass_grid2t        opposite of fill_mass_grid2
!   prgmmr: parrish          org: np22                date: 2004-06-22
!
! abstract: This is the same as unfill_mass_grid2t, 
!           but only add analysis increment over land. 
!
! program history log:
!   2015-01-15  Hu      - apply the land/sea mask here for soil adjustment
!                          fields
!
!   input argument list:
!     gout     - input A-grid (reorganized for distibution to local domains)
!     gin      - preexisting input values to be added to on C-grid
!     nx,ny    - input grid dimensions
!     deltaT   - delta T between atmosphere and tsk/tslb(1)
!     i_snowT_check - input option for snow Temperature adjustment
!                     =0: input gin is not temperature
!                     =1: tsk make sure surface temperature over snow is below 0C
!                     =2: input gin is soil mositure, don't adjust over seaice
!                     =3: soil temperature
!                     =4: soilt1
!
!   output argument list:
!     gin      - output result on C grid
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_single,i_kind,r_kind
  use constants, only: partialSnowThreshold

  implicit none

  integer(i_kind), intent(in   ) :: nx,ny
  integer(i_kind), intent(in   ) :: i_snowt_check
  real(r_kind) , intent(inout) :: b(nx,ny)
  real(r_kind) , intent(in)    :: landmask(nx,ny)
  real(r_kind) , intent(in)    :: sncovr(nx,ny)
  real(r_kind) , intent(in)    :: deltaT(nx,ny)
  
  integer(i_kind) i,j
  real(r_kind) :: DTsTmax

  DTsTmax = 20.0_r_kind                             ! maximum allowed difference between Ts and T 1st level

! only add analysis increment over land
!tgs--- Landmask in UFS: 0 - water, 1 - land, 2 - ice
     do j=1,ny
        do i=1,nx
           if(abs(landmask(i,j)) < 1.0e-10_r_kind)  b(i,j)=0.0_r_kind ! water
           if(i_snowT_check==2 .and. abs(landmask(i,j)-2.0_r_kind)< 1.0e-10_r_kind)  b(i,j)=0.0_r_kind ! ice
!  don't change soil T (TSBL) under thick snow on land: snow cover > 0.99
           if(i_snowT_check==3 .and. ( abs(landmask(i,j)-1.0_r_kind)<1.0e-10_r_kind .and. sncovr(i,j) > partialSnowThreshold)) b(i,j)=0.0_r_kind
! Limit application of soil temp nudging in fine grid as follows:
!  - If cooling is indicated, apply locally only
!        if deltaT = Tskin - T(k=1) > -20K. for TSK and SOILT1  
!        if deltaT = TSLB(1) - T(k=1) > -20K. for TSLB  
! Idea:  If skin temp is already much colder than atmos temp,
!        it's useless to cool off the soil any more
!        As we also know, the repeated application will created
!          unrealistic values.
!  - Do similar for indicated soil warming, apply locally only:
!        if deltaT = Tskin - T(k=1) < 20K. for TSK and SOILT1  
!        if deltaT = TSLB(1) - T(k=1) < 20K. for TSLB  
!
           if( (i_snowT_check==1 .or. i_snowT_check==3 .or. i_snowT_check==4) ) then
              if(deltaT(i,j) < -DTsTmax .and. b(i,j) < 0.0_r_kind) & 
                  b(i,j)=0.0_r_kind
              if(deltaT(i,j) >  DTsTmax .and. b(i,j) > 0.0_r_kind) &
                  b(i,j)=0.0_r_kind
           endif
        end do
     end do
  
end subroutine unfill_fv3_grid2t_ldmk
