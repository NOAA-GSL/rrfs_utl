subroutine genqsat_2m(qsat,tsen,prsl,lat2,lon2,nsig,ice)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    genqsat
!   prgmmr: derber           org: np23                date: 1998-01-14
!
! abstract: obtain saturation specific humidity for given temperature.
!
! program history log:
!
!   input argument list:
!     tsen      - input sensibile temperature field (lat2,lon2,nsig)
!     prsl      - input layer mean pressure field (lat2,lon2,nsig)
!     lat2      - number of latitudes                              
!     lon2      - number of longitudes                             
!     nsig      - number of levels                              
!     ice       - logical flag:  T=include ice and ice-water effects,
!                 depending on t, in qsat calcuations.
!                 otherwise, compute qsat with respect to water surface
!
!   output argument list:
!     qsat      - saturation specific humidity (output)
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: xai,tmix,xb,omeps,eps,xbi,one,zero,&
       xa,psat,ttp,half,one_tenth,qmin
  implicit none

  logical                               ,intent(in   ) :: ice
  real(r_kind),dimension(lat2,lon2,nsig),intent(inout) :: qsat
  real(r_kind),dimension(lat2,lon2,nsig),intent(in   ) :: tsen,prsl
  integer(i_kind)                       ,intent(in   ) :: lat2,lon2,nsig


  integer(i_kind) k,j,i
  real(r_kind) pw,tdry,tr,es,es2
  real(r_kind) w,onep3,esmax
  real(r_kind) desidt,deswdt,dwdt,desdt,esi,esw
  real(r_kind),dimension(lat2):: mint,estmax
  integer(i_kind),dimension(lat2):: lmint
  logical:: idtupdate,idpupdate
  real(r_kind) :: qsat_local

! Declare local parameters
  real(r_kind),parameter:: r015 = 0.15_r_kind

  onep3 = 1.e3_r_kind
!  qsat=zero

  write(*,*) lat2,lon2,nsig
  write(6,*) 'tsen=',maxval(tsen),minval(tsen)
  write(6,*) 'prsl=',maxval(prsl),minval(prsl)

  do j=1,lon2
     do i=1,lat2
        mint(i)=340._r_kind
        lmint(i)=1
     end do
     do k=1,nsig
        do i=1,lat2
           if((prsl(i,j,k) < 30._r_kind .and.  &
               prsl(i,j,k) > 2._r_kind) .and.  &
               tsen(i,j,k) < mint(i))then
              lmint(i)=k
              mint(i)=tsen(i,j,k)
           end if
        end do
     end do
     do i=1,lat2
        tdry = mint(i)
        tr = ttp/tdry
        if (tdry >= ttp .or. .not. ice) then
           estmax(i) = psat * (tr**xa) * exp(xb*(one-tr))
        elseif (tdry < tmix) then
           estmax(i) = psat * (tr**xai) * exp(xbi*(one-tr))
        else
           w  = (tdry - tmix) / (ttp - tmix)
           estmax(i) =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                   + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))
        endif
     end do

     do k = 1,nsig
        do i = 1,lat2
           tdry = tsen(i,j,k)
           tr = ttp/tdry
           if (tdry >= ttp .or. .not. ice) then
              es = psat * (tr**xa) * exp(xb*(one-tr))
           elseif (tdry < tmix) then
              es = psat * (tr**xai) * exp(xbi*(one-tr))
           else
              esw = psat * (tr**xa) * exp(xb*(one-tr)) 
              esi = psat * (tr**xai) * exp(xbi*(one-tr)) 
              w  = (tdry - tmix) / (ttp - tmix)
!             es =  w * esw + (one-w) * esi
              es =  w * psat * (tr**xa) * exp(xb*(one-tr)) &
                       + (one-w) * psat * (tr**xai) * exp(xbi*(one-tr))

           endif

           pw = onep3*prsl(i,j,k)
           esmax = es
           if(lmint(i) < k)then
              esmax=0.1_r_kind*pw
              esmax=min(esmax,estmax(i))
           end if
           es2=min(es,esmax)
           qsat_local = eps * es2 / (pw - omeps * es2)
           qsat(i,j,k) = max(qmin,qsat_local)
           !qsat_local = max(qmin,qsat_local)

        end do
     end do
  end do
  return
end subroutine genqsat_2m

