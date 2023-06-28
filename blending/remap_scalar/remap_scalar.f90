 subroutine main(km, npz, ncnst, ak0, bk0, Atm_ak, Atm_bk, psc, qa, zh, omga, t_in, &
                 is, ie, js, je, Atm_pt, Atm_q, Atm_delp, Atm_phis, Atm_ps)
 use ISO_FORTRAN_ENV
 use omp_lib
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
 integer,                          intent(IN)    ::  is, ie, js, je
 integer,                          intent(IN)    ::  km       ! 128
 integer,                          intent(IN)    ::  npz      ! 127
 integer,                          intent(IN)    ::  ncnst    !   7
 real(kind=8), dimension(:),       intent(IN)    ::  ak0      ! (129,)
 real(kind=8), dimension(:),       intent(IN)    ::  bk0      ! (129,)
 real(kind=8), dimension(:,:),     intent(IN)    ::  psc      ! (768, 768)
 real(kind=8), dimension(:,:,:),   intent(IN)    ::  zh       ! (768, 768, 129)
 real(kind=8), dimension(:,:,:),   intent(IN)    ::  omga     ! (768, 768, 128)
 real(kind=8), dimension(:,:,:),   intent(IN)    ::  t_in     ! (768, 768, 128)
 real(kind=8), dimension(:,:,:,:), intent(IN)    ::  qa       ! (768, 768, 128, 7)
 real(kind=8), dimension(:),       intent(IN)    ::  Atm_ak   ! (128,)
 real(kind=8), dimension(:),       intent(IN)    ::  Atm_bk   ! (128,)
 real(kind=8), dimension(:,:),     intent(IN)    ::  Atm_phis ! (768, 768)
 real(kind=8), dimension(:,:),     intent(INOUT) ::  Atm_ps   ! (768, 768)
 real(kind=8), dimension(:,:,:),   intent(INOUT) ::  Atm_delp ! (768, 768, 127)
 real(kind=8), dimension(:,:,:),   intent(INOUT) ::  Atm_pt   ! (768, 768, 127)
 real(kind=8), dimension(:,:,:,:), intent(INOUT) ::  Atm_q    ! (768, 768, 128, 7)
 real(kind=8)                                    ::  Atm_ptop
 real(kind=8)                                    ::  pst
 real(kind=8), dimension(2*km+1)                 ::  pn,gz
 real(kind=8), dimension(npz+1)                  ::  gz_fv
 real(kind=8), dimension(is:ie,js:je,npz)        ::  Atm_delz
 real(kind=8), dimension(is:ie,js:je,npz)        ::  Atm_w
 real(kind=8), dimension(is:ie,npz+1,js:je)      ::  Atm_peln
 !real(kind=8), dimension(:,:),     intent(IN)    ::  Atm_phis ! (768, 768)
 real(kind=8), dimension(is:ie,npz)              ::  dp2
 real(kind=8), dimension(is:ie,npz)              ::  qn1
 real(kind=8), dimension(is:ie,km+1)             ::  pe0
 real(kind=8), dimension(is:ie,npz+1)            ::  pe1
 real(kind=8), dimension(is:ie,km+1)             ::  pn0
 real(kind=8), dimension(is:ie,npz+1)            ::  pn1
 real(kind=8), dimension(is:ie,km)               ::  qp
 real(kind=8), dimension(is:ie,js:je)            ::  z500

 integer :: sphum,liq_wat,o3mr,ice_wat,rainwat,snowwat,graupel
 integer :: i,j,k,iq
 integer :: k2,l,itoa,m
 integer :: nwat = 6
 logical :: data_source_fv3gfs = .true.
 logical :: hydrostatic = .true.
 logical :: nggps_ic = .true.
 logical :: ncep_ic = .false.
 logical :: USE_ISOTHERMO=.false.

!Constants
 !real(kind=8), parameter:: grav= 9.80616   !< acceleration due to gravity (m/s2)
 real(kind=8), parameter:: grav  =9.80665_r8_kind   !< acceleration due to gravity (m/s2)
 real(kind=8), parameter:: rdgas = 287.05_r8_kind    !< gfs: gas constant for dry air
 real(kind=8), parameter :: rvgas = 461.5_r8_kind
 real(kind=8), parameter:: zvir =  rvgas/rdgas - 1.0_r8_kind


 k2 = max(10, km/2)
 itoa = km - npz + 1
 Atm_ptop = Atm_ak(1)

 ! This is the order in my python code
 sphum   = 1
 liq_wat = 2
 o3mr    = 3
 ice_wat = 4
 rainwat = 5
 snowwat = 6
 graupel = 7

!$OMP parallel do default(none) &
!$OMP             shared(ncnst,npz,is,ie,js,je,km,k2,ak0,bk0,psc,zh,omga,qa,z500,t_in, &
!$OMP                    sphum,liq_wat,ice_wat,rainwat,snowwat,graupel, &
!$OMP                    Atm_phis,Atm_ak,Atm_bk,Atm_ptop, &
!$OMP                    Atm_ps,Atm_delp,Atm_w,Atm_q, &
!$OMP                    Atm_pt,Atm_peln,Atm_delz, &
!$OMP                    data_source_fv3gfs,hydrostatic,nwat,ncep_ic,nggps_ic,USE_ISOTHERMO) &
!$OMP             private(l,m,pst,pn,gz,pe0,pn0,pe1,pn1,dp2,qp,qn1,gz_fv)


  do 5000 j=js,je
     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
           pn0(i,k) = log(pe0(i,k))
        enddo
     enddo

     iloop: do i=is,ie !i-loop start
        do k=1,km+1
           pn(k) = pn0(i,k)
           gz(k) = zh(i,j,k)*grav
        enddo
! Use log-p for interpolation/extrapolation
! mirror image method:
        do k=km+2, km+k2
               l = 2*(km+1) - k
           gz(k) = 2.*gz(km+1) - gz(l)
           pn(k) = 2.*pn(km+1) - pn(l)
        enddo
        do k=km+k2-1, 2, -1
           if( Atm_phis(i,j).le.gz(k) .and. Atm_phis(i,j).ge.gz(k+1) ) then
              pst = pn(k) + (pn(k+1)-pn(k))*(gz(k)-Atm_phis(i,j))/(gz(k)-gz(k+1))
              go to 123
           endif
        enddo
123     Atm_ps(i,j) = exp(pst)

 ! ------------------
 ! Find 500-mb height
 ! ------------------
        pst = log(500.e2)
        do k=km+k2-1, 2, -1
           if( pst.le.pn(k+1) .and. pst.ge.pn(k) ) then
              z500(i,j) = (gz(k+1) + (gz(k)-gz(k+1))*(pn(k+1)-pst)/(pn(k+1)-pn(k)))/grav
              go to 124
           endif
        enddo
124     continue
     enddo iloop  ! i-loop

     do i=is,ie
        pe1(i,1) = Atm_ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
        do i=is,ie
           pe1(i,k) = Atm_ak(k) + Atm_bk(k)*Atm_ps(i,j)
           pn1(i,k) = log(pe1(i,k))
        enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           dp2(i,k) = pe1(i,k+1) - pe1(i,k)
           Atm_delp(i,j,k) = dp2(i,k)
        enddo
     enddo

! map tracers
     tracers: do iq=1,1 !ncnst
        if (floor(qa(is,j,1,iq)) > -999) then !skip missing scalars
           do k=1,km
              do i=is,ie
                 qp(i,k) = qa(i,j,k,iq)
              enddo
           enddo
           call mappm(km, pe0, qp, npz, pe1,  qn1, is,ie, 0, 8, Atm_ptop)
           if ( iq==sphum ) then
              call fillq(ie-is+1, npz, 1, qn1, dp2)
           else
              call fillz(ie-is+1, npz, 1, qn1, dp2)
           endif
           do k=1,npz
              do i=is,ie
                 Atm_q(i,j,k,iq) = qn1(i,k)
              enddo
           enddo
        endif
     enddo tracers

!---------------------------------------------------
! Retrive temperature using  geopotential height from external data
!---------------------------------------------------
   iloop2: do i=is,ie
! Make sure FV3 top is lower than GFS; can not do extrapolation above the top at this point
      if ( pn1(i,1) .lt. pn0(i,1) ) then
           write(*,*) 'pn1,pn0',pn1(i,1),pn0(i,1)
           write(*,*) 'FV3 top higher than external data'
           stop
      endif

      do k=1,km+1
         pn(k) = pn0(i,k)
         gz(k) = zh(i,j,k)*grav
      enddo
!-------------------------------------------------
      do k=km+2, km+k2
         l = 2*(km+1) - k
         gz(k) = 2.*gz(km+1) - gz(l)
         pn(k) = 2.*pn(km+1) - pn(l)
      enddo
!-------------------------------------------------

      gz_fv(npz+1) = Atm_phis(i,j)

      m = 1

      do k=1,npz
! Searching using FV3 log(pe): pn1
!#ifdef USE_ISOTHERMO
if(USE_ISOTHERMO) then
         do l=m,km
            if ( (pn1(i,k).le.pn(l+1)) .and. (pn1(i,k).ge.pn(l)) ) then
                gz_fv(k) = gz(l) + (gz(l+1)-gz(l))*(pn1(i,k)-pn(l))/(pn(l+1)-pn(l))
                goto 555
            elseif ( pn1(i,k) .gt. pn(km+1) ) then
! Isothermal under ground; linear in log-p extra-polation
                gz_fv(k) = gz(km+1) + (gz_fv(npz+1)-gz(km+1))*(pn1(i,k)-pn(km+1))/(pn1(i,npz+1)-pn(km+1))
                goto 555
            endif
         enddo
else
         do l=m,km+k2-1
            if ( (pn1(i,k).le.pn(l+1)) .and. (pn1(i,k).ge.pn(l)) ) then
                gz_fv(k) = gz(l) + (gz(l+1)-gz(l))*(pn1(i,k)-pn(l))/(pn(l+1)-pn(l))
                goto 555
            endif
         enddo
endif
555   m = l
      enddo

      do k=1,npz+1
         Atm_peln(i,k,j) = pn1(i,k)
      enddo

!----------------------------------------------------
! Compute true temperature using hydrostatic balance
!----------------------------------------------------
      !if (.not. data_source_fv3gfs .or. .not. present(t_in)) then
      !if (.not. data_source_fv3gfs ) then
      if (data_source_fv3gfs) then
        do k=1,npz
           Atm_pt(i,j,k) = (gz_fv(k)-gz_fv(k+1))/( rdgas*(pn1(i,k+1)-pn1(i,k))*(1.+zvir*Atm_q(i,j,k,sphum)) )
        enddo
!------------------------------
! Remap input T logarithmically in p.
!------------------------------
      else
        do k=1,km
            qp(i,k) = t_in(i,j,k)
        enddo

        call mappm(km, log(pe0), qp, npz, log(pe1), qn1, is,ie, 2, 4, Atm_ptop) ! pn0 and pn1 are higher-precision
                                                                                ! and cannot be passed to mappm
        do k=1,npz
            Atm_pt(i,j,k) = qn1(i,k)
        enddo
      endif
      if ( .not. hydrostatic ) then
         do k=1,npz
            Atm_delz(i,j,k) = (gz_fv(k+1) - gz_fv(k)) / grav
         enddo
      endif

   enddo  iloop2  ! i-loop

!-----------------------------------------------------------------------
! seperate cloud water and cloud ice from Jan-Huey Chen's HiRAM code
! only use for NCEP IC and GFDL microphy
!-----------------------------------------------------------------------
   if (.not. data_source_fv3gfs) then
      if ((nwat .eq. 3 .or. nwat .eq. 6) .and. (ncep_ic .or. nggps_ic)) then
         do k=1,npz
            do i=is,ie

               qn1(i,k) = Atm_q(i,j,k,liq_wat)
               !if (cld_amt .gt. 0) Atm_q(i,j,k,cld_amt) = 0.

               if ( Atm_pt(i,j,k) > 273.16 ) then       ! > 0C all liq_wat
                  Atm_q(i,j,k,liq_wat) = qn1(i,k)
                  Atm_q(i,j,k,ice_wat) = 0.
!#ifdef ORIG_CLOUDS_PART
               else if ( Atm_pt(i,j,k) < 258.16 ) then  ! < -15C all ice_wat
                  Atm_q(i,j,k,liq_wat) = 0.
                  Atm_q(i,j,k,ice_wat) = qn1(i,k)
               else                                     ! between -15~0C: linear interpolation
                  Atm_q(i,j,k,liq_wat) = qn1(i,k)*((Atm_pt(i,j,k)-258.16)/15.)
                  Atm_q(i,j,k,ice_wat) = qn1(i,k) - Atm_q(i,j,k,liq_wat)
               endif
!#else
!               else if ( Atm_pt(i,j,k) < 233.16 ) then  ! < -40C all ice_wat  Atm_q(i,j,k,liq_wat) = 0.
!                  Atm_q(i,j,k,ice_wat) = qn1(i,k)
!               else
!                  if ( k.eq.1 ) then  ! between [-40,0]: linear interpolation
!                     Atm_q(i,j,k,liq_wat) = qn1(i,k)*((Atm_pt(i,j,k)-233.16)/40.)
!                     Atm_q(i,j,k,ice_wat) = qn1(i,k) - Atm_q(i,j,k,liq_wat)
!                  else
!                     if (Atm_pt(i,j,k)<258.16 .and. Atm_q(i,j,k-1,ice_wat)>1.e-5 ) then
!                        Atm_q(i,j,k,liq_wat) = 0.
!                        Atm_q(i,j,k,ice_wat) = qn1(i,k)
!                     else  ! between [-40,0]: linear interpolation
!                        Atm_q(i,j,k,liq_wat) = qn1(i,k)*((Atm_pt(i,j,k)-233.16)/40.)
!                        Atm_q(i,j,k,ice_wat) = qn1(i,k) - Atm_q(i,j,k,liq_wat)
!                     endif
!                  endif
!               endif
!
!#endif
               if (nwat .eq. 6) then ! no need to check for nwat=7 (hail) since only nwat=3,6 treated here
                  Atm_q(i,j,k,rainwat) = 0.
                  Atm_q(i,j,k,snowwat) = 0.
                  Atm_q(i,j,k,graupel) = 0.
                  call mp_auto_conversion(Atm_q(i,j,k,liq_wat), Atm_q(i,j,k,rainwat),  &
                       Atm_q(i,j,k,ice_wat), Atm_q(i,j,k,snowwat) )
               endif
            enddo
         enddo
      endif
  endif ! data source /= FV3GFS GAUSSIAN NEMSIO/NETCDF and GRIB2 FILE

! For GFS spectral input, omega in pa/sec is stored as w in the input data so actual w(m/s) is calculated
! For GFS nemsio input, omega is 0, so best not to use for input since boundary data will not exist for w
! For FV3GFS NEMSIO input, w is already in m/s (but the code reads in as omga) and just needs to be remapped
!-------------------------------------------------------------
! map omega or w
!------- ------------------------------------------------------
   if ( (.not. hydrostatic) .and. (.not. ncep_ic) ) then
      do k=1,km
         do i=is,ie
            qp(i,k) = omga(i,j,k)
         enddo
      enddo
      call mappm(km, pe0, qp, npz, pe1, qn1, is,ie, -1, 4, Atm_ptop)
    if (data_source_fv3gfs) then
      do k=1,npz
         do i=is,ie
            Atm_w(i,j,k) = qn1(i,k)
         enddo
      enddo
    else
      do k=1,npz
         do i=is,ie
            Atm_w(i,j,k) = qn1(i,k)/Atm_delp(i,j,k)*Atm_delz(i,j,k)
         enddo
      enddo
     endif
   endif

   5000 continue

 end subroutine main

!-------------------------------------------------------------------------------------------------
!>@brief The subroutine 'mappm' is a general-purpose routine for remapping
!! one set of vertical levels to another.
 subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds

! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real(kind=8), intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1) !< pe1: pressure at layer edges from model top to bottom
                                                      !!      surface in the ORIGINAL vertical coordinate
                                                      !< pe2: pressure at layer edges from model top to bottom
                                                      !!      surface in the NEW vertical coordinate
! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
 real(kind=8), intent(in )::  q1(i1:i2,km)
 real(kind=8), intent(out)::  q2(i1:i2,kn)
 real(kind=8), intent(IN) ::  ptop
! local
      real(kind=8) qs(i1:i2)
      real(kind=8) dp1(i1:i2,km)
      real(kind=8) a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real(kind=8) pl, pr, tt, delp, qsum, dpsum, esl, r3, r23
      logical :: NGGPS_SUBMITTED=.true.

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
!#ifdef NGGPS_SUBMITTED
if(NGGPS_SUBMITTED) then
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo
!#endif
endif

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k) .le. pe1(i,1)) then
! above old ptop
            q2(i,k) = q1(i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
!#ifdef NGGPS_SUBMITTED
if(NGGPS_SUBMITTED) then
            q2(i,k) = a4(3,i,km)   ! this is not good.
!#else
else
            q2(i,k) = q1(i,km)
!#endif
endif
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
!#ifdef NGGPS_SUBMITTED
if(NGGPS_SUBMITTED) then
           qsum = qsum + delp * a4(3,i,km)    ! not good.
!#else
else
           qsum = qsum + delp * q1(i,km)
!#endif
endif
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm

 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: iv      !< iv =-1: winds
                               !< iv = 0: positive definite scalars
                               !< iv = 1: others
 integer, intent(in):: kord
 real(kind=8), intent(in)   ::   qs(i1:i2)
 real(kind=8), intent(in)   :: delp(i1:i2,km)     !< layer pressure thickness
 real(kind=8), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real(kind=8) gam(i1:i2,km)
 real(kind=8)   q(i1:i2,km+1)
 real(kind=8)  d4(i1:i2)
 real(kind=8)  bet, a_bot, grat
 real(kind=8)  pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im


 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
       do i=i1,i2
          if ( ext5(i,k) ) then  ! c90_mp122
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
                  a4(2,i,k) = a4(1,i,k)
                  a4(3,i,k) = a4(1,i,k)
             endif
          elseif( ext6(i,k) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile

 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)
 use ISO_FORTRAN_ENV

! !INPUT PARAMETERS:
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
 integer, intent(in):: iv      !< iv =-1: winds
                               !! iv = 0: positive definite scalars
                               !! iv = 1: others
                               !! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      !< Starting longitude
 integer, intent(in):: i2      !< Finishing longitude
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: kord    !< Order (or more accurately method no.):
                               !!
 real(kind=8), intent(in):: delp(i1:i2,km)     !< layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real(kind=8), intent(inout):: a4(4,i1:i2,km)  !< Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real(kind=8)   dc(i1:i2,km)
      real(kind=8)   h2(i1:i2,km)
      real(kind=8) delq(i1:i2,km)
      real(kind=8)  df2(i1:i2,km)
      real(kind=8)   d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real(kind=8) fac
      real(kind=8) a1, a2, c1, c2, c3, d1, d2
      real(kind=8) qm, dq, lac, qmp, pmp
      logical :: BOT_MONO=.true.

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

!#ifdef BOT_MONO
if(BOT_MONO) then
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile

 subroutine cs_limiters(im, extm, a4, iv)
 use ISO_FORTRAN_ENV
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real(kind=8), intent(inout) :: a4(4,im)   !< PPM array
! LOCAL VARIABLES:
 real(kind=8) da1, da2, a6da, r12
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters

 subroutine ppm_limiters(dm, a4, itot, lmt)
 use ISO_FORTRAN_ENV

! INPUT PARAMETERS:
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
      real(kind=8), intent(in):: dm(*)     !< Linear slope
      integer, intent(in) :: itot      !< Total Longitudes
      integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
                                       !< (Lin) 2: Positive definite constraint
                                       !< 3: do nothing (return immediately)
! INPUT/OUTPUT PARAMETERS:
      real(kind=8), intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
! LOCAL VARIABLES:
      real(kind=8) qmp, r12
      real(kind=8) da1, da2, a6da
      real(kind=8) fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters

 subroutine mp_auto_conversion(ql, qr, qi, qs)
 use ISO_FORTRAN_ENV
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
 real(kind=8), intent(inout):: ql, qr, qi, qs
 real(kind=8), parameter:: qi0_max = 2.0e-3
 real(kind=8), parameter:: ql0_max = 2.5e-3

! Convert excess cloud water into rain:
  if ( ql > ql0_max ) then
       qr = ql - ql0_max
       ql = ql0_max
  endif 
! Convert excess cloud ice into snow:
  if ( qi > qi0_max ) then
       qs = qi - qi0_max
       qi = qi0_max
  endif

 end subroutine mp_auto_conversion

 subroutine fillq(im, km, nq, q, dp)
 use ISO_FORTRAN_ENV
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
   integer,  intent(in):: im            !< No. of longitudes
   integer,  intent(in):: km            !< No. of levels
   integer,  intent(in):: nq            !< Total number of tracers
   real(kind=8), intent(in)::  dp(im,km)       !< pressure thickness
   real(kind=8), intent(inout) :: q(im,km,nq)  !< tracer mixing ratio
! !LOCAL VARIABLES:
   integer i, k, ic, k1

   do ic=1,nq
! Bottom up:
      do k=km,2,-1
         k1 = k-1
         do i=1,im
           if( q(i,k,ic) < 0. ) then
               q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
               q(i,k ,ic) = 0.
           endif
         enddo
      enddo
! Top down:
      do k=1,km-1
         k1 = k+1
         do i=1,im
            if( q(i,k,ic) < 0. ) then
                q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                q(i,k ,ic) = 0.
            endif
         enddo
      enddo

   enddo

 end subroutine fillq

!>@brief The subroutine 'fillz' is for mass-conservative filling of nonphysical negative values in the tracers. 
!>@details This routine takes mass from adjacent cells in the same column to fill negatives, if possible.
 subroutine fillz(im, km, nq, q, dp)
 use ISO_FORTRAN_ENV
 integer, parameter :: r8_kind = selected_real_kind(15) ! 15 decimal digits
   integer,  intent(in):: im                !< No. of longitudes
   integer,  intent(in):: km                !< No. of levels
   integer,  intent(in):: nq                !< Total number of tracers
   real(kind=8), intent(in)::  dp(im,km)           !< pressure thickness
   real(kind=8), intent(inout) :: q(im,km,nq)      !< tracer mixing ratio
! LOCAL VARIABLES:
   logical:: zfix(im)
   real(kind=8)::  dm(km)
   integer i, k, ic !, k1
   real(kind=8) qup, qly, dup, dq, sum0, sum1, fac
   logical :: DEV_GFS_PHYS=.true.

   do ic=1,nq
!#ifdef DEV_GFS_PHYS
if(DEV_GFS_PHYS) then
! Bottom up:
      do k=km,2,-1
         k1 = k-1
         do i=1,im
           if( q(i,k,ic) < 0. ) then
               q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
               q(i,k ,ic) = 0.
           endif
         enddo
      enddo
! Top down:
      do k=1,km-1
         k1 = k+1
         do i=1,im
            if( q(i,k,ic) < 0. ) then
                q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                q(i,k ,ic) = 0.
            endif
         enddo
      enddo
else
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0. ) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      zfix(:) = .false.
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
             zfix(i) = .true.
             if ( q(i,k-1,ic) > 0. ) then
! Borrow from above
                dq = min ( q(i,k-1,ic)*dp(i,k-1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k-1,ic) = q(i,k-1,ic) - dq/dp(i,k-1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
             if ( q(i,k,ic)<0.0 .and. q(i,k+1,ic)>0. ) then
! Borrow from below:
                dq = min ( q(i,k+1,ic)*dp(i,k+1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k+1,ic) = q(i,k+1,ic) - dq/dp(i,k+1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
          endif
         enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic)<0. .and. q(i,k-1,ic)>0.) then
             zfix(i) = .true.
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min(qly, qup)
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
          endif
      enddo

! Perform final check and non-local fix if needed
      do i=1,im
         if ( zfix(i) ) then

           sum0 = 0.
           do k=2,km
              dm(k) = q(i,k,ic)*dp(i,k)
              sum0 = sum0 + dm(k)
           enddo

           if ( sum0 > 0. ) then
             sum1 = 0.
             do k=2,km
                sum1 = sum1 + max(0., dm(k))
             enddo
             fac = sum0 / sum1
             do k=2,km
                q(i,k,ic) = max(0., fac*dm(k)/dp(i,k))
             enddo
           endif

         endif
      enddo
endif

   enddo
 end subroutine fillz
