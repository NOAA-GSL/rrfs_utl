module module_update_bc
!$$$   module documentation block
!                .      .    .                                       .
! module:   update boundary condition based on the analysis results
! from GSD for RR
!   prgmmr: Hu              org: gsd                date: 2022-03-18
!
! abstract: module for updating surface, soil, moisture from GSD for RR
!
! program history log:
!   2022-03-18  Hu
!
! subroutines included:
!
! Variable Definitions:

  implicit none

! set default to private
  private
! set subroutines to public
  public :: update_bc_4side
! set passed variables to public

contains

subroutine update_bc_4side(fv3bdy,fv3bk,varname)

  use kinds, only: r_kind,i_kind
  use module_io_fv3lam_bdy , only : io_fv3lam_bdy
  use module_io_fv3lam_bk , only : io_fv3lam_bk

  implicit none

  character(len=*), intent(in) :: varname
  type(io_fv3lam_bdy) :: fv3bdy
  type(io_fv3lam_bk) :: fv3bk

!
  integer :: nlon,nlat,nz
  integer :: nxi,nyj
  integer :: i,j,k
  integer :: ii,jj,kk
  integer :: is,ie,js,je
  integer :: isb,ieb,jsb,jeb
  integer :: ist,iet,jst,jet
  integer :: isl,iel,jsl,jel
  integer :: isr,ier,jsr,jer
  real,allocatable :: r2d4b(:,:)
  real,allocatable :: r2dwind(:,:)
  real :: delta
!
  nlon=fv3bk%nlon
  nlat=fv3bk%nlat
  nz=fv3bk%nlvl
!
  if(trim(varname)=='u_s') then
     nlat=nlat+1
  else if(trim(varname)=='u_w') then
     nlon=nlon+1
  elseif(trim(varname)=='v_s') then
     nlat=nlat+1
  elseif(trim(varname)=='v_w') then
     nlon=nlon+1
  elseif(trim(varname)=='zh') then
     nz=nz+1
  elseif(trim(varname)=='ps') then
     nz=1
  endif

  if(trim(varname) == 'u_w' .or. trim(varname) == 'v_w') then
!
! bottom
     nxi=fv3bdy%iw_top_bottom_size
     nyj=fv3bdy%jw_top_bottom_size
     jeb=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_w_bottom(ii)
           j=fv3bdy%j_w_bottom(jj)
           if( i==1 .and. j==1) then
              isb=ii
              jsb=jj
           endif
           if( i==nlon .and. j==1) then
              ieb=ii
           endif
        enddo
     enddo
     write(6,'(a50,10I5)') 'update u/v w bottom is,ie,js,je,nxi,nyj ',isb,ieb,jsb,jeb,nxi,nyj
! top
     nxi=fv3bdy%iw_top_bottom_size
     nyj=fv3bdy%jw_top_bottom_size
     jet=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_w_top(ii)
           j=fv3bdy%j_w_top(jj)
           if( i==1 .and. j==nlat) then
              ist=ii
              jst=jj
           endif
           if( i==nlon .and. j==nlat) then
              iet=ii
           endif
        enddo
     enddo
     write(6,'(a50,10I5)') 'update u/v w top is,ie,js,je,nxi,nyj ',ist,iet,jst,jet,nxi,nyj
! left
     nxi=fv3bdy%iw_right_left_size
     nyj=fv3bdy%jw_right_left_size
     jsl=1
     jel=nyj
     iel=nxi
     do ii=1,nxi
        i=fv3bdy%i_w_left(ii)
        if( i==1) isl=ii
     enddo
     write(6,'(a50,10I5)') 'update u/v w left is,ie,js,je,nxi,nyj ',isl,iel,jsl,jel,nxi,nyj
! right
     nxi=fv3bdy%iw_right_left_size
     nyj=fv3bdy%jw_right_left_size
     jsr=1
     jer=nyj
     ier=nxi
     do ii=1,nxi 
        i=fv3bdy%i_w_right(ii)
        if(i==nlon) isr=ii
     enddo 
     write(6,'(a50,10I5)') 'update u/v s right is,ie,js,je,nxi,nyj',isr,ier,jsr,jer,nxi,nyj

!
!  process each level
!
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'v_w') then
           r2dwind=fv3bk%field3d(:,:,k)
        elseif(trim(varname) == 'u_w') then
           allocate(r2d4b(nlon-1,nlat+1))
           r2d4b=fv3bk%field3d(:,:,k)
           call dgrid_to_cgrid_u(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
!!!!!!!!!!!!!!!!!!!!
! bottom u/v w
!
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        is=isb
        ie=ieb
        js=jsb
        je=jeb
        allocate(r2d4b(nxi,nyj))

        r2d4b=fv3bdy%bdy_bottom(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_w_bottom(ii)
              j=fv3bdy%j_w_bottom(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
! get bottom halo update
        do i=is,ie
           delta=r2d4b(i,js) !-fv3bdy%bdy_bottom(i,js,kk)
           do j=1,js-1
              r2d4b(i,j)=delta !+ r2d4b(i,j)
           enddo
        enddo
!
! the boundary wind is earch grid but the background is on native grid
!  can not caculate difference. simply replace the wind in halo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_bottom(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+ r2d4b(i,j)
           enddo
! get right halo update
           delta=r2d4b(ie,j) !-fv3bdy%bdy_bottom(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=delta !+ r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_bottom(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!!!
! top u/v w
        nxi=fv3bdy%iw_top_bottom_size
        nyj=fv3bdy%jw_top_bottom_size
        is=ist
        ie=iet
        js=jst
        je=jet
        allocate(r2d4b(nxi,nyj))

        r2d4b=fv3bdy%bdy_top(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_w_top(ii)
              j=fv3bdy%j_w_top(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
! get top halo update
        do i=is,ie
           delta=r2d4b(i,js) !-fv3bdy%bdy_top(i,js,kk)
           do j=js+1,je
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j)!-fv3bdy%bdy_top(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
! get right halo update
           delta=r2d4b(ie,j) !-fv3bdy%bdy_top(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_top(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!!
! left u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        is=isl
        ie=iel
        js=jsl
        je=jel
        allocate(r2d4b(nxi,nyj))
        r2d4b=fv3bdy%bdy_left(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_w_left(ii)
              j=fv3bdy%j_w_left(jj)
              if( (i >=1 .and. i<=nlon) .and.    &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_left(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_left(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!
! right u/v w
        nxi=fv3bdy%iw_right_left_size
        nyj=fv3bdy%jw_right_left_size
        is=isr
        ie=ier
        js=jsr
        je=jer
        allocate(r2d4b(nxi,nyj))
        r2d4b=fv3bdy%bdy_right(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_w_right(ii)
              j=fv3bdy%j_w_right(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get right halo update
           delta=r2d4b(is,j)!-fv3bdy%bdy_right(is,j,kk)
           do i=is+1,ie
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_right(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)

     enddo

     deallocate(r2dwind)

  elseif(trim(varname) == 'u_s' .or. trim(varname) == 'v_s') then
!
! bottom
     nxi=fv3bdy%is_top_bottom_size
     nyj=fv3bdy%js_top_bottom_size
     jeb=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_s_bottom(ii)
           j=fv3bdy%j_s_bottom(jj)
           if( i==1 .and. j==1) then
              isb=ii
              jsb=jj
           endif
           if( i==nlon .and. j==1) then
              ieb=ii
           endif
        enddo
     enddo
     write(6,'(a50,10I5)') 'update u/v s bottom is,ie,js,je,nxi,nyj ',isb,ieb,jsb,jeb,nxi,nyj
! top
     nxi=fv3bdy%is_top_bottom_size
     nyj=fv3bdy%js_top_bottom_size
     jet=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_s_top(ii)
           j=fv3bdy%j_s_top(jj)
           if( i==1 .and. j==nlat) then
              ist=ii
              jst=jj
           endif
           if( i==nlon .and. j==nlat) then
              iet=ii
           endif
        enddo
     enddo
     write(6,'(a50,10I5)') 'update u/v s top is,ie,js,je,nxi,nyj',ist,iet,jst,jet,nxi,nyj

! left
     nxi=fv3bdy%is_right_left_size
     nyj=fv3bdy%js_right_left_size
     jsl=1
     jel=nyj
     iel=nxi
     do ii=1,nxi
        i=fv3bdy%i_s_left(ii)
        if( i==1) isl=ii
     enddo
     write(6,'(a50,10I5)') 'update u/v s left is,ie,js,je,nxi,nyj',isl,iel,jsl,jel,nxi,nyj

! right
     nxi=fv3bdy%is_right_left_size
     nyj=fv3bdy%js_right_left_size
     jsr=1
     jer=nyj
     ier=nxi
     do ii=1,nxi
        i=fv3bdy%i_s_right(ii)
        if(i==nlon) isr=ii
     enddo
     write(6,'(a50,10I5)') 'update u/v s right is,ie,js,je,nxi,nyj',isr,ier,jsr,jer,nxi,nyj
!
!  process each level
!    
     allocate(r2dwind(nlon,nlat))

     do k=1,nz
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        if(trim(varname) == 'u_s') then
           r2dwind=fv3bk%field3d(:,:,k)
        elseif(trim(varname) == 'v_s') then
           allocate(r2d4b(nlon+1,nlat-1))
           r2d4b=fv3bk%field3d(:,:,k)
           call dgrid_to_cgrid_v(nlon,nlat,r2d4b,r2dwind)
           deallocate(r2d4b)
        else
           stop 1234
        endif
!!!!!!!!!!!!!!!!!
! bottom u/v s
!
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        is=isb
        ie=ieb
        js=jsb
        je=jeb
        allocate(r2d4b(nxi,nyj))

        r2d4b=fv3bdy%bdy_bottom(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_s_bottom(ii)
              j=fv3bdy%j_s_bottom(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
! get bottom halo update
        do i=is,ie
           delta=r2d4b(i,js) !-fv3bdy%bdy_bottom(i,js,kk)
           do j=1,js-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_bottom(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
! get right halo update
           delta=r2d4b(ie,j) !-fv3bdy%bdy_bottom(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_bottom(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!!!!!!!!
! top u/v s
        nxi=fv3bdy%is_top_bottom_size
        nyj=fv3bdy%js_top_bottom_size
        is=ist
        ie=iet
        js=jst
        je=jet
        allocate(r2d4b(nxi,nyj))
        r2d4b=fv3bdy%bdy_top(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_s_top(ii)
              j=fv3bdy%j_s_top(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
! get top halo update
        do i=is,ie
           delta=r2d4b(i,js)!-fv3bdy%bdy_top(i,js,kk)
           do j=js+1,je
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_top(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
! get right halo update
           delta=r2d4b(ie,j)!-fv3bdy%bdy_top(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_top(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!!!!!
! left u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        is=isl
        ie=iel
        js=jsl
        je=jel
        allocate(r2d4b(nxi,nyj))

        r2d4b=fv3bdy%bdy_left(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_s_left(ii)
              j=fv3bdy%j_s_left(jj)
              if( (i >=1 .and. i<=nlon) .and.    &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_left(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_left(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)
!!!!!!!!!!!!!!!
! right u/v s
        nxi=fv3bdy%is_right_left_size
        nyj=fv3bdy%js_right_left_size
        is=isr
        ie=ier
        js=jsr
        je=jer
        allocate(r2d4b(nxi,nyj))
        r2d4b=fv3bdy%bdy_right(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_s_right(ii)
              j=fv3bdy%j_s_right(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=r2dwind(i,j)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get right halo update
           delta=r2d4b(is,j) !-fv3bdy%bdy_right(is,j,kk)
           do i=is+1,ie
              r2d4b(i,j)=delta !+r2d4b(i,j)
           enddo
        enddo
!
        fv3bdy%bdy_right(:,:,kk)=r2d4b(:,:)

        deallocate(r2d4b)

     enddo

     deallocate(r2dwind)

  else
!!!!!!!!!!!!!!!!
! bottom
! start grid of the background in boundary area
     nxi=fv3bdy%i_top_bottom_size
     nyj=fv3bdy%j_top_bottom_size
     je=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_bottom(ii)
           j=fv3bdy%j_bottom(jj)
           if( i==1 .and. j==1) then
              is=ii
              js=jj
           endif
           if( i==nlon .and. j==1) then
              ie=ii
           endif
        enddo
     enddo

     write(6,'(a50,10I5)') 'update bottom is,ie,js,je,nxi,nyj ',is,ie,js,je,nxi,nyj
     allocate(r2d4b(nxi,nyj))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bdy%bdy_bottom(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_bottom(ii)
              j=fv3bdy%j_bottom(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=fv3bk%field3d(i,j,k)
              endif
           enddo
        enddo
! get bottom halo update
        do i=is,ie
           delta=r2d4b(i,js)-fv3bdy%bdy_bottom(i,js,kk)
           do j=1,js-1
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j)-fv3bdy%bdy_bottom(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
! get right halo update
           delta=r2d4b(ie,j)-fv3bdy%bdy_bottom(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        fv3bdy%bdy_bottom(:,:,kk)=r2d4b(:,:)
     enddo ! k

     deallocate(r2d4b)

!!!!!!!!!!!!!!!!
! top
! start grid of the background in boundary area
     nxi=fv3bdy%i_top_bottom_size
     nyj=fv3bdy%j_top_bottom_size
     je=nyj
     do jj=1,nyj
        do ii=1,nxi
           i=fv3bdy%i_top(ii)
           j=fv3bdy%j_top(jj)
           if( i==1 .and. j==nlat) then
              is=ii
              js=jj
           endif
           if( i==nlon .and. j==nlat) then
              ie=ii
           endif
        enddo
     enddo
     write(6,'(a50,10I5)') 'update top:',is,ie,js,je,nxi,nyj
     allocate(r2d4b(nxi,nyj))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bdy%bdy_top(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_top(ii)
              j=fv3bdy%j_top(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=fv3bk%field3d(i,j,k)
              endif
           enddo
        enddo
! get top halo update
        do i=is,ie
           delta=r2d4b(i,js)-fv3bdy%bdy_top(i,js,kk)
           do j=js+1,je
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j)-fv3bdy%bdy_top(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
! get right halo update
           delta=r2d4b(ie,j)-fv3bdy%bdy_top(ie,j,kk)
           do i=ie+1,nxi
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        fv3bdy%bdy_top(:,:,kk)=r2d4b(:,:)
     enddo ! k
     deallocate(r2d4b)

!!!!!!!!!!!!!!!!
! left
! start grid of the background in boundary area
     nxi=fv3bdy%i_right_left_size
     nyj=fv3bdy%j_right_left_size
     js=1
     je=nyj
     ie=nxi
     do ii=1,nxi
        i=fv3bdy%i_left(ii)
        if( i==1) is=ii
     enddo

     write(6,'(a50,10I5)') 'update left:',is,ie,js,je,nxi,nyj
     allocate(r2d4b(nxi,nyj))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bdy%bdy_left(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_left(ii)
              j=fv3bdy%j_left(jj)
              if( (i >=1 .and. i<=nlon) .and.    &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=fv3bk%field3d(i,j,k)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get left halo update
           delta=r2d4b(is,j)-fv3bdy%bdy_left(is,j,kk)
           do i=1,is-1
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        fv3bdy%bdy_left(:,:,kk)=r2d4b(:,:)
     enddo ! k

     deallocate(r2d4b)
!!!!!!!!!!!!!!!!
! right
! start grid of the background in boundary area
     nxi=fv3bdy%i_right_left_size
     nyj=fv3bdy%j_right_left_size
     js=1
     je=nyj
     ie=nxi
     do ii=1,nxi
        i=fv3bdy%i_right(ii)
        if(i==nlon) is=ii
     enddo
     write(6,'(a50,10I5)') 'update right:',is,ie,js,je,nxi,nyj

     allocate(r2d4b(nxi,nyj))
     do k=1,nz
! fill in all boundary area with background
        if( nz > 1) then 
           kk=k+1
        else
           kk=k
        endif
        r2d4b=fv3bdy%bdy_right(:,:,kk)
        do jj=1,nyj
           do ii=1,nxi
              i=fv3bdy%i_right(ii)
              j=fv3bdy%j_right(jj)
              if( (i >=1 .and. i<=nlon) .and.   &
                  (j >=1 .and. j<=nlat) ) then
                 r2d4b(ii,jj)=fv3bk%field3d(i,j,k)
              endif
           enddo
        enddo
!
        do j=1,nyj
! get right halo update
           delta=r2d4b(is,j)-fv3bdy%bdy_right(is,j,kk)
           do i=is+1,ie
              r2d4b(i,j)=r2d4b(i,j)+delta
           enddo
        enddo
!
        fv3bdy%bdy_right(:,:,kk)=r2d4b(:,:)
     enddo ! k

     deallocate(r2d4b)
!!!!!!!!!!
!
  endif
end subroutine update_bc_4side

subroutine dgrid_to_cgrid_u(nlon,nlat,ud,uc)
!
!-----------------------------------------------------------------------
!***  The GSI updates the D-grid winds but the BC file also needs the
!***  C-grid winds.   Compute the C-grid winds in the boundary rows
!***  by interpolating from the D-grid winds and write them into the 
!***  BC file.
!-----------------------------------------------------------------------
!
      integer,intent(in)  :: nlon,nlat
      real, intent(in)    :: ud(nlon-1,nlat+1)
      real, intent(inout) :: uc(nlon,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

      do j=1,nlat
          do i=2,nlon-1
            uc(i,j)=0.25*(ud(i-1,j)+ud(i-1,j+1)      &
                         +ud(i  ,j)+ud(i  ,j+1))
          enddo
          uc(1,j)=2.*uc(2,j)-uc(3,j)
          uc(nlon,j)=2.*uc(nlon-1,j)-uc(nlon-2,j)
      enddo

!
end subroutine dgrid_to_cgrid_u

subroutine dgrid_to_cgrid_v(nlon,nlat,vd,vc)
!
!-----------------------------------------------------------------------
!***  The GSI updates the D-grid winds but the BC file also needs the
!***  C-grid winds.   Compute the C-grid winds in the boundary rows
!***  by interpolating from the D-grid winds and write them into the 
!***  BC file.
!-----------------------------------------------------------------------
!
      integer,intent(in)  :: nlon,nlat
      real, intent(in)    :: vd(nlon+1,nlat-1)
      real, intent(inout) :: vc(nlon,nlat)
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx

      do j=2,nlat-1
          do i=1,nlon
            vc(i,j)=0.25*(vd(i  ,j-1)+vd(i  ,j)    &
                         +vd(i+1,j-1)+vd(i+1,j))
          enddo
      enddo
      do i=1,nlon
        vc(i,1)=2.*vc(i,2)-vc(i,3)
        vc(i,nlat)=2.*vc(i,nlat-1)-vc(i,nlat-2)
      enddo

!
end subroutine dgrid_to_cgrid_v


end module module_update_bc
