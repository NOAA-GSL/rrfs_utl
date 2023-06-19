module general_sub2grid_simple_mod
!$$$   module documentation block
!                .      .    .                                       .
! module:    generic_sub2grid_simple_mod  generalized sub2grid and grid2sub
!                                         style routines based on GSI
!                                         general_sub2grid_mod
!   prgmmr: Ming Hu          org: gsl                 date: 2023-06-07
!
! abstract: 
!
! program history log:
!   2023-06-07  Hu, initial documentation
!$$$ end documentation block

   use kinds, only: r_double,i_kind,i_long,r_single,r_kind
   use mpi

   implicit none

! set default to private
   private
! set subroutines to public
   public :: general_sub2grid
   public :: general_grid2sub
   public :: general_sub2grid_create_info
   public :: general_sub2grid_destroy_info
! set passed variables to public
   public :: sub2grid_info

   interface general_sub2grid
     module procedure general_sub2grid_r_single_rank3
   end interface

   interface general_grid2sub
     module procedure general_grid2sub_r_single_rank3
   end interface

   type sub2grid_info

      integer(i_kind):: lat1=0          ! no. of lats on subdomain (no buffer)
      integer(i_kind):: lon1=0          ! no. of lons on subdomain (no buffer)
      integer(i_kind):: lat2=0          ! no. of lats on subdomain (buffer)
      integer(i_kind):: lon2=0          ! no. of lons on subdomain (buffer)
      integer(i_kind):: latlon11=0      ! no. of points on subdomain (including buffer)
      integer(i_kind):: latlon1n=0      ! latlon11*nsig
      integer(i_kind):: nlat=0          ! no. of latitudes
      integer(i_kind):: nlon=0          ! no. of longitudes
      integer(i_kind):: nsig=0          ! no. of vertical levels
      integer(i_kind):: num_fields=0    ! total number of fields/levels
      integer(i_kind):: iglobal=0       ! number of horizontal points on global grid
      integer(i_kind):: itotsub=0       ! number of horizontal points of all subdomains combined
      integer(i_kind):: kbegin_loc=0    ! starting slab index for local processor
      integer(i_kind):: kend_loc=0      ! ending slab index for local processor
      integer(i_kind),pointer :: nlevs_loc(:) => null()  ! number of active local levels( = kend_loc-kbegin_loc+1)
      integer(i_kind):: npe=0           ! total number of processors
      integer(i_kind):: mype=-1         ! local processor
      integer(i_kind):: nskip=0         ! # of processors skipped between full horizontal fields in grid mode.
      logical,pointer :: vector(:)     => null()    ! logical flag, true for vector variables
      integer(i_kind),pointer :: ilat1(:)       => null()    !  no. of lats for each subdomain (no buffer)
      integer(i_kind),pointer :: jlon1(:)       => null()    !  no. of lons for each subdomain (no buffer)
      integer(i_kind),pointer :: istart(:)      => null()    !  start lat of the whole array on each pe
      integer(i_kind),pointer :: jstart(:)      => null()    !  start lon of the whole array on each pe
      integer(i_kind),pointer :: recvcounts(:)  => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::  displs_g(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: rdispls(:)     => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: sendcounts(:)  => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: sdispls(:)     => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: ijn(:)         => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: ltosj(:)       => null()    !  lat index for reordering slab
      integer(i_kind),pointer :: ltosi(:)       => null()    !  lon index for reordering slab
      integer(i_kind),pointer :: recvcounts_s(:)=> null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::     irc_s(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::     ird_s(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::     isc_g(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::     isd_g(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer ::  displs_s(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: rdispls_s(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: sendcounts_s(:)=> null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: sdispls_s(:)   => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: ijn_s(:)       => null()    !  for mpi_alltoallv (sub2grid)
      integer(i_kind),pointer :: ltosj_s(:)     => null()    !  lat index for reordering slab
      integer(i_kind),pointer :: ltosi_s(:)     => null()    !  lon index for reordering slab
      integer(i_kind),pointer :: kbegin(:)      => null()    !  starting slab index for each processor
      integer(i_kind),pointer :: kend(:)        => null()    !  ending slab index for each processor

      logical:: lallocated = .false.

   end type sub2grid_info

   logical :: print_verbose=.false.

!  other declarations  ...

   contains

   subroutine general_sub2grid_create_info(s,mype,npe,nlat,nlon,nsig,num_fields,kbegin,kend)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    general_sub2grid_create_info populate info variable s
!   prgmmr: parrish          org: np22                date: 2010-02-12
!
! abstract: given dimensions of horizontal domain and various other 
!              information, obtain all required constants to allow
!              use of general_sub2grid_ and general_grid2sub_ and store them
!              in structure variable s.
!
! program history log:
!   2023-06-07  Hu, initial documentation
!
!   input argument list:
!     s          - structure variable, waiting for all necessary information for
!                    use with general_sub2grid and general_grid2sub.
!     inner_vars - inner index, reserved for eventually putting all ensemble
!     members 
!                    on 1st (most rapidly varying) array index.
!     nlat       - number of horizontal grid points in "latitude" direction
!     nlon       - number of horizontal grid points in "longitude"
!     nsig       - number of vertical levels for 1 3d variable.
!     num_fields - total number of 2d fields to be processed.

!   output argument list:
!     s          - structure variable, contains all necessary information for
!                    moving this set of subdomain variables sub_vars to
!                    the corresponding set of full horizontal grid variables.


!$$$
      implicit none

      type(sub2grid_info),     intent(inout) :: s
      integer(i_kind),         intent(in   ) :: mype,npe
      integer(i_kind),         intent(in   ) :: nlat,nlon,nsig,num_fields
      integer(i_kind),         intent(in   ) :: kbegin(npe),kend(npe)

      integer(i_kind) i,ierror,j,k,num_loc_groups,nextra,mm1,n,ns,npe_used,iadd
      integer(i_kind),allocatable:: idoit(:)

      s%npe=npe
      s%mype=mype
      s%nlat=nlat
      s%nlon=nlon
      s%iglobal=nlat*nlon
      s%nsig=nsig
      s%num_fields=num_fields
      if(s%lallocated) then
         call general_sub2grid_destroy_info(s)
      end if
      allocate(s%jstart(s%npe),s%istart(s%npe),s%ilat1(s%npe),s%jlon1(s%npe))
      allocate(s%ijn(s%npe),s%ijn_s(s%npe))

!      first determine subdomains
      call general_deter_subdomain(s%npe,s%mype,s%nlat,s%nlon, &
            s%lon1,s%lon2,s%lat1,s%lat2,s%ilat1,s%istart,s%jlon1,s%jstart)
      s%latlon11=s%lat2*s%lon2
      s%latlon1n=s%latlon11*s%nsig

      allocate(s%isc_g(s%npe),s%isd_g(s%npe),s%displs_g(s%npe),s%displs_s(s%npe))
      allocate(s%ird_s(s%npe),s%irc_s(s%npe))

      s%ijn=s%ilat1*s%jlon1
      s%ijn_s=(s%ilat1+2)*(s%jlon1+2)
      mm1=s%mype+1
      do i=1,s%npe
         s%irc_s(i)=s%ijn_s(mm1)
         s%isc_g(i)=s%ijn(mm1)
      end do

!        obtain ltosi,ltosj
      allocate(s%ltosi(s%nlat*s%nlon),s%ltosj(s%nlat*s%nlon))
      do i=1,s%nlat*s%nlon
         s%ltosi(i)=0
         s%ltosj(i)=0
      end do

!                       load arrays dealing with global grids
      s%isd_g(1)=0
      s%displs_g(1)=0
      do n=1,s%npe
         if(n/=1) then
            s%isd_g(n)=s%isd_g(n-1)+s%isc_g(n-1)
            s%displs_g(n)=s%displs_g(n-1)+s%ijn(n-1)
         end if
         do j=1,s%jlon1(n)
            ns=s%displs_g(n)+(j-1)*s%ilat1(n)
            do i=1,s%ilat1(n)
               ns=ns+1
               s%ltosi(ns)=s%istart(n)+i-1
               s%ltosj(ns)=s%jstart(n)+j-1
            end do
         end do
      end do

! Load arrays dealing with subdomain grids
      s%ird_s(1)=0
      s%displs_s(1)=0
      do n=1,s%npe
         if(n/=1) then
            s%ird_s(n)=s%ird_s(n-1)+s%irc_s(n-1)
            s%displs_s(n)=s%displs_s(n-1)+s%ijn_s(n-1)
         end if
      end do
! set total number of points from all subdomain grids
      s%itotsub=s%displs_s(s%npe)+s%ijn_s(s%npe)

!        obtain ltosi_s,ltosj_s
      allocate(s%ltosi_s(s%itotsub),s%ltosj_s(s%itotsub))
      do i=1,s%itotsub
         s%ltosi_s(i)=0
         s%ltosj_s(i)=0
      end do

      do n=1,s%npe
         do j=1,s%jlon1(n)+2
            ns=s%displs_s(n)+(j-1)*(s%ilat1(n)+2)
            do i=1,s%ilat1(n)+2
               ns=ns+1
               s%ltosi_s(ns)=s%istart(n)+i-2
               s%ltosj_s(ns)=s%jstart(n)+j-2
               if(s%ltosi_s(ns)==0) s%ltosi_s(ns)=1
               if(s%ltosi_s(ns)==nlat+1) s%ltosi_s(ns)=s%nlat
               if(s%ltosj_s(ns)==0) s%ltosj_s(ns)=1
               if(s%ltosj_s(ns)==nlon+1) s%ltosj_s(ns)=s%nlon
            end do
         end do
      end do  ! end do over npe

!      next, determine vertical layout:
      allocate(idoit(0:s%npe-1))
      if(s%num_fields<s%npe) then
         idoit=1
         idoit(s%num_fields-1:s%npe-1)=0
         npe_used=s%num_fields
         if(s%mype==0.and.print_verbose) &
           write(6,'(a,4I7)')'npe,num_fields,npe_used,idoit=',s%npe,s%num_fields,npe_used,idoit
      else
         idoit=0
         npe_used=0
         do n=0,s%npe-1
            npe_used=npe_used+1
            idoit(n)=1
         end do
      end if
      allocate(s%kbegin(0:s%npe),s%kend(0:s%npe-1))
      allocate(s%nlevs_loc(0:s%npe))
      s%kbegin=0
      s%kend=0
      s%nlevs_loc=0
      do n=1,s%npe
         s%kbegin(n-1)=kbegin(n)
         s%kend(n-1)=kend(n)
         if(s%kend(n-1) > 0) then
            s%nlevs_loc(n-1)=s%kend(n-1)-s%kbegin(n-1)+1
         else
            s%nlevs_loc(n-1)=0
         endif
      end do
      if(print_verbose) then
         do k=0,s%npe-1
            write(6,'(a,4I7)')' in general_sub2grid_create_info,k,kbegin,kend,nlevs_loc=', &
               k,s%kbegin(k),s%kend(k),s%nlevs_loc(k)
         end do
      end if

      s%kbegin_loc=s%kbegin(s%mype)
      s%kend_loc=s%kend(s%mype)

!         get alltoallv indices for sub2grid
      allocate(s%sendcounts(0:s%npe-1),s%sdispls(0:s%npe))
      allocate(s%recvcounts(0:s%npe-1),s%rdispls(0:s%npe))
      s%sdispls(0)=0
      do k=0,s%npe-1
         s%sendcounts(k)=s%ijn(k+1)*s%nlevs_loc(s%mype)
         s%sdispls(k+1)=s%sdispls(k)+s%sendcounts(k)
      end do
      s%rdispls(0)=0
      do k=0,s%npe-1
         s%recvcounts(k)=s%ijn(s%mype+1)*s%nlevs_loc(k)
         s%rdispls(k+1)=s%rdispls(k)+s%recvcounts(k)
      end do

!         get alltoallv indices for grid2sub
      allocate(s%sendcounts_s(0:s%npe-1),s%sdispls_s(0:s%npe))
      allocate(s%recvcounts_s(0:s%npe-1),s%rdispls_s(0:s%npe))
      s%sdispls_s(0)=0
      do k=0,s%npe-1
         s%sendcounts_s(k)=s%ijn_s(k+1)*s%nlevs_loc(s%mype)
         s%sdispls_s(k+1)=s%sdispls_s(k)+s%sendcounts_s(k)
      end do
      s%rdispls_s(0)=0
      do k=0,s%npe-1
         s%recvcounts_s(k)=s%ijn_s(s%mype+1)*s%nlevs_loc(k)
         s%rdispls_s(k+1)=s%rdispls_s(k)+s%recvcounts_s(k)
      end do

      s%lallocated=.true.
      if(print_verbose) then
         write(6,'(a)')' in general_sub2grid_create_info'
         write(6,'(a)')'   k, sendcounts_s, sdispls_s, recvcounts_s, rdispls_s'
         do k=0,s%npe-1
            write(6,'(I5,4I12)') k,s%sendcounts_s(k),s%sdispls_s(k),s%recvcounts_s(k),s%rdispls_s(k)
         end do
      end if

   end subroutine general_sub2grid_create_info


   subroutine general_sub2grid_destroy_info(s)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    general_sub2grid_destroy_info deallocate variable s
!   prgmmr: parrish          org: np22                date: 2010-02-12
!
! abstract: deallocate all components of type(sub2grid_info)  variable s.
!
! program history log:
!   2023-06-07  Hu, initial documentation
!   input argument list:
!     s          - structure variable, contains all necessary information for
!                    moving this set of subdomain variables sub_vars to
!                    the corresponding set of full horizontal grid variables.
!
!   output argument list:
!     s          - returned with all allocatable pointers pointed to NULL and
!                   structure variable s%lallocated set to .false.
!
!$$$
      implicit none

      type(sub2grid_info),     intent(inout) :: s

      if(s%lallocated) then
         deallocate(s%ilat1,s%jlon1,s%istart,s%jstart,s%recvcounts,s%displs_g)
         deallocate(s%rdispls,s%sendcounts,s%sdispls,s%ijn,s%recvcounts_s)
         deallocate(s%irc_s,s%ird_s,s%isc_g,s%isd_g,s%displs_s,s%rdispls_s,s%sendcounts_s,s%sdispls_s)
         deallocate(s%ijn_s,s%kbegin,s%kend)
         deallocate(s%ltosj,s%ltosi,s%ltosj_s,s%ltosi_s)
         s%lallocated=.false.
      end if

   end subroutine general_sub2grid_destroy_info


   subroutine general_deter_subdomain(npe,mype,nlat,nlon, &
                    lon1,lon2,lat1,lat2,ilat1,istart,jlon1,jstart)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    general_deter_subdomain_nolayout   perform domain decomposition
!   prgmmr: weiyu yang       org: np20                date: 1998-05-14
!
! abstract: Given an array of the observation computation load and
!           the number of available mpi tasks (npe), this routine 
!           decomposes the total analysis grid into npe subdomains
!
! program history log:
!   2023-06-07  weiyu yang
!
!   input argument list:
!     mype      - mpi task number
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
      implicit none

!     Declare passed variables
      integer(i_kind),intent(in   ) :: npe,mype,nlat,nlon
      integer(i_kind),intent(  out) :: lon1,lon2,lat1,lat2
      integer(i_kind),intent(  out) :: ilat1(npe),istart(npe),jlon1(npe),jstart(npe)

!     Declare local variables
      integer(i_kind) npts,nrnc,iinum,iileft,jrows,jleft,k,i,jjnum
      integer(i_kind) j,mm1,iicnt,ipts,jjleft
      integer(i_kind),dimension(npe+1):: iiend,jjend,iistart
      real(r_kind):: anperpe

!************************************************************************
!     Compute number of points on full grid and target number of
!     point per mpi task (pe)
      npts=nlat*nlon
      anperpe=float(npts)/float(npe)

!     Start with square subdomains
      nrnc=sqrt(anperpe)
      iinum=nlon/nrnc
      if(iinum==0) iinum=1
      iicnt=nlon/iinum
      iileft=nlon-iicnt*iinum
      jrows=npe/iinum
      jleft=npe-jrows*iinum

!     Adjust subdomain boundaries
      k=0
      istart=1
      jstart=1
      iistart(1)=1
      do i=1,iinum
         ipts = iicnt
         if(i <= iileft)ipts=ipts+1
         iiend(i)=iistart(i)+ipts-1
         iistart(i+1)=iiend(i)+1
         jjnum=jrows
         if(i <= jleft)jjnum=jrows+1
         do j=1,jjnum
            k=k+1
            jlon1(k)=ipts
            jstart(k)= iistart(i)
            ilat1(k)=nlat/jjnum
            jjleft=nlat-ilat1(k)*jjnum
            if(j <= jjleft)ilat1(k)=ilat1(k)+1
            if(j > 1)istart(k)=jjend(j-1)+1
            jjend(j)=istart(k)+ilat1(k)-1

            if(mype==0) &
            write(6,100) k-1,istart(k),jstart(k),ilat1(k),jlon1(k)
         end do
      end do
    100 format('general_DETER_SUBDOMAIN: task,istart,jstart,ilat1,jlon1=',6(i6,1x))


! Set number of latitude and longitude for given subdomain
      mm1=mype+1
      lat1=ilat1(mm1)
      lon1=jlon1(mm1)
      lat2=lat1+2
      lon2=lon1+2

      return

   end subroutine general_deter_subdomain

   subroutine general_grid2sub_r_single_rank3(s,grid_vars,sub_vars)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    general_sub2grid  convert from subdomains to full horizontal
! grid
!   prgmmr: parrish          org: np22                date: 2010-02-11
!
! abstract: generalized version of grid2sub--uses only gsi module kinds.
!              All information needed is contained in the structure variable
!              "s", instead of various modules.  This allows
!              for easy adaptation for any collection/ordering of variables
!              defined on subdomains, which need to be made available on
!              full horizontal grid for horizontal operations.
!              The structure variable is specified by subroutine
!              general_sub2grid_setup.
!              This version works with single precision (4-byte) real variables.
!              Output sub_vars, the desired arrays on horizontal subdomains, has
!              one 
!              halo row, for now, which is filled with zero, since for ensemble
!              use,
!              there is no need for a halo, but is easiest for now to keep it.
!              A later version will have variable number of halo rows, filled
!              with proper values.
!
! program history log:
!   2023-06-08  Hu, initial documentation
!
!   input argument list:
!     s          - structure variable, contains all necessary information for
!                    moving this set of subdomain variables sub_vars to
!                    the corresponding set of full horizontal grid variables.
!     grid_vars  - input grid values in horizontal slab mode.
!
!   output argument list:
!     sub_vars   - output grid values in vertical subdomain mode
!
! attributes:
!
!$$$

      use constants, only: zero
      implicit none

      type(sub2grid_info),intent(in   ) :: s
      real(r_single),     intent(in   ) :: grid_vars(s%nlon,s%nlat,s%kbegin_loc:s%kend_loc)
      real(r_single),     intent(  out) :: sub_vars(s%lat2,s%lon2,s%num_fields)

      real(r_single) :: temp(s%itotsub*(s%kend_loc-s%kbegin_loc+1))
      integer(i_kind) iloc,i,ii,k,n,ilat,jlon,ierror,icount
      integer(i_kind),dimension(s%npe) ::iskip
      integer(i_long) mpi_string

!     reorganize for eventual distribution to local domains
      iskip(1)=0
      do n=2,s%npe
        iskip(n)=iskip(n-1)+s%ijn_s(n-1)*(s%kend_loc-s%kbegin_loc+1)
      end do
!$omp parallel do  schedule(dynamic,1) private(n,k,i,jlon,ii,ilat,iloc,icount)
      do k=s%kbegin_loc,s%kend_loc
         icount=0
         do n=1,s%npe
            iloc=iskip(n)+(k-s%kbegin_loc)*s%ijn_s(n)
            do i=1,s%ijn_s(n)
               iloc=iloc+1
               icount=icount+1
               ilat=s%ltosi_s(icount)
               jlon=s%ltosj_s(icount)
               temp(iloc)=grid_vars(jlon,ilat,k)
            end do
         end do
      end do


      call mpi_type_contiguous(1,mpi_real4,mpi_string,ierror)
      call mpi_type_commit(mpi_string,ierror)

      call mpi_alltoallv(temp,s%sendcounts_s,s%sdispls_s,mpi_string, &
                        sub_vars,s%recvcounts_s,s%rdispls_s,mpi_string,mpi_comm_world,ierror)
      call mpi_type_free(mpi_string,ierror)

   end subroutine general_grid2sub_r_single_rank3

   subroutine general_sub2grid_r_single_rank3(s,sub_vars,grid_vars)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    general_sub2grid_r_single_rank4  convert from subdomains to
! full horizontal grid
!   prgmmr: parrish          org: np22                date: 2010-02-11
!
! abstract: generalized version of sub2grid--uses only gsi module kinds.
!              All information needed is contained in the structure variable
!              "s", instead of various modules.  This allows
!              for easy adaptation for any collection/ordering of variables
!              defined on subdomains, which need to be made available on
!              full horizontal grid for horizontal operations.
!              The structure variable is specified by subroutine
!              general_sub2grid_setup.
!              This version works with single precision (4-byte) real variables.
!              Input sub_vars, the desired arrays on horizontal subdomains, has
!              one
!              halo row, for now, which is filled with zero, since for ensemble
!              use,
!              there is no need for a halo, but is easiest for now to keep it.
!              A later version will have variable number of halo rows, filled
!              with proper values.
!
! program history log:
!   2010-02-11  parrish, initial documentation
!
!   input argument list:
!     s          - structure variable, contains all necessary information for
!                    moving this set of subdomain variables sub_vars to
!                    the corresponding set of full horizontal grid variables.
!     sub_vars   - input grid values in vertical subdomain mode (contains one
!     halo row)
!
!   output argument list:
!     grid_vars  - output grid values in horizontal slab mode.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

      implicit none

      type(sub2grid_info),intent(in   ) :: s
      real(r_single),     intent(in   ) :: sub_vars(s%lat2,s%lon2,s%num_fields)
      real(r_single),     intent(  out) :: grid_vars(s%nlon,s%nlat,s%kbegin_loc:s%kend_loc)

      real(r_single) :: sub_vars0(s%lat1,s%lon1,s%num_fields)
      real(r_single) :: work(s%itotsub*(s%kend_loc-s%kbegin_loc+1))
      integer(i_kind) iloc,iskip,i,i0,ii,j,j0,k,n,k_in,ilat,jlon,ierror,ioffset
      integer(i_long) mpi_string

!    remove halo row
!$omp parallel do  schedule(dynamic,1) private(k,j,j0,i0,i,ii)
      do k=1,s%num_fields
         do j=2,s%lon2-1
            j0=j-1
            do i=2,s%lat2-1
               i0=i-1
               sub_vars0(i0,j0,k)=sub_vars(i,j,k)
            end do
         end do
      end do
      call mpi_type_contiguous(1,mpi_real4,mpi_string,ierror)
      call mpi_type_commit(mpi_string,ierror)

      call mpi_alltoallv(sub_vars0,s%recvcounts,s%rdispls,mpi_string, &
                        work,s%sendcounts,s%sdispls,mpi_string,mpi_comm_world,ierror)

      call mpi_type_free(mpi_string,ierror)

      k_in=s%kend_loc-s%kbegin_loc+1

! Load grid_vars array in desired order
!$omp parallel do  schedule(dynamic,1) private(k,iskip,iloc,n,i,ilat,jlon,ii,ioffset)
      do k=s%kbegin_loc,s%kend_loc
         iskip=0
         iloc=0
         do n=1,s%npe
            if (n/=1) then
               iskip=iskip+s%ijn(n-1)*k_in
            end if
            ioffset=iskip+(k-s%kbegin_loc)*s%ijn(n)
            do i=1,s%ijn(n)
               iloc=iloc+1
               ilat=s%ltosi(iloc)
               jlon=s%ltosj(iloc)
               grid_vars(jlon,ilat,k)=work(i + ioffset)
            end do
         end do
      end do

   end subroutine general_sub2grid_r_single_rank3


end module general_sub2grid_simple_mod
