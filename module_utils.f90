! =============================================================================
! Modules holding variables needed for running the CRTM
! Written by Man-Yau Chan
! =============================================================================
! Module index:
! 1) namelist_mod:  Variables and subroutines to load namelist
! 2) constants:     Some useful constants
! 2) mpi_module:    Variables and subroutines for parallelization.
! 3) wrf_utils:     Variables and subroutines to read essential WRF fields for
!                   running CRTM. Also contains subroutines to divvy up the
!                   WRF fields into "slabs" among the processes.
! 4) crtm_utils:    Variables and subroutines to evoke CRTM.


! =============================================================================
! Module to read namelist variables
! =============================================================================
module namelist_mod

   implicit none
   !--- Parameters associated running the IR sensors
   ! Total number of channels to run CRTM for (up to 99)
   ! Suppose 2 satellites, A and B, and we want channels 1 and 2 from A, and
   ! channels 6 and 9 from B. Then, channel_count = 2 + 2 = 4
   integer                            :: channel_count
   ! Name of the sensor corresponding to each channel.
   ! Suppose sensor on A is SA, and sensor on B is SB, then
   ! sensor_list = 'SA', 'SB'
   character(len=50), dimension( 99 ) :: sensor_list
   ! Channel numbers for each of the specified channels
   ! For given example, channel_id = 1,2,6,9
   integer          , dimension( 99 ) :: channel_id
   ! Physical heights of the channels in quesiton.
   ! sensor_heights = 35786000.0, 35786000.0, 35786000.0, 35786000.0
   real             , dimension( 99 ) :: sensor_heights
   ! Satellite longitudes in question
   real             , dimension( 99 ) :: sensor_lons
   ! Geographical limits of CRTM calculation for each channel. Must be s.t.
   ! the maximum abs zenith angle is 80
   real             , dimension( 99 ) :: sensor_lon_min, sensor_lon_max
   real             , dimension( 99 ) :: sensor_lat_min, sensor_lat_max
   ! Directory to the CRTM coefficients
   CHARACTER(256):: crtm_coeff_path

  ! --- Parameters associated with the WRF file
  ! Thinning of the WRF grid before doing CRTM
  integer :: grid_thinning
  ! Starting index of point to run the CRTM
  integer :: first_index


  !-- Namelist contents :
  namelist / crtm_parameters / channel_count, sensor_list, channel_id, &
                               sensor_heights, sensor_lons, sensor_lon_min, &
                               sensor_lon_max, sensor_lat_min, sensor_lat_max, &
                               crtm_coeff_path
  namelist / wrf_parameters / grid_thinning, first_index


  CONTAINS

  ! Namelist reader subroutine
  ! ---------------------------
  subroutine  read_namelist()

    implicit none
    !  Local scalars:
    character(len=80)   :: namelist_file       ! Input namelist filename.
    integer, parameter  :: namelist_unit=7     ! Input namelist unit.
    integer             :: iost                 ! Error code.
    integer             :: i

    !-- initialize crtm parameters
    channel_count = 99
    sensor_list = '     '
    channel_id = -99
    sensor_lons = -999
    sensor_heights = 35786000.0
    sensor_lon_min = -999
    sensor_lon_max = -999
    sensor_lat_min = -999
    sensor_lat_max = -999
    crtm_coeff_path = './coefficients'
    ! -- initialize wrf parameters
    grid_thinning = 1
    first_index = 1


    !  read namelist
    !-----------------
    namelist_file = 'namelist.crtm'
    iost = 0
    ! Open namelist file
    open ( file = namelist_file, unit = namelist_unit,                   &
           status = 'old', access = 'sequential', iostat = iost )
    if( iost .ne. 0 ) then
        write(*,*)'namelist.crtm does not exist, please check it.'
        stop 'read_namelist'
    endif

    ! Load crtm parameters
    iost = 0
    read ( unit = namelist_unit, nml = crtm_parameters, iostat = iost )
    if( iost .ne. 0 ) then
        write(*,*)'crtm_parameter, please check it.'
        stop 'read_namelist crtm_parameter'
    endif

    ! Load wrf parameters
    iost = 0
    read ( unit = namelist_unit, nml = wrf_parameters, iostat = iost )
    if( iost .ne. 0 ) then
        write(*,*)'wrf_parameter, please check it.'
        stop 'read_namelist wrf_parameter'
    endif

    close ( unit = namelist_unit )


    ! Convert sensor lon to radians
    sensor_lons(1:channel_count) &
    = sensor_lons(1:channel_count) /180*3.14159

  end subroutine  read_namelist



end module namelist_mod
! =============================================================================


! =============================================================================
! Module for constants
! =============================================================================

module constants
  ! Constants
  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re      = 6378000.0
  end module constants






! =============================================================================
! Module for parallelization
! =============================================================================
module mpi_module

   use MPI

   integer, parameter :: IO_NODE = 0
   integer :: nprocs, my_proc_id, comm, ierr, request, nprocs_mod, tag
   double precision   :: time_start, time_end
   integer, dimension(MPI_STATUS_SIZE)   :: status

   real                                  :: mpisend
   real, allocatable, dimension(:  )     :: mpirecv

   real, allocatable, dimension(:  )     :: mpisend1d
   real, allocatable, dimension(:,:)     :: mpirecv1d

   integer :: prev, next

   contains

  !-----------------------------------------------------------------------------
  !  PURPOSE: to basically set up a communicator for a rectangular mesh.
  !------------------------------------------------------------------------------
   subroutine parallel_start()

      implicit none

      ! Local variables
      integer :: mpi_rank, mpi_size
      integer :: ierr

      ! Find out our rank and the total number of processors
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)
      time_start = MPI_Wtime()

      comm = MPI_COMM_WORLD
      nprocs = mpi_size
      my_proc_id = mpi_rank
   end subroutine parallel_start


  !----------------------------------------------------------------------------
  !  PURPOSE: Free up, deallocate, and for MPI, finalize.
  !----------------------------------------------------------------------------
   subroutine parallel_finish()

      implicit none

      ! Local variables
      integer :: ierr

      time_end = MPI_Wtime()
      call MPI_Finalize(ierr)

   end subroutine parallel_finish

end module mpi_module
! =============================================================================


! =============================================================================
! Module of WRF variables and subroutines needed for CRTM
! =============================================================================
module wrf_utils

  use constants


  implicit none


  ! Allocatable arrays for root process
  REAL, allocatable, dimension(:,:) :: xlat, xlong, lat, lon, psfc, hgt, tsk, landmask
  REAL, allocatable, dimension(:,:,:) :: p, pb, pres, ph, phb, t, tk
  REAL, allocatable, dimension(:,:,:) :: qvapor, qcloud, qrain, qice, qsnow, qgraup


  ! WRF domain scalars
  INTEGER        :: xmax, ymax, zmax

  ! Allocatable arrays for all processes. Used to hold split up wrf fields
  REAL, allocatable, dimension(:) :: delz
  REAL, allocatable, dimension(:) :: xlat_sub, xlong_sub, lat_sub, lon_sub
  REAL, allocatable, dimension(:) :: psfc_sub, hgt_sub, tsk_sub, landmask_sub
  REAL, allocatable, dimension(:,:) :: p_sub, pb_sub, pres_sub, ph_sub, phb_sub
  REAL, allocatable, dimension(:,:) :: t_sub, tk_sub
  REAL, allocatable, dimension(:,:) :: qvapor_sub, qcloud_sub, qrain_sub
  REAL, allocatable, dimension(:,:) :: qice_sub, qsnow_sub, qgraup_sub
  INTEGER, allocatable, dimension(:) :: xpos, ypos
  INTEGER, allocatable, dimension(:,:) :: xmesh, ymesh

  ! Variables used for dividing up the domain
  INTEGER :: n_wrf_cols, total_wrf_cols, total_wrf_cols_regular


  CONTAINS


  ! --------------------------------------------------------------------------
  ! Subroutines to load WRF variables into root process
  ! --------------------------------------------------------------------------

  ! Allocate and read wrf variables into root process
  subroutine root_allocate_load_wrf_variables(wrf_fname)
    use mpi_module
    implicit none
    CHARACTER(256), intent(in) :: wrf_fname

    if ( my_proc_id == 0) then
      call root_allocate_wrf_variables(wrf_fname)
      call root_load_wrf_variables(wrf_fname)
    endif
  end subroutine root_allocate_load_wrf_variables

  ! Allocate all wrf variables in root process
  subroutine root_allocate_wrf_variables( wrf_fname )

    use netcdf
    implicit none

    CHARACTER(256), intent(in) :: wrf_fname

    ! Determine domain dimensions
    call get_ij( wrf_fname, xmax, ymax, zmax)

    ! Allocate arrays
    allocate(  xlat(xmax,ymax)  )  ! latitude
    allocate(  xlong(xmax,ymax) )  ! longitude
    allocate(  lat(xmax,ymax)   )  ! in radian
    allocate(  lon(xmax,ymax)   )  ! in radian
    allocate(  p(xmax,ymax,zmax)  )
    allocate(  pb(xmax,ymax,zmax)  )
    allocate(  pres(xmax,ymax,zmax)  )
    allocate(  ph(xmax,ymax,zmax+1)  )
    allocate(  phb(xmax,ymax,zmax+1)  )
    allocate(  t(xmax,ymax,zmax)  )
    allocate(  tk(xmax,ymax,zmax)  )
    allocate(  qvapor(xmax,ymax,zmax)  )
    allocate(  qcloud(xmax,ymax,zmax)  )
    allocate(  qrain(xmax,ymax,zmax)  )
    allocate(  qice(xmax,ymax,zmax)  )
    allocate(  qsnow(xmax,ymax,zmax)  )
    allocate(  qgraup(xmax,ymax,zmax)  )
    allocate(  psfc(xmax,ymax)  )
    allocate(  hgt(xmax,ymax)  )
    allocate(  tsk(xmax,ymax)  )
    allocate(  landmask(xmax,ymax)  )
  end subroutine root_allocate_wrf_variables


  ! Allocate all wrf variables into root
  subroutine root_load_wrf_variables( wrf_fname )
    use netcdf
    implicit none
    CHARACTER(256), intent(in) :: wrf_fname

    call get_variable2d(wrf_fname,'XLAT',xmax,ymax,1,xlat)
    call get_variable2d(wrf_fname,'XLONG',xmax,ymax,1,xlong)
    call get_variable3d(wrf_fname,'P',xmax,ymax,zmax,1,p)
    call get_variable3d(wrf_fname,'PB',xmax,ymax,zmax,1,pb)
    call get_variable3d(wrf_fname,'PH',xmax,ymax,zmax+1,1,ph)
    call get_variable3d(wrf_fname,'PHB',xmax,ymax,zmax+1,1,phb)
    call get_variable3d(wrf_fname,'T',xmax,ymax,zmax,1,t)
    call get_variable3d(wrf_fname,'QVAPOR',xmax,ymax,zmax,1,qvapor)
    call get_variable2d(wrf_fname,'PSFC',xmax,ymax,1,psfc)
    call get_variable2d(wrf_fname,'TSK',xmax,ymax,1,tsk)
    call get_variable2d(wrf_fname,'HGT',xmax,ymax,1,hgt)
    call get_variable3d(wrf_fname,'QCLOUD',xmax,ymax,zmax,1,qcloud)
    call get_variable3d(wrf_fname,'QRAIN',xmax,ymax,zmax,1,qrain)
    call get_variable3d(wrf_fname,'QICE',xmax,ymax,zmax,1,qice)
    call get_variable3d(wrf_fname,'QSNOW',xmax,ymax,zmax,1,qsnow)
    call get_variable3d(wrf_fname,'QGRAUP',xmax,ymax,zmax,1,qgraup)

    ! Processing WRF data
    lat = xlat/180.0*3.14159
    lon = xlong/180.0*3.14159
    pres = P + PB
    tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )

    ! Removing spurious values
    where(qvapor.lt.0.0) qvapor=1.0e-8
    where(qcloud.lt.0.0) qcloud=0.0
    where(qice.lt.0.0)   qice=0.0
    where(qrain.lt.0.0)  qrain=0.0
    where(qsnow.lt.0.0)  qsnow=0.0
    where(qgraup.lt.0.0) qgraup=0.0

    !qcloud = 0.5
    !qrain  = 0.5
    !qice   = 0.5
    !qsnow  = 0.5
    !qgraup = 0.5

    !qcloud = qcloud*0 !1000
    !qrain  = qrain* 0 !1000
    !qice   = qice * 0 !1000
    !qsnow  = qsnow* 0 !1000
    !qgraup = qgraup*0 !1000
    write(*,*) 'Max value for qcloud', maxval( qcloud )
    write(*,*) 'Max value for  qrain', maxval(  qrain )
    write(*,*) 'Max value for   qice', maxval(   qice )
    write(*,*) 'Max value for  qsnow', maxval(  qsnow )
    write(*,*) 'Max value for qgraup', maxval( qgraup )

  end subroutine root_load_wrf_variables



  ! ---------------------------------------------------------------------------
  ! Subroutines to pass information from root process to all processes
  ! ---------------------------------------------------------------------------
  ! Subroutine to break up wrf domain
  subroutine wrf_domain_breakup()

    call init_domain_breakup()
    call mpi_scatter2d_domain()
    call mpi_scatter3d_domain()

  end subroutine wrf_domain_breakup


  ! Subroutine to figure out how much wrf fields should be divided among WRF
  ! processes. Also allocates memory for subdomains
  subroutine init_domain_breakup()

    use mpi_module
    implicit none
    ! Counter variables
    integer :: ii,jj

    ! Pass dimensional info from root process to all processes
    call MPI_BCAST( xmax, 1, MPI_INTEGER, 0, comm, ierr )
    call MPI_BCAST( ymax, 1, MPI_INTEGER, 0, comm, ierr )
    call MPI_BCAST( zmax, 1, MPI_INTEGER, 0, comm, ierr )
    total_wrf_cols = xmax*ymax

    ! Determine number of atmospheric columns going to each process
    n_wrf_cols = CEILING( (xmax * ymax * 1.)/nprocs )
    total_wrf_cols_regular = n_wrf_cols * nprocs

    ! Allocate memory for subdomains
    allocate(  xlat_sub       (n_wrf_cols)  )  ! latitude
    allocate(  xlong_sub      (n_wrf_cols) )  ! longitude
    allocate(  lat_sub        (n_wrf_cols)   )  ! in radian
    allocate(  lon_sub        (n_wrf_cols)   )  ! in radian
    allocate(  p_sub          (n_wrf_cols,zmax)  )
    allocate(  pb_sub         (n_wrf_cols,zmax)  )
    allocate(  pres_sub       (n_wrf_cols,zmax)  )
    allocate(  ph_sub         (n_wrf_cols,zmax+1)  )
    allocate(  phb_sub        (n_wrf_cols,zmax+1)  )
    allocate(  delz           (zmax)  )
    allocate(  t_sub          (n_wrf_cols,zmax)  )
    allocate(  tk_sub         (n_wrf_cols,zmax)  )
    allocate(  qvapor_sub     (n_wrf_cols,zmax)  )
    allocate(  qcloud_sub     (n_wrf_cols,zmax)  )
    allocate(  qrain_sub      (n_wrf_cols,zmax)  )
    allocate(  qice_sub       (n_wrf_cols,zmax)  )
    allocate(  qsnow_sub      (n_wrf_cols,zmax)  )
    allocate(  qgraup_sub     (n_wrf_cols,zmax)  )
    allocate(  psfc_sub       (n_wrf_cols)  )
    allocate(  hgt_sub        (n_wrf_cols)  )
    allocate(  tsk_sub        (n_wrf_cols)  )
    allocate(  landmask_sub   (n_wrf_cols)  )
    allocate(  xpos           (n_wrf_cols)  )
    allocate(  ypos           (n_wrf_cols)  )

    ! Setting up meshgrid of indices
    allocate( xmesh(xmax,ymax) ); allocate( ymesh(xmax,ymax) )
    fill_xmesh: do ii = 1, xmax
      xmesh(ii,:) = ii
    enddo fill_xmesh
    fill_ymesh: do ii = 1, ymax
      ymesh(:,ii) = ii
    enddo fill_ymesh

  end subroutine init_domain_breakup



  ! Subroutine to mpi_scatter wrf domain 2d fields into processes' subdomains
  subroutine mpi_scatter2d_domain()
    use mpi_module

    implicit none
    ! Define buffer variables
    real, dimension( total_wrf_cols_regular )   :: fbuffer
    integer, dimension( total_wrf_cols_regular ):: ibuffer

    ! Initialize buffers
    ibuffer = -99999
    fbuffer = -99999.

    ! Scatter xmesh and ymesh
    if (my_proc_id == 0) write(*,*) size( pack(xmesh,.true.) )


    ibuffer( 1:total_wrf_cols ) = pack(xmesh, .true.)
    call MPI_SCATTER( ibuffer, n_wrf_cols, MPI_INTEGER, xpos, n_wrf_cols, &
                      MPI_INTEGER, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    ibuffer( 1:total_wrf_cols ) = pack(ymesh, .true.)
    call MPI_SCATTER( ibuffer, n_wrf_cols, MPI_INTEGER, ypos, n_wrf_cols, &
                      MPI_INTEGER, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )


    ! Scatter 2d variables
    fbuffer( 1:total_wrf_cols ) = pack(xlat, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, xlat_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(xlong, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, xlong_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(lat, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, lat_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(lon, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, lon_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(psfc, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, psfc_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(hgt, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, hgt_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(tsk, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, tsk_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(landmask, .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, landmask_sub, &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

  end subroutine mpi_scatter2d_domain




  ! Subroutine to mpi_scatter wrf domain 3d fields into processes' subdomains
  subroutine mpi_scatter3d_domain()
    use mpi_module

    implicit none
    ! Define buffer variables
    real, dimension( total_wrf_cols_regular )   :: fbuffer
    ! Define counter variables
    integer :: kk

    ! Initialize buffer
    fbuffer = -99999.

    ! Scatter 3d variables with zmax dimensions
    scatter_zmax_levels: do kk = 1, zmax

      fbuffer( 1:total_wrf_cols ) = pack(p(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, p_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(pb(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, pb_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(pres(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, pres_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(ph(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, ph_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(phb(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, phb_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(t(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, t_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(tk(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, tk_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qvapor(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qvapor_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qcloud(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qcloud_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qrain(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qrain_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qice(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qice_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qsnow(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qsnow_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

      fbuffer( 1:total_wrf_cols ) = pack(qgraup(:,:,kk), .true.)
      call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, qgraup_sub(:,kk), &
                        n_wrf_cols, MPI_REAL, 0, comm, ierr  )
      call MPI_BARRIER( comm, ierr )

    enddo scatter_zmax_levels

    ! Dealing with variables zmax+1 levels
    fbuffer( 1:total_wrf_cols ) = pack(ph(:,:,zmax+1), .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, ph_sub(:,zmax+1), &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

    fbuffer( 1:total_wrf_cols ) = pack(phb(:,:,zmax+1), .true.)
    call MPI_SCATTER( fbuffer, n_wrf_cols, MPI_REAL, phb_sub(:,zmax+1), &
                      n_wrf_cols, MPI_REAL, 0, comm, ierr  )
    call MPI_BARRIER( comm, ierr )

  end subroutine mpi_scatter3d_domain


end module wrf_utils
! =============================================================================




! =============================================================================
! CRTM module for variable definitions and running radiative transfer
! =============================================================================
module crtm_utils

  use CRTM_Module
  use constants

  implicit none

  ! CRTM profile dimensions
  INTEGER :: N_LAYERS, N_CLOUDS
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1

  ! Other CRTM variables
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: CRTM_COEFF_PATH
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels

  real, allocatable, dimension(:,:,:) :: BTsend, BT

  ! CRTM structures
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)


  CONTAINS



  ! Subroutine to initialize atmospheric structure for crtm
  subroutine crtm_init_atm()

    use wrf_utils

    CALL CRTM_Atmosphere_Create( Atm, zmax, N_ABSORBERS, zmax*5, N_AEROSOLS)
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
      Message = 'Error allocating CRTM Atmosphere structures'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF

  end subroutine crtm_init_atm




  ! Subroutine to initialize crtm for the sensor
  subroutine crtm_init_sensor( Sensor_Id )

    use namelist_mod

    implicit none
    CHARACTER(len=50), intent(in) :: Sensor_Id

    Error_Status = CRTM_Init(   &
                        (/trim(Sensor_Id)/), &  ! Input... must be an array, hencethe (/../)
                        ChannelInfo  , &  ! Output
                        Load_CloudCoeff = .True., &
                        IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                        IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                        File_Path= './coefficients/', & !CRTM_COEFF_PATH, &
                        Quiet=.true. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF

  end subroutine crtm_init_sensor




  ! Subroutine to initialize crtm for one specific sensor channel
  subroutine crtm_init_single_channel( ch_id )

    use wrf_utils
    implicit none
    integer, intent(in) :: ch_id

    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), &
                                            Channel_Subset = (/ch_id/) )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error initializing CRTM sensor channels'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF

    n_Channels = 1

    ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
      Message = 'Error allocating structure arrays'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF


  end subroutine crtm_init_single_channel



  ! Subroutine to prepare crtm structured types for one wrf column
  subroutine crtm_input_wrf_col( sat_lon, sat_h, ll )

    use mpi_module
    use wrf_utils

    implicit none
    real, intent(in) :: sat_lon, sat_h
    integer, intent(in) :: ll  ! Index position along the flattened sub domain

    INTEGER :: l, m, irec, pt_end, pt_start, n_pts
    integer :: ncid,ncrcode
    integer :: tt,v,z,n,reci,ens,n_ec
    integer :: pt_id, tot_pts
    INTEGER :: ncl,icl,k1,k2


    ! Satellite viewing geometry
    sat_dis = Re**2.0 + (Re+sat_h)**2.0   &
              - 2.0*Re*(Re+sat_h)*cos(lon_sub(ll)-sat_lon)*cos(lat_sub(ll))
    sat_dis = sqrt( sat_dis )
    SCAN_ANGLE = sqrt(1-(cos(lon_sub(ll)-sat_lon)*cos(lat_sub(ll)))**2)
    SCAN_ANGLE = Re/sat_dis * SCAN_ANGLE
    SCAN_ANGLE = (180.0/3.14159) * asin( SCAN_ANGLE )
    ZENITH_ANGLE = cos(lon_sub(ll)-sat_lon) * cos(lat_sub(ll))
    ZENITH_ANGLE = SCAN_ANGLE + 180.0/3.14159*acos(ZENITH_ANGLE)
    if ( abs( ZENITH_ANGLE ) >= 80 ) then
      write(*,*) my_proc_id, lon_sub(ll), lat_sub(ll)
      write(*,*) sat_dis, SCAN_ANGLE, ZENITH_ANGLE
      stop
    endif

    ! calculating delz
    do z=1,zmax
      if(z.eq.1) then
        delz(z) = (ph_sub(ll,z+1) + phb_sub(ll,z+1)) / 9.806 - hgt_sub(ll)
      else
        delz(z) = ( (ph_sub(ll,z+1) + phb_sub(ll,z+1) ) &
                     -(ph_sub(ll,z) + phb_sub(ll,z)) &
                  )/9.806
      endif
    enddo
    if (delz(1) <= 0.) delz(1) = delz(2)

    !---Atmospheric Profile
    atm(1)%Climatology         = TROPICAL
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)
    atm(1)%Level_Pressure(0) = (pres_sub(ll,zmax)*3.0/2.0 &
                               - pres_sub(ll,zmax-1)/2.0)/100.0  ! convert from Pa to hPA

    do z=zmax,1,-1
      if(z.eq.1) then
        atm(1)%Level_Pressure(zmax-z+1) = psfc_sub(ll)/100.0    ! Pa -> hPa
        !max(psfc_sub(ll), pres_sub(ll,1)*3.0/2.0-pres_sub(ll,2)/2.0)/100.0
      else
        atm(1)%Level_Pressure(zmax-z+1) = ((pres_sub(ll,z-1)+pres_sub(ll,z))/2.0)/100.0  ! convert from Pa to hPA
      endif
      atm(1)%Pressure(zmax-z+1)       = pres_sub(ll,z) / 100.0
      atm(1)%Temperature(zmax-z+1)    = tk_sub(ll,z)
      atm(1)%Absorber(zmax-z+1,1)     = qvapor_sub(ll,z)*1000.0
    enddo
    atm(1)%Absorber(:,2) = 5.0E-02
    !---Cloud Profile
    do z=1,zmax*5
     atm(1)%Cloud(z)%Type = 0
     atm(1)%Cloud(z)%Effective_Radius = 0.0
     atm(1)%Cloud(z)%Water_Content = 0.0
    enddo
    ncl = 0
    icl = 0
    !--calculating # of clouds (cloud and rain)
    do z=zmax,1,-1
      if(qcloud_sub(ll,z).gt.1e-6) then
        ncl = ncl + 1
      endif
      if(qrain_sub(ll,z).gt.1e-6) then
        ncl = ncl + 1
      endif
      if(qice_sub(ll,z).gt.1e-6) then
        ncl = ncl + 1
      endif
      if(qsnow_sub(ll,z).gt. 1e-6) then
        ncl = ncl + 1
      endif
      if(qgraup_sub(ll,z).gt. 1e-6) then
        ncl = ncl + 1
      endif
    enddo


    !--Data for cloud
    atm(1)%n_Clouds         = ncl
    IF ( atm(1)%n_Clouds > 0 ) THEN
    do z=zmax,1,-1
      if(qcloud_sub(ll,z).gt. 1e-6) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = WATER_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qcloud_sub(ll,z)*pres_sub(ll,z)/287.2/(tk_sub(ll,z)+0.61*(qvapor_sub(ll,z)/(1+qvapor_sub(ll,z))))*delz(z)
      endif
    enddo

    do z=zmax,1,-1
      if(qrain_sub(ll,z).gt. 1e-6) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = RAIN_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qrain_sub(ll,z)*pres_sub(ll,z)/287.2/(tk_sub(ll,z)+0.61*(qvapor_sub(ll,z)/(1+qvapor_sub(ll,z))))*delz(z)
      endif
    enddo

    do z=zmax,1,-1
      if(qice_sub(ll,z).gt. 1e-6) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = ICE_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qice_sub(ll,z)*pres_sub(ll,z)/287.2/(tk_sub(ll,z)+0.61*(qvapor_sub(ll,z)/(1+qvapor_sub(ll,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qsnow_sub(ll,z).gt. 1e-6) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = SNOW_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qsnow_sub(ll,z)*pres_sub(ll,z)/287.2/(tk_sub(ll,z)+0.61*(qvapor_sub(ll,z)/(1+qvapor_sub(ll,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qgraup_sub(ll,z).gt. 1e-6) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qgraup_sub(ll,z)*pres_sub(ll,z)/287.2/(tk_sub(ll,z)+0.61*(qvapor_sub(ll,z)/(1+qvapor_sub(ll,z))))*delz(z)
      endif
    enddo
    ENDIF


    !---Surface data
    if( landmask_sub(ll) .eq.1.0) then
     sfc(1)%Water_Coverage = 0.0_fp
     sfc(1)%Land_Coverage = 1.0_fp
     sfc(1)%Land_Temperature = tsk_sub(ll)
     sfc(1)%Soil_Temperature = tsk_sub(ll)
    else
     sfc(1)%Water_Coverage = 1.0_fp
     sfc(1)%Land_Coverage = 0.0_fp
     sfc(1)%Water_Type = 1  ! Sea water
     sfc(1)%Water_Temperature = tsk_sub(ll)
    endif


    ! Input satellite viewing geometry
    ! The Sensor_Scan_Angle is optional.
    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE )

    write(*,*) ncl

  end subroutine crtm_input_wrf_col



  ! Subroutine to run CRTM for single column, single sensor, single channel
  subroutine crtm_single_solve( bt )

    implicit none
    real, intent(out) :: bt

    ! Preemptive print stuff
     write(*,*) 'Pre-CRTM cloud count', atm(1)%n_Clouds


    ! Run CRTM forward model
    Options%RT_Algorithm_ID = RT_SOI
    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution , &
                                 Options = Options )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF

    ! Output solution
    bt = real( RTSolution(1,1)%Brightness_Temperature )

    ! Print out stuff
    write(*,*) 'Cloudy BT', bt, 'Post-CRTM cloud count', atm(1)%n_Clouds
    if ( atm(1)%n_Clouds > 5) STOP

  end subroutine crtm_single_solve


  ! Subroutine to kill crtm for the current channel
  subroutine kill_crtm()
    implicit none

    Error_Status = CRTM_Destroy( ChannelInfo )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error destroying CRTM'
      CALL Display_Message( 'crtm', Message, FAILURE )
      STOP
    END IF

    deallocate( RTSolution )
  end subroutine kill_crtm


end module crtm_utils
! =============================================================================



! ============================================================================
! Module to post-process the outputs of CRTM
! ============================================================================
module post_process

  use mpi_module
  use wrf_utils
  use namelist_mod

  ! 3D BT field
  real, allocatable, dimension(:,:,:) :: BT_field3d


  contains

  ! Subroutine to gather all BT back together to root process
  subroutine gather_BT_to_root( BT_sub )

    implicit none

    ! Input BT subdomain values
    real, dimension( channel_count, n_wrf_cols ), intent(in) :: BT_sub

    ! Gathered x, y positions and BT
    integer, allocatable, dimension(:) :: xpos_all, ypos_all
    real, allocatable, dimension(:,:) :: BT_all
    logical, allocatable, dimension(:,:,:) :: BT_mask
    real :: t0,t1

    ! Counter variables
    integer :: cc, ii

    ! Allocate BT_field3d at root
    if ( my_proc_id == 0 ) then
      allocate( BT_field3d( channel_count, xmax, ymax ) )
      allocate( xpos_all ( total_wrf_cols_regular ) )
      allocate( ypos_all ( total_wrf_cols_regular ) )
      allocate( BT_all   ( channel_count, total_wrf_cols_regular ) )
      allocate( BT_mask  ( channel_count, xmax, ymax ) ) 
      BT_mask(:,:,:) = .true.
    endif

    ! Gather indices into root process
    call MPI_GATHER( xpos, n_wrf_cols, MPI_INTEGER, xpos_all, n_wrf_cols,  &
                     MPI_INTEGER, 0, comm, ierr )
    call MPI_GATHER( ypos, n_wrf_cols, MPI_INTEGER, ypos_all, n_wrf_cols,  &
                     MPI_INTEGER, 0, comm, ierr )

    ! Gather BT field into root process
    if (my_proc_id == 0) call cpu_time( t0 )
    gather_bt_ch_loop: do cc= 1, channel_count
      call MPI_GATHER( BT_sub(cc,:), n_wrf_cols, MPI_REAL, BT_all(cc,:), n_wrf_cols, &
                       MPI_REAL, 0, comm, ierr )
    enddo gather_bt_ch_loop
    if (my_proc_id == 0) call cpu_time( t1 )
    if (my_proc_id == 0) write(*,'(a,x,f5.1,x,a)') 'BT MPI_Gather took', t1 -t0, 'seconds'

    ! Using unpack to construct the 3d bt field
    if (my_proc_id == 0) then
      call cpu_time(t0)
      BT_field3d = -999.
      unpack_channel_loop: do cc = 1, channel_count
        BT_field3d(cc,:,:) &
        = unpack( BT_all(cc,1:total_wrf_cols), BT_mask(cc,:,:), & !BT_mask(cc,:), &
                  BT_field3d(cc,:,:) )
      enddo unpack_channel_loop
    endif
    if (my_proc_id == 0) call cpu_time(t1)
    if (my_proc_id == 0) write(*,'(a,x,f5.1,x,a)') 'Unpack BT reshaping', t1 -t0, 'seconds'



    ! Barrier to make sure all processes wait for root proc to clear
    call MPI_BARRIER( comm, ierr )

  end subroutine gather_BT_to_root


  ! Subroutine to write BT fields into netcdf file
  ! Will only perform writing using the root process
  ! Note that BT_field3d to be written exists within this module!
  ! Based off https://www.unidata.ucar.edu/software/netcdf/examples/programs/pres_temp_4D_wr.f90
  subroutine root_write_bt_ncfile( out_fname )

    ! Important modules
    use namelist_mod   ! Supplies channel and sensor info
    use netcdf         ! Supplies routines to write ncfile
    use wrf_utils      ! Supplies domain info and variables

    implicit none
    character(256), intent(in) :: out_fname
    integer :: ncid                         ! ncfile id
    character(256) :: sensor_channel_string ! ncfile varname for the channel
    integer :: lat_dimid, lon_dimid, bt_dimid, lat_varid, lon_varid
    integer, dimension(99) :: ch_varid
    integer, dimension(2):: dimid2d
    integer :: cc     ! Channel counter variable
    integer :: rcode  ! Error code handler
    integer :: sensor_name_len

    ! Create ncfile
    rcode = nf_create( trim(out_fname), nf_clobber, ncid )
    write(*,*) 'created ncfile'

    ! ---- DEFINE DIMENSIONS AND VARIABLE ------ !
    ! Define the geo dimensions
    rcode = nf_def_dim(ncid, 'latitude' , ymax, lat_dimid)
    rcode = nf_def_dim(ncid, 'longitude', xmax, lon_dimid)
    write(*,*) 'Defined lat lon dimensions'

    ! Construct 2d array dimension ids
    dimid2d = (/ lon_dimid, lat_dimid /)
    write(*,*) 'Defined 2d Ids ', dimid2d

    ! Define lat lon variables and add attributes
    rcode = nf_def_var( ncid, 'latitude' , NF_REAL, 2, dimid2d, lat_varid )
    rcode = nf_def_var( ncid, 'longitude', NF_REAL, 2, dimid2d, lon_varid )
    write(*,*) 'Defined lat lon variables'

    rcode = nf_put_att_text( ncid, lat_varid, 'units', &
                             len('degrees north'), 'degrees north' )
    rcode = nf_put_att_text( ncid, lon_varid, 'units', &
                             len('degrees east' ), 'degrees east' )
    write(*,*) 'Defined lat lon attributes'


    ! For each channel, define variable and add attributes
    ncf_ch_define: do cc=1, channel_count
      write( sensor_channel_string , '(a,a,I0.3)' )  &
                    trim(sensor_list(cc)), '_ch', channel_id(cc)
      write(*,*) trim(sensor_channel_string)
      rcode = nf_def_var( ncid, trim(sensor_channel_string),  &
                          NF_REAL, 2, dimid2d, ch_varid(cc) )
      rcode = nf_put_att_text( ncid, ch_varid(cc), 'units', &
                               len('Kelvins'), 'Kelvins' )
      rcode = nf_put_att_text( ncid, ch_varid(cc), 'masked_value', &
                               len('-999.'), '-999.' )
    enddo ncf_ch_define
    write(*,*) 'Defined BT variables'


    ! End definition mode
    rcode = nf_enddef( ncid )
    write(*,*) 'Closed definition mode'

    ! ----- INSERT DATA INTO THE NETCDF FILE ------ !
    ! Insert the coordinate meshgrid
    rcode = nf_put_var( ncid, lat_varid, xlat )
    rcode = nf_put_var( ncid, lon_varid, xlong )
    write(*,*) 'Inserted lat lon meshgrid into file'

    ! Insert the computed BT fields
    ncf_ch_write: do cc=1, channel_count
      rcode = nf_put_var( ncid, ch_varid(cc), BT_field3d(cc,:,:) )
    enddo ncf_ch_write
    write(*,*) 'Inserted BT field'


    ! ----- SAVE AND CLOSE NETCDF FILE ------ !
    rcode = nf_close(ncid)
    write(*,*) 'Closed and flushed file into hard drive'


  end subroutine root_write_bt_ncfile


end module post_process
! ============================================================================
