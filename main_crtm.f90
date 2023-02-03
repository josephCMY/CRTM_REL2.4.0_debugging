! =========================================================================
!              PARALLELIZED OFFLINE GEOSTATIONARY IR CRTM CODE
! =========================================================================
! Written by M.-Y. Chan
! Based on xb_to_radiance subroutine in EnSRF/xb.f
! Generic. Works for any geostationary IR.
! =========================================================================
! IMPORTANT NOTES:
! 1) Always check if the desired satellite and channels are set up properly.
! 2) To run after compilation (crtm.exe):
!    >>> PARALLEL_RUN crtm.exe wrf_file crtm_output_name.bin
!    where PARALLEL_RUN can be mpirun, srun, ibrun (whichever suitable)
! =========================================================================



PROGRAM crtm

  ! -----------------------------------------------------------------------
  ! 1. Modules used in this program
  ! -----------------------------------------------------------------------
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use wrf_utils
  use crtm_utils
  use namelist_mod
  use post_process
  use constants


  ! -----------------------------------------------------------------------
  ! 2. Predefining variables (no implicit variables)
  ! -----------------------------------------------------------------------
  implicit none

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'

  ! Input and output file names
  CHARACTER(256) :: wrf_fname
  CHARACTER(256) :: out_fname

  ! Counter variables
  integer :: cc, ll

  ! BT field (dims: channel, column)
  real, allocatable, dimension(:,:) :: BT_sub




  ! -----------------------------------------------------------------------
  ! 3. Read in namelist
  ! -----------------------------------------------------------------------
  call read_namelist()


  ! -----------------------------------------------------------------------
  ! 4. Read in wrf file name and output file name
  ! -----------------------------------------------------------------------
  ! Read in file names
  call getarg( 1, wrf_fname )
  call getarg( 2, out_fname )


  ! -----------------------------------------------------------------------
  ! 5. Initialize parallelization and load variables into processes
  ! -----------------------------------------------------------------------
  ! Initialize parallization
  call parallel_start()

  ! Load variables into root process
  call root_allocate_load_wrf_variables(wrf_fname)

  ! Scatter wrf fields among processes
  call wrf_domain_breakup()

  ! Checking out root process stuff
  if ( my_proc_id == 0) write(*,*) 'length of zone',n_wrf_cols

  ! Initializing number of layers and clouds
  N_LAYERS = zmax; N_CLOUDS = zmax*5

  ! Initialize BT field to hold stuff
  allocate( BT_sub (channel_count, total_wrf_cols_regular) )
  BT_sub = -999.


  ! Prepare CRTM atmospheric structure
  call crtm_init_atm()

  ! -----------------------------------------------------------------------
  ! 6. Iterate over each channel specified in namelist
  ! -----------------------------------------------------------------------
  channel_loop: do cc = 1, channel_count

    ! Initialize CRTM for sensor in question
    if( my_proc_id == 0) write(*,*) trim(sensor_list(cc) ), channel_id(cc)
    call crtm_init_sensor( sensor_list(cc) )

    ! Prepare CRTM for one channel
    call crtm_init_single_channel( channel_id(cc) )

    ! --------------------------------------------------------------------
    ! 7. Iterate over all columns held by the process
    ! --------------------------------------------------------------------
    wrf_column_loop: do ll = 1, n_wrf_cols


      ! Check if selected column falls within the desired location
      if( ( xpos(ll) < first_index ) .or. ( ypos(ll) < first_index ) ) &
        cycle wrf_column_loop
      if ( MODULO( ( xpos(ll) - first_index ), grid_thinning ) .ne. 0 ) &
        cycle wrf_column_loop
      if ( MODULO( ( ypos(ll) - first_index ), grid_thinning ) .ne. 0 ) &
        cycle wrf_column_loop

      ! Check if column falls outside of sensor's zone
      if ( ( xlat_sub(ll) < sensor_lat_min(cc) ) .or. &
           ( xlat_sub(ll) > sensor_lat_max(cc) )         ) &
          cycle wrf_column_loop
      if ( ( xlong_sub(ll) < sensor_lon_min(cc) ) .or. &
           ( xlong_sub(ll) > sensor_lon_max(cc) )         ) &
          cycle wrf_column_loop


      !! Debug mode. Skipping over all points but one.
      !if ( ( xpos(ll) .ne. 101 ) .or. (ypos(ll) .ne. 101) ) then
      !  cycle wrf_column_loop
      !endif


      ! Load wrf column into the crtm structure
      call crtm_input_wrf_col( sensor_lons(cc), sensor_heights(cc), ll )


      ! Run CRTM forward calculation
      call crtm_single_solve( BT_sub(cc,ll) )

    enddo wrf_column_loop

    ! Kill active CRTM before moving to the next channel
    ! Will restart CRTM for the next channel
    call kill_crtm()

  enddo channel_loop

  ! -----------------------------------------------------------------------
  ! 8. Gather BT to root process and generate 3d array of BT values
  ! -----------------------------------------------------------------------
  call MPI_Barrier(comm, ierr)
  call gather_BT_to_root( BT_sub )

  ! -----------------------------------------------------------------------
  ! 9. Output BT values as a netcdf file
  ! -----------------------------------------------------------------------
  ! Only activate this outputting in root
  if (my_proc_id == 0) call root_write_bt_ncfile( out_fname )

  ! Barrier for safety reasons
  call MPI_Barrier( comm, ierr )

  ! -----------------------------------------------------------------------
  ! 10. Kill MPI processes
  ! -----------------------------------------------------------------------
  call parallel_finish()



END PROGRAM crtm
