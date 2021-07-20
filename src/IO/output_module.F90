module output_module
  !/**********************************************************************\
  !| OUTPUT MODULE                                                        |
  !| -------------                                                        |
  !| This module contains the interface to write to an output file from a |
  !| state_type derived type. The output module supports spectral and     |
  !| physical output files in 3D in NetCDF format and vertical profiles   |
  !| in NetCDF format. It also supports a diagnostic file in CSV format.  |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type
  implicit none
  private
  public nc_ids_type,controllers_type
  public open_files,close_files,write_output_files
  public check

  !/**********************************************************************\
  !| FLAGS                                                                |
  !\**********************************************************************/
  integer, parameter :: io_spectral = 1
  integer, parameter :: io_vertical = 2
  integer, parameter :: io_single_record = 4

  !/**********************************************************************\
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type nc_ids_type
    !/********************************************************************\
    !| NC_IDS_TYPE                                                        |
    !| -----------                                                        |
    !| This derived type contains the NetCDF IDs for the variables output |
    !| in any NetCDF file. These are not intended to be manipulated by    |
    !| the user but rather controlled by open_file and used by            |
    !| write_file.                                                        |
    !\********************************************************************/
    character(len=100) :: file_name ! The name of the file
    integer :: ncid ! The netcdf ID of the file
    integer :: time_dimid !| The netcdf ID of the time dimension
    integer :: x_dimid,y_dimid,z_dimid ! The netcdf IDs of the dimensions
    integer :: x,y,z ! The netcdf IDs of the x, y, and z variables
    integer :: kx,ky,kz ! The netcdf IDs of the wavenumber variables
    integer :: time,timestep,dt ! The netcdf IDs of the time variables
    integer :: t,c,u,v,w ! The netcdf IDs of the basic arrays
    integer :: taux,caux,uaux,vaux,waux ! The netcdf IDs of the aux arrays
    logical :: stationary_output ! Whether to include stationary outputs
    integer :: x_remaps,y_remaps ! The number of remaps in x and y
    integer :: flags ! The io flags of the output file
    integer :: step ! The current step of the output
    integer :: dims ! The number of dimensions in the output
    integer :: dimids(4) ! The array of netcdf dimension IDs
  end type nc_ids_type

  type controllers_type
    !/********************************************************************\
    !| CONTROLLERS_TYPE                                                   |
    !| ----------------                                                   |
    !| This derived type contains the parameters governing the IO         |
    !| routines. It contains the handles to the relevant data files and   |
    !| the information on the file names and steps in between outputs.    |
    !| This derived type is typically read by open_files in a namelist    |
    !| called io, but the object can be edited manually prior to calling  |
    !| open_files instead.                                                |
    !\********************************************************************/
    type(nc_ids_type) :: restart ! The handle on the restart file
    type(nc_ids_type) :: data ! The handle on the data file
    type(nc_ids_type) :: profile ! The handle on the profile file

    integer :: steps_diag = 1 ! The number of steps between diag files
    integer :: steps_data = 100 ! The number of steps between data files
    integer :: steps_restart = 100 ! The number of steps between restart files
    integer :: steps_profile = 100 ! The number of steps between profile files

    logical :: stationary_output = .false. ! Whether to include stationary outputs
    
    character(len=100) :: file_diag = "diag.csv" ! The name of the diag file
    character(len=100) :: file_data = "data.nc" ! The name of the data file
    character(len=100) :: file_profile = "profiles.nc" ! The name of the profile file
    character(len=100) :: file_restart = "restart.nc" ! The name of the restart file

    character(len=100) :: file_input = "" ! The name of the input file

    integer :: unit_diag = 79 ! The unit of the diag file
  end type controllers_type

contains

  !/**********************************************************************\
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function mean_flux(grid,x,u,x_dx_mean,x_dz_mean,vel_mean)
    !/********************************************************************\
    !| MEAN_FLUX                                                          |
    !| --------                                                           |
    !| This function returns the mean turbulent flux of a quantity x,     |
    !| given the velocity u and the mean background gradients of x in the |
    !| x and z directions. If the velocity has a nonzero mean, that can   |
    !| be optionally specified to be subtracted.                          |
    !| It is assumed that the total field associated with x is given by   |
    !|     x_total = x + x_dx_mean * x + x_dz_mean * z                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   x: the physical array for which the turbulent mean is desired    |
    !|   u: the physical array for the component of velocity              |
    !|   x_dx_mean: the mean x-gradient for the array x                   |
    !|   x_dz_mean: the mean z-gradient for the array x                   |
    !|   vel_mean: the mean velocity to be subtracted from u              |
    !|   returns: the average turbulent flux of x in the domain           |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type) :: grid
    real(kind=dp) :: x(:,:,:)
    real(kind=dp) :: u(:,:,:)
    real(kind=dp) :: x_dx_mean,x_dz_mean
    real(kind=dp), optional :: vel_mean
    real(kind=dp) :: mean_flux

    real(kind=dp) :: rtmp,vel_meanin,vel
    integer :: i,j,k,ierror

    vel_meanin = 0.0_dp
    if (present(vel_mean)) vel_meanin = vel_mean

    mean_flux = 0.0_dp
    do k = 1,grid%z_p_count
      do j = 1,grid%y_p_count
        do i = 1,grid%x_n
          vel = u(i,j,k) - vel_meanin
          rtmp = x(i,j,k) + x_dx_mean * grid%lx(i) + x_dz_mean * grid%lz(k)
          mean_flux = mean_flux + grid%dv * vel * rtmp
        end do
      end do
    end do

    mean_flux = mean_flux / grid%volume 

    mean_flux = collect_sum(mean_flux)
  end function

  function dissipation(grid,time,x)
    !/********************************************************************\
    !| DISSIPATION                                                        |
    !| -----------                                                        |
    !| This function returns the mean dissipation a quantity x, that is,  |
    !| it returns (gradient(x))^2, where the square indicates the inner   |
    !| product of a vector with itself. This is useful in understanding   |
    !| the turbulence in the system and can serve as a proxy for the flux |
    !| in some circumstances.                                             |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   time: the current time of the simulation                         |
    !|   x: the physical array for which the turbulent mean is desired    |
    !|   returns: the average dissipation of x in the domain              |
    !\********************************************************************/
    use grid_module, only: allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use fields_module, only: sum_squares
    implicit none

    type(grid_type) :: grid
    real(kind=dp) :: time
    complex(kind=dp) :: x(:,:,:)

    complex(kind=dp), pointer :: work_s(:,:,:)
    real(kind=dp) :: dissipation

    call allocate_s(grid,work_s)
    
    call deriv_x(grid,x,work_s)
    dissipation = sum_squares(grid,work_s)

    call deriv_y(grid,x,work_s)
    dissipation = dissipation + sum_squares(grid,work_s)

    call deriv_z(grid,time,x,work_s)
    dissipation = dissipation + sum_squares(grid,work_s)

    dissipation = dissipation / grid%xyz_n

    deallocate(work_s)
  end function

  !/**********************************************************************\
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine open_files(grid,params,controls)
    !/********************************************************************\
    !| OPEN_FILES                                                         |
    !| ----------                                                         |
    !| This subroutine reads in the io parameters and opens all output    |
    !| files.                                                             |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   controls: a controllers_type object that handles io parameters   |
    !\********************************************************************/
    use comm_module, only: get_id,abort_comm
    implicit none
    
    type(grid_type) :: grid
    type(params_type) :: params
    type(controllers_type) :: controls

    integer :: ierror

    namelist /io/ controls

    read(params%unit,nml=io,iostat=ierror)  

    if (ierror .ne. 0) then 
      print*,"WARNING: While reading IO namelist, encountered error (",ierror,")"
      print*,"WARNING: Using defaults for IO."
    end if
    close(params%unit)
    
    ! Open the diagnostic file and write the header
    if (get_id() .eq. 0) then 
      open(unit=controls%unit_diag,FILE=controls%file_diag,action="write")
      write(controls%unit_diag,*) "# timestep,time,dt,rms_u,rms_v,rms_w,prod_uw,prod_vw,", &
        & "rms_t,rms_c,rms_cmt,flux_t,flux_c,diss_t,diss_c"
    end if

    ! Open the restart file
    call open_file(grid,params,controls%restart,controls%file_restart,ior(io_single_record,io_spectral))

    ! If the simulation is shearing, check if the user wants stationary outputs
    controls%data%stationary_output = controls%stationary_output
    if (grid%x_shear .eq. 0.0_dp .and. grid%y_shear .eq. 0.0_dp) then
      controls%data%stationary_output = .false.
    end if

    ! Open the data file
    call open_file(grid,params,controls%data,controls%file_data)
    
    ! Open the profile file
    call open_file(grid,params,controls%profile,controls%file_profile,io_vertical)
  end subroutine

  subroutine close_files(controls)
    !/********************************************************************\
    !| CLOSE_FILES                                                        |
    !| -----------                                                        |
    !| This subroutine closes all output files.                           |
    !|   controls: a controllers_type object that handles io parameters   |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(controllers_type) :: controls
    integer :: i,ierror
  
    if (get_id() .eq. 0) then
        CLOSE(UNIT=controls%unit_diag)
    end if
  
    call close_file(controls%restart)
    call close_file(controls%data)
    call close_file(controls%profile)
  end subroutine

  subroutine write_output_files(grid,params,controls,state)
    !/********************************************************************\
    !| WRITE_FILES                                                        |
    !| ----------                                                         |
    !| This subroutine evaluates whether a step should be output to the   |
    !| output files and then calls the relevant subroutine to write that  |
    !| step to file.                                                      |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   controls: a controllers_type object that handles io parameters   |
    !|   state: a state_type object that contains the simulation data     |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(controllers_type), intent(in) :: controls
    type(state_type), intent(in) :: state
  
    ! Write to the diafnostics file
    if (MOD(state%timestep,controls%steps_diag) .eq. 0) then
      call write_diagnostics(grid,params,controls%unit_diag,state)
    end if

    ! Write to the profile file
    if (MOD(state%timestep,controls%steps_profile) .eq. 0) then
      call write_file(grid,params,controls%profile,state)
    end if

    ! Write to the restart file
    if (MOD(state%timestep,controls%steps_restart) .eq. 0) then
      call write_file(grid,params,controls%restart,state)
    end if

    ! Write to the data file
    if (MOD(state%timestep,controls%steps_data) .eq. 0) then
      call write_file(grid,params,controls%data,state)
    end if
  end subroutine

  subroutine write_diagnostics(grid,params,unit,state)
    !/********************************************************************\
    !| WRITE_DIAGNOSTICS                                                  |
    !| -----------------                                                  |
    !| This subroutine computes various diagnostics regarding the state   |
    !| of the simulation and then writes those to a CSV file.             |
    !|   grid: a grid_type object that describes the grid                 |
    !|   params: a params_type object that specifies the parameters       | 
    !|   unit: the unit of the output file to write to                    |
    !|   state: a state_type object that contains the simulation data     |
    !\********************************************************************/
    use comm_module, only: get_id
    use grid_ops_module, only: div
    use fields_module, only: sum_squares,sum_squares_vector,sum_product_p
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(state_type) :: state
    integer :: unit

    real(kind=dp) :: flux_t,flux_c,hflux_t,hflux_c
    real(kind=dp) :: diss_t,diss_c
    real(kind=dp) :: rms_u,rms_v,rms_w,rms_t,rms_c,rms_cmt,prod_uw,prod_vw
  
    ! Calculate the vertical turbulent fluxes of the temperature and concentration
    flux_t = mean_flux(grid,state%t%p,state%w%p,params%horiz_t,params%strat_t,grid%z_vel)
    flux_c = mean_flux(grid,state%c%p,state%w%p,params%horiz_c,params%strat_c,grid%z_vel)

    ! Calculate the horizontal turbulent fluxes of the temperature and concentration
    hflux_t = mean_flux(grid,state%t%p,state%u%p,params%horiz_t,params%strat_t)
    hflux_c = mean_flux(grid,state%c%p,state%u%p,params%horiz_c,params%strat_c)

    ! Calculate the dissipation of the temperature and concentration
    diss_t = dissipation(grid,state%time,state%t%s)
    diss_c = dissipation(grid,state%time,state%c%s)

    ! Calculate the rms of velocity, temperature (T), concentration (C), and C - T
    rms_u = sqrt(sum_squares(grid,state%u%s))
    rms_v = sqrt(sum_squares(grid,state%v%s))
    rms_w = sqrt(sum_squares(grid,state%w%s))
    prod_uw = sum_product_p(grid,state%u%p,state%w%p) / grid%xyz_n
    prod_vw = sum_product_p(grid,state%v%p,state%w%p) / grid%xyz_n
    rms_t = sqrt(sum_squares(grid,state%t%s))
    rms_c = sqrt(sum_squares(grid,state%c%s))
    rms_cmt = sqrt(sum_squares(grid,state%c%s - state%t%s))
    
    ! Write to the diagnostics file
    if (get_id() .eq. 0) then
      write(unit,'(i7,14(",",1E21.15))') &
        & state%timestep,state%time,state%dt,rms_u,rms_v,rms_w,prod_uw,prod_vw, &
        & rms_t,rms_c,rms_cmt,flux_t,flux_c,diss_t,diss_c
    end if
  end subroutine

#include "output_netcdf.F90"

end module