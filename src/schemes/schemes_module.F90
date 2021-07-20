module schemes_module
  !/**********************************************************************\
  !| SCHEMES MODULE                                                       |
  !| --------------                                                       |
  !| This module contains the interface and selectors to run the various  |
  !| numerical schemes made available by the program. The primary         |
  !| available interface is the update subroutine, which steps the        |
  !| program forward one timestep.                                        |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type,scheme_rkcn2,scheme_abbdf3
  use fields_module, only: state_type,rhs_type,field_type
  implicit none
  private
  public update

contains

  !/**********************************************************************\
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function get_dt(grid,params,state,scheme)
    !/********************************************************************\
    !| GET_DT                                                             |
    !| ------                                                             |
    !| Calculates the length of the next timestep, using the limit        |
    !| appropriate to the chosen scheme.                                  |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   scheme: the numerical scheme                                     |
    !|   returns: the length of the next timestep                         |
    !\********************************************************************/
    use rkcn2_module, only: rkcn2_dt_limit
    use abbdf3_module, only: abbdf3_dt_limit
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    integer, intent(in) :: scheme
    real(kind=dp) :: get_dt

    real(kind=dp) :: diff

    diff = min(params%diff_t,params%diff_c,params%diff_vel)
    get_dt = state%dt

    ! Calculate the limit appropriate to scheme
    select case (scheme)
    case (scheme_rkcn2)
      get_dt = params%limit_mult_dt * rkcn2_dt_limit(grid,state,diff)
    case (scheme_abbdf3)
      get_dt = params%limit_mult_dt * abbdf3_dt_limit(grid,state,diff)
    end select

    ! The timestep can only increase by a fixed factor each step
    if (get_dt .le. params%max_increase_dt * state%dt) then
      get_dt = min(get_dt,params%max_dt)
    else
      get_dt = min(params%max_increase_dt * state%dt,params%max_dt)
    end if
  end function

  !/**********************************************************************\
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine update(grid,params,state,prior,rhs)
    !/********************************************************************\
    !| UPDATE                                                             |
    !| ------                                                             |
    !| Update the simulation to the next simulation state by choosing the |
    !| appropriate scheme from the parameters and executing run_scheme.   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous times       |
    !|   rhs: an array of rhs_type objects for prior times                |
    !\********************************************************************/
    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)

    state%timestep = state%timestep + 1

    ! Choose the appropriate scheme.
    if (state%timestep - state%restart .lt. 3) then
      call run_scheme(grid,params,state,prior,rhs,params%start_schemein)
    else 
      call run_scheme(grid,params,state,prior,rhs,params%schemein)
    end if
  end subroutine

  subroutine run_scheme(grid,params,state,prior,rhs,scheme)
    !/********************************************************************\
    !| RUN_SCHEME                                                         |
    !| ----------                                                         |
    !| Update the simulation to the next simulation state. If necessary,  |
    !| perform a remap operation.                                         |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous times       |
    !|   rhs: an array of rhs_type objects for prior times                |
    !|   scheme: the scheme to be run                                     |
    !\********************************************************************/
    use comm_module, only: get_id,abort_comm
    use grid_ops_module, only: get_incline
    use fields_module, only: step_time,remap,transform_all_to_p
    use rkcn2_module, only: rkcn2_run
    use abbdf3_module, only: abbdf3_run
    implicit none

    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)
    integer, intent(in) :: scheme

    real(kind=dp) :: incline(2),rtmp

    ! If the timestep hasn't been set yet, set it now
    if (state%dt .eq. 0.0_dp) state%dt = params%start_dt

    ! Cycle the arrays to the next time step and copy 
    call step_time(state,prior,rhs)

    ! Calculate the timestep limit, based on the scheme stability region
    state%dt = get_dt(grid,params,state,scheme)
    state%time = prior(1)%time + state%dt

    ! Check if the grid has tilted past the remap condition
    incline = get_incline(grid,state%time)
    
    if (abs(incline(1) * grid%z_length) .gt. grid%x_length / 2.0_dp) then
      grid%remap = .true.
    end if
    if (abs(incline(2) * grid%z_length) .gt. grid%y_length / 2.0_dp) then
      grid%remap = .true.
    end if

    ! Update the simulation in spectral space
    select case (scheme)
    case (scheme_rkcn2)
      call rkcn2_run(grid,params,state,prior,rhs)
    case (scheme_abbdf3)
      call abbdf3_run(grid,params,state,prior,rhs)
    end select

    ! Check if the simulation has crashed
    rtmp = sum(abs(state%all_s))
    if (rtmp .ne. rtmp) then
      print*,"FATAL: NaN detected, ending simulation."
      call abort_comm
    end if
    if (rtmp .gt. huge(1.0_dp) .or. rtmp .lt. -huge(1.0_dp)) then
      print*,"FATAL: +/- infinity detected, ending simulation."
      call abort_comm
    end if

    ! If necessary, remap the simulation
    if (grid%code_remap .and. grid%remap) then
      call remap(grid,state)
      grid%remap = .false.
      state%restart = state%timestep
    end if

    ! Update the physical space representation of the state
    call transform_all_to_p(grid,state)
  end subroutine

end module
