module params_module
  use grid_module,only: dp,li,pi,grid_type
  implicit none 
  private
  public params_type,read_params
  public scheme_rkcn2,scheme_abbdf3

  !/**********************************************************************\
  !| PARAMETERS                                                           |
  !\**********************************************************************/
  integer, parameter :: scheme_rkcn2 = 1
  integer, parameter :: scheme_abbdf3 = 2

  !/**********************************************************************\
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type params_type
    !/********************************************************************\
    !| PARAMS_TYPE                                                        |
    !| -----------                                                        |
    !| The params type contains the parameters needed to describe the     |
    !| physical parameters of the simulation. It also contains the        |
    !| relevant time stepping information. These can be specified         |
    !| manually, but it is recommended to use the read_params subroutine. |
    !\********************************************************************/
    integer :: unit ! The file unit used for the parameter file

    integer :: max_steps = huge(1) ! The maximum number of timesteps taken
    real(kind=dp) :: max_time = huge(1.0_dp) ! The maximum time to reach
    real(kind=dp) :: start_dt = 1.0e-7_dp ! Starting dt
    real(kind=dp) :: max_dt = 1.0_dp ! Maximum dt
    real(kind=dp) :: limit_mult_dt = 0.5_dp ! Multiplier on dt limit
    real(kind=dp) :: max_increase_dt = 1.02_dp ! Factor dt can increase per step

    ! Like Botella and Peyret (2001), default RK-CN2 scheme to start
    character(len=10) :: start_scheme = "" ! The scheme used at start
    integer :: start_schemein = scheme_rkcn2 ! The scheme used at start
    character(len=10) :: scheme = "" ! The scheme used in the bulk run
    integer :: schemein = scheme_abbdf3 ! The scheme used in the bulk run
    
    real(kind=dp) :: diff_vel ! The viscosity

    real(kind=dp) :: strat_t = 0.0_dp ! The t stratification
    real(kind=dp) :: horiz_t = 0.0_dp ! The horizontal t gradient
    real(kind=dp) :: buoy_t ! The buoyancy coefficient of t
    real(kind=dp) :: diff_t ! The diffusivity of t

    real(kind=dp) :: strat_c = 0.0_dp ! The c stratification
    real(kind=dp) :: horiz_c = 0.0_dp ! The horizontal c gradient
    real(kind=dp) :: buoy_c ! The buoyancy coefficient of c
    real(kind=dp) :: diff_c ! The diffusivity of c

    ! This parameter is essentially the ratio of g_x/g_z
    ! Ordinarily, this is 0 unless the domain is inclined
    real(kind=dp) :: buoy_horiz_mult = 0.0_dp ! Horizontal buoyancy multiplier
  end type params_type

contains

  !/**********************************************************************\
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function get_scheme(scheme,default)
    !/********************************************************************\
    !| GET_SCHEME                                                         |
    !| ----------                                                         |
    !| Given the string representation of a scheme, return the integer    |
    !| representation of that scheme. Options are currently "abbdf3" and  |
    !| "rkcn2". If the scheme is an empty string, use the default scheme  |
    !| instead.                                                           |
    !|   scheme: the string representation of the scheme                  |
    !|     "abbdf3": Adams--Bashforth/Backward Difference Formula (3)     |
    !|     "rkcn2": Runge--Kutta/Crank--Nicolson (2)                      | 
    !|     "": use the scheme specified by default                        |
    !|   default: the integer representation of the default scheme        |
    !|   returns: the integer representation of the chosen scheme         |
    !\********************************************************************/
    implicit none

    character(len=10), intent(in) :: scheme
    integer, intent(in) :: default
    integer :: get_scheme

    select case (trim(scheme))
    case ("abbdf3")
      get_scheme = scheme_abbdf3
    case ("rkcn2")
      get_scheme = scheme_rkcn2
    case ("")
      get_scheme = default
    case default
      get_scheme = 0
    end select
  end function

  !/**********************************************************************\
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  ! TODO: Move argument reader to main
  subroutine read_params(grid,params)
    !/********************************************************************\
    !| READ_PARAMS                                                        |
    !| -----------                                                        |
    !| Read in the grid and parameter information from the parameters     |
    !| file. This subroutine will check to see if an alternate file name  |
    !| was specified as an argument at the command line; otherwise, it    |
    !| will look for a file named "parameters". The file should contain   |
    !| three namelists: grid, parameters, and io. Each should specify     |
    !| variables contained within the grid_type, params_type, and         |
    !| controllers_type respectively. For the grid namelist, it is        |
    !| expected that these variables are bare, but in parameters, they    |
    !| should appear as params%variable = value and in io, as             |
    !| controls%variable = value.                                         |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !\********************************************************************/
    use comm_module,only : get_id,abort_comm
    implicit none

    type(grid_type), intent(inout) :: grid
    type(params_type), intent(inout) :: params

    character(len=100) :: file
    integer :: ierror

    namelist /parameters/ params

    ! Check to see if a file name was specified at the command line
    call getarg(1,file)

    ! If no file was specified, use "parameters"
    if (trim(file) .eq. "") file = "parameters"

    ! Open the file and read in the grid information
    params%unit = 80 + get_id()
    open(unit=params%unit,file=trim(file),action="read")
    read(params%unit,fmt="(DT)",iostat=ierror) grid

    if (ierror .ne. 0) then 
      print*,"FATAL: While reading grid namelist, encountered error (",ierror,")"
      call abort_comm
    end if

    ! Read in the parameter information
    read(params%unit,nml=parameters,iostat=ierror)  

    if (ierror .ne. 0) then 
      print*,"FATAL: While reading input namelist, encountered error (",ierror,")"
      call abort_comm
    end if

    ! Convert the string representations of the schemes to integers
    params%start_schemein = get_scheme(params%start_scheme,params%start_schemein)
    params%schemein = get_scheme(params%scheme,params%schemein)

    ! Check that the chosen schemes are appropriate
    if (params%start_schemein .eq. scheme_abbdf3) then
      print*,"FATAL: Cannot use multi-step scheme as start_scheme."
      call abort_comm
    end if

    if (params%start_schemein .eq. 0) then
      print*,"FATAL: No start_scheme selected."
      call abort_comm
    end if

    if (params%schemein .eq. 0) then
      print*,"FATAL: No scheme selected."
      call abort_comm
    end if
  end subroutine

endmodule params_module

