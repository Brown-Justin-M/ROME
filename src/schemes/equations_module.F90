module equations_module
  !/**********************************************************************\
  !| EQUATIONS MODULE                                                     |
  !| ----------------                                                     |
  !| This module contains the subroutines specific to the equations. This |
  !| includes the calculation of the explicit and implicit terms of the   |
  !| equations and the general form for the timestepping limit.           |
  !\**********************************************************************/
    use grid_module, only: dp,li,grid_type
    use params_module, only: params_type
    use fields_module, only: state_type,rhs_type,field_type
    implicit none
    private
    public dt_limit,explicit_rhs,explicit_implicit_rhs,implicit_update
  
contains

  !/**********************************************************************\
  !| FUNCTIONS                                                            |
  !\**********************************************************************/ 
  function dt_limit(grid,state,diff,a,b)
    !/********************************************************************\
    !| DT_LIMIT                                                           |
    !| --------                                                           |
    !| Calculates the limit on dt to ensure numerical stability. This is  |
    !| adapted from the 1D version described in 4.1.1(e) of Peyret (2002) |
    !| Spectral Methods for Incompressible Viscous Flow, Equation 4.97.   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the simulation data        |
    !|   diff: the smallest diffusion coeffient for the equations         |
    !|   returns: the length of the next timestep                         |
    !\********************************************************************/
    use comm_module, only: collect_max,collect_array_max
    use grid_ops_module, only: get_incline
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(in) :: state
    real(kind=dp), intent(in) :: diff
    real(kind=dp), intent(in) :: a,b
    real(kind=dp) :: dt_limit
    
    real(kind=dp) :: hkx,hky,hkz,ksquared,incline(2)
    real(kind=dp) :: umax(3),advect,rtmp
    integer :: i,j,k

    ! Get the maximum velocity on the grid
    umax(1) = collect_array_max(abs(state%u%p))
    umax(2) = collect_array_max(abs(state%v%p))
    umax(3) = collect_array_max(abs(state%w%p))

    incline = get_incline(grid,state%time)
    
    ! Calculate the maximum denominator of Equation 4.97
    rtmp = epsilon(0.0_dp)
    do j = 1,grid%y_s_count
      hky = abs(grid%lky(j))
      do i = 1,grid%x_s_count
        hkx = abs(grid%lkx(i))
        do k = 1,grid%z_s_count
          hkz = abs(grid%lkz(k) - hkx * incline(1) - hky * incline(2))
          ksquared = hkx**2 + hky**2 + hkz**2
          advect = hkx * umax(1) + hky * umax(2) + hkz * umax(3)
          rtmp = max(rtmp, a * advect - diff * ksquared)
        end do
      end do
    end do
    rtmp = collect_max(rtmp)

    if (rtmp .le. 0.0_dp) then
      dt_limit = huge(0.0_dp)
    else
      dt_limit = a * b / rtmp
    end if
  end function

  !/**********************************************************************\
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine explicit_rhs(grid,params,prior,state,rhs,b_rhs)
    !/********************************************************************\
    !| EXPLICICT_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the explicit parts of the right hand sides of the        |
    !| equations. Store them in rhs. Note the warning in explicit_rhs_vel |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous time        |
    !|   rhs: an rhs_type object for the current time                     |
    !|   b_rhs: the multiplier on the rhs in the final scheme             |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    type(rhs_type), intent(inout) :: rhs
    real(kind=dp), optional :: b_rhs

    real(kind=dp) :: b_rhsin

    b_rhsin = state%dt
    if (present(b_rhs)) b_rhsin = b_rhs

    rhs%all = 0.0_dp

    call explicit_rhs_vel(grid,params,prior,state,rhs%vel,b_rhsin)
    call explicit_rhs_t(grid,params,prior,rhs%t)
    call explicit_rhs_c(grid,params,prior,rhs%c)
  end subroutine

  subroutine explicit_implicit_rhs(grid,params,prior,state,rhs)
    !/********************************************************************\
    !| EXPLICIT_IMPLICIT_RHS                                              |
    !| ---------------------                                              |
    !| Calculate the implicit part of the equations, but treat them as    |
    !| though they were explicit. This is useful for the Crank--Nicolson  |
    !| scheme.                                                            |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous time        |
    !|   rhs: an rhs_type object for the current time                     |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    type(rhs_type), intent(out) :: rhs

    rhs%all = 0.0_dp

    call diffusion_rhs(grid,params,state,prior%u,rhs%u,params%diff_vel)
    call diffusion_rhs(grid,params,state,prior%v,rhs%v,params%diff_vel)
    call diffusion_rhs(grid,params,state,prior%w,rhs%w,params%diff_vel)
    call diffusion_rhs(grid,params,state,prior%t,rhs%t,params%diff_t)
    call diffusion_rhs(grid,params,state,prior%c,rhs%c,params%diff_c)
  end subroutine
  
  subroutine explicit_rhs_t(grid,params,state,out)
    !/********************************************************************\
    !| EXPLICICT_RHS_T                                                    |
    !| ---------------                                                    |
    !| Calculate the explicit parts of the right hand side of the         |
    !| temperature equation. Store it in out. Note that the state is the  |
    !| state for which the right hand side should be evaluated            |
    !| (typically the previous state).                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) :: out(:,:,:)

    call advection_rhs(grid,params,state,state%t,out)

    if (params%strat_t .ne. 0.0_dp) then
      out = out - params%strat_t * state%w%s
    end if

    if (params%horiz_t .ne. 0.0_dp) then
      out = out - params%horiz_t * state%u%s
    end if

    ! In the absence of a net velocity, the mean doesn't change
    if (get_id() .eq. 0 .and. grid%z_vel .eq. 0.0_dp) then
      out(1,1,1) = (0.0_dp,0.0_dp)
    end if
  end subroutine

  subroutine explicit_rhs_c(grid,params,state,out)
    !/********************************************************************\
    !| EXPLICICT_RHS_C                                                    |
    !| ---------------                                                    |
    !| Calculate the explicit parts of the right hand side of the         |
    !| concentration equation. Store it in out. Note that the state is    |
    !| the state for which the right hand side should be evaluated        |
    !| (typically the previous state).                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) :: out(:,:,:)

    call advection_rhs(grid,params,state,state%c,out)

    if (params%strat_c .ne. 0.0_dp) then
      out = out - params%strat_c * state%w%s
    end if

    if (params%horiz_t .ne. 0.0_dp) then
      out = out - params%horiz_c * state%u%s
    end if

    ! In the absence of a net velocity, the mean doesn't change
    if (get_id() .eq. 0 .and. grid%z_vel .eq. 0.0_dp) then
      out(1,1,1) = (0.0_dp,0.0_dp)
    end if
  end subroutine

  subroutine explicit_rhs_vel(grid,params,prior,state,out,b_rhs)
    !/********************************************************************\
    !| EXPLICIT_RHS_VEL                                                   |
    !| ----------------                                                   |
    !| **WARNING** For now, the velocity in "state%vel%s" is expected to  |
    !| already include the scheme estimates from prior times.             |
    !| Calculate the explicit parts of the right hand side of the         |
    !| velocity equation. Store it in out. Note that the state is the     |
    !| state for which the right hand side should be evaluated (typically |
    !| the previous state). This also ensures that the pressure component |
    !| of the right hand side makes state%vel + b_rhs * out               |
    !| incompressible.                                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !|   b_rhs: the multiplier on the rhs in the final scheme             |
    !\********************************************************************/
    use grid_module, only: allocate_p
    use grid_ops_module, only: cross,remove_div,get_shear
    use comm_module, only: get_id
    use transform_module, only: transform_p_to_s,transform_s_to_p 
    use mpi
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) ::  out(:,:,:,:)
    real(kind=dp) :: b_rhs

    real(kind=dp) :: shear(2)
    real(kind=dp), pointer :: work_p(:,:,:)

    call allocate_p(grid,work_p)

    ! Evaluate the RHS of the u-component
    call cross(prior%vel%p,prior%crl%p,work_p,1)

    ! Add in the mean shear component of advection
    shear = get_shear(grid,prior%time)
    work_p = work_p - prior%w%p * shear(1)

    ! Add in buoyancy due to horizontal gravity component
    if (params%buoy_horiz_mult .ne. 0.0_dp) then
      work_p = work_p + params%buoy_t * params%buoy_horiz_mult * state%t%p 
      work_p = work_p - params%buoy_c * params%buoy_horiz_mult * state%c%p
    end if
    call transform_p_to_s(work_p,out(:,:,:,1))

    ! Evaluate the RHS of the v-component
    call cross(prior%vel%p,prior%crl%p,work_p,2)

    ! Add in the mean shear component of advection
    work_p = work_p - prior%w%p * shear(2)
    call transform_p_to_s(work_p,out(:,:,:,2))

    ! Evaluate the RHS of the w-component
    call cross(prior%vel%p,prior%crl%p,work_p,3)

    ! Add in buoyancy due to vertical gravity component
    work_p = work_p + params%buoy_t*prior%t%p
    work_p = work_p - params%buoy_c*prior%c%p

    call transform_p_to_s(work_p,out(:,:,:,3))
    deallocate(work_p)

    ! Add pressure term to rhs by removing the divergence of the next step
    call remove_div(grid,state%time,b_rhs * out + state%vel%s,out,b_rhs)

    ! The change in the mean velocity is 0 by construction
    if (get_id() .eq. 0) out(1,1,1,:) = (0.0_dp,0.0_dp)
  end subroutine

  subroutine advection_rhs(grid,params,state,field,out)
    !/********************************************************************\
    !| ADVECTION_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the advection term in the right hand side of an equation.|
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use transform_module, only: transform_p_to_s
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)

    call allocate_p(grid,work_p) 
    call allocate_s(grid,work_s)

    work_p = field%p * state%u%p
    call transform_p_to_s(work_p,work_s)
    call deriv_x(grid,work_s,work_s)
    out = out - work_s

    work_p = field%p * state%v%p
    call transform_p_to_s(work_p,work_s)
    call deriv_y(grid,work_s,work_s)
    out = out - work_s

    work_p = field%p * state%w%p
    call transform_p_to_s(work_p,work_s)
    call deriv_z(grid,state%time,work_s,work_s)
    out = out - work_s
    
    deallocate(work_p,work_s)
  end subroutine

  subroutine diffusion_rhs(grid,params,state,field,out,diff)
    !/********************************************************************\
    !| DIFFUSION_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the diffusion term in the right hand side of an equation.|
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   out: an array used to contain the rhs information                |
    !|   diff: the diffusion coefficient                                  |
    !\********************************************************************/
    use grid_ops_module, only: get_incline
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)
    real(kind=dp), intent(in) :: diff

    real(kind=dp)  :: incline(2)
    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k

    incline = get_incline(grid,state%time)

    do j = 1,grid%y_s_count
      hky = grid%lky(j) 
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k) - hky * incline(2) - hkx * incline(1)
          ksquared = hkx**2 + hky**2 + hkz**2
          out(k,i,j) = out(k,i,j) - ksquared * diff * field%s(k,i,j)
        end do
      end do
    end do
  end subroutine

  subroutine implicit_update(grid,params,state,a0)
    !/********************************************************************\
    !| IMPLICIT_UPDATE                                                    |
    !| ---------------                                                    |
    !| Once the right hand side has been calculated and the effects of    |
    !| the prior steps have been included in state, use this to perform   |
    !| the final update to the new time with                              |
    !|     a0 * u(n+1) + dt * I(n+1) = b_rhs * RHS + PRIOR_TERMS          |
    !| This subroutine updates every variable in state.                   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   a0: the multiplier on the new value in the final scheme          |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    real(kind=dp), optional :: a0

    real(kind=dp) a0in

    a0in = 1.0_dp
    if (present(a0)) a0in = a0

    call diffuse(grid,state,state%u%s,state%u%s,params%diff_vel,a0in)
    call diffuse(grid,state,state%v%s,state%v%s,params%diff_vel,a0in)
    call diffuse(grid,state,state%w%s,state%w%s,params%diff_vel,a0in)
    call diffuse(grid,state,state%t%s,state%t%s,params%diff_t,a0in)
    call diffuse(grid,state,state%c%s,state%c%s,params%diff_c,a0in)
  end subroutine

  subroutine diffuse(grid,state,in,out,diff,a0)
    !/********************************************************************\
    !| DIFFUSE                                                            |
    !| -------                                                            |
    !| Once the right hand side has been calculated and the effects of    |
    !| the prior steps have been included in state, use this to perform   |
    !| the final update to the new time with                              |
    !|     a0 * u(n+1) + dt * I(n+1) = b_rhs * RHS + PRIOR_TERMS          |
    !| This subroutine updates one variable.                              |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the simulation data        |
    !|   in: the input array containing the prior time information        |
    !|   out: the output array to contain the new variable state          |
    !|   diff: the diffusion coefficient                                  |
    !|   a0: the multiplier on the new value in the final scheme          |
    !\********************************************************************/
    use grid_ops_module, only: get_incline
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)
    real(kind=dp), intent(in) :: diff
    real(kind=dp), intent(in) :: a0

    real(kind=dp)  :: incline(2)
    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k

    incline = get_incline(grid,state%time)

    do j = 1,grid%y_s_count
      hky = grid%lky(j) 
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k) - hky * incline(2) - hkx * incline(1)
          ksquared = hkx**2 + hky**2 + hkz**2
          rtmp = a0 + state%dt * ksquared * diff
          out(k,i,j) = in(k,i,j) / rtmp
        end do
      end do
    end do
  end subroutine

end module equations_module