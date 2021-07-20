module input_module
  !/**********************************************************************\
  !| INPUT MODULE                                                         |
  !| ------------                                                         |
  !| This module contains the interface to read an input file into a      |
  !| state_type derived type. The input module supports inputs in both    |
  !| physical and spectral space.                                         |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type
  use output_module, only: nc_ids_type,check
  implicit none
  private
  public read_input

contains

#include "input_netcdf.F90"

end module