/***********************************************************************
   vprofile.c                                                          
   UDF for specifying steady-state velocity profile boundary condition 
************************************************************************/
#include "udf.h"

DEFINE_PROFILE(inlet_velocity_1, thread, position) 
{
  real x[ND_ND];		/* this will hold the position vector */
  face_t f;

  begin_f_loop(f, thread)
    {
      F_CENTROID(x,f,thread);
      F_PROFILE(f, thread, position) = 2.0*1.2732*(1-pow((x[1]-0)/5.0e-4,2.0)-pow((x[2]-0.005712)/5.0e-4,2.0));
    }
  end_f_loop(f, thread)
}
