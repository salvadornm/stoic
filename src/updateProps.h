#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, particleset  & vdmean, int p, double dt)
{
    //update density

    //update velocity
    vd.template getProp<i_velocity>(p)[0] += (vdmean.template getProp<i_velocity>(p)[0])*dt*.5;
    vd.template getProp<i_velocity>(p)[1] += (vdmean.template getProp<i_velocity>(p)[1])*dt*.5;
    vd.template getProp<i_velocity>(p)[2] += (vdmean.template getProp<i_velocity>(p)[2])*dt*.5;

}
#endif // _updateProps_h