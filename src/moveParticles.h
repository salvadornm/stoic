#ifndef _moveParticles_h
#define _moveParticles_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void moveParticles(particleset  & vd, particleset  & vdmean, int p, double dt)
{
    vd.getPos(p)[0] += vd.template getProp<i_velocity>(p)[0]*dt;
    vd.getPos(p)[1] += vd.template getProp<i_velocity>(p)[1]*dt;
    vd.getPos(p)[2] += vd.template getProp<i_velocity>(p)[2]*dt;

    //check boundary conditions
    
}
#endif // _moveParticles_h
