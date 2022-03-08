#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, particleset  & vdmean, int p, double dt)
{
    //update density
    //calculate pressure
    //calculate pressure velocity term...
    double time_turb = 0.1;
    double eps_turb = 0.1;

    //update velocity transport term
    // drift term
    double du = (vdmean.template getProp<i_velocity>(p)[0] -  vd.template getProp<i_velocity>(p)[0]);
    // stoic
    double dW = 0.0; // eps_turb*distribution(generator)*sqrt(dt);

    vd.template getProp<i_velocity>(p)[0] += du*dt/(time_turb + dt) + dW;
    vd.template getProp<i_velocity>(p)[1] += (vdmean.template getProp<i_velocity>(p)[1])*dt;
    vd.template getProp<i_velocity>(p)[2] += (vdmean.template getProp<i_velocity>(p)[2])*dt;

    //check continuity?
    //update velocity

}
#endif // _updateProps_h