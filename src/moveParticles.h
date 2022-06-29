#ifndef _moveParticles_h
#define _moveParticles_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void moveParticles(particleset  & vd, int p, double dt, engine eng)
{
    //vd.getPos(p)[0] += vd.template getProp<i_velocity>(p)[0]*dt;
    //vd.getPos(p)[1] += vd.template getProp<i_velocity>(p)[1]*dt;
    //vd.getPos(p)[2] += vd.template getProp<i_velocity>(p)[2]*dt;

    double pos_x, pos_y, pos_z;
    double vel_x, vel_y, vel_z;

    // update temporary particle position
    pos_x = vd.getPos(p)[0] + vd.template getProp<i_velocity>(p)[0]*dt;
    pos_y = vd.getPos(p)[1] + vd.template getProp<i_velocity>(p)[1]*dt;
    pos_z = vd.getPos(p)[2] + vd.template getProp<i_velocity>(p)[2]*dt;

    //check boundary conditions
    //if 1 < p < 0 then the position of the particle is out of the box
    if (pos_x < 0 || pos_x > eng.bore)
        vel_x = -vd.template getProp<i_velocity>(p)[0];
        // reflect position over the boundary
        



    //update to new particle position
    vd.getPos(p)[0] = pos_x;
    vd.getPos(p)[1] = pos_y;
    vd.getPos(p)[2] = pos_z;
}
#endif // _moveParticles_h
