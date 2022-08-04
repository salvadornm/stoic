#ifndef _moveParticles_h
#define _moveParticles_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void moveParticles(particleset  & vd, int p, double dt, engine eng)
{
    Point <3,double> psi {0.0,0.0,0.0};  //distance out of bounds

    // update temporary particle position
    Point <3,double> pos_zero {vd.getPos(p)[0],vd.getPos(p)[1],vd.getPos(p)[2]};
    Point <3,double> pos_new {pos_zero[0] + vd.template getProp<i_velocity>(p)[0]*dt, pos_zero[1] + vd.template getProp<i_velocity>(p)[1]*dt, pos_zero[2] + vd.template getProp<i_velocity>(p)[2]*dt};
    vector <double> vel {vd.template getProp<i_velocity>(p)[0],vd.template getProp<i_velocity>(p)[1],vd.template getProp<i_velocity>(p)[2]};

    //check boundary conditions
    int bc_flag = inCylinder(vd, pos_zero, pos_new, p, eng, psi);
    
    //update to new particle position
    if (bc_flag == 1){
        vd.getPos(p)[0] = pos_new[0];
        vd.getPos(p)[1] = pos_new[1];
        vd.getPos(p)[2] = pos_new[2];
        cout << "IN BOUNDS!" << endl;
        return;
    }

    while (bc_flag == 0){

        topBound(pos_new, vel, eng);
        pistonBound(pos_new, vel, eng);
        sideBC(vd, pos_zero, pos_new, vel, p, eng, dt);

        //need to determine which side it hits first... use psi...
        /*
        if (psi[0] > psi[2] || psi[1] > psi[2]){
            sideBC(vd, pos_zero, pos_new, vel, p, eng, dt);
            topBound(pos_new, vel, eng);
            pistonBound(pos_new, vel, eng);
        }
        else{
            topBound(pos_new, vel, eng);
            pistonBound(pos_new, vel, eng);
            sideBC(vd, pos_zero, pos_new, vel, p, eng, dt);
        }
        */

        bc_flag = inCylinder(vd, pos_zero, pos_new, p, eng, psi);
        if (bc_flag == 1){
            cout << "IN BOUNDS!" << endl;
        }
    }
    
    vd.getPos(p)[0] = pos_new[0];
    vd.getPos(p)[1] = pos_new[1];
    vd.getPos(p)[2] = pos_new[2];
    vd.template getProp<i_velocity>(p)[0] = vel[0];
    vd.template getProp<i_velocity>(p)[1] = vel[1];
    vd.template getProp<i_velocity>(p)[2] = vel[2];   
}
#endif // _moveParticles_h
