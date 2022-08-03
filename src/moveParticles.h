#ifndef _moveParticles_h
#define _moveParticles_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"

void moveParticles(particleset  & vd, int p, double dt, engine eng)
{
    double pos_x, pos_y, pos_z;
    double vel_x, vel_y, vel_z;
    vector <double> psi {0.0,0.0,0.0};  //distance out of bounds

    // update temporary particle position
    pos_x = vd.getPos(p)[0] + vd.template getProp<i_velocity>(p)[0]*dt;
    pos_y = vd.getPos(p)[1] + vd.template getProp<i_velocity>(p)[1]*dt;
    pos_z = vd.getPos(p)[2] + vd.template getProp<i_velocity>(p)[2]*dt;
    vector <double> pos {pos_x,pos_y,pos_z};
    vel_x = vd.template getProp<i_velocity>(p)[0];
    vel_y = vd.template getProp<i_velocity>(p)[1];
    vel_z = vd.template getProp<i_velocity>(p)[2];
    vector <double> vel {vel_x,vel_y,vel_z};

    //check boundary conditions
    int bc_flag = inCylinder(vd, pos_x, pos_y, pos_z, p, eng, psi);
    
    //update to new particle position
    if (bc_flag == 1){
        vd.getPos(p)[0] = pos_x;
        vd.getPos(p)[1] = pos_y;
        vd.getPos(p)[2] = pos_z;
        cout << "IN BOUNDS!" << endl;
        return;
    }

    cout << "Vel_z: " << vel_z;
    while (bc_flag == 0){
        topBound(vd, vel_z, pos_z, p, eng);
        pistonBound(vd, vel_z, pos_z, p, eng);
        cout << " updated vel_z: " << vel_z << endl;
        sideBC(vd, vel_x, vel_y, pos_x, pos_y, p, eng, dt);

        bc_flag = inCylinder(vd, pos_x, pos_y, pos_z, p, eng, psi);
        if (bc_flag == 1){
            vd.getPos(p)[0] = pos_x;
            vd.getPos(p)[1] = pos_y;
            vd.getPos(p)[2] = pos_z;
            vd.template getProp<i_velocity>(p)[0] = vel_x;
            vd.template getProp<i_velocity>(p)[1] = vel_y;
            vd.template getProp<i_velocity>(p)[2] = vel_z;
            cout << "IN BOUNDS!" << endl;
        return;
        }
    }
}
#endif // _moveParticles_h
