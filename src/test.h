#ifndef _test_h
#define _test_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include "kernel.h"
#include "calculations.h"

using namespace std;

void kernel_test(double H, Point<3,double> dr);
void output_kernel(double r, double h);
void output_vd(particleset  & vd, int p);
void output_properties(double mom_p, double drho, double rho_new, double Au_p, double dWien);
void output_bc_props(vector<double> vel, Point <3,double> pos, Point <3,double> pos_new, Point <3,double> pos_wall, engine eng, Point<3,double> & psi);
void output_energy_props(particleset &vd, int p, double dh, double dvisc, double edensity_p, double edensity_new, double energy_new);
void vary_initialization(particleset &vd, Cfd simulation, int key);
void limit_velocity(particleset &vd, int key, int i);

template<typename CellList> int stateOfNeighbors(particleset  & vd, CellList & NN)
{    auto part = vd.getDomainIterator();

    while(part.isNext())
    {
        auto a = part.get();

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a))); 
        cout << a.getKey() << " position " << vd.getPos(a)[0] << ", " << vd.getPos(a)[1] << ", "<< vd.getPos(a)[2] << ", " << endl;   

        while (Np.isNext() == true)
        {
            auto b = Np.get();

            Point<3,double> xa = vd.getPos(a);
            Point<3,double> xb = vd.getPos(b);
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);
            double r = sqrt(r2);

            if (r < H) {
                cout << b << ":";
                cout << "r = " << r << " , ";                 
            }
            // r/simulation.H

            ++Np;
        }
        cout << "}" << endl;
        ++part;
    }
    return 0;
}
    
#endif