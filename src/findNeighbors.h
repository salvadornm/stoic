
#include "Vector/vector_dist.hpp"
#include <math.h>
#include <iostream>
#include "kernel.h"
#include "global.h"



template<typename CellList> inline void find_neighbors(particleset  & vd, CellList & NN, double & max_visc, double H)
{
    const double Eta2 = 0.01 * H*H;// Eta in the formulas
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);

    double W_dap = 1.0/Wab(H/1.5);

    const int i_velocity = 4;

    // For each particle ...
    while (part.isNext())
    {
        auto a = part.get();

        // Get all properties of the particle
        Point<3,double> xa = vd.getPos(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            if (a.getKey() == b)    {++Np; continue;};// if (p == q) skip this particle
                
            // Get all properties of the particle
            Point<3,double> vb = vd.template getProp<i_velocity>(b);
            // Get the distance between p and q
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);

            // If the particles interact ...
            if (r2 < 4.0*H*H)
            {
                // ... calculate delta rho
                double r = sqrt(r2);
                Point<3,double> dv = va - vb;
                Point<3,double> DW;
                DWab(dr,DW,r,false);
                const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
                const double dot_rr2 = dot/(r2+Eta2);                
            }
            ++Np;
        }
        ++part;
    }
}