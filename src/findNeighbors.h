#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"

const double dp = 0.0085;//initial spacing between particles dp in the formulas
double mass_p = 1.0; //m_tot/N; 

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, double & max_visc, double H)
{
    cout << " start nearest neighbor  " << endl;
    const double Eta2 = 0.01 * H*H;// Eta in the formulas
    const double W_dap = 1.0/Wab(H/1.5);

    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);


    // For each particle ...
    while (part.isNext())
    {
        auto a = part.get();

        // Get all properties of the particle
        Point<3,double> xa = vd.getPos(a);
        double massa = mass_p; //(vd.getProp<type>(a) == FLUID)?MassFluid:MassBound;
        double rhoa = vd.getProp<i_rho>(a);
        double Pa = vd.getProp<i_pressure>(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

        //reset force counters
        double forcex = 0.0;
        double forcey = 0.0;
        double forcez = -9.81;
        double drho = 0.0;

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle, (assuming all particles are fluid particles)
        while (Np.isNext() == true)
        {
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            
            if (a.getKey() == b)    {++Np; continue;};// if (p == q) skip this particle
                
            // Get all properties of the particle
            double massb = mass_p; //(vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;
            Point<3,double> vb = vd.template getProp<i_velocity>(b);
            double Pb = vd.getProp<i_pressure>(b);
            double rhob = vd.getProp<i_rho>(b);

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

                // //if it is a boundary particle: calculate drho
                // const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
                // const double dot_rr2 = dot/(r2+Eta2);
                // max_visc=std::max(dot_rr2,max_visc);
                // drho += massb*(dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));

                // else if it is a fluid particle: 
                double factor = - massb*((vd.getProp<i_pressure>(a) + vd.getProp<i_pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb,W_dap) + viscous(dr,r2,dv,rhoa,rhob,massb,max_visc));
                forcex += factor * DW.get(0);
                forcey += factor * DW.get(1);
                forcez += factor * DW.get(2);
                drho += massb*(dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));
                
            }
            ++Np;
        }
        ++part;
    }
    cout << " end nearest neighbor  " << endl;
}

#endif // _neighbors_h