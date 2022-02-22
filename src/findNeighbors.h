#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"

const double dp = 0.0085;//initial spacing between particles dp in the formulas
double W_dap = 0.0;

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, double & max_visc, double H)
{
    cout << " start nearest neighbor  " << endl;
    const double Eta2 = 0.01 * H*H;// Eta in the formulas
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);

    W_dap = 1.0/Wab(H/1.5);

    // For each particle ...
    while (part.isNext())
    {
        auto a = part.get();

        // Get all properties of the particle
        Point<3,double> xa = vd.getPos(a);
        //double massa = (vd.getProp<type>(a) == FLUID)?MassFluid:MassBound;
        //double rhoa = vd.getProp<rho>(a);
        //double Pa = vd.getProp<Pressure>(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            if (a.getKey() == b)    {++Np; continue;};// if (p == q) skip this particle
                
            // Get all properties of the particle
            //double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;
            Point<3,double> vb = vd.template getProp<i_velocity>(b);
            //double Pb = vd.getProp<Pressure>(b);
            //double rhob = vd.getProp<rho>(b);

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
                
                //max_visc=std::max(dot_rr2,max_visc);
                //vd.getProp<drho>(a) += massb*(dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));
            }
            ++Np;
        }
        ++part;
    }
    cout << " end nearest neighbor  " << endl;
}

#endif // _neighbors_h