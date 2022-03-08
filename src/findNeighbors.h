#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"

const double dp = 0.0085;//initial spacing between particles dp in the formulas
double mass_p = 1.0; //m_tot/N; 

template<typename CellList> void find_neighbors(particleset  & vd, particleset  & vdmean, CellList & NN, double H)
{
    //cout << " start nearest neighbor  " << endl;
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
        double tempa = vd.getProp<i_temperature>(a); 
        double Pa = vd.getProp<i_pressure>(a);
        double rhoa = vd.getProp<i_rho>(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

        //reset counters
        vdmean.template getProp<i_rho>(a)         = 0.0;
        vdmean.template getProp<i_temperature>(a) = 0.0;
        vdmean.template getProp<i_pressure>(a)    = 0.0;
        vdmean.template getProp<i_velocity>(a)[0] = 0.0;
        vdmean.template getProp<i_velocity>(a)[1] = 0.0;
        vdmean.template getProp<i_velocity>(a)[2] = 0.0;

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle, (assuming all particles are fluid particles)
        while (Np.isNext() == true)
        {
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            
            if (a.getKey() == b)    {++Np; continue;};// if (p == q) skip this particle
                
            // Get all properties of the particle
            
            Point<3,double> vb = vd.template getProp<i_velocity>(b);
            double Pb = vd.getProp<i_pressure>(b);
            double rhob = vd.getProp<i_rho>(b);
            double temp_b = vd.getProp<i_temperature>(b); 

            // Get the distance between p and q
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);

            // If the particles interact ... (if r2 is less than the cutoff radius)
            if (r2 < 4.0*H*H)
            //if (r2 < 0.1)
            {
                double r = sqrt(r2);
                Point<3,double> dv = va - vb;
                Point<3,double> DW;
                DWab(dr,DW,r,false);    //compute gradient but do not use gradient...
                double W = Wab(r);  //weighting
                
                // else if it is a fluid particle: 
                // collect the mean of the properties of each negihbor, weighted
                vdmean.template getProp<i_rho>(a)         += W*rhob;
                vdmean.template getProp<i_temperature>(a) += W*temp_b;
                vdmean.template getProp<i_pressure>(a)    += W*Pb;
                vdmean.template getProp<i_velocity>(a)[0] += DW.get(0)* vb[0];
                vdmean.template getProp<i_velocity>(a)[1] += DW.get(1)* vb[1];
                vdmean.template getProp<i_velocity>(a)[2] += DW.get(2)* vb[2];

            }
            ++Np;
        }
        ++part;
    }
}

#endif // _neighbors_h