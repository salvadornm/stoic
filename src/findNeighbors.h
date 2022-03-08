#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"

template<typename CellList> void find_neighbors(particleset  & vd, particleset &vdmean, particleset &dvdmeanx, CellList & NN, double H)
{
    int n,ingh,ip;
    float iavg=0;

    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);

  //  cout << " star looping  " << endl;

    // For each particle ...
    ip=1;
    while (part.isNext())
    {
        auto a = part.get();
        ip++;

        // Get all properties of the particle a
        Point<3,double> xa = vd.getPos(a);

        // reset counters
        vdmean.template getProp<i_rho>(a)         = 0.0;
        vdmean.template getProp<i_temperature>(a) = 0.0;
        vdmean.template getProp<i_pressure>(a)    = 0.0;
        vdmean.template getProp<i_velocity>(a)[0] = 0.0;
        vdmean.template getProp<i_velocity>(a)[1] = 0.0;
        vdmean.template getProp<i_velocity>(a)[2] = 0.0;

        dvdmeanx.template getProp<i_velocity>(a)[0] = 0.0;
        dvdmeanx.template getProp<i_velocity>(a)[1] = 0.0;
        dvdmeanx.template getProp<i_velocity>(a)[2] = 0.0;

        
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle, (assuming all particles are fluid particles)
        ingh =0;
        while (Np.isNext() == true)
        {

            ingh++; 
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            
            if (a.getKey() == b)    {++Np; continue;};// if (p == q) skip this particle
                
            // Get the distance between p and q
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);
            // ... calculate delta rho
            double r = sqrt(r2);
            Point<3,double> DW;
            DWab(dr,DW,r,false); // gradient kernel //
            double W = Wab(r); //kernel

            vdmean.template getProp<i_rho>(a)         += W*vd.getProp<i_rho>(b);
            vdmean.template getProp<i_temperature>(a) += W*vd.getProp<i_temperature>(b);
            vdmean.template getProp<i_pressure>(a)    += W*vd.getProp<i_pressure>(b);
            vdmean.template getProp<i_velocity>(a)[0] += W*vd.getProp<i_velocity>(b)[0];
            vdmean.template getProp<i_velocity>(a)[1] += W*vd.getProp<i_velocity>(b)[1];
            vdmean.template getProp<i_velocity>(a)[2] += W*vd.getProp<i_velocity>(b)[2];

            // grad x of partciel property at a particle position
            dvdmeanx.template getProp<i_velocity>(a)[0] += DW.get(0)*vd.getProp<i_velocity>(b)[0];
            dvdmeanx.template getProp<i_velocity>(a)[1] += DW.get(0)*vd.getProp<i_velocity>(b)[1];
            dvdmeanx.template getProp<i_velocity>(a)[2] += DW.get(0)*vd.getProp<i_velocity>(b)[2];

            ++Np;
        }
        iavg = iavg + ingh;
        ++part;
    }
    //cout << "i = " << ip << " neigh= "<< ingh  << endl;
    cout << " Avg number of neighbours=" << iavg/ip << endl;
}

#endif // _neighbors_h