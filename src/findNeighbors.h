#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"

template<typename CellList> void find_neighbors(particleset  & vd, particleset &vdmean, gradientset &dvdmean, gradientset &dvdmeanx, gradientset &dvdmeany, gradientset &dvdmeanz, CellList & NN)
{
    int n,ingh,ip;
    float iavg=0;

    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);

    // For each particle ...
    ip=1;
    //timer tsim;
    //tsim.start();
    
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

        dvdmean.template getProp<i_momentum>(a) = 0.0; 
        dvdmean.template getProp<i_rho>(a)  = 0.0;
        dvdmean.template getProp<i_energy>(a)   = 0.0;
        dvdmean.template getProp<i_temperature>(a) = 0.0;
        dvdmean.template getProp<i_pressure>(a)    = 0.0;
       
        dvdmeanx = dvdmean;
        dvdmeany = dvdmean;
        dvdmeanz = dvdmean;

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        
        // For each neighborhood particle
        ingh =0;
        while (Np.isNext() == true)
        {

            ingh++; 
            auto b = Np.get();
              
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            
            if (a.getKey() == b)    {++Np; continue;};// if (a == b) skip this particle
                
            // Get the distance between a and b
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
            vdmean.template getProp<i_energy>(a)      += W*vd.getProp<i_energy>(b);
            vdmean.template getProp<i_velocity>(a)[0] += W*vd.getProp<i_velocity>(b)[0];
            vdmean.template getProp<i_velocity>(a)[1] += W*vd.getProp<i_velocity>(b)[1];
            vdmean.template getProp<i_velocity>(a)[2] += W*vd.getProp<i_velocity>(b)[2];

            for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
            {
                // grad of particle property at a particle position
                //dvdmean.template getProp<i_momentum>(a) += DW.get(i)*vd.getProp<i_momentum>(b);
                //dvdmean.template getProp<i_density>(a)  += DW.get(i)*vd.getProp<i_density>(b);
                //dvdmean.template getProp<i_energy>(a)   += DW.get(i)*vd.getProp<i_energy>(b);
                dvdmean.template getProp<i_rho>(a)         += DW.get(i)*vd.getProp<i_rho>(b);
                dvdmean.template getProp<i_pressure>(a)    += DW.get(i)*vd.getProp<i_pressure>(b);

                switch(i){
                    case 0: dvdmeanx = dvdmean; break;
                    case 1: dvdmeany = dvdmean; break;
                    case 2: dvdmeanz = dvdmean; break;
                }
            }

            ++Np;
        }
        iavg = iavg + ingh;
        ++part;
    }

    cout << "i = " << ip << " neigh= "<< ingh  << endl;
    cout << " Avg number of neighbours=" << iavg/ip << endl;
    //tsim.stop();
    //std::cout << "Time: " << tsim.getwct() << std::endl;
}

#endif // _neighbors_h