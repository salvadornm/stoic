#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"
#include"test.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN)
{
    int n,ingh,ip;
    float iavg=0;
    
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);

    // For each particle ...
    ip=0;   
    //timer tsim;
    //tsim.start();
    
    while (part.isNext())
    {
        auto a = part.get();
        ip++;

        // Get all properties of the particle a
        Point<3,double> xa = vd.getPos(a);

       // reset counters        
        for (size_t j = 0; j < 6.0 ; j++)
        { vd.template getProp<i_vdmean>(a)[j] = 0.0; }
        for (size_t j = 0; j < 3.0 ; j++)
        { 
            vd.template getProp<i_dvdmean>(a)[j][i_momentum]  = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_rho]  = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_energy]  = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_pressure]  = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_temperature]  = 0.0;
        }

        // Get an iterator of all the particles neighborhood of p
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
            double r = sqrt(r2);
            // r/simulation.H
            //cout << "r = " << r << endl;    // r = 0.028 is the cutoff. > = no impact
            
            Point<3,double> DW;
            double factor = DWab(dr,DW,r,false); // gradient kernel //
            double W = Wab(r); //kernel

            vd.template getProp<i_vdmean>(a)[i_rho] += W*vd.getProp<i_rho>(b);
            vd.template getProp<i_vdmean>(a)[i_temperature] += W*vd.getProp<i_temperature>(b);
            vd.template getProp<i_vdmean>(a)[i_pressure]    += W*vd.getProp<i_pressure>(b);
            vd.template getProp<i_vdmean>(a)[i_energy]      += W*vd.getProp<i_energy>(b);
            vd.template getProp<i_vdmean>(a)[i_velx] += W*vd.getProp<i_velocity>(b)[0];
            vd.template getProp<i_vdmean>(a)[i_vely] += W*vd.getProp<i_velocity>(b)[1];
            vd.template getProp<i_vdmean>(a)[i_velz] += W*vd.getProp<i_velocity>(b)[2];

            for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
            {
                // grad of particle property at a particle position
                vd.template getProp<i_dvdmean>(a)[i][i_rho]         += DW.get(i)*vd.getProp<i_rho>(b);
                vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += DW.get(i)*vd.getProp<i_pressure>(b);
                vd.template getProp<i_dvdmean>(a)[i][i_temperature] += DW.get(i)*vd.getProp<i_temperature>(b);
                vd.template getProp<i_dvdmean>(a)[i][i_momentum]    += DW.get(i)*(vd.getProp<i_velocity>(b)[i]*vd.getProp<i_rho>(b));
                //vd.template getProp<i_dvdmean>(a)[i][i_energy]      += W*vd.getProp<i_energy>(b);

            }
            ++Np;
        }
        
        //vd.template getProp<i_vdmean>(a)[i_velocity] /= ingh;
        //vd.template getProp<i_vdmean>(a)[i_vely] /= ingh;
        //vd.template getProp<i_vdmean>(a)[i_velz] /= ingh;
        cout << "particle = " << ip << " neigh= "<< ingh  << endl;
        iavg += ingh;
        ++part;
    }

    cout << " Avg number of neighbours=" << iavg/ip << endl;
    //tsim.stop();
    //std::cout << "Time: " << tsim.getwct() << std::endl;
}

#endif // _neighbors_h