#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"
#include"test.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN){
    int n,ingh,ip;
    float iavg=0;
    
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    //stateOfNeighbors(vd, NN);

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
        Point<3,double> va = vd.getProp<i_velocity>(a);

       // reset counters        
        for (size_t j = 0; j < 7.0 ; j++)
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
        double tot_W = 0;
        
        // For each neighborhood particle
        ingh =0;
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            Point<3,double> vb = vd.getProp<i_velocity>(b);
            
            if (a.getKey() == b)    {++Np; continue;};// if (a == b) skip this particle
                
            // Get the distance between a and b
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);  //norm2 = (sum of the squares)
            double r = sqrt(r2);
                        
            //if the particles interact...
            if (r < H) {
                Point<3,double> DW;
                double factor = DWab(dr,DW,r,false); // gradient kernel //
                double W = Wab(r); //kernel
                tot_W += W;            
                
                vd.template getProp<i_vdmean>(a)[i_rho] += W*vd.getProp<i_rho>(b);
                vd.template getProp<i_vdmean>(a)[i_temperature] += W*vd.getProp<i_temperature>(b);
                vd.template getProp<i_vdmean>(a)[i_pressure]    += W*vd.getProp<i_pressure>(b);
                vd.template getProp<i_vdmean>(a)[i_energy]      += W*vd.getProp<i_energy>(b);
                vd.template getProp<i_vdmean>(a)[i_velx] += W*vd.getProp<i_velocity>(b)[0];
                vd.template getProp<i_vdmean>(a)[i_vely] += W*vd.getProp<i_velocity>(b)[1];
                vd.template getProp<i_vdmean>(a)[i_velz] += W*vd.getProp<i_velocity>(b)[2];

                //gradient of particle property
                Point<3,double> dv = va - vb;
                Point<3,double> dP = vd.getProp<i_pressure>(a) - vd.getProp<i_pressure>(b);
                Point<3,double> dT = vd.getProp<i_temperature>(a)-vd.getProp<i_temperature>(b);
                Point<3,double> drho =  vd.getProp<i_rho>(a)*va-vd.getProp<i_rho>(b)*vb; //0.1*(va - vb);
                double visc = viscous(dr, r2, dv,vd.getProp<i_rho>(a), vd.getProp<i_rho>(b), 1, 0);

                for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
                {
                    // grad of particle property at a particle position
                    vd.template getProp<i_dvdmean>(a)[i][i_momentum] += (drho.get(i)*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += (dP.get(i)*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += (dT.get(i)*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_visc] += (visc*DW.get(i));
                }
                ingh++;
            }
            ++Np;
        }
        
        //divide by total kernel weight to get mean
        vd.template getProp<i_vdmean>(a)[i_rho] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_temperature] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_pressure] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_energy] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_velx] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_vely] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_velz] /= tot_W;
        
        //cout << "particle = " << ip << " neigh= "<< ingh  << endl;
        iavg += ingh;
        ++part;
    }

    cout << " Avg number of neighbours=" << iavg/ip << endl;

    //tsim.stop();
    //std::cout << "Time: " << tsim.getwct() << std::endl;
}

#endif // _neighbors_h