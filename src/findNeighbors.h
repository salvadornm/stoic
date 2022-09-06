#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"
#include"test.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, Cfd sim, engine eng){
    int n,ingh,ip;
    float iavg=0;
    
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    //stateOfNeighbors(vd, NN);     //output

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
        double tot_factor = 0;
        
        // For each neighborhood particle
        ingh =0;
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            Point<3,double> xb = vd.getPos(b);  // position xp of the particle
            Point<3,double> vb = vd.getProp<i_velocity>(b);
            
            //if (a.getKey() == b)    {++Np; continue;};// if (a == b) skip this particle
                
            // Get the distance between a and b
            Point<3,double> dr = xa - xb;
            double r2 = norm2(dr);  //norm2 = (sum of the squares)
            double r = sqrt(r2);
                        
            //if the particles interact...
            if (r < 2.0*sim.H) {
                Point<3,double> DW {0.0,0.0,0.0};
                double factor = DWab(dr,DW,r,sim.H); // gradient kernel //
                double W = Wab(r, sim.H); //kernel
                tot_W += W;
                tot_factor += factor;          
                
                vd.template getProp<i_vdmean>(a)[i_rho] += W*vd.getProp<i_rho>(b);
                vd.template getProp<i_vdmean>(a)[i_temperature] += W*vd.getProp<i_temperature>(b);
                vd.template getProp<i_vdmean>(a)[i_pressure]    += W*vd.getProp<i_pressure>(b);
                vd.template getProp<i_vdmean>(a)[i_energy]      += W*vd.getProp<i_energy>(b);
                vd.template getProp<i_vdmean>(a)[i_velx] += W*vd.getProp<i_velocity>(b)[0];
                vd.template getProp<i_vdmean>(a)[i_vely] += W*vd.getProp<i_velocity>(b)[1];
                vd.template getProp<i_vdmean>(a)[i_velz] += W*vd.getProp<i_velocity>(b)[2];

                //gradient of particle property
                Point<3,double> dv = va - vb;
                double dP = vd.getProp<i_pressure>(a) - vd.getProp<i_pressure>(b);
                double dT = vd.getProp<i_temperature>(a)-vd.getProp<i_temperature>(b);
                Point<3,double> drho =  vd.getProp<i_rho>(a)*va-vd.getProp<i_rho>(b)*vb; //0.1*(va - vb);
                //drho =  vd.getProp<i_rho>(a)-vd.getProp<i_rho>(b); //0.1*(va - vb);
                double visc = viscous(dr, r2, dv,vd.getProp<i_rho>(a), vd.getProp<i_rho>(b), 1, 0);
                double dviscP =  visc*(vd.getProp<i_pressure>(a) - vd.getProp<i_pressure>(b));
                //Point<3,double> dviscP =  va*vd.getProp<i_pressure>(a) - vb*vd.getProp<i_pressure>(b);

                for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
                {
                    // grad of particle property at a particle position
                    /*cout << "i= " << i << "------------" << endl;
                    cout << "va: " << va.get(i) << endl;
                    cout << "vb: " << vb.get(i) << endl;
                    cout << "rhoa: " << vd.getProp<i_rho>(a) << endl;
                    cout << "rhob: " << vd.getProp<i_rho>(b) << endl;
                    cout << "drho: " << drho.get(i) << endl;
                    double drho_double =  vd.getProp<i_rho>(a)*va.get(i)-vd.getProp<i_rho>(b)*vb.get(i); 
                    cout << "drho: " << drho_double << endl;
                    */
                    //vd.template getProp<i_dvdmean>(a)[i][i_momentum] += (vd.getProp<i_rho>(b)/eng.volumeC)*dv.get(i)*DW.get(i);
                    
                    //cout << "DW: " << DW.get(i) << endl;

                    vd.template getProp<i_dvdmean>(a)[i][i_momentum] += (drho.get(i)*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += (dP*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += (dT*DW.get(i));
                    vd.template getProp<i_dvdmean>(a)[i][i_visc] += (dviscP*DW.get(i));
                    /*
                    vd.template getProp<i_dvdmean>(a)[i][i_momentum] += (drho.get(i)*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += (dP*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += (dT*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_visc] += (dviscP*factor);
                    */
                }
                ingh++;
            }
            ++Np;
        }

        //cout << "ingh: " << ingh << endl;
        //W ~ G/V
        //to find the gradients need to divide by total kernel weight
        for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
        {
            vd.template getProp<i_dvdmean>(a)[i][i_momentum] /= tot_W; //ingh;
            vd.template getProp<i_dvdmean>(a)[i][i_pressure]    /= tot_W; //ingh;
            vd.template getProp<i_dvdmean>(a)[i][i_temperature] /= tot_W; //ingh;
            vd.template getProp<i_dvdmean>(a)[i][i_visc] /= tot_W; //ingh;
        }

        //divide by total kernel weight to get mean
        tot_W = max(tot_W,1.0);
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