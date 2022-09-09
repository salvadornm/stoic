#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"
#include"test.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, Cfd sim){
    int n,ingh,ip;
    float iavg=0;
    double avg_error, avg_grad, avg_mean;
    double *dW; double W;
    
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    
    //stateOfNeighbors(vd, NN);     //output for testing
    timer tsim;
    tsim.start();

    // For each particle ...
    ip = 0;   
    avg_error = 0;
    avg_grad = 0;
    avg_mean = 0;
  
    while (part.isNext())
    {
        auto a = part.get();
        ip++;

        // Position and Velocity of Particle a
        Point<3,double> xa = vd.getPos(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

       // reset counters        
        for (size_t j = 0; j < 7.0 ; j++)
        { vd.template getProp<i_vdmean>(a)[j] = 0.0; }
        for (size_t j = 0; j < NDIM ; j++)
        { 
            vd.template getProp<i_dvdmean>(a)[j][i_momentum]    = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_rho]         = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_energy]      = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_pressure]    = 0.0;
            vd.template getProp<i_dvdmean>(a)[j][i_temperature] = 0.0;
        }

        // Get an iterator of all the particles neighborhood of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        double tot_W = 0;
        double tot_dW = 0;
        
        // For each neighborhood particle
        ingh =0;
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            // Position and velocity of particle b
            Point<3,double> xb = vd.getPos(b);
            Point<3,double> vb = vd.getProp<i_velocity>(b);
                
            // Get the distance between a and b
            Point<3,double> dr = xb - xa;
            double r2   = norm2(dr);
            double r    = abs(sqrt(r2));
                        
            //if the particles interact...
            if (r < 2.0*sim.H) {
                
                W = Wab(r, sim.H);          // kernel
                dW = DWab(dr,r,sim.H);      // gradient kernel
                
                if (a.getKey() == b){tot_dW -= W;}         //if particle is itself 
                tot_W += W;
                tot_dW += W;
                
                vd.template getProp<i_vdmean>(a)[i_rho]         += W*vd.getProp<i_rho>(b);
                vd.template getProp<i_vdmean>(a)[i_temperature] += W*vd.getProp<i_temperature>(b);
                vd.template getProp<i_vdmean>(a)[i_pressure]    += W*vd.getProp<i_pressure>(b);
                vd.template getProp<i_vdmean>(a)[i_energy]      += W*vd.getProp<i_energy>(b);
                vd.template getProp<i_vdmean>(a)[i_velx]        += W*vd.getProp<i_velocity>(b)[0];
                vd.template getProp<i_vdmean>(a)[i_vely]        += W*vd.getProp<i_velocity>(b)[1];
                vd.template getProp<i_vdmean>(a)[i_velz]        += W*vd.getProp<i_velocity>(b)[2];

                //gradient of particle property
                Point<3,double> dv      = vb - va;
                Point<3,double> drho    = vd.getProp<i_rho>(b)*vb-vd.getProp<i_rho>(a)*va;
                double dP               = vd.getProp<i_pressure>(b) - vd.getProp<i_pressure>(a);
                double dT               = vd.getProp<i_temperature>(b)-vd.getProp<i_temperature>(a);
                //double visc             = viscous(dr, r2, dv,vd.getProp<i_rho>(b), vd.getProp<i_rho>(a), 1, 0);
                //double dviscP           = visc*(vd.getProp<i_pressure>(b) - vd.getProp<i_pressure>(a));

                for (size_t i = 0; i < NDIM ; i++) //loop through x,y,z directions
                {   
                    double dx = 1/(dr.get(i) + 1e-8);
                    //cout << "dx: " << dx << endl;
                    //cout << "W: " << W << endl;
                    
                    //Using Approach 1
                    vd.template getProp<i_dvdmean>(a)[i][i_momentum]    += (W*drho.get(i)*dx);
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += (W*dP*dx);
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += (W*dT*dx);
                    //vd.template getProp<i_dvdmean>(a)[i][i_visc]        += (W*dviscP*DW.get(i));
                                        
                    /*
                    vd.template getProp<i_dvdmean>(a)[i][i_momentum]    += (drho.get(i)*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += (dP*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += (dT*factor);
                    vd.template getProp<i_dvdmean>(a)[i][i_visc]        += (dviscP*factor);
                    */
                }
                ingh++;
            }
            ++Np;
        }
        //cout << "ingh: " << ingh << endl;
        //cout << "tot_W: " << tot_W << endl;

        //divide by total kernel weight to get mean
        vd.template getProp<i_vdmean>(a)[i_rho]         /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_temperature] /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_pressure]    /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_energy]      /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_velx]        /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_vely]        /= tot_W;
        vd.template getProp<i_vdmean>(a)[i_velz]        /= tot_W;

        // divide by total kernel weight to find the gradients.....
        for (size_t i = 0; i < NDIM ; i++)
        {
            vd.template getProp<i_dvdmean>(a)[i][i_momentum]    /= tot_dW;
            vd.template getProp<i_dvdmean>(a)[i][i_pressure]    /= tot_dW;
            vd.template getProp<i_dvdmean>(a)[i][i_temperature] /= tot_dW;
            vd.template getProp<i_dvdmean>(a)[i][i_visc]        /= tot_dW;
        }
        
        iavg += ingh;
        ++part;

        // temp snm ******
        //cout << "totW= " << tot_W;
        //cout << "totdW= " << tot_dW;
        //cout << ", dp/dx= " << vd.template getProp<i_dvdmean>(a)[0][i_pressure] << endl;
        
        avg_error += abs(vd.template getProp<i_dvdmean>(a)[0][i_pressure] - 101300)/101300;
        avg_grad += vd.template getProp<i_dvdmean>(a)[0][i_pressure];
        avg_mean += vd.template getProp<i_vdmean>(a)[i_pressure];
    }

    cout << " Avg number of neighbours=" << iavg/ip << endl;
    cout << " Avg error=" << avg_error*100/ip << " %" << endl;
    cout << " Avg grad =" << avg_grad/ip << endl;
    cout << " Avg mean =" << avg_mean/ip << endl;

    tsim.stop();
    std::cout << "Time: " << tsim.getwct() << std::endl;
}

#endif // _neighbors_h