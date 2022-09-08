#ifndef _neighbors_h
#define _neighbors_h

#include "global.h"
#include "kernel.h"
#include "calculations.h"
#include"test.h"

template<typename CellList> void find_neighbors(particleset  & vd, CellList & NN, Cfd sim){
    int n,ingh,ip;
    float iavg=0;
    double avg_error,avg_grad;
    double *dW;
     
    // ------------------------------------------------------------------ 
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    //stateOfNeighbors(vd, NN);     //output

    // For each particle ...
    ip=0;  avg_error = 0; avg_grad=0;
  
    while (part.isNext())
    {
        auto a = part.get();
        ip++;

        // position and velocity particle a
        Point<3,double> xa = vd.getPos(a);
        Point<3,double> va = vd.getProp<i_velocity>(a);

       // reset counters        
        for (int nv = 0; nv < NVARSOLVE; nv++)
        { 
            vd.template getProp<i_vdmean>(a)[nv] = 0.0; 
            for (int i = 0; i < NDIM ; i++)
                vd.template getProp<i_dvdmean>(a)[i][nv]    = 0.0;  
        }

        // Get an iterator of all the particles neighborhood of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        double tot_W = 0.0; double tot_dW = 0.0;
        
        // For each neighborhood particle
        ingh =0;
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            // position and velocity particle b
            Point<3,double> xb = vd.getPos(b);          
            Point<3,double> vb = vd.getProp<i_velocity>(b);
                
            // Get the distance between a and b
            Point<3,double> dr = xb - xa;
            double r2 = norm2(dr);     
            double r = abs(sqrt(r2));
                        
            //if the particles interact...
            if (r < 2.0*sim.H) 
            {    
                double W = Wab(r, sim.H);               //kernel 
                dW = DWab(r,sim.H,dr);    // gradient kernel       
                tot_W += W;

                // q<- getProp<*>(b)
                
                vd.template getProp<i_vdmean>(a)[i_rho]         += W*vd.getProp<i_rho>(b);
                vd.template getProp<i_vdmean>(a)[i_temperature] += W*vd.getProp<i_temperature>(b);
                vd.template getProp<i_vdmean>(a)[i_pressure]    += W*vd.getProp<i_pressure>(b);
                vd.template getProp<i_vdmean>(a)[i_energy]      += W*vd.getProp<i_energy>(b);
                vd.template getProp<i_vdmean>(a)[i_velx]        += W*vd.getProp<i_velocity>(b)[0];
                vd.template getProp<i_vdmean>(a)[i_vely]        += W*vd.getProp<i_velocity>(b)[1];
                vd.template getProp<i_vdmean>(a)[i_velz]        += W*vd.getProp<i_velocity>(b)[2];

                //gradient of particle property
                Point<3,double> drho    = vd.getProp<i_rho>(b)*vb      -vd.getProp<i_rho>(a)*va;
                double dP               = vd.getProp<i_pressure>(b)    -vd.getProp<i_pressure>(a);
                double dT               = vd.getProp<i_temperature>(b) -vd.getProp<i_temperature>(a);

                
                for (int i = 0; i < NDIM ; i++) 
                {                       
                    vd.template getProp<i_dvdmean>(a)[i][i_momentum]    += W*drho.get(i)*dW[i];
                    vd.template getProp<i_dvdmean>(a)[i][i_pressure]    += W*dP*dW[i];
                    vd.template getProp<i_dvdmean>(a)[i][i_temperature] += W*dT*dW[i];           
                }
                ingh++;

                tot_dW += W ;
                if (a.getKey() == b)  
                  { tot_dW -= W;} 
 
            // cout << "tot_W= " <<  tot_W << " 1/W=" << 1/W << endl;
            // cout << "r= " <<  r << " h= " << sim.H << " r/h= " << r/sim.H << endl;           
            // cout << "Pb= " << Pb  << "Pa= " << Pa << endl;
            // cout << "xb =" << vd.getPos(b)[0]  << "xa= " << vd.getPos(a)[0] << endl;
            // cout << "dP/dx (est) " << (Pb - Pa)/(vd.getPos(b)[0] -vd.getPos(a)[0] )<< endl;
            // cout << "dP*dW[0] " <<  dP*dW[0]<< endl;
            // cout << "dW" << dW[0] << " " << dW[1] << " " << dW[2] << endl;
            //             exit(0);//temp snm

                
            }
            ++Np;
        }
        //divide by total kernel weight to get mean
        for (int nv = 0; nv < NVARSOLVE; nv++)
        {
            vd.template getProp<i_vdmean>(a)[nv]         /= tot_W;
        }
        
        //to find the gradients need to divide by total kernel weight.....
        for (size_t i = 0; i < 3 ; i++)                 
        {
            vd.template getProp<i_dvdmean>(a)[i][i_momentum]    /= tot_dW;
            vd.template getProp<i_dvdmean>(a)[i][i_pressure]    /= tot_dW;
            vd.template getProp<i_dvdmean>(a)[i][i_temperature] /= tot_dW;
        }

        
        iavg += ingh;
        ++part;

       // temp snm ******
        cout << "dp/dx= " << vd.template getProp<i_dvdmean>(a)[0][i_pressure] << endl;
        cout << "totW= " << tot_W;
        avg_error += abs(vd.template getProp<i_dvdmean>(a)[0][i_pressure] - 101300)/101300;

        avg_grad += vd.template getProp<i_dvdmean>(a)[0][i_pressure];

    }

    cout << " Avg number of neighbours=" << iavg/ip << endl;
    cout << " Avg error=" << avg_error*100/ip << " %" << endl;
    cout << " Avg grad =" << avg_grad/ip << endl;
    



    //tsim.stop();
    //std::cout << "Time: " << tsim.getwct() << std::endl;
}

#endif // _neighbors_h