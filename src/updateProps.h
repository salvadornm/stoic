//input:
//output:

#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, int p, double dt, double l, turbulence turb)
{
    // stoic
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    double time_turb = 0.1;
    double eps_turb = 0.1;
    double dWien = distribution(generator); // eps_turb*distribution(generator)*sqrt(dt);
    
    //initial particle properties (vd(p))
    double Tp = vd.template getProp<i_temperature>(p);
    double Pp = vd.template getProp<i_pressure>(p);
    double rho_p = vd.template getProp<i_rho>(p);
    vector <double> mom_p {0.0,0.0,0.0};
    double energy_p = vd.template getProp<i_energy>(p);
    vector <double> u_p {vd.template getProp<i_velocity>(p)[0],vd.template getProp<i_velocity>(p)[1],vd.template getProp<i_velocity>(p)[2]};
    double vel_p = sqrt(u_p[0]*u_p[0]+u_p[1]*u_p[1]+u_p[2]*u_p[2]);
    double mass_p, edensity_p;


    //mean particle properties (vdmean(p))
    double Tmean = vd.template getProp<i_vdmean>(p)[i_temperature];
    double Pmean = vd.template getProp<i_vdmean>(p)[i_pressure];
    double rho_mean = vd.template getProp<i_vdmean>(p)[i_rho];
    double energy_mean = vd.template getProp<i_vdmean>(p)[i_temperature];
    vector <double> u_mean {vd.template getProp<i_vdmean>(p)[i_velocity],vd.template getProp<i_vdmean>(p)[i_vely],vd.template getProp<i_vdmean>(p)[i_velz]};
    //double vel_mean = sqrt(u_mean[0]*u_mean[0]+u_mean[1]*u_mean[1]+u_mean[2]*u_mean[2]);
    
    //gradient particle properties (dvdmean(p))
    double energy_grad, T_grad;
    vector <double> P_grad {vd.template getProp<i_dvdmean>(p)[0][i_pressure],vd.template getProp<i_dvdmean>(p)[1][i_pressure],vd.template getProp<i_dvdmean>(p)[2][i_pressure]}; 
    vector <double> mom_grad {vd.template getProp<i_dvdmean>(p)[0][i_rho],vd.template getProp<i_dvdmean>(p)[1][i_rho],vd.template getProp<i_dvdmean>(p)[2][i_rho]};

    //initialize deltas
    double drho, de, dm, dYk;
    double du, dh; //drift term
    
    //intialize timescales
    double tau_sgs, tau_mol, tau_eq, tau_eq_energy;

    // initialize coefficients
    double Au_turb, Au_mol, Au_p, Au_p_alt, B;
    double Ae_p, h_p, h_mean, D_therm;

    // initialize new particles
    double edensiy_new, energy_new, u_new, m_new;
    double rho_new; 

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    turb.k_sgs = 0.0;
    for (size_t i = 0; i < 3 ; i++) turb.k_sgs += (u_p[i] - u_mean[i])*(u_p[i] - u_mean[i]);
    turb.Eps_sgs = turb.k_sgs/l;
    
    //time scales
    double Cu = 2.1; //Kolmogorov constant
    double C0 = 1;
    tau_sgs = turb.k_sgs/turb.Eps_sgs;
    tau_mol = 0.001;//(l*l)/(rho_p*u_p[i]);    //l = kernel width (H)
    tau_eq = 1/((1/tau_mol)+(Cu/tau_sgs));

    // (0) check continuity

    // (1) update density ---------------------------------
    drho = - (mom_grad[0] + mom_grad[1] + mom_grad[2])*dt;
    rho_new = rho_p + drho;
    //std::cout << "drho = " << drho << " rho_new = " << rho_new << endl;
    vd.template getProp<i_rho>(p) =  rho_new;

    // (2) update momentum ---------------------
    Au_turb = ((rho_p*Cu)/tau_sgs);
    Au_mol = (rho_p/tau_mol);
    Au_p = Au_mol + Au_turb;
    Au_p_alt = (rho_p/(tau_eq+dt));   //confirmed Au_p and Au_p_alt (wo +dt) are equivalent
    B = C0*sqrt(turb.Eps_sgs);        //turbulent diffusion
        
    for (size_t i = 0; i < 3 ; i++) 
    {
        //test velocity! in each direction
        //du = (u_mean[i] -  u_p[i]); 
        //vd.template getProp<i_velocity>(p)[i] += (du*dt);
    }
    /*
        // drift term
        du = (u_mean[i] -  u_p[i]);
        //solve momentum equation
        mom_p[i]  = rho_p * u_p[i];
        mom_p[i] +=  P_grad[i] * dt + Au_p_alt * du * dt + B * dWien * sqrt(dt);
        
        //std::cout << "New momentum = " << mom_p[i] << endl;
        //std::cout << "rho_p = " << rho_p << endl;
        
        // (3) update velocity
        vd.template getProp<i_velocity>(p)[i] = mom_p[i] / rho_new;
        //std::cout << "New u velocity = " << vd.template getProp<i_velocity>(p)[i] << endl;
    }

    // (4) solve energy density
    // (5) find specific enthalpy
    h_p =  energy_p      + (Pp / rho_p);
    h_mean = energy_mean - (Pmean / rho_mean);
    dh = h_mean - h_p;

    D_therm = (1);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

    tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));
    Ae_p = (rho_p/(tau_eq_energy + dt))*dh;

    //energy density equation
    edensity_p = rho_p * energy_p;
    //  edensity_p += - (P_grad * u_grad[i] * dt) + (Ae_p * dt) + (B * dWien);

    // vd.template getProp<i_energy>(p)= edensity_p/rho_new;

    //check continuity again?

    //UPDATE OVERALL PROPERTIES
    //std::cout << "New specific energy = " << vd.template getProp<i_energy>(p) << " old energy = "<< energy_p  << endl;
    //std::cout << "New density = " << vd.template getProp<i_rho>(p) << " old density = "<< rho_p  << endl; 

    //std::cout << "--------------------" << endl; 
    */
}
#endif // _updateProps_h