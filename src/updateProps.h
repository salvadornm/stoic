//input:
//output:

#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, int p, double dt, double l, turbulence turb, Cfd sim )
{
    // stoic    
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(0.0,sim.lx);
    double time_turb = 0.1;
    double eps_turb = 0.1;
    vector <double> dWien {distribution(generator)*sqrt(dt),distribution(generator)*sqrt(dt),distribution(generator)*sqrt(dt)}; // eps_turb*distribution(generator)*sqrt(dt);
    
    //initial particle properties (vd(p))
    double Tp = vd.template getProp<i_temperature>(p);
    double Pp = vd.template getProp<i_pressure>(p);
    double rho_p = vd.template getProp<i_rho>(p);
    vector <double> mom_p {0.0,0.0,0.0};
    double energy_p = vd.template getProp<i_energy>(p);
    vector <double> u_p {vd.template getProp<i_velocity>(p)[0],vd.template getProp<i_velocity>(p)[1],vd.template getProp<i_velocity>(p)[2]};
    double mass_p, edensity_p;


    //mean particle properties (vdmean(p))
    double Tmean = vd.template getProp<i_vdmean>(p)[i_temperature];
    double Pmean = vd.template getProp<i_vdmean>(p)[i_pressure];
    double rho_mean = vd.template getProp<i_vdmean>(p)[i_rho];
    double energy_mean = vd.template getProp<i_vdmean>(p)[i_temperature];
    vector <double> u_mean {vd.template getProp<i_vdmean>(p)[i_velx],vd.template getProp<i_vdmean>(p)[i_vely],vd.template getProp<i_vdmean>(p)[i_velz]};
        
    //gradient particle properties (dvdmean(p))
    double energy_grad;
    vector <double> T_grad {vd.template getProp<i_dvdmean>(p)[0][i_temperature],vd.template getProp<i_dvdmean>(p)[1][i_temperature],vd.template getProp<i_dvdmean>(p)[2][i_temperature]}; 
    vector <double> P_grad {vd.template getProp<i_dvdmean>(p)[0][i_pressure],vd.template getProp<i_dvdmean>(p)[1][i_pressure],vd.template getProp<i_dvdmean>(p)[2][i_pressure]}; 
    vector <double> mom_grad {vd.template getProp<i_dvdmean>(p)[0][i_momentum],vd.template getProp<i_dvdmean>(p)[1][i_momentum],vd.template getProp<i_dvdmean>(p)[2][i_momentum]};
    vector <double> visc_grad {vd.template getProp<i_dvdmean>(p)[0][i_visc],vd.template getProp<i_dvdmean>(p)[1][i_visc],vd.template getProp<i_dvdmean>(p)[2][i_visc]};

    //initialize deltas
    double drho, de, dm, dYk, dvisc;
    double du, dh; //drift term
    
    //intialize timescales
    double tau_sgs, tau_mol, tau_eq, tau_eq_energy;

    // initialize coefficients
    double Au_turb, Au_mol, Au_p, Au_p_alt, B;
    double Ae_p, h_p, h_mean, D_therm;

    // initialize new particles
    double edensity_new, energy_new, u_new, m_new;
    double rho_new; 

    // aux
    double vel_i,kinetic_energy,freq_l,freq_sgs,freq_mol,freq_eq;
    double viscp= 1.3e-5;
    double k = .025;    //[W/m K] thermal conductivity

    double Cu = 2.1; //Kolmogorov constant
    double C0 = 1;

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    // calculate  tubulence kinetic energy
    turb.k_sgs = 1.e-8;
    for (size_t i = 0; i < 3 ; i++) turb.k_sgs += .5*(u_p[i] - u_mean[i])*(u_p[i] - u_mean[i]);
    turb.Eps_sgs = pow(turb.k_sgs,1.5)/l; //energy dissipation
    freq_sgs = turb.Eps_sgs/turb.k_sgs; 
    tau_sgs = 1.0/freq_sgs;
    
    
    // std::cout << "turb.k_sgs= " << turb.k_sgs <<  std::endl;
    // for (size_t i = 0; i < 3 ; i++) {
    // std::cout << "u_mean[i]= " << u_mean[i] <<  std::endl;
    // std::cout << "u_p[i]= " << u_p[i] <<  std::endl;
    // }
    // std::cout << "tau_sgs = " << tau_sgs <<  std::endl;


    // compute velocity partcile (absolute value)
    vel_i = 0.0;
    for (size_t i = 0; i < 3 ; i++) vel_i += u_p[i]*u_p[i]; 
    vel_i = sqrt(vel_i);
    freq_l = vel_i/l; // unit [1/s]
    // moleculart frequency and time scale
    freq_mol = viscp/(rho_p*l*l);
    tau_mol = 1.0/freq_mol;
    
    //std::cout << "tau_mol = " << tau_mol <<  std::endl;


    //time scales
    freq_eq = Cu*freq_sgs + freq_mol;
    tau_eq = dt + 1.0/freq_eq;    // dt added for stability
    tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));
    tau_eq_energy = tau_eq ;
    
    D_therm = k/(rho_p*cp_global);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

   // std::cout << "tau_eq = " << tau_eq <<  std::endl;


    // (0) check continuity
    drho = - (mom_grad[0] + mom_grad[1] + mom_grad[2])*dt;

    // (1) update density ---------------------------------
    rho_new = rho_p + drho;
    vd.template getProp<i_rho>(p) =  rho_new;

    // (2) update momentum ---------------------
    Au_p = (0.5 + 0.75*Cu)*rho_p/tau_eq;      // Au > 0 
    B = C0*sqrt(turb.Eps_sgs)* sqrt(dt);   //turbulent diffusion
    
    // std::cout << "Au = " << Au_p <<  std::endl;
    // std::cout << "B = " << B <<  std::endl;
    // std::cout << "Esgs= " << turb.Eps_sgs <<  std::endl;
 
     //Au_p = 0.1;
    //B = 0;

    for (size_t i = 0; i < 3 ; i++) 
    {    
        du = u_mean[i] -  u_p[i];
        //solve momentum 
        mom_p[i]  = rho_p * u_p[i];
        mom_p[i] +=  P_grad[i] * dt + Au_p * du * dt + B * dWien[i];
                
        // (3) update velocity ---------------------
        vel_i = mom_p[i] / rho_new;

        // clipping 
        vel_i = min(vel_i,10.0);
        vel_i = max(vel_i,-10.0);

        // snm
        // mom_p[i]  = rho_p * u_p[i] + 0.1 * du * dt;
        // vel_i  = mom_p[i] / rho_p;

        kinetic_energy += 0.5*vel_i*vel_i;
        vd.template getProp<i_velocity>(p)[i] = vel_i;
        // limit_velocity(vd, p, i);
        
    }

    // (4) find specific enthalpy ---------------------
    h_p =  energy_p      + (Pp / rho_p);
    h_mean = energy_mean + (Pmean / rho_mean);
    dh = h_mean - h_p;

    // (5) solve energy density ---------------------
    Ae_p = rho_p/tau_eq_energy;

    edensity_p = rho_p * energy_p;
    dvisc = (visc_grad[0] + visc_grad[1] + visc_grad[2])*dt; //CHECK

    edensity_new = edensity_p - (dvisc) + Ae_p * dh* dt;    //check viscosity term
    
    energy_new = edensity_new/rho_new;    //specific total energy

    // SNM: total energy not updated 
   // energy_new = energy_p;

    // SNM:  U = Etotal - Ekin
    vd.template getProp<i_energy>(p)= energy_new - kinetic_energy; 

    // cout << "dh: " << dh << " dviscP: " << dvisc << endl;
    // cout << "edensity_p: " << edensity_p << " edensity_new: " << edensity_new << endl;
    // cout << "energy: " << energy_new << " internal energy: " << vd.template getProp<i_energy>(p) << endl;
    
     //---------------------------

    // (6): update temperature and pressure
    updateThermalProperties1(vd, p);
    // cout << "New temp: " << vd.template getProp<i_temperature>(p);
    // cout << " New pressure: " << vd.template getProp<i_pressure>(p) << endl;



    //UPDATE OVERALL PROPERTIES

}
#endif // _updateProps_h