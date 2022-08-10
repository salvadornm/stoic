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
    double rho_new, vel_new; 

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    turb.k_sgs = 0.0;
    for (size_t i = 0; i < 3 ; i++) turb.k_sgs += .5*(u_p[i] - u_mean[i])*(u_p[i] - u_mean[i]);
    turb.Eps_sgs = turb.k_sgs/l;
    
    //time scales
    double Cu = 2.1; //Kolmogorov constant
    double C0 = 1;
    double k = .025;    //[W/m K] thermal conductivity
    D_therm = k/(rho_p*cp_global);    //placeholder for thermal diffusivity (dependent on equivalence ratio)
    tau_sgs = turb.k_sgs/turb.Eps_sgs;
    for (size_t i = 0; i < 3 ; i++) tau_mol += (l*l)/(rho_p*u_p[i]);    //l = kernel width (H)
    tau_eq = 1/((1/tau_mol)+(Cu/tau_sgs));    
    tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));

    // (0) check continuity
    drho = - (mom_grad[0] + mom_grad[1] + mom_grad[2])*dt;

    // (1) update density ---------------------------------
    rho_new = rho_p + drho;
    vd.template getProp<i_rho>(p) =  rho_new;

    // (2) update momentum ---------------------
    Au_turb = ((rho_p*Cu)/tau_sgs);
    Au_mol = (rho_p/tau_mol);
    Au_p = Au_mol + Au_turb;
    Au_p_alt = (rho_p/(tau_eq+dt));   //confirmed Au_p and Au_p_alt (wo +dt) are equivalent
    B = C0*sqrt(turb.Eps_sgs);        //turbulent diffusion

    //SNM
    Au_p = 0.1;
    turb.k_sgs = 0;

    for (size_t i = 0; i < 3 ; i++) 
    {    
        du = (u_mean[i] -  u_p[i]);
        //solve momentum 
        mom_p[i]  = rho_p * u_p[i];
        mom_p[i] +=  P_grad[i] * dt + Au_p * du * dt + B * dWien[i] * sqrt(dt);
        
        vel_new =  mom_p[i] / rho_new;
        // clipping 
        vel_new = min(vel_new,50.0);
        vel_new = max(vel_new,-50.0);

        // (3) update velocity ---------------------
        vd.template getProp<i_velocity>(p)[i] = vel_new;

        turb.k_sgs += .5*(vel_new*vel_new);
    }

    // (4) find specific enthalpy ---------------------
    h_p =  energy_p      + (Pp / rho_p);
    h_mean = energy_mean + (Pmean / rho_mean);
    dh = h_mean - h_p;

    // (5) solve energy density ---------------------
    Ae_p = rho_p/(tau_eq_energy + dt);
    Ae_p = 0.1; //SNM

    edensity_p = rho_p * energy_p;
    dvisc = (visc_grad[0] + visc_grad[1] + visc_grad[2])*dt;

    edensity_new = edensity_p - (dvisc) + (Ae_p * dt * dh);    //check viscosity term
    energy_new = edensity_new/rho_new;    //specific total energy

    //SNM: U = Etotal - Ekin
    vd.template getProp<i_energy>(p)= energy_new - turb.k_sgs; //internal energy
    
    // (6) check continuity again - update temperature and pressure
    updateThermalProperties1(vd, p);
    
    //output_energy_props(vd, p, dh, dvisc, edensity_p, edensity_new, energy_new);

    //UPDATE OVERALL PROPERTIES

}
#endif // _updateProps_h