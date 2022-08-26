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
    double energy_mean = vd.template getProp<i_vdmean>(p)[i_energy];
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
    double Au_p, Ae_p, B;
    double h_p, h_mean, D_therm;

    // initialize new particles
    double edensity_new, energy_new, u_new, m_new;
    double rho_new, vel_new; 

    //intialize temp values
    double freq_sgs, freq_mol, freq_eq, freq_eq_energy;

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    //approximate viscosity (air,400K,1bar)
    double viscp = 2.59E-5;
    turb.k_sgs = 1e-8;

    //time scales
    double Cu = 2.1; //Kolmogorov constant
    turb.C0 = 1.0;
    turb.delta = 9.92e-4;
    double k = .025;    //[W/m K] thermal conductivity

    D_therm = k/(rho_p*cp_global);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

    tau_mol = 0.0; tau_sgs = 0; tau_eq = 0; tau_eq_energy = 0;
    
    
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    //turbulence time scales
    for (size_t i = 0; i < 3 ; i++) turb.k_sgs += 0.5*(u_p[i] - u_mean[i])*(u_p[i] - u_mean[i]);
    turb.Eps_sgs = pow(turb.k_sgs,1.5)/l; 
    
    freq_sgs = turb.Eps_sgs/turb.k_sgs; 
    tau_sgs = 1/freq_sgs;

    //molecular times scale
    freq_mol = viscp/(rho_p*l*l);
    tau_mol = 1.0/freq_mol;

    //cout << "tau_mol: " << tau_mol << " tau_sgs: " << tau_sgs << endl; 
    //cout << "cu: " << Cu << " freq_sgs: " << freq_sgs << " freq_mol: " << freq_mol << endl; 

    //equivalent time scale
    freq_eq = Cu*freq_sgs + freq_mol;
    tau_eq = dt + 1/freq_eq;

    freq_eq_energy = (D_therm/(l*l))+ (Cu*freq_sgs);    
    tau_eq_energy = dt + 1/freq_eq_energy;

    //cout << "tau_eq: " << tau_eq << " tau_eq_energy: " << tau_eq_energy << endl; 

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // (0) check continuity
    drho = - (mom_grad[0] + mom_grad[1] + mom_grad[2])*dt;

    // (1) update density ---------------------------------
    rho_new = rho_p + drho;
    vd.template getProp<i_rho>(p) =  rho_new;
    cout << "drho: " << drho << " rhoNew: " << rho_new << endl;

    //tau_eq = 10; tau_eq_energy = tau_eq;
   
    // (2) update momentum ---------------------
    Au_p = (0.5 + 0.75*Cu) * rho_p / tau_eq;
    B = sqrt(turb.C0*turb.Eps_sgs) * sqrt(dt);        //turbulent diffusion
    
    //cout << "Au_P solved: " << Au_p << " tau_eq: " << tau_eq << endl;

    //SNM
    double ke_p = 0;
    double ke = 0;
    double ke_mean = 0;

    for (size_t i = 0; i < 3 ; i++) 
    {    
        du = (u_mean[i] -  u_p[i]);
        //solve momentum 
        mom_p[i]  = rho_p * u_p[i];
        mom_p[i] +=  P_grad[i] * dt + Au_p * du * dt + B * dWien[i];
        
        vel_new =  mom_p[i] / rho_new;
        // clipping 
        vel_new = min(vel_new,50.0);
        vel_new = max(vel_new,-50.0);

        // (3) update velocity ---------------------
        vd.template getProp<i_velocity>(p)[i] = vel_new;

        ke_p += .5*(u_p[i]*u_p[i]);
        ke += .5*(vel_new*vel_new);
        ke_mean += .5*(u_mean[i]*u_mean[i]);
    }

    // (4) find specific enthalpy ---------------------
    h_p =  (energy_p) + (Pp / rho_p);
    h_mean = (energy_mean) + (Pmean / rho_mean);
    dh = h_mean - h_p; //dh = cp_global*(Tmean - Tp);

    // (5) solve energy density ---------------------
    Ae_p = (rho_p * dh)/(tau_eq_energy);
    
    edensity_p = rho_p * (energy_p);
    dvisc = (visc_grad[0] + visc_grad[1] + visc_grad[2])*dt;    //CHECK ... should this be velocity * P instead of visc*P?
    //dvisc = 0;
    //Ae_p = 0;

    edensity_new = edensity_p - (dvisc) + (Ae_p * dt);
    energy_new = edensity_new/rho_new;    //specific total energy

    //SNM: U = Etotal - Ekin
    vd.template getProp<i_energy>(p)= energy_new - ke; //internal energy
    
    // (6) check continuity again - update temperature and pressure
    updateThermalProperties1(vd, p);
    //output_energy_props(vd, p, dh, dvisc, edensity_p, edensity_new, energy_new);

}
#endif // _updateProps_h