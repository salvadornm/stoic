//input:
//output:

#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, particleset  & vdmean, gradientset  & dvdmeanx, gradientset  & dvdmeany, gradientset  & dvdmeanz, int p, double dt, double l, turbulence turb)
{
    // stoic
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    double time_turb = 0.1;
    double eps_turb = 0.1;
    double dWien = distribution(generator); // eps_turb*distribution(generator)*sqrt(dt);
    
    //initial particle properties
    double Tp = vd.template getProp<i_temperature>(p);
    double Pp = vd.template getProp<i_pressure>(p);
    double rho_p = vd.template getProp<i_rho>(p);
    double energy_p = vd.template getProp<i_energy>(p);
    vector <double> up {vd.template getProp<i_velocity>(p)[0],vd.template getProp<i_velocity>(p)[1],vd.template getProp<i_velocity>(p)[2]};
    double vel_p = sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
    double m_p, mass_p, edensity_p;

    //mean particle properties
    double Tmean = vdmean.template getProp<i_temperature>(p);
    double Pmean = vdmean.template getProp<i_pressure>(p);
    double rho_mean = vdmean.template getProp<i_rho>(p);
    double energy_mean = vdmean.template getProp<i_energy>(p);
    vector <double> u_mean {vdmean.template getProp<i_velocity>(p)[0],vdmean.template getProp<i_velocity>(p)[1],vdmean.template getProp<i_velocity>(p)[2]};
    double vel_mean = sqrt(u_mean[0]*u_mean[0]+u_mean[1]*u_mean[1]+u_mean[2]*u_mean[2]);
    
    //gradient particle properties
    double P_grad, rho_grad, energy_grad, m_grad, T_grad;

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
    auto dvdmean = dvdmeanx;

    for (size_t i = 0; i < 3 ; i++) //loop through x,y,z directions
    {
        switch(i){
            case 0: dvdmean = dvdmeanx; break;
            case 1: dvdmean = dvdmeany; break;
            case 2: dvdmean = dvdmeanz; break;
        }

        //gradient particle properties
        T_grad = dvdmean.template getProp<i_temperature>(p);  //not used
        P_grad = dvdmean.template getProp<i_pressure>(p); 
        rho_grad = dvdmean.template getProp<i_rho>(p);
        energy_grad = dvdmean.template getProp<i_energy>(p);
        m_grad = dvdmean.template getProp<i_momentum>(p);
        // (0) check continuity

        // (1) update density
        drho = -(rho_grad)*dt;
        rho_new += drho;
        cout << "drho = " << drho << " rho_grad, rho_mean = " << rho_grad << ", " << rho_mean << endl;
 
        // (2) update momentum
        // drift term
        du = (u_mean[i] -  up[i]);
    
        //move to a class to manipulate outside
        //turbulent kinetic energy is half the sum of the variances of the velocities
        turb.k_sgs = .5 * ((up[i] - u_mean[i])*(up[i] - u_mean[i]));
        turb.Eps_sgs = turb.k_sgs/l;

        //time scales
        double Cu = 2.1; //Kolmogorov constant
        double C0 = 1;
        tau_sgs = turb.k_sgs/turb.Eps_sgs;
        tau_mol = (l*l)/(rho_p*up[i]);    //l = kernel width (H)
        tau_eq = 1/((1/tau_mol)+(Cu/tau_sgs));

        //solve momentum equation
        // dm = -dP/dx|dt + Audt + BdW
        Au_turb = ((rho_p*Cu)/tau_sgs)*du;
        Au_mol = (rho_p/tau_mol)*du;
        Au_p = Au_mol + Au_turb;
        Au_p_alt = (rho_p/(tau_eq+dt))*du;   //confirmed Au_p and Au_p_alt (wo +dt) are equivalent
        B = C0*sqrt(turb.Eps_sgs);   //turbulent diffusion

        //momentum equation
        m_p = rho_p * up[i];
        m_new = m_p - (P_grad * dt) + (Au_p_alt * dt) + (B * dWien);

        cout << "New momentum = " << m_new << " old momentum = "<< m_p  << endl;

        // (3) update velocity
        u_new = m_new / rho_p;
        vd.template getProp<i_velocity>(p)[i] = u_new;

        cout << "New u velocity = " << u_new << endl;

        // (4) solve energy density
        // (5) find specific enthalpy
        h_p = up[i] - (Pp / rho_p);
        h_mean = u_mean[i] - (Pmean / rho_mean);
        dh = h_mean - h_p;

        D_therm = (1);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

        tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));
        Ae_p = (rho_p/(tau_eq_energy + dt))*dh;

        //energy density equation
        edensity_p = rho_p * energy_p;
        //edensity_new += edensity_p - (P_grad * u_grad[i] * dt) + (Ae_p * dt) + (B * dWien);

        //energy_new = edensity_new / rho_p;
    }

    //check continuity again?

    //UPDATE OVERALL PROPERTIES
    vd.template getProp<i_rho>(p) = drho;
    vd.template getProp<i_energy>(p) = edensity_p / rho_p;// edensity_new / rho_p;
    cout << "New specific energy = " << vd.template getProp<i_energy>(p) << " old energy = "<< energy_p  << endl;
    cout << "New density = " << vd.template getProp<i_rho>(p) << " old density = "<< rho_p  << endl; 

    cout << "--------------------" << endl; 
}
#endif // _updateProps_h