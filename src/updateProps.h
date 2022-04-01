//input:
//output:

#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, particleset  & vdmean, particleset  & dvdmean, int p, double dt, double l)
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
    double up = vd.template getProp<i_velocity>(p)[0];
    double vp = vd.template getProp<i_velocity>(p)[1];
    double wp = vd.template getProp<i_velocity>(p)[2];
    double vel_p = sqrt(up*up+vp*vp+wp*wp);
    double m_p, mass_p, edensity_p;

    //mean particle properties
    double Tmean = vdmean.template getProp<i_temperature>(p);
    double Pmean = vdmean.template getProp<i_pressure>(p);
    double rho_mean = vdmean.template getProp<i_rho>(p);
    double energy_mean = vdmean.template getProp<i_energy>(p);
    double u2 = vdmean.template getProp<i_velocity>(p)[0];
    double v2 = vdmean.template getProp<i_velocity>(p)[1];
    double w2 = vdmean.template getProp<i_velocity>(p)[2];
    double vel_mean = sqrt(u2*u2+v2*v2+w2*w2);

    //gradient particle properties
    double Tgrad = dvdmean.template getProp<i_temperature>(p);  //not used
    double Pgrad = dvdmean.template getProp<i_pressure>(p); 
    double rho_grad = dvdmean.template getProp<i_rho>(p);
    double energy_grad = dvdmean.template getProp<i_energy>(p);
    double u3 = dvdmean.template getProp<i_velocity>(p)[0];
    double v3 = dvdmean.template getProp<i_velocity>(p)[1];
    double w3 = dvdmean.template getProp<i_velocity>(p)[2];
    double vel_grad = sqrt(u3*u3+v3*v3+w3*w3);

    //initialize deltas
    double drho, de, dm, dYk;
    double du, dh; //drift term
    //initialize new properties <-- dont need this bc should be += ?
    double rho_new, edensity_new, energy_new, Yk_new, T_new, P_new;
    double m_new, vel_new;
    //intialize timescales
    double tau_sgs, tau_mol, tau_eq, tau_eq_energy;
    // initialize coefficients
    double Au_turb, Au_mol, Au_p, Au_p_alt, B;
    double Ae_p, h_p, h_mean, D_therm;

    //check continuity to get density: 
    //drho/dt = -dmi/dxi -> drho = -dmi/dxi|dt; drho = -(sph gradient d(rho*u))/dx)*dt;
    drho = -(rho_grad * vel_grad)*dt;
    vd.template getProp<i_rho>(p) += drho;
    cout << "drho = " << drho << " rho_grad, rho_mean = " << rho_grad << ", " << rho_mean << endl;
    cout << "New density = " << vd.template getProp<i_rho>(p) << " old density = "<< rho_p  << endl;
        
    // drift term
    du = (vel_mean -  vel_p);
    cout << "vel_grad = "<< vel_grad  << " vel_mean = " << vel_mean << endl;
      
    //turbulent kinetic energy is half the sum of the variances of the velocities
    double k_sgs = .5 * ((vel_p - vel_mean)*(vel_p - vel_mean));
    double Eps_sgs = k_sgs/l;

    //time scales
    double Cu = 2.1; //Kolmogorov constant <-- check this
    double C0 = 1;
    tau_sgs = k_sgs/Eps_sgs;
    tau_mol = (l*l)/(rho_p*vel_p);    //l = kernel width (H)
    tau_eq = 1/((1/tau_mol)+(Cu/tau_sgs));

    //solve momentum equation
    // dm = -dP/dx|dt + Audt + BdW
    Au_turb = ((rho_p*Cu)/tau_sgs)*du;
    Au_mol = (rho_p/tau_mol)*du;
    Au_p = Au_mol + Au_turb;
    Au_p_alt = (rho_p/(tau_eq+dt))*du;   //Au_p and Au_p_alt (wo +dt) are equivalent
    B = C0*sqrt(Eps_sgs);   //turbulent diffusion

    //cout << "Au_p = " << Au_p << " Au_p_alt= "<< Au_p_alt  << endl;

    //momentum equation
    m_p = rho_p * vel_p;
    m_new = m_p - (Pgrad * dt) + (Au_p_alt * dt) + (B * dWien);

    cout << "New momentum = " << m_new << " old momentum = "<< m_p  << endl;

    //update velocity
    vel_new = m_new / vd.template getProp<i_rho>(p);

    cout << "New velocity = " << vel_new << " old velocity = "<< vel_p  << endl;

    //solve energy density
    //find specific enthalpy
    h_p = vel_p - (Pp / rho_p);
    h_mean = vel_mean - (Pmean / rho_mean);
    dh = h_mean - h_p;

    D_therm = (1);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

    tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));
    Ae_p = (rho_p/(tau_eq_energy + dt))*dh;

    //energy density equation
    edensity_p = rho_p * energy_p;
    edensity_new = edensity_p - (Pgrad * vel_grad * dt) + (Ae_p * dt) + (B * dWien);

    vd.template getProp<i_energy>(p) = edensity_new / vd.template getProp<i_rho>(p);

    cout << "New specific energy = " << vd.template getProp<i_energy>(p) << " old energy = "<< energy_p  << endl;




    //check continuity?
    
    //update velocity transport term
    //du*dt/(time_turb + dt) + dW;
    //vd.template getProp<i_velocity>(p)[0] += (vdmean.template getProp<i_velocity>(p)[0])*dt;
    //vd.template getProp<i_velocity>(p)[1] += (vdmean.template getProp<i_velocity>(p)[1])*dt;
    //vd.template getProp<i_velocity>(p)[2] += (vdmean.template getProp<i_velocity>(p)[2])*dt;

    cout << "--------------------" << endl;
    
}
#endif // _updateProps_h