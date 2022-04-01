//input:
//output:

#ifndef _updateProps_h
#define _updateProps_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"


const double Rideal = 8.3144621;
void updateParticleProperties(particleset  & vd, particleset  & vdmean, particleset  & dvdmeanx, int p, double dt, double l)
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
    double u_mean = vdmean.template getProp<i_velocity>(p)[0];
    double v_mean = vdmean.template getProp<i_velocity>(p)[1];
    double w_mean = vdmean.template getProp<i_velocity>(p)[2];
    double vel_mean = sqrt(u_mean*u_mean+v_mean*v_mean+w_mean*w_mean);

    //gradient particle properties
    double Tgrad = dvdmeanx.template getProp<i_temperature>(p);  //not used
    double Pgrad = dvdmeanx.template getProp<i_pressure>(p); 
    double rho_grad = dvdmeanx.template getProp<i_rho>(p);
    double energy_grad = dvdmeanx.template getProp<i_energy>(p);
    double u_grad = dvdmeanx.template getProp<i_velocity>(p)[0];
    double v_grad = dvdmeanx.template getProp<i_velocity>(p)[1];
    double w_grad = dvdmeanx.template getProp<i_velocity>(p)[2];
    double vel_grad = sqrt(u_grad*u_grad+v_grad*v_grad+w_grad*w_grad);

    //initialize deltas
    double drho, de, dm, dYk;
    double du, dh; //drift term
    //initialize new properties <-- dont need this bc should be += ?
    double rho_new, edensity_new, energy_new, Yk_new, T_new, P_new;
    double m_new, u_new;
    //intialize timescales
    double tau_sgs, tau_mol, tau_eq, tau_eq_energy;
    // initialize coefficients
    double Au_turb, Au_mol, Au_p, Au_p_alt, B;
    double Ae_p, h_p, h_mean, D_therm;

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //TESTING IN X DIRECTION ONLY//

    //check continuity to get density: 
    //drho/dt = -dmi/dxi -> drho = -dmi/dxi|dt; drho = -(sph gradient d(rho*u))/dx)*dt;
    drho = -(rho_grad * u_grad)*dt;
    vd.template getProp<i_rho>(p) += drho;
    cout << "drho = " << drho << " rho_grad, rho_mean = " << rho_grad << ", " << rho_mean << endl;
    cout << "New density = " << vd.template getProp<i_rho>(p) << " old density = "<< rho_p  << endl;
        
    // drift term
    du = (u_mean -  up);
    cout << "velx_grad = "<< u_grad  << " velx_mean = " << u_mean << endl;
      
    //turbulent kinetic energy is half the sum of the variances of the velocities
    double k_sgs = .5 * ((up - u_mean)*(up - u_mean));
    double Eps_sgs = k_sgs/l;

    //time scales
    double Cu = 2.1; //Kolmogorov constant <-- check this
    double C0 = 1;
    tau_sgs = k_sgs/Eps_sgs;
    tau_mol = (l*l)/(rho_p*up);    //l = kernel width (H)
    tau_eq = 1/((1/tau_mol)+(Cu/tau_sgs));

    //solve momentum equation
    // dm = -dP/dx|dt + Audt + BdW
    Au_turb = ((rho_p*Cu)/tau_sgs)*du;
    Au_mol = (rho_p/tau_mol)*du;
    Au_p = Au_mol + Au_turb;
    Au_p_alt = (rho_p/(tau_eq+dt))*du;   //Au_p and Au_p_alt (wo +dt) are equivalent
    B = C0*sqrt(Eps_sgs);   //turbulent diffusion

    //cout << "Aup = " << Au_p << " Au_p_alt= "<< Au_p_alt  << endl;

    //momentum equation
    m_p = rho_p * up;
    m_new = m_p - (Pgrad * dt) + (Au_p_alt * dt) + (B * dWien);

    cout << "New momentum = " << m_new << " old momentum = "<< m_p  << endl;

    //update velocity
    u_new = m_new / vd.template getProp<i_rho>(p);
    vd.template getProp<i_velocity>(p)[0] = u_new;

    cout << "New u velocity = " << u_new << " old u velocity = "<< up  << endl;

    //solve energy density
    //find specific enthalpy
    h_p = up - (Pp / rho_p);
    h_mean = u_mean - (Pmean / rho_mean);
    dh = h_mean - h_p;

    D_therm = (1);    //placeholder for thermal diffusivity (dependent on equivalence ratio)

    tau_eq_energy = 1/((D_therm/(l*l))+(Cu/tau_sgs));
    Ae_p = (rho_p/(tau_eq_energy + dt))*dh;

    //energy density equation
    edensity_p = rho_p * energy_p;
    edensity_new = edensity_p - (Pgrad * u_grad * dt) + (Ae_p * dt) + (B * dWien);

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