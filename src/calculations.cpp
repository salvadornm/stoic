#include "calculations.h"

// calculate the property means
void calculateMeans(particleset & vd, Cfd &simulation)
{   
    double Tmean = 0;
    double Pmean = 0;
    double Rhomean = 0;
    double count = 0;
    //calculate mean temperature from previous iteration
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        Tmean += vd.template getProp<i_temperature>(a);
        Pmean += vd.template getProp<i_pressure>(a);
        Rhomean += vd.template getProp<i_rho>(a);
        
        ++it; ++count;
    }

    simulation.Pmean = Pmean / count;
    simulation.Tmean = Tmean / count;
    simulation.Rhomean = Rhomean / count;
}


// update initial cylinder properties for means
void updateInitialProps(particleset & vd, Cfd &simulation)
{
    calculateMeans(vd, simulation);

    simulation.P0 = simulation.Pmean;
    simulation.T0 = simulation.Tmean;
}


//solve for density and energy from Pressure and Temperature
// rho = (P)/(R*T)
void updateThermalProperties2(particleset & vd, int a)
{
    double P = vd.template getProp<i_pressure>(a);
    double T = vd.template getProp<i_temperature>(a);

    vd.template getProp<i_rho>(a) = P/(R_air*T);        //ideal gas law: assume dry air for now
    vd.template getProp<i_energy>(a) = T*cv_global;
}


//solve for Pressure and Temperature from density and energy 
// P = (rho*R*T)/M
void updateThermalProperties1(particleset & vd, int a)
{
    double rho = vd.template getProp<i_rho>(a);
    double e_int = vd.template getProp<i_energy>(a);
    double k = R_air/cp_global;

    //vd.template getProp<i_temperature>(a) = vd.template getProp<i_vdmean>(a)[i_temperature] - k * vd.template getProp<i_energy>(a);
    vd.template getProp<i_temperature>(a) = e_int/cv_global;
    vd.template getProp<i_pressure>(a) = rho*(R_air*vd.template getProp<i_temperature>(a));        //assume dry air for now
}


// calculate mass of the volume in the cylinder
double calculateMass(particleset & vd, engine eng)
{
    auto it = vd.getDomainIterator();
    double rho_tot;
    int count = 0;
    while (it.isNext())
    {
        auto a = it.get();
        double rho_a = vd.template getProp<i_rho>(a);
        rho_tot += rho_a;
        
        ++it; ++count;
    }
    rho_tot = rho_tot/count;    //find the average density
    double mass_mix = rho_tot/eng.volumeC;
    return mass_mix;
}

// calculate heat loss from interaction with walls of cylinder
double calculateHeatloss(particleset & vd, Cfd simulation, engine eng)
{
    double rho = simulation.m_tot / eng.volumeC;
    double P = rho*(R_air*simulation.Tmean);

    // Monitored pressure  (from old code)(in vgas would be P-Pm )
    /*
    Pm_old=Pm
    Pm=Pref*(eng.VBDC/eng.volumeC)**GAMA
    dPm=Pm-Pm_old
    */

    // Compute gas velocity
    double v_gas = 2.28*eng.smp + (0.00324/6.0)*simulation.T0*(eng.Vdisp/eng.VBDC)*(P)/simulation.P0;

    // Heat Transfer coefficient; modified Woschni correlation by Chang (2004)
    double h_coeff = 3.4 * pow(P,0.80) * pow(v_gas,0.80) * pow((simulation.ly-eng.dStime),-0.20) * pow(simulation.Tmean,-0.73) ; 
    // P(t), T(t), do not depend on position h_coeff is not local, global. [W/m²K

    double Tdelta = 0; double heat_transfer = 0;
    double SA = 0; double Q = 0; int count2 = 0;

    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        Point<3,double> xa = vd.getPos(a);
        // set up geometries
        double r_cyl = eng.bore/2;
        Point<3,double> cyl_center {r_cyl, r_cyl, xa[2]};
        double r_pos = sqrt(norm2(xa - cyl_center));

        //if @ walls, then:
        if ( r_cyl - r_pos < 1e-8){
        SA = (2*pi*r_cyl) * eng.height;
        } else if((eng.height - xa[2] < 1e-8) || ((xa[2] - eng.s_inst) < 1e-8) ){
        // else if @ top/bottom
        SA = pow(r_cyl,2) * (pi);
        }

        Tdelta = vd.template getProp<i_temperature>(a) - eng.Twall;
        Q = h_coeff * Tdelta  * simulation.dt;
        vd.template getProp<i_energy>(a) = vd.template getProp<i_energy>(a) - (Q* SA) -  (vd.template getProp<i_pressure>(a) / vd.template getProp<i_rho>(a));
    
        //update temperature based on new enthalpy/energy
        vd.template getProp<i_temperature>(a) = vd.template getProp<i_energy>(a)/cv_global;

        heat_transfer += Q;
        count2++; ++it;
    }

    // mean heat transferred of a particle per unit area [J/m2]
    heat_transfer = heat_transfer/count2;
    double heat_flux = heat_transfer/(simulation.dt*10e5);    //[MW/m2]
    
    return heat_transfer;
}


//calculate pressure based on local density (EqState in OpenFPM)
// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;
double B = 20.0*20.0*9.81*1.0*rho_zero / 7.0; //(coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;

void updateEqtnState(particleset & vd)
{
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        double rho_a = vd.template getProp<i_rho>(a);
        double rho_frac = rho_a / rho_zero;
        vd.template getProp<i_pressure>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);
        ++it;
    }
}

//implements the viscous term (from OpenFPM)
double cbar = 20.0 * sqrt(9.81*1.0);    //cbar in formula //coeff_sound * sqrt(gravity * h_swl);
const double visco = 0.1;               // alpha in the formula

double viscous(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double visc, Cfd simulation) //double & visc)
{
    double Eta2 = 0.01 * simulation.H * simulation.H;
    const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
    const double dot_rr2 = dot/(rr2+Eta2);
    visc=std::max(dot_rr2,visc);
    if(dot < 0)
    {
        const float amubar=simulation.H*dot_rr2;
        const float robar=(rhoa+rhob)*0.5f;
        const float pi_visc=(-visco*cbar*amubar/robar);
        return pi_visc;
    }
    else
        return 0.0;
}
