#include "calculations.h"

// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;
double B = 20.0*20.0*9.81*1.0*rho_zero / 7.0; //(coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;

//calculate pressure based on local density (called EqState in openfpm)
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

//solve for density and energy from Pressure and Temperature
// rho = (M*P)/(R*T)
void updateDensity(particleset & vd)
{
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        double P = vd.template getProp<i_pressure>(a);
        double T = vd.template getProp<i_temperature>(a);
        vd.template getProp<i_rho>(a) = P/(R_air*T);        //ideal gas law: assume dry air for now
        ++it;
    }
}

//solve for Pressure and Temperature from density and energy 
// P = (rho*R*T)/M
void updateThermalProperties1(particleset & vd, particleset & vdmean, int a)
{
    double rho = vd.template getProp<i_rho>(a);
    //k = vd.template getProp<i_species>(a);
    double k = 1;   //temporary
    vd.template getProp<i_temperature>(a) = vdmean.template getProp<i_temperature>(a) - k * vd.template getProp<i_energy>(a);
    vd.template getProp<i_pressure>(a) = rho*(R_air*vd.template getProp<i_temperature>(a));        //assume dry air for now
}

double cbar = 20.0 * sqrt(9.81*1.0); //cbar in formula //coeff_sound * sqrt(gravity * h_swl);
const double visco = 0.1; // alpha in the formula

//implements the viscous term
double viscous(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc)
{
    const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
    const double dot_rr2 = dot/(rr2+Eta2);
    visc=std::max(dot_rr2,visc);
    if(dot < 0)
    {
        const float amubar=H*dot_rr2;
        const float robar=(rhoa+rhob)*0.5f;
        const float pi_visc=(-visco*cbar*amubar/robar);
        return pi_visc;
    }
    else
        return 0.0;
}
