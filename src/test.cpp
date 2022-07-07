#include "test.h"

void kernel_test( double H, Point<3,double> dr)
{
    // test kernel function
    std::cout << " kernel(0.0 H)  " << Wab(0.0) << std::endl;
    std::cout << " kernel(0.99 H) " << Wab(0.99*H)<< std::endl;
    std::cout << " kernel(1.01 H) " << Wab(1.01*H)<< std::endl;
    std::cout << " kernel(1.9 H)  " << Wab(1.9*H) << std::endl;

    double r = -0.1;
                Point<3,double> DW;
                double factor = DWab(dr,DW,r,false); // gradient kernel //
                double W = Wab(r); //kernel
    while (r < 2)
    {
        r += 0.1;
        W = Wab(r*H); //kernel
        //cout << " W = " << W << endl;
        factor = DWab(dr,DW,r*H,false);
        cout << "r= " << r << endl;
        cout << "dW = " << factor << ", W/a = " << W << endl;
        cout << "dWab = " << DW.get(0) << " , " << DW.get(1) << " , " << DW.get(2) << endl;
    }
}

void output_vd(particleset  & vd, int p)
{
    //test outputs
    double a1 = vd.template getProp<i_rho>(p);
    double a2 = vd.getProp<i_temperature>(p);
    double a3 = vd.getProp<i_pressure>(p);
    double a4 = vd.getProp<i_energy>(p);
    double a5 = vd.getProp<i_velocity>(p)[0];
    double a6 = vd.getProp<i_velocity>(p)[1];
    double a7 = vd.getProp<i_velocity>(p)[2];

    double b1 = vd.getProp<i_vdmean>(p)[i_rho];
    double b2 = vd.getProp<i_vdmean>(p)[i_temperature]; 
    double b3 = vd.getProp<i_vdmean>(p)[i_pressure];
    double b4 = vd.getProp<i_vdmean>(p)[i_energy];
    double b5 = vd.getProp<i_vdmean>(p)[i_velx];
    double b6 = vd.getProp<i_vdmean>(p)[i_vely];
    double b7 = vd.getProp<i_vdmean>(p)[i_velz];

    double c1x = vd.getProp<i_dvdmean>(p)[0][i_rho];
    double c2x = vd.getProp<i_dvdmean>(p)[0][i_temperature]; 
    double c3x = vd.getProp<i_dvdmean>(p)[0][i_pressure];
    double c4x = vd.getProp<i_dvdmean>(p)[0][i_energy];

    double c1y = vd.getProp<i_dvdmean>(p)[1][i_rho];
    double c2y = vd.getProp<i_dvdmean>(p)[1][i_temperature]; 
    double c3y = vd.getProp<i_dvdmean>(p)[1][i_pressure];
    double c4y = vd.getProp<i_dvdmean>(p)[1][i_energy];

    double c1z = vd.getProp<i_dvdmean>(p)[2][i_rho];
    double c2z = vd.getProp<i_dvdmean>(p)[2][i_temperature]; 
    double c3z = vd.getProp<i_dvdmean>(p)[2][i_pressure];
    double c4z = vd.getProp<i_dvdmean>(p)[2][i_energy];

    //print
    std::cout << "(vd particle) --------" << std::endl;
    std::cout << " temp = " << a2 << " p = " << a3 << std::endl;
    std::cout << "density = " << a1 << " energy = " << a4 << std::endl;
    std::cout << " vx = " << a5 << " vy = " << a6 << " vz = " << a7 << std::endl;

    std::cout << "(vdmean particle) --------" << std::endl;
    std::cout << " temp = " << b2 << " p = " << b3 << std::endl;
    std::cout << "density = " << b1 << " energy = " << b4 << std::endl;
    std::cout << " vx = " << b5 << " vy = " << b6 << " vz = " << b7 << std::endl;

    std::cout << "(dvdmean particle) --------" << std::endl;
    std::cout << " temp = " << c2x << " p = " << c3x << std::endl;
    std::cout << "density = " << c1x << " energy = " << c4x << std::endl;

/*
    std::cout << "(dvdmean_y particle) --------" << std::endl;
    std::cout << " temp = " << c2y << " p = " << c3y << std::endl;
    std::cout << "density = " << c1y << " energy = " << c4y << std::endl;

    std::cout << "(dvdmean_z particle) --------" << std::endl;
    std::cout << " temp = " << c2z << " p = " << c3z << std::endl;
    std::cout << "density = " << c1z << " energy = " << c4z << std::endl;
    */
}