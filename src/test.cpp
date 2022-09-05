#include "test.h"

void kernel_test( double H, Point<3,double> dr)
{
    // test kernel function
    std::cout << " kernel(0.0 H)  " << Wab(0.0,H) << std::endl;
    std::cout << " kernel(0.99 H) " << Wab(0.99*H,H)<< std::endl;
    std::cout << " kernel(1.01 H) " << Wab(1.01*H,H)<< std::endl;
    std::cout << " kernel(1.9 H)  " << Wab(1.9*H,H) << std::endl;

    double r = -0.1;
    Point<3,double> DW;
    double factor = DWab(dr,DW,r,H); // gradient kernel //
    double W = Wab(r,H); //kernel
    
    double a = 1.0/(pi*pow(H,3));
    double a_h = a/H;

    while (r < 2)
    {
        r += 0.1;
        W = Wab(r*H,H); //kernel
        cout << " H = " << H << endl;
        factor = DWab(dr,DW,r*H,H);
        cout << "r= " << r << endl;
        cout << "dW = " << factor << ", W = " << W << endl;
        cout << "dW/(a/h) = " << factor/a_h << ", W/a = " << W/a << endl;
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

    double c1x = vd.getProp<i_dvdmean>(p)[0][i_momentum];
    double c2x = vd.getProp<i_dvdmean>(p)[0][i_temperature]; 
    double c3x = vd.getProp<i_dvdmean>(p)[0][i_pressure];
    double c4x = vd.getProp<i_dvdmean>(p)[0][i_energy];

    double c1y = vd.getProp<i_dvdmean>(p)[1][i_momentum];
    double c2y = vd.getProp<i_dvdmean>(p)[1][i_temperature]; 
    double c3y = vd.getProp<i_dvdmean>(p)[1][i_pressure];
    double c4y = vd.getProp<i_dvdmean>(p)[1][i_energy];

    double c1z = vd.getProp<i_dvdmean>(p)[2][i_momentum];
    double c2z = vd.getProp<i_dvdmean>(p)[2][i_temperature]; 
    double c3z = vd.getProp<i_dvdmean>(p)[2][i_pressure];
    double c4z = vd.getProp<i_dvdmean>(p)[2][i_energy];

    //print
    std::cout << "(vd particle) --------" << std::endl;
    std::cout << " temp = " << a2 << " p = " << a3 << std::endl;
    std::cout << "density = " << a1 << " energy = " << a4 << std::endl;
    std::cout << " vx = " << a5 << " vy = " << a6 << " vz = " << a7 << std::endl;
    std::cout << " px = " << vd.getPos(p)[0] << " py = " << vd.getPos(p)[1] << " pz = " << vd.getPos(p)[2] << std::endl;
    
    
    std::cout << "(vdmean particle) --------" << std::endl;
    std::cout << " temp = " << b2 << " p = " << b3 << std::endl;
    std::cout << "density = " << b1 << " energy = " << b4 << std::endl;
    std::cout << " vx = " << b5 << " vy = " << b6 << " vz = " << b7 << std::endl;

    std::cout << "(dvdmean particle) --------" << std::endl;
    std::cout << " temp = " << c2x << " p = " << c3x << std::endl;
    std::cout << " mom = " << c1x << " visc*P = " << c4x << std::endl;

    std::cout << "(dvdmean_y particle) --------" << std::endl;
    std::cout << " temp = " << c2y << " p = " << c3y << std::endl;
    std::cout << " mom = " << c1y << " visc*P = " << c4y << std::endl;

    std::cout << "(dvdmean_z particle) --------" << std::endl;
    std::cout << " temp = " << c2z << " p = " << c3z << std::endl;
    std::cout << " mom = " << c1z << " visc*P = " << c4z << std::endl;
    
    
}

void output_properties(double mom_p, double drho, double rho_new, double Au_p, double dWien)
{
    std::cout << "momentum = " << mom_p << endl;
    std::cout << "drho = " << drho << " rho_new = " << rho_new << endl;
    std::cout << "Au_p: " << Au_p << endl;
    std::cout << "dwien: " << dWien << endl;
    //std::cout << "New u velocity = " << vd.template getProp<i_velocity>(p)[i] << endl
}

void output_bc_props(vector<double> vel, Point <3,double> pos, Point <3,double> pos_new, Point <3,double> pos_wall, engine eng, Point<3,double> & psi){
    std::cout << "(bc / movement) --------" << std::endl;
    cout << "point og: " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;
    cout << "point new: " << pos_new[0] << ", " << pos_new[1] << ", " << pos_new[2] << endl;
    cout << "wall p: " << pos_wall[0] << ", " << pos_wall[1] << ", " << pos_wall[2] << endl;
    cout << "psi: " << psi[0] << ", " << psi[1] << ", " << psi[2] << endl;
    cout << "vel x: " << vel[0] << " vel y: " << vel[1] << " vel z: " << vel[2] << endl;       
}

void output_energy_props(particleset &vd, int p, double dh, double dvisc, double edensity_p, double edensity_new, double energy_new){
    cout << "dh: " << dh << " dviscP: " << dvisc << endl;
    cout << "edensity_p: " << edensity_p << " edensity_new: " << edensity_new << endl;
    cout << "energy: " << energy_new << " internal energy: " << vd.template getProp<i_energy>(p) << endl;
    cout << "New temp: " << vd.template getProp<i_temperature>(p);
    cout << " New pressure: " << vd.template getProp<i_pressure>(p) << endl;
}


void vary_initialization(particleset &vd, Cfd simulation, int key)
{  
    //update the inputs as desired      
    if (vd.getPos(key)[0] < (0.8*simulation.lx) && vd.getPos(key)[1] < (0.8*simulation.ly) && vd.getPos(key)[2] < (0.8*simulation.lz))
    {
        vd.template getProp<i_velocity>(key)[0] = 0.1;
        vd.template getProp<i_velocity>(key)[1] = 0.1;
        vd.template getProp<i_velocity>(key)[2] = 0.1;
        vd.template getProp<i_rho>(key) = .9; 
    }
}
void limit_velocity(particleset &vd, int key, int i)
{  
    double maxVelComp = 1.0;
    //update the inputs as desired      
    if (vd.template getProp<i_velocity>(key)[i] > maxVelComp) 
    {
        vd.template getProp<i_velocity>(key)[i] = maxVelComp;
        cout << "velocity limited to 0.5m/s" << endl;
    } else if (vd.template getProp<i_velocity>(key)[i] < -maxVelComp) {
      vd.template getProp<i_velocity>(key)[i] = -maxVelComp;
        cout << "velocity limited to -0.5m/s" << endl;
    }
}

void outputdata_to_csv(particleset vd, int step){
    ofstream outputfile("stoic_output.csv", std::ios::app);
    //outputfile.open("stoic_output.csv");
    outputfile << "Step,ID,Posx,Posy,Posz,Vx,Vy,Vz,P,T,rho,U_in,Vx_avg,Vy_avg,Vz_avg,Pavg,Tavg,rho_avg,Uin_avg,Mom_gradx,Mom_grady,Mom_gradz,P_gradx,P_grady,P_gradz,T_gradx,T_grady,T_gradz,\n";
    
    auto part = vd.getDomainIterator();
    while(part.isNext())
    {
        auto p = part.get();
        int id = p.getKey();

        double pos0 = vd.getPos(p)[0];
        double pos1 = vd.getPos(p)[1];
        double pos2 = vd.getPos(p)[2];

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

    double c1x = vd.getProp<i_dvdmean>(p)[0][i_momentum];
    double c2x = vd.getProp<i_dvdmean>(p)[0][i_temperature]; 
    double c3x = vd.getProp<i_dvdmean>(p)[0][i_pressure];
    double c4x = vd.getProp<i_dvdmean>(p)[0][i_energy];

    double c1y = vd.getProp<i_dvdmean>(p)[1][i_momentum];
    double c2y = vd.getProp<i_dvdmean>(p)[1][i_temperature]; 
    double c3y = vd.getProp<i_dvdmean>(p)[1][i_pressure];
    double c4y = vd.getProp<i_dvdmean>(p)[1][i_energy];

    double c1z = vd.getProp<i_dvdmean>(p)[2][i_momentum];
    double c2z = vd.getProp<i_dvdmean>(p)[2][i_temperature]; 
    double c3z = vd.getProp<i_dvdmean>(p)[2][i_pressure];
    double c4z = vd.getProp<i_dvdmean>(p)[2][i_energy];

    outputfile << step << "," << id << ",";
    outputfile << pos0 << "," << pos1 << "," << pos2 << "," << a5 << "," << a6 << "," << a7 << "," << a3 << "," << a2 << "," << a1 << "," << a4 << ",";
    outputfile << b5 << "," << b6 << "," << b7 << "," << b3 << "," << b2 << "," << b1 << "," << b4 << ",";
    outputfile << c1x << "," << c1y << "," << c1z << "," << c3x << "," << c3y << "," << c3z << "," << c2x << "," << c2y << "," << c2z << "\n";
    
    ++part;
    }

}