#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

constexpr int velocity = 0;
constexpr int force = 1;

/*
 * Calculate the forces between particles. Requires the vector of particles
 * Cell list and scaling factor for the Lennard-Jhones potential.
 */
template<typename CellList> void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6, double r_cut2)
{
	vd.updateCellList(NN);
	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();
		Point<3,double> xp = vd.getPos(p);  // Get position xp of the particle

		// Reset the force counter
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Iterate over the neighborhood particles of p
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			auto q = Np.get();
			if (q == p.getKey())	{++Np; continue;};

			Point<3,double> xq = vd.getPos(q);
			Point<3,double> r = xp - xq;    // Get the distance between p and q

			double rn = norm2(r);   // take the norm of this vector

			if (rn > r_cut2)    {++Np; continue;};

			// Calculate the force
			Point<3,double> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

			// sum the force produced by q on p
			vd.template getProp<force>(p)[0] += f.get(0);
			vd.template getProp<force>(p)[1] += f.get(1);
			vd.template getProp<force>(p)[2] += f.get(2);

			++Np;
		}
		++it2;
	}
}


//Calculate energy. Requires the same parameter as calculate forces
template<typename CellList> double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6, double r_cut2)
{
    double rc = r_cut2;
    double shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	double E = 0.0; //zero energy

	/*!
	 * Iterate over the particles and get its position. 
	 * Calculate the energy based on the Lennard-Jhones potential:
	 *
	 * \f$ V(x_p,x_q) = 4(\frac{1}{r^{12}} - \frac{1}{r^{6}}) \f$
	 */

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,double> xp = vd.getPos(p);
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			auto q = Np.get();
			if (q == p.getKey())	{++Np; continue;};

			Point<3,double> xq = vd.getPos(q);// Get position of the particle q
			double rn = norm2(xp - xq);// take the normalized direction

			if (rn > r_cut2)    {++Np;continue;}

			// Calculate potential energy
			E += 2.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

			++Np;
		}

		// Kinetic energy of the particle given by its actual speed
		E +=   (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
				vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
				vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;

		++it2;
	}

	return E;   // Calculated energy
}


int main(int argc, char* argv[])
{
    cout << "Hello World \n" << endl;
    cout << rand() % 10;
    
    //** VARIABLES **//
	double sigma = 0.1;
	double r_cut = 3.0*sigma;   //cutoff radius
    //time step integration:
	double dt = 0.0005;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

    //** INITIALIZATION **//
    openfpm_init(&argc,&argv);

    Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});	//define 3D box
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};    // boundary conditions
    Ghost<3,double> ghost(r_cut);   // extended to contain interaction radius
    
	size_t sz[3] = {10,10,10};// place particles on a 10x10x10 Grid like

    //** VECTOR INSTANTIATION **//
    //contains two vectoral properties
	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,domain,bc,ghost);
    

    //**ASSIGN POSITION**//
    //should generate all values
    auto it = vd.getGridIterator(sz);
    while (it.isNext())
    {
        vd.add();   //???????
        auto key = it.get();    //contains (i,j,k) index of grid
        
        // get position of the last particle added
		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

        // set the property values of the last particle we added
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;
        
        // next particle
        ++it;
    }

    //** Map Particles to grid **/
    vd.map();           //??????
	vd.ghost_get<>();   //syncs the ghost with the newly mapped particles


    //**ASSIGNING VALUES**//
    //here we do 10000 MD steps using verlet integrator

	timer tsim;
	tsim.start();

	auto NN = vd.getCellList<CELL_MEMBAL(3,double)>(r_cut); //cell list structure


	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);
	unsigned long int f = 0;

// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		auto it3 = vd.getDomainIterator();

		// (1) - integrate velocity and space based on the calculated forces
		while (it3.isNext())
		{
			auto p = it3.get();

			// calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			// calculate x(tn + 1)
			vd.getPos(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
			vd.getPos(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
			vd.getPos(p)[2] += vd.template getProp<velocity>(p)[2]*dt;

			++it3;
		}

		// Map particles and re-sync the ghost
		vd.map();
		vd.template ghost_get<>();

		// (2) - calculate forces or a(tn + 1)
		calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);


		// (3) - Integrate the velocity
		auto it4 = vd.getDomainIterator();

		while (it4.isNext())
		{
			auto p = it4.get();

			//calculate v(tn + 1)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			++it4;
		}

        // collect some statistic about the configuration
        if (i % 100 == 0)
        {
            // Write the particle position for visualization (Without ghost)
            vd.deleteGhost();
            vd.write_frame("particles_",f);

            // resync the ghost
            vd.ghost_get<>();

            // calculate the energy
            double energy = calc_energy(vd,NN,sigma12,sigma6,r_cut*r_cut);
            auto & vcl = create_vcluster();
            vcl.sum(energy);
            vcl.execute();

            // save the energy calculated at time step i c contain the time-step y contain the energy
            x.add(i);
            y.add({energy});

            // We also print on terminal the value of the energy
            // only one processor (master) write on terminal
            if (vcl.getProcessUnitID() == 0)
                    std::cout << "Energy: " << energy << std::endl;

            f++;
        }
	}

	tsim.stop();
	std::cout << "Time: " << tsim.getwct() << std::endl;

    
    //**REDUCE**//
    

    //**VISUALIZATION**// <--- NEEDS UPDATED STILL
    //A VTK file contains information about particle position and properties

    openfpm::vector<std::string> names({"scalar","vector","tensor"});
    vd.setPropNames(names);
    // save vtk format (vtk is always the default)
    vd.write("particles");
    // save in vtk format with time
    vd.write("particles_with_time","time=1.234");
    // save in vtk format with time
    vd.write("particles_with_time_bin","time=1.234",VTK_WRITER | FORMAT_BINARY);
    
    //De - initialize openfpm
    openfpm_finalize();
}
