#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include "timer.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

constexpr int velocity = 0;
constexpr int force = 0;

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


int main(int argc, char* argv[])
{
    cout << "Hello World \n" << endl;
    
    //** VARIABLES **//
	double sigma = 0.1;
	double r_cut = 3.0*sigma;   //cutoff radius
    //time step integration:
	double dt = 0.01;
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
	//vector_dist<3,double, aggregate<float,double[3],double[3]> > vd(0,domain,bc,ghost);
    //typedef vector_dist<3,double,aggregate<size_t,double,double,double,double,double[3],double[3], double[3]>> particles;
	vector_dist<3,float, aggregate<float,float[3],float[3][3]> > vd(4096,domain,bc,g); //each processor runs with 4096 particles
    // the scalar is the element at position 0 in the aggregate
    const int scalar = 0;
    // the vector is the element at position 1 in the aggregate
    const int vector = 1;
    // the tensor is the element at position 2 in the aggregate
    const int tensor = 2;

    size_t cnt = 0; //used later

    //**ASSIGN POSITION**//
    //should generate all values
    auto it = vd.getGridIterator(sz);
    while (it.isNext())
    {
        vd.add();   //adds particle?
        auto key = it.get();    //contains (i,j,k) index of grid
        
        // get position of the last particle added
		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

        // set the property values of the last particle we added
		vd.template getLastProp<velocity>()[0] = 01.0;
		vd.template getLastProp<velocity>()[1] = 1.0;
		vd.template getLastProp<velocity>()[2] = 1.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;
        
        // next particle
        ++it;
        cnt++;
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
	//calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);
	unsigned long int f = 0;

// MD time stepping
	for (size_t i = 0; i < 100 ; i++)
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
		//calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);


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
        if (i % 10 == 0)
        {
            // Write the particle position for visualization (Without ghost)
            vd.deleteGhost();
            vd.write_frame("particles_",f);

            // resync the ghost
            vd.ghost_get<>();

            //**Reduce**//
            auto & v_cl = create_vcluster();
            v_cl.sum(cnt);
            v_cl.execute();

            f++;
        }
	}

	tsim.stop();
	std::cout << "Time: " << tsim.getwct() << std::endl;


    //**VISUALIZATION**// <--- NEEDS UPDATED STILL
    //A VTK file contains information about particle position and properties

    openfpm::vector<std::string> names({"velocity","force"});
    vd.setPropNames(names);
    // save vtk format (vtk is always the default)
    vd.write("particles_moving");
    // save in vtk format with time
    //vd.write("particles_moving_with_time","time=1.234");
    // save in vtk format with time
    vd.write("particles_moving_with_time_bin","time=1.234",VTK_WRITER | FORMAT_BINARY);
    
    //De - initialize openfpm
    openfpm_finalize();
}
