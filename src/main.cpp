#include <stddef.h>
#include "Vector/vector_dist.hpp"
int main(int argc, char* argv[])
{
    // initialize the library
    openfpm_init(&argc,&argv);
    // Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
    Box<2,float> domain({0.0,0.0},{1.0,1.0});
    // Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<2,float> g(0.01);
    
    vector_dist<2,float, aggregate<float,float[3],float[3][3]> > vd(4096,domain,bc,g);
    // the scalar is the element at position 0 in the aggregate
    const int scalar = 0;
    // the vector is the element at position 1 in the aggregate
    const int vector = 1;
    // the tensor is the element at position 2 in the aggregate
    const int tensor = 2;
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();
        // we define x, assign a random position between 0.0 and 1.0
        vd.getPos(key)[0] = (float)rand() / RAND_MAX;
        // we define y, assign a random position between 0.0 and 1.0
        vd.getPos(key)[1] = (float)rand() / RAND_MAX;
        // next particle
        ++it;
    }
    vd.map();
    //Counter we use it later
    size_t cnt = 0;
    // Get a particle iterator
    it = vd.getDomainIterator();
    // For each particle ...
    while (it.isNext())
    {
        // ... p
        auto p = it.get();
        // we set the properties of the particle p
        
         // the scalar property
        vd.template getProp<scalar>(p) = 1.0;
        vd.template getProp<vector>(p)[0] = 1.0;
        vd.template getProp<vector>(p)[1] = 1.0;
        vd.template getProp<vector>(p)[2] = 1.0;
        vd.template getProp<tensor>(p)[0][0] = 1.0;
        vd.template getProp<tensor>(p)[0][1] = 1.0;
        vd.template getProp<tensor>(p)[0][2] = 1.0;
        vd.template getProp<tensor>(p)[1][0] = 1.0;
        vd.template getProp<tensor>(p)[1][1] = 1.0;
        vd.template getProp<tensor>(p)[1][2] = 1.0;
        vd.template getProp<tensor>(p)[2][0] = 1.0;
        vd.template getProp<tensor>(p)[2][1] = 1.0;
        vd.template getProp<tensor>(p)[2][2] = 1.0;
        // increment the counter
        cnt++;
        // next particle
        ++it;
    }
    
    auto & v_cl = create_vcluster();
    v_cl.sum(cnt);
    v_cl.execute();
    
    openfpm::vector<std::string> names({"scalar","vector","tensor"});
    vd.setPropNames(names);
    // save vtk format (vtk is always the default)
    vd.write("particles");
    // save in vtk binary format
    vd.write("particles_bin",VTK_WRITER | FORMAT_BINARY);
    // save in vtk format with time
    vd.write("particles_with_time","time=1.234");
    // save in vtk format with time
    vd.write("particles_with_time_bin","time=1.234",VTK_WRITER | FORMAT_BINARY);
    openfpm_finalize();
}
