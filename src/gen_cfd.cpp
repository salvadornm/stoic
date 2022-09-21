#include "class_cfd.hpp"

   
///// TESTING /////////////////////////////////
int main()
{
   Cfd sim;
   string filename("input.txt");

   cout << " hello ....  \n";
   sim.displayCfd();

   sim.inputCfd(filename);

   sim.displayCfd();

  return 0;
}
