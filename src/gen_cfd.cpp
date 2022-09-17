#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

class Cfd
    {
    public:
        Cfd();
        ~Cfd();
        int nsteps;
        int nparticles;
        int frame;
        double dt,dx;
        double lx, ly, lz;  
        //engine eng; (outer class) 
        void displayCfd();
        void inputCfd(std::string filename);
    private:
        int NKEYINT;
        int  cfd_int[3];
        double cfd_double[3];    
    };
    Cfd::Cfd() // constructor
    {
    nsteps = 0; nparticles = 0; frame = 0;
    dt = 0.0;dx = 0.0;
    lx = 1.0; ly=1.0;lz=1.0;
    NKEYINT = 3;
    //cfd_int= new int [NKEYINT];
    for (int nv = 0; nv < NKEYINT ; nv++) cfd_int[nv]=0;    
    }
    Cfd::~Cfd() // destructor (take sno action)
    {    
    }
    // definition of display
    void Cfd::displayCfd()
    {
    cout << "------ CFD display  ------ \n"; 
    cout << " nsteps = " << nsteps << " \n"; 
    cout << " nparts = " << nparticles << " \n";  
    cout << " nframe = " << frame << " \n";         
    cout << "-------------------------- \n";     
    }  

    // read file
    void Cfd::inputCfd(std::string filename)
    {
    cout << " ** reading  input file:" << filename << "  \n"; 
    string item_name;
    ifstream nameFileout;
    nameFileout.open(filename);
    int NCHAR=6; //length of key workds 
    int NUMLENGTH=8; //length of number
    char cword[NUMLENGTH];
    string line;
    string cfd_keyint[3]={".nsteps",".nparts",".nframe"};
    int nline=0,nt;
    while(getline(nameFileout, line))
        {    
        // std::cout << "line  " << nline << " :" << line << std::endl; 
        nline++;    
       
        if (!line.empty() )
        { 
            if (line[0]=='%')
            {
            cout << " comment detected \n";
            }
            else
            {
            //cout << "line  " << nline << " :" << line << "\n";                      
            // look for key integer             
            for (int nv = 0; nv < NKEYINT ; nv++)
                {
                size_t found = line.find(cfd_keyint[nv]);
                if (found!=string::npos)
                    {
                    int pos = found; 
                    cout << cfd_keyint[nv] << " found at: " << pos << "\n";                  
                    for (int j = 0; j < NUMLENGTH ; j++) {cword[j]=line[pos+NCHAR+2+j];}
                    string sword(cword);                            
                    cfd_int[nv]= stoi(sword);
                    }
                }
            // look for key doubles                            
            }        
        } 
        }
        // assign values
        nparticles = cfd_int[0];
        nsteps     = cfd_int[1];
        frame      = cfd_int[2];
    }

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
