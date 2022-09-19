#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

constexpr int NKEYINT=3;    // number of key word with integer
constexpr int NKEYFLOAT=4;  // number of key word with float
string cfd_keyint[NKEYINT]={".nsteps",".nparts",".nframe"};
string cfd_keydbl[NKEYFLOAT]={".dt",".lx",".ly",".lz"};    


class Cfd
    {
    public:
        Cfd();
        ~Cfd();
        int nsteps;
        int nparticles;
        int frame;
        double dt;
        double lx, ly, lz;  
        //engine eng; (outer class) 
        void displayCfd();
        void inputCfd(std::string filename);
    private:
        int    cfd_int[NKEYINT];
        double cfd_dbl[NKEYFLOAT];  
        void assignCfd();  
    };
    Cfd::Cfd() // constructor
    {        
        for (int nv = 0; nv < NKEYINT   ; nv++) cfd_int[nv]=0;    
        for (int nv = 0; nv < NKEYFLOAT ; nv++) cfd_dbl[nv]=0;   
        assignCfd();     
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
    cout << " dt = " << dt << " \n";             
    cout << " lx = " << lx << "  ly = " << ly << "  lz = " << lz << "\n";             
    cout << "-------------------------- \n";     
    }  

    // fill public variables  (change if you add more keywords)
    void Cfd::assignCfd()
    {
        nparticles = cfd_int[0];
        nsteps     = cfd_int[1];
        frame      = cfd_int[2];
        dt  = cfd_dbl[0];
        lx = cfd_dbl[1];ly = cfd_dbl[2];lz = cfd_dbl[3];   
    }

    // read file
    void Cfd::inputCfd(std::string filename)
    {
    cout << " ** reading  input file:" << filename << "  \n"; 
    string item_name;
    ifstream nameFileout;
    nameFileout.open(filename);
    int NCHAR=6,NCHAR2=2; //length of key workds (integer and float)
    int NUMLENGTH=8; //length of number
    char cword[NUMLENGTH];
    string line;
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
            for (int nv = 0; nv < NKEYFLOAT ; nv++)
                {
                size_t found = line.find(cfd_keydbl[nv]);
                if (found!=string::npos)
                    {
                    int pos = found; 
                    cout << cfd_keydbl[nv] << " found at: " << pos << "\n";                  
                    for (int j = 0; j < NUMLENGTH ; j++) {cword[j]=line[pos+NCHAR2+2+j];}
                    string sword(cword);                            
                    cfd_dbl[nv]= stod(sword);
                    }
                }

            }        
        } 
        }        
        assignCfd(); // assign values
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
