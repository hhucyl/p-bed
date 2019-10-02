#include "./lbm/Domain.h"
#include <time.h>






double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 10;
    double nu = 5e-4;

    //rwP
    int RWP = 5e4;
    if(argc>=2) Nproc = atoi(argv[0]);
    
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double Dm = 0.01;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    
    dom.Nproc = Nproc;       

    //initial
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,0.0,0.0);
    dom.IsF = true;
    
    dom.Initial(1.0,v0,g0);
   
    //RWParticles
    if(!dom.IsRWContinue)
    {
        int nn = 0;
        
        for(int ip=0; ip<RWP; ++ip)
        {
            Vec3_t xt(49,49,0);
            int xx = std::round(xt(0));
            int yy = std::round(xt(1));
            if(dom.Gamma[xx][yy][0]<1e-12)
            {
                dom.RWParticles.push_back(RW::Particle(xt,2,Dm));
                dom.RWParticles.back().Tag = nn;
                nn++;
            }
        } 
        
        std::cout<<"RW particles complete "<<std::endl;
        std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;
    }else{
        for(size_t ip=0; ip<dom.RWParticles.size(); ++ip)
        {
            dom.RWParticles[ip].Dm = Dm;
        }
        std::cout<<"Assign Dm"<<std::endl;
        
    }
    


    double Tf = 1e4;
    double dtout = 100;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolvePRW( Tf, dtout, "test_rw2", NULL, NULL);
    
    return 0;
}MECHSYS_CATCH
