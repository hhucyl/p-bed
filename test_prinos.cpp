#include "./lbm/Domain.h"
#include <time.h>

struct myUserData
{
    std::ofstream oss_ss;
    double y_start;
    double vmax;
    double H;
    Vec3_t g0;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    // size_t nx = dom.Ndim(0);
    // size_t ny = dom.Ndim(1);
    // size_t nz = dom.Ndim(2);
    // double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","mvbed");
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Collision"<<Util::_8s<<"GhostNum\n";
    }else{
        
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.ListofContacts.size()<<Util::_8s<<dom.GhostParticles[0].Ghost<<std::endl;
        
    }
}

void Setup(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0; ix<nx; ++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        f[4] = f[2];
        f[7] = f[6];
        f[8] = f[5];
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
    }
    
    // double V = 0;
    // double num = 1;    
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iy=0; iy<ny; ++iy)
    // {
    //     if(dom.Gamma[ix][iy][0]<1) 
    //     {
    //         V += dom.Vel[ix][iy][0](0);
    //         num += 1; 
    //     }     
    // }
    // V /= (num);
    // // std::cout<<dat.V - V<<std::endl;
    // Vec3_t g(dat.V - V,0,0);
    // #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iy=0; iy<ny; ++iy)
    // {
    //     dom.BForce[ix][iy][0] = g;

    // }
    
    // std::cout<<g(0)<<std::endl;
}


void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);

    double py = dat.y_start;
    double P = 2*dat.vmax/(dat.H*dat.H);
    std::cout<<"py = "<<py<<std::endl;
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        if(iy<py)
        {
            Vec3_t vtemp(0, 0, 0);
            // Vec3_t vtemp((double) dat.vb, 0, 0);
            dom.Rho[ix][iy][0] = 1.0;
            dom.Vel[ix][iy][0] = vtemp;
            dom.BForce[ix][iy][0] = dat.g0;
            for(size_t k=0; k<dom.Nneigh; ++k)
            {
                dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
                dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            }
        }else{
            double yy = (double) iy;
            double Y = yy-py;
            // double uy = dat.g/(2.0*dat.nu)*(H*(yy-py) - (yy-py)*(yy-py)); 
            double uy = P*(dat.H*Y - 0.5*Y*Y); 
            Vec3_t vtemp(uy, 0, 0);
            // Vec3_t vtemp((double) dat.vb, 0, 0);
            dom.Rho[ix][iy][0] = 1.0;
            dom.Vel[ix][iy][0] = vtemp;
            dom.BForce[ix][iy][0] = dat.g0;
            for(size_t k=0; k<dom.Nneigh; ++k)
            {
                dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
                dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            }
        }
        

    }
}



double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 4;
    size_t Dn = 20;
    size_t Npx = 3;
    size_t Npy = 3;
    size_t H = 5*Dn;
    size_t nx = std::ceil(2.5*Dn*3);
    size_t ny = H + std::ceil(5.5*Dn);
    std::cout<<"nx = "<<nx<<" ny = "<<ny<<" H = "<<H<<std::endl;
    double nu = 5e-4;
    double ratio_l = (double) Dn/0.01; //Ll/Lr
    double ratio_nu = nu/1e-6; // nul/nur
    double ratio_t = ratio_l*ratio_l/ratio_nu; // Tl/Tr
    
    std::cout<<"ratio_l = "<<ratio_l<<" ratio_t = "<<ratio_t<<std::endl;
    
    double Re = 1.437e4;
    double vmax = nu*Re/H*1.5;
    double gx = 19.6e-3*ratio_l/(ratio_t*ratio_t);
    std::cout<<"vmax = "<<vmax<<" gx = "<<gx<<std::endl;

    double gapx = std::ceil(1.5*Dn);
    double gapy = 1*Dn;
    std::cout<<"gapx = "<<gapx<<std::endl;
    std::cout<<"gapy = "<<gapy<<std::endl;

    if(argc>=2) Nproc = atoi(argv[1]);
    
    double R = (double) Dn * 0.5;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    
    double Ga = 0.0;
    double rho = 1.0;
    double rhos = 2.0;
    double gy = 0;
    std::cout<<"gy = "<<gy<<std::endl;

    std::cout<<"R = "<<R<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;

    my_dat.y_start = (double) (ny - H);
    Vec3_t g0(gx,0.0,0.0);
    my_dat.g0 = g0;
    my_dat.vmax = vmax;
    my_dat.H = (double) H;
    dom.Nproc = Nproc;       

    // //initial
    
    
    
    Vec3_t pos(0,0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    double rand[3] = {0, -1, 1};
    //DEM
    dom.dtdem = 1*dt;
    for(size_t ipx=0; ipx<Npx; ++ipx)
    {
        pos = 0.5*gapx+R + (gapx+R*2)*ipx, gapy*0.5+R + (gapy+R*2)*0 + rand[0], 0;
        dom.Particles.push_back(DEM::Disk(0, pos, v, w, rhos, R, dom.dtdem));

        pos = 0.5*gapx+R + (gapx+R*2)*ipx, gapy*0.5+R + (gapy+R*2)*1 + rand[1], 0;
        dom.Particles.push_back(DEM::Disk(0, pos, v, w, rhos, R, dom.dtdem));

        pos = 0.5*gapx+R + (gapx+R*2)*ipx, gapy*0.5+R + (gapy+R*2)*2 + rand[2], 0;
        dom.Particles.push_back(DEM::Disk(0, pos, v, w, rhos, R, dom.dtdem));

    }

  

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, -M_PI*R*R*(rhos/rho-1)*gy, 0.0;
        dom.Particles[ip].Kn = 0;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = R;
        dom.Particles[ip].FixVeloc();

    }
    
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        // dom.IsSolid[ix][ny-1][0] = true;
    }

    // Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    
    // // dom.Initial(rho,v0,g0);
    // // dom.InitialFromH5("test_pbed1_0017.h5",g0);
    Initial(dom, dom.UserData);
    // // dom.InitialFromH5("test_pbed4_0113.h5",g0);


    double Tf = 2;
    double dtout = 1;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    dom.Sc = 0.17； 
    //solving
    dom.SolveP( Tf, dtout, "test_pbed1", Setup, NULL);
    
    return 0;
}MECHSYS_CATCH
