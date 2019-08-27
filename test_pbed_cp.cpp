#include "./lbm/Domain.h"
#include <time.h>

struct myUserData
{
    std::ofstream oss_ss;
    double nu;
    double R;
    double rhos;
    Vec3_t g0;
    double V;
    double K;
    double Tf;
    std::vector<double> Y;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    // double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","permeability");
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"U"<<Util::_8s<<"K\n";
    }else{
        double U = 0;
        double num = 0;
        for(size_t ix=0;ix<nx;++ix)
        for(size_t iy=0;iy<ny;++iy)
        for(size_t iz=0;iz<nz;++iz)
        {
            if(dom.Gamma[ix][iy][iz]>1e-9) continue;
            num += 1.0;
            U += dot(dom.Vel[ix][iy][iz],dat.g0/norm(dat.g0));
        }
        std::cout<<dat.g0/norm(dat.g0)<<std::endl;
        U /= num;
        double K = (U*dat.nu)/norm(dat.g0);

        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<U<<Util::_8s<<K<<std::endl;
        
        double r = std::fabs(K-dat.K)/dat.K;

        if(r<1e-5||U<0||std::fabs(U)>1.0)
        {
            dom.Time = dat.Tf;
        }
        dat.K = K;
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

    // double py = dat.Ny*2*dat.R + (dat.Ny-1)*dat.gap;
    double py = *std::max_element(dat.Y.begin(),dat.Y.end())+dat.R;
    double H = ((double) ny-1) - py;
    std::cout<<"py = "<<py<<" H = "<<H<<std::endl;
    std::cout<<"max vel "<<norm(dat.g0)/(2.0*dat.nu)*H*H/4.0<<std::endl;
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
            double uy = norm(dat.g0)/(2.0*dat.nu)*(-1.0*((Y-H)*(Y-H)-H*H)); 
            Vec3_t vtemp(uy, 0, 0);
            // Vec3_t vtemp((double) dat.vb, 0, 0);
            dom.Rho[ix][iy][0] = 1.0;
            dom.Vel[ix][iy][0] = vtemp;
            dom.BForce[ix][iy][0] = dat.g0;
            dat.V += uy;
            for(size_t k=0; k<dom.Nneigh; ++k)
            {
                dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
                dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            }
        }
        

    }
    dat.V /= (nx*ny);
    std::cout<<dat.V<<std::endl;
}



double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 5;
    size_t Rn = 10;
    double nu = 0.1;
    std::cout<<"nu "<<nu<<std::endl;

    if(argc>=2) Nproc = atoi(argv[1]);
    
    size_t nx = 250;//caution may vary
    size_t ny = 250+1;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = (double) Rn;
    double rho = 1.0;
    double rhos = 2.0;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.R = R;
    Vec3_t g0(0.0,-1e-8,0.0);
    my_dat.g0 = g0;
    my_dat.V = 0;
    std::cout<<"g = "<<g0<<std::endl;
    dom.Nproc = Nproc;       

    //initial
    
    my_dat.rhos = rhos;
    
    
    Vec3_t pos(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.01*dt;
    //fixed
    char const *infilename = "circle_mono.txt";
    String fn;
    fn.Printf("%s",infilename);
    std::fstream ifile(fn.CStr(),std::ios::in);
    int N = 0;
    double fi = 0;
    if(!ifile.fail()){
        ifile>>fi;
        ifile>>N;
        std::cout<<"porosity "<<fi<<"Number of Circle "<<N<<std::endl;
        for(size_t i=0; i<N; ++i)
        {
            double x;
            double y;
            double r;
            ifile>>x>>y>>r;
            pos = x,y,0;
            dom.Particles.push_back(DEM::Disk(i, pos, v, w, rhos, r, dom.dtdem));
            dom.Particles[i].FixVeloc();
            my_dat.Y.push_back(y);
        }
    }else{
        std::cout<<"READ TXT ERRO!!!!"<<std::endl;
    }


    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, 0.0, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = -0.3;
        dom.Particles[ip].Kt = 2.5;
        dom.Particles[ip].Mu = 0.4;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
        dom.Particles[ip].FixVeloc();

    }
    
    

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    
    dom.Initial(rho,v0,g0);
    // Initial(dom, dom.UserData);
    // dom.InitialFromH5("test_pbed_0999.h5",g0);


    double Tf = 3;
    my_dat.Tf = Tf;
    double dtout = 1;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    dom.IsOutH5 = false;
    //solving
    dom.SolveP( Tf, dtout, "test_pbed_cp", NULL, Report);
    
    return 0;
}MECHSYS_CATCH
