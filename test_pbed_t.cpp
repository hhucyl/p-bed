#include "./lbm/Domain.h"
#include <time.h>

struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    int Ny;
    double gap;
    double rhos;
    Vec3_t g0;
    double V;
    std::vector<double> Y;
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
    double ppy = *std::max_element(dat.Y.begin(),dat.Y.end())+dat.R;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for(int i=0;i<(int) dom.RWParticles.size();++i)
    {
        RW::Particle *RWP = &dom.RWParticles[i]; 
        // double px = RWP->X(0);
        double py = RWP->X(1);
        if(py<ppy)
        {
            RWP->Time += 1;
            RWP->Tag = -1;
        }else{
            if(RWP->Time>RWP->Tlimt){
                RWP->TimeArray.push_back(RWP->Time);
                RWP->OutTimeArray.push_back(dom.Time);
            }
            RWP->Tag = 1;
            RWP->Time = 0;
        }
    }

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
    std::cout<<"max vel "<<dat.g/(2.0*dat.nu)*H*H/4.0<<std::endl;
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
            double uy = dat.g/(2.0*dat.nu)*(-0.25*((Y-H)*(Y-H)-H*H)); 
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
    size_t Nproc = 12;
    int Nx = 10;
    int Ny = 5;
    size_t Rn = 10;
    double Re = 1e4;
    size_t H = 50;
    double vmax = 0.15;
    double nu = 2.0/3.0*vmax*H/Re;
    std::cout<<"nu "<<nu<<std::endl;

    double gap = 2;
    double mag = 1;
    double Dm0 = 2.3e-9;
    double Dm = nu/1e-6*Dm0;
    std::cout<<"Dm0 = "<<Dm0<<" Dm = "<<Dm<<std::endl;
    //rwP
    int RWP = 20;
    if(argc>=2) Nproc = atoi(argv[0]);
    
    int gapn = std::ceil(gap*Nx);
    int gapny = std::ceil(gap*Ny);
    std::cout<<"extra gap n "<<gapn<<std::endl;
    size_t nx = 2*Rn*Nx+gapn;
    size_t ny = 2*Rn*Ny+gapny+H;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = (double) Rn;
    double Ga = 0.0;
    double rho = 1.0;
    double rhos = 2.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;

    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 8.0*nu*vmax/((double)H*(double)H);
    my_dat.R = R;
    my_dat.gap = gap;
    my_dat.Ny = Ny;
    Vec3_t g0(my_dat.g,0.0,0.0);
    my_dat.g0 = g0;
    my_dat.V = 0;
    std::cout<<"gx = "<<my_dat.g<<std::endl;
    dom.Nproc = Nproc;       

    //initial
    
    my_dat.rhos = rhos;
    
    
    Vec3_t pos(R+gap,R,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.01*dt;
    int pnum = 0;
    //fixed
    for(int ip=0; ip<Nx; ++ip)
    {
        // std::cout<<pos<<std::endl;
        dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, R, dom.dtdem));
        dom.Particles[ip].FixVeloc();
        dxp = 2.0*R+gap,0.0,0.0;
        pos = pos+dxp;
        pnum++;
        my_dat.Y.push_back(R);
    }
    //move
    for(int ipy=0; ipy<Ny-1; ++ipy)
    {
        pos = R+gap,(2*ipy+3)*R + (ipy+1)*gap,0.0;
        for(int ipx=0; ipx<Nx; ++ipx)
        {
            Vec3_t dxr(random(-mag,mag),random(-mag,mag),0.0);
            // Vec3_t dxr(0.0,0.0,0.0);
            dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            dxp = 2.0*R+gap,0.0,0.0;
            pos = pos+dxp;
            pnum++;
        }
        my_dat.Y.push_back(pos(1));

    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, -M_PI*R*R*(rhos/rho-1)*gy, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = -0.3;
        dom.Particles[ip].Kt = 2.5;
        dom.Particles[ip].Mu = 0.4;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
        dom.Particles[ip].FixVeloc();

    }
    
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        // dom.IsSolid[ix][ny-1][0] = true;
    }

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    
    // dom.Initial(rho,v0,g0);
    Initial(dom, dom.UserData);
    dom.InitialFromH5("test_pbed2_0999.h5",g0);

    //RWParticles
    double py = *std::max_element(my_dat.Y.begin(),my_dat.Y.end())+my_dat.R;
    int starty = std::ceil(py);
    // std::vector<int> startx{22,65,109,154,199};
    // for(int i=0; i<startx.size(); ++i)
    for(int i=0; i<nx-1; ++i)
    for(int iy=starty+1; iy<starty+21; ++iy)
    {
        // int x1 = startx[i];
        int x1 = i;
        int y1 = iy;
        int x2 = x1+1;
        int y2 = y1+1;
        for(int ip=0; ip<RWP; ++ip)
        {
            Vec3_t xt(random(x1,x2),random(y1,y2),0);
            int xx = std::round(xt(0));
            int yy = std::round(xt(1));
            if(dom.Gamma[xx][yy][0]<1e-9)
                dom.RWParticles.push_back(RW::Particle(xt,2,Dm));
        } 
    }
    for(int i=0;i<(int) dom.RWParticles.size();++i)
    {
        dom.RWParticles[i].Tag = 1;
        dom.RWParticles[i].Tlimt = 1e3;
    }
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;

    


    double Tf = 1e6;
    double dtout = 1e3;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolvePRW( Tf, dtout, "test_pbed_t2", Setup, NULL);
    

    //last time to caculate always in
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
    for(int i=0;i<(int) dom.RWParticles.size();++i)
    {
        RW::Particle *RWP = &dom.RWParticles[i]; 
        // double px = RWP->X(0);
        double rpy = RWP->X(1);
        if(rpy<py)
        {
            RWP->Time += 1;
            RWP->Tag = -1;
            if(RWP->Time>RWP->Tlimt){
                RWP->TimeArray.push_back(RWP->Time);
            }
        }else{
            if(RWP->Time>RWP->Tlimt){
                RWP->TimeArray.push_back(RWP->Time);
                RWP->OutTimeArray.push_back(dom.Time);
            }
            RWP->Tag = 1;
            RWP->Time = 0;
        }
    }

    //output time
    std::ofstream ofile;
    ofile.open("time.txt",std::ios::out);
    for(int i=0;i<(int) dom.RWParticles.size();++i)
    {
        RW::Particle *RWP = &dom.RWParticles[i]; 
        ofile<<Util::_10_6<<i;
        for(size_t it=0; it<RWP->TimeArray.size(); ++it)
        {
            ofile<<Util::_8s<<dom.RWParticles[i].TimeArray[it];
        }
        ofile<<std::endl;
        ofile<<Util::_10_6<<i;
        for(size_t it=0; it<RWP->OutTimeArray.size(); ++it)
        {
            ofile<<Util::_8s<<dom.RWParticles[i].OutTimeArray[it];
        }
        ofile<<std::endl;
    }
    
    return 0;
}MECHSYS_CATCH
