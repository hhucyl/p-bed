#ifndef LBM_SOLVE_H
#define LBM_SOLVE_H

inline void Domain::Solve(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveP(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            if(IsOutH5)
            {
                String fn;
                fn.Printf("%s_%04d", TheFileKey, idx_out);
            
                WriteXDMF(fn.CStr());
            }
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        


        
         
        //set fluid force
        if(std::fabs(Time)<1e-6 && (!IsContinue)){
            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
            AddDisksG();
        }

        //collide and streaming
        // CollideSRTGamma();
        CollideTRTGamma();
        Stream();
        BounceBack(false);
        CalcProps();

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveIBM(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
            
        
         
            //set fluid force
            if(flag){
                AddDisksIBM();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
                // UpdateParticlesContactsVL();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
            // GhostParticles.clear();
            // GhostParticles.assign(Particles.begin(), Particles.end());
            // GhostPeriodic();

            //test for friction
            // for(size_t ip=0; ip<Particles.size(); ++ip)
            // {
            //     if(Particles[ip].IsFree()) continue;

            //     Particles[ip].X = Particles[ip+1].X;
            //     Particles[ip].X(1) = Particles[ip].Xb(1);
            //     Particles[ip].V = 0.0,0.0,0.0; 
            // }

        }
        

        //collide and streaming
        CollideMRTIBM();
        Stream();
        BounceBack(false);
        CalcProps();
        // std::cout<<std::boolalpha<<GhostParticles[0].Ghost<<std::endl;
        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolvePRW(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        


        
         
        //set fluid force
        if(std::fabs(Time)<1e-6 && (!IsContinue)){
            AddDisksG();
        }
        if(std::fabs(Time)<1e-6){
            ApplyDisksCheck();
        }

        //collide and streaming
        // CollideSRTGamma();
        CollideTRTGamma();
        Stream();
        BounceBack(false);
        CalcProps();

        //trace particle
        rwsolve_sub(dt);

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveRW(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        //trace particle
        rwsolve_sub(dt);

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::rwsolve_sub1(double dt)
{
    // std::cout<<1<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<Ndim(0); ++ix)
    for(size_t iy=0; iy<Ndim(1); ++iy)
    {
        Con[ix][iy][0] = 0;
    }


    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<(int) RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        int x1 = std::floor(RWP->X(0));
        int y1 = std::floor(RWP->X(1));
        int x2 = x1+1;
        int y2 = y1+1;
        if(x1<0 || y1<0 || x2>(int) Ndim(0)-1 || y2>(int) Ndim(1)-1)
        {   
            if(y2>(int) Ndim(1)-1)
            {
                RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
            }else{
                RWP->X = -1000,-1000,0;
            }
        }else{
            std::vector<Vec3_t> VV{Vel[x1][y1][0],Vel[x2][y1][0],Vel[x1][y2][0],Vel[x2][y2][0]};
            std::vector<int> idx{x1,x2,y1,y2};
            RWP->Move(VV,idx,dt);
            if(x1<0 || y1<0 || x2>(int) Ndim(0)-1 || y2>(int) Ndim(1)-1)
            {   
                if(y2>(int) Ndim(1)-1)
                {
                    RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
                }else{
                    RWP->X = -1000,-1000,0;
                }
            }
            if(RWP->X(0)<0) continue;
            int ix = std::round(RWP->X(0));
            int iy = std::round(RWP->X(1));
            if(Gamma[ix][iy][0]>1e-9 && Check[ix][iy][0]>0) 
            {
                int ip = Check[ix][iy][0];
                if(Norm(Particles[ip].X-RWP->X)<=Particles[ip].Rh && Norm(Particles[ip].Xb-RWP->Xb)>Particles[ip].Rh)
                {
                    RWP->Reflect(Particles[ip].X,Particles[ip].Rh,Time);
                    if(x1<0 || y1<0 || x2>(int) Ndim(0)-1 || y2>(int) Ndim(1)-1)
                    {   
                        if(y2>(int) Ndim(1)-1)
                        {
                            RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
                        }else{
                            RWP->X = -1000,-1000,0;
                        }
                    }

                }
            }

            ix = std::round(RWP->X(0));
            iy = std::round(RWP->X(1));
            #pragma omp atomic
                Con[ix][iy][0] += 1;
        }
        
    }

}


inline void Domain::rwsolve_sub(double dt)
{
    // std::cout<<1<<std::endl;
    if(Time<1) std::cout<<"--- rwsub ---"<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<Ndim(0); ++ix)
    for(size_t iy=0; iy<Ndim(1); ++iy)
    {
        Con[ix][iy][0] = 0;
    }

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<(int) RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        int x1 = std::floor(RWP->X(0));
        int y1 = std::floor(RWP->X(1));
        int x2 = x1+1;
        int y2 = y1+1;
        
        // if(x1<0 || y1<0 || x2>(int) Ndim(0)-1 || y2>(int) Ndim(1)-1)
        // {
        //     RWP->Leave(modexy,Box);
        //     if(std::floor(RWP->X(1))+1>Ndim(1)-1)
        //     {
        //         RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
        //     }
        //     if(std::floor(RWP->X(1))<0)
        //     {
        //         RWP->X(1) = - RWP->X(1);
        //     }

        //     x1 = std::floor(RWP->X(0));
        //     y1 = std::floor(RWP->X(1));
        //     x2 = x1+1;
        //     y2 = y1+1;
        // }
        std::vector<Vec3_t> VV{Vel[x1][y1][0],Vel[x2][y1][0],Vel[x1][y2][0],Vel[x2][y2][0]};
        std::vector<int> idx{x1,x2,y1,y2};
        RWP->Move(VV,idx,dt);
        RWP->Leave(modexy,Box);
        if(std::floor(RWP->X(1))+1>(int) Ndim(1)-1)
        {
            RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
        }
        if(std::floor(RWP->X(1))<0)
        {
            RWP->X(1) = - RWP->X(1);
        }
        int ix = std::round(RWP->X(0));
        int iy = std::round(RWP->X(1));
        if(Gamma[ix][iy][0]>1e-9 && Check[ix][iy][0]>-1e-12) 
        {
            int ip = Check[ix][iy][0];
            // if(Norm(Particles[ip].X-RWP->X)<=Particles[ip].Rh && Norm(Particles[ip].X-RWP->Xb)>Particles[ip].Rh)
            if(Norm(Particles[ip].X-RWP->X) < Particles[ip].Rh)
            {
                if(Norm(Particles[ip].X-RWP->Xb)>Particles[ip].Rh)
                {
                    RWP->Reflect(Particles[ip].X,Particles[ip].Rh,Time);
                }else{
                    RWP->X = RWP->Xb;
                }
                // std::cout<<2<<std::endl;
                RWP->Leave(modexy,Box);
                if(std::floor(RWP->X(1))+1>(int) Ndim(1)-1)
                {
                    RWP->X(1) = 2*(Ndim(1)-1) - RWP->X(1);
                }
                if(std::floor(RWP->X(1))<0)
                {
                    RWP->X(1) = - RWP->X(1);
                }
            }
            if(Norm(Particles[ip].X-RWP->X) < Particles[ip].Rh)
            {
                RWP->X = RWP->Xb;
            }

        }
        ix = std::round(RWP->X(0));
        iy = std::round(RWP->X(1));
        #pragma omp atomic
            Con[ix][iy][0] += 1;
    }

}

inline void Domain::CheckInside()
{
    for(size_t i=0; i<RWParticles.size(); ++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        for(size_t ip=0; ip<Particles.size(); ++ip)
        {
            DEM::Disk *Pa = &Particles[ip];
            if(Norm(Pa->X-RWP->X)<Pa->Rh)
            {
                RWParticles.erase(RWParticles.begin()+i);
                break;
            }
        }

    }
}

#endif