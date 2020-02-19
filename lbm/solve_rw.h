#ifndef LBM_SOLVE_RW_H
#define LBM_SOLVE_RW_H

inline void Domain::rwsolve_sub2(double dt)
{
    // std::cout<<1<<std::endl;
    if(Time<1e-6) std::cout<<"--- rwsolve_sub2 ---"<<std::endl;
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
        std::vector<int> idc{-1,-1};
       

        
        Vec3_t v0(0,0,0); 
        std::vector<Vec3_t> VV{v0,v0,v0,v0};
        std::vector<int> idx{-1,-1,-1,-1};
        // std::cout<<111<<std::endl;
        Pa2GridV(RWP,idx,VV);
        // std::cout<<222<<std::endl;
        RWP->Move(VV,idx,dt);
        RWP->Leave1(modexy,Box);
        RWP->LeaveReflect(modexy1,Box1);

        
        Pa2Grid(RWP,idc);
        int ix = idc[0];
        int iy = idc[1];
        // std::cerr<<RWP->X(0)<<" "<<RWP->X(1)<<" "<<ix<<" "<<iy<<std::endl;
        // if(ix<0 || ix>Ndim(0)-1) std::cout<<"x "<<ix<<std::endl;
        // if(iy<0 || iy>Ndim(1)-1) std::cout<<"y "<<iy<<std::endl;
        if(Check[ix][iy][0]>-0.5) 
        {


            int ip = Check[ix][iy][0];
            DEM::Disk *P = &Particles[ip];
            DEM::Disk *GP = &GhostParticles[ip];
            if(Norm(P->X-RWP->X) < P->Rh)
            {
                if(Norm(P->X-RWP->Xb)>P->Rh)
                {
                    RWP->Reflect(P->X,P->Rh,Time);
                }else{
                    RWP->X = RWP->Xb;
                }
                // std::cout<<2<<std::endl;
                RWP->Leave1(modexy,Box);
                RWP->LeaveReflect(modexy1,Box1);
            }
            if(Norm(GP->X-RWP->X) < GP->Rh)
            {
                if(Norm(GP->X-RWP->Xb)>GP->Rh)
                {
                    RWP->Reflect(GP->X,GP->Rh,Time);
                }else{
                    RWP->X = RWP->Xb;
                }
                // std::cout<<2<<std::endl;
                RWP->Leave1(modexy,Box);
                RWP->LeaveReflect(modexy1,Box1);
            }
            if(Norm(P->X-RWP->X) < P->Rh)
            {
                RWP->X = RWP->Xb;
            }
            
            if(Norm(GP->X-RWP->X) < GP->Rh)
            {
                RWP->X = RWP->Xb;
            }

            
        }
        
        // std::cout<<444<<std::endl;
        Pa2Grid(RWP,idc);
        ix = idc[0];
        iy = idc[1];
        #pragma omp atomic
            Con[ix][iy][0] += 1;

    }

}



inline void Domain::Pa2Grid(RW::Particle *RWP,std::vector<int> &idc)
{
    int xx = std::round(RWP->X(0));
    int yy = std::round(RWP->X(1));
    if(xx<0) xx = (int)Ndim(0) + xx;
    if(xx>(int)Ndim(0)-1) xx = xx - (int)Ndim(0);
    if(yy<0) yy = (int)Ndim(1) + yy;
    if(yy>(int)Ndim(0)-1) yy = yy - (int)Ndim(1);
    idc[0] = xx;
    idc[1] = yy;  
}
inline void Domain::Pa2GridV(RW::Particle *RWP, std::vector<int> &idx, std::vector<Vec3_t> &VV)
{
    //VV 0=>11 1=>21 2=>12 3=>22
    //*12   *22
    //*11   *21
    //idx 0=>x1 1=>x2 2=>y1 3=>y2
    //!! The other dimension of periodic dimension is considered as no flux boundary using specular reflection
    // the relative location is not change only the velocity was changed and idc
    int x1 = std::floor(RWP->X(0));
    int y1 = std::floor(RWP->X(1));
    int x2 = x1+1;
    int y2 = y1+1;
    idx[0] = x1;
    idx[1] = x2;
    idx[2] = y1;
    idx[3] = y2;
    if(modexy ==0)
    {
        if(x1<0)
        {
            int x11 = (int)Ndim(modexy)+x1;
            VV[0] = Vel[x11][y1][0];
            VV[1] = Vel[x2][y1][0];
            VV[2] = Vel[x11][y2][0];
            VV[3] = Vel[x2][y2][0];
        }else{
            if(x2>(int)Ndim(modexy)-1)
            {
                int x22 = x2 - (int)Ndim(modexy);
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x22][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x22][y2][0];
            }else{
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x2][y2][0];
            }
            
        }
        

    }else{
        if(y1<0)
        {
            int y11 = (int)Ndim(modexy)+y1;
            VV[0] = Vel[x1][y11][0];
            VV[1] = Vel[x2][y11][0];
            VV[2] = Vel[x1][y2][0];
            VV[3] = Vel[x2][y2][0];
        }else{
            if(y2>(int)Ndim(modexy)-1)
            {
                int y22 = y2 - (int)Ndim(modexy);
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y22][0];
                VV[3] = Vel[x2][y22][0];
            }else{
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x2][y2][0];
            }
            
        }

    }
}


//OLD CODE



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

#endif