#ifndef RW_PARTICLE_H
#define RW_PARTICLE_H

#include <mechsys/linalg/matvec.h>
#include<random>

namespace RW
{
class Particle
{
    public:
    Particle(Vec3_t &X, int d, double dm);
    int Tag;
    double Time;
    double Tlimt;

    Vec3_t X;
    Vec3_t Xb;
    Vec3_t Xb1;
    int D;
    double Dm;
    std::vector<double> TimeArray;
    std::vector<double> OutTimeArray;
    void Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt);
    void Leave(int modexy, Vec3_t &Box);
    void Leave1(int modexy, Vec3_t &Box);
    void Leave2(int modexy, Vec3_t &Box);//leave and set as -200
    void LeaveReflect(int modexy1, Vec3_t &Box1);
    void Reflect(Vec3_t &C, double R, double Time);
};
inline Particle::Particle(Vec3_t &X0, int d, double dm)
{
    X = X0;
    Xb = X;
    Xb1 = X;

    D = d;
    Dm = dm;
    Time = 0;
    Tlimt = 1e4;
    Tag = 1;
}

inline void Particle::Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normal(0,1);
    Vec3_t e(normal(gen),normal(gen), 0.0);
    e /= Norm(e);
    // std::cout<<e<<std::endl;
    //VV 0=>11 1=>21 2=>12 3=>22
    //*12   *22
    //*11   *21
    //idx 0=>x1 1=>x2 2=>y1 3=>y2
    double x1 = (double) idx[0];
    double x2 = (double) idx[1];
    double y1 = (double) idx[2];
    double y2 = (double) idx[3];
    double x = X(0);
    double y = X(1);
    Vec3_t V(0,0,0);
    V = 1.0/((x2-x1)*(y2-y1))*(VV[0]*(x2-x)*(y2-y)+VV[1]*(x-x1)*(y2-y)+VV[2]*(x2-x)*(y-y1)+VV[3]*(x-x1)*(y-y1));
    // std::cout<<VV[0]<<" "<<VV[1]<<" "<<VV[2]<<" "<<VV[3]<<std::endl;
    // std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
    // std::cout<<V(0)<<"    "<<V(1)<<std::endl;
    Xb1 = Xb;
    Xb = X;
    X = X+V*dt + std::sqrt(2*D*Dm*dt)*e;    
}

inline void Particle::Leave(int modexy, Vec3_t &Box)
{
    if(1000*X(modexy)<1000*Box(0))
    {   
        double dist = std::fabs(X(modexy)-Box(0));
        X(modexy) = Box(1)-dist;
        double distb = Xb(modexy)-Box(0);
        Xb(modexy) = Box(1)+distb;

        // double dist = std::fabs(X(modexy)-Box(0));
        // X(modexy) = Box(1)-dist+1.0;
        // double distb = Xb(modexy)-Box(0);
        // Xb(modexy) = Box(1)+distb+1.0;   

    }
    if(X(modexy)>Box(1))
    {
        double dist = std::fabs(X(modexy)-Box(1));
        X(modexy) = Box(0)+dist;
        double distb = Xb(modexy)-Box(1);
        Xb(modexy) = Box(0)+distb;

        // double dist = std::fabs(X(modexy)-Box(1));
        // X(modexy) = Box(0)+dist-1.0;
        // double distb = Xb(modexy)-Box(1);
        // Xb(modexy) = Box(0)+distb-1.0;

    }
    


}

inline void Particle::Leave1(int modexy, Vec3_t &Box)
{
    if(1000*X(modexy)<1000*Box(0))
    {   
        // double dist = std::fabs(X(modexy)-Box(0));
        // X(modexy) = Box(1)-dist;
        // double distb = Xb(modexy)-Box(0);
        // Xb(modexy) = Box(1)+distb;

        double dist = std::fabs(X(modexy)-Box(0));
        X(modexy) = Box(1)-dist+1.0;
        double distb = Xb(modexy)-Box(0);
        Xb(modexy) = Box(1)+distb+1.0;   

    }
    if(X(modexy)>Box(1))
    {
        // double dist = std::fabs(X(modexy)-Box(1));
        // X(modexy) = Box(0)+dist;
        // double distb = Xb(modexy)-Box(1);
        // Xb(modexy) = Box(0)+distb;

        double dist = std::fabs(X(modexy)-Box(1));
        X(modexy) = Box(0)+dist-1.0;
        double distb = Xb(modexy)-Box(1);
        Xb(modexy) = Box(0)+distb-1.0;

    }
    


}

inline void Particle::Leave2(int modexy, Vec3_t &Box)
{
    if(1000*X(modexy)<1000*Box(0))
    {   
        X(modexy) = -200;
        Xb(modexy) = -200;   

    }
    if(X(modexy)>Box(1))
    {

        X(modexy) = -200;
        Xb(modexy) = -200;

    }
    


}

inline void Particle::LeaveReflect(int modexy1, Vec3_t &Box1)
{
    if(X(modexy1)<Box1(0))
    {
        X(modexy1) = Box1(0) + std::fabs(X(modexy1) - Box1(0));
    }
    if(X(modexy1)>Box1(1))
    {
        X(modexy1) = Box1(1) - std::fabs(X(modexy1) - Box1(1));
    }
}

inline void Particle::Reflect(Vec3_t &C, double R, double Time)
{
    double L1 = Norm(X-Xb);
    double L2 = Norm(Xb-C);
    // std::cout<<"L1 = "<<L1<<" L2 ="<<L2<<std::endl;
    double AA = L1*L1;
    double BB = -2.0*((Xb(1)-C(1))*(Xb(1)-X(1))+(Xb(0)-C(0))*(Xb(0)-X(0)));
    double CC = L2*L2 - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    // std::cout<<"Delta = "<<Delta<<std::endl;
    
    if(Delta>0)
    {
        double q = 0;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0 && q1-1<1e-9;
        bool flag2 = q2>=0 && q2-1<1e-9;
        bool flag = true;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                flag = false;
                std::cout<<q1<<" "<<q2<<" RW Now "<<X<<" RW Pre "<<Xb<<" Circle Center "<<C<<" ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
                std::cout<<Tag<<" Time "<<Time<<" L "<<Norm(X-C)<<" Lb "<<Norm(Xb-C)<<" Lb1 "<<Norm(Xb1-C)<<" Circle R "<<R<<" ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        if(flag)
        {
            Vec3_t Xi(Xb(0)-q*(Xb(0)-X(0)),Xb(1)-q*(Xb(1)-X(1)),0.0);
            // std::cout<<"D = "<<Xi<<std::endl;
            double y = -((X(0)-2.*Xi(0))*(C(0)-Xi(0))*(C(1)-Xi(1)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1))*(C(1)-Xi(1)) + X(0)*(C(0)-Xi(0))*(C(1)-Xi(1)) - X(1)*(C(0)-Xi(0))*(C(0)-Xi(0)))/(R*R);
            double x;
            if(std::fabs(C(1)-Xi(1))>0)
            {
                x = ((C(0)-Xi(0))*(y-X(1)))/(C(1)-Xi(1)) + X(0);
            }else{
                x = -((C(1)-Xi(1))*y + (X(0)-2.*Xi(0))*(C(0)-Xi(0)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1)))/(C(0)-Xi(0));
            }
            X = x,y,0.0;
        }
        
        
        
    }else{
        std::cout<<Delta<<std::endl;
        std::cout<<" RW Now "<<X<<" RW Pre "<<Xb<<" Circle Center "<<C<<" ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
        std::cout<<Tag<<" Time "<<Time<<" L "<<Norm(X-C)<<" Lb "<<Norm(Xb-C)<<" Lb1 "<<Norm(Xb1-C)<<" Circle R "<<R<<" ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
    }
    
    
}


}


#endif