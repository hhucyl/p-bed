#include <iostream>
#include "../lbm/Domain.h"


int main()
{
    char const * TheFileKey = "/home/user/p-bed/test_pbed_r2";
    char const * FileKey = "turbulence";
    String fn;
    size_t fileNum = 1000;
    size_t Nneigh = 9;
    double nu = 1e-3;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(1,1,1),1.0,1.0);
    std::vector<double *> VV;
    std::vector<double *> Vvt;
    std::vector<double *> GGa;
    int nx,ny,nz;

    for(size_t i=0;i<fileNum; ++i)
    {
        fn.Printf("%s_%04d", TheFileKey, i);
        fn.append(".h5");
        // std::cout<<fn.CStr()<<std::endl;
        hid_t file_id = H5Fopen(fn.CStr(),H5F_ACC_RDONLY,H5P_DEFAULT);
        int N[1];
        H5LTread_dataset_int(file_id,"/Nx",N);
        nx = N[0];
        H5LTread_dataset_int(file_id,"/Ny",N);
        ny = N[0];
        H5LTread_dataset_int(file_id,"/Nz",N);
        nz = N[0];
        double *Ff = new double[Nneigh*nx*ny*nz];
        double *Ga = new double[nx*ny*nz];
        double *vt = new double[nx*ny*nz];
        double *Vvel = new double[3*nx*ny*nz];
        H5LTread_dataset_double(file_id,"/F_0",Ff);
        H5LTread_dataset_double(file_id,"/Gamma",Ga);
        size_t nn=0;
        for(size_t ix=0; ix<nx; ++ix)
        for(size_t iy=0; iy<ny; ++iy)
        for(size_t iz=0; iz<nz; ++iz)
        {
            Vec3_t vel(0.0,0.0,0.0);
            double rho = 0.0;
            double Gamma  = Ga[nn];
		    for(size_t k=0; k<Nneigh; k++)
		    {
			
			    double f = Ff[Nneigh*nn + k];
                rho += f;
                vel += f*dom.C[k];
		    }
            vel *= dom.Cs/rho;
            Vvel[3*nn] = vel(0); 
            Vvel[3*nn+1] = vel(1); 
            Vvel[3*nn+2] = vel(2); 
            if(Gamma>1.0)
            {
                vel = OrthoSys::O;
                rho = 1.0;
            }
            double tau = dom.Tau;
            double VdotV = dot(vel,vel);
            double NonEq[Nneigh];
            double Q = 0.0;
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,dom.C[k]);
                double Feq   = dom.W[k]*rho*(1.0 + 3.0*VdotC + 4.5*VdotC*VdotC - 1.5*VdotV);
                double f = Ff[Nneigh*nn + k];
                NonEq[k] = f - Feq;
                Q +=  NonEq[k]*NonEq[k]*dom.EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*dom.Sc/rho));
            vt[nn] = 1.0/3.0*(tau - dom.Tau);
            nn++;
        }
        
        VV.push_back(Vvel);
        Vvt.push_back(vt);
        GGa.push_back(Ga);

        



        delete[] Ff;
        H5Fflush(file_id,H5F_SCOPE_GLOBAL);
        H5Fclose(file_id);

        
        

    
    }
    
    double *Va = new double[3*nx*ny*nz];
    memset(Va,0,3*nx*ny*nz);
    for(size_t i=0;i<fileNum; ++i)
    {
        for(size_t ii=0; ii<3*nx*ny*nz; ++ii)
        {
            Va[ii] += VV[i][ii];
        }

    }
    for(size_t ii=0; ii<3*nx*ny*nz; ++ii)
    {
        Va[ii] /= fileNum;
    }
    // std::cout<<1<<std::endl;
    for(size_t i=0;i<fileNum; ++i)
    {

        double *Vhas = new double[3*nx*ny*nz];
        for(size_t ii=0; ii<3*nx*ny*nz; ++ii)
        {
            Vhas[ii] = VV[i][ii]-Va[ii];
        }

        String fn1,fn2;
        fn1.Printf("%s_%04d", FileKey, i);
        fn2.Printf("%s_%04d", FileKey, i);
        fn1.append(".h5");

        hid_t     file_id1;
        file_id1 = H5Fcreate(fn1.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hsize_t dims[1];
        String dsname;

        int N[1];
        dims[0] = 1;
        N[0] = nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id1,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id1,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id1,dsname.CStr(),1,dims,N);
        
        dims[0] = nx*ny*nz;
        dsname.Printf("vt");
        H5LTmake_dataset_double(file_id1,dsname.CStr(),1,dims,Vvt[i]);
        dsname.Printf("Gamma");
        H5LTmake_dataset_double(file_id1,dsname.CStr(),1,dims,GGa[i]);
        
        dims[0] = 3*nx*ny*nz;
        dsname.Printf("Vhas");
        H5LTmake_dataset_double(file_id1,dsname.CStr(),1,dims,Vhas);
        dsname.Printf("Vave");
        H5LTmake_dataset_double(file_id1,dsname.CStr(),1,dims,Va);
       
        delete[] Vhas;
        delete[] VV[i];
        delete[] Vvt[i];
        delete[] GGa[i];
        H5Fflush(file_id1,H5F_SCOPE_GLOBAL);
        H5Fclose(file_id1);

        std::ostringstream oss;


   
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ny << " " << nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn1.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"vt\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn1.CStr() <<":/vt\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Vhas"<< "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn1.CStr() <<":/Vhas" << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Vave"<< "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn1.CStr() <<":/Vave" << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
        
        fn2.append(".xmf");
        std::ofstream of(fn2.CStr(), std::ios::out);
        of << oss.str();
        of.close();


        

    }
    delete[] Va;

}
