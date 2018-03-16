#include"swfd.h"

#define PI 3.14159
using namespace std;
static void Malloc2D(float **&a, int Nx ,int Nz)
{
	int i,j;
	a=new float *[Nx];
	for(i=0;i<Nx;i++)
	{
		a[i]=new float [Nz];
	}
}
static void Init_0(float **&a, int Nx ,int Nz)
{
	int i,j;
	for (i=0;i<Nx;i++)
	{
		for (j=0;j<Nz;j++)
		{
			a[i][j]=float(0);
		}
	}
}
static void CreatResultFolder(const char *Dir ) 
{
	if (_access(Dir, 0) == 0)  
    {  
        cout<<Dir<<" exists"<<endl;  
		/*cout<<"Are you sure to cover the file in "<<Dir<<" by the following result data?"<<endl;
		system("pause");*/
    }  
	if( (_access( Dir, 0 )) == -1 ) 
	{
		int flag = _mkdir(Dir);
		if (flag==0)
			cout<<Dir<<" Creat directory successfully"<<endl;
		else
			cout<<Dir<<" Fail to creat directory!"<<endl;
	}
}
static void outputfile(float **a,int start_X,int end_X, int start_Z,int end_Z, string filename)//文件输出函数
{	
	int i,j;
	ofstream outf(filename,ios::binary);
	if (!outf)
	{
		cout<<"cannot open output file"<<endl;
		abort();
	}
	for(i=start_X;i<end_X;i++)// i表示列，即按列输出
	{
		for(j=start_Z;j<end_Z;j++)
			{outf.write((char*)&a[i][j],sizeof(a[i][j]));
		/*cout<<j<<endl;
		system("pause");*/}
	}
	outf.close();
}
static float Gausswavelet(float t)
	{
		float wavelet,tmp;
		float pai,f0;
		f0=30;
		pai=(float)PI;
		tmp=pai*f0*t;
		wavelet=(float)(-2000.0*pai*pai*f0*f0*t*exp(-1*tmp*tmp)*(3-2*tmp*tmp));
		return wavelet;
		
	}
static float d_pml(float x,float Vp,float Dx, int N_pml)
	{
		float pi=3.14159;
		float attenuation,R,pml_width;
		R=pow(10.,-6);
		pml_width=Dx*N_pml;
		attenuation=4.*(float)(3*Vp*log10(1.0/R)*pow(x/pml_width,2)/(2*pml_width));
		return attenuation;
	}
void SetGridParmeters(FD2DP & a,float dt,float dx,float dz,int nx,int nz,int nt,int pml,\
				   char *FileName_snap,char *FileName_seismicrecord,char *FilName_model)
{
	a.dt=dt,a.dx=dx,a.dz=dz;    //网格步长和时间步长
	a.nx=nx;a.nz=nz;a.nt=nt;    //模型区域大小和时间记录长度	
	a.pml=pml;                  //吸收边界网格点数
	//a.isx=isx; a.isz=isz;       //震源位置
	a.FileName_snap=FileName_snap;
	a.FileName_seismicrecord=FileName_seismicrecord;
	a.FilName_model=FilName_model;
}
void SetModelParmeters(FD2DP &a , int createmodelway,vector<float>&Source_X,vector<float>&Source_Z,\
					   vector<float>&h,vector<float>&Vp,vector<float>&Vs,vector<float>&rou)
{
	//a.n_source=n_source;
	a.Source_X=Source_X;
	a.Source_Z=Source_Z;
	a.createmodelway=createmodelway;
	a.h=h;
	a.Vp=Vp;
	a.Vs=Vs;
	a.rou=rou;
}
FD2D::FD2D_openmp(const FD2DP &par)
{
	//NN=8;                    //差分阶数
	dt=par.dt, dx=par.dx, dz=par.dz; //网格步长和时间步长
	nx=par.nx, nz=par.nz, nt=par.nt;      //模型区域大小和时间记录长度	
	pml=par.pml;                      //吸收边界网格点数
	n_pml=pml+NN/2;          //实际吸收边界网格点数
	nx_pml=nx+2*n_pml;		 //加上PML后的网格大小
	nz_pml=nz+NN/2+n_pml;
	NX=nx;					 //NX:接受点数
	//isx=nx/2+n_pml; isz=0+NN/2;    //震源位置
	//n_source=par.n_source;
	Source_X=par.Source_X;
	Source_Z=par.Source_Z;
	for (int i = 0; i < Source_X.size(); i++)
	{
		Source_X[i] = Source_X[i]+n_pml;
		Source_Z[i] = Source_Z[i]+NN/2;
	}
	createmodelway=par.createmodelway;
	h=par.h;
	vp=par.Vp;
	vs=par.Vs;
	rou=par.rou;
	FileName_snap=par.FileName_snap;
	FileName_seismicrecord=par.FileName_seismicrecord;
	FilName_model=par.FilName_model;

	CreatResultFolder(FileName_seismicrecord);
	CreatResultFolder(FileName_snap);
	CreatResultFolder(FilName_model);

	Malloc2D(U ,nx_pml,nz_pml);	   Init_0(U,nx_pml,nz_pml);	
	Malloc2D(V ,nx_pml,nz_pml);	    Init_0(V,nx_pml,nz_pml);
	Malloc2D(Ux ,nx_pml,nz_pml);	Init_0(Ux,nx_pml,nz_pml);
	Malloc2D(Uz ,nx_pml,nz_pml);	Init_0(Uz, nx_pml,nz_pml);
	Malloc2D(Vx ,nx_pml,nz_pml);	Init_0(Vx,nx_pml,nz_pml);
	Malloc2D(Vz ,nx_pml,nz_pml);	Init_0(Vz,nx_pml,nz_pml);
	Malloc2D(R ,nx_pml,nz_pml);	    Init_0(R,nx_pml,nz_pml);
	Malloc2D(T ,nx_pml,nz_pml);		Init_0(T,nx_pml,nz_pml);
	Malloc2D(H ,nx_pml,nz_pml);		Init_0(H,nx_pml,nz_pml);
	Malloc2D(Rx ,nx_pml,nz_pml);	Init_0(Rx,nx_pml,nz_pml);
	Malloc2D(Rz ,nx_pml,nz_pml);	Init_0(Rz,nx_pml,nz_pml);
	Malloc2D(Tx ,nx_pml,nz_pml);	Init_0(Tx,nx_pml,nz_pml);
	Malloc2D(Tz ,nx_pml,nz_pml);	Init_0(Tz,nx_pml,nz_pml);
	Malloc2D(Hx ,nx_pml,nz_pml);	Init_0(Hx,nx_pml,nz_pml);
	Malloc2D(Hz ,nx_pml,nz_pml);	Init_0(Hz,nx_pml,nz_pml);
	Malloc2D(Rou ,nx_pml,nz_pml);	Init_0(Rou,nx_pml,nz_pml);
	Malloc2D(Vp ,nx_pml,nz_pml);	Init_0(Vp,nx_pml,nz_pml);
	Malloc2D(Vs ,nx_pml,nz_pml);	Init_0(Vs,nx_pml,nz_pml);
	Malloc2D(mu ,nx_pml,nz_pml);	Init_0(mu,nx_pml,nz_pml);
	Malloc2D(numd ,nx_pml,nz_pml);	Init_0(numd,nx_pml,nz_pml);
	Malloc2D(record ,nx_pml,nt);	Init_0(record,nx_pml,nt);
}
void FD2D::CreateModel()
{
	int i,j,k;
	if(createmodelway==0)
	{
		if(h.size()==1) 
		{
			for(i=0;i<nx_pml;i++)
			{
				for(j=0;j<nz_pml;j++)
				{
					Vp[i][j]=vp.back();
					Vs[i][j]=vs.back();
					Rou[i][j]=rou.back();
				}
			}
		}
		else
		{
			for(i=0;i<h.size();i++)
			{h[i]=h[i]/dz;}
			vector<float> av(h.size()+1, 0.0);					//auxiliary vector
			copy(h.begin(),h.end(),av.begin()+1);
			partial_sum(av.begin(),av.end(), av.begin());	
			
			for(i=0;i<nx_pml;i++)
			{
				for(j=0;j<nz_pml;j++)
				{
					Vp[i][j]=vp.back();
					Vs[i][j]=vs.back();
					Rou[i][j]=rou.back();
					for(k=0;k<h.size()-1;k++)
					{
						if(j>=NN/2+av[k]&&j<NN/2+av[k]+h[k])
						{		
							Vp[i][j]=vp[k];
							Vs[i][j]=vs[k];
							Rou[i][j]=rou[k];
						}
					}
				}
			}
		}
		//for(i=0;i<nx_pml;i++)//空洞模拟
		//{
		//	for(j=0;j<nz_pml;j++)
		//	{
		//		if(i>=139+n_pml&&i<=149+n_pml&&j>=10&&j<20)
		//		{
		//			Rou[i][j]=400.;
		//			Vp[i][j]=340.;
		//			Vs[i][j]=17.;
		//		}
		//	}
		//}
		//Vp[isx][isz]=1.0;		
	}
	if(createmodelway==1)
	{
		float **vp = NULL;
		float **vs = NULL;
		float **den = NULL;
		vp = new float *[nx];
		vs = new float *[nx];
		den = new float *[nx];
		for(i=0;i<nx;i++)
		{
			vp[i] = new float [nz];
			vs[i] = new float [nz];
			den[i]= new float [nz];
		}
		ifstream inputVp("Velocity_Vp.bin",ios::binary);
		ifstream inputVs("Velocity_Vs.bin",ios::binary);
		ifstream inputDen("Den.bin",ios::binary);
		if (!inputVp)
		{
			cout<<"cannot open Vp file data file"<<endl;
			getchar();
		}
		if (!inputVs)
		{
			cout<<"cannot open Vs file data file"<<endl;
			getchar();
		}
		if (!inputVp)
		{
			cout<<"cannot open Density file data file"<<endl;
			getchar();
		}
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				inputVp.read((char*)&vp[i][j],sizeof(vp[i][j]));
				inputVs.read((char*)&vs[i][j],sizeof(vs[i][j]));
				inputDen.read((char*)&den[i][j],sizeof(den[i][j]));
			}
		for(i=0;i<nx;i++)
		{
			for(j=0;j<nz;j++)
			{	Vp[i+n_pml][j+NN/2] = vp[i][j];
				Vs[i+n_pml][j+NN/2] = vs[i][j];
				Rou[i+n_pml][j+NN/2]= den[i][j];
			}
		}
		for(i=0;i<nx_pml;i++)//上
		{
			for(j=0;j<NN/2;j++)
			{	Vp[i][j] = 0.0;
				Vs[i][j] = 0.0;
				Rou[i][j]= 0.0;
			}
		}
		for(i=0;i<n_pml;i++)//左
		{
			for(j=NN/2;j<nz+NN/2;j++)
			{	
				Vp[i][j] = Vp[n_pml][j];
				Vs[i][j] = Vs[n_pml][j];
				Rou[i][j]= Rou[n_pml][j];
			}
		}
		for(i=nx+n_pml;i<nx_pml;i++)//右
		{
			for(j=NN/2;j<nz+NN/2;j++)
			{	
				Vp[i][j] = Vp[nx+n_pml-1][j];
				Vs[i][j] = Vs[nx+n_pml-1][j];
				Rou[i][j]= Rou[nx+n_pml-1][j];
			}
		}
		for(i=0;i<nx_pml;i++)//下
		{
			for(j=nz+NN/2;j<nz_pml;j++)
			{	Vp[i][j] = Vp[i][nz+NN/2-1];
				Vs[i][j] = Vs[i][nz+NN/2-1];
				Rou[i][j]= Rou[i][nz+NN/2-1];
			}
		}
	}
	for(i=0;i<nx_pml;i++)
	{
		for(j=0;j<nz_pml;j++)
		{
			mu[i][j]=Rou[i][j]*pow(Vs[i][j],2);
			numd[i][j]=Rou[i][j]*(pow(Vp[i][j],2)-2.0*pow(Vs[i][j],2));
		}
	}

	string FilePath = "./"+string(FilName_model)+"/";
	char vpF[100];strcpy_s(vpF,(FilePath+"Model_Vp").c_str());
	char vsF[100];strcpy_s(vsF,(FilePath+"Model_Vs").c_str());
	char rouF[100];strcpy_s(rouF,(FilePath+"Model_rou").c_str());
	outputfile(Vp,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,vpF);
	outputfile(Vs,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,vsF);
	outputfile(Rou,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,rouF);
}
void FD2D::pre_pmlabc(int i,int j,float dUx,float dUz,float dVz,float dVx)
{
	float ad_Tx,ad_Rx,ad_Hx,ad_Tz,ad_Rz,ad_Hz;
	if(i<n_pml)					//左吸收边界
	{
		ad_Tx=d_pml(dx*((float)(n_pml-i)),Vp[i][j],dx,n_pml);
		ad_Rx=ad_Tx;
		ad_Hx=d_pml(dx*((float)(n_pml-i)-dx/2.0),Vp[i][j],dx,n_pml);
	}
	else if(i>=nx+n_pml)		//右吸收边界
	{
		ad_Tx=d_pml(dx*((float)(i-nx-n_pml)),Vp[i][j],dx,n_pml);
		ad_Rx=ad_Tx;
		ad_Hx=d_pml(dx*((float)(i-nx-n_pml)+dx/2.0),Vp[i][j],dx,n_pml);
	}
	else
	{
		ad_Tx=0.0;
		ad_Rx=0.0;
		ad_Hx=0.0;
	}
	if(j>=nz+NN/2)  //下吸收边界
	{
		ad_Tz=d_pml(dz*(j-nz-NN/2),Vp[i][j],dz,n_pml);
		ad_Rz=ad_Tz;
		ad_Hz=d_pml(dz*((float)(j-nz-NN/2)+dz/2.0),Vp[i][j],dz,n_pml);
	}
	else
	{
		ad_Tz=0.0;
		ad_Rz=0.0;
		ad_Hz=0.0;
	}
	Rx[i][j]=Rx[i][j]*(2.0-dt*ad_Rx)/(2.0+dt*ad_Rx)+2.0*dt*\
			   (numd[i][j]+2*mu[i][j])/((2.0+dt*ad_Rx)*dx)*dUx;

	Tx[i][j]=Tx[i][j]*(2.0-dt*ad_Tx)/(2.0+dt*ad_Tx)+2.0*dt*\
			   numd[i][j]/((2.0+dt*ad_Tx)*dx)*dUx;	

	Rz[i][j]=Rz[i][j]*(2.0-dt*ad_Rz)/(2.0+dt*ad_Rz)+2.0*dt*\
		       numd[i][j]/((2.0+dt*ad_Rz)*dz)*dVz;

	Tz[i][j]=Tz[i][j]*(2.0-dt*ad_Tz)/(2.0+dt*ad_Tz)+2.0*dt*\
			   (numd[i][j]+2*mu[i][j])/((2.0+dt*ad_Tz)*dz)*dVz;

	Hx[i][j]=Hx[i][j]*(2.0-dt*ad_Hx)/(2.0+dt*ad_Hx)+2.0*dt*\
			   mu[i][j]/((2.0+dt*ad_Hx)*dx)*dVx;	

	Hz[i][j]=Hz[i][j]*(2.0-dt*ad_Hz)/(2.0+dt*ad_Hz)+2.0*dt*\
			   mu[i][j]/((2.0+dt*ad_Hz)*dz)*dUz;

	R[i][j]=Rx[i][j]+Rz[i][j];
	T[i][j]=Tx[i][j]+Tz[i][j];
	H[i][j]=Hx[i][j]+Hz[i][j];
}
void FD2D::vel_pmlabc(int i,int j,float dRx,float dHx_U,float dHx_V,float dTz)
{
	float ad_Ux,ad_Vx,ad_Uz,ad_Vz;
	if(i<n_pml)				//左吸收边界
	{
		ad_Ux=d_pml(dx*((float)(n_pml-i)-dx/2),Vp[i][j],dx,n_pml);
		ad_Vx=d_pml(dx*((float)(n_pml-i)),Vp[i][j],dx,n_pml);
	}
	else if(i>=nx+n_pml)    //右吸收边界
	{
		ad_Ux=d_pml(dx*((float)(i-nx-n_pml)+dx/2),Vp[i][j],dx,n_pml);
		ad_Vx=d_pml(dx*((float)(i-nx-n_pml)),Vp[i][j],dx,n_pml);
	}
	else
	{
		ad_Ux=0.0;
		ad_Vx=0.0;
	}
	if(j>=nz+NN/2)  //下吸收边界
	{
		ad_Uz=d_pml(dz*(j-nz-NN/2),Vp[i][j],dz,n_pml);
		ad_Vz=d_pml(dz*((float)(j-nz-NN/2)+0.5),Vp[i][j],dz,n_pml);
	}
	else
	{
		ad_Uz=0.0;
		ad_Vz=0.0;
	}
	Ux[i][j]=Ux[i][j]*(2.0-dt*ad_Ux)/(2.0+dt*ad_Ux)\
				+2.0*dt/(Rou[i][j]*dx*(2.0+dt*ad_Ux))*dRx;
	Vx[i][j]=Vx[i][j]*(2.0-dt*ad_Vx)/(2.0+dt*ad_Vx)\
				+2.0*dt/(Rou[i][j]*dx*(2.0+dt*ad_Vx))*dHx_V;
	Uz[i][j]=Uz[i][j]*(2.0-dt*ad_Uz)/(2.0+dt*ad_Uz)\
				+2.0*dt/(Rou[i][j]*dz*(2.0+dt*ad_Uz))*dHx_U;
	Vz[i][j]=Vz[i][j]*(2.0-dt*ad_Vz)/(2.0+dt*ad_Vz)\
				+2.0*dt/(Rou[i][j]*dz*(2.0+dt*ad_Vz))*dTz;
	U[i][j]=Ux[i][j]+Uz[i][j];
	V[i][j]=Vx[i][j]+Vz[i][j];
}
void FD2D::fd2d_hoeizontal_surface()
{
	float C1=1.196289;
	float C2= -0.0797526;
	float C3= 0.009570313;
	float C4= -0.0006975447;
	float bb= 0.1;	        //震源衰减系数
	int i,j,k;
	float fz,dRx,dHx_U,dHx_V,dTz,dUx,dUz,dVz,dVx;
	char Snap_Vx_Record[256],Snap_Vz_Record[256];
	/*---------------自由边界条件------------------------*/
	for(j=0;j<=NN/2;j++)//for(j=0;j<NN/2;j++)
	{
		for(i=0;i<nx_pml;i++)
		{
			numd[i][j]=0.0;
			mu[i][j]=mu[i][j]*0.5;
			//mu[i][j]=mu[i][j];
			Rou[i][j]=Rou[i][j]*0.5;
		}
	}
/*-------------------------------------------------------------------------------------------------------------------------------*/
/*                  R:Txx    T:Tzz   H:Txz     U:Vx    V:Vz                                                                      */
/*--------------------------------------------------根据差分格式由U/V计算R/H/T-----------------------------------------------------*/
	for(k=0;k<nt;k++)
	{
		//if (k%100==0)
			cout<<k<<endl;
/*-------------------------------------------自由边界条件--------------------------------------------------------------------------------------*/
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			T[i][NN/2]=0.0;
		}
/*---------------------------------------------根据差分格式由R/H/T计算U/V--------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for \
	schedule(dynamic) \
	private(i,j,dRx,dHx_U,dHx_V,dTz,fz) /*\
	shared(Rou, Vp, numd, mu, U, V, \
	R, T, H, Ux, Uz, Vx, Vz)*/
#endif
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			for(j=NN/2;j<nz_pml-NN/2;j++)
			{
				dRx=C1*(R[i+1][j]-R[i][j])+C2*(R[i+2][j]-R[i-1][j])+\
					C3*(R[i+3][j]-R[i-2][j])+C4*(R[i+4][j]-R[i-3][j]);
				dHx_U=C1*(H[i][j]-H[i][j-1])+C2*(H[i][j+1]-H[i][j-2])+\
					C3*(H[i][j+2]-H[i][j-3])+C4*(H[i][j+3]-H[i][j-4]);
				dHx_V=C1*(H[i][j]-H[i-1][j])+C2*(H[i+1][j]-H[i-2][j])+\
					C3*(H[i+2][j]-H[i-3][j])+C4*(H[i+3][j]-H[i-4][j]);
				dTz=C1*(T[i][j+1]-T[i][j])+C2*(T[i][j+2]-T[i][j-1])+\
					C3*(T[i][j+3]-T[i][j-2])+C4*(T[i][j+4]-T[i][j-3]);
				//if(i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)           // 模型区域
				//{
				//	fz=exp(-1.0*bb*(pow(float((i-isx)),2)+pow(float((j-isz)),2)))*Gausswavelet(dt*(k+0.5));
				//	U[i][j]=U[i][j]+dt/Rou[i][j]*(dRx/dx+dHx_U/dz);
				//	V[i][j]=V[i][j]+dt/Rou[i][j]*(dHx_V/dx+dTz/dz)+fz;
				//}
				
				if(i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)           // 模型区域
				{
					fz=0.0;
					for(int ii=0; ii<Source_X.size(); ii++)
					{
						fz=fz+exp(-1.0*bb*(pow(float((i-Source_X[ii])),2)+pow(float((j-Source_Z[ii])),2)))\
							*Gausswavelet(dt*(k+0.5));
					}
					U[i][j]=U[i][j]+dt/Rou[i][j]*(dRx/dx+dHx_U/dz);
					V[i][j]=V[i][j]+dt/Rou[i][j]*(dHx_V/dx+dTz/dz)+fz;
				}
				else                                          //pml吸收边界条件     
				{
					vel_pmlabc(i,j, dRx, dHx_U, dHx_V, dTz);
				}
			}
		}
#ifdef _OPENMP
#pragma omp parallel for \
	schedule(dynamic) \
	private(i,j,dUx,dUz,dVz,dVx) /*\
	shared(Rou, Vp, mu, numd, U, V,\
	R, T, H, Rx, Rz, Tx, Tz, Hx, Hz)*/
#endif
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			for(j=NN/2;j<nz_pml-NN/2;j++)
			{
				/*差分格式*/
				dUx=C1*(U[i][j]-U[i-1][j])+C2*(U[i+1][j]-U[i-2][j])+\
					C3*(U[i+2][j]-U[i-3][j])+C4*(U[i+3][j]-U[i-4][j]);
				dUz=C1*(U[i][j+1]-U[i][j])+C2*(U[i][j+2]-U[i][j-1])+\
					C3*(U[i][j+3]-U[i][j-2])+C4*(U[i][j+4]-U[i][j-3]);
				dVz=C1*(V[i][j]-V[i][j-1])+C2*(V[i][j+1]-V[i][j-2])+\
					C3*(V[i][j+2]-V[i][j-3])+C4*(V[i][j+3]-V[i][j-4]);
				dVx=C1*(V[i+1][j]-V[i][j])+C2*(V[i+2][j]-V[i-1][j])+\
					C3*(V[i+3][j]-V[i-2][j])+C4*(V[i+4][j]-V[i-3][j]);
				if (i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)                  //模型区域
				{
					R[i][j]=R[i][j]+dt*((numd[i][j]+2*mu[i][j])*\
						dUx/dx+numd[i][j]*dVz/dz);
				    T[i][j]=T[i][j]+dt*(numd[i][j]*dUx/dx+\
						(numd[i][j]+2*mu[i][j])*dVz/dz);
				    H[i][j]=H[i][j]+dt*(mu[i][j]*dVx/dx+mu[i][j]*dUz/dz);
				}
				else                                                 //pml吸收边界条件
				{
					pre_pmlabc(i,j,dUx,dUz,dVz,dVx);
				}
			}
		}
/*------------------------------------------------------------------获取波场快照---------------------------------------------------------*/				
       if(k%50==0)
		{
			string FilePath = "./"+string(FileName_snap)+"/";
			char Snap_Vx_Record[256],Snap_Vz_Record[256];
			sprintf_s(Snap_Vx_Record,"Snap_Vx_%d",k);
			sprintf_s(Snap_Vz_Record,"Snap_Vz_%d",k);
			strcpy_s(Snap_Vx_Record,(FilePath+string(Snap_Vx_Record)).c_str());
			strcpy_s(Snap_Vz_Record,(FilePath+string(Snap_Vz_Record)).c_str());
			outputfile(U,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,Snap_Vx_Record);
			outputfile(V,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,Snap_Vz_Record);
		}
/*-----------------获取地震记录------------------------*/
	   for(i=0;i<NX;i++)
	   {
		   record[i][k]=V[i+n_pml][NN/2];
	   }
/*---------------------------------------------------------------------------------------------------------------------------------------*/
	}
	string RecordPath = "./"+string(FileName_seismicrecord)+"/";
	char seis[100];
	strcpy_s(seis,(RecordPath+"Seismic Record").c_str());
	outputfile(record,0,NX,0,nt,seis);
}
void FD2D::fd2d_fluctuant_surface()
{
	float C1=1.196289;
	float C2= -0.0797526;
	float C3= 0.009570313;
	float C4= -0.0006975447;
	float bb= 0.1;	                           //震源衰减系数
	int i,j,k;
	float fz,dRx,dHx_U,dHx_V,dTz,dUx,dUz,dVz,dVx;
	char Snap_Vx_Record[256],Snap_Vz_Record[256];
	/*---------------自由边界条件AEA------------------------*/
	for(j=0;j<=NN/2;j++)
	{
		for(i=0;i<nx_pml;i++)
		{
			numd[i][j]=0.0;
			mu[i][j]=mu[i][j]*0.5;
			//mu[i][j]=mu[i][j];
			Rou[i][j]=Rou[i][j]*0.5;
		}
	}
/*-------------------------------------------------------------------------------------------------------------------------------*/
/*                  R:Txx    T:Tzz   H:Txz     U:Vx    V:Vz                                                                      */
/*-------------------------------------------------------------------------------------------------------------------------------*/
	for(k=0;k<nt;k++)
	{
		//if (k%100==0)
			cout<<k<<endl;
		/*---------------自由边界条件AEA------------------------*/
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			T[i][NN/2]=0.0;
		}
/*---------------------------------------------根据差分格式由R/H/T计算U/V--------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for \
	schedule(dynamic) \
	private(i,j,dRx,dHx_U,dHx_V,dTz,fz) /*\
	shared(Rou, Vp, numd, mu, U, V, \
	R, T, H, Ux, Uz, Vx, Vz)*/
#endif
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			for(j=NN/2;j<nz_pml-NN/2;j++)
			{
				dRx=C1*(R[i+1][j]-R[i][j])+C2*(R[i+2][j]-R[i-1][j])+\
					C3*(R[i+3][j]-R[i-2][j])+C4*(R[i+4][j]-R[i-3][j]);
				dHx_U=C1*(H[i][j]-H[i][j-1])+C2*(H[i][j+1]-H[i][j-2])+\
					C3*(H[i][j+2]-H[i][j-3])+C4*(H[i][j+3]-H[i][j-4]);
				dHx_V=C1*(H[i][j]-H[i-1][j])+C2*(H[i+1][j]-H[i-2][j])+\
					C3*(H[i+2][j]-H[i-3][j])+C4*(H[i+3][j]-H[i-4][j]);
				dTz=C1*(T[i][j+1]-T[i][j])+C2*(T[i][j+2]-T[i][j-1])+\
					C3*(T[i][j+3]-T[i][j-2])+C4*(T[i][j+4]-T[i][j-3]);
				//if(i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)           // 模型区域
				//{
				//	fz=exp(-1.0*bb*(pow(float((i-isx)),2)+pow(float((j-isz)),2)))*Gausswavelet(dt*(k+0.5));
				//	U[i][j]=U[i][j]+dt/Rou[i][j]*(dRx/dx+dHx_U/dz);
				//	V[i][j]=V[i][j]+dt/Rou[i][j]*(dHx_V/dx+dTz/dz)+fz;
				//}
				
				if(i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)           // 模型区域
				{
					fz=0.0;
					for(int ii=0; ii<Source_X.size(); ii++)
					{
						fz=fz+exp(-1.0*bb*(pow(float((i-Source_X[ii])),2)+pow(float((j-Source_Z[ii])),2)))\
							*Gausswavelet(dt*(k+0.5));
					}
					U[i][j]=U[i][j]+dt/Rou[i][j]*(dRx/dx+dHx_U/dz);
					V[i][j]=V[i][j]+dt/Rou[i][j]*(dHx_V/dx+dTz/dz)+fz;
				}
				else                                          //pml吸收边界条件     
				{
					vel_pmlabc(i,j, dRx, dHx_U, dHx_V, dTz);
				}
			}
		}
/*--------------------------------------------------根据差分格式由U/V计算R/H/T-----------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for \
	schedule(dynamic) \
	private(i,j,dUx,dUz,dVz,dVx) /*\
	shared(Rou, Vp, mu, numd, U, V,\
	R, T, H, Rx, Rz, Tx, Tz, Hx, Hz)*/
#endif
		for(i=NN/2;i<nx_pml-NN/2;i++)
		{
			for(j=NN/2;j<nz_pml-NN/2;j++)
			{
				/*差分格式*/
				dUx=C1*(U[i][j]-U[i-1][j])+C2*(U[i+1][j]-U[i-2][j])+\
					C3*(U[i+2][j]-U[i-3][j])+C4*(U[i+3][j]-U[i-4][j]);
				dUz=C1*(U[i][j+1]-U[i][j])+C2*(U[i][j+2]-U[i][j-1])+\
					C3*(U[i][j+3]-U[i][j-2])+C4*(U[i][j+4]-U[i][j-3]);
				dVz=C1*(V[i][j]-V[i][j-1])+C2*(V[i][j+1]-V[i][j-2])+\
					C3*(V[i][j+2]-V[i][j-3])+C4*(V[i][j+3]-V[i][j-4]);
				dVx=C1*(V[i+1][j]-V[i][j])+C2*(V[i+2][j]-V[i-1][j])+\
					C3*(V[i+3][j]-V[i-2][j])+C4*(V[i+4][j]-V[i-3][j]);
				if (i>=n_pml&&i<nx+n_pml&&j<nz+NN/2)                  //模型区域
				{
					R[i][j]=R[i][j]+dt*((numd[i][j]+2*mu[i][j])*dUx/dx+numd[i][j]*dVz/dz);
				    T[i][j]=T[i][j]+dt*(numd[i][j]*dUx/dx+(numd[i][j]+2*mu[i][j])*dVz/dz);
				    H[i][j]=H[i][j]+dt*(mu[i][j]*dVx/dx+mu[i][j]*dUz/dz);
				}
				else                                                 //pml吸收边界条件
				{
					pre_pmlabc(i,j,dUx,dUz,dVz,dVx);
				}
			}
		}
/*------------------------------------------------------------------获取波场快照---------------------------------------------------------*/				
       if(k%50==0)
		{
			string FilePath = "./"+string(FileName_snap)+"/";
			char Snap_Vx_Record[256],Snap_Vz_Record[256];
			sprintf_s(Snap_Vx_Record,"Snap_Vx_%d",k);
			sprintf_s(Snap_Vz_Record,"Snap_Vz_%d",k);
			strcpy_s(Snap_Vx_Record,(FilePath+string(Snap_Vx_Record)).c_str());
			strcpy_s(Snap_Vz_Record,(FilePath+string(Snap_Vz_Record)).c_str());
			outputfile(U,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,Snap_Vx_Record);
			outputfile(V,n_pml,nx_pml-n_pml,NN/2,nz_pml-n_pml,Snap_Vz_Record);
		}
/*-----------------获取地震记录------------------------*/
	   for(i=0;i<NX;i++)
	   {
		   record[i][k]=V[i+n_pml][NN/2];
	   }
/*---------------------------------------------------------------------------------------------------------------------------------------*/
	}
	string RecordPath = "./"+string(FileName_seismicrecord)+"/";
	char seis[100];
	strcpy_s(seis,(RecordPath+"Seismic Record").c_str());
	outputfile(record,0,NX,0,nt,seis);
}
