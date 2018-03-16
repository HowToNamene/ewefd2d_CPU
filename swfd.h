#ifndef  _Vel_Pre_H_
#define _Vel_Pre_H_
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<ctime>
#include <algorithm>
#include<numeric>
#include <io.h> 
#include <direct.h> 
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef struct FD2DParemter
{
	int nx_pml,nz_pml;//,isx,isz;
	int pml;
	int nx,nz,nt;
	float dt,dx,dz;
	int createmodelway;
	//int n_source;
	vector<float> Source_X;
	vector<float> Source_Z;
	int n_layer;
	vector<float> h;
	vector<float> Vp;
	vector<float> Vs;
	vector<float> rou;
	char *FileName_snap;
	char *FileName_seismicrecord;
	char *FilName_model;
}FD2DP;
void SetGridParmeters(FD2DP &,float dt,float dx,float dz,int nx,int nz,int nt,int pml,\
				   char *FileName_snap,\
				  char *FileName_seismicrecord,char *FilName_model);
void SetModelParmeters(FD2DP &,int createmodelway,vector<float>&Source_X,vector<float>&Source_Z,\
					   vector<float>&h,vector<float>&Vp,vector<float>&Vs,vector<float>&rou);
typedef class FD2D_openmp
{
public:
	FD2D_openmp(const FD2DP &par);
	void CreateModel();
	void pre_pmlabc( const int i,const int j,float dUx,float dUz,float dVz,float dVx);
	void vel_pmlabc(const int i,const int j,float dRx,float dHx_U,float dHx_V,float dTz);
	void fd2d_hoeizontal_surface();
	void fd2d_fluctuant_surface();
private:
	enum FDorder{NN=8};
	int nx_pml,nz_pml;//,isx,isz;
	int pml,n_pml;
	int nx,nz,nt,NX;
	float dt,dx,dz;
	//int n_source;
	vector<float> Source_X;
	vector<float> Source_Z;
	int createmodelway;
	vector<float> h;
	vector<float> vp;
	vector<float> vs;
	vector<float> rou;
	char *FileName_snap;
	char *FileName_seismicrecord;
	char *FilName_model;
private:
	float **U,**V,**Ux,**Uz,**Vx,**Vz;
	float **R,**T,**H,**Rx,**Rz,**Tx,**Tz,**Hx,**Hz;
	float **Rou,**Vp,**Vs,**mu,**numd;
	float **record;

}FD2D;
#endif