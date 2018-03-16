#include"swfd.h"

/*------------------------速度等模型参数在swfd.cpp中设置--------------------------------------------------*/
int main(int argc, char*argv[])
{
	clock_t start,finish;
	double totaltime;
	start=clock();
/*--------------------网格参数设置----------------------*/
	float dt=0.0005, dx=1.0, dz=1.0; //网格步长和时间步长
	int nx=500, nz=500, nt=1500;      //模型区域大小和时间记录长度	
	int pml=50;                      //吸收边界网格点数
/*--------------------模型参数设置----------------------*/
	int createmodelway=0;					 //0或1,0:输入下面参数创建，1：程序读入文件，则下面参数可以不更改
	const int n_layer=1;					 //地层层数,1:均匀介质
	float h[n_layer]={100};	//层厚，单位米，最好为网格步长的整数倍,最后一层不管，除开前面层厚，剩下的就是最后一层
	float vp[n_layer]={1000.0};
	float vs[n_layer]={vp[0]/1.73};
	float rou[n_layer]={1800.0};
/*--------------------震源参数设置----------------------*/
	const int n_source=1;			     		 //震源个数
	int source_X[n_source]={nx/2};				 //震源位置横坐标
	int source_Z[n_source]={0};				     //震源位置纵坐标
/*--------------------输出文件----------------------*/
	char *FileName_snap="snap";
	char *FileName_seismicrecord="seismic record";
	char *FilName_model="model";
/*----------------------------------------------------*/
/*----------------------------------------------------*/
	vector<float> Source_X(source_X,source_X+n_source);
	vector<float> Source_Z(source_Z,source_Z+n_source);
	vector<float> H(h,h+n_layer);
	vector<float> Vp(vp,vp+n_layer);
	vector<float> Vs(vs,vs+n_layer);
	vector<float> Rou(rou,rou+n_layer);
	FD2DP par; 
	SetGridParmeters(par, dt, dx, dz, nx, nz, nt, pml, \
					FileName_snap,FileName_seismicrecord, FilName_model);
	SetModelParmeters(par,createmodelway,Source_X,Source_Z,H,Vp,Vs,Rou);
	FD2D fd2d(par);
	fd2d.CreateModel();
	fd2d.fd2d_fluctuant_surface();

	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"\n     此程序的运行时间为:  "<<totaltime<<"秒"<<endl;
	system("pause");
	return(0);
}