#include"swfd.h"

/*------------------------�ٶȵ�ģ�Ͳ�����swfd.cpp������--------------------------------------------------*/
int main(int argc, char*argv[])
{
	clock_t start,finish;
	double totaltime;
	start=clock();
/*--------------------�����������----------------------*/
	float dt=0.0005, dx=1.0, dz=1.0; //���񲽳���ʱ�䲽��
	int nx=500, nz=500, nt=1500;      //ģ�������С��ʱ���¼����	
	int pml=50;                      //���ձ߽��������
/*--------------------ģ�Ͳ�������----------------------*/
	int createmodelway=0;					 //0��1,0:�����������������1����������ļ���������������Բ�����
	const int n_layer=1;					 //�ز����,1:���Ƚ���
	float h[n_layer]={100};	//��񣬵�λ�ף����Ϊ���񲽳���������,���һ�㲻�ܣ�����ǰ����ʣ�µľ������һ��
	float vp[n_layer]={1000.0};
	float vs[n_layer]={vp[0]/1.73};
	float rou[n_layer]={1800.0};
/*--------------------��Դ��������----------------------*/
	const int n_source=1;			     		 //��Դ����
	int source_X[n_source]={nx/2};				 //��Դλ�ú�����
	int source_Z[n_source]={0};				     //��Դλ��������
/*--------------------����ļ�----------------------*/
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
	cout<<"\n     �˳��������ʱ��Ϊ:  "<<totaltime<<"��"<<endl;
	system("pause");
	return(0);
}