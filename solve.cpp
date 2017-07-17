#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cctype>
#include<algorithm>
#include<cmath>
#include<string>
#define ll long long
#define ull unsigned long long
#define inf 0x3f3f3f3f
using namespace std;
double mat[6][6]={
    4,-1,0,-1,0,0,
    -1,4,-1,0,-1,0,
    0,-1,4,-1,0,-1,
    -1,0,-1,4,-1,0,
    0,-1,0,-1,4,-1,
    0,0,-1,0,-1,4
};
double bXJ[6]={0,0,0,0,0,0};
double bGJ[6]={0,0,0,0,0,0};
double bSJ[6]={0,0,0,0,0,0};
double XJ[6]={0,0,0,0,0,0};
double GJ[6]={0,0,0,0,0,0};
double SJ[6]={0,0,0,0,0,0};
double B[6]={0,5,-2,5,-2,6};
bool judge(int type)
{
    if(type==1)
    {
        double res=0;
        for(int i=1;i<6;++i)
        {
            res+=(bXJ[i]-XJ[i])*(bXJ[i]-XJ[i]);
        }
        res=sqrt(res);
        if(res<=0.0001) return true;
        else return false;
    }
    else if(type==2)
    {
        double res=0;
        for(int i=1;i<6;++i)
        {
            res+=(bGJ[i]-GJ[i])*(bGJ[i]-GJ[i]);
        }
        res=sqrt(res);
        if(res<=0.0001) return true;
        else return false;
    }
    else
    {
        double res=0;
        for(int i=1;i<6;++i)
        {
            res+=(bSJ[i]-SJ[i])*(bSJ[i]-SJ[i]);
        }
        res=sqrt(res);
        if(res<=0.0001) return true;
        else return false;
    }
}
void Jacobi()
{
    int k=1;
    int type=1;
    while(!judge(type)||k==1)
    {
        for(int i=0;i<6;++i)
        {
            bXJ[i]=XJ[i];
        }
        printf("No.%d: ",k);
        for(int i=0;i<6;++i)
        {
            double sum=0;
            for(int j=0;j<6;++j)
            {
                if(i!=j) sum+=mat[i][j]*bXJ[j];
            }
            XJ[i]=(B[i]-sum)/mat[i][i];
        }
        for(int i=0;i<6;++i) printf("%f ",XJ[i]);
        cout<<endl;
        k++;
    }
}
void GS()
{
    int k=1;
    int type=2;
    while(!judge(type)||k==1)
    {
        for(int i=0;i<6;++i)
        {
            bGJ[i]=GJ[i];
        }
        printf("No.%d: ",k);
        for(int i=0;i<6;++i)
        {
            double sum=0;
            for(int j=0;j<i;++j)
            {
                sum+=mat[i][j]*GJ[j];
            }
            for(int j=i+1;j<6;++j)
            {
                sum+=mat[i][j]*bGJ[j];
            }
            GJ[i]=(B[i]-sum)/mat[i][i];
        }
        for(int i=0;i<6;++i) printf("%f ",GJ[i]);
        cout<<endl;
        k++;
    }
}
void SOR(double val)
{
    memset(SJ,0,sizeof(SJ));
    memset(bSJ,0,sizeof(bSJ));
    int k=1;
    int type=3;
    while(!judge(type)||k==1)
    {
        for(int i=0;i<6;++i)
        {
            bSJ[i]=SJ[i];
        }
        printf("No.%d: ",k);
        for(int i=0;i<6;++i)
        {
            double sum=0;
            for(int j=0;j<i;++j)
            {
                sum+=mat[i][j]*SJ[j];
            }
            for(int j=i;j<6;++j)
            {
                sum+=mat[i][j]*bSJ[j];
            }
            SJ[i]=bSJ[i]+val/mat[i][i]*(B[i]-sum);
        }
        for(int i=0;i<6;++i) printf("%f ",SJ[i]);
        cout<<endl;
        k++;
    }
}
int main() 
{ 
    printf("Jacobi\n");
    Jacobi();
    printf("\nGauss-Seidel\n");
    GS();
    printf("\nSOR w=1.334\n");
    SOR(1.334);
    printf("\nSOR w=1.95\n");
    SOR(1.95);
    printf("\nSOR w=0.95\n");
    SOR(0.95);
    getchar();
    return 0;  
}