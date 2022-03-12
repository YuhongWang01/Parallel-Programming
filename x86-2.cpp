#include <iostream>
#include<windows.h>
#include<stdlib.h>
#include<ctime>
#include<cstdlib>
using namespace std;






int main()
{

    int n=4096*4;
    cout<<"��������ģ��"<<n<<endl;
    int*a=new int[n];
    int sum=0;

    long long head,tail,freq,times=0;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );


    for(int i=0;i<n;i++)
    {
        srand(time(NULL));

        a[i]=rand()%100+1;
    }
    int c=100;//ѭ������


    for(int k=c;k>0;k--)
    {
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for (int i = 0; i < n; i++)
            sum += a[i];
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        times+=tail-head;
    }


    cout<<"sum:"<<sum/c<<endl;
    cout<<"ƽ���㷨ʱ�䣺"<<times*1000.0/freq/c<<"ms"<<endl;

    int sum1 = 0, sum2 = 0;
    times=0;
    for(int k=c;k>0;k--)
    {
        QueryPerformanceCounter((LARGE_INTEGER *)&head);

        for (int i = 0;i < n; i += 2) {
            sum1 += a[i];
            sum2 += a[i + 1];
        }
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        times+=tail-head;

    }

    sum = sum1 + sum2;
    cout<<"sum:"<<sum/c<<endl;
    cout<<"����·ʽ�㷨ʱ�䣺"<<times*1000.0/freq/c<<"ms"<<endl;






    times=0;

    for(int k=c;k>0;k--)
    {
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for (int m = n; m > 1; m /= 2) // log(n)������
            for (int i = 0; i < m / 2; i++)
                a[i] = a[i * 2] + a[i * 2 + 1];
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        times+=tail-head;
    }

    cout<<"sum:"<<a[0]<<endl;
    cout<<"����ѭ���㷨ʱ�䣺"<<times*1000.0/freq/c<<"ms"<<endl;
    return 0;
}
