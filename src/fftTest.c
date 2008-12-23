#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include"fftw3.h"

#define N 100
//#define T 300
//#define DF -0.000102
int main()
{
	int i;
	fftw_plan     fftwPlan;
	fftw_complex *fftTableInB;
	fftw_complex *fftTableOut;
		
	fftTableInB = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftTableOut = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	
	for(i=0;i<N;i++)
	{
		fftTableInB[i][0] = cos(0.2*i);
		fftTableInB[i][1] = 0;//sin(0.2*i);		
	}
	
	fftwPlan = fftw_plan_dft_1d(N,fftTableInB,fftTableOut,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(fftwPlan);
		
	fftw_free(fftTableInB);
	fftw_free(fftTableOut);
	

	return 0;
}

