#include "time.h"
#include "nmblas.h"


#define SIZE 1024

#ifdef __GNUC__ //  NMC-GCC C++ compilier 
double buffer_a[SIZE] __attribute__ ((section (".data_imu1")));
double buffer_b[SIZE] __attribute__ ((section (".data_imu2")));
double buffer_c[SIZE] __attribute__ ((section (".data_imu3")));
double buffer_d[SIZE] __attribute__ ((section (".data_imu4")));
#else 		// NMSDK C++ compiler
#pragma data_section ".data_imu1"
double buffer_a[SIZE];
#pragma data_section ".data_imu2"
double buffer_b[SIZE];
#pragma data_section ".data_imu3"
double buffer_c[SIZE];
#pragma data_section ".data_imu4"
double buffer_d[SIZE];

#endif 

int main(){

	clock_t t0,t1;
	int i;
	for (i=0; i<SIZE; i++){
		buffer_a[i]=i+1;
		buffer_b[i]=2;
		buffer_c[i]=1000;
	}
	//d =accmul_cpp(a, b, c);
	t0=clock();
	accmul32_64f(buffer_a, buffer_b, buffer_c, buffer_d, SIZE);
	t1=clock();
	return t1-t0;
}

// SIZE CLOCKS
// 1024 4343  4343/1024 = 4,2412109375 (a,b,c,d,code in the same bank)
// 1024 4179  4179/1024 = 4,0810546875 (a,d,b,d - in one bank & code in other bnak)
// 1024 1639 mc12101run: 1639/1024 = 1,6005859375 (a,b,c,d,code - in different banks)
// 1024 1639 verilog   : 1674/1024 = 1,6305859375 (a,b,c,d,code - in different banks)
 
