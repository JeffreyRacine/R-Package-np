/* A sleep utility so that R can be started at different time 
Copyright 2005 Hao Yu */

#include <stdlib.h>
// #include <stdio.h>
#include <math.h>
#include <WinSock.h>


int main(int argc, char *argv[])
{
	//int ms=(int)(500atof(argv[1])/32767); 
	int ms=atoi(argv[1]);
	
	// ms = (int)floor(ran/step);
	// printf("ms = %d\n",ms);
	
	Sleep(ms);

	return 0;
}
