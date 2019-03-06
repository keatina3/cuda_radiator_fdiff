#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/time.h>

int main(int argc, char **argv){
	int m = 32, n = 32, p = 10;
	int t = 0, w = 0, option = 0;
	struct timeval start, end;

	while((option=getopt(argc,argv,"n:m:b:rtwp"))!=-1){
		switch(option){
			case 'n': n = atoi(optarg);
				break;
			case 'm': m = atoi(optarg);
				break;
			case 'p': p = atoi(optarg);
				break;
			case 't': t = 1;
				break;
			case 'w': w = 1;	// to write results to file //
				break;
			default:
				printf("Incorrect options entered!\n");
				return 1;
		}
	}	
	if(argc != optind){
		printf("Too many arguments provided, exiting!\n");
		return 1;
	}

	return 0;
}
