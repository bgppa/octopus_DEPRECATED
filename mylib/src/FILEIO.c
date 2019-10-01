/* Library for managing file IO interaction.
 * TO CHECK AND IMPROVE */
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

/* Read a matrix of numbers from a file.
 * Arguments are:
 * file_name, then a target array that will contain the data
 * (so list_points has to be at leat of dimension dime_each times lines)
 * dim_each : how many numbers are in each line
 * lines : how many lines */
int readPoints(char* file_name, double* list_points, int dim_each, int lines){
	assert(file_name != NULL);
	assert(list_points != NULL);

	if(dim_each < 1 || lines < 1){
		printf("*error*: readPoints: negative values (%d, %d)\n", dim_each, lines);
		return 0;
	}

        FILE* file = fopen(file_name, "r");
        if(file == NULL){
                printf("*error* unable to open %s\n", file_name);
                return 0;
        }
        else{
                int n=0;
                for(n=0; n < dim_each*lines; ++n){
			if(fscanf(file, "%lf", list_points+n) == EOF){
				printf("*err* unable to read %d-th data\n", n);
				return 0;
			}
                }
		/* All the data have been successfully saved */
		printf("readPoints: %d data stored.\n", n);
                fclose(file);
                return 1;
        }
}


/* Return the number of read data */
int dataFromFile(char* file_name, double** target){
        assert(file_name != NULL);
        assert(target != NULL);
        FILE* f = fopen(file_name, "r");
        if(f == NULL){
                fprintf(stderr, "unable to open %s\n", file_name);
                return -1;
        }

        int size = 2; /* size must be bigger that 1, because 
                        * of the if(i==size-1) coparison later */
        double* tmp;
        tmp = realloc(*target, sizeof(double)*size);
        if(tmp != NULL){
                *target = tmp;
                tmp = NULL;
        }
        else{
                fprintf(stderr, "Unable to realloc\n");
		fclose(f);
                return -1;
        }
        
        int i=0;
        while( fscanf(f, "%lf", (*target)+i ) != EOF ){
//                printf("Read: %f\n", (*target)[i]);
                ++i;
                if(i == size-1){
                        size *= 2;
                        tmp = realloc(*target, sizeof(double)*size);
                        if(tmp == NULL){
                                fprintf(stderr, "Unable to realloc, too many points! > %d\n", size);
                                fprintf(stderr, "(but %d have been successfully read)\n", i);
				fclose(f);
                                return i;
                        }
                        else{
                                /* Successful reallocation */
                                *target = tmp;
                                tmp = NULL;
                        }
                }
        }
	fclose(f);
        /* Re-schrank target to precisely i elements, freeing the left memory */
        /* (i is surely smaller that size) */
        tmp = realloc(*target, sizeof(double) * i);
        if(tmp != NULL){
                *target = tmp;
                tmp = NULL;
                printf("dataFromFile: read %d data\n", i);
                return i;
        }
        else{
                printf("Read %d data, but returning array of dimension %d\n", i, size);
                return size;
        }
}
