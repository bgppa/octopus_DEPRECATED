/* Header for fileio.h, file interaction **to improve**!!! */
#ifndef _FILEIO_H_
#define _FILEIO_H_

/* Read a matrix of numbers from a file.
 * Arguments are:
 * file_name, then a target array that will contain the data
 * (so list_points has to be at leat of dimension dime_each times lines)
 * dim_each : how many numbers are in each line
 * lines : how many lines */
int readPoints(char *file_name, double *list_points, int dim_each, int lines);

/* Return the number of read data */
int dataFromFile(char *file_name, double **target);
#endif
