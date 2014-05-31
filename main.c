#include <mpi.h>
#include <stdio.h>
#include "limits.h"
#define MESSAGE 1

int* route_array;
int* best_route_array;
int** city_matrix;
int start = 0;
int weight_best_route = INT_MAX;
int dimension;
int index = 1;
int myrank;
int size;
int devide_depth = 3;
int devide_depth_index = 0;
int solution_in_thread=0;
MPI_Status* status;
MPI_Request* request;

int** readFileMatrix(char* file_name, int* dimm){
	int i,j;
	FILE* fp;
	int** matrix;
	fp = fopen(file_name,"r");
	if (fp == NULL) {
		perror("Error while opening the file.\n");
		exit(1);
	}
	fscanf(fp, "%i", dimm);
	matrix = (int**) malloc(sizeof(int*) * *dimm);
	for (i = 0; i < *dimm; i++) {
		matrix[i] = (int*) malloc(sizeof(int) * *dimm);
	}
	for (i = 0; i < *dimm; i++) {
		for (j = 0; j < *dimm; j++) {
			fscanf(fp, "%i", &matrix[i][j]);
		}
	}
	return matrix;
}

void freeMainVariables(){
	freeMatrix(city_matrix,dimension);
	free(status);
	free(request);
	free(route_array);
	free(best_route_array);
}


int main(int argc, char* argv[]){
	int i,j;
	double t1,t2;
	MPI_Init(&argc, &argv);

	//t1=MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	status = (MPI_Status*) malloc(sizeof(MPI_Status));
	request = (MPI_Request*) malloc(sizeof(MPI_Request));

	if(myrank ==0){
		city_matrix = readFileMatrix(argv[1], &dimension);
		route_array = (int*) malloc(sizeof(int) * dimension);
		best_route_array = (int*) malloc(sizeof(int) * dimension);

		//distribute matrix over processes
		//first distribute dimension through broadcast
		MPI_Bcast(&dimension, 1, MPI_INT, myrank, MPI_COMM_WORLD);

		for (i = 0; i < dimension; i++) {	//send each row from matrix
			for (j = 1; j < size; j++) {
				MPI_Isend(city_matrix[i], dimension, MPI_INT, j, i, MPI_COMM_WORLD,	status);
			}
		}
	} else {
		MPI_Bcast(&dimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
		route_array = (int*) malloc(sizeof(int) * dimension);
		best_route_array = (int*) malloc(sizeof(int) * dimension);

		//allocate matrix
		city_matrix = (int**) malloc(sizeof(int*) * dimension);
		for (i = 0; i < dimension; i++) {
			city_matrix[i] = (int*) malloc(sizeof(int) * dimension);
		}

		request = malloc(sizeof(MPI_Request) * dimension);

		for (i = 0; i < dimension; i++) {
			MPI_Irecv(city_matrix[i], dimension, MPI_INT, 0, i, MPI_COMM_WORLD, &request[i]);
		}
		MPI_Waitall(dimension, request, NULL);
		free(request);
		request = (MPI_Request*) malloc(sizeof(MPI_Request));	//dit en vorige free vervangen door realloc
	}

	for (i = 0; i < dimension; i++) {
		route_array[i] = i;
	}
	//printf("size %i\n", size);
	//printf("[%2i] done initialising\n", myrank);

	if (myrank == 0){
		greedyMinions();
		//printf("weight best route %i\n",weight_best_route);
		for (i = 0; i < dimension; i++) {
			route_array[i] = i;
		}
		index=1;
		srand (time(NULL));
		for(i=0;i<200;i++)
			minionDominationPathAlgorithm();
		for (i = 0; i < dimension; i++) {
			route_array[i] = i;
		}
		index=1;
	}

	
	WDSPA(0,start);

	MPI_Barrier(MPI_COMM_WORLD);
	ListenToYourNewLeader();

	if(solution_in_thread){
		while(best_route_array[0] != 0){
			int temp = best_route_array[0];
			for(i=0;i<dimension-1;i++){
				best_route_array[i] = best_route_array[i+1];
			}
			best_route_array[dimension-1] = temp;
		}
	}

	if(solution_in_thread && myrank !=0){	//can occur multiple times if different sollutions in different threads have the same weight
		// send
		MPI_Isend(best_route_array, dimension, MPI_INT, 0, 123, MPI_COMM_WORLD,	status);
	} else if (solution_in_thread && myrank == 0) {
		// print
		printf("%i\n",weight_best_route);
		for(i=0;i<dimension;i++){
			printf("%i ", best_route_array[i]);
		}
		printf("0\n");
	} else if (myrank ==0) {	//only takes the first in the receive queue
		//recieve & print
		MPI_Recv(best_route_array, dimension, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("%i\n",weight_best_route);
		for(i=0;i<dimension;i++){
			printf("%i ", best_route_array[i]);
		}
		printf("0\n");
	}

	/*
	t2=MPI_Wtime();
	if (myrank == 0)
		printf("time spent is: %F\n", t2-t1);
	*/

	MPI_Finalize();
	freeMainVariables();
	return 0;
}

void WDSPA(int weight, int city){
	int i,j;
	int temp;
	//check if new upper bound
	ListenToYourNewLeader();

	//printf("weight %i, city %i\n",weight,city);
	if(index==dimension){
		//all cities are traversed
		if((weight+city_matrix[city][start]) < weight_best_route){
			weight_best_route=weight+city_matrix[city][start];
			memcpy(best_route_array,route_array, dimension*sizeof(int));
			// broadcast other processes new upper bound
			PreachWorldDomination();
		}
	} else {
		for(i=index;i<dimension;i++){
			temp = route_array[i];
			//printf("i %i, index %i, temp %i\n",i, index, temp);
			if( city_matrix[city][temp]+weight < weight_best_route){	//check if you need to enter this city
				if (index == devide_depth){					//if you're on the devidedepth, increment devide_depth_index and check if you can continue
					devide_depth_index++;
					if (devide_depth_index%size==myrank){
						//printf("[%i/%i] %i\n",myrank,size,devide_depth_index);
						//xorSwap(&route_array[i],&route_array[index]);
						swap(i,index);
						index++;
						WDSPA(weight+city_matrix[city][temp], temp);
						index--;
						//xorSwap(&route_array[index],&route_array[i]);
						swap(index, i);
					}
				} else {			//you're not on the devide_depth and you want to continue
					swap(i,index);
					//xorSwap(&route_array[i],&route_array[index]);
					index++;
					WDSPA(weight+city_matrix[city][temp], temp);
					index--;
					swap(index, i);
					//xorSwap(&route_array[index],&route_array[i]);
				}
			} else if (index < devide_depth) {		//you're above the devide_depth and you do not want to continue with this city, so you need to increment devide_depth_index appropriately
				temp=1;
				for(j=index;j<devide_depth;j++){
					temp=temp*(dimension-j-1);
				}
				devide_depth_index+=temp;
			} else if (index == devide_depth){
				devide_depth_index++;
			}
		}
	}
}

void greedyMinions(){	// Shortest Path First
	int i,temp,indexsmallestcost;
	weight_best_route=0;
	while(index<dimension){
		//get index smallest cost
		temp=INT_MAX;
		for(i=index;i<dimension;i++){
			if(temp > city_matrix[route_array[index-1]][route_array[i]]){
				temp=city_matrix[route_array[index-1]][route_array[i]];
				indexsmallestcost=i;
			}
		}
		weight_best_route+=temp;
		swap(indexsmallestcost,index);
		index++;
	}
	weight_best_route+=city_matrix[route_array[dimension-1]][0];

	solution_in_thread = 1;
	ListenToYourNewLeader();

	//if our found weight is better than an already found path, broadcast this to other processes
	if (solution_in_thread) {
		memcpy(best_route_array, route_array, dimension * sizeof(int));
		PreachWorldDomination();
	}
}

void printArray(int dimension, int* array){
	int i;
	printf("[");
	for(i=0;i<dimension-1;i++){
		printf("%i, ", array[i]);
	}
	printf("%i]\n", array[dimension-1]);
}

inline void ListenToYourNewLeader(){
	int temp=0;
	int tempnewroute=0;
	MPI_Iprobe(MPI_ANY_SOURCE, MESSAGE, MPI_COMM_WORLD, &temp, status);
	while(temp){
		//printf("recvtry\n");
		MPI_Recv(&tempnewroute, 1, MPI_INT, MPI_ANY_SOURCE, MESSAGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(tempnewroute<weight_best_route){	//only if the sollution is better will this be used
			solution_in_thread=0;			//best solution is not in this thread anymore
			//printf("[%2i] old: %i, ", myrank, weight_best_route);
			weight_best_route = tempnewroute;
			//printf("new: %i\n", weight_best_route);
		}
		MPI_Iprobe(MPI_ANY_SOURCE, MESSAGE, MPI_COMM_WORLD, &temp, status);
	}
}

inline void PreachWorldDomination(){
	int i;
	//printf("[%2i] preaching %i\n",myrank, weight_best_route);
	solution_in_thread=1;									//best solution is in this thread
	for(i=0;i<size;i++){
		if(i!=myrank){
			//int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
			MPI_Isend(&weight_best_route, 1, MPI_INT, i, MESSAGE, MPI_COMM_WORLD, request);
			//printf("sent %i %i\n", weight_best_route, &weight_best_route);
		}
	}
}

inline void swap(int index1, int index2){
	if(index1 != index2){
		int temp = route_array[index1];
		route_array[index1] = route_array[index2];
		route_array[index2]=temp;
	}
}

inline void xorSwap (int *y, int *x) {
	if (x != y) {
		*x ^= *y;
		*y ^= *x;
		*x ^= *y;
	}
}

void freeMatrix(int** matrix, int dimension){
	int i;
	for(i=0;i<dimension;i++){
		free(matrix[i]);
	}
	free(matrix);
}


void minionDominationPathAlgorithm(){	//2opt
	int i, j, improvement, biggest_improvement, swap_i, swap_j,weight,temp;

	//create random path
	for (i = 0; i < dimension; i++) {
		index = (rand() % (dimension - i)) + i;
		temp = route_array[index];
		route_array[index] = route_array[i];
		route_array[i] = temp;
	}

	//2-opt algorithm on random path choosing/applying best improvement each iteration until there are no more improvements
	do{
		biggest_improvement=0;
		for(i=0;i<dimension-2;i++){
			for(j=i+2;j<dimension;j++){
				improvement = getImprovement(i,j);
				if (improvement < biggest_improvement){
					biggest_improvement = improvement;
					swap_i = i;
					swap_j = j;
				}
			}
		}

		if (biggest_improvement<0){		//only swap if there is an improvement
			SwapMinionsSwap(swap_i,swap_j);
		}
	} while (biggest_improvement<0);

	//calculate weight of found path and check if this is better than an already found path
	weight = calculateWeight();
	if (weight < weight_best_route){	//if weight is better than the current best weight see if there has been broadcast a better boundry if not broadcast my boundry
		weight_best_route = weight;
		solution_in_thread=1;
		ListenToYourNewLeader();
		//printf("weight: %i    bestweight: %i\n",weight, weight_best_route);

		//if our found weight is better than an already found path, broadcast this to other processes
		if(solution_in_thread){
			//printf("2opt sent weight\n");
			memcpy(best_route_array,route_array, dimension*sizeof(int));
			PreachWorldDomination();
		}
	}
}

int calculateWeight() {
	int weight = 0;
	int i;
	for (i = 0; i < dimension - 1; i++) {
		weight += city_matrix[route_array[i]][route_array[(i + 1)]];
	}
	return weight + city_matrix[route_array[dimension - 1]][route_array[0]];
}

void SwapMinionsSwap(int i, int j) {
	int k;
	i++;
	swap(i, j);
	for (k = 1; k <= ((j - i) / 2); k++) {
		swap(i + k, j - k);
	}
}

int getImprovement(int i, int j) {
	return city_matrix[route_array[i]][route_array[j]]
			+ city_matrix[route_array[(i + 1) % dimension]][route_array[(j + 1) % dimension]]
			- city_matrix[route_array[i]][route_array[(i + 1) % dimension]]
			- city_matrix[route_array[j]][route_array[(j + 1) % dimension]];
}

