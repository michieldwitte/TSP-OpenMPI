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
int devide_depth = 2;
int devide_depth_index = 0;
int solution_in_thread=0;

/*
 * heuristiek uitleggen adhv het boek (2opt) en hoe dit is opgebouwd ect
 * uitvoeringstijden zonder en met heuristieken
 * shortest path vergelijken met 2opt in snelheid en in oplossing
 */

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
	//free(status);
	//free(request);
	free(route_array);
	free(best_route_array);
}


// MPI_STATUS_IGNORE	als status gebruiken als niets nodig
// isend en recv voor matrix door te sturen vervangen door bcast ?
int main(int argc, char* argv[]){
	int i,j;
	//double t1,t2;
	//t1=time(0);
	city_matrix = readFileMatrix(argv[1], &dimension);
	route_array = (int*) malloc(sizeof(int) * dimension);
	best_route_array = (int*) malloc(sizeof(int) * dimension);

	for (i = 0; i < dimension; i++) {
		route_array[i] = i;
	}
	
	zoekzoek();
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
		
	WDSPA(0,start);


	printf("%i\n", weight_best_route);

	while(best_route_array[0] != 0){
		int temp = best_route_array[0];
		for(i=0;i<dimension-1;i++){
			best_route_array[i] = best_route_array[i+1];
		}
		best_route_array[dimension-1] = temp;
	}
	for(i=0;i<dimension;i++)
		printf("%i ", best_route_array[i]);
	printf("\n");
	//printArray(dimension,best_route_array);

	//t2=time(0);
	//printf("time spent is: %F\n", t2-t1);
	
	freeMainVariables();
	return 0;
}

void WDSPA(int weight, int city){
	int i,j;
	int temp;
	//check if new upper bound

	//printf("weight %i, city %i\n",weight,city);
	if(index==dimension){
		//alle nodes zijn traversed
		if((weight+city_matrix[city][start]) < weight_best_route){
			weight_best_route=weight+city_matrix[city][start];
			memcpy(best_route_array,route_array, dimension*sizeof(int));
			// broadcast other processes new upper bound
		}
	} else {
		for(i=index;i<dimension;i++){
			temp = route_array[i];
			//printf("i %i, index %i, temp %i\n",i, index, temp);
			if( city_matrix[city][temp]+weight < weight_best_route){	//kijken als je deze node moet binnen gaan
				
						//printf("[%i/%i] %i\n",myrank,size,devide_depth_index);
						//xorSwap(&route_array[i],&route_array[index]);
						swap(i,index);
						index++;
						WDSPA(weight+city_matrix[city][temp], temp);
						index--;
						//xorSwap(&route_array[index],&route_array[i]);
						swap(index, i);		//??????? dit zorgt ervoor dat de lijst bij elk process in dezelfde volgorde staat en dus kan zo geteld worden welke nodes nog moeten gedaan worden ??????
											// {zet in verslag} bij elke node op die diepte teller++, dit houdt ook in dat je niet mag bounden voor de splitsdiepte
											// of hier net na een elseif (index<splitsdiepte){teller=+(index*......)} daje toch alles kan overslaan	
					}
		}
	}
}

void zoekzoek(){	//beta-code
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
	memcpy(best_route_array, route_array, dimension * sizeof(int));
}

void printArray(int dimension, int* array){
	int i;
	printf("[");
	for(i=0;i<dimension-1;i++){
		printf("%i, ", array[i]);
	}
	printf("%i]\n", array[dimension-1]);
}

inline void swap(index1, index2){
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
		memcpy(best_route_array,route_array, dimension*sizeof(int));
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

void SwapMinionsSwap( i, j) {
	int k;
	i++;
	swap(i, j);
	for (k = 1; k <= ((j - i) / 2); k++) {
		swap(i + k, j - k);
	}
}

int getImprovement( i, j) {
	return city_matrix[route_array[i]][route_array[j]]
			+ city_matrix[route_array[(i + 1) % dimension]][route_array[(j + 1) % dimension]]
			- city_matrix[route_array[i]][route_array[(i + 1) % dimension]]
			- city_matrix[route_array[j]][route_array[(j + 1) % dimension]];
}

