/* 
    CA - practical work OpenMP
    gengroups_s.c SERIAL VERSION

    Processing genetic characteristics to discover information about diseases
    Classify in NGROUPS groups, elements of NFEAT features, according to "distances"    

    Input:  dbgen.dat 	   input file with genetic information
            dbdise.dat     input file with information about diseases
    Output: results_s.out  centroids, number of group members and compactness, and diseases

    Compile with module fungg_s.c and include option -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "../shared/definegg.h"
#include "../shared/fungg.h"

float **elems;				  // matrix to keep information about every element
struct ginfo iingrs[NGROUPS]; // vector to store information about each group: members and size

float **dise;					   // probabilities of diseases (from dbdise.dat)
struct analysis disepro[TDISEASE]; // vector to store information about each disease (max, min, group...)

// Main program
// ============
void main(int argc, char *argv[])
{
	float cent[NGROUPS][NFEAT], newcent[NGROUPS][NFEAT]; // centroids and new centroids
	double additions[NGROUPS][NFEAT + 1];
	float compact[NGROUPS]; // compactness of each group or cluster

	int i, j;
	int nelems, group;
	int *grind; // group assigned to each element
	int finish = 0, niter = 0;
	double discent;

	FILE *f1, *f2;
	struct timespec t1, t2, t3, t4, t5, t6, t7;
	double t_read, t_clus, t_org, t_compact, t_anal, t_write;

	if ((argc < 3) || (argc > 4))
	{
		printf("ATTENTION:  progr file1 (elems) file2 (dise) [num elems])\n");
		exit(-1);
	}

	printf("\n >> Serial execution\n");
	clock_gettime(CLOCK_REALTIME, &t1);

	// read data from files: elems[i][j] and dise[i][j]
	// ===============================================
	f1 = fopen(argv[1], "r");
	if (f1 == NULL)
	{
		printf("Error opening file %s \n", argv[1]);
		exit(-1);
	}

	fscanf(f1, "%d", &nelems);
	if (argc == 4)
		nelems = atoi(argv[3]);

	// Assign memory dynamically to elems, dise and grind
	elems = (float **)malloc(nelems * sizeof(float *));
	dise = (float **)malloc(nelems * sizeof(float *));
	grind = (int *)malloc(nelems * sizeof(int));
	for (i = 0; i < nelems; i++)
	{
		elems[i] = (float *)malloc(NFEAT * sizeof(float));
		dise[i] = (float *)malloc(TDISEASE * sizeof(float));
	}

	for (i = 0; i < nelems; i++)
		for (j = 0; j < NFEAT; j++)
			fscanf(f1, "%f", &(elems[i][j]));

	fclose(f1);

	f1 = fopen(argv[2], "r");
	if (f1 == NULL)
	{
		printf("Error opening file %s \n", argv[1]);
		exit(-1);
	}

	for (i = 0; i < nelems; i++)
		for (j = 0; j < TDISEASE; j++)
			fscanf(f1, "%f", &dise[i][j]);

	fclose(f1);

	clock_gettime(CLOCK_REALTIME, &t2);

	// select randomly the first centroids
	// ===================================
	srand(147);
	for (i = 0; i < NGROUPS; i++)
		for (j = 0; j < NFEAT / 2; j++)
		{
			cent[i][j] = (rand() % 10000) / 100.0;
			cent[i][j + NFEAT / 2] = cent[i][j];
		}

	// Phase 1: classify elements and calculate new centroids
	// ======================================================
	niter = 0;
	finish = 0;
	while ((finish == 0) && (niter < MAXIT))
	{

		// Obtain the closest group or cluster for each element
		closestgroup(nelems, elems, cent, grind);

		// calculate new centroids for each group:  average of each dimension or feature
		// additions: to accumulate the values for each feature and cluster. Last value: number of elements in the group
		for (i = 0; i < NGROUPS; i++)
			for (j = 0; j < NFEAT + 1; j++)
				additions[i][j] = 0.0;

		for (i = 0; i < nelems; i++)
		{
			for (j = 0; j < NFEAT; j++)
				additions[grind[i]][j] += elems[i][j];
			additions[grind[i]][NFEAT]++;
		}

		// Calculate new centroids and decide to finish or not depending on DELTA
		finish = 1;
		for (i = 0; i < NGROUPS; i++)
		{
			if (additions[i][NFEAT] > 0)
			{ // the group is not empty
				for (j = 0; j < NFEAT; j++)
					newcent[i][j] = additions[i][j] / additions[i][NFEAT];

				// decide if the process needs to be finished
				discent = geneticdistance(&newcent[i][0], &cent[i][0]);
				if (discent > DELTA)
					finish = 0; // there is change at least in one of the dimensions; continue with the process

				// copy new centroids
				for (j = 0; j < NFEAT; j++)
					cent[i][j] = newcent[i][j];
			}
		}

		niter++;
	} // while

	clock_gettime(CLOCK_REALTIME, &t3);

	// Phase 2: count the number of elements of each group and calculate the "compactness" of the group
	// and analyse diseases
	// ================================================================================================
	for (i = 0; i < NGROUPS; i++)
		iingrs[i].size = 0;

	// number of elements and classification
	for (i = 0; i < nelems; i++)
	{
		group = grind[i];
		iingrs[group].members[iingrs[group].size] = i;
		iingrs[group].size++;
	}

	// free the memory
	free(grind);

	clock_gettime(CLOCK_REALTIME, &t4);

	// compactness of each group: average distance between elements
	groupcompactness(elems, iingrs, compact);

	clock_gettime(CLOCK_REALTIME, &t5);

	// diseases analysis

	diseases(nelems, iingrs, dise, disepro);

	// Free the memory
	free(elems);
	free(dise);

	clock_gettime(CLOCK_REALTIME, &t6);

	// write results in a file
	// =======================

	f2 = fopen("results_s.out", "w");
	if (f2 == NULL)
	{
		printf("Error when opening file results_s.outs \n");
		exit(-1);
	}

	fprintf(f2, " Centroids of groups \n\n");
	for (i = 0; i < NGROUPS; i++)
	{
		for (j = 0; j < NFEAT; j++)
			fprintf(f2, "%7.3f", cent[i][j]);
		fprintf(f2, "\n");
	}

	fprintf(f2, "\n >> Size of the groups \n\n");
	for (i = 0; i < NGROUPS / 10; i++)
	{
		for (j = 0; j < 10; j++)
			fprintf(f2, "%6d", iingrs[10 * i + j].size);
		fprintf(f2, "\n");
	}

	fprintf(f2, "\n >> Group compactness \n\n");
	for (i = 0; i < NGROUPS / 10; i++)
	{
		for (j = 0; j < 10; j++)
			fprintf(f2, "%9.2f", compact[10 * i + j]);
		fprintf(f2, "\n");
	}

	fprintf(f2, "\n\n Analysis of deseases (medians)\n\n");
	fprintf(f2, "\n Dise.  M_max - Group   M_min - Group");
	fprintf(f2, "\n ==================================\n");
	for (i = 0; i < TDISEASE; i++)
		fprintf(f2, "  %2d     %4.2f - %2d      %4.2f - %2d\n", i, disepro[i].mmax,
				disepro[i].gmax, disepro[i].mmin, disepro[i].gmin);

	fclose(f2);

	clock_gettime(CLOCK_REALTIME, &t7);

	// print some results
	// ===================

	t_read = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec) / (double)1e9;
	t_clus = (t3.tv_sec - t2.tv_sec) + (t3.tv_nsec - t2.tv_nsec) / (double)1e9;
	t_org = (t4.tv_sec - t3.tv_sec) + (t4.tv_nsec - t3.tv_nsec) / (double)1e9;
	t_compact = (t5.tv_sec - t4.tv_sec) + (t5.tv_nsec - t4.tv_nsec) / (double)1e9;
	t_anal = (t6.tv_sec - t5.tv_sec) + (t6.tv_nsec - t5.tv_nsec) / (double)1e9;
	t_write = (t7.tv_sec - t6.tv_sec) + (t7.tv_nsec - t6.tv_nsec) / (double)1e9;

	printf("\n    Number of iterations: %d", niter);
	printf("\n    T_read:    %6.3f s", t_read);
	printf("\n    T_clus:    %6.3f s", t_clus);
	printf("\n    T_org:     %6.3f s", t_org);
	printf("\n    T_compact: %6.3f s", t_compact);
	printf("\n    T_anal:    %6.3f s", t_anal);
	printf("\n    T_write:   %6.3f s", t_write);
	printf("\n    ========================");
	printf("\n    T_total:  %6.3f s\n\n", t_read + t_clus + t_org + t_compact + t_anal + t_write);

	printf("\n centroids 0, 40 and 80 and the compactness of their group\n ");
	for (i = 0; i < NGROUPS; i += 40)
	{
		printf("\n  z%2d -- ", i);
		for (j = 0; j < NFEAT; j++)
			printf("%5.1f", cent[i][j]);
		printf("\n          %5.6f\n", compact[i]);
	}

	printf("\n >> Size of the groups \n\n");
	for (i = 0; i < NGROUPS / 10; i++)
	{
		for (j = 0; j < 10; j++)
			printf("%6d", iingrs[10 * i + j].size);
		printf("\n");
	}

	printf("\n >> Group compactness \n\n");
	for (i = 0; i < NGROUPS / 10; i++)
	{
		for (j = 0; j < 10; j++)
			printf("%9.2f", compact[10 * i + j]);
		printf("\n");
	}

	printf("\n\n Analysis of diseases (medians)\n\n");
	printf("\n Dise.  M_max - Group   M_min - Group");
	printf("\n ==================================\n");
	for (i = 0; i < TDISEASE; i++)
		printf("  %2d     %4.2f - %2d      %4.2f - %2d\n", i, disepro[i].mmax,
			   disepro[i].gmax, disepro[i].mmin, disepro[i].gmin);

	printf("\n");
}
