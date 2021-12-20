/*
CA - OpenMP
fungg_s.c
Routines used in gengroups_s.c program
*/

#include <math.h>
#include <float.h>
#include "../shared/definegg.h" // definition of constants
#include <stdio.h>
#include <omp.h>

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(float arr[], int l, int m, int r);

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(float arr[], int l, int r);

/* 1 - Function to calculate the genetic distance; Euclidean distance between two elements.
       Input:   two elements of NFEAT characteristics (by reference)
       Output:  distance (double)
***************************************************************************************************/
double geneticdistance(float *elem1, float *elem2)
{
   double distance = 0;

   #pragma omp parallel for reduction(+:distance)
   for (int i = 0; i < NFEAT; i++)
      // Get euclidean distance of specific feature
      distance += pow(elem1[i] - elem2[i], 2);

   return sqrt(distance);
}

/* 2 - Function to calculate the closest group (closest centroid) for each element.
   Input:   nelems   number of elements, int
            elems    matrix, with the information of the elements, of size nelems x NFEAT, by reference
            cent    matrix, with the centroids, of size NGROUPS x NFEAT, by reference
   Output:  grind   vector of size nelems, by reference, closest group for each element
***************************************************************************************************/
void closestgroup(int nelems, float **elems, float cent[][NFEAT], int *grind)
{
   double min_d; // Auxiliary variable to store the closest centroid value
   int min_d_i;  // Auxiliary variable to store the closest centroid index
   double aux_d; // Auxiliary variable to store the output of geneticdistance

   // Iterate over all elements
   for (int i = 0; i < nelems; i++)
   {
      // Initialize the minimum distance
      min_d = FLT_MAX;
      min_d_i = -1;

      // For each element, get the distance with every centroid, and store the
      // index of the centroid with closest distance
      for (int j = 0; j < NGROUPS; j++)
      {
         aux_d = geneticdistance(elems[i], cent[j]);
         if (aux_d < min_d)
         {
            min_d = aux_d;
            min_d_i = j;
         }
      }

      // Assign to the element i the closest centroid, stored in min
      grind[i] = min_d_i;
   }
}

/* 3 - Function to calculate the compactness of each group (average distance between all the elements in the group) 
   Input:  elems     elements (matrix of size nelems x NFEAT, by reference)
           iingrs   indices of the elements in each group (vector of size NGROUPS with information for each group)
   Output: compact  compactness of each group (vector of size NGROUPS, by reference) 
***************************************************************************************************/
void groupcompactness(float **elems, struct ginfo *iingrs, float *compact)
{
   // We need this variable because compact variable points to an average of distances
   // if we try to calculate the sum of all distances in this variable, if the group is so big,
   // it may not fit and the compiler will try to round it, so the final result is not going to be
   // completely accurate
   double comp_aux;
   int gsize;

   // Iterate over each group
   for (int i = 0; i < NGROUPS; i++)
   {
      gsize = iingrs[i].size;
      if (gsize <= 1)
         compact[i] = 0.0;
      else
      {
         comp_aux = 0.0;
         // Iterate over each element of the group
         for (int j = 0; j < gsize; j++)
            // Get the distance of the element j with respect the rest of the elements of the group
            for (int k = j + 1; k < gsize; k++)
               comp_aux += geneticdistance(elems[iingrs[i].members[j]], elems[iingrs[i].members[k]]);
         compact[i] = (double)(comp_aux / ((gsize * (gsize - 1)) / 2));
      }
   }
}

/* 4 - Function to analyse diseases 
   Input:  iingrs   indices of the elements in each group (matrix of size NGROUPS x nelems, by reference)
           dise     information about the diseases (NGROUPS x TDISEASE)
   Output: disepro  analysis of the diseases: maximum, minimum of the medians and groups
***************************************************************************************************/
void diseases(int nelems, struct ginfo *iingrs, float **dise, struct analysis *disepro)
{
	float *diseaseList;
	float median;
	int gsize;

	// Intialize disepro struct for the current disease
	for (int i = 0; i < TDISEASE; i++)
	{
		disepro[i].mmax = FLT_MIN;
		disepro[i].mmin = FLT_MAX;
	}

	for (int i = 0; i < NGROUPS; i++)
	{
		gsize = iingrs[i].size;

		// Allocate memory for diseaseList
		diseaseList = (float *)malloc(gsize * sizeof(float));

		if (gsize > 0)
		{
			for (int j = 0; j < TDISEASE; j++)
			{
				// Fill the auxiliary array with all the value for current disease of current group
				for (int k = 0; k < gsize; k++)
					diseaseList[k] = dise[iingrs[i].members[k]][j];

				mergeSort(diseaseList, 0, gsize - 1);

				// Get the median
				if (gsize % 2 == 0) median = diseaseList[(gsize + 1) / 2];
				else median = diseaseList[(gsize) / 2];

				// Check if it is the new maximum / minimum of the current disease´
				if (median < disepro[j].mmin)
				{
					disepro[j].mmin = median;
					disepro[j].gmin = i;
				}
				else if (median > disepro[j].mmax)
				{
					disepro[j].mmax = median;
					disepro[j].gmax = i;
				}
			}
		}

		// Free the memory
		free(diseaseList);
	}
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(float arr[], int l, int m, int r)
{
   int i, j, k;
   int n1 = m - l + 1;
   int n2 = r - m;

   /* create temp arrays */
   float L[n1], R[n2];

   /* Copy data to temp arrays L[] and R[] */
   for (i = 0; i < n1; i++)
      L[i] = arr[l + i];
   for (j = 0; j < n2; j++)
      R[j] = arr[m + 1 + j];

   /* Merge the temp arrays back into arr[l..r]*/
   i = 0; // Initial index of first subarray
   j = 0; // Initial index of second subarray
   k = l; // Initial index of merged subarray
   while (i < n1 && j < n2)
   {
      if (L[i] <= R[j])
      {
         arr[k] = L[i];
         i++;
      }
      else
      {
         arr[k] = R[j];
         j++;
      }
      k++;
   }

   /* Copy the remaining elements of L[], if there
    are any */
   while (i < n1)
   {
      arr[k] = L[i];
      i++;
      k++;
   }

   /* Copy the remaining elements of R[], if there
    are any */
   while (j < n2)
   {
      arr[k] = R[j];
      j++;
      k++;
   }
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(float arr[], int l, int r)
{
   if (l < r)
   {
      // Same as (l+r)/2, but avoids overflow for
      // large l and h
      int m = l + (r - l) / 2;

      // Sort first and second halves
      mergeSort(arr, l, m);
      mergeSort(arr, m + 1, r);

      merge(arr, l, m, r);
   }
}