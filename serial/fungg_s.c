/*
CA - OpenMP
fungg_s.c
Routines used in gengroups_s.c program
*/

#include <math.h>
#include <float.h>
#include "../shared/definegg.h" // definition of constants
#include <stdio.h>

void swap(float *xp, float *yp);
void bubbleSort(float arr[], int n);

/* 1 - Function to calculate the genetic distance; Euclidean distance between two elements.
       Input:   two elements of NFEAT characteristics (by reference)
       Output:  distance (double)
***************************************************************************************************/
double geneticdistance(float *elem1, float *elem2)
{
   double distance = 0;
   for (int i = 0; i < NFEAT; i++)
      // Get euclidean distance of specific feature
      distance += pow(elem1[i] - elem2[i], 2);

   return sqrt(distance);
}

/* 2 - Function to calculate the closest group (closest centroid) for each element.
   Input:   nelems   number of elements, int
            elems    matrix, with the information of the elements, of size MAXELE x NFEAT, by reference
            cent    matrix, with the centroids, of size NGROUPS x NFEAT, by reference
   Output:  grind   vector of size MAXELE, by reference, closest group for each element
***************************************************************************************************/
void closestgroup(int nelems, float elems[][NFEAT], float cent[][NFEAT], int *grind)
{
   float min_d; // Auxiliary variable to store the closest centroid value
   int min_d_i; // Auxiliary variable to store the closest centroid index
   float aux_d; // Auxiliary variable to store the output of geneticdistance

   // Iterave over all elements
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
   Input:  elems     elements (matrix of size MAXELE x NFEAT, by reference)
           iingrs   indices of the elements in each group (vector of size NGROUPS with information for each group)
   Output: compact  compactness of each group (vector of size NGROUPS, by reference) 
***************************************************************************************************/
void groupcompactness(float elems[][NFEAT], struct ginfo *iingrs, float *compact)
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
   Input:  iingrs   indices of the elements in each group (matrix of size NGROUPS x MAXELE, by reference)
           dise     information about the diseases (NGROUPS x TDISEASE)
   Output: disepro  analysis of the diseases: maximum, minimum of the medians and groups
***************************************************************************************************/
void diseases(struct ginfo *iingrs, float dise[][TDISEASE], struct analysis *disepro)
{
   float diseaseList[MAXELE];
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

      if (gsize > 0)
      {
         for (int j = 0; j < TDISEASE; j++)
         {
            // Fill the auxiliary array with all the value for current disease of current group
            for (int k = 0; k < gsize; k++)
               diseaseList[k] = dise[iingrs[i].members[k]][j];

            bubbleSort(diseaseList, gsize);

            // Get the median
            if (gsize % 2 == 0)
               median = diseaseList[(gsize + 1) / 2];
            else
               median = diseaseList[(gsize) / 2];

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
   }
}

void swap(float *xp, float *yp)
{
   float temp = *xp;
   *xp = *yp;
   *yp = temp;
}

// A function to implement bubble sort
void bubbleSort(float arr[], int n)
{
   int i, j;
   for (i = 0; i < n - 1; i++)

      // Last i elements are already in place
      for (j = 0; j < n - i - 1; j++)
         if (arr[j] > arr[j + 1])
            swap(&arr[j], &arr[j + 1]);
}
