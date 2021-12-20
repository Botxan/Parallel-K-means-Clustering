/*
   fungg.h
   headers of the functions used in gengroups_s.c
***************************************************************/

extern double geneticdistance(float *elem1, float *elem2);
extern void closestgroup(int nelem, float **elem, float cent[][NFEAT], int *grind);
extern void groupcompactness(float **elem, struct ginfo *iingrs, float *compact);
extern void diseases(int nelems, struct ginfo *iingrs, float **dise, struct analysis *disepro);
