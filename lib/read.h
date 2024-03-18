#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../lib/atom.h"
#ifndef read_h
#define read_h
#endif
int length(FILE *ipt){
int rnr;
    while (!feof(ipt)) 
        {
            fscanf(ipt, "%*s %*d %*s %*c %*s %d %*f %*f %*f %*f %*f %*c",&rnr);
        }
    if(rnr==0)
    {
        exit(1);
    }
    return rnr;
}
int * positionO2(int length, FILE *ipt){
    int *positions = malloc(sizeof(int) * length);
    while (!feof(ipt)) 
    {
        char type[10]; int rnr, anr,pos;
        if(fscanf(ipt, "%*s %d %s %*c %*s %d %*f %*f %*f %*f %*f %*c", &anr, type, &rnr)>0)
        {
            if(strcmp(type,"O2")==0)
            {
                pos=rnr-1;
                *(positions+pos)=anr;
            }
        }
    }
    return (positions);
}
int * positionC5(int length, FILE *ipt){
    int *positions = malloc(sizeof(int) * length);
    while (!feof(ipt)) 
    {
        char type[10]; int rnr, anr,pos;
        if(fscanf(ipt, "%*s %d %s %*c %*s %d %*f %*f %*f %*f %*f %*c", &anr, type, &rnr)>0)
        {
            if(strcmp(type,"C5")==0)
            {
                pos=rnr-1;
                *(positions+pos)=anr;
            }
        }
    }
    return (positions);
}
