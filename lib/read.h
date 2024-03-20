#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef read_h
#define read_h
#endif
int length(FILE *ipt){
rewind(ipt);
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
    rewind(ipt);
    int *positionO2 = malloc(sizeof(int) * length);
    while (!feof(ipt)) 
    {
        char type[10]; int rnr, anr,pos;
        if(fscanf(ipt, "%*s %d %s %*c %*s %d %*f %*f %*f %*f %*f %*c", &anr, type, &rnr)>0)
        {
            if(strcmp(type,"O2")==0)
            {
                pos=rnr-1;
                *(positionO2+pos)=anr;
            }
        }
    }
    return (positionO2);
}
int * positionC5(int length, FILE *ipt){
    rewind(ipt);
    int *positionC5 = malloc(sizeof(int) * length);
    while (!feof(ipt)) 
    {
        char type[10]; int rnr, anr,pos;
        if(fscanf(ipt, "%*s %d %s %*c %*s %d %*f %*f %*f %*f %*f %*c", &anr, type, &rnr)>0)
        {
            if(strcmp(type,"C5")==0)
            {
                pos=rnr-1;
                *(positionC5+pos)=anr;
            }
        }
    }
    return (positionC5);
}
float * getpos(FILE *ipt, int pos){
rewind(ipt);
float * posit=malloc(sizeof(float)*3);
while (!feof(ipt)) 
    {
        char type[10];int anr;
        if(fscanf(ipt, "%*s %d %*s %*c %*s %*d %f %f %f %*f %*f %*c",&anr, (posit+0),(posit+1),(posit+2))>0)
        {
            if(anr==pos){
                break;
            }
        }
    }
    return ( posit );
}


