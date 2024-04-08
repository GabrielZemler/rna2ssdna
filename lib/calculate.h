#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef calculate_h
#define calculate_h
#endif
float *getvec(float * vec1, float *vec2)
{
    float *vec=(float *)malloc(sizeof(float) * 3);
    *(vec+0)=*(vec2+0)-*(vec1+0);
    *(vec+1)=*(vec2+1)-*(vec1+1);
    *(vec+2)=*(vec2+2)-*(vec1+2);
    return ( vec );
}
//SPHERICAL CALCULATION

float phi(float* vec)
{
    float phi;
    float mag_vec1=sqrtf(powf(*(vec+0),2)+powf(*(vec+1),2));
    float mag_vec2=1;
    float abs_scal=*(vec+0);
    if(abs_scal<0){
        abs_scal=-abs_scal;
    }
    phi=acosf(abs_scal/(mag_vec1*mag_vec2));
    if(*(vec+0)<0 &&*(vec+1)>0)
    {
        phi=3.142-phi;
    }
    else if(*(vec+0)<0 &&*(vec+1)<0)
    {
        phi=4.712-phi;
    }
    else if(*(vec+0)>0 &&*(vec+1)<0){
        phi=6.283-phi;
    }
    return phi;
}
float theta(float* vec)
{
    float theta;
    float mag_vec1=sqrtf(powf(*(vec+0),2)+powf(*(vec+1),2)+powf(*(vec+2),2));
    //printf("mag:%f\n",mag_vec1);
    float mag_vec2=1;
    float abs_scal=*(vec+2);
    if(abs_scal<0){
        abs_scal=-abs_scal;
    }
    theta=acosf(abs_scal/(mag_vec1*mag_vec2));
    //printf("theta%f\n",theta);
    if(*(vec+2)<0){
        theta=3.142-theta;
    }
    return theta;
}
float r(float* vec)
{
    float r;
    r=sqrtf(powf(*(vec+0),2)+powf(*(vec+1),2)+powf(*(vec+2),2));
    return r;
}
float * revert(float r, float theta, float phi)
{
    float * xyz=(float *)malloc(sizeof(float)*3);
    *(xyz+0)=r*sinf(theta)*cosf(phi);
    *(xyz+1)=r*sinf(theta)*sinf(phi);
    *(xyz+2)=r*cos(theta);
    return (xyz);
}