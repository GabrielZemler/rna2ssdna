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
float *getangle(float* vec1, float *vec2)
{
    float *angle=(float *)malloc(sizeof(float));
    float mag_vec1=*(vec1+0)**(vec1+1)**(vec1+2);
    float mag_vec2=sqrtf(powf(*(vec2+0),2)+powf(*(vec2+1),2)+powf(*(vec2+2),2));
    float abs_scal=*(vec1+0)**(vec2+0)+*(vec1+1)**(vec2+1)+*(vec1+2)**(vec2+2);
    if(abs_scal<0){
        abs_scal=-abs_scal;
    }
    *angle=acosf(abs_scal/(mag_vec1+mag_vec2));
    return ( angle );
}
float *rotate(float * vec, float alpha, float beta, float gamma){
    float * rot=(float *)malloc(sizeof(float)*3);
    float Rxyz[3][3];
    Rxyz[0][0]=cosf(alpha)*cosf(beta);
    Rxyz[0][1]=cosf(alpha)*sinf(beta)*sinf(gamma)-sinf(alpha)*sinf(beta);
    Rxyz[0][2]=cosf(alpha)*sinf(beta)*cosf(gamma)+sinf(alpha)*sinf(gamma);
    Rxyz[1][0]=sinf(alpha)*cosf(beta);
    Rxyz[1][1]=sinf(alpha)*sinf(beta)*sinf(gamma)+cosf(alpha)*cosf(gamma);
    Rxyz[1][2]=sinf(alpha)*sinf(beta)*cosf(gamma)-cosf(alpha)*sinf(gamma);
    Rxyz[2][0]=-sinf(beta);
    Rxyz[2][1]=cosf(beta)*sinf(gamma);
    Rxyz[2][2]=cosf(beta)*cosf(gamma);
    printf("%f %f %f\n",alpha, beta, gamma );
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%f ", Rxyz[i][j]);
        }
        printf("\n");
    }
    
    for (int i = 0; i < 3; i++)
    {
       for (int j = 0; j < 3; j++)
       {
       *(rot+i)=*(rot+i)+*(vec+j)*Rxyz[i][j];
       }
       
    }
    return (rot);
}