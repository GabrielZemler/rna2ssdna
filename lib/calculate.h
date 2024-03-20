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
    float mag_vec1=sqrtf(powf(*(vec1+0),2)+powf(*(vec1+1),2)+powf(*(vec1+2),2));
    float mag_vec2=sqrtf(powf(*(vec2+0),2)+powf(*(vec2+1),2)+powf(*(vec2+2),2));
    float abs_scal=*(vec1+0)**(vec2+0)+*(vec1+1)**(vec2+1)+*(vec1+2)**(vec2+2);
    if(abs_scal<0){
        abs_scal=-abs_scal;
    }
    *angle=acosf(abs_scal/(mag_vec1+mag_vec2));
    return ( angle );
}
float * allangles(float* vec){
    float **abg=(float**)malloc(sizeof(float)*3);
    float v1[3]={1.00,0.00,0.00};
    float v2[3]={0.00,1.00,0.00}; 
    float v3[3]={0.00,0.00,1.00};

    return *((abg));
}
float *rotatexyz(float * vec, float alpha, float beta, float gamma){
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
    //printf("%f %f %f\n",alpha, beta, gamma );
    for (int i = 0; i < 3; i++)
    {
       for (int j = 0; j < 3; j++)
       {
       *(rot+i)=*(rot+i)+*(vec+j)*Rxyz[i][j];
       }
       
    }
    return (rot);
}
float roll(float * vec)
{
    //vec_y=0 vecr=100
    float mag_nvec=sqrtf(powf(*(vec+0),2)+powf(*(vec+2),2));
    float abs_scal=sqrtf(powf(*(vec+0),2));
    float roll=0.0;
    if(mag_nvec!=0)
    {
        roll= acosf(abs_scal/(mag_nvec));
    }
    return roll;
}
float pitch(float * vec)
{
    //vec_x=0 vecr=0 0 1
    float mag_nvec=sqrtf(powf(*(vec+1),2)+powf(*(vec+2),2));
    float abs_scal=sqrtf(powf(*(vec+2),2));
    float pitch=0.0;
    if(mag_nvec!=0)
    {
        pitch= acosf(abs_scal/(mag_nvec));
    }
    return pitch;
}
float yaw(float * vec)
{
    //vec_z=0 vecr=0 1 0
    float mag_nvec=sqrtf(powf(*(vec+0),2)+powf(*(vec+1),2));
    float abs_scal=sqrtf(powf(*(vec+1),2));
    float yaw=0.0;
    if(mag_nvec!=0)
    {
        yaw= acosf(abs_scal/(mag_nvec));
    }
    return yaw;
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
    printf("used values: r=%.3f theta=%.3f phi=%.3f\n",r,theta,phi);
    float * xyz=(float *)malloc(sizeof(float)*3);
    *(xyz+0)=r*sinf(theta)*cosf(phi);
    *(xyz+1)=r*sinf(theta)*sinf(phi);
    *(xyz+2)=r*cos(theta);
    return (xyz);
}
float * rotate(float * xyz, float theta, float phi)
{
    float * nxyz= (float *)malloc(sizeof(float)*3);
    *(nxyz+0)=*(xyz+0)*cosf(phi)*cosf(theta)-*(xyz+1)*sinf(theta)+*(xyz+2)*sinf(phi)*cosf(theta);
    *(nxyz+1)=*(xyz+0)*sinf(theta)*cosf(phi)+*(xyz+1)*cosf(theta)+*(xyz+2)*sinf(phi)*sinf(theta);
    *(nxyz+2)=-*(xyz)*sinf(phi)+*(xyz+2)*cosf(theta);
    return (nxyz);
}