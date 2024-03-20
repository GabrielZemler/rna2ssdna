#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef write_h
#define write_h
#endif
#include "../lib/read.h"
#include "../lib/calculate.h"
void fasta(FILE *ipt, char *oname)
{
    FILE *opt;
    FILE *log;
    opt = fopen(oname, "wr");
    log = fopen("rna2ssdna.log", "a");
    if (log == NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to append to log file rna2ssdna.log\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    if (opt == NULL)
    {
        fprintf(log, "Fatal error: Unable to write to %s\nExiting\n", oname);
        fprintf(log, "----------------------------------------------------------\n");
        exit(1);
    }
    fprintf(opt, ">sequence\n ; generated using rna2ssdna\n");
    int x = 0;
    fprintf(log, "Writing Sequenz in fasta code to %s\n", oname);
    fprintf(log, "----------------------------------------------------------\n");
    while (!feof(ipt))
    {

        int rnr = 1;
        char rtype;
        if (fscanf(ipt, "%*s %*d %*s %c %*s %d %*f %*f %*f %*f %*f %*c", &rtype, &rnr) > 0)
        {
            if (rnr > x)
            {
                if (rtype == 'U')
                {
                    rtype = 'T';
                }
                fprintf(opt, "%c", rtype);
                x++;
            }
        }
    }
    fclose(log);
    fclose(opt);
}
void writelog()
{
    FILE *log;
    log = fopen("rna2ssdna.log", "w");
    fprintf(log, "----------------------------------------------------------\n");
    fprintf(log, "WELCOME TO RNA2SSDNA\nPlease cite:\n");
    fprintf(log, "----------------------------------------------------------\n");
    fclose(log);
}
float *gethydrogen()
{
    FILE *hydro = fopen("../src/rpl/hydrogen.cord", "r");
    FILE *log = fopen("rna2ssdna.log", "a");
    if (hydro == NULL)
    {
        fprintf(log, "Fatal error: Unable to read coordinates from hydrogen.cord\n");
        fprintf(log, "----------------------------------------------------------\n");
        exit(1);
    }
    float *hp = malloc(sizeof(float) * 3);
    if (fscanf(hydro, "%*s %f %f %f", (hp + 0), (hp + 1), (hp + 2)) > 2)
    {
        fprintf(log, "Setting relative coordinates of H2'' to:\n x=%.3f\n z=%.3f \n z=%.3f\n", hp[0], hp[1], hp[2]);
        fprintf(log, "----------------------------------------------------------\n");
    }
    fclose(hydro);
    fclose(log);
    return (hp);
}
void convert(FILE *ipt)
{
    int l = length(ipt);
    int *posC5 = positionC5(l, ipt);
    int *posO2 = positionO2(l, ipt);
    FILE *log = fopen("rna2ssdna.log", "a");
    FILE *opt = fopen("output.pdb", "w");
    if (opt == NULL)
    {
        fprintf(log, "Fatal error: Unable to write to output.pdb\n");
        fprintf(log, "----------------------------------------------------------\n");
        exit(1);
    }
    fprintf(opt, "#COMMENT gernerated with rna2ssdna\n");
    rewind(ipt);
    int j = 1;
    int line_of_O = -3;
    int curr_line = 1;
    float *hp = gethydrogen();
    FILE *methyl;
    methyl = fopen("../src/rpl/methyl.cord", "r");
    if (methyl == NULL)
    {
        fprintf(log, "Fatal error: Unable to read coordinates from methyl.cord\n");
        fprintf(log, "----------------------------------------------------------\n");
        exit(1);
    }
    float C7[3], H71[3], H72[3], H73[3];
    if (fscanf(methyl, "%*s %f %f %f", &C7[0], &C7[1], &C7[2]) > 0)
    {
        fprintf(log, "Reading parameters of C7:\n r=%.3f\n theta=%.3f \n phi=%.3f\n", C7[0], C7[1], C7[2]);
        fprintf(log, "----------------------------------------------------------\n");
    }
    if (fscanf(methyl, "%*s %f %f %f", &H71[0], &H71[1], &H71[2]) > 0)
    {
        fprintf(log, "Reading parameters of H71:\n r=%.3f\n theta=%.3f \n phi=%.3f\n", H71[0], H71[1], H71[2]);
        fprintf(log, "----------------------------------------------------------\n");
    }
    if (fscanf(methyl, "%*s %f %f %f", &H72[0], &H72[1], &H72[2]) > 0)
    {
        fprintf(log, "Reading parameters of H72:\n r=%.3f\n theta=%.3f \n phi=%.3f\n", H72[0], H72[1], H72[2]);
        fprintf(log, "----------------------------------------------------------\n");
    }
    if (fscanf(methyl, "%*s %f %f %f", &H73[0], &H73[1], &H73[2]) > 0)
    {
        fprintf(log, "Reading parameters of H73:\n r=%.3f\n theta=%.3f \n phi=%.3f\n", H73[0], H73[1], H73[2]);
        fprintf(log, "----------------------------------------------------------\n");
    }
    while (!feof(ipt))
    {
        int len;
        int pos = 0;
        char rtype;
        int rnr = 1;
        char ident[20] = " ";
        char type[10] = " ";
        char chain[3] = " ";
        float coord_x;
        float coord_y;
        float coord_z;
        float occ;
        float temp_f;
        char elem;
        if (fscanf(ipt, "%s %*d %s %c %s %d %f %f %f %f %f %c", ident, type, &rtype, chain, &rnr, &coord_x, &coord_y, &coord_z, &occ, &temp_f, &elem) > 0)
        {
            if (strcmp(ident, "ATOM") == 0)
            {
                if (rtype == 'U')
                {
                    rtype = 'T';
                }
                if (rtype == 'T' && (strcmp(type, "H5") == 0))
                {
                    fprintf(log, "Deleting H5 in residue %d\n", rnr);
                }
                else if (rtype == 'T' && (strcmp(type, "C5") == 0))
                {
                    len = ftell(ipt);
                    pos = rnr - 1;
                    // printf("Atom pos of C5: %d \n O2: %d\n",*(posC5+pos), *(posO2+pos));
                    float *positC5 = getpos(ipt, *(posC5 + pos));
                    float *positO2 = getpos(ipt, *(posO2 + pos));
                    printf("pos of C5 and 02 in res %d:\n %f %f %f\n %f %f %f\n", rnr, *(positC5 + 0), *(positC5 + 1), *(positC5 + 2), *(positO2 + 0), *(positO2 + 1), *(positO2 + 2));
                    float *vec = getvec(positO2, positC5);
                    printf("vector o2c5: %.3f %.3f %.3f\n", *(vec + 0), *(vec + 1), *(vec + 2));
                    // get theta and phi of c5o2
                    float c5o2_theta = theta(vec);
                    float c5o2_phi = phi(vec);
                    if(c5o2_theta>3.15){
                        c5o2_theta=c5o2_theta-3.14;
                    }
                    printf("Theta  and phi of o2c5:%.3f %.3f\n",c5o2_theta,c5o2_phi);
                    //rotate and give back in cartesian
                    float *nc7 = revert(C7[0], c5o2_theta+C7[1], c5o2_phi+C7[2]);
                    float *nh1 = revert(H71[0], c5o2_theta+H71[1], c5o2_phi+H71[2]);
                    float *nh2 = revert(H72[0], c5o2_theta+H72[1], c5o2_phi+H72[2]);
                    float *nh3 = revert(H73[0], c5o2_theta+H73[1], c5o2_phi+H73[2]);
                    printf("releative coordinates of h71: %.3f %.3f %.3f\n", *(nh1+0), *(nh1+1), *(nh1+2));
                    printf("coordinates coorx: %.3f coory %.3f coorz %.3f\n", coord_x, coord_y, coord_z);
                    fseek(ipt, len, SEEK_SET);
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, type, rtype, chain, rnr, coord_x, coord_y, coord_z, occ, temp_f, elem);
                    curr_line++;
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "C7", rtype, chain, rnr, coord_x + *(nc7 + 0), coord_y + *(nc7 + 1), coord_z + *(nc7 + 2), 1.00, 0.00, 'C');
                    curr_line++;
                    fprintf(log, "Placing C7 in resudue %d at coordinates %.3f %.3f %.3f\n", rnr, coord_x+*(nc7 + 0), coord_y + (*nc7 + 1), coord_z + *(nc7 + 2));
                    printf("coordinates coorx: %.3f coory %.3f coorz %.3f\n", coord_x, coord_y, coord_z);
                    printf("calculation h71y=%.3f\n", *(nh1+1)+coord_y);
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "H71", rtype, chain, rnr, coord_x + *(nh1 + 0), coord_y + *(nh1 + 1), coord_z + *(nh1 + 2), 1.00, 0.00, 'H');
                    curr_line++;
                    fprintf(log, "Placing H71 in resudue %d at coordinates %.3f %.3f %.3f\n", rnr, coord_x+*(nh1 + 0), coord_y + *(nh1 + 1), coord_z + *(nh1 + 2));
                    
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "H72", rtype, chain, rnr, coord_x + *(nh2 + 0), coord_y + *(nh2 + 1), coord_z + *(nh1 + 2), 1.00, 0.00, 'H');
                    curr_line++;
                    fprintf(log, "Placing H72 in resudue %d at coordinates %.3f %.3f %.3f\n", rnr, coord_x+*(nh2 + 0), coord_y + (*nh2 + 1), coord_z + *(nh2 + 2));
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "H73", rtype, chain, rnr, coord_x + *(nh3 + 0), coord_y + *(nh3 + 1), coord_z + *(nh3 + 2), 1.00, 0.00, 'H');
                    curr_line++;
                    fprintf(log, "Placing H73 in resudue %d at coordinates %.3f %.3f %.3f\n", rnr, coord_x+*(nh3 + 0), coord_y + (*nh3 + 1), coord_z + *(nh3 + 2));
                    free(nh1);
                    free(nh2);
                    free(nh3);
                    free(vec);
                    free(nc7);
                    free(positC5);
                    free(positO2);
                }
                else if (strcmp(type, "O2'") == 0)
                {
                    fprintf(log, "Repacing O2' in residue %d with H2'' and correcting atom positons\n", rnr);
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "H2''", rtype, chain, rnr, coord_x + *(hp + 0), coord_y + *(hp + 1), coord_z + *(hp + 2), 1.00, 0.00, 'H');
                    line_of_O = j;
                    curr_line++;
                }
                else if (strcmp(type, "HO2'") == 0)
                {
                    fprintf(log, "Deleting atom HO2' in residue %d\n", rnr);
                }
                else
                {
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c", curr_line, type, rtype, chain, rnr, coord_x, coord_y, coord_z, occ, temp_f, elem);
                    curr_line++;
                    if (strcmp(type, "OP2") == 0)
                    {
                        fprintf(opt, "1-\n");
                    }
                    else
                    {
                        fprintf(opt, "\n");
                    }
                }
            }
            else if (strcmp(ident, "TER") == 0)
            {
                fprintf(log, "Reached end of file by TER (line %d)\n", j);
                fprintf(log, "----------------------------------------------------------\n");
                fprintf(opt, "TER\n");
            }
            else if (strcmp(ident, "END") == 0)
            {
                fprintf(log, "Found end of structure by END (line %d)\n", j);
                fprintf(log, "----------------------------------------------------------\n");
                fprintf(opt, "TER\n");
                break;
            }
        }
        j++;
    }
    free(hp);
    fclose(log);
    fclose(opt);
}
void test()
{
    float v1[3] = {1.00, 0.00, 0.00};
    float v2[3] = {0.00, 1.00, 0.00};
    float v3[3] = {0.00, 0.00, 1.00};
    float vecC5[3] = {-3.645, 2.504, 23.333};
    float vec02[3] = {-3.230, -1.489, 23.423};
    float *vecO2C5 = getvec(vec02, vecC5);
    float c7[3] = {-0.131, 1.494, -0.031};
    float h1[3] = {-0.743, 1.823, 0.809};
    float h2[3] = {0.858, 1.947, 0.040};
    float h3[3] = {-0.603, 1.799, -0.965};
    float thta = theta(vecO2C5);
    float ph = phi(vecO2C5);
    float len_of_o2c5 = r(vecO2C5);
    printf("Vector c5O2: %.3f %.3f %.3f\n", *(vecO2C5 + 0), *(vecO2C5 + 1), *(vecO2C5 + 2));
    printf("theta: %.3f\n", thta);
    printf("phi: %.3f\n", ph);
    /*printf("theta: %.3f\n", thta);
    printf("phi: %.3f\n", ph);
    float* nh1=rotate(h1,-thta,-ph);
    float* nh2=rotate(h2,-thta,-ph);
    float* nh3=rotate(h3,-thta,-ph);
    float* nc7=rotate(c7,-thta,-ph);
    float* no2c5=rotate(vecO2C5,0,-ph);
    no2c5=rotate(vecO2C5,-thta,0);
     printf("new coord of vec o2c5: %.3f %.3f %.3f\n", *(no2c5+0), *(no2c5+1),*(no2c5+2));
    printf("new pos of h1: %.3f %.3f %.3f\n", *(nh1+0), *(nh1+1),*(nh1+2));
    printf("new pos of h2: %.3f %.3f %.3f\n", *(nh2+0), *(nh2+1),*(nh2+2));
    printf("new pos of h3: %.3f %.3f %.3f\n", *(nh3+0), *(nh3+1),*(nh3+2));
    printf("new pos of c7: %.3f %.3f %.3f\n", *(nc7+0), *(nc7+1),*(nc7+2));
    free(nh1);
    free(nh2);
    free(nh3);
    free(nc7);
    free(vecO2C5);*/
    float c7_phi = phi(c7);
    float c7_theta = theta(c7);
    float c7_r = r(c7);
    printf("%.3f %.3f %.3f\n", c7_r, c7_theta, c7_phi);
    float h1_phi = phi(h1);
    float h1_theta = theta(h1);
    float h1_r = r(h1);
    printf("%.3f %.3f %.3f\n", h1_r, h1_theta, h1_phi);
    float h2_phi = phi(h2);
    float h2_theta = theta(h2);
    float h2_r = r(h2);
    printf("%.3f %.3f %.3f\n", h2_r, h2_theta, h2_phi);
    float h3_phi = phi(h3);
    float h3_theta = theta(h3);
    float h3_r = r(h3);
    printf("%.3f %.3f %.3f\n", h3_r, h3_theta, h3_phi);
    printf("Rotation by theta=0 phi=3.14\n");
    float *h1_xyz = revert(h1_r, h1_theta - thta, h1_phi);
    printf("new coords of h71: %.3f %.3f %.3f\n", *(h1_xyz + 0), *(h1_xyz + 1), *(h1_xyz + 2));
    float *h2_xyz = revert(h2_r, h2_theta - thta, h2_phi);
    printf("new coords of h72: %.3f %.3f %.3f\n", *(h2_xyz + 0), *(h2_xyz + 1), *(h2_xyz + 2));
    free(h2_xyz);
    float *h3_xyz = revert(h3_r, h3_theta - thta, h3_phi);
    printf("new coords of h73: %.3f %.3f %.3f\n", *(h3_xyz + 0), *(h3_xyz + 1), *(h3_xyz + 2));
    free(h3_xyz);
    float *c7_xyz = revert(c7_r, c7_theta - thta, c7_phi);
    printf("new coords of c7: %.3f %.3f %.3f\n", *(c7_xyz + 0), *(c7_xyz + 1), *(c7_xyz + 2));
    free(c7_xyz);
    h1_xyz = revert(len_of_o2c5, 0, ph);
    printf("renew coords of h71: %.3f %.3f %.3f\n", *(h1_xyz + 0), *(h1_xyz + 1), *(h1_xyz + 2));
    free(h1_xyz);
    printf("%.3f %.3f %.3f\n", c7_r, c7_theta-thta, c7_phi-ph);
    printf("%.3f %.3f %.3f\n", h1_r, h1_theta-thta, h1_phi-ph);
    printf("%.3f %.3f %.3f\n", h2_r, h2_theta - thta, h2_phi-ph);
    printf("%.3f %.3f %.3f\n", h3_r, h3_theta-thta, h3_phi-ph);
    /*
    float hp1=yaw(vecO2C5);
    float hp2=roll(vecO2C5);
    float hp3=pitch(vecO2C5);
    for (int i = 0; i < 3; i++)
    {
        printf("x%d:%f\n",i+1,*(vecO2C5+i));
    }
    printf("alpha:%f\n",hp1);
    printf("beta:%f\n",hp2);
    printf("gamma:%f\n",hp3);
    float *nc7=rotate(c7, -hp1,-hp2,-hp3);
    float *nh1=rotate(h1, -hp1,-hp2,-hp3);
    float *nh2=rotate(h2, -hp1,-hp2,-hp3);
    float *nh3=rotate(h3, -hp1,-hp2,-hp3);
    printf("h71:%.3f %.3f %.3f\n",*(nh1+0) ,*(nh1+1),*(nh1+2));
    printf("h72:%.3f %.3f %.3f\n", *(nh2+0),*(nh2+1),*(nh2+2));
    printf("h73:%.3f %.3f %.3f\n", *(nh3+0),*(nh3+1),*(nh3+2));
    printf("c7:%.3f %.3f %.3f\n", *(nc7+0),*(nc7+1),*(nc7+2));*/
}