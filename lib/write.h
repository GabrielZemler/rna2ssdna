#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef write_h
#define write_h
#endif
#include "../lib/read.h"
#include "../lib/calculate.h"
void fasta(FILE *ipt, FILE *opt)
{
    FILE *log;
    log = fopen("rna2ssdna.log", "a");
    if (log == NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to append to log file rna2ssdna.log\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    fprintf(opt, ">sequence\n ; generated using rna2ssdna\n");
    int x = 0;
    fprintf(log, "Writing Sequenz in fasta code\n");
    fprintf(log, "----------------------------------------------------------\n");
    while (!feof(ipt))
    {

        int rnr = 1;
        char rtype;
        if (fscanf(ipt, "%*s %*d %*s %c %*s %d %*f %*f %*f %*f %*f %*c", &rtype, &rnr) > 1)
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
void writesupporitve(FILE * ipt)
{
rewind(ipt);
FILE * supp=fopen("supportive.log", "wr");
fasta(ipt, supp);
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
float *gethydrogen(FILE *hydro)
{
    if(hydro==NULL){
        exit(1);
    }
    FILE *log = fopen("rna2ssdna.log", "a");
    if (hydro == NULL)
    {
        fprintf(log, "Fatal error: Unable to read coordinates from hydrogen.cord\n");
        fprintf(log, "----------------------------------------------------------\n");
        exit(1);
    }
    float *hp = (float *)malloc(sizeof(float) * 3);
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
    float hp[3];

    fprintf(log, "Using standard coordinates for H2''\n");
    *(hp+0)=0.000;
    *(hp+1)=0.000;
    *(hp+2)=0.000;
    fprintf(log, "H2'': x=%.3f y=%.3f z=%.3f\n",*(hp+0), *(hp+1), *(hp+2));
    fprintf(log, "----------------------------------------------------------\n");
    /*
    else{
        printf("called h\n");
        hp = gethydrogen(hydro);
    }*/
    float C7[3], H71[3], H72[3], H73[3];
         fprintf(log, "Using standard parameters for C7, H71, H72 and H73\n");
        *(C7+0)=1.500;
        *(C7+1)=-0.002;
        *(C7+2)=-0.016;
        fprintf(log, "C7: r=%.3f theta=%.3f phi=%.3f\n",*(C7+0), *(C7+1), *(C7+2));
        *(H71+0)=2.128;
        *(H71+1)=-0.113;
        *(H71+2)=0.283;
        fprintf(log, "H71: r=%.3f theta=%.3f phi=%.3f\n",*(H71+0), *(H71+1), *(H71+2));
        *(H72+0)=2.128;
        *(H72+1)=-0.042;
        *(H72+2)=-0.519;
        fprintf(log, "H72: r=%.3f theta=%.3f phi=%.3f\n",*(H72+0), *(H72+1), *(H72+2));
        *(H73+0)=2.129 ;
        *(H73+1)=0.448 ;
        *(H73+2)=0.220;
        fprintf(log, "H73: r=%.3f theta=%.3f phi=%.3f\n",*(H73+0), *(H73+1), *(H73+2));
        fprintf(log, "----------------------------------------------------------\n");
        /*else{
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
        } */
    
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
                    fprintf(log,"Position of C5 and 02 in residue %d:\n %f %f %f\n %f %f %f\n", rnr, *(positC5 + 0), *(positC5 + 1), *(positC5 + 2), *(positO2 + 0), *(positO2 + 1), *(positO2 + 2));
                    float *vec = getvec(positO2, positC5);
                    fprintf(log,"Vector O2C5 in residue %d: %.3f %.3f %.3f\n",rnr, *(vec + 0), *(vec + 1), *(vec + 2));
                    // get theta and phi of c5o2
                    float c5o2_theta = theta(vec);
                    float c5o2_phi = phi(vec);
                    if(c5o2_theta>3.15){
                        c5o2_theta=c5o2_theta-3.14;
                    }
                    fprintf(log,"Values of theta and phi of vector O2C5 in residue %d:%.3f %.3f\n",rnr,c5o2_theta,c5o2_phi);
                    //rotate and give back in cartesian
                    float *nc7 = revert(C7[0], c5o2_theta+C7[1], c5o2_phi+C7[2]);
                    float *nh1 = revert(H71[0], c5o2_theta+H71[1], c5o2_phi+H71[2]);
                    float *nh2 = revert(H72[0], c5o2_theta+H72[1], c5o2_phi+H72[2]);
                    float *nh3 = revert(H73[0], c5o2_theta+H73[1], c5o2_phi+H73[2]);
                    fseek(ipt, len, SEEK_SET);
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, type, rtype, chain, rnr, coord_x, coord_y, coord_z, occ, temp_f, elem);
                    curr_line++;
                    fprintf(opt, "ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n", curr_line, "C7", rtype, chain, rnr, coord_x + *(nc7 + 0), coord_y + *(nc7 + 1), coord_z + *(nc7 + 2), 1.00, 0.00, 'C');
                    curr_line++;
                    fprintf(log, "Placing C7 in resudue %d at coordinates %.3f %.3f %.3f\n", rnr, coord_x+*(nc7 + 0), coord_y + (*nc7 + 1), coord_z + *(nc7 + 2));
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
    fclose(log);
    fclose(opt);
}
