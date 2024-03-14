#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef reader_h
#define reader_h
#endif
void read(char *inpt)
{
    FILE *opt;
    FILE *ipt;
    opt=fopen("output.pdb","wr");
    ipt=fopen(inpt,"r");
    if(ipt==NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to read from input file\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    FILE *hydro;
     hydro=fopen("../src/rpl/hydrogen.cord", "r");
    if(hydro==NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to read coordinates from hydrogen.cord\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    float hx,hy,hz;     
    if(fscanf(hydro, "%*s %f %f %f", &hx,&hy,&hz) >2)
    {
        printf("----------------------------------------------------------\n");        
        printf("Setting relative coordinates of H01 to:\n x=%.3f\n z=%.3f \n z=%.3f\n",hx,hy,hz);
        printf("----------------------------------------------------------\n");
    }  
    FILE *methyl;
    methyl=fopen("../src/rpl/methyl.cord", "r");
    if(methyl==NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to read coordinates from methyl.cord\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    float C7[3],H71[3], H72[3], H73[3];
    if(fscanf(methyl, "%*s %f %f %f",&C7[0],&C7[1],&C7[2]) > 0)
    {
        printf("----------------------------------------------------------\n");        
        printf("Setting relative coordinates of C7 to:\n x=%.3f\n z=%.3f \n z=%.3f\n",C7[0],C7[1],C7[2]);
        printf("----------------------------------------------------------\n");
    }
    fclose(hydro);
    if(fscanf(methyl, "%*s %f %f %f",&H71[0],&H71[1],&H71[2]) > 0)
    {
        printf("----------------------------------------------------------\n");        
        printf("Setting relative coordinates of H71 to:\n x=%.3f\n z=%.3f \n z=%.3f\n",H71[0],H71[1],H71[2]);
        printf("----------------------------------------------------------\n");
    }
    if(fscanf(methyl, "%*s %f %f %f",&H72[0],&H72[1],&H72[2]) > 0)
    {
        printf("----------------------------------------------------------\n");        
        printf("Setting relative coordinates of H72 to:\n x=%.3f\n z=%.3f \n z=%.3f\n",H72[0],H72[1],H72[2]);
        printf("----------------------------------------------------------\n");
    }
    if(fscanf(methyl, "%*s %f %f %f",&H73[0],&H73[1],&H73[2]) > 0)
    {
        printf("----------------------------------------------------------\n");        
        printf("Setting relative coordinates of H73 to:\n x=%.3f\n z=%.3f \n z=%.3f\n",H73[0],H73[1],H73[2]);
        printf("----------------------------------------------------------\n");
    }
    fclose(methyl);
    if(opt==NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to create output file\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    else
    {   printf("----------------------------------------------------------\n");
        printf("Reading input file\n");
        FILE *fasta;
        int x=0;
        fasta=fopen("dna.fasta", "wr");
        if(fasta==NULL)
        {
            printf("----------------------------------------------------------\n");
            printf("Fatal error: Unable to create output file dna.fasta\n");
            printf("----------------------------------------------------------\n");
            exit(1);
        }
        fprintf(fasta, ">sequence\n ; generated using rna2ssdna\n");
        int j=1;
        int line_of_O=-3;
        int curr_line=1;
        while (!feof(ipt)) 
        {
            if(j==1)
            {
                fprintf(opt,"#COMMENT generated using rna2ssdna\n");
                fprintf(opt,"#COMMENT Please cite:\n");
            }
            char rtype;int rnr=1; char ident[20]=" "; char type[10]= " "; char chain[3]=" "; float coord_x; float coord_y; float coord_z; float occ; float temp_f; char elem;
            if (fscanf(ipt, "%s %*d %s %c %s %d %f %f %f %f %f %c",ident,type,&rtype,chain, &rnr,&coord_x,&coord_y,&coord_z,&occ,&temp_f,&elem) > 0){
                   if(strcmp(ident, "ATOM")==0)
                   {
                        if(rtype=='U'){
                            rtype= 'T';
                            if(strcmp(type,"H5")==0)
                            {
                                fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"C7",rtype,chain,rnr,coord_x+C7[0],coord_y+C7[1],coord_z+C7[2],1.00,0.00,'C');
                                curr_line++;
                                fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"H71",rtype,chain,rnr,coord_x+H71[0],coord_y+H71[1],coord_z+H71[2],1.00,0.00,'H');
                                curr_line++;
                                fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"H72",rtype,chain,rnr,coord_x+H72[0],coord_y+H72[1],coord_z+H72[2],1.00,0.00,'H');
                                curr_line++;
                                fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"H73",rtype,chain,rnr,coord_x+H73[0],coord_y+H73[1],coord_z+H73[2],1.00,0.00,'H');
                                curr_line++;
                            }
                            printf("Residue %d is Uracil\n", rnr );
                        }
                        if(rnr>x)
                        {
                            fprintf(fasta,"%c", rtype);
                        x++;
                        }  
                        if(strcmp(type,"O2'")==0)
                        {
                                fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"H01",rtype,chain,rnr,coord_x+hx,coord_y+hy,coord_z+hz,1.00,0.00,'H');
                                line_of_O=j;
                                curr_line++;
                        }
                        else if(j==line_of_O+1)
                        {

                        }
                        else
                        {
                            fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c",curr_line,type,rtype,chain,rnr,coord_x,coord_y,coord_z,occ,temp_f,elem);
                            curr_line++;
                            if(strcmp(type,"OP2")==0){
                                fprintf(opt, "1-\n");
                            }
                            else{
                                fprintf(opt,"\n");
                            }

                        }
                    }
                    else if(strcmp(ident, "TER")==0)  
                    {
                        printf("Found end of stem by TER (line %d)\n", j);
                        printf("----------------------------------------------------------\n");
                        fprintf(opt,"TER\n");
                        fprintf(opt,"END");
                        break;
                    }
                    else if(strcmp(ident, "END")==0)  
                    {
                        printf("Found end of structure by END (line %d)\n",j);
                        printf("----------------------------------------------------------\n");
                        fprintf(opt,"TER\n");
                        fprintf(opt,"END");
                        break;
                    }
                    else
                    {
                        if(j==0)
                        {
                            printf("Skipping line(s) at top of input file ");
                        }
                    }  
            }
            j++;
        }
        fclose(fasta);
        fclose(ipt);
    }
}