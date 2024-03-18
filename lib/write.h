#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef write_h
#define write_h
#endif
void fasta(FILE *ipt, char *oname){
    FILE *opt;
    FILE *log;
    opt=fopen(oname, "wr");
    log=fopen("rna2ssdna.log", "a");
    if(log==NULL)
        {
            printf("----------------------------------------------------------\n");
            printf("Fatal error: Unable to append to log file rna2ssdna.log\n");
            printf("----------------------------------------------------------\n");
            exit(1);
        }
    if(opt==NULL)
        {
            fprintf(log, "Fatal error: Unable to write to %s\nExiting\n", oname);
            fprintf(log,"----------------------------------------------------------\n");
            exit(1);
        }
    fprintf(opt, ">sequence\n ; generated using rna2ssdna\n");
    int x=0;
    fprintf(log, "Writing Sequenz in fasta code to %s\n", oname);
    fprintf(log,"----------------------------------------------------------\n");
    while (!feof(ipt)) 
        {
           
            int rnr=1; char rtype;
            if (fscanf(ipt, "%*s %*d %*s %c %*s %d %*f %*f %*f %*f %*f %*c",&rtype, &rnr) > 0)
                {
                    if(rnr>x)
                        {
                        if(rtype=='U'){
                            rtype='T';
                        }
                        fprintf(opt,"%c", rtype);
                        x++;
                        }   
                } 
        }  
    fclose(log);
    fclose(opt);
}
void writelog(){
    FILE *log;
    log=fopen("rna2ssdna.log", "w");
    fprintf(log,"----------------------------------------------------------\n");
    fprintf(log,"WELCOME TO RNA2SSDNA\nPlease cite:\n");
    fprintf(log,"----------------------------------------------------------\n");
    fclose(log);
}
float *gethydrogen(){
    FILE *hydro=fopen("../src/rpl/hydrogen.cord", "r");
    FILE *log=fopen("rna2ssdna.log", "a");
    if(hydro==NULL)
    {
        fprintf(log,"Fatal error: Unable to read coordinates from hydrogen.cord\n");
        fprintf(log,"----------------------------------------------------------\n");
        exit(1);
    }
    float *hp=malloc(sizeof(float) * 3);
        if(fscanf(hydro, "%*s %f %f %f", (hp+0),(hp+1),(hp+2)) >2)
        {       
        fprintf(log,"Setting relative coordinates of H2'' to:\n x=%.3f\n z=%.3f \n z=%.3f\n",hp[0],hp[1],hp[2]);
        fprintf(log,"----------------------------------------------------------\n");
        }
        fclose(hydro);
        fclose(log);
    return ( hp );
}
void convert(FILE *ipt){
    FILE *opt=fopen("output.pdb", "w");
    rewind(ipt);
    int j=1;
    int line_of_O=-3;
    int curr_line=1;
    float *hp=gethydrogen();
    FILE *log=fopen("rna2ssdna.log", "a");
    while (!feof(ipt))
    {
        char rtype;int rnr=1; char ident[20]=" "; char type[10]= " "; char chain[3]=" "; float coord_x; float coord_y; float coord_z; float occ; float temp_f; char elem;
        if (fscanf(ipt, "%s %*d %s %c %s %d %f %f %f %f %f %c",ident,type,&rtype,chain, &rnr,&coord_x,&coord_y,&coord_z,&occ,&temp_f,&elem) > 0)
        {
            if(strcmp(ident,"ATOM")==0)
            {
                if(strcmp(type,"O2'")==0)
                {
                    fprintf(log,"Repacing O2' in residue %d with H2'' and correcting atom positons\n", rnr);
                    fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c\n",curr_line,"H2''",rtype,chain,rnr,coord_x+*(hp+0),coord_y+*(hp+1),coord_z+*(hp+2),1.00,0.00,'H');
                    line_of_O=j;
                    curr_line++;
                }
                else if(strcmp(type,"HO2'")==0)
                {
                    fprintf(log,"Deleting atom HO2' in residue %d\n", rnr);
                }
                else{
                    fprintf(opt,"ATOM%7d%5s  D%c%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12c",curr_line,type,rtype,chain,rnr,coord_x,coord_y,coord_z,occ,temp_f,elem);
                    curr_line++;
                    if(strcmp(type,"OP2")==0)
                    {
                        fprintf(opt, "1-\n");
                    }
                    else
                    {
                        fprintf(opt,"\n");
                    } 
                }
            }
            else if (strcmp(ident,"TER")==0)
            {
                fprintf(log, "Reached end of file by TER (line %d)\n", j);
                fprintf(log,"----------------------------------------------------------\n");
                fprintf(opt,"TER\n");
            }
            else if(strcmp(ident, "END")==0)  
            {
                fprintf(log,"Found end of structure by END (line %d)\n",j);
                fprintf(log,"----------------------------------------------------------\n");
                fprintf(opt,"TER\n");
                break;
            }
        }
        j++;
    }
    free(hp);
    fclose(log);
    fclose(opt);   
}