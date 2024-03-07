#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef reader_h
#define reader_h
#endif
void read(char *inpt)
{
    FILE *ipt;
    char * line = NULL;
    size_t len = 0;
    size_t read;
    ipt=fopen(inpt,"r");
    if(ipt==NULL)
    {
        printf("----------------------------------------------------------\n");
        printf("Fatal error: Unable to read from input file\n");
        printf("----------------------------------------------------------\n");
        exit(1);
    }
    else
    {
        printf("Reading \n");
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
        while (!feof(ipt)) 
        {
            char d;int n=1; char a[100]="ATOM";
            char *ip;
            ip=a;
            if (fscanf(ipt, "%*s %*d %*s %c %*c %d %*f %*f %*f %*f %*f %*c",&d, &n) == 2){
                    if(n>x)
                    {
                        fprintf(fasta,"%c", d );
                        x++;
                    }           
                }
                else if(fscanf(ipt,"%s%*[^\n]", a))
                {
                    printf("%s",a);
                    if(strcmp(a,"TER")==0)
                    {
                        printf("ITS TER\n");
                       
                    }
                }
                else
                {
                    fscanf(ipt, "%*[^\n]\n");
                    printf("called");
                }
        }
        fclose(fasta);
        fclose(ipt);
    }
}
void flag(){

}