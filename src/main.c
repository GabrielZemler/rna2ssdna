/* This is the main file of the programm which reads the user input and calls the needed functions in ./lib
*/
#include<stdio.h>
#include <stdlib.h>
#include "../lib/write.h"
int main(int argc, char **argv)
{
    writelog();
    printf("----------------------------------------------------------\n");
    printf("                WELLCOME TO RNA2ssDNA \n");
    printf("----------------------------------------------------------\n");
    int i;
    FILE *ipt;
    char *inpt; 
    for (i = 1; i < argc; i++) 
    {
        if (argv[i][0] == '-') 
        {
            switch (argv[i][1]) 
            {
 
                case 'f':
                        if(argv[i+1]!=NULL)
                        {
                            printf("Beginning conversion of %s\n", argv[i+1]);
                            inpt=argv[i+1];
                            ipt=fopen(inpt, "r");
                             if(ipt==NULL)
                             {
                                printf("----------------------------------------------------------\n");
                                printf("Fatal error: Unable to read from input file\n");
                                printf("----------------------------------------------------------\n");
                                exit(1);
                            }
                            printf("Appending information to rna2ssdna.log\n");
                            printf("Writing output to output.pdb\n");
                            convert(ipt);
                            //test();
                            fclose(ipt);
                            printf("Done\n");
                            printf("----------------------------------------------------------\n");
                        }
                        else
                        {
                            printf("----------------------------------------------------------\n");
                            printf("Fatal error: Input file does not exist or is unaccesible\n");
                            printf("----------------------------------------------------------\n");
                            exit(1);
                        }
                        break;
                case 'h':
                        printf("gerneric usage: rna2ssdna -f rna-file.pdb\n");
                        printf("-f input.pdb:   declare input file for conversion\n");
                        printf("-o output.pdb:  declare output file (default: dna.pdb)\n");
                        printf("-s dna.fasta:   get supportive information: output dna in fasta code\n");
                        printf("-h:             call help\n ");
                        break;
                case 's': 
                        if(inpt==NULL){
                            printf("Fatal error: No input file given");
                            printf("----------------------------------------------------------\n");
                            exit(-1);
                        }
                        if(argv[i+1]!=NULL)
                        {
                            char *opt=argv[i+1];
                            ipt=fopen(inpt, "r");
                             if(ipt==NULL)
                             {
                                printf("Fatal error: Unable to read from input file\n");
                                printf("----------------------------------------------------------\n");
                                exit(1);
                            }
                            fasta(ipt,opt);
                        }
                        break;
            }
        }
    }
}
