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
                            printf("INFO Beginning conversion of %s\n", argv[i+1]);
                            inpt=argv[i+1];
                            ipt=fopen(inpt, "r");
                             if(ipt==NULL)
                             {
                                printf("----------------------------------------------------------\n");
                                printf("Fatal error: Unable to read from input file\n");
                                printf("----------------------------------------------------------\n");
                                exit(1);
                            }
                            else{
                                printf("INFO Reading from input file %s\n", inpt);
                            }
                            printf("INFO Appending information to rna2ssdna.log\n");
                            printf("INFO Writing output to output.pdb\n");
                            convert(ipt);
                            
                            fclose(ipt);
                            
                            printf("INFO Done\n");
                            printf("----------------------------------------------------------\n");
                        }
                        else
                        {
                            printf("----------------------------------------------------------\n");
                            printf("ERROR Input file does not exist or is unaccesible\n");
                            printf("----------------------------------------------------------\n");
                            exit(1);
                        }
                        break;
                case 'h':
                        printf("gerneric usage: rna2ssdna -f rna-file.pdb\n");
                        printf("-f input.pdb:               declare input file for conversion\n");
                        //printf("-i insert.pdb aindx a1 a2:  insert structure at atom index aindx\n                            and rotate by axis (a1a2)\n");
                        printf("-s:                         get supportive information\n");
                        printf("-h:                         call help\n");
                        //printf("-o methyl.cord hydro.cord:  use own paramters for atom positions\n");
                        printf("----------------------------------------------------------\n");
                        break;
                case 's': 
                        if(inpt==NULL){
                            printf("ERROR No input file given");
                            printf("----------------------------------------------------------\n");
                            exit(-1);
                        }
                        
                            char *opt=argv[i+1];
                            ipt=fopen(inpt, "r");
                             if(ipt==NULL)
                             {
                                printf("ERROR Unable to read from input file\n");
                                printf("----------------------------------------------------------\n");
                                exit(1);
                            }
                            writesupporitve(ipt);
                        
                        break;
            }
        }
    }
}
