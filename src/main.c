/* This is the main file of the programm which reads the user input and calls the needed functions in ./lib
*/
#include<stdio.h>
#include <stdlib.h>
#include "../lib/reader.h"
int main(int argc, char **argv)
{
    printf("----------------------------------------------------------\n");
    printf("                WELLCOME TO RNA2ssDNA \n");
    printf("----------------------------------------------------------\n");
    int i;
    for (i = 1; i < argc; i++) 
    {
        if (argv[i][0] == '-') 
        {
            switch (argv[i][1]) 
            {
 
                case 'f':
                        if(argv[i+1]!=NULL)
                        {
                            printf("Trying to read input from file %s\n", argv[i+1]);
                            char *inpt=argv[i+1];
                            read(inpt);
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
                        printf("gerneric usage: rna2ssdna -f rna-file.pdb\n -o output.pdb:  declare output file (default: dna.pdb)\n -h:             call help\n");
                        break;
            }
        }
    }
}
