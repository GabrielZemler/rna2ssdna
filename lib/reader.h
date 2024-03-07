#include <stdio.h>
#include <stdlib.h>

#ifndef reader_h
#define reader_h
#endif
void read(char *inpt){
FILE *ipt;
char * line = NULL;
size_t len = 0;
size_t read;
ipt=fopen(inpt,"r");
if(ipt==NULL){
    printf("----------------------------------------------------------\n");
    printf("Fatal error: Unable to read from input file\n");
    printf("----------------------------------------------------------\n");
}
else{
    printf("Reading \n");
    /*while ((read = getline(&line, &len, ipt)) != -1) {
        printf("Retrieved line of length %zu:\n", read);
        printf("%s", line);
    }*/
    while (!feof(ipt)) {
    char c[4]; char d;int n=0; int x=0; int p=0;
    if (fscanf(ipt, "%*s %*d %*s %c %*c %*d %*f %*f %*f %*f %*f %*c",  &d) < 1)
      break;
    printf("%c\n", d);
 }
}
}
void flag(){

}