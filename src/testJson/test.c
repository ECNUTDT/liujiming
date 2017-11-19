#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 


int main(int argc, char **argv) {

//clain element
pairing_t pairing; 
element_t g,g1;
pbc_demo_pairing_init(pairing, argc, argv);

element_init_G1(g,pairing);
element_init_G1(g1,pairing);
element_random(g);

element_printf("g = %B\n", g);

int size = pairing_length_in_bytes_compressed_G1(pairing);
element_printf("size = %d\n", size);
int size2 = element_length_in_bytes(g);
element_printf("size2 = %d\n", size2);

//from elements
unsigned char *str = malloc(size2);
int readB=element_to_bytes(str,g);
int writeB=element_from_bytes(g1,str);


element_printf("g1 = %B\n", g1);

//write file
FILE *fc;
fc=fopen("test.data","w+");
fwrite( str, sizeof( unsigned char ), size2, fc );
fclose(fc);

//read file
fc = fopen("test.data", "rb");
unsigned char buf[size2];
fread(buf,sizeof(unsigned char), size2,fc);
fclose(fc);

//to element_t
element_from_bytes(g1,buf);
element_printf("g1 = %B\n", g1);





return 0;

}
