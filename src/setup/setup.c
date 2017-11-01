
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <stdio.h>


//gcc -o setup setup.c -L. -lgmp -lpbc
//./setup <../../data/param/a.param 


int main(int argc, char **argv) {

//clain element
pairing_t pairing; 
element_t g,g1,g2;
element_t A;
element_t v;
element_t z,y;


//int n,m;
//scanf("%d %d", &n,&m);
int n=10,m=10;

element_t T[n];
element_t Z[m];
element_t V[m];

pbc_demo_pairing_init(pairing, argc, argv);

element_init_G1(g,pairing);
element_init_G1(g1,pairing);
element_init_G1(g2,pairing);
element_init_G1(v,pairing);
element_init_Zr(y,pairing);
element_init_Zr(z,pairing);
element_init_GT(A,pairing);

element_random(g);
 element_printf("g = %B\n", g);

element_random(g2);
element_printf("g2 = %B\n", g2);

element_random(y);
element_printf("y = %B\n", y);


element_random(z);
element_printf("z = %B\n", z);


for (int i=0;i<n;i++){
	element_init_G1(T[i],pairing);
	element_random(T[i]);
	element_printf("T[%d] = %B\n",i, T[i]);
}

for (int i=0;i<m;i++){
	element_init_Zr(Z[i],pairing);
	element_random(Z[i]);
}

for (int i=0;i<m;i++){
	element_init_G1(V[i],pairing);
	element_pow_zn(V[i],g,Z[i]);
	element_printf("V[%d] = %B\n",i, V[i]);
}

element_pow_zn(g1,g,y);
element_printf("g1 = %B\n", g1);

element_pow_zn(v,g,z);
element_printf("v = %B\n", v);

element_pairing(A,g1,g2);
element_printf("A = %B\n", A);

//write file PP : g,g1,g2,v,v1...vm,t1...tn,A
FILE *fc;
fc=fopen("../..//data/setup_data/PP","w+");
element_fprintf(fc,"{\"g\":\"%B\"",g);
element_fprintf(fc,",\r\"g1\":\"%B\"",g1);
element_fprintf(fc,",\r\"g2\":\"%B\"",g2);
element_fprintf(fc,",\r\"z\":\"%B\"",z);
element_fprintf(fc,",\r\"v\":\"%B\"",v);
element_fprintf(fc,",\r\"A\":\"%B\"",A);
for (int i=0;i<m;i++){
	
	element_fprintf(fc,",\r\"V-%d\":\"%B\"",i,V[i]);
}
for (int i=0;i<n;i++){

	element_fprintf(fc,",\r\"T-%d\":\"%B\"",i,T[i]);
	
}
element_fprintf(fc,"}");
fclose(fc);

//write file MK :y
fc=fopen("../../data/setup_data/MK","w+");
element_fprintf(fc,"{\"y\":\"%B\"",y);
fclose(fc);


//clear element
element_clear(g);
element_clear(g1);
element_clear(g2);
element_clear(y);
element_clear(z);
element_clear(v);
pairing_clear(pairing);

return 0;


}
