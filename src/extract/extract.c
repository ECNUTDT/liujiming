
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <stdio.h>
#include <math.h>

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


element_t t[n];
element_t Z[m];
element_t V[m];
element_t d[n+1];
element_t D[n+1];

element_t w[n+1];
element_t r[n+1];
element_t T[n+1];

element_t element_long_i;
element_t element_long_n;

element_t temp_neg_ri;
element_t temp_pow_delta;
element_t temp_delta;
element_t temp_delta_pai;
element_t temp_xn;
element_t temp_sub1;
element_t temp_sub2;
element_t temp_div;

element_t polynomial_q;


pbc_demo_pairing_init(pairing, argc, argv);

element_init_Zr(element_long_i,pairing);
element_init_Zr(element_long_n,pairing);


element_init_G1(g,pairing);
element_init_G1(g1,pairing);
element_init_G1(g2,pairing);
element_init_G1(v,pairing);
element_init_Zr(y,pairing);
element_init_Zr(z,pairing);
element_init_GT(A,pairing);

element_random(g);
 //element_printf("g = %B\n", g);

element_random(g2);
//element_printf("g2 = %B\n", g2);

element_random(y);
//element_printf("y = %B\n", y);


element_random(z);
//element_printf("z = %B\n", z);


for (int i=0;i<n;i++){
	element_init_G1(t[i],pairing);
	element_random(t[i]);
	//element_printf("T[%d] = %B\n",i, T[i]);
}

for (int i=0;i<m;i++){
	element_init_Zr(Z[i],pairing);
	//element_random(Z[i]);
}

for (int i=0;i<m;i++){
	element_init_G1(V[i],pairing);
	element_pow_zn(V[i],g,Z[i]);
	//element_printf("V[%d] = %B\n",i, V[i]);
}

/*
element_pow_zn(g1,g,y);
element_printf("g1 = %B\n", g1);

element_pow_zn(v,g,z);
element_printf("v = %B\n", v);

element_pairing(A,g1,g2);
element_printf("A = %B\n", A);
*/
element_printf("===============start===============");


for(int i=0;i<n+1;i++){
	
	
	element_init_Zr(w[i],pairing);
	element_random(w[i]);

	element_init_Zr(r[i],pairing);
	element_random(r[i]);
}

for(int i=0;i<n+1;i++){

	

	//di=element_pow_zn(di,g,-ri) 
	element_init_G1(d[i],pairing);
	element_init_Zr(temp_neg_ri,pairing);
	element_neg(temp_neg_ri,r[i]);
	
	element_pow_zn(d[i],g,temp_neg_ri);
	element_printf("d-%d = %B\n",i, d[i]);

	//delta
	element_init_Zr(temp_delta,pairing);
	element_init_G1(temp_pow_delta,pairing);
	element_init_G1(temp_delta_pai,pairing);
	element_init_Zr(temp_sub1,pairing);
	element_init_Zr(temp_sub2,pairing);
	element_init_Zr(temp_div,pairing);

	element_set1(temp_delta_pai);
	for(int j=0;j<n+1;j++){
		element_set1(temp_delta);
		//element_printf("temp_delta-%d = %B\n",j, temp_delta);
		for(int k=0;k<n+1;k++){
			if(k==j) continue;
			else if(k==i) continue;
			else {
				element_sub(temp_sub1,w[i],w[k]);
				element_sub(temp_sub2,w[j],w[k]);
				element_div(temp_div,temp_sub1,temp_sub2);
				element_mul(temp_delta,temp_delta,temp_div);
			}
		}
		//element_printf("temp_delta-%d = %B\n",j, temp_delta);
		element_pow_zn(temp_pow_delta,t[i],temp_delta);
		element_mul(temp_delta_pai,temp_delta_pai,temp_pow_delta);
		
		
	}	
	//element_printf("temp_delta_pai-%d = %B\n",i, temp_delta_pai);
	
	//T(X)=element_pow_zn(g2^(x^n)) ti^*
	element_init_G1(T[i],pairing);
	element_init_Zr(temp_xn,pairing);
	element_set_si(element_long_i,(signed long int)i);
	element_set_si(element_long_n,(signed long int)n);

	element_pow_zn(temp_xn,w[i],element_long_n);
	element_pow_zn(T[i],g2,temp_xn);
	element_mul(T[i],temp_delta_pai,T[i]);
	//element_printf("T-%d = %B\n",i, T[i]);

	//Di=element_pow2_zn(Di,g2,q(i),T(i),ri)
	element_init_G1(D[i],pairing);
	element_init_Zr(polynomial_q,pairing);
	element_set(polynomial_q,element_long_i);
	element_pow2_zn(D[i],g2,polynomial_q,T[i],r[i]);
	element_printf("D-%d = %B\n",i, D[i]);
}

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

