
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include "json/json.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

using namespace std;

const int N=10;
const int M=10;


//g++ -o setup setup.c -L. -lgmp -lpbc -ljson
//./setup <../../data/param/a.param 

unsigned char* transfer(element_t t){

  int leng = element_length_in_bytes(t);
  unsigned char *p = new unsigned char[leng + 1];
  int writeLength = element_to_bytes_compressed(p,t);
  return p;
}


int main(int argc, char **argv) {

//clain element
pairing_t pairing; 
element_t g,g1,g2;
element_t A;
element_t v;
element_t z,y;

element_t T[N+1];
element_t Z[M];
element_t V[M];

pbc_demo_pairing_init(pairing, argc, argv);

element_init_G1(g,pairing);
element_init_G1(g1,pairing);
element_init_G1(g2,pairing);
element_init_G1(v,pairing);
element_init_Zr(y,pairing);
element_init_Zr(z,pairing);
element_init_GT(A,pairing);

element_random(g);
element_random(g2);
element_random(y);
element_random(z);


element_pow_zn(g1,g,y);
element_pow_zn(v,g,z);
element_pairing(A,g1,g2);

for (int i=0;i<N+1;i++){
	element_init_G1(T[i],pairing);
	element_random(T[i]);
	element_printf("T[%d] = %B\n",i, T[i]);
}

for (int i=0;i<M;i++){
	element_init_Zr(Z[i],pairing);
	element_random(Z[i]);
}

for (int i=0;i<M;i++){
	element_init_G1(V[i],pairing);
	element_pow_zn(V[i],g,Z[i]);
	element_printf("V[%d] = %B\n",i, V[i]);
}

element_printf("A= %B\n",A);

//write file PP : g,g1,g2,v,v1...vm,t1...tn,A
/*
FILE *fc;
fc=fopen("../..//data/setup_data/PP","w+");
element_fprintf(fc,"{\"g\":\"%B\"",g);
element_fprintf(fc,",\r\"g1\":\"%B\"",g1);
element_fprintf(fc,",\r\"g2\":\"%B\"",g2);
element_fprintf(fc,",\r\"z\":\"%B\"",z);
element_fprintf(fc,",\r\"v\":\"%B\"",v);
element_fprintf(fc,",\r\"A\":\"%B\"",A);
for (int i=0;i<M;i++){
	element_fprintf(fc,",\r\"V-%d\":\"%B\"",i,V[i]);
}
for (int i=0;i<N;i++){
	element_fprintf(fc,",\r\"T-%d\":\"%B\"",i,T[i]);
}
element_fprintf(fc,"\r}");
fclose(fc);
*/

//jsoncpp
Json::Value root;
Json::StyledWriter sw;
ofstream os;
os.open("../../data/setup_data/PP");

root["g"] = Json::Value((char*) transfer(g));
root["g1"] = Json::Value((char*) transfer(g1));
root["g2"] = Json::Value((char*) transfer(g2));


    

root["v"] = Json::Value((char*) transfer(v));
//root["z"] = Json::Value((char*) transfer(z));
//root["A"] = Json::Value((char*) transfer(A));

    for (int i = 0; i < M; i++) {
        string index = "";
        stringstream st;
        st << (i + 1);
        st >> index;
        index = "v-" + index;
        root[index] = Json::Value((char * ) transfer(V[i]));
    }
    for (int i = 0; i < N+1; i++) {
        string index = "";
        stringstream st;
        st << (i + 1);
        st >> index;
        index = "t-" + index;
        root[index] = Json::Value((char * ) transfer(T[i]));
    }
    os << sw.write(root);
    os.close();

    //write file MK :y
    /*
	fc = fopen("../../data/setup_data/MK", "w+");
    element_fprintf(fc, "{\"y\":\"%B\"}", y);
    fclose(fc);
	*/

    Json::Value root1;
    Json::StyledWriter sw1;
    ofstream os1;
    

    //create file MK to keep the variables 
    if (freopen("../../data/setup_data/MK", "w", stdout) == NULL) {
        fprintf(stderr, "error2\n");
    }
    printf("{\n");
    element_printf("\"y\":\"%B\"\n", y);
    //end writing data to file MK 
    printf("}\n");

    fclose(stdout);

os1.open("../../data/config/config");
root1["N"] = Json::Value(N);
root1["M"] = Json::Value(M);
os1 << sw.write(root1);
os1.close();


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
