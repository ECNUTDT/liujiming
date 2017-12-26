
#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include "json/json.h"
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iterator>

using namespace std;


const string configPath = "../../data/config/config";
const string PPPath = "../../data/setup_data/PP";
const string MKPath = "../../data/setup_data/MK";
const int degree = 5;

int N,M;
element_t g,g1,g2;
element_t v;
element_t y;


//g++ -o extract extract.cpp -L. -lgmp -lpbc -ljson
//./extract <../../data/param/a.param 

//
unsigned char* transfer(element_t t){

  int leng = element_length_in_bytes(t);
  unsigned char *p = new unsigned char[leng + 1];
  int writeLength = element_to_bytes_compressed(p,t);
  return p;
}
void getJsonValueNKey(Json::Value value,int tag,pairing_t pairing,element_t *arrN,element_t *arrM){

Json::Value::Members members;
members = value.getMemberNames();
int i = 0;
for (Json::Value::Members::iterator iterMember = members.begin(); iterMember != members.end(); iterMember++){  
  string strKey = *iterMember;
	cout<<"strKey:"+strKey+"\n";

    if (tag == 0) {
        element_t temp;
        element_init_G1(temp, pairing);
        string str = value[strKey.c_str()].asString();
        char * p;
        int len = str.length();
        p = new char[len + 1];
        strcpy(p, str.c_str());
        unsigned char * sb = (unsigned char * ) p;
        element_from_bytes_compressed(temp, sb);
        int pos = strKey.find("-");
        if (pos == -1) {
            if (strKey == "g") {
                element_init_G1(g, pairing);
                element_from_bytes_compressed(g, sb);
            } else if (strKey == "g1") {
                element_init_G1(g1, pairing);
                element_from_bytes_compressed(g1, sb);
            } else if (strKey == "g2") {
                element_init_G1(g2, pairing);
                element_from_bytes_compressed(g2, sb);
            } else if (strKey == "v") {
                element_init_G1(v, pairing);
                element_from_bytes_compressed(v, sb);
            }
        } else {
            string type = strKey.substr(0, pos);
            string num = strKey.substr(pos + 1, strKey.size());
            stringstream ss;
            ss << num;
            int index;
            ss >> index;
            if (type == "t") {
                element_init_G1(arrN[index - 1], pairing);
                element_from_bytes_compressed(arrN[index - 1], sb);
            } else if (type == "v") {
                element_init_G1(arrM[index - 1], pairing);
                element_from_bytes_compressed(arrM[index - 1], sb);
            }
        }
    }
    if (tag == 1) {
        element_init_Zr(y, pairing);
        string strVal = value[strKey.c_str()].asString();
        mpz_t gmp_y;
        const char * chargmp_y = strVal.c_str();
        mpz_init_set_str(gmp_y, chargmp_y, 10);
        element_set_mpz(y, gmp_y);
    }
    if (tag == 2) {
        int iVal = value[strKey.c_str()].asInt();
        if (strKey == "N") {
            N = iVal;
        } else if (strKey == "M") {
            M = iVal;
        }
    }
    if (tag == 3) {
       	element_t temp;
        element_init_G1(temp, pairing);
        string str = value[strKey.c_str()].asString();
        char * p;
        int len = str.length();
        p = new char[len + 1];
        strcpy(p, str.c_str());
        unsigned char * sb = (unsigned char * ) p;
        element_from_bytes_compressed(temp, sb);
        int pos = strKey.find("-");
        if (pos == -1) {
            if (strKey == "g") {
                element_init_G1(g, pairing);
                element_from_bytes_compressed(g, sb);
            } else if (strKey == "g1") {
                element_init_G1(g1, pairing);
                element_from_bytes_compressed(g1, sb);
            } else if (strKey == "g2") {
                element_init_G1(g2, pairing);
                element_from_bytes_compressed(g2, sb);
            } else if (strKey == "v") {
                element_init_G1(v, pairing);
                element_from_bytes_compressed(v, sb);
            }
        } else {
            string type = strKey.substr(0, pos);
            string num = strKey.substr(pos + 1, strKey.size());
            stringstream ss;
            ss << num;
            int index;
            ss >> index;
            if (type == "D") {
                element_init_G1(arrN[index - 1], pairing);
                element_from_bytes_compressed(arrN[index - 1], sb);
            } else if (type == "d") {
                element_init_G1(arrM[index - 1], pairing);
                element_from_bytes_compressed(arrM[index - 1], sb);
            }
        }
    }

}

}
string readFile(string path){
	stringstream ss;
	fstream readPathData(path.c_str());
	ss << readPathData.rdbuf();
	string data = ss.str();
	ss.clear();
	ss.str("");
	readPathData.close();
	return data;
}




int main(int argc, char **argv) {

//clain element
element_t A;
element_t z;


//int n,m,d;
//scanf("%d %d", &n,&m);
//int n=10,m=10;



//--------READ FILE--------------
//=============read data from file end===================
string config = readFile(configPath);
string PP = readFile(PPPath);
string MK = readFile(MKPath);


pairing_t pairing;
pbc_demo_pairing_init(pairing, argc, argv);

Json::Reader reader;
Json::Reader reader1;
Json::Reader reader2;

Json::Value value;
Json::Value value1;
Json::Value value2;


if(!reader2.parse(config,value2)){
	cout<<"cannot read config file!\n";
	return 0;
}
getJsonValueNKey(value2,2,pairing,NULL,NULL);
element_t *t = new element_t[N];
element_t *V = new element_t[M];


if(!reader.parse(PP,value)){
	cout<<"cannot read PP file!\n";
	return 0;
}
getJsonValueNKey(value,0,pairing,t,V);


if(!reader1.parse(MK,value1)){
	cout<<"cannot read MK file!\n";
	return 0;
}
getJsonValueNKey(value1,1,pairing,NULL,NULL);




//--------EXTRACT--------------


element_t Z[M];

element_t d[N+1];
element_t D[N+1];

element_t w[N+1];
element_t r[N+1];
element_t T[N+1];
element_t q[degree];

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

element_init_Zr(element_long_i,pairing);
element_init_Zr(element_long_n,pairing);

for(int i=0;i<N+1;i++){
	
	
	element_init_Zr(w[i],pairing);
	element_random(w[i]);

	element_init_Zr(r[i],pairing);
	element_random(r[i]);
}



//delta
element_init_Zr(temp_delta,pairing);
element_init_G1(temp_pow_delta,pairing);
element_init_G1(temp_delta_pai,pairing);
element_init_Zr(temp_sub1,pairing);
element_init_Zr(temp_sub2,pairing);
element_init_Zr(temp_div,pairing);



for(int i=0;i<N+1;i++){
	//di=element_pow_zn(di,g,-ri) 
	element_init_G1(d[i],pairing);
	element_init_Zr(temp_neg_ri,pairing);
	element_neg(temp_neg_ri,r[i]);
	element_pow_zn(d[i],g,temp_neg_ri);
	element_printf("d-%d = %B\n",i, d[i]);

	//π△(x) 
	element_set1(temp_delta_pai);
	for(int j=0;j<N+1;j++){
		element_set1(temp_delta);
		element_set0(temp_sub1);
		element_set0(temp_sub2);
		//element_printf("temp_delta-%d = %B\n",j, temp_delta);
		for(int k=0;k<N+1;k++){
			
			
			if(k == j) continue;
			else if(k==i) continue;
			else {
				element_sub(temp_sub1,w[i],w[k]);
				element_sub(temp_sub2,w[j],w[k]);
				element_div(temp_div,temp_sub1,temp_sub2);
				element_mul(temp_delta,temp_delta,temp_div);
			}
			
		}
		element_pow_zn(temp_pow_delta,t[i],temp_delta);//t[i]<i:1-n+1>
		element_mul(temp_delta_pai,temp_delta_pai,temp_pow_delta);
				
	}	

	//T(X)=g2^(x^n)*π△(x) 
	element_init_G1(T[i],pairing);
	element_init_Zr(temp_xn,pairing);
	element_set_si(element_long_i,(signed long int)i);
	element_set_si(element_long_n,(signed long int)N);
	element_pow_zn(temp_xn,w[i],element_long_n);
	element_pow_zn(T[i],g2,temp_xn);
	element_mul(T[i],temp_delta_pai,T[i]);

	//Di=element_pow2_zn(Di,g2,q(i),T(i),ri)
	element_init_G1(D[i],pairing);
	element_init_Zr(polynomial_q,pairing);
	element_set0(polynomial_q);

	//q(i)=r1*i^(d-1)+r2*(d-2)+.....+y
	for(int l=0;l<degree;l++){
		element_init_Zr(q[l],pairing);
		element_random(q[l]);
		if(l==0) element_set(q[l],y);
		if(l==degree-1) {
			while(element_is0(q[l])) element_random(q[l]);
		}
		
		element_pow_zn(temp_xn,w[i],element_long_n);
		element_mul(temp_xn,temp_xn,q[l]);
		element_add(polynomial_q,polynomial_q,temp_xn);
	}
	element_pow2_zn(D[i],g2,polynomial_q,T[i],r[i]);
	element_printf("D-%d = %B\n",i, D[i]);
	
}



/*
//write file Dd	
FILE *fc;
fc=fopen("../..//data/extract_data/Dd","w+");
element_fprintf(fc,"{");
for(int i=0;i<N+1;i++){
	if(i>0) element_fprintf(fc,",");
	element_fprintf(fc,"\r\"d-%d\":\"%B\"",i, d[i]);
	element_fprintf(fc,",\r\"D-%d\":\"%B\"",i, D[i]);
}
element_fprintf(fc,"\r}");
fclose(fc);
*/

Json::Value root;
Json::StyledWriter sw;
ofstream os;
os.open("../../data/extract_data/Dd");
for (int i = 0; i < N+1; i++) {
        string index = "";
        stringstream st;
        st << (i + 1);
        st >> index;
        index = "d-" + index;
        root[index] = Json::Value((char * ) transfer(d[i]));
}

for (int i = 0; i < N+1; i++) {
        string index = "";
        stringstream st;
        st << (i + 1);
        st >> index;
        index = "D-" + index;
        root[index] = Json::Value((char * ) transfer(D[i]));
}
os << sw.write(root);
os.close();





//clear element
element_clear(g);
element_clear(g1);
element_clear(g2);
element_clear(y);
element_clear(v);


element_clear(element_long_i);
element_clear(element_long_n);
element_clear(temp_neg_ri);
element_clear(temp_pow_delta);
element_clear(temp_delta);
element_clear(temp_delta_pai);
element_clear(temp_xn);
element_clear(temp_sub1);
element_clear(temp_sub2);
element_clear(temp_div);

delete []t;
delete []V;


pairing_clear(pairing);

return 0;

}


