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


//g++ -o sign sign.cpp -L. -lgmp -lpbc -ljson
//./sign <../../data/param/a.param 

const string configPath = "../../data/config/config";
const string PPPath = "../../data/setup_data/PP";
const string MKPath = "../../data/setup_data/MK";
const int degree = 5;

int N,M;

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
pairing_t pairing; 
element_t g,g1,g2;
element_t A;
element_t v;
element_t z,y;



//int n,m,d;
//scanf("%d %d", &n,&m);
int n=10,m=10,degree=5;


element_t t[n+1];
element_t Z[m];
element_t V[m];
element_t d[n+1];
element_t D[n+1];

element_t w[n+1];
element_t r[n+1];
element_t T[n+1];
element_t q[degree];

element_t miu[m];
element_t s[n+1];
element_t S[n+1][3];
element_t temp_G1;
element_t pai_v;

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

char* str="76794694781183289556915698583624103921063413964367385662225670703628416913618214459476707216239711143815363248997679867517448662092880568570217328577056715082897505846747235950749109596558737919561125054842757776548508220404079200400120171875332425702500370972342317374668041488772875181617059537446057415909";
char* stry="479669737328012410214860217391068436330049605306";







int readB=element_from_bytes(g,str);

int readY=element_from_bytes(y,stry);

int readB=element_to_bytes(str,g);


element_printf("====read======%d====%d=====\n", readB,readY);

element_printf("g = %B\n", g);

element_printf("y = %B\n", y);



return 0;


//--------SETUP--------------
element_random(g);
element_printf("g = %B\n", g);


element_random(g2);
//element_printf("g2 = %B\n", g2);

element_random(y);
//element_printf("y = %B\n", y);


element_random(z);
//element_printf("z = %B\n", z);


for (int i=0;i<n+1;i++){
	element_init_G1(t[i],pairing);
	element_random(t[i]);
	//element_printf("t[%d] = %B\n",i, t[i]);
}

for (int i=0;i<m;i++){
	element_init_Zr(Z[i],pairing);
	element_random(Z[i]);
}

for (int i=0;i<m;i++){
	element_init_G1(V[i],pairing);
	element_pow_zn(V[i],g,Z[i]);
	//element_printf("V[%d] = %B\n",i, V[i]);
}


element_pow_zn(g1,g,y);
//element_printf("g1 = %B\n", g1);

element_pow_zn(v,g,z);
//element_printf("v = %B\n", v);

element_pairing(A,g1,g2);
//element_printf("A = %B\n", A);


//--------EXTRACT--------------
for(int i=0;i<n+1;i++){
	
	
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



for(int i=0;i<n+1;i++){
	//di=element_pow_zn(di,g,-ri) 
	element_init_G1(d[i],pairing);
	element_init_Zr(temp_neg_ri,pairing);
	element_neg(temp_neg_ri,r[i]);
	element_pow_zn(d[i],g,temp_neg_ri);
	element_printf("d-%d = %B\n",i, d[i]);

	//π△(x) 
	element_set1(temp_delta_pai);
	for(int j=0;j<n+1;j++){
		element_set1(temp_delta);
		element_set0(temp_sub1);
		element_set0(temp_sub2);
		//element_printf("temp_delta-%d = %B\n",j, temp_delta);
		for(int k=0;k<n+1;k++){
			
			
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
	element_set_si(element_long_n,(signed long int)n);
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

//--------EXTRACT--------------
element_init_G1(temp_G1,pairing);
element_init_G1(pai_v,pairing);
element_set1(pai_v);
for(int i=0;i<m;i++){
	element_init_Zr(miu[i],pairing);
	int num=rand()%2;
	element_set_si(miu[i],(signed long int)num);
	element_pow_zn(temp_G1,V[i],miu[i]);
	element_mul(pai_v,pai_v,temp_G1);
	element_printf("\n=====%d=====\n%B",i,pai_v);
	//if(i==m-3) continue;
	
}



element_t temp_neg_si;
element_init_Zr(temp_neg_si,pairing);
for (int i=0;i<n+1;i++){
	element_init_Zr(s[i],pairing);
	element_random(s[i]);
	element_neg(temp_neg_si,s[i]);
	for(int j=0;j<3;j++){
		element_init_G1(S[i][j],pairing);
	}
		element_mul(temp_G1,pai_v,v);
		element_pow_zn(temp_G1,temp_G1,s[i]);
		element_mul(temp_G1,temp_G1,D[i]);
		element_set(S[i][0],temp_G1);
		element_set(S[i][1],d[i]);
		element_pow_zn(temp_G1,g,temp_neg_si);
		element_set(S[i][2],temp_G1);
}



//write file S	
FILE *fc;
fc=fopen("../../data/sign_data/Sign","w+");
element_fprintf(fc,"{");
for(int i=0;i<n+1;i++){
	element_fprintf(fc,"\r\"S-%d-1\":\"%B\"",i, S[i][0]);
	element_fprintf(fc,"\r\"S-%d-2\":\"%B\"",i, S[i][1]);
	element_fprintf(fc,"\r\"S-%d-3\":\"%B\"",i, S[i][2]);
}
element_fprintf(fc,"\r}");
fclose(fc);



//clear element
element_clear(g);
element_clear(g1);
element_clear(g2);
element_clear(y);
element_clear(z);
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



pairing_clear(pairing);

return 0;

}

