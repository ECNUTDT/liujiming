
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
const string DdPath = "../../data/extract_data/Dd";
const int degree = 5;

int N,M;
element_t g,g1,g2;
element_t v;
element_t y;


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


//--------READ FILE--------------
//=============read data from file end===================
string config = readFile(configPath);
string PP = readFile(PPPath);
string MK = readFile(MKPath);
string Dd = readFile(DdPath);

pairing_t pairing;
pbc_demo_pairing_init(pairing, argc, argv);

Json::Reader reader;
Json::Reader reader1;
Json::Reader reader2;
Json::Reader reader3;
Json::Value value;
Json::Value value1;
Json::Value value2;
Json::Value value3;

if(!reader2.parse(config,value2)){
	cout<<"cannot read config file!\n";
	return 0;
}
getJsonValueNKey(value2,2,pairing,NULL,NULL);
element_t *t = new element_t[N];
element_t *V = new element_t[M];
element_t *D = new element_t[N+1];
element_t *d = new element_t[N+1];

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

if(!reader3.parse(Dd,value3)){
	cout<<"cannot read Dd file!\n";
	return 0;
}

getJsonValueNKey(value3,3,pairing,D,d);




element_t w[N+1];
element_t r[N+1];
element_t T[N+1];
element_t q[degree];

element_t miu[M];
element_t s[N+1];
element_t S[N+1][3];
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



//--------SIGN--------------
element_init_G1(temp_G1,pairing);
element_init_G1(pai_v,pairing);
element_set1(pai_v);
for(int i=0;i<M;i++){
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
for (int i=0;i<N+1;i++){
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
for(int i=0;i<N+1;i++){
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

