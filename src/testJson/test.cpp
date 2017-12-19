
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include<json/json.h>
using namespace std;
int main(){
	Json::Value root;
	Json::StyledWriter fast;
	root["DateTime"]=("23");
	root["Name"]=("hello");
	root["class"]=("12");
	cout<<fast.write(root)<<endl;

	ofstream os;
	os.open("test_data");

os << fast.write(root);  
  os.close();  


	return 0;

}
