#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#include "utils.h"
using namespace std;


int main(){

    double L_S_min = 6.296; //nm
    string line;
    ifstream infile("checking.gro");
    if (!infile.is_open()) {
        cerr << "Error opening file" << endl;
        return 1;
    }

    vector<float> X1;
    while (getline(infile, line)) {
        istringstream iss(line);
        vector<string> tokens;
        string token;

        while (iss >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() == 6) {
            X1.push_back(stof(tokens[3]));
        } else if (tokens.size() == 5) {
            X1.push_back(stof(tokens[2]));
        }
    }

    infile.close();


    int H2Onumber = 0;
    for(int i = 0; i < X1.size();i+=4){

        if(X1[i] < L_S_min){
            H2Onumber++;

        }
    }

    cout << H2Onumber << endl;


    return 0;
}