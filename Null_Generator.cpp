#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <string>
//#include <boost/multiprecision/cpp_int.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
//#include "Null_Generator.h"

using namespace std;
using namespace boost;
// This function extracts the input file information into a vector of pairs!
vector<pair<complex<double>, vector<int>>> data_extract(const string& fileName){
    vector<pair<complex<double>, vector<int>>> data;

    ifstream inputFile(fileName);
    if (!inputFile) {
        cout << "Failed to open the input file!" << endl;
        return data;
    }

    string line;
    while (getline(inputFile, line)) {
        // The first non-empty element is the coefficient
        istringstream iss(line);
        pair<complex<double>,vector<int>> linedata;

        // Extracting the complex coefficient:
        double realpart, imagpart=0;
        char sign;
        string complexPart;
        iss >> complexPart; 

        istringstream complexIss(complexPart);
        complexIss >> realpart >> imagpart;
        linedata.first = complex<double> (realpart , imagpart);

        //Extracting the integer vectors of qubits and paulis:
        string token;
        vector<int> integers;
        while (iss >> token){
            try {
                int num = boost::lexical_cast<int>(token);
                integers.push_back(num);
            } catch (const boost::bad_lexical_cast&) {
                // Ignore non-integer tokens
            }
        }
        linedata.second = integers;

        data.push_back(linedata);
        }
        inputFile.close();
    return data;
}
/*
using boost::multiprecision::cpp_int;
// ********* Computing the Zs and Ps and coefficients ******** //
typedef vector<complex<double>> Coeffs;
typedef vector<int> ZsVec;
struct PZdata {
    vector<cpp_int> Ps;
    Coeffs coeffs;
    ZsVec<int> Zs;
    ZsVec<int> Z_track;
};

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    int l = data.size();
    int z_count = 1, no_qubit = 0;
    vector<complex<double>> coefficients(l);
    vector<cpp_int> Ps;
    Coeffs coeffs;
    Indexed_Array<int> Zs(l); // This size is a max size. Some of Z[i]s might not be used!
    for (int i = 0; i<l; i++){
        coefficients[i] = data[i].first;
        vector<int> data_i = data[i].second;
        int num = 0;
        ZsVec zs;
        for (size_t j = 0; j < data_i.size() / 2; j++) {
            int pauli_j = data_i[2 * j];
            int qubit = data_i[2 * j - 1];
            if (qubit > no_qubit)
                no_qubit = qubit;
            if (pauli_j == 1) {
                num += 1 << (qubit - 1);
            } else if (pauli_j == 2) {
                num += 1 << (qubit - 1);
                coefficients[i] *= complex<double>(0, 1);
                zs.push_back(qubit);
            } else if (pauli_j == 3) {
                zs.push_back(qubit);
            }
        }
        auto it = find(Ps.begin(), Ps.end(), num);

        if(it == Ps.end()) {
            Ps.push_back(num);
            coeffs.push_back(coefficients[i]);
            // Zs.push_back(zs);
            Zs.add(z_count , zs);
            z_count++;
        } else {
            int common_p = distance(Ps.begin(), it);
            int lhz = Zs.get(common_p).size();
            int common_l = 0;
            for (int l = 0; l < lhz; l++) {
                if (Zs.get(common_p)[l] == zs) {
                    common_l = l;
                    break;
                }
            }
            if (common_l == 0) {
                Zs.add(common_p , zs);
                coeffs.push_back(coefficients[i]);
            } else {
                coeffs[common_l] += coefficients[i];
            }
        }
    }
}

vector<int> index_keep(Coeffs coeffs){
    int s = coeffs.size();
    vector<int> indices(s);
    iota(indices.begin() , indices.end() , 0);
    for (int i = 0; i < s; i++) {
        if (abs(coeffs[i]) < 1E-7) {
            indices.erase(indices.begin() + i);
        }
    }
    return indices;
}*/

// The PZComp structure is as follows:
// 1 - We have a vector of ints for Ps
// 2 - We have a vector of vector of int for Zs
// 3 - We have a vector of complexes for the coefficients
// 4 - We have a vector of vector of ints to map each Ps index to the corresponding indices of Zs and coefficients (z_tracker)!
//      More comments on this: Since Zs and coefficients get stacked on top of each other, we need to keep track of 
//          which element of Zs belongs to which Ps. The z_tracker keeps track of which Zs belong to which Ps. It is a one-to-many
//          mapping. 
using namespace boost;
int main(){
    string fileName = "test.txt";  // Replace "input.txt" with the name of your input file
    vector<pair<complex<double>, vector<int>>> data = data_extract(fileName);

    //vector<string> result = extractLinesFromFile(fileName);

    // Printing the contents of the cell array
    /*for (const auto& cell : result) {
        cout << cell << endl;
    }*/

/*     vector<int> z1 = {1,2,3};
    vector<int> z2 = {2,3};
    vector<vector<int>> Z;
    pair<complex<double>,vector<int>> data_1;
    data_1.first = complex<double> (1.0, -1.5);
    data_1.second = vector<int> {1,5,9,11};
    vector<pair<complex<double>,vector<int>>> data;
    data.push_back(data_1);
    Z.push_back(z1);
    Z.push_back(z2);
    int l = Z.size();
    for(size_t i = 0; i<l; i++){
        cout << i << " th index of Z is ";
        for(int element : Z[i]){
            cout << element << " ";
        }
        cout << endl;
    }
    cout << "That was the end of Z" << endl;

    cout << "first element of the data is " << data[0].first << endl;
    cout << "second element of the data is ";
    for(int element : data[0].second){
            cout << element << " ";
        }
    cout << endl; */
    
    
/*     string input = "1.240 1 4 5 2";
    istringstream iss(input);

    double realpart, complexpart=0;
    char sign;
    string complexPart;
    //getline(iss, complexPart, ' ');
    iss >> complexPart; 
    istringstream complexIss(complexPart);
    complexIss >> realpart >> complexpart;

    //Extracting the vectors:

    string token;
    vector<int> integers;
    while (iss >> token){
        try {
            int num = boost::lexical_cast<int>(token);
            integers.push_back(num);
        } catch (const boost::bad_lexical_cast&) {
            // Ignore non-integer tokens
        }
    }
    cout << "the real part is: " << realpart << endl; 

    std::cout << "Integers: ";
    for (const auto& num : integers) {
        std::cout << num << " ";
    }
    std::cout << std::endl; */


    return 0;

}

