#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>
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

using boost::multiprecision::cpp_int;
// ********* Computing the Zs and Ps and coefficients ******** //
typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
struct PZdata {
    vector<cpp_int> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track;
};

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    PZdata PZ_data;
    int l = data.size(),z_count = 0, no_qubit = 0;
    vector<cpp_int> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track; //This vector maps the Zs to Ps it is a many to one mapping!

    for (int i = 0; i<l; i++){
        complex<double> coeff_i = data[i].first;
        vector<int> zs_i; // For every line zs extracts the qubits on which a pauli Z acts! 
        vector<int> data_i = data[i].second; // Extracts the array of qubits and paulis for every line of input!
        int num = 0; // This variable keeps track of the index of the permutation matrix we get for each line of data!

        for (size_t j = 0; j < data_i.size() / 2; j++) {
            // Format of the input file: The 1st, 3rd, 5th, ... indicate the qubits
            int qubit = data_i[2 * j];
            // Format of the input file: The 2nd, 4th, 6th, ... indicate the paulis
            int pauli_j = data_i[2 * j + 1];

            if (qubit > no_qubit)
                no_qubit = qubit;
            if (pauli_j == 1) {
                num += 1 << (qubit - 1);
            } else if (pauli_j == 2) {
                num += 1 << (qubit - 1);
                coeff_i *= complex<double>(0, 1);
                zs_i.push_back(qubit);
            } else if (pauli_j == 3) {
                zs_i.push_back(qubit);
            }
        }
        coeffs.push_back(coeff_i);
        Zs.push_back(zs_i); // If zs is empty then no non-trivial diagonal components! (All identity operators)
        
        auto it = find(Ps.begin(), Ps.end(), num); // Look for num in the previous list of Ps (permutations)
        
        if(it == Ps.end()) {
            // num was not found in Ps, thus a new permutation matrix!
            Ps.push_back(num);
            // We will add i-th element of coeffs and Zs to be associated with the current Ps!
            Z_track.push_back(vector<int> {i});
        } else {
            // num was found in Ps, and P_index will be the index of Ps that matched num.
            int P_index = std::distance(Ps.begin(), it);

            vector<int> z_indices = Z_track[P_index];

            bool z_found = false;
            for (int k = 0; k < z_indices.size(); k++){
                if (Zs[z_indices[k]] == zs_i){
                    coeffs[P_index] += coeff_i;
                    z_found = true;
                    break;
                }
            }
            // If the z array is new, we add it to the Z_track for the Ps associated with num! 
            if (!z_found){
                Z_track[P_index].push_back(i);
            }
        }
    }
    // Throw away the zero coefficients:
    vector<cpp_int> Ps_kept;
    ZVecs Z_track_kept;

    for(int k = 0; k < Z_track.size(); k++){
        vector<int> ztrack_k;
        for(int l = 0; l < Z_track[k].size(); l++){
            int k_l_index = Z_track[k][l];
            complex<double> coeff_k_l = coeffs[k_l_index];
            if (abs(coeff_k_l) > 1e-8){
                //zero_coeffs_for_k.push_back(l);
                ztrack_k.push_back(k_l_index);
            }
        }
        if(ztrack_k.size() > 0){
            Z_track_kept.push_back(ztrack_k);
            Ps_kept.push_back(Ps[k]);
        }
    }

    // Sorting everything based on indices of Ps:
    vector<pair<cpp_int, int>> indexedP;
    for(int i = 0; i < Ps_kept.size(); i++) {
        indexedP.emplace_back(Ps_kept[i], i);
    }

    vector<cpp_int> Ps_sorted;
    ZVecs Z_track_sorted;

    sort(indexedP.begin() , indexedP.end());
    for(int i = 0; i < Ps_kept.size(); i++){
            int index = indexedP[i].second;
            Ps_sorted.push_back(indexedP[i].first);
            Z_track_sorted.push_back(Z_track_kept[index]);
    }

    PZ_data.coeffs = coeffs; //Coeffs and Zs are kept as the original data
    PZ_data.Ps = Ps_sorted;
    PZ_data.Zs = Zs;
    PZ_data.Z_track = Z_track_sorted;

    return PZ_data;
}

/*
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
    PZdata PZ_data = PZcomp(data);
    vector<cpp_int> Ps = PZ_data.Ps;
    vector<complex<double>> coefficients = PZ_data.coeffs;
    vector<vector<int>> Z_track = PZ_data.Z_track;
    vector<vector<int>> Zs = PZ_data.Zs;

    return 0;
}

