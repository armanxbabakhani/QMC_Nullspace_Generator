#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>
//#include "Null_Generator.h"

using namespace std;
using namespace boost;
using boost::multiprecision::cpp_int;

// *****************************    Functions ******************************************* //

void printMatrix(const vector<vector<int>>& matrix) {
    int m = matrix.size();
    int n = matrix[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

pair<vector<int> , vector<int>> int_to_bin(cpp_int num , int no_qubits) {
    vector<int> binaryVector , binaryVector_padded;
    pair<vector<int> , vector<int>> binary_combo;
    // The binary vectors are created from left to right. For example, the number 6 --> [1 1 0] , 14 --> [1 1 1 0]
    while (num > 0) {
        int remainder = static_cast<int> (num % 2);
        binaryVector.push_back(remainder);
        num /= 2;
    }

    int no_of_padding = no_qubits - binaryVector.size();
    reverse(binaryVector.begin() , binaryVector.end());
    binaryVector_padded = binaryVector;
    // Padding zeros from the left
    for(int i=0; i < no_of_padding; i++){
        binaryVector_padded.insert(binaryVector_padded.begin() , 0);
    }
    binary_combo.first = binaryVector_padded;
    binary_combo.second = binaryVector;
    return binary_combo;
}

// This function is for binary addition (XOR)
vector<int> GF2_add(vector<int> vec1 , vector<int> vec2){
    vector<int> result;
    if (vec1.size() != vec2.size()){
        cerr << "The size of the added vectors are differen!" << endl;
    }

    for(int i = 0; i < vec1.size(); i++){
        result.push_back((vec1[i] + vec2[i])%2);
    }
    return result;
}

// Function to compute the mod 2 nullspace
vector<vector<int>> Null2(const vector<vector<int>>& matrix) {
    int numRows = matrix.size();
    int numCols = matrix[0].size();
    vector<vector<int>> matrix_RE = matrix; // Copy the original matrix
    vector<int> nullspaceVector(numCols, 0);
    vector<int> marked_rows;

    // Gaussian Elimination (GF(2)) to obtain matrix_GE
    for(int j = 0; j < numCols; j++){
        for(int i=0; i < numRows; i++){
            if(matrix_RE[i][j] == 1){
                marked_rows.push_back(i);
                for(int k = 0; k < numCols; k++ ){
                    if(k != j){
                        if(matrix_RE[i][k] == 1){
                            // Adding column j to k!
                            for(int l = 0; l < matrix_RE.size(); l++){
                                matrix_RE[l][k] = matrix_RE[l][j] ^ matrix_RE[l][k];
                            }
                        }
                    }
                }
                break;
            }
        }
    }
    // sort the marked row
    std::sort(marked_rows.begin(), marked_rows.end());

    // Remove consecutive duplicates (taking only unieque marked rows!)
    auto it = unique(marked_rows.begin(), marked_rows.end());
    marked_rows.erase(it, marked_rows.end());

    vector<vector<int>> marked_matrix;
    // Make the marked matrix:
    for(int i = 0 ; i < marked_rows.size() ; i++){
        marked_matrix.push_back(matrix_RE[marked_rows[i]]);
    }
    // Go over the marked rows and find the nullspace!
    // The rows that are unmarked will define the nullspace
    vector<vector<int>> nullspaceBasis;
    for(int i = 0; i < numRows; i++){
        auto iter = find(marked_rows.begin(), marked_rows.end(), i); 
        if(iter == marked_rows.end()){ // Finding unmarked (independent) rows
            // Step 1 - Then find the ones in that row!
            // 2 - Go over each one found in that row and find the corresponding
            // row which has that that index of one in that column

            // Step 1 : 
            vector<int> row_i = matrix_RE[i] , nullvector_j(numRows, 0);
            nullvector_j[i] = 1;
            for(int j = 0; j < numCols; j++){
                if(row_i[j] == 1){
                    // Go over the marked rows and extract the null space:
                    for(int k = 0; k < marked_matrix.size(); k++){
                        if(marked_matrix[k][j] == 1){
                            nullvector_j[marked_rows[k]] = 1;
                        } 
                    }
                }
            }
            nullspaceBasis.push_back(nullvector_j);
        }
    }
    return nullspaceBasis;
}

// This function minimizes the size of the cycles (nullspace basis vectors with least 1s):
int int_sum(vector<int> vec){
    int sum = 0;
    for(int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}

void eig_minimize(vector<vector<int>>& null_eigs){
    int null_size = null_eigs.size();
    vector<int> null_n(null_size);
    for(int i = 0; i < null_size; i++){
        null_n[i] = int_sum(null_eigs[i]);
    }
    vector<vector<int>> null_eigs_min , null_eigs_high;
    vector<int> null_eigs_min_ind , null_eigs_high_ind;
    for(int i = 0; i < null_size; i++){
        if(null_n[i] == 3){
            null_eigs_min.push_back(null_eigs[i]);
            null_eigs_min_ind.push_back(i);
        }
        else{
            null_eigs_high.push_back(null_eigs[i]);
            null_eigs_high_ind.push_back(i);
        }
    }

    int high_size = null_eigs_high_ind.size();
    for(int k = 0; k < high_size; k++){
        vector<int> high_eig_k = null_eigs[null_eigs_high_ind[k]];
        int null_k = int_sum(high_eig_k);
        for(int l = 0; l < null_eigs_min_ind.size(); l++){
            vector<int> high_to_min = GF2_add(high_eig_k , null_eigs[null_eigs_min_ind[l]]);
            int low_k = int_sum(high_to_min);
            if(low_k == 3){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                null_eigs_min_ind.push_back(null_eigs_high_ind[k]);
                std::sort(null_eigs_min_ind.begin() , null_eigs_min_ind.end());

                //null_eigs_high_ind.erase(null_eigs_high_ind.begin() + k);
                break;
            }
            else if(low_k < null_k){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                high_eig_k = high_to_min;
            }
        }
    }
}

// This function converts the array of integers into a corresponding binary string (used to convert the indices of Z into string of bitsets)
string int_to_str(vector<int> Z){
    string Z_string = "";
    if(Z.size() < 1){
        return "0";
    }
    std::sort(Z.begin() , Z.end());
    int Z_size = Z.size();
    int count = 1, ind_z_count = 0, max_z = Z[Z_size-1];

    while(ind_z_count < Z_size){
        if(Z[ind_z_count] == count){
            Z_string = "1" + Z_string;
            ind_z_count++;
        }
        else{
            Z_string = "0" + Z_string;
        }
        count++;
    }
    return Z_string;
}

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

// ------------- Functions to compute the Zs and Ps and coefficients ----------
typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
struct PZdata {
    vector<cpp_int> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track;
    int no_qubit;
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
        cpp_int num("0"); // This variable keeps track of the index of the permutation matrix we get for each line of data!

        for (size_t j = 0; j < data_i.size() / 2; j++) {
            // Format of the input file: The 1st, 3rd, 5th, ... indicate the qubits
            int qubit = data_i[2 * j];
            // Format of the input file: The 2nd, 4th, 6th, ... indicate the paulis
            int pauli_j = data_i[2 * j + 1];

            if (qubit > no_qubit)
                no_qubit = qubit;
            if (pauli_j == 1) {
                num += boost::multiprecision::pow(cpp_int(2), qubit-1);
            } else if (pauli_j == 2) {
                //cpp_int perm_term = 1 << (qubit - 1);
                //num += perm_term;
                num += boost::multiprecision::pow(cpp_int(2), qubit-1);
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

    std::sort(indexedP.begin() , indexedP.end());
    for(int i = 0; i < Ps_kept.size(); i++){
            int index = indexedP[i].second;
            Ps_sorted.push_back(indexedP[i].first);
            Z_track_sorted.push_back(Z_track_kept[index]);
    }

    PZ_data.coeffs = coeffs; //Coeffs and Zs are kept as the original data
    PZ_data.Ps = Ps_sorted;
    PZ_data.Zs = Zs;
    PZ_data.Z_track = Z_track_sorted;
    PZ_data.no_qubit = no_qubit;

    return PZ_data;
}

// ************************************************************************************************************** //
// -------------------------------------------------------------------------------------------------------------- //


using namespace boost;
int main(int argc , char* argv[]){
    string fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> data = data_extract(fileName);

    // Unpacking the data from the input file "fileName"
    PZdata PZ_data = PZcomp(data);
    vector<cpp_int> Ps = PZ_data.Ps;
    vector<complex<double>> coefficients = PZ_data.coeffs;
    vector<vector<int>> Z_track = PZ_data.Z_track;
    vector<vector<int>> Zs = PZ_data.Zs;
    vector<string> Zs_string;
    int no_qubit = PZ_data.no_qubit;

    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs.size();i++){
        Zs_string.push_back(int_to_str(Zs[i]));
    }

    // Convert Ps indices into vector of ints
    vector<cpp_int> Ps_nontrivial = Ps;
    if(Ps[0] == 0){
        Ps_nontrivial.erase(Ps_nontrivial.begin());
    }
    vector<vector<int>> Ps_binary , Ps_binary_unpad;
    for(int i = 0; i < Ps_nontrivial.size(); i++){
        pair<vector<int>, vector<int>> bin_i_combo = int_to_bin(Ps_nontrivial[i], no_qubit);
        vector<int> bin_i_pad = bin_i_combo.first;
        vector<int> bin_i_unpad = bin_i_combo.second;
        Ps_binary.push_back(bin_i_pad);
        Ps_binary_unpad.push_back(bin_i_unpad);
    }

    // Minimizing the size of the fundamental cycless
    vector<vector<int>> nullspace = Null2(Ps_binary);
    eig_minimize(nullspace); 
    cout << "Calculations done!" << endl;
    cout << "Making the output files ... " << endl;

    int no_ps = Ps_nontrivial.size();
    int nullity = nullspace.size();
    string output_h = fileName.substr(0, fileName.find_last_of(".")) + ".h" , output_cycle = fileName.substr(0, fileName.find_last_of(".")) + "_cycles.txt";
    ofstream output(output_h) , output_c(output_cycle);

    // ********************************************************************** //
    // -------------------- Creating the the .h file ------------------------ //
    if(output.is_open()){
        output << "#define N        " << no_qubit << endl;
        output << "#define Nop      " << no_ps << endl;
        output << "#define Ncycles  " << nullity << endl;
        output << endl;
        
        // ---------------- Permutation matrices and cycles --------------------- //
        // The permutation bitsets
        output << "std::bitset<N> P_matrix[Nop] = {";
        for(int i = 0; i < no_ps; i++){
            output << "std::bitset<N>(\"";
            for(int j=0; j < Ps_binary_unpad[i].size(); j++){
                output << Ps_binary_unpad[i][j];
            }
            output << "\")";
            if(i < no_ps - 1){
                output << ", ";
            }
        }
        output << "};" << endl;

        // The cycle bitsets
        output << "std::bitset<Nop> cycles[Ncycles] = {";
        for(int i = 0; i < nullity; i++){
            output << "std::bitset<Nop>(\"";
            for(int j=0; j < no_ps; j++){
                output << nullspace[i][j];
            }
            output << "\")";
            if(i < nullity-1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << endl;

        // ---------------------------------------------------------------------- //
        // -------------------------- Diagonal terms ---------------------------- //
        output << "const int D0_size = " << Z_track[0].size() << ";" << endl;
        output << "double D0_coeff[D0_size] = {";
        for(int i=0;i<Z_track[0].size();i++){
            complex<double> c_ij = coefficients[Z_track[0][i]];
            if(abs(c_ij.imag()) > 1e-7 ){
                if(c_ij.imag() < 0){
                output << c_ij.real() << c_ij.imag() << "i";
                }
                else{
                    output << c_ij.real() << "+" << c_ij.imag() << "i";
                }
            }
            else{
                output << c_ij.real();
            }
            if(i < Z_track[0].size() - 1){
                output << ", "; 
            }
        }
        output << "};" << endl;
        output << "std::bitset<N> D0_product[D0_size] = {";
        for(int i = 0; i < Z_track[0].size(); i++){
            vector<int> Zs_i = Zs[Z_track[0][i]];
            output << "std::bitset<N>(\"" << int_to_str(Zs_i) << "\")";
            if(i < Z_track[0].size() - 1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << endl;

        // ------------------------------------------------------------------------ //
        // --------------------------- Off-Diagonal terms ------------------------- //
        // Finding D_maxsize:
        int D_max = 0 , D_start;
        vector<int> D_size;
        if(Ps_nontrivial.size() == Ps.size()) // In case of no diagonal term (i.e. D_0 = 0)
        {
            D_start = 0;
        }else{
            D_start = 1;
        }
        for(int i = D_start; i < Z_track.size(); i++){
            int z_size_i = Z_track[i].size();
            D_size.push_back(z_size_i);
            if(z_size_i > D_max){
                D_max = z_size_i;
            }
        }
        output << "const int D_maxsize = " << D_max << " ;" << endl;
        output << "int D_size[Nop] = {";
        for(int i=0;i<no_ps;i++){
            output << D_size[i];
            if(i < no_ps - 1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << "double D_coeff[Nop][D_maxsize] = {";
        for(int i = D_start; i < Z_track.size() ; i++){
            output << "{";
            for(int j = 0; j < Z_track[i].size(); j++){
                complex<double> c_ij = coefficients[Z_track[i][j]];
                if(abs(c_ij.imag()) > 1e-7 ){
                    if(c_ij.imag() < 0){
                        output << c_ij.real() << c_ij.imag() << "i";
                    }
                    else{
                        output << c_ij.real() << "+" << c_ij.imag() <<"i";
                    }
                }
                else{
                    output << c_ij.real();
                }
                if(j < Z_track[i].size() - 1){
                    output << ", "; 
                }
            }
            output << "}";
            if(i < Z_track.size()-1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << "std::bitset<N> D_product[Nop][D_maxsize] = {";
        for(int i = D_start; i < Z_track.size() ; i++){
            output << "{";
            for(int j = 0; j < Z_track[i].size(); j++){
                vector<int> z_ij = Zs[Z_track[i][j]];
                output << "std::bitset<N>(\"" << int_to_str(z_ij) << "\")";
                if(j < Z_track[i].size() - 1){
                    output << ", "; 
                }
            }
            output << "}";
            if(i < Z_track.size()-1){
                output << ", ";
            }
        }
        output << "};" << endl;

        output.close();
    }

    // ******************************************************************************//
    // ------------------------ Creating the cycles text file -----------------------//
    
    if(output_c.is_open()){
        for(int i = 0; i < Ps.size(); i++){
            output_c << "   " << boost::lexical_cast<std::string>(Ps[i]) << endl;
            output_c << "------------------------------" << endl;
            for(int j = 0; j < Z_track[i].size(); j++){
                complex<double> c_ij = coefficients[Z_track[i][j]];
                if(c_ij.imag() < 0){
                    output_c << "   " << c_ij.real() << c_ij.imag() << "i   ";
                }
                else{
                    output_c << "   " << c_ij.real() << "+" <<c_ij.imag() << "i   ";
                }
                vector<int> Zs_ij = Zs[Z_track[i][j]];
                for(int k=0; k < Zs_ij.size() ; k++){
                    output_c << Zs_ij[k] << "  ";
                }
                output_c << endl;
            }
            output_c << endl;
        }
        output_c << " The nullity is: " << endl;
        output_c << "   " << nullity << endl;
        output_c << endl;
        output_c << " The nullspace basis vectors (fundamental cycles) are: " << endl;
        for(int i=0; i < nullity; i++){
            output_c << endl;
            output_c << " The " << i+1 << " th fundamental cycles: " << endl;
            for(int j = 0; j < no_ps ; j++){
                output_c << "  " << nullspace[i][j] << endl;
            } 
        }
        output_c.close();
    }
    
    cout << "All done!";
    return 0;
}

