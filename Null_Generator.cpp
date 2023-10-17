#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>

using namespace std;

int no_qubit = 0; // number of qubits is a global variable!

// *****************************    Functions ******************************************* //
template<typename T>
void Print_matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    int n = matrix[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
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

int Int_sum(vector<int> vec){
    int sum = 0;
    for(int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}

// This function minimizes the size of the cycles (nullspace basis vectors with least 1s):
int Cycle_minimize(vector<vector<int>>& null_eigs){
    int nullsize = null_eigs.size() , number_minimized = 0;

    vector<int> null_n(nullsize);
    for(int i = 0; i < nullsize; i++){
        null_n[i] = Int_sum(null_eigs[i]);
    }
    vector<int> null_eigs_high_ind;
    auto min_cyc = min_element(null_n.begin(), null_n.end());
    
    for(int i = 0; i < nullsize; i++){
        if(null_n[i] > *min_cyc){
            null_eigs_high_ind.push_back(i);
        }
    }
    
    int high_size = null_eigs_high_ind.size();
    for(int k = 0; k < high_size; k++){
        vector<int> high_eig_k = null_eigs[null_eigs_high_ind[k]];
        int null_k = Int_sum(high_eig_k);
        for(int l = 0; l < nullsize; l++){
            if( l != null_eigs_high_ind[k]){
                vector<int> high_to_min = GF2_add(high_eig_k , null_eigs[l]);
                int low_k = Int_sum(high_to_min);
                if(low_k == *min_cyc){
                    number_minimized++;
                    null_eigs[null_eigs_high_ind[k]] = high_to_min;
                    break;
                }
                else if(low_k < null_k){
                    null_eigs[null_eigs_high_ind[k]] = high_to_min;
                    high_eig_k = high_to_min;
                }
            }
        }
    }
    return high_size - number_minimized;
}

pair<int, int> Cycle_minimize_0(vector<vector<int>>& null_eigs , int min_cyc_length){
    int nullsize = null_eigs.size();
    pair<int, int> cycdata;
    vector<int> null_n(nullsize);
    for(int i = 0; i < nullsize; i++){
        null_n[i] = Int_sum(null_eigs[i]);
    }
    auto min_cyc = min_element(null_n.begin(), null_n.end());
    bool min_cycle_satisfied = true;

    vector<vector<int>> null_eigs_min , null_eigs_high;
    vector<int> null_eigs_min_ind , null_eigs_high_ind;
    
    for(int i = 0; i < nullsize; i++){
        if(null_n[i] == *min_cyc & min_cycle_satisfied){
            null_eigs_min.push_back(null_eigs[i]);
            null_eigs_min_ind.push_back(i);
            if(*min_cyc > min_cyc_length){
                // this is to ensure that if no cycles are of length three or less, only one minimum is taken to be
                // in the null_eigs_min s. This is to make sure that all higher than three cycles get a chance to 
                // be minimized.
                min_cycle_satisfied = false;
            }
        }
        else{
            null_eigs_high.push_back(null_eigs[i]);
            null_eigs_high_ind.push_back(i);
        }
    }
    
    int high_size = null_eigs_high_ind.size() , number_of_minimized=0;
    for(int k = 0; k < high_size; k++){
        vector<int> high_eig_k = null_eigs[null_eigs_high_ind[k]];
        int null_k = Int_sum(high_eig_k);
        for(int l = 0; l < null_eigs_min_ind.size(); l++){
            vector<int> high_to_min = GF2_add(high_eig_k , null_eigs[null_eigs_min_ind[l]]);
            int low_k = Int_sum(high_to_min);
            if(low_k <= min_cyc_length){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                null_eigs_min_ind.push_back(null_eigs_high_ind[k]);
                std::sort(null_eigs_min_ind.begin() , null_eigs_min_ind.end());
                number_of_minimized++;
                break;
            }
            else if(low_k < null_k){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                high_eig_k = high_to_min;
            }
        }
    }
    cycdata.first = high_size - number_of_minimized;
    cycdata.second = Int_sum(null_eigs[null_eigs_min_ind[0]]);
    return cycdata;
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
vector<pair<complex<double>, vector<int>>> Data_extract(const string& fileName){
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
            int num = std::stoi(token);
            integers.push_back(num);
        }
        linedata.second = integers;

        data.push_back(linedata);
        }
        inputFile.close();
    return data;
}

// Finding bitset in a vector of bitsets:
// The max bitset will be 64 (the maximum number of qubits), and we will cut off accordingly
// ..... when the number of qubits of the system is less than this number.
pair<bool , int> Bit_is_in_set(bitset<5000> bitstring, vector<bitset<5000>> bitsetVector) {
    bool found = false;
    int found_indx = 0, ind_count = 0;
    pair<bool, int> output;
    for (const auto& bitset : bitsetVector) {
            if (bitset == bitstring) {
                found = true;
                found_indx = ind_count;
                break;
            }
            ind_count++;
        }
    output.first = found;
    output.second = found_indx;
    return output;
}

// This function downsizes the vector of bitsets from bitset<64> to bitset<no_qubit> to avoid 
// ..... having redundant zeros.
vector<vector<bool>> Downsize_bitset(vector<bitset<5000>> bitsetVector){
    extern int no_qubit;
    vector<vector<bool>> bitset_down;
    for(const auto& bitset : bitsetVector){
        vector<bool> bitset_vector;
        for(int i = 0; i < no_qubit; i++){
            bitset_vector.push_back(bitset[i]);
        }
        reverse(bitset_vector.begin() , bitset_vector.end())
;        bitset_down.push_back(bitset_vector);
    }
    return bitset_down;
}

vector<vector<int>> Bit_to_intvec(vector<vector<bool>> bitsetVec){
    extern int no_qubit;
    vector<vector<int>> bit_int;
    for(auto const& bitset : bitsetVec){
        vector<int> bit_int_i;
        for(int i=0; i < no_qubit; i++){
            bit_int_i.push_back(bitset[i]);
        }
        reverse(bit_int_i.begin() , bit_int_i.end());
        bit_int.push_back(bit_int_i);
    }

    return bit_int;
}

// ------------- Functions to compute the Zs and Ps and coefficients ----------
struct BitsetComparator {
    bool operator()(const std::bitset<5000>& lhs, const std::bitset<5000>& rhs) const {
        return lhs.to_ulong() < rhs.to_ulong();
    }
}; // This struct is for bitset comparison (in order to sort the bitsets)

typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
struct PZdata {
    vector<bitset<5000>> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track;
};

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    PZdata PZ_data;
    int l = data.size(),z_count = 0;
    extern int no_qubit;
    vector<bitset<5000>> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track; //This vector maps the Zs to Ps it is a many to one mapping!

    for (int i = 0; i<l; i++){
        complex<double> coeff_i = data[i].first;
        vector<int> zs_i; // For every line zs extracts the qubits on which a pauli Z acts! 
        vector<int> data_i = data[i].second; // Extracts the array of qubits and paulis for every line of input!
        bitset<5000> bit_num; // This variable keeps track of the index of the permutation matrix we get for each line of data!

        for (size_t j = 0; j < data_i.size() / 2; j++) {
            // Format of the input file: The 1st, 3rd, 5th, ... indicate the qubits
            int qubit = data_i[2 * j];
            // Format of the input file: The 2nd, 4th, 6th, ... indicate the paulis
            int pauli_j = data_i[2 * j + 1];

            if (qubit > no_qubit)
                no_qubit = qubit;
            if (pauli_j == 1) {
                bit_num.set(qubit-1, true);
                // Add a 1 in the num-th position from the right of the bit string 
            } else if (pauli_j == 2) {
                //cpp_int perm_term = 1 << (qubit - 1);
                //num += perm_term;
                bit_num.set(qubit-1, true);
                coeff_i *= complex<double>(0, 1);
                zs_i.push_back(qubit);
            } else if (pauli_j == 3) {
                zs_i.push_back(qubit);
            }
        }
        coeffs.push_back(coeff_i);
        Zs.push_back(zs_i); // If zs is empty then no non-trivial diagonal components! (All identity operators)
            
        // Look for num in the previous list of Ps (permutations)
        pair<bool , int> bit_in_set = Bit_is_in_set(bit_num , Ps);
        if(!bit_in_set.first){ 
            // num was not found in Ps, thus a new permutation matrix!
            Ps.push_back(bit_num);
            // We will add i-th element of coeffs and Zs to be associated with the current Ps!
            Z_track.push_back(vector<int> {i});
        }else {
            // num was found in Ps, and P_index will be the index of Ps that matched num.
            int P_index = bit_in_set.second;

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
    vector<bitset<5000>> Ps_kept;
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
    vector<int> indices;
    for(int i = 0; i < Ps.size(); i++){
        indices.push_back(i);
    }
    sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return BitsetComparator()(Ps_kept[a], Ps_kept[b]);
    });


    vector<bitset<5000>> Ps_sorted;
    ZVecs Z_track_sorted;
    for(int i = 0; i < Ps_kept.size(); i++){
        Ps_sorted.push_back(Ps_kept[indices[i]]);
        Z_track_sorted.push_back(Z_track_kept[indices[i]]);
    }

    PZ_data.coeffs = coeffs; //Coeffs and Zs are kept as the original data
    PZ_data.Ps = Ps_sorted;
    PZ_data.Zs = Zs;
    PZ_data.Z_track = Z_track_sorted;

    return PZ_data;
}

// ************************************************************************************************************** //
// -------------------------------------------------------------------------------------------------------------- //


//using namespace boost;
int main(int argc , char* argv[]){
    string fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> data = Data_extract(fileName);

    // Unpacking the data from the input file "fileName"
    PZdata PZ_data = PZcomp(data);
    vector<bitset<5000>> Ps_bit = PZ_data.Ps;
    vector<vector<bool>> Ps = Downsize_bitset(Ps_bit);
    vector<complex<double>> coefficients = PZ_data.coeffs;
    vector<vector<int>> Z_track = PZ_data.Z_track;
    vector<vector<int>> Zs = PZ_data.Zs;
    vector<string> Zs_string;
    bool D0_exists = false;
    //int no_qubit = PZ_data.no_qubit;

    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs.size();i++){
        Zs_string.push_back(int_to_str(Zs[i]));
    }

    // Convert Ps indices into vector of ints
    vector<vector<bool>> Ps_nontrivial = Ps;
    if(Ps_bit[0].to_ullong() < 1e-6){
        Ps_nontrivial.erase(Ps_nontrivial.begin());
        D0_exists = true;
    }

    // Minimizing the size of the fundamental cycless
    vector<vector<int>> Ps_binary = Bit_to_intvec(Ps_nontrivial);
    vector<vector<int>> nullspace = Null2(Ps_binary);
    // Create a for loop to minimize the cycles!
    int high_size = 1000;
    int count = 0 , min_cyc_len = 3;
    pair<int,int> cycdata;
    cycdata = Cycle_minimize_0(nullspace , min_cyc_len);
    high_size = cycdata.first;
    min_cyc_len = cycdata.second;
    while(high_size > 0){
        cycdata = Cycle_minimize_0(nullspace , min_cyc_len);
        high_size = cycdata.first;
        min_cyc_len = cycdata.second;
        count++;
    }
    cout << "the number of times minimization was run: " << count << endl;
    int no_ps = Ps_binary.size();
    int nullity = nullspace.size();

    cout << "Calculations done!" << endl;
    cout << "Making the output files ... " << endl;

    string output_h = fileName.substr(0, fileName.find_last_of(".")) + ".h";
    ofstream output(output_h);

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
            for(int j=0; j < Ps_nontrivial[i].size(); j++){
                output << Ps_nontrivial[i][j];
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
        if(D0_exists)
            output << "const int D0_size = " << Z_track[0].size() << ";" << endl;  
        else
            output << "const int D0_size = 0;" << endl;
            
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
    
    cout << "All done!";
    return 0;
}

