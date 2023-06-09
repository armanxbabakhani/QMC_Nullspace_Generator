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

// ***************************** Functions ******************************************* //
template<typename T>
void printMatrix(const vector<vector<T>>& matrix) {
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

    // Remove consecutive duplicates (taking only unique marked rows!)
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
// This Eig_minimize is specifically designed for operator permutations, and functions differently
//      than the one in Null_Generator, which is used to find closed cycles of minimum length 3.
void Eig_minimize(vector<vector<int>>& nullEigs){
    int nullsize = nullEigs.size() , noops = nullEigs[0].size(); //noops is the number of operators

    vector<int> nulln(nullsize);
    for(int i = 0; i < nullsize; i++){
        nulln[i] = Int_sum(nullEigs[i]);
    }
    auto mincyc = min_element(nulln.begin(), nulln.end());

    vector<int> nullEigsMinind , nullEigsHighind, nullEigsOnesind; // the indices of eigenvectors ending with 1
    
    for(int i = 0; i < nullsize; i++){
        if(nulln[i] == *mincyc){
            nullEigsMinind.push_back(i);
        }
        else if(nullEigs[i][noops-1] == 1){
            nullEigsOnesind.push_back(i);
        }
        else{
            nullEigsHighind.push_back(i);
        }
    }
    int highsize = nullEigsHighind.size();
    int minsize = nullEigsMinind.size();
    
    for(int i=0; i < minsize; i++){
        int vecisum = Int_sum(nullEigs[nullEigsMinind[i]]);
        vector<int> bestVec = nullEigs[nullEigsMinind[i]]; 
        for(int j=0; j < minsize; j++){
            if(j!= i){
                vector<int> tryVecij = GF2_add(nullEigs[nullEigsMinind[i]] , nullEigs[nullEigsMinind[j]]);
                int tryvecijsum = Int_sum(tryVecij);
                if(tryvecijsum < vecisum){
                    vecisum = tryvecijsum;
                    bestVec = tryVecij;
                }
            }
        }
        nullEigs[nullEigsMinind[i]] = bestVec;
    }

    for(int k = 0; k < nullEigsOnesind.size(); k++){
        int oneskind = nullEigsOnesind[k];
        vector<int> highEigk = nullEigs[oneskind];
        vector<int> bestEigk = highEigk;
        int bestnullk = Int_sum(highEigk);
        for(int l = 0; l < nullsize; l++){
            if(l != oneskind){
                vector<int> highTomin = GF2_add(highEigk , nullEigs[l]);
                int lowk = Int_sum(highTomin);
                if(lowk < 3){
                    bestEigk = highTomin;
                    bestnullk = lowk;
                    break;
                }
                else if(lowk < bestnullk){
                    bestEigk = highTomin;
                    bestnullk = lowk;
                }
            }
        }
        nullEigs[oneskind] = bestEigk;
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
pair<bool , int> bit_is_in_set(bitset<64> bitstring, vector<bitset<64>> bitsetVector) {
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
vector<vector<bool>> downsize_bitset(vector<bitset<64>> bitsetVector){
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

vector<vector<int>> bit_to_intvec(vector<vector<bool>> bitsetVec){
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
    bool operator()(const std::bitset<64>& lhs, const std::bitset<64>& rhs) const {
        return lhs.to_ulong() < rhs.to_ulong();
    }
}; // This struct is for bitset comparison (in order to sort the bitsets)

typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
struct PZdata {
    vector<bitset<64>> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track;
};

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    PZdata PZ_data;
    int l = data.size(),z_count = 0;
    extern int no_qubit;
    vector<bitset<64>> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track; //This vector maps the Zs to Ps: it is a many to one mapping!

    for (int i = 0; i<l; i++){
        complex<double> coeff_i = data[i].first;
        vector<int> zs_i; // For every line zs extracts the qubits on which a pauli Z acts! 
        vector<int> data_i = data[i].second; // Extracts the array of qubits and paulis for every line of input!
        bitset<64> bit_num; // This variable keeps track of the index of the permutation matrix we get for each line of data!

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
        pair<bool , int> bit_in_set = bit_is_in_set(bit_num , Ps);
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
    vector<bitset<64>> Ps_kept;
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


    vector<bitset<64>> Ps_sorted;
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
    string Ham_fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> Ham_data = data_extract(Ham_fileName);

    string Op_fileName(argv[2]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> Op_data = data_extract(Op_fileName);

    // Unpacking the Hamiltonian data from the input file 
    PZdata PZ_data_Ham = PZcomp(Ham_data);
    vector<bitset<64>> Ps_bit_Ham = PZ_data_Ham.Ps;
    vector<vector<bool>> Ps_Ham = downsize_bitset(Ps_bit_Ham);
    vector<complex<double>> coefficients_Ham = PZ_data_Ham.coeffs;
    vector<vector<int>> Z_track_Ham = PZ_data_Ham.Z_track;
    vector<vector<int>> Zs_Ham = PZ_data_Ham.Zs;
    vector<string> Zs_string_Ham;

    // Unpacking the Operator data from the input file 
    PZdata PZ_data_Op = PZcomp(Op_data);
    vector<bitset<64>> Ps_bit_Op = PZ_data_Op.Ps;
    vector<vector<bool>> Ps_Op = downsize_bitset(Ps_bit_Op);
    vector<complex<double>> coefficients_Op = PZ_data_Op.coeffs;
    vector<vector<int>> Z_track_Op = PZ_data_Op.Z_track;
    vector<vector<int>> Zs_Op = PZ_data_Op.Zs;
    vector<string> Zs_string_Op;
    bool D0_exists = false, P0_exists = false;

    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs_Ham.size();i++){
        Zs_string_Ham.push_back(int_to_str(Zs_Ham[i]));
    }

    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs_Op.size();i++){
        Zs_string_Op.push_back(int_to_str(Zs_Op[i]));
    }

    // Convert Ps indices into vector of bools
    vector<vector<bool>> Ps_Ham_nontrivial = Ps_Ham;
    if(Ps_bit_Ham[0].to_ullong() < 1e-6){
        Ps_Ham_nontrivial.erase(Ps_Ham_nontrivial.begin());
        D0_exists = true;
    }

    // Convert Ps indices into vector of bools
    vector<vector<bool>> Ps_Op_nontrivial = Ps_Op;
    if(Ps_bit_Op[0].to_ullong() < 1e-6){
        Ps_Op_nontrivial.erase(Ps_Op_nontrivial.begin());
        P0_exists = true;
    }

    vector<vector<int>> Ps_binary_Ham = bit_to_intvec(Ps_Ham_nontrivial);
    vector<vector<int>> Ps_binary_Op = bit_to_intvec(Ps_Op_nontrivial);

    // Combine the permutations of Hamiltonians with each permutation from the operator and compute the nullspace
    int no_ps = Ps_Ham_nontrivial.size() , no_ops = Ps_Op_nontrivial.size();
    vector<vector<int>> Op_to_Ham;
    vector<int> zero_vector;
    for(int i = 0; i < no_ps; i++){
        zero_vector.push_back(0);
    }

    for(int k=0; k < Ps_binary_Op.size(); k++){
        vector<vector<int>> Ps_bin_total = Ps_binary_Ham, nullspace_k; 
        bool permutation_found = false;
        Ps_bin_total.push_back(Ps_binary_Op[k]);
        cout << "k = " << k << endl;
        cout << "The permutation set is: " << endl;
        printMatrix(Ps_bin_total);
        cout << "The nullspace is: " << endl;
        nullspace_k = Null2(Ps_bin_total);
        Eig_minimize(nullspace_k);
        int min_ops=100000 , min_index;
        for(int i = 0; i< nullspace_k.size(); i++){
            int sum_ops = 0;
            if(nullspace_k[i][no_ps] == 1){
                permutation_found = true;
                sum_ops = Int_sum(nullspace_k[i]);
                if(sum_ops < min_ops){
                    min_ops = sum_ops;
                    min_index = i;
                }
            }
        }
        cout << "min_index is " << min_index << endl; 
        if(permutation_found){
            Op_to_Ham.push_back(nullspace_k[min_index]);
        }
        else{
            Op_to_Ham.push_back(zero_vector);
        }
        cout << endl;
    }

    cout << "Calculations done!" << endl;
    cout << "Making the output files ... " << endl;

    string output_operator = Op_fileName.substr(0, Op_fileName.find_last_of(".")) + "_permutations.h";
    ofstream output(output_operator);


    // ------- Operator Measurement permutations in terms of Hamiltonian permutations ----- //
    // Printing only the relevant permutations for the measurement operator (any operator that can be
    // written as a mod 2 addition of the permutations in the Hamiltonian.)
    vector<int> rel_perms;
    bool non_triv_offdiags_exists = false;
    for(int i = 0; i < no_ops; i++){
        if(Int_sum(Op_to_Ham[i]) > 0){
            rel_perms.push_back(i);
        }
    }
    // --------------------------------------------------------
    // Delete the Z_tracks of the trivial permutations:
    vector<vector<int>> Z_track_Op_kept;
    int counter = 0;
    for(int i = 0 ; i < no_ops; i++){
        if(rel_perms[counter] == i){
            Z_track_Op_kept.push_back(Z_track_Op[i + int(P0_exists)]);
            counter++;
        }
    }
    no_ops = rel_perms.size();
    // ********************************************************************** //
    // -------------------- Creating the the .h file ------------------------ //
    if(output.is_open()){
        output << "#define N        " << no_qubit << "  // This defines the number of qubits " << endl;
        output << "#define Nop      " << no_ps << "  // This defines the number of operators in the Hamiltonian" << endl; // Number of permutations of the measurement operator!
        output << "#define Ncycles  " << no_ops << "  // This defines the number of non-trivial permutations of the measurement operator" << endl;
        output << endl;
        
        // ---------------- Permutation matrices and cycles --------------------- //
        // The permutation matrices of the Hamiltonian
        output << "// P_matrix is the set of permutation matrices in the Hamiltonian" <<endl;
        output << "std::bitset<N> P_matrix[Nop] = {";
        for(int i = 0; i < no_ps; i++){
            output << "std::bitset<N>(\"";
            for(int j=0; j < Ps_Ham_nontrivial[i].size(); j++){
                output << Ps_Ham_nontrivial[i][j];
            }
            output << "\")";
            if(i < no_ps - 1){
                output << ", ";
            }
        }
        output << "};" << endl;

        if(rel_perms.size() > 0){
            non_triv_offdiags_exists = true;
            output << "// Meas_Op[Ncycles] is a set of bitsets representing each permutation of the operator" <<endl;
            output << "//       as a product of permutations available in the Hamiltonian." <<endl;
            output << "std::bitset<Nop> Meas_Op[Ncycles] = {";
            for(int i = 0; i < no_ops; i++){
                output << "std::bitset<Nop>(\"";
                for(int j=0; j < no_ps; j++){
                    output << Op_to_Ham[rel_perms[i]][j];
                }
                output << "\")";
                if(i < no_ops-1){
                    output << ", ";
                }
            }
            output << "};" << endl;
            output << endl; 
        }else{
            output << "No relevant permutation operators. The expectation value of this operator is zero." << endl;
        }
        

        // ---------- Permutation Matrix Representation of the measurement operator -------- //
        // --------------------------------- Diagonal terms -------------------------------- //
        if(P0_exists){
            output << "const int D0_size = " << Z_track_Op[0].size() << ";" << endl;
            output << "double D0_coeff[D0_size] = {";
            for(int i=0;i<Z_track_Op[0].size();i++){
                complex<double> c_ij = coefficients_Op[Z_track_Op[0][i]];
                if(abs(c_ij.real()) > 1e-7){
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
                }else{
                    if(abs(c_ij.imag()) > 1e-7 ){
                        if(c_ij.imag() < 0){
                        output << c_ij.imag() << "i";
                        }
                        else{
                            output  << c_ij.imag() << "i";
                        }
                    }
                }
                if(i < Z_track_Op[0].size() - 1){
                    output << ", "; 
                }
            }
            output << "};" << endl;
            output << "std::bitset<N> D0_product[D0_size] = {";
            for(int i = 0; i < Z_track_Op[0].size(); i++){
                //vector<int> Zs_i = Zs_Op[Z_track_Op[0][i]];
                //output << "std::bitset<N>(\"" << int_to_str(Zs_i) << "\")";
                output << "std::bitset<N>(\"" << Zs_string_Op[Z_track_Op[0][i]] << "\")";
                if(i < Z_track_Op[0].size() - 1){
                    output << ", ";
                }
            }
            output << "};" << endl;
            output << endl; 
        }else{
            output << "const int D0_size = 0;" << endl;
            output << endl;
        }

        // ------------------------------------------------------------------------ //
        // --------------------------- Off-Diagonal terms ------------------------- //
        // Finding D_maxsize:
        if(non_triv_offdiags_exists){
            int D_max = 0;
            vector<int> D_size;
            for(int i = 0; i < no_ops; i++){
                int z_size_i = Z_track_Op_kept[i].size();
                D_size.push_back(z_size_i);
                if(z_size_i > D_max){
                    D_max = z_size_i;
                }
            }
            output << "const int D_maxsize = " << D_max << " ;" << endl;
            output << "int D_size[Ncycles] = {";
            for(int i=0;i<no_ops;i++){
                output << D_size[i];
                if(i < no_ops - 1){
                    output << ", ";
                }
            }

            output << "};" << endl;
            output << "double D_coeff[Ncycles][D_maxsize] = {";
            for(int i = 0; i < no_ops; i++){
                output << "{";
                for(int j = 0; j < Z_track_Op_kept[i].size(); j++){
                    complex<double> c_ij = coefficients_Op[Z_track_Op_kept[i][j]];
                    if(abs(c_ij.real()) > 1e-7){
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
                    }else{
                        if(abs(c_ij.imag()) > 1e-7 ){
                            if(c_ij.imag() < 0){
                            output << c_ij.imag() << "i";
                            }
                            else{
                                output  << c_ij.imag() << "i";
                            }
                        }
                    }
                    if(j < Z_track_Op_kept[i].size() - 1){
                        output << ", "; 
                    }
                }
                output << "}";
                if(i < no_ops-1){
                    output << ", ";
                }
            }
            output << "};" << endl;
            output << "std::bitset<N> D_product[Ncycles][D_maxsize] = {";
            for(int i = 0; i < Z_track_Op_kept.size() ; i++){
                output << "{";
                for(int j = 0; j < Z_track_Op_kept[i].size(); j++){
                    output << "std::bitset<N>(\"" << Zs_string_Op[Z_track_Op_kept[i][j]] << "\")";
                    if(j < Z_track_Op_kept[i].size() - 1){
                        output << ", "; 
                    }
                }
                output << "}";
                if(i < Z_track_Op_kept.size()-1){
                    output << ", ";
                }
            }
            output << "};" << endl;

            output.close();
        }
    }
    
    cout << "All done!";
    return 0;
}

