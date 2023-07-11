#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>

using namespace std;

struct BitsetComparator {
    bool operator()(const std::bitset<64>& lhs, const std::bitset<64>& rhs) const {
        return lhs.to_ulong() < rhs.to_ulong();
    }
};

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

int Int_sum(vector<int> vec){
    int sum = 0;
    for(int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}

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

void Cycle_minimize_0(vector<vector<int>>& null_eigs){
    int nullsize = null_eigs.size();

    vector<int> null_n(nullsize);
    for(int i = 0; i < nullsize; i++){
        null_n[i] = Int_sum(null_eigs[i]);
    }
    auto min_cyc = min_element(null_n.begin(), null_n.end());
    bool min_equal_three = true;

    vector<vector<int>> null_eigs_min , null_eigs_high;
    vector<int> null_eigs_min_ind , null_eigs_high_ind;
    
    for(int i = 0; i < nullsize; i++){
        if(null_n[i] == *min_cyc & min_equal_three){
            null_eigs_min.push_back(null_eigs[i]);
            null_eigs_min_ind.push_back(i);
            if(*min_cyc > 3){
                // this is to ensure that if no cycles are of length three or less, only one minimum is taken to be
                // in the null_eigs_min s. This is to make sure that all higher than three cycles get a chance to 
                // be minimized.
                min_equal_three = false;
            }
        }
        else{
            null_eigs_high.push_back(null_eigs[i]);
            null_eigs_high_ind.push_back(i);
        }
    }
    
    int high_size = null_eigs_high_ind.size();
    for(int k = 0; k < high_size; k++){
        vector<int> high_eig_k = null_eigs[null_eigs_high_ind[k]];
        int null_k = Int_sum(high_eig_k);
        for(int l = 0; l < null_eigs_min_ind.size(); l++){
            vector<int> high_to_min = GF2_add(high_eig_k , null_eigs[null_eigs_min_ind[l]]);
            int low_k = Int_sum(high_to_min);
            if(low_k <= 3){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                null_eigs_min_ind.push_back(null_eigs_high_ind[k]);
                std::sort(null_eigs_min_ind.begin() , null_eigs_min_ind.end());
                break;
            }
            else if(low_k < null_k){
                null_eigs[null_eigs_high_ind[k]] = high_to_min;
                high_eig_k = high_to_min;
            }
        }
    }
}

void Cycle_minimize(vector<vector<int>>& null_eigs){
    int nullsize = null_eigs.size();

    vector<int> null_n(nullsize);
    for(int i = 0; i < nullsize; i++){
        null_n[i] = Int_sum(null_eigs[i]);
    }
    vector<int> null_eigs_min_ind , null_eigs_high_ind;
    
    for(int i = 0; i < nullsize; i++){
        if(null_n[i] > 3){
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
                if(low_k == 3){
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
}

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

int main(){
    // Learning bitsets
    int no_qubit =2;
    bitset<64> bit0;
    //cout << "first bitset is: " << bit0;

    // sorting a vector of bitsets

    vector<vector<int>> vec;
    vector<int> vec1 = {1,1,0,0,0,1,0} , vec2={0,1,1,1,0,0,1} , vec3 = {0,1,1,0,1,1,0} , vec4 = {1,0,0,0,1,1,1};
    vec.push_back(vec1);
    vec.push_back(vec2);
    vec.push_back(vec3);
    vec.push_back(vec4);
    vector<vector<int>> vec_cop = vec;
    Cycle_minimize(vec);
    Cycle_minimize_0(vec_cop);
    
    cout << "previous version: " << endl;
    Print_matrix(vec_cop);
    cout << "current version: " << endl;
    Print_matrix(vec);

    Cycle_minimize(vec);
    cout << "current version after twice: " << endl;
    Print_matrix(vec);

    return 0;
}