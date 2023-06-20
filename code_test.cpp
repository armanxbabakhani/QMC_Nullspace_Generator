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
    vector<int> vec1 = {0,1,1,1,1,0,0,0} , vec2={1,1,1,0,0,1,0,0} , vec3 = {1,1,1,0,1,1,0,1} , vec4 = {1,0,1,0,1,0,0,0};
    vec.push_back(vec1);
    vec.push_back(vec2);
    vec.push_back(vec3);
    vec.push_back(vec4);
    Eig_minimize(vec);
    Print_matrix(vec);

    return 0;
}