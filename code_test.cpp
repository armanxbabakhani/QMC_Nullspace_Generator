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

int main(){
    // Learning bitsets
    int no_qubit =2;
    bitset<64> bit0;
    bitset<64> bit1;
    bitset<64> bit2;
    //cout << "first bitset is: " << bit0;

    // sorting a vector of bitsets
    vector<bitset<64>> bitsetvec;
    bitsetvec.push_back(bit0);
    bit1.set(4, true);
    bitsetvec.push_back(bit1);
    bit2.set(2, true);
    bitsetvec.push_back(bit2);
    bit1.set(12, true);
    bitsetvec.push_back(bit1);
    bit2.set(57, true);
    bitsetvec.push_back(bit2);
    bitsetvec.push_back(bit0);
    bit2.set(62, true);
    bitsetvec.push_back(bit2);
    bitsetvec.push_back(bit0);
    vector<int> indices;
    for(int i = 0; i < bitsetvec.size(); i++){
        indices.push_back(i);
    }
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return BitsetComparator()(bitsetvec[a], bitsetvec[b]);
    });

    no_qubit++;
    for(int i=0; i<indices.size();i++){
        cout << bitsetvec[indices[i]].to_ullong() << endl;
    }
    cout << endl;
    for(int i=0; i < bitsetvec.size(); i++){
        cout << bitsetvec[indices[i]] << endl;
    }

    return 0;
}