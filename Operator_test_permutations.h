#define N        6  // This defines the number of qubits 
#define Nop      10  // This defines the number of operators in the Hamiltonian
#define Ncycles  2  // This defines the number of non-trivial permutations of the measurement operator

// P_matrix is the set of permutation matrices in the Hamiltonian
std::bitset<N> P_matrix[Nop] = {std::bitset<N>("000011"), std::bitset<N>("000101"), std::bitset<N>("000110"), std::bitset<N>("001010"), std::bitset<N>("001100"), std::bitset<N>("010010"), std::bitset<N>("010100"), std::bitset<N>("011000"), std::bitset<N>("101000"), std::bitset<N>("110000")};
// Meas_Op[Ncycles] is a set of bitsets representing each permutation of the operator
//       as a product of permutations available in the Hamiltonian.
std::bitset<Nop> Meas_Op[Ncycles] = {std::bitset<Nop>("0001000000"), std::bitset<Nop>("0000100000")};

const int D0_size = 1;
double D0_coeff[D0_size] = {3.91};
std::bitset<N> D0_product[D0_size] = {std::bitset<N>("101")};

const int D_maxsize = 2 ;
int D_size[Ncycles] = {1, 2};
double D_coeff[Ncycles][D_maxsize] = {{-1.86i}, {-2.5+1.5i, -1.9}};
std::bitset<N> D_product[Ncycles][D_maxsize] = {{std::bitset<N>("10")}, {std::bitset<N>("1000"), std::bitset<N>("1100")}};
