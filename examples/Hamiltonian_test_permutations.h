#define N        6  // This defines the number of qubits 
#define Nop      10  // This defines the number of operators in the Hamiltonian
#define Ncycles  10  // This defines the number of non-trivial permutations of the measurement operator

// P_matrix is the set of permutation matrices in the Hamiltonian
std::bitset<N> P_matrix[Nop] = {std::bitset<N>("000011"), std::bitset<N>("000101"), std::bitset<N>("000110"), std::bitset<N>("001010"), std::bitset<N>("001100"), std::bitset<N>("010010"), std::bitset<N>("010100"), std::bitset<N>("011000"), std::bitset<N>("101000"), std::bitset<N>("110000")};
// Meas_Op[Ncycles] is a set of bitsets representing each permutation of the operator
//       as a product of permutations available in the Hamiltonian.
std::bitset<Nop> Meas_Op[Ncycles] = {std::bitset<Nop>("1000000000"), std::bitset<Nop>("0100000000"), std::bitset<Nop>("0010000000"), std::bitset<Nop>("0001000000"), std::bitset<Nop>("0000100000"), std::bitset<Nop>("0000010000"), std::bitset<Nop>("0000001000"), std::bitset<Nop>("0000000100"), std::bitset<Nop>("0000000010"), std::bitset<Nop>("0000000001")};

const int D0_size = 12;
double D0_coeff[D0_size] = {1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, -1};
std::bitset<N> D0_product[D0_size] = {std::bitset<N>("11"), std::bitset<N>("101"), std::bitset<N>("110"), std::bitset<N>("1010"), std::bitset<N>("1100"), std::bitset<N>("10100"), std::bitset<N>("10010"), std::bitset<N>("11000"), std::bitset<N>("101000"), std::bitset<N>("110000"), std::bitset<N>("100"), std::bitset<N>("10")};

const int D_maxsize = 4 ;
int D_size[Ncycles] = {4, 4, 2, 4, 4, 1, 1, 2, 2, 2};
double D_coeff[Ncycles][D_maxsize] = {{0.707107, -0.707107, 0.707107, -0.707107}, {-0.707107, 0.707107, 0.707107, -0.707107}, {-0.5, 0.5}, {0.707107, -0.707107, 0.707107, -1.41421}, {-0.707107, 0.707107, 0.707107, -1.41421}, {0.707107}, {0.707107}, {1, -1}, {1, -1}, {1, -1}};
std::bitset<N> D_product[Ncycles][D_maxsize] = {{std::bitset<N>("0"), std::bitset<N>("11"), std::bitset<N>("100"), std::bitset<N>("111")}, {std::bitset<N>("10"), std::bitset<N>("111"), std::bitset<N>("0"), std::bitset<N>("101")}, {std::bitset<N>("10000"), std::bitset<N>("10110")}, {std::bitset<N>("0"), std::bitset<N>("1010"), std::bitset<N>("100"), std::bitset<N>("1110")}, {std::bitset<N>("10"), std::bitset<N>("1110"), std::bitset<N>("0"), std::bitset<N>("1100")}, {std::bitset<N>("100")}, {std::bitset<N>("0")}, {std::bitset<N>("0"), std::bitset<N>("11000")}, {std::bitset<N>("0"), std::bitset<N>("101000")}, {std::bitset<N>("0"), std::bitset<N>("110000")}};
