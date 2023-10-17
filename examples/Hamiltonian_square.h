#define N        9
#define Nop      12
#define Ncycles  4

std::bitset<N> P_matrix[Nop] = {std::bitset<N>("000000011"), std::bitset<N>("000000110"), std::bitset<N>("000001001"), std::bitset<N>("000010010"), std::bitset<N>("000011000"), std::bitset<N>("000100100"), std::bitset<N>("000110000"), std::bitset<N>("001001000"), std::bitset<N>("010010000"), std::bitset<N>("011000000"), std::bitset<N>("100100000"), std::bitset<N>("110000000")};
std::bitset<Nop> cycles[Ncycles] = {std::bitset<Nop>("101110000000"), std::bitset<Nop>("010101100000"), std::bitset<Nop>("000010011100"), std::bitset<Nop>("010101001011")};

const int D0_size = 1;
double D0_coeff[D0_size] = {1};
std::bitset<N> D0_product[D0_size] = {std::bitset<N>("0")};

const int D_maxsize = 1 ;
int D_size[Nop] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
double D_coeff[Nop][D_maxsize] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
std::bitset<N> D_product[Nop][D_maxsize] = {{std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}};
