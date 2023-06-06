#define N        6
#define Nop      11
#define Ncycles  5

std::bitset<N> P_matrix[Nop] = {std::bitset<N>("000011"), std::bitset<N>("000101"), std::bitset<N>("000110"), std::bitset<N>("001010"), std::bitset<N>("001100"), std::bitset<N>("010010"), std::bitset<N>("010100"), std::bitset<N>("011000"), std::bitset<N>("101000"), std::bitset<N>("110000"), std::bitset<N>("")};
std::bitset<Nop> cycles[Ncycles] = {std::bitset<Nop>("11100000000"), std::bitset<Nop>("00111000000"), std::bitset<Nop>("00100110000"), std::bitset<Nop>("00010101000"), std::bitset<Nop>("00000001110")};

const int D0_size = 12;
double D0_coeff[D0_size] = {1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, -1};
std::bitset<N> D0_product[D0_size] = {std::bitset<N>("11"), std::bitset<N>("101"), std::bitset<N>("110"), std::bitset<N>("1010"), std::bitset<N>("1100"), std::bitset<N>("10100"), std::bitset<N>("10010"), std::bitset<N>("11000"), std::bitset<N>("101000"), std::bitset<N>("110000"), std::bitset<N>("100"), std::bitset<N>("10")};

const int D_maxsize = 4 ;
int D_size[Nop] = {4, 4, 2, 4, 4, 1, 1, 2, 2, 2, 1};
double D_coeff[Nop][D_maxsize] = {{0.707107, -0.707107, 0.707107, -0.707107}, {-0.707107, 0.707107, 0.707107, -0.707107}, {-0.5, 0.5}, {0.707107, -0.707107, 0.707107, -1.41421}, {-0.707107, 0.707107, 0.707107, -1.41421}, {0.707107}, {0.707107}, {1, -1}, {1, -1}, {1, -1}};
std::bitset<N> D_product[Nop][D_maxsize] = {{std::bitset<N>("0"), std::bitset<N>("11"), std::bitset<N>("100"), std::bitset<N>("111")}, {std::bitset<N>("10"), std::bitset<N>("111"), std::bitset<N>("0"), std::bitset<N>("101")}, {std::bitset<N>("10000"), std::bitset<N>("10110")}, {std::bitset<N>("0"), std::bitset<N>("1010"), std::bitset<N>("100"), std::bitset<N>("1110")}, {std::bitset<N>("10"), std::bitset<N>("1110"), std::bitset<N>("0"), std::bitset<N>("1100")}, {std::bitset<N>("100")}, {std::bitset<N>("0")}, {std::bitset<N>("0"), std::bitset<N>("11000")}, {std::bitset<N>("0"), std::bitset<N>("101000")}, {std::bitset<N>("0"), std::bitset<N>("110000")}};
