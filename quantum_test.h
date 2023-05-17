#define N        9
#define Nop      25
#define Ncycles  16

std::bitset<N> P_matrix[Nop] = {std::bitset<N>("1"), std::bitset<N>("10"), std::bitset<N>("11"), std::bitset<N>("100"), std::bitset<N>("110"), std::bitset<N>("1000"), std::bitset<N>("1001"), std::bitset<N>("1010"), std::bitset<N>("10000"), std::bitset<N>("10010"), std::bitset<N>("10100"), std::bitset<N>("11000"), std::bitset<N>("100000"), std::bitset<N>("100100"), std::bitset<N>("110000"), std::bitset<N>("1000000"), std::bitset<N>("1001000"), std::bitset<N>("1010000"), std::bitset<N>("10000000"), std::bitset<N>("10010000"), std::bitset<N>("10100000"), std::bitset<N>("11000000"), std::bitset<N>("100000000"), std::bitset<N>("100100000"), std::bitset<N>("110000000")};
std::bitset<Nop> cycles[Ncycles] = {std::bitset<Nop>("1110000000000000000000000"), std::bitset<Nop>("0101100000000000000000000"), std::bitset<Nop>("1000011000000000000000000"), std::bitset<Nop>("0100010100000000000000000"), std::bitset<Nop>("0100000011000000000000000"), std::bitset<Nop>("0001000010100000000000000"), std::bitset<Nop>("0000010010010000000000000"), std::bitset<Nop>("0001000000001100000000000"), std::bitset<Nop>("0000000010001010000000000"), std::bitset<Nop>("0000010000000001100000000"), std::bitset<Nop>("0000000010000001010000000"), std::bitset<Nop>("0000000010000000001100000"), std::bitset<Nop>("0000000000001000001010000"), std::bitset<Nop>("0000000000000001001001000"), std::bitset<Nop>("0000000000001000000000110"), std::bitset<Nop>("0000000000000000001000101")};

const int D0_size = 16;
double D0_coeff[D0_size] = {-1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1};
std::bitset<N> D0_product[D0_size] = {std::bitset<N>("11"), std::bitset<N>("1001"), std::bitset<N>("110"), std::bitset<N>("1010"), std::bitset<N>("10010"), std::bitset<N>("10100"), std::bitset<N>("100100"), std::bitset<N>("11000"), std::bitset<N>("1001000"), std::bitset<N>("110000"), std::bitset<N>("1010000"), std::bitset<N>("10010000"), std::bitset<N>("10100000"), std::bitset<N>("100100000"), std::bitset<N>("11000000"), std::bitset<N>("110000000");

const int D_maxsize = 1 ;
int D_size[Nop] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
double D_coeff[Nop][D_maxsize] = {{1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}, {1+0}};
std::bitset<N> D_product[Nop][D_maxsize] = {{std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}, {std::bitset<N>("0")}};