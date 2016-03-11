#include "converter.h"

char Converter::i2v [] = {'A','E','I','O','U'};
char Converter::i2k [] = {'B','C','D','F','G','H','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'};

void Converter::get_file_name(char buffer[], size_t idx){
	size_t tmp = idx;
	buffer[0]  = i2k[tmp%20];
	tmp       /= 20;
	buffer[1]  = i2v[tmp%5];
	tmp       /= 5;
	buffer[2]  = i2k[tmp%20];
	tmp       /= 20;
	buffer[3]  = i2k[tmp%20];
	tmp       /= 20;
	buffer[4]  = i2v[tmp%5];
	tmp       /= 5;
	buffer[5]  = i2k[tmp%20];
	tmp       /= 20;
	buffer[6]  = i2v[tmp%5];
	tmp       /= 5;
	buffer[7]  = i2k[tmp%20];
	tmp       /= 20;
	buffer[8]  = i2v[tmp%5];
	tmp       /= 5;
}
//                          A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P,  Q,  R,  S,  T,  U,  V,  W, X,  Y, Z
int Converter::char2i [] = {0, 0, 1, 2, 1, 3, 4, 5, 2, 0, 6, 7, 8, 9, 3, 10, 11, 12, 13, 14, 4, 15, 16, 17, 18, 19};

int Converter::get_idx_for_file(const char buffer[]){
    int ret  = char2i[buffer[0] - 65];
    ret += char2i[buffer[1] - 65] * 20;
    ret += char2i[buffer[2] - 65] * 100;
    ret += char2i[buffer[3] - 65] * 2000;
    ret += char2i[buffer[4] - 65] * 40000;
    ret += char2i[buffer[5] - 65] * 200000;
    ret += char2i[buffer[6] - 65] * 4000000;
    ret += char2i[buffer[7] - 65] * 20000000;
    ret += char2i[buffer[8] - 65] * 400000000;

    return ret;
}
