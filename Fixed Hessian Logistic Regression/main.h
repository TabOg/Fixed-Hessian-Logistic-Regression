#pragma once
#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "logregtools.h"

using namespace std;
using namespace seal;

int Fixed_Hessian_Compact(bool ringdim);
int Fixed_Hessian_Chebyshev(const unsigned int precision, const unsigned number_of_threads);

int Fixed_Hessian_Taylor();
int Fixed_Hessian_IDASH();
int Nesterov_GD_split(const unsigned int precision);
int Nesterov_GD(const unsigned int precision);

int GD(const unsigned int precision, const unsigned int number_of_threads);
int Plaintext_LR(string filename, int iternum);
int Plaintext_LR_NAG(string filename, int iternum);
int Plaintext_LR_lowdeg(string filename, int iternum);
