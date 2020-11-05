#pragma once
#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "logregtools.h"

using namespace std;
using namespace seal;

int Fixed_Hessian_Chebyshev(bool ringdim);
int Fixed_Hessian_Compact(bool ringdim);
int Fixed_Hessian_Taylor();
int Fixed_Hessian_IDASH();
int Nesterov_GD(bool ringbool);
int GD(bool ringdim);
int Plaintext_LR(string filename, int iternum);
int Plaintext_LR_NAG(string filename, int iternum);



