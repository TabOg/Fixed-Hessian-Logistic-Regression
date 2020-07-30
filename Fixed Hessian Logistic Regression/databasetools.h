#pragma once
#ifndef DATABASE_TOOLS
#define DATABASE_TOOLS

#include <vector>
#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <cstdlib>
#include "math.h"
#include"seal/seal.h"

using namespace std;
using namespace seal;

typedef vector<double> dVec;
typedef vector<vector<double>> dMat;
typedef vector<Ciphertext> cVec;
typedef vector<Plaintext> pVec;

int ImportData(dMat& Z, string filename);
int ImportDataLR(dMat& Z, string filename);
double inner_prod(dVec v, dVec u, int start = 0);
void CVRandomSampling(dMat*& train, dMat*& test, dMat data);
void AllSum(Ciphertext encrypted, Ciphertext& allsum, int slot_count, shared_ptr<SEALContext> context, GaloisKeys gal_keys);
int ImportDataLR_eighth(dMat& Matrix, string filename);
int ImportDataLR_half(dMat& Matrix, string filename);
int ImportDataLR_yfirst(dMat& Z, string filename);
int ImportDataLR_half_yfirst(dMat& Matrix, string filename);

#endif
