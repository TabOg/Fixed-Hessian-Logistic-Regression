#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "math.h"
#include "seal/seal.h"
#include "time.h"
#include "databasetools.h"

using namespace std;
using namespace seal;

int ImportData(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int long ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		entry1.push_back(stod(record[ncolumns - 1]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]));

		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		entryi.push_back(stod(record[ncolumns - 1]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]));
		//add it to the matrix
		Matrix.push_back(entryi);
	}

	return 0;
}
int ImportDataLR(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		entry1.push_back(stod(record[1.*ncolumns - 1]) * 2 - 1);
		//preprocessing for logistic regression
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]) * entry1[0]);
		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		entryi.push_back(stod(record[ncolumns - 1]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]) * entryi[0]);
		//add it to the matrix
		Matrix.push_back(entryi);
	}

	return 0;

}
double inner_prod(dVec v, dVec u, int start) {
	if (v.size() != u.size()) {
		cout << "error - vectors are different sizes";
		return 0;
	}
	else {
		double prod = 0.0;
		for (int i = start; i < u.size(); i++) prod += v[i] * u[i];
		return prod;
	}

}

void CVRandomSampling(dMat*& train, dMat*& test, dMat data) {

	int n = data.size();
	int n_test[5];
	//generate data fold sizes, equal sizes but losing up to 4 records
	int m = floor(n / 5);
	n_test[0] = m;
	n_test[1] = m;
	n_test[2] = m;
	n_test[3] = m;
	n_test[4] = n - 4 * m;
	//label all pieces of vector as "unchosen"
	vector<bool> randombits(n, false);

	//form test[i] for i = 0,...,3
	for (int i = 0; i < 4; i++) {
		//start a counter
		int k = 0;
		while (k < n_test[i]) {
			//sample a random number from [data.size()]
			int j = rand() % data.size();
			//if it's unchosen, add it to the fold & change to true
			if (randombits[j] == false) {
				randombits[j] = true;
				test[i].push_back(data[j]);
				k++;
			}
		}
	}
	//add the remaining records to the fifth fold
	for (int i = 0; i < n; i++) {
		if (randombits[i] == false)test[4].push_back(data[i]);
	}

	//generate the training sets train[0],...,train[4]
	for (int m = 0; m < 5; ++m) {
		for (int l = m + 1; l < 5; ++l) {
			for (int i = 0; i < n_test[l]; ++i) {
				train[m].push_back(test[l][i]);
			}
		}

		for (int l = 0; l < m; ++l) {
			for (int i = 0; i < n_test[l]; ++i) {
				train[m].push_back(test[l][i]);
			}
		}
	}
}

void AllSum(Ciphertext encrypted, Ciphertext& allsum, int slot_count, shared_ptr<SEALContext> context, GaloisKeys gal_keys) {
	Evaluator evaluator(context);
	allsum = encrypted;
	Ciphertext temp = encrypted;
	for (int j = 1; j < slot_count; j++) {
		evaluator.rotate_vector(temp, 1, gal_keys, temp);
		evaluator.add(allsum, temp, allsum);
	}
}

int ImportDataLR_eighth(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		//divide first entry by eight, preprocessing
		entry1.push_back((stod(record[ncolumns - 1]) * 2 - 1) / 8);
		//preprocessing for logistic regression
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]) * entry1[0]);
		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		//dividing by 8 again
		entryi.push_back((stod(record[ncolumns - 1]) * 2 - 1) / 8);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]) * entryi[0]);
		//add it to the matrix
		Matrix.push_back(entryi);
	}
	return 0;

}


int ImportDataLR_half(dMat& Matrix, string filename){
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		//divide first entry by two, preprocessing
		entry1.push_back((stod(record[ncolumns - 1]) * 2 - 1) / 2);
		//preprocessing for logistic regression
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]) * entry1[0]);
		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		//dividing by 2 again
		entryi.push_back((stod(record[ncolumns - 1]) * 2 - 1) / 2);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]) * entryi[0]);
		//add it to the matrix
		Matrix.push_back(entryi);
	}
	return 0;
}
int ImportDataLR_yfirst(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		entry1.push_back(stod(record[0]) * 2 - 1);
		//preprocessing for logistic regression
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]) * entry1[0]);
		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		entryi.push_back(stod(record[0]) * 2 - 1);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]) * entryi[0]);
		//add it to the matrix
		Matrix.push_back(entryi);
	}

	return 0;

}

int ImportDataLR_half_yfirst(dMat& Matrix, string filename) {
	//open file
	ifstream inFile;
	inFile.open(filename);
	//check file is open
	if (!inFile) {
		cout << "unable to open file";
		return 0;
	}

	string line;
	char split_char = '\t';
	int ncolumns;
	//process first row, moving class to the front and extracting number of columns
	if (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		ncolumns = record.size();
		vector<double> entry1;
		//divide first entry by two, preprocessing
		entry1.push_back((stod(record[0]) * 2 - 1) / 2);
		//preprocessing for logistic regression
		for (int i = 0; i < ncolumns - 1; i++) entry1.push_back(stod(record[i]) * entry1[0]);
		//add to matrix
		Matrix.push_back(entry1);
	}
	else {
		cout << "could not read file" << exit;
	}
	//process rest of the data
	while (getline(inFile, line)) {
		istringstream split(line);
		vector<string> record;
		for (string entry; getline(split, entry, split_char); record.push_back(entry));
		//record should have the same number of features
		if (record.size() != ncolumns) {
			cout << "database dimension error" << exit;
		}
		//define a new entry
		vector<double> entryi;
		//dividing by 2 again
		entryi.push_back((stod(record[0]) * 2 - 1) / 2);
		for (int i = 0; i < ncolumns - 1; i++) entryi.push_back(stod(record[i]) * entryi[0]);
		//add it to the matrix
		Matrix.push_back(entryi);
	}
	return 0;
}
	



