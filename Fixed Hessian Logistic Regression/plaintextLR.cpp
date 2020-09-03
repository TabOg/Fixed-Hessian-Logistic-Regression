#include "logregtools.h"
#include "databasetools.h"
#include <iostream>

int Plaintext_LR(string filename,int iternum) {
	dMat data;	
	ImportDataLR(data, "edin.txt", false);
	dVec weights(data[0].size(), 0.0);
	for (int k = 0; k < iternum; k++) {
		LR_iteration(data, weights, 5, data.size(), data[0].size());
		cout << "iteration " << k + 1 << " accuracy: " << accuracy_LR(weights, data) << "\n";
		cout << "iteration " << k + 1 << " AUC: " << getAUC(weights, data) << "\n";
	}
	return 0;
}