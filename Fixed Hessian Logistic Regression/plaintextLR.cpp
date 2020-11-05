#include "logregtools.h"
#include "databasetools.h"
#include <iostream>

int Plaintext_LR(string filename,int iternum) {
	dMat data;	
	ImportDataLR(data, "edin.txt", false);
	dMatMat cvtrain, cvtest;
	CVrandomSampling(cvtrain, cvtest, data);
	cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrain[0].size() << "\n";
	cout << "2 -- " << cvtrain[1].size() << "--" << cvtrain[1].size() << "\n";
	cout << "3 -- " << cvtrain[2].size() << "--" << cvtrain[2].size() << "\n";
	cout << "4 -- " << cvtrain[3].size() << "--" << cvtrain[3].size() << "\n";
	cout << "5 -- " << cvtrain[4].size() << "--" << cvtrain[4].size() << "\n";

	cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtest[0].size() << "\n";
	cout << "2 -- " << cvtest[1].size() << "--" << cvtest[1].size() << "\n";
	cout << "3 -- " << cvtest[2].size() << "--" << cvtest[2].size() << "\n";
	cout << "4 -- " << cvtest[3].size() << "--" << cvtest[3].size() << "\n";
	cout << "5 -- " << cvtest[4].size() << "--" << cvtest[4].size() << "\n";
	dVec weights(data[0].size(), 0.0);
	double auc = 0;
	double accuracy = 0;
	for (int j = 0; j < 5; j++) {
		fill(weights.begin(), weights.end(), 0.0);
		for (int k = 0; k < iternum; k++) {
			LR_iteration(cvtrain[j], weights, 10 / (k + 2), cvtrain[j].size(), cvtrain[j][0].size());
			cout << "iteration " << k + 1 << " accuracy: " << accuracy_LR(weights, cvtrain[j]) << "\n";
			cout << "iteration " << k + 1 << " AUC: " << getAUC(weights, cvtrain[j]) << "\n";
		}
		cout << "fold " << j + 1 << " CV accuracy: " << accuracy_LR(weights, cvtest[j]) << "%\n";
		cout << "fold " << j + 1 << " CV AUC: " << getAUC(weights, cvtest[j]) << "\n";
		accuracy += accuracy_LR(weights, cvtest[j]);
		auc += getAUC(weights, cvtest[j]);
	}
	cout << "average CV accuracy: " << 0.2 * accuracy << "%\n";
	cout << "average CV AUC: " << 0.2* auc << "\n";
	return 0;
}

int Plaintext_LR_NAG(string filename, int iternum) {
	dMat data;
	ImportDataLR(data, "edin.txt", false);
	dMatMat cvtrain, cvtest;
	CVrandomSampling(cvtrain, cvtest, data);
	cout << "fold sizes are:\n1 -- " << cvtrain[0].size() << "--" << cvtrain[0].size() << "\n";
	cout << "2 -- " << cvtrain[1].size() << "--" << cvtrain[1].size() << "\n";
	cout << "3 -- " << cvtrain[2].size() << "--" << cvtrain[2].size() << "\n";
	cout << "4 -- " << cvtrain[3].size() << "--" << cvtrain[3].size() << "\n";
	cout << "5 -- " << cvtrain[4].size() << "--" << cvtrain[4].size() << "\n";

	cout << "test sizes are:\n1 -- " << cvtest[0].size() << "--" << cvtest[0].size() << "\n";
	cout << "2 -- " << cvtest[1].size() << "--" << cvtest[1].size() << "\n";
	cout << "3 -- " << cvtest[2].size() << "--" << cvtest[2].size() << "\n";
	cout << "4 -- " << cvtest[3].size() << "--" << cvtest[3].size() << "\n";
	cout << "5 -- " << cvtest[4].size() << "--" << cvtest[4].size() << "\n";
	dVec beta(data[0].size(), 0.0);
	dVec v(data[0].size(), 0.0);
	double auc = 0;
	double accuracy = 0;
	double T, t, theta;
	for (int j = 0; j < 5; j++) {
		t = 1;
		fill(beta.begin(), beta.end(), 0.0);
		fill(v.begin(), v.end(), 0.0);
		int n = cvtrain[j].size();
		int nfeatures = cvtrain[j][0].size();
		for (int k = 0; k < iternum; k++) {
			T = updatet(t);
			theta = -(t - 1) / T;
			LR_NV_iteration(cvtrain[j], beta, v, 10 / (k + 2), theta, n, nfeatures);
			cout << "iteration " << k + 1 << " accuracy: " << accuracy_LR(v, cvtrain[j]) << "\n";
			cout << "iteration " << k + 1 << " AUC: " << getAUC(v, cvtrain[j]) << "\n";
		}
		cout << "fold " << j + 1 << " CV accuracy: " << accuracy_LR(v, cvtest[j]) << "%\n";
		cout << "fold " << j + 1 << " CV AUC: " << getAUC(v, cvtest[j]) << "\n";
		accuracy += accuracy_LR(v, cvtest[j]);
		auc += getAUC(v, cvtest[j]);
	}
	cout << "average CV accuracy: " << 0.2 * accuracy << "%\n";
	cout << "average CV AUC: " << 0.2 * auc << "\n";
	return 0;
}