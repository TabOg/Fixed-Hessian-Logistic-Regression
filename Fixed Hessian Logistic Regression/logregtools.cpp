#include "logregtools.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "math.h"
#include "seal/seal.h"
#include "databasetools.h"


double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}
//define one iteration of logistic regression with preprocessing
int LR_iteration(dMat Matrix, dVec& weights, double learning_rate, int n, int nfeatures) {
    //temp grad vector
    dVec grad(nfeatures, 0);
    ////loop over the entries
    for (int i = 0; i < n; i++) {
        //find the value of the sigmoid function
        double sig = sigmoid(-inner_prod(Matrix[i], weights));
        //loop over the features, adding to the grad vector:
        for (int j = 0; j < nfeatures; j++) grad[j] += sig * Matrix[i][j];
    }
    //add to the weight vector
    for (int i = 0; i < nfeatures; i++) weights[i] += learning_rate * grad[i] / (1.0 * n);
    return 0;
}
//function combining import & plaintext logistic regression, writing resulting weights to a vector
int LR(dMat Matrix, dVec& weights, int max_iter, double learning_rate) {
    //extract number of records
    int n = Matrix.size();
    //extract number of weights needed (remember we need a constant weight as well)
    int nfeatures = Matrix[0].size();
    //iterate 
    for (int i = 0; i < max_iter; i++) {
        LR_iteration(Matrix, weights, learning_rate, n, nfeatures);
    }
    return 0;
}


int predict_LR(dVec weights, dVec sample, double threshold) {
    
    //compute probability
    double prob = sigmoid(sample[0]*inner_prod(weights, sample));
    //compare to threshold
    if (prob > threshold) return 1;
    else return -1;

}

double accuracy_LR(dVec weights, dMat test, double threshold) {
    double score = 0;
    for (int i = 0; i < test.size(); i++) {
        if (predict_LR(weights, test[i], threshold) == test[i][0])score++;
    }
    return 100 * score / test.size();
}

double getAUC(dVec theta, dMat zTest) {
    //initialise counters
    int n_fail_y1 = 0;
    int n_fail_y0 = 0;

    dVec xtheta_y1;
    dVec xtheta_y0;

    for (int i = 0; i < zTest.size(); ++i) {
        if (zTest[i][0] == 1.0) {
            if (inner_prod(zTest[i], theta) < 0) n_fail_y1++;
            xtheta_y1.push_back(zTest[i][0] * inner_prod(zTest[i], theta, 1));
        }
        else {
            if (inner_prod(zTest[i], theta) < 0) n_fail_y0++;
            xtheta_y0.push_back(zTest[i][0] * inner_prod(zTest[i], theta, 1));
        }
    }

    /*double correctness = 100.0 - (100.0 * (n_fail_y0 + n_fail_y1) / zTest.size());
    cout << "Failure rate: (y = 1) " << n_fail_y1 << "/" << xtheta_y1.size() << " + (y = 0) " << n_fail_y0 << "/";
    cout << xtheta_y0.size() << " = " << (100.0 * (n_fail_y0 + n_fail_y1) / zTest.size()) << " %." << endl;
    cout << "Correctness: " << correctness << " %." << endl;*/


    if (xtheta_y0.size() == 0 || xtheta_y1.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        return 0.0;
    }
    else {
        double auc = 0.0;
        for (int i = 0; i < xtheta_y1.size(); ++i) {
            for (int j = 0; j < xtheta_y0.size(); ++j) {
                if (xtheta_y0[j] <= xtheta_y1[i]) auc++;
            }
        }
        /*auc /= xtheta_y1.size() * xtheta_y0.size();*/
        cout << "AUC: " << auc << endl;
        return auc;

    }
}

double AUC(dVec weights, dMat test) {
    dVec D0;
    dVec D1;
    //sort into positive a negative samples, only pushing back the relevant dot product:
    for (int i = 0; i < test.size(); i++) {
        if (test[i][0] == -1)D0.push_back(test[i][0] * inner_prod(test[i], weights));
        else D1.push_back(test[i][0] * inner_prod(test[i], weights));
    }
    double auc = 0;
    for (int i = 0; i < D0.size(); i++) {
        for (int j = 0; j < D1.size(); j++) {
            if (D0[i] < D1[j]) auc++;
        }
    }
    auc /= D0.size();
    auc /= D1.size();
    return auc;
}

//one iteration of Nesterov's gradient descent: not great to be extracting dimensions each time? 
int LR_NV_iteration(dMat train, dVec& beta, dVec& v, double alpha, double gamma, int n, int nfeatures) {
    //first find J(v), since this is the only time we need the whole training set:
    dVec J(nfeatures, 0);
    double sig;
    for (int i = 0; i < n; i++) {
        //find the value of the sigmoid function
        sig = /*sigmoid(-inner_prod(train[i], v));*/1 / 2 - 1.20096 * inner_prod(train[i], v) / 8 + 0.81562 * pow(inner_prod(train[i], v) / 8, 3);
        sig /= n;
        //loop over the features, adding to the grad vector:
        for (int j = 0; j < nfeatures; j++) J[j] += sig * train[i][j];
    }
    dVec temp(nfeatures, 0);
    for (int i = 0; i < nfeatures; i++) {

        temp[i] = v[i] + alpha * J[i];
        v[i] = (1 - gamma) * temp[i] + gamma * beta[i];

    }
    beta = temp;
    return 0;
}

int LR_NV(dMat train, dVec& beta, dVec& v, int max_iter) {
    int n = train.size();
    int nfeatures = train[0].size();
    double t = 1.;
    double T, theta;
    for (int i = 0; i < max_iter; i++) {
        T = updatet(t);
        theta = -(t - 1) / T;
        LR_NV_iteration(train, beta, v, 10 / (i + 1), theta, n, nfeatures);
    }
    return 0;
}

//this function calculates t(k+1) from tk. Recall gammat = -(t(k)-1)/t(k+1)
double updatet(double t) {
    double T = (1. + sqrt(1. + 4 * t * t)) / 2.;
    return T;
}

double T1(double X) {
    return 8 * (X + 1) / (X * X + 6 * X + 1);
}

double T2(double X) {
    return (-1*8) / (X * X + 6 * X + 1);
}
