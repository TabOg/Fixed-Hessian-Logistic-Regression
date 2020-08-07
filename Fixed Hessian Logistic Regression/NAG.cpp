#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "logregtools.h"

using namespace std;
using namespace seal;

int Nesterov_GD() {
    cout << "Running Nesterov Accelerated Gradient Descent:\n";
    auto begin = chrono::steady_clock::now();
    dMat Matrix;
    dVec weights;
    ImportDataLR(Matrix, "edin.txt",false,8);
    int n = Matrix.size();
    int nfeatures = Matrix[0].size();
    //make the first eight columns one matrix
    dVec train1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 8; j++)train1.push_back(Matrix[i][j]);
    }
    //and the last 2 columns another
    dVec train2;
    for (int i = 0; i < n; i++) {
        for (int j = 8; j < 10; j++)train2.push_back(Matrix[i][j]);
        //fill the rows up with zeroes so it matches with above
        for (int j = 0; j < 6; j++)train2.push_back(0);
    }

    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 32768;

    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 40,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,40 }));
    cout << "Generating context...";
    auto start = chrono::steady_clock::now();
    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    RelinKeys relin_keys = keygen.relin_keys_local();
    GaloisKeys gal_keys = keygen.galois_keys_local();

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "KeyGen time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
    CKKSEncoder encoder(context);

    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << "\n";
    //ncol is the number of rows our matrix will have: we need to know this to do all sum on each column
    int nrow = slot_count / 8;
    cout << "Encoding...";
    start = chrono::steady_clock::now();
    double scale = pow(2.0, 30);
    
    Plaintext data1, data2;
    encoder.encode(train1, scale, data1);
    encoder.encode(train2, scale, data2);
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    Ciphertext dataenc1, dataenc2;
    start = chrono::steady_clock::now();
    cout << "Encrypting...";

    encryptor.encrypt(data1, dataenc1);
    encryptor.encrypt(data2, dataenc2);
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    double a0 = 0.5;
    double a1 = -1.20096;
    double a3 = 0.81562;
    double a4 = a1 / a3;
    double sc = 4.0 / (1.0 * n);
    double t = 1.;
    double T;
    Plaintext gammap, mgammap;

    cout << "encoding polynomial coefficients...";
    Plaintext coeff0, coeff1, coeff3, coeff4, coeff4temp, scaler;
    encoder.encode(a0, scale, coeff0);
    encoder.encode(a1, scale, coeff1);
    encoder.encode(2 * a3, scale, coeff3);
    encoder.encode(a4, scale, coeff4);
    encoder.encode(sc, scale, scaler);

    cout << "done \n";

    cout << "creating the matrix C...";
    //creating the matrix C -- start with the message version
    dVec Cm(n * 8.0, 0);
    for (int i = 0; i < n; i++)Cm[8.0 * i] += 1;
    //now plaintext version
    Plaintext Cp;
    encoder.encode(Cm, scale, Cp);

    cout << "done \n";
    cout << "creating 2a3zi/8...";
    Ciphertext poly1, poly2;
    poly1 = dataenc1;
    poly2 = dataenc2;
    evaluator.multiply_plain_inplace(poly1, coeff3);
    evaluator.multiply_plain_inplace(poly2, coeff3);
    evaluator.rescale_to_next_inplace(poly1);
    evaluator.rescale_to_next_inplace(poly2);
    cout << "done. \n";

    cout << "creating allsum matrices...";
    //this plays an analogous role to ct.sum. In column i, all entries are sum(z_ij/8)
    Ciphertext temp, ctsum1, ctsum2;
    ctsum1 = dataenc1;
    ctsum2 = dataenc2;
    
    for (int i = 0; i < log2(nrow); i++) {
        temp = ctsum1;
        evaluator.rotate_vector(temp, 8 * pow(2, i), gal_keys, temp);
        evaluator.add_inplace(ctsum1, temp);
    }
    
    for (int i = 0; i < log2(nrow); i++) {
        temp = ctsum2;
        evaluator.rotate_vector(temp, 8 * pow(2, i), gal_keys, temp);
        evaluator.add_inplace(ctsum2, temp);
    }
    

    //store a copy for the first iteration
    Ciphertext Beta1 = ctsum1;
    Ciphertext Beta2 = ctsum2;
    //scale by 4/n
    
    evaluator.multiply_plain_inplace(ctsum1, scaler);
    evaluator.rescale_to_next_inplace(ctsum1);
    evaluator.multiply_plain_inplace(ctsum2, scaler);
    evaluator.rescale_to_next_inplace(ctsum2);
    
    //This is the first iteration

    Ciphertext v1, v2;

    //learning rate is 10/2*4/n for the first iteration. 
    double sc1 = 5 * sc;
    Plaintext scaler1;
    encoder.encode(sc1, scale, scaler1);

    //Beta is updated to 5 ctsum
    evaluator.multiply_plain_inplace(Beta1, scaler1);
    evaluator.rescale_to_next_inplace(Beta1);
    evaluator.multiply_plain_inplace(Beta2, scaler1);
    evaluator.rescale_to_next_inplace(Beta2);

    //now update the v, recalling beta_0 is zero: 
    T = (1. + sqrt(1. + 4 * t * t)) / 2.;
    double gamma = -(t - 1) / T;
    double mgamma = 1 - gamma;
    encoder.encode(mgamma, scale, mgammap);
    t = T;
    evaluator.mod_switch_to_next_inplace(mgammap);

    evaluator.multiply_plain(Beta1, mgammap, v1);
    evaluator.multiply_plain(Beta2, mgammap, v2);
    evaluator.rescale_to_next_inplace(v1);
    evaluator.rescale_to_next_inplace(v2);

    dMat weightsmat;
    Plaintext p1, p2;
    dVec w1, w2;
    decryptor.decrypt(Beta1, p1);
    decryptor.decrypt(Beta2, p2);
    encoder.decode(p1, w1);
    encoder.decode(p2, w2);
    weights.clear();
    
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < 8; j++)weights.push_back(w1[8.0 * i + 1.0 * j]);
        for (int j = 0; j < 2; j++)weights.push_back(w2[8.0 * i + 1.0 * j]);
        weightsmat.push_back(weights);
        weights.clear();
    }

    weights = weightsmat[0];
    
    for (int i = 1; i < nrow; i++) {
        for (int j = 0; j < 10; j++)weights[j] += weightsmat[i][j];        
    }
    for (int i = 0; i < weights.size(); i++)weights[i] /= (1. * nrow);
    cout << "weights vector has "<<weights.size()<<" elements";
    Matrix.clear();
    ImportDataLR(Matrix, "edin.txt", false,1.0,'\t');
    cout << "1st iteration AUC is " << 100 * getAUC(weights, Matrix) << "%\n";
    cout << "1st accuracy is " << accuracy_LR(weights, Matrix) << "%";

    weightsmat.clear();
    weights.clear();

    //create all the temp variables we will need each iteration

    Plaintext alphap, Ctemp;
    Ciphertext allsum1, allsum2, poly1temp, poly2temp, weighttemp1, weighttemp2, datatemp;
    double alpha;
    //beginning of an iteration:
    for (int i = 2; i < 6; i++) {

        //encode learning rate alpha
        alpha = 10 / (1.0 * i + 1);
        encoder.encode(alpha, scale, alphap);
        //calculate and encode the weights gamma and 1 - gamma
        T = (1. + sqrt(1. + 4 * t * t)) / 2.;
        gamma = -(t - 1) / T;
        mgamma = 1 - gamma;
        encoder.encode(gamma, scale, gammap);
        encoder.encode(mgamma, scale, mgammap);

        t = T;
        
        weighttemp1 = v1;
        datatemp = dataenc1;
        //switch down encrypted data and multiply with the weight vector v

        evaluator.mod_switch_to_inplace(datatemp, weighttemp1.parms_id());
        evaluator.multiply_inplace(weighttemp1, datatemp);
        evaluator.relinearize_inplace(weighttemp1, relin_keys);
        evaluator.rescale_to_next_inplace(weighttemp1);

        //now for the last 2 columns
        weighttemp2 = v2;
        datatemp = dataenc2;

        evaluator.mod_switch_to_inplace(datatemp, weighttemp2.parms_id());
        evaluator.multiply_inplace(weighttemp2, datatemp);
        evaluator.relinearize_inplace(weighttemp2, relin_keys);
        evaluator.rescale_to_next_inplace(weighttemp2);

        //allsum to create inner products

        allsum1 = weighttemp1;
        allsum2 = weighttemp2;

        //all sum on the first vector
        
        for (int i = 0; i < 3; i++) {
            temp = allsum1;

            evaluator.rotate_vector_inplace(temp, pow(2, i), gal_keys);

            evaluator.add_inplace(allsum1, temp);
        }
        //all sum on the second
        
        temp = allsum2;

        evaluator.rotate_vector_inplace(temp, 1, gal_keys);

        evaluator.add_inplace(allsum2, temp);

        //adding together to get inner products in the first column
        evaluator.add_inplace(allsum1, allsum2);

        //multiply by C, first changing C's parameters 
        Ctemp = Cp;

        evaluator.mod_switch_to_inplace(Ctemp, allsum1.parms_id());
        evaluator.multiply_plain_inplace(allsum1, Ctemp);

        evaluator.rescale_to_next_inplace(allsum1);

        //replicating the inner product across the columns of an 8 x nrow matrix

        for (int i = 0; i < 3; i++) {
            temp = allsum1;
            evaluator.rotate_vector(temp, -pow(2, i), gal_keys, temp);
            evaluator.add_inplace(allsum1, temp);
        }

        //evaluating the polynomial. Start by squaring wi:

        allsum2 = allsum1;
        evaluator.square_inplace(allsum2);
        evaluator.relinearize_inplace(allsum2, relin_keys);
        evaluator.rescale_to_next_inplace(allsum2);
        
        //modify scale, mod switch a1/a3, and add:
        coeff4temp = coeff4;
        coeff4temp.scale() = allsum2.scale();
        evaluator.mod_switch_to_inplace(coeff4temp, allsum2.parms_id());
        evaluator.add_plain_inplace(allsum2, coeff4temp);
        
        //switch down (a copy of) 2a3zi/8, multiply by wi:
        poly1temp = poly1;

        evaluator.mod_switch_to_inplace(poly1temp, allsum1.parms_id());
        evaluator.multiply_inplace(poly1temp, allsum1);
        evaluator.relinearize_inplace(poly1temp, relin_keys);
        evaluator.rescale_to_next_inplace(poly1temp);
        //repeat for second half of matrix;
        poly2temp = poly2;

        evaluator.mod_switch_to_inplace(poly2temp, allsum1.parms_id());
        evaluator.multiply_inplace(poly2temp, allsum1);
        evaluator.relinearize_inplace(poly2temp, relin_keys);
        evaluator.rescale_to_next_inplace(poly2temp);

        //complete polynomial evaluation of (2g(x)-1)zij/8:

        evaluator.multiply_inplace(poly1temp, allsum2);
        evaluator.relinearize_inplace(poly1temp, relin_keys);
        evaluator.rescale_to_next_inplace(poly1temp);
        evaluator.multiply_inplace(poly2temp, allsum2);
        evaluator.relinearize_inplace(poly2temp, relin_keys);
        evaluator.rescale_to_next_inplace(poly2temp);

        //multiply ctsum1 and 2 by the learning rate
        allsum1 = ctsum1;

        evaluator.mod_switch_to_next_inplace(alphap);
        evaluator.multiply_plain_inplace(allsum1, alphap);
        evaluator.rescale_to_next_inplace(allsum1);
        allsum2 = ctsum2;
        evaluator.multiply_plain_inplace(allsum2, alphap);
        evaluator.rescale_to_next_inplace(allsum2);

        //time to allsum the columns.
        for (int i = 0; i < log2(nrow); i++) {
            temp = poly1temp;
            evaluator.rotate_vector(temp, 8 * pow(2, i), gal_keys, temp);
            evaluator.add_inplace(poly1temp, temp);
        }

        for (int i = 0; i < log2(nrow); i++) {
            temp = poly2temp;
            evaluator.rotate_vector(temp, 8 * pow(2, i), gal_keys, temp);
            evaluator.add_inplace(poly2temp, temp);
        }

        //multiply polytemp 1 and 2 by 4alpha/n:
        alpha = alpha * sc;
        encoder.encode(alpha, scale, alphap);

        evaluator.mod_switch_to_inplace(alphap, poly1temp.parms_id());
        evaluator.multiply_plain_inplace(poly1temp, alphap);
        evaluator.rescale_to_next_inplace(poly1temp);
        evaluator.multiply_plain_inplace(poly2temp, alphap);
        evaluator.rescale_to_next_inplace(poly2temp);

        //finally complete polynomial evaluation:
        
        allsum1.scale() = poly1temp.scale();
        allsum2.scale() = poly2temp.scale();
        evaluator.mod_switch_to_inplace(allsum1, poly1temp.parms_id());
        evaluator.mod_switch_to_inplace(allsum2, poly2temp.parms_id());
        evaluator.add_inplace(poly1temp, allsum1);
        evaluator.add_inplace(poly2temp, allsum2);

        weighttemp1 = Beta1;    
        weighttemp2 = Beta2;
        Beta1 = v1;
        Beta2 = v2;

        //switch v1 and v2 down and add:
        Beta1.scale() = poly1temp.scale();

        Beta2.scale() = poly2temp.scale();

        evaluator.mod_switch_to_inplace(Beta1, poly1temp.parms_id());
        evaluator.mod_switch_to_inplace(Beta2, poly2temp.parms_id());

        evaluator.add_inplace(Beta1, poly1temp);
        evaluator.add_inplace(Beta2, poly2temp);

        //first line of update done. Now for second
        //start by switching down the learning rates 1 - gamma and gamma       

        evaluator.mod_switch_to(mgammap, Beta1.parms_id(), mgammap);

        evaluator.mod_switch_to(gammap, weighttemp1.parms_id(), gammap);

        //multiply by 1 - gamma and gamma
        v1 = Beta1;
        v2 = Beta2;
        evaluator.multiply_plain_inplace(v1, mgammap);
        evaluator.rescale_to_next_inplace(v1);
        evaluator.multiply_plain_inplace(v2, mgammap);
        evaluator.rescale_to_next_inplace(v2);
        evaluator.multiply_plain_inplace(weighttemp1, gammap);
        evaluator.rescale_to_next_inplace(weighttemp1);
        evaluator.multiply_plain_inplace(weighttemp2, gammap);
        evaluator.rescale_to_next_inplace(weighttemp2);

        //switch down the "old" weights, so they can be added.
        weighttemp1.scale() = v1.scale();
        weighttemp2.scale() = v2.scale();
        evaluator.mod_switch_to(weighttemp1, v1.parms_id(), weighttemp1);
        evaluator.mod_switch_to(weighttemp2, v2.parms_id(), weighttemp2);

        //now update the vectors v1 and v2:
        evaluator.add_inplace(v1, weighttemp1);
        evaluator.add_inplace(v2, weighttemp2);


        decryptor.decrypt(v1, p1);
        decryptor.decrypt(v2, p2);
        encoder.decode(p1, w1);
        encoder.decode(p2, w2);
        weights.clear();
        
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < 8; j++)weights.push_back(w1[8.0 * i + 1.0 * j]);
            for (int j = 0; j < 2; j++)weights.push_back(w2[8.0 * i + 1.0 * j]);
            weightsmat.push_back(weights);
            weights.clear();
        }

        weights = weightsmat[0];
        
        for (int i = 1; i < nrow; i++) {
            for (int j = 0; j < nfeatures; j++)weights[j] += weightsmat[i][j];

        }
        for (int i = 0; i < weights.size(); i++)weights[i] /= (1. * nrow);
        cout << "weights vector has " << weights.size() << " elements\n";
        cout << "iteration "<<i<<" AUC is " << 100 * getAUC(weights, Matrix) << "%\n";
        cout << "iteration " << i<< " accuracy is " << accuracy_LR(weights, Matrix) << "%\n";

        weightsmat.clear();
        weights.clear();

    }
    auto ending = chrono::steady_clock::now();
    auto total = ending - begin;
    cout << "done. Time =" << chrono::duration <double, milli>(total).count() / 1000.0 << "s. \n";


    decryptor.decrypt(v1, p1);
    decryptor.decrypt(v2, p2);
    encoder.decode(p1, w1);
    encoder.decode(p2, w2);
    
    weights.clear();
    
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < 8; j++)weights.push_back(w1[8.0 * i + 1.0 * j]);
        for (int j = 0; j < 2; j++)weights.push_back(w2[8.0 * i + 1.0 * j]);
        weightsmat.push_back(weights);
        weights.clear();
    }

    weights = weightsmat[0];
    
    for (int i = 1; i < nrow; i++) {
        for (int j = 0; j < 10; j++)weights[j] += weightsmat[i][j];
    }
    for (int i = 0; i < weights.size(); i++)weights[i] /= (1. * nrow);
    cout << "final AUC is " << 100 * getAUC(weights, Matrix) << "%\n";
    cout << "final accuracy is " << accuracy_LR(weights, Matrix) << "%";
    return 0;
}