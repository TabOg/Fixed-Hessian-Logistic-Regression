#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "logregtools.h"
#include "algorithm"

using namespace std;
using namespace seal;

int Fixed_Hessian_Chebyshev(bool ringdim) {
    cout << "Running Fixed Hessian with Chebyshev sigmoid approximation:\n";
    dMat Matrix;
    ImportDataLR(Matrix, "edin.txt", false, 2);

    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = ringdim ? 32768 : 65536;

    parms.set_poly_modulus_degree(poly_modulus_degree);
    vector<int> mod;
    int x = ringdim ? 17 : 41;
    mod.push_back(50);
    for (int i = 0; i < x; i++)mod.push_back(40);
    mod.push_back(50);
    ringdim ? 5 : 12;
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, mod));
    cout << "Generating context...";
    auto start = chrono::steady_clock::now();
    auto context = SEALContext::Create(parms, true, sec_level_type::none);
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
    dMatMat cvtrain, cvtest;
    CVrandomSampling(cvtrain, cvtest, Matrix);
    Matrix.clear();

    //set 40 bits of precision
    double scale = pow(2.0, 40);
    pVec dataplain;
    Plaintext plain;
    dVec input;
    Ciphertext ctemp;
    cVec dataenc;
    cVec H;
    Plaintext P_T1, P_T2;
    double t1, t2;
    Ciphertext Htemp, allsumtemp;
    cVec AllSum;
    cVec Beta;
    cVec dataencscale;
    dVec weights;
    double accuracy, auc;
    accuracy = 0;
    auc = 0;
    for (int l = 0; l < 5; l++) {
        cout << "starting fold " << l + 1 << "...\n";
        int n = cvtrain[l].size();
        int nfeatures = cvtrain[l][0].size();
        cout << "Encoding...";
        start = chrono::steady_clock::now();

        //encode the data, a column at a time. The ith element of data corresponds to the ith feature
        dataplain.clear();
        for (int i = 0; i < nfeatures; i++) {
            input.clear();
            for (int j = 0; j < cvtrain[l].size(); j++)input.push_back(cvtrain[l][j][i]);
            encoder.encode(input, scale, plain);
            dataplain.push_back(plain);
        }
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        cout << "Encrypting...";
        start = chrono::steady_clock::now();
        //encrypt these plaintexts        
        dataenc.clear();

        for (int i = 0; i < dataplain.size(); i++) {
            encryptor.encrypt(dataplain[i], ctemp);
            dataenc.push_back(ctemp);
        }
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        cout << "Creating H...";
        start = chrono::steady_clock::now();
        //creating H: first add all ciphertexts
        H.clear();
        for (int i = 1; i < nfeatures; i++)evaluator.add_inplace(ctemp, dataenc[1. * nfeatures - 1 - i]);
        H = dataenc;
        //now create H(i,i)
        for (int i = 0; i < nfeatures; i++) {
            evaluator.multiply_inplace(H[i], ctemp);
            evaluator.relinearize_inplace(H[i], relin_keys);
            evaluator.rescale_to_next_inplace(H[i]);
            //allsum H[i]
            allsumtemp = H[i];
            for (int j = 0; j < log2(slot_count); j++) {
                H[i] = allsumtemp;
                evaluator.rotate_vector(H[i], pow(2, j), gal_keys, H[i]);
                evaluator.add_inplace(allsumtemp, H[i]);
            }
            H[i] = allsumtemp;
        }
        //time to calculate 1/H(i,i) -- first we need our starting point, T1 + T2D
        cout << "Calculating 1/H(i)...";
        t1 = T1(0.25 * (nfeatures * n));
        t2 = T2(0.25 * (nfeatures * n));
        encoder.encode(t1, scale, P_T1);
        encoder.encode(t2, scale, P_T2);

        //T2 needs to be multiplied by something 1 level down, so:
        evaluator.mod_switch_to_next_inplace(P_T2);

        //T1 will need to be added to something 2 levels down, so:
        evaluator.mod_switch_to_next_inplace(P_T1);
        evaluator.mod_switch_to_next_inplace(P_T1);

        for (int i = 0; i < nfeatures; i++) {
            //negate and store a copy of H(i,i) for later:
            Htemp = H[i];
            evaluator.negate_inplace(Htemp);

            //find v0 = T1 + T2H(i,i):
            evaluator.multiply_plain_inplace(H[i], P_T2);
            evaluator.rescale_to_next_inplace(H[i]);
            P_T1.scale() = H[i].scale();
            evaluator.add_plain_inplace(H[i], P_T1);

            //now iterate Newton Raphson: each update is 2v - H[i,i]v^2
            for (int j = 0; j < 3; j++) {

                //first double and store the result
                evaluator.add(H[i], H[i], ctemp);

                //now square the current value, relin and rescale

                evaluator.square_inplace(H[i]);
                evaluator.relinearize_inplace(H[i], relin_keys);
                evaluator.rescale_to_next_inplace(H[i]);

                //now mod switch down our stored value Htemp, and multiply
                evaluator.mod_switch_to_inplace(Htemp, H[i].parms_id());
                evaluator.multiply_inplace(H[i], Htemp);
                evaluator.relinearize_inplace(H[i], relin_keys);
                evaluator.rescale_to_next_inplace(H[i]);

                //modify scale of ctemp = 2v, mod switch down, and add
                ctemp.scale() = H[i].scale();
                evaluator.mod_switch_to_inplace(ctemp, H[i].parms_id());
                evaluator.add_inplace(H[i], ctemp);
            }
        }

        //precompute AllSum: AllSum[i] = sum(zji/2)
        AllSum.clear();
        AllSum.reserve(nfeatures);
        for (int i = 0; i < nfeatures; i++) {

            allsumtemp = dataenc[i];
            for (int j = 0; j < log2(slot_count); j++) {
                ctemp = allsumtemp;
                evaluator.rotate_vector_inplace(ctemp, pow(2, j), gal_keys);
                evaluator.add_inplace(allsumtemp, ctemp);
            }
            AllSum.push_back(allsumtemp);
        }

        //compute first iteration: beta[i]=H[i]AllSum[i]
        Beta.clear();
        Beta.reserve(nfeatures);

        for (int i = 0; i < nfeatures; i++) {
            allsumtemp = AllSum[i];
            evaluator.mod_switch_to_inplace(allsumtemp, H[i].parms_id());
            evaluator.multiply_inplace(allsumtemp, H[i]);
            evaluator.relinearize_inplace(allsumtemp, relin_keys);
            evaluator.rescale_to_next_inplace(allsumtemp);
            Beta.push_back(allsumtemp);
        }


        //we will also need (a copy of) each data vector to be multiplied by 5/8: we do this now
        dataencscale = dataenc;
        encoder.encode(0.625, scale, plain);
        for (int i = 0; i < nfeatures; i++) {
            evaluator.multiply_plain_inplace(dataencscale[i], plain);
            evaluator.rescale_to_next_inplace(dataencscale[i]);
        }
        weights.clear();

        for (int i = 0; i < nfeatures; i++) {
            decryptor.decrypt(Beta[i], plain);
            encoder.decode(plain, input);
            weights.push_back(input[0]);
        }
        cout << weights.size();
        cout << "," << cvtrain[l][0].size() << "\n";
        cout << "1st iteration accuracy: " << accuracy_LR(weights, cvtrain[l], 2) << "%\n";
        cout << "1st iteration AUC: " << 100 * getAUC(weights, cvtrain[l], 2) << "%\n";
        //start of an iteration: we are performing the update beta[i] <- beta[i] + H[i](AllSum[i] -5/8sum Beta.z(j)/2.z(ji)/2

        for (int k = 2; k < x; k++) {

            //calculate the inner product: this is the only calculation where we can't go feature by feature
            evaluator.mod_switch_to_inplace(dataenc[0], Beta[0].parms_id());
            evaluator.multiply(dataenc[0], Beta[0], Htemp);
            evaluator.relinearize_inplace(Htemp, relin_keys);
            evaluator.rescale_to_next_inplace(Htemp);

            for (int i = 1; i < nfeatures; i++) {
                evaluator.mod_switch_to_inplace(dataenc[i], Beta[i].parms_id());
                evaluator.multiply(dataenc[i], Beta[i], ctemp);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);
                evaluator.add_inplace(Htemp, ctemp);
            }
            /*cout << "6\n";*/
            //now we have a vector Htemp {beta.zi/2} i=1,...,n. so we can evaluate the update circuit
            //feature by feature
            for (int i = 0; i < nfeatures; i++) {
                /*cout << "7\n";*/
                //first modify 5/8zji/2 so that it can be multiplied by the inner product
                evaluator.mod_switch_to_inplace(dataencscale[i], Htemp.parms_id());
                //now multiply, relin, rescale
                evaluator.multiply(dataencscale[i], Htemp, ctemp);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);

                //now we need to allsum this vector:
                for (int j = 0; j < log2(slot_count); j++) {
                    allsumtemp = ctemp;
                    evaluator.rotate_vector_inplace(allsumtemp, pow(2, j), gal_keys);
                    evaluator.add_inplace(ctemp, allsumtemp);
                }

                //subtract this from AllSum[i], first switching down & modifying scale
                allsumtemp = AllSum[i];
                evaluator.mod_switch_to_inplace(allsumtemp, ctemp.parms_id());
                allsumtemp.scale() = ctemp.scale();
                evaluator.negate_inplace(ctemp);
                evaluator.add_inplace(ctemp, allsumtemp);


                //now multiply by H[i]
                evaluator.mod_switch_to_inplace(H[i], ctemp.parms_id());
                evaluator.multiply_inplace(ctemp, H[i]);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);

                //and finally update Beta[i]
                evaluator.mod_switch_to_inplace(Beta[i], ctemp.parms_id());
                Beta[i].scale() = ctemp.scale();
                evaluator.add_inplace(Beta[i], ctemp);
            }
            weights.clear();
            for (int i = 0; i < nfeatures; i++) {
                decryptor.decrypt(Beta[i], plain);
                encoder.decode(plain, input);
                weights.push_back(input[0]);

            }
            cout << "iteration " << k << " accuracy is: " << accuracy_LR(weights, cvtrain[l], 2) << "%\n";
            cout << "AUC is: " << 100 * getAUC(weights, cvtrain[l], 2) << "%\n";
        }
        cout << "final level is: " << context->get_context_data(Beta[0].parms_id())->chain_index();

        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "fold " << l + 1 << " training done. Time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
        cout << "CV accuracy is " << accuracy_LR(weights, cvtest[l], 2) << "%\n";
        cout << "CV AUC is " << getAUC(weights, cvtest[l], 2) << "\n";
        accuracy += accuracy_LR(weights, cvtest[l], 2);
        auc += getAUC(weights, cvtest[l], 2);
    }
    cout << "average CV accuracy is: " << accuracy / 5 << "%\n";
    cout << "average CV AUC is: " << auc / 5 << "%\n";
    return 0;   
}

int Fixed_Hessian_Taylor() {
    cout << "Running Fixed Hessian with Taylor sigmoid approximation:\n";
    dMat Matrix;
    ImportDataLR(Matrix, "edin.txt",false,2);
    int n = Matrix.size();
    int nfeatures = Matrix[0].size();

    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 32768;

    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 50,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,50 }));
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
    cout << "Encoding...";
    start = chrono::steady_clock::now();

    //set 40 bits of precision
    double scale = pow(2.0, 40);

    //encode the data, a column at a time. The ith element of data corresponds to the ith feature
    pVec dataplain;
    dataplain.reserve(nfeatures);
    Plaintext plain;
    dVec input;
    for (int i = 0; i < nfeatures; i++) {
        input.clear();
        for (int j = 0; j < Matrix.size(); j++)input.push_back(Matrix[j][i]);
        encoder.encode(input, scale, plain);
        dataplain.push_back(plain);
    }

    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    cout << "Encrypting...";
    start = chrono::steady_clock::now();
    //encrypt these plaintexts
    cVec dataenc;
    dataenc.reserve(nfeatures);
    Ciphertext ctemp;
    for (int i = 0; i < dataplain.size(); i++) {
        encryptor.encrypt(dataplain[i], ctemp);
        dataenc.push_back(ctemp);
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";
    cout << "Creating H...";
    start = chrono::steady_clock::now();
    //creating H: first add all ciphertexts
    cVec H;

    H.reserve(nfeatures);

    for (int i = 1; i < nfeatures; i++)evaluator.add_inplace(ctemp, dataenc[1. * nfeatures - 1 - i]);

    //now create H(i,i)
    Ciphertext Htemp, allsumtemp;
    for (int i = 0; i < nfeatures; i++) {
        evaluator.multiply(dataenc[i], ctemp, Htemp);
        evaluator.relinearize_inplace(Htemp, relin_keys);
        evaluator.rescale_to_next_inplace(Htemp);
        //allsum Htemp
        allsumtemp = Htemp;
        for (int j = 0; j < log2(slot_count); j++) {
            Htemp = allsumtemp;
            evaluator.rotate_vector(Htemp, pow(2, j), gal_keys, Htemp);
            evaluator.add_inplace(allsumtemp, Htemp);
        }
        //push back H(i,i) to the vector H
        H.push_back(allsumtemp);
    }
    //time to calculate 1/H(i,i) -- first we need our starting point, T1 + T2D
    cout << "Calculating 1/H(i)...";

    Plaintext P_T1, P_T2;
    double t1, t2;
    t1 = T1(2.5 * n);
    t2 = T2(2.5 * n);
    encoder.encode(t1, scale, P_T1);
    encoder.encode(t2, scale, P_T2);

    //T2 needs to be multiplied by something 1 level down, so:
    evaluator.mod_switch_to_next_inplace(P_T2);

    //T1 will need to be added to something 2 levels down, so:
    evaluator.mod_switch_to_next_inplace(P_T1);
    evaluator.mod_switch_to_next_inplace(P_T1);

    for (int i = 0; i < nfeatures; i++) {

        //negate and store a copy of H(i,i) for later:
        Htemp = H[i];
        evaluator.negate_inplace(Htemp);

        //find v0 = T1 + T2H(i,i):
        evaluator.multiply_plain_inplace(H[i], P_T2);
        evaluator.rescale_to_next_inplace(H[i]);
        P_T1.scale() = H[i].scale();
        evaluator.add_plain_inplace(H[i], P_T1);

        //now iterate Newton Raphson: each update is 2v - H[i,i]v^2
        for (int j = 0; j < 3; j++) {

            //first double and store the result
            evaluator.add(H[i], H[i], ctemp);

            //now square the current value, relin and rescale

            evaluator.square_inplace(H[i]);
            evaluator.relinearize_inplace(H[i], relin_keys);
            evaluator.rescale_to_next_inplace(H[i]);

            //now mod switch down our stored value Htemp, and multiply
            evaluator.mod_switch_to_inplace(Htemp, H[i].parms_id());
            evaluator.multiply_inplace(H[i], Htemp);
            evaluator.relinearize_inplace(H[i], relin_keys);
            evaluator.rescale_to_next_inplace(H[i]);

            //modify scale of ctemp = 2v, mod switch down, and add
            ctemp.scale() = H[i].scale();
            evaluator.mod_switch_to_inplace(ctemp, H[i].parms_id());
            evaluator.add_inplace(H[i], ctemp);
        }
    }
    //precompute AllSum: AllSum[i] = sum(zji/2)
    cVec AllSum;
    AllSum.reserve(nfeatures);
    for (int i = 0; i < nfeatures; i++) {

        allsumtemp = dataenc[i];
        for (int j = 0; j < log2(slot_count); j++) {
            ctemp = allsumtemp;
            evaluator.rotate_vector_inplace(ctemp, pow(2, j), gal_keys);
            evaluator.add_inplace(allsumtemp, ctemp);
        }
        AllSum.push_back(allsumtemp);
    }
    
    //compute first iteration: beta[i]=H[i]AllSum[i]
    cVec Beta;

    Beta.reserve(nfeatures);

    for (int i = 0; i < nfeatures; i++) {
        allsumtemp = AllSum[i];
        evaluator.mod_switch_to_inplace(allsumtemp, H[i].parms_id());
        evaluator.multiply_inplace(allsumtemp, H[i]);
        evaluator.relinearize_inplace(allsumtemp, relin_keys);
        evaluator.rescale_to_next_inplace(allsumtemp);
        Beta.push_back(allsumtemp);
    }
    
    //for the Taylor polynomial, each update only uses AllSum multiplied by H(i)
    AllSum = Beta;

    //we will also need (a copy of) each data vector to be multiplied by H(i): we do this now
    cVec dataencscale;
    dataencscale = dataenc;
    
    for (int i = 0; i < nfeatures; i++) {
        evaluator.mod_switch_to_inplace(dataencscale[i], H[i].parms_id());
        evaluator.multiply_inplace(dataencscale[i], H[i]);
        evaluator.relinearize_inplace(dataencscale[i], relin_keys);
        evaluator.rescale_to_next_inplace(dataenc[i]);

    }
    dVec weights(nfeatures, 0.0);
    
    for (int i = 0; i < nfeatures; i++) {
        decryptor.decrypt(Beta[i], plain);
        encoder.decode(plain, input);
        for (int j = 0; j < input.size(); j++)weights[i] += input[j];
        weights[i] /= input.size();
    }
    Matrix.clear();
    ImportDataLR(Matrix, "edin.txt",false);
    cout << "accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
    cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
    //start of an iteration: we are performing the update beta[i] <- beta[i] + H[i]AllSum[i] -sum Beta.z(j)/2.z(ji)/2

    for (int k = 2; k < 5; k++) {

        //calculate the inner product: this is the only calculation where we can't go feature by feature
        evaluator.mod_switch_to_inplace(dataenc[0], Beta[0].parms_id());
        evaluator.multiply(dataenc[0], Beta[0], Htemp);
        evaluator.relinearize_inplace(Htemp, relin_keys);
        evaluator.rescale_to_next_inplace(Htemp);
        
        for (int i = 1; i < nfeatures; i++) {
            evaluator.mod_switch_to_inplace(dataenc[i], Beta[i].parms_id());
            evaluator.multiply(dataenc[i], Beta[i], ctemp);
            evaluator.relinearize_inplace(ctemp, relin_keys);
            evaluator.rescale_to_next_inplace(ctemp);
            evaluator.add_inplace(Htemp, ctemp);
        }
        
        //now we have a vector Htemp {beta.zi/2} i=1,...,n. so we can evaluate the update circuit
        //feature by feature
        for (int i = 0; i < nfeatures; i++) {
            
            //first modify H(i)zji/2 so that it can be multiplied by the inner product
            evaluator.mod_switch_to_inplace(dataencscale[i], Htemp.parms_id());
            //now multiply, relin, rescale
            evaluator.multiply(dataencscale[i], Htemp, ctemp);
            evaluator.relinearize_inplace(ctemp, relin_keys);
            evaluator.rescale_to_next_inplace(ctemp);
            
            //now we need to allsum this vector:
            for (int j = 0; j < log2(slot_count); j++) {
                allsumtemp = ctemp;
                evaluator.rotate_vector_inplace(allsumtemp, pow(2, j), gal_keys);
                evaluator.add_inplace(ctemp, allsumtemp);
            }
            
            //subtract this from AllSum[i], first switching down & modifying scale
            allsumtemp = AllSum[i];
            evaluator.mod_switch_to_inplace(allsumtemp, ctemp.parms_id());
            allsumtemp.scale() = ctemp.scale();
            evaluator.negate_inplace(ctemp);
            evaluator.add_inplace(ctemp, allsumtemp);

            //and finally update Beta[i]
            evaluator.mod_switch_to_inplace(Beta[i], ctemp.parms_id());
            Beta[i].scale() = ctemp.scale();
            evaluator.add_inplace(Beta[i], ctemp);
        }
        for (int i = 0; i < nfeatures; i++) {
            decryptor.decrypt(Beta[i], plain);
            encoder.decode(plain, input);
            for (int j = 0; j < input.size(); j++)weights[i] += input[j];
            weights[i] /= input.size();
        }
        Matrix.clear();
        ImportDataLR(Matrix, "edin.txt",false);
        cout << "iteration " << k << " accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
        cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
    }
    cout << "final level is: " << context->get_context_data(Beta[0].parms_id())->chain_index();
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "done. Time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    for (int i = 0; i < nfeatures; i++) {
        decryptor.decrypt(Beta[i], plain);
        encoder.decode(plain, input);
        for (int j = 0; j < input.size(); j++)weights[i] += input[j];
        weights[i] /= input.size();
    }
    Matrix.clear();
    ImportDataLR(Matrix, "edin.txt",false);
    cout << "accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
    cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
    return 0;
}