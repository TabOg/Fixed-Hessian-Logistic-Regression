#include "seal/seal.h"
#include "databasetools.h"
#include "logregtools.h"
#include <iostream>
#include "threadpool.hpp"

using namespace std;
using namespace seal;

int GD(bool ringdim) {
    thread_pool::thread_pool tp(std::thread::hardware_concurrency());
    double accuracy = 0;
    double auc = 0;
    dMat Matrix;
    ImportDataLR(Matrix, "edin.txt", false, 8);
    dMatMat cvtrain, cvtest;
    CVrandomSampling(cvtrain, cvtest, Matrix);
    Matrix.clear();
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = ringdim? 32768:65536;
    vector<int> mod;
    mod.push_back(38);
    for (int i = 0; i < ringdim ? 25 : 57; i++)mod.push_back(28);
    mod.push_back(38);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, mod));
    cout << "Generating context..."<<endl;
    auto start = chrono::steady_clock::now();
    auto context = SEALContext::Create(parms, true,sec_level_type::none);
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
    cout << "KeyGen time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n"<<endl;
    CKKSEncoder encoder(context);

    size_t slot_count = encoder.slot_count();
    double a1 = -1.20096;
    double a3 = 0.81562;
    double a4 = a1 / a3;
    double sc;
    double scale = pow(2.0, 28);
    cout << "encoding polynomial coefficients..."<<endl;
    Plaintext coeff3, coeff4;
    encoder.encode(2 * a3, scale, coeff3);
    encoder.encode(a4, scale, coeff4);
    cout << "done \n";
    cout << "Number of slots: " << slot_count << "\n";

    int n, nfeatures;

    pVec data;
    Plaintext plain, scaler, scaler1,coeff4temp;
    dVec input,weights;
    cVec dataenc, ctsum;
    Ciphertext datatemp;
    Ciphertext allsum, temp;
    cVec Beta;
    for (int l = 0; l < 5; l++) {
        cout << "starting fold " << l + 1 << "...\n";
        
        Beta.clear();
        ctsum.clear();

        data.clear();
        n = cvtrain[l].size();
        nfeatures = cvtrain[l][0].size();
        sc = 4.0 / (1.0 * n);
        cout << "Encoding..."<<endl;
        start = chrono::steady_clock::now();
        for (int i = 0; i < nfeatures; i++) {
            input.clear();
            for (int j = 0; j < cvtrain[l].size(); j++)input.push_back(cvtrain[l][j][i]);
            encoder.encode(input, scale, plain);
            data.push_back(plain);
        }
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        start = chrono::steady_clock::now();
        cout << "Encrypting..."<<endl;
        dataenc.clear();
        for (int i = 0; i < nfeatures; i++) {
            encryptor.encrypt(data[i], datatemp);
            dataenc.push_back(datatemp);
        }
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

        auto begin = chrono::steady_clock::now();
        cout << "calculating ct.sum...\n";
        start = chrono::steady_clock::now();

        //creating the vector ctsum, ctsum[j] = 4/nsum(z_{ij}/8)
        encoder.encode(sc, scale, scaler);
        encoder.encode(5 * sc, scale, scaler1);
        
        std::mutex mutex;
        Beta.resize(nfeatures);
        ctsum.resize(nfeatures);
        
    
        for (int i = 0; i < nfeatures; i++) {
            tp.push([i, &Beta, &ctsum, &mutex, &evaluator, &scaler, &scaler1,  &dataenc, slot_count, &gal_keys] {   
            
                auto allsum = dataenc[i];

                for (int j = 0; j < log2(slot_count); j++) {
                    auto temp = allsum;
                    evaluator.rotate_vector(temp, pow(2, j), gal_keys, temp);
                    evaluator.add_inplace(allsum, temp);
                }

                auto temp = allsum;
                // Writes to allsum and temp
                evaluator.multiply_plain_inplace(allsum, scaler);
                evaluator.multiply_plain_inplace(temp, scaler1);
                evaluator.rescale_to_next_inplace(allsum);
                evaluator.rescale_to_next_inplace(temp);
                std::lock_guard<std::mutex> lock(mutex);
                ctsum[i] = allsum;
                Beta[i] = temp;
            });
        }

        tp.wait_work();

        weights.clear();
        for (int k = 0; k < nfeatures; k++) {
            decryptor.decrypt(Beta[k], plain);
            encoder.decode(plain, input);
            weights.push_back(input[0]);
        }
        cout << weights.size() << ", " << cvtrain[l][0].size() << "\n";
        cout << "fold " << l + 1 << " 1st iteration AUC is " << 100 * getAUC(weights, cvtrain[l], 8) << "%,";
        cout << "1st iteration accuracy is " << accuracy_LR(weights, cvtrain[l], 8) << "%\n";
        cVec dataencscale = dataenc;
        for (int i = 0; i < dataenc.size(); i++) {
            evaluator.multiply_plain_inplace(dataencscale[i], coeff3);
            evaluator.rescale_to_next_inplace(dataencscale[i]);
        }        
        Plaintext alphap, p;
        Ciphertext mult, innerprod, square;
        double alpha;
        //iterations
        for (int k = 2; k < ringdim ? 8 : 16; k++)
        {
            alpha = 10 / (k + 1);
            encoder.encode(alpha * sc, scale, scaler);
            encoder.encode(alpha, scale, alphap);
            evaluator.mod_switch_to_next_inplace(alphap);
            //creating the inner product vector
            evaluator.mod_switch_to(dataenc[0], Beta[0].parms_id(), innerprod);
            //start with b0*z[i][0]
            evaluator.multiply_inplace(innerprod, Beta[0]);
            evaluator.relinearize_inplace(innerprod, relin_keys);
            evaluator.rescale_to_next_inplace(innerprod);
            std::vector<Ciphertext> mults(nfeatures - 1);
            for (int i = 1; i < nfeatures; i++) {
                tp.push([&Beta, &mults, &mutex, &dataenc, i, &evaluator, &relin_keys] {
                    Ciphertext mult;
                    evaluator.mod_switch_to(dataenc[i], Beta[i].parms_id(), mult);
                    evaluator.multiply_inplace(mult, Beta[i]);
                    evaluator.relinearize_inplace(mult, relin_keys);
                    evaluator.rescale_to_next_inplace(mult);
                    std::lock_guard<std::mutex> lock(mutex);
                    mults[i - 1] = mult;
                    });
            }

            tp.wait_work();

            for (const auto& mult : mults)
            {
                evaluator.add_inplace(innerprod, mult);
            }

            //time to evaluate the polynomial!
            square = innerprod;
            evaluator.square_inplace(square);
            evaluator.relinearize_inplace(square, relin_keys);
            evaluator.rescale_to_next_inplace(square);
            evaluator.mod_switch_to(coeff4, square.parms_id(), coeff4temp);
            coeff4temp.scale() = square.scale();
            evaluator.add_plain_inplace(square, coeff4temp);
            //we now have wi^2 + a1/a3 -- calculate the value of the polynomial for each j
            mults.resize(nfeatures);
            // Note: this loop is done in parallel: we don't have any dependencies, so we can use a temporary ciphertext
            // (called mult) in each thread.
            for (int j = 0; j < nfeatures; j++) {
                tp.push([&evaluator, &mutex, &dataencscale, j, innerprod, &relin_keys, &square, &mults] {
                    Ciphertext mult;
                    // Writes to first
                    evaluator.mod_switch_to_inplace(dataencscale[j], innerprod.parms_id());
                    // Writes to last
                    evaluator.multiply(innerprod, dataencscale[j], mult);
                    // Writes to first
                    evaluator.relinearize_inplace(mult, relin_keys);
                    // Writes to first
                    evaluator.rescale_to_next_inplace(mult);
                    /*cout << context->get_context_data(square.parms_id())->chain_index() << "\n";
                    cout << context->get_context_data(squarea4.parms_id())->chain_index() << "\n";*/
                    evaluator.multiply_inplace(mult, square);
                    evaluator.relinearize_inplace(mult, relin_keys);
                    evaluator.rescale_to_next_inplace(mult);

                    std::lock_guard<std::mutex> lock(mutex);
                    mults[j] = mult;
                    });
            }

            tp.wait_work();

            // This loop has a dependency on scaler (mod_switch_to_implace writes to it), and 
            // I couldn't see if there were any interactions (i.e if scaler is recycled across iterations: it shouldn't but 
            // you never know.)
            // So, for safety we can only do this bit in serial
            for (auto& mult : mults) {
                evaluator.mod_switch_to_inplace(scaler, mult.parms_id());
                evaluator.multiply_plain_inplace(mult, scaler);
                evaluator.rescale_to_next_inplace(mult);
            }

            // We do this loop in parallel too: again we use a mutex to control
            // writes to mults
            for (int j = 0; j < nfeatures; j++) {
                tp.push([&mults, &mutex, &evaluator, j, &gal_keys, slot_count]() {
                    Ciphertext mult = mults[j];
                    for (int i = 0; i < log2(slot_count); i++) {
                        Ciphertext temp = mult;
                        evaluator.rotate_vector_inplace(temp, pow(2, i), gal_keys);
                        evaluator.add_inplace(mult, temp);
                    }
                    std::lock_guard<std::mutex> lock(mutex);
                    mults[j] = mult;
                    });
            }

            // Wait for this to finish before continuing
            tp.wait_work();


            for (unsigned int j = 0; j < nfeatures; j++) {
                tp.push([j, &Beta, &ctsum, &mults, &alphap, &evaluator, &mutex]() {

                    Ciphertext temp;
                    auto mult = mults[j];
                    auto beta = Beta[j];
                    auto ctsum_j = ctsum[j];
                    //switch down a copy of ctsum[j], add to all sum, update poly
                    evaluator.multiply_plain(ctsum_j, alphap, temp);
                    // Writes to temp
                    evaluator.rescale_to_next_inplace(temp);
                    temp.scale() = mult.scale();
                    evaluator.mod_switch_to_inplace(temp, mult.parms_id());
                    // Writes to previous mult
                    evaluator.add_inplace(mult, temp);
                    //update weights vector
                    mult.scale() = beta.scale();
                    // Writes to both Betas, but mult is unique per iteration

                    evaluator.mod_switch_to_inplace(beta, mult.parms_id());
                    evaluator.add_inplace(beta, mult);
                    // Synchronise to beta and ctsum
                    std::lock_guard<std::mutex> lock(mutex);
                    // These should be implicit moves in C++17.
                    beta.scale() = pow(2, 28);
                    Beta[j] = beta;
                    ctsum[j] = ctsum_j;
                    });
            }

            tp.wait_work();

            weights.clear();
            for (int i = 0; i < nfeatures; i++) {
                decryptor.decrypt(Beta[i], plain);
                encoder.decode(plain, input);
                weights.push_back(input[0]);
            }
            cout << "iteration " << k << " AUC is " << 100 * getAUC(weights, cvtrain[l], 8) << "%,";
            cout << "iteration " << k << " accuracy is " << accuracy_LR(weights, cvtrain[l], 8) << "%\n";
        }
        cout << "fold " << l + 1 << " CV accuracy is " << accuracy_LR(weights, cvtest[l], 8) << "%\n";
        cout << "fold " << l + 1 << " CV AUC is " << getAUC(weights, cvtest[l], 8) << "\n";
        accuracy += accuracy_LR(weights, cvtest[l], 8);
        auc += getAUC(weights, cvtest[l], 8);
        end = chrono::steady_clock::now();
        diff = end - start;

        cout << "fold " << l + 1 << " training done. Time =" << chrono::duration <double, milli>(diff).count() / 1000.0 << "s. \n";
    }
    cout << "CV average accuracy: " << accuracy / 5 << "%\n";
    cout << "CV average auc: " << auc / 5 << "\n";
    return 0;
}
