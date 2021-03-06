#include "algorithm"
#include "databasetools.h"
#include "logregtools.h"
#include "seal/seal.h"
#include "threadpool.hpp"
#include <iostream>

using namespace std;
using namespace seal;

// Create the set of parameters that are to be used for this file.
// Note: const so that mis-using the [] operator will cause an error.
static const std::unordered_map<unsigned int, ParameterPack> supported_parameters = {
    {39, ParameterPack(5, 49, 39, 22)},
    {40, ParameterPack(6, 50, 40, 18)},
    {45, ParameterPack(4, 55, 45, 19)}};

int Fixed_Hessian_Chebyshev(const unsigned precision, const unsigned number_of_threads){
    const auto parameter_iter = supported_parameters.find(precision);
    if (parameter_iter == supported_parameters.cend()) {
        cerr << "Error: unsupported parameter pack requested. Please try with one of "
          "the supported parameter sets.\n";
        exit(1);
    }

    // Else, grab the parameters and treat them as an rvalue for nicer calling
    const auto &parameter_set = parameter_iter->second;

    cout << "Running Fixed Hessian with Chebyshev sigmoid approximation:\n";
    dMat Matrix;
    ImportDataLR(Matrix, "edin.txt", false, 2);
    thread_pool::thread_pool tp(number_of_threads);
    EncryptionParameters parms(scheme_type::CKKS);
    constexpr size_t poly_modulus_degree = 32768;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    // Add 2 to count for the first and last prime in the vector
    vector<int> mod(parameter_set.number_of_primes + 2,
                    parameter_set.middle_prime_length);
    mod[0] = parameter_set.sentinel_prime_length;
    mod[parameter_set.number_of_primes + 1] = parameter_set.sentinel_prime_length;

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
    cout << "KeyGen time = "
         << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << "\n";
    dMatMat cvtrain, cvtest;
    CVrandomSampling(cvtrain, cvtest, Matrix);
    Matrix.clear();

    // set the desired bits of precision
    double scale = pow(2.0, precision);
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
    double accuracy{0.0}, auc{0.0};
    
    for (unsigned l = 0; l < 5; l++) {
        cout << "Starting fold " << l + 1 << "...\n";
        const unsigned int n = cvtrain[l].size();
        const unsigned int nfeatures = cvtrain[l][0].size();
        cout << "Encoding...";
        start = chrono::steady_clock::now();

        // encode the data, a column at a time. The ith element of data corresponds to
        // the ith feature
        dataplain.clear();
        dataplain.reserve(nfeatures);
        for (unsigned int i = 0; i < nfeatures; i++) {
            input.clear();
            for (unsigned int j = 0; j < cvtrain[l].size(); j++) {
                input.push_back(cvtrain[l][j][i]);
            }
            
            encoder.encode(input, scale, plain);
            dataplain.emplace_back(plain);
        }
        
        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encoding time = "
            << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";

        cout << "Encrypting...";
        start = chrono::steady_clock::now();
        
        // encrypt these plaintexts
        dataenc.clear();
        dataenc.reserve(dataplain.size());

        for(unsigned i = 0; i < dataplain.size(); i++) {
            encryptor.encrypt(dataplain[i], ctemp);
            dataenc.emplace_back(ctemp);
        }

        end = chrono::steady_clock::now();
        diff = end - start;
        cout << "Encrypting time = "
             << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
        cout << "Creating H...";
        start = chrono::steady_clock::now();

        // creating H: first add all ciphertexts
        H.clear();
        for (unsigned int i = 1; i < nfeatures; i++) {
            evaluator.add_inplace(ctemp, dataenc[1. * nfeatures - 1 - i]);
        }
        
        H = dataenc;
        // now create H(i,i)
        // We do this in parallel, treating each H[i] independently.
        std::mutex H_mutex;
        for (unsigned int i = 0; i < nfeatures; i++) {
            tp.push([&H, i, &H_mutex, &evaluator, &ctemp, &relin_keys, &gal_keys,
                    slot_count]() {
                auto local_temp = H[i];

                evaluator.multiply_inplace(local_temp, ctemp);
                evaluator.relinearize_inplace(local_temp, relin_keys);
                evaluator.rescale_to_next_inplace(local_temp);
                
                // allsum H[i]
                Ciphertext allsumtemp = local_temp;
                for (unsigned j = 0; j < log2(slot_count); j++) {
                    local_temp = allsumtemp;
                    evaluator.rotate_vector(local_temp, pow(2, j), gal_keys, local_temp);
                    evaluator.add_inplace(allsumtemp, local_temp);
                }
                
                std::lock_guard<std::mutex> lock(H_mutex);
                H[i] = allsumtemp;
            });
       }

        // Wait for the previous loop to finish
        tp.wait_work();
        // time to calculate 1/H(i,i) -- first we need our starting point, T1 + T2D
        cout << "Calculating 1/H(i)...";
        t1 = T1(0.25 * (nfeatures * n));
        t2 = T2(0.25 * (nfeatures * n));
        encoder.encode(t1, scale, P_T1);
        encoder.encode(t2, scale, P_T2);
        // T2 needs to be multiplied by something 1 level down, so:
        evaluator.mod_switch_to_next_inplace(P_T2);

        // T1 will need to be added to something 2 levels down, so:
        evaluator.mod_switch_to_next_inplace(P_T1);
        evaluator.mod_switch_to_next_inplace(P_T1);

        // Similarly to above we can parallelise this loop
        // We reuse h_mutex here too to make sure we write back safely
        for (unsigned i = 0; i < nfeatures; i++) {
            tp.push([&H, i, &H_mutex, &evaluator, &relin_keys, &P_T1, &P_T2]() {
                // negate and store a copy of H(i,i) for later:
                auto Htemp = H[i];
                auto local_hi = H[i];
                evaluator.negate_inplace(Htemp);

                // find v0 = T1 + T2H(i,i):
                evaluator.multiply_plain_inplace(local_hi, P_T2);
                evaluator.rescale_to_next_inplace(local_hi);
                P_T1.scale() = local_hi.scale();
                evaluator.add_plain_inplace(local_hi, P_T1);

                Ciphertext ctemp;
                // now iterate Newton Raphson: each update is 2v - H[i,i]v^2
                for (int j = 0; j < 3; j++) {
                    // first double and store the result
                    evaluator.add(local_hi, local_hi, ctemp);

                    // now square the current value, relin and rescale
                    evaluator.square_inplace(local_hi);
                    evaluator.relinearize_inplace(local_hi, relin_keys);
                    evaluator.rescale_to_next_inplace(local_hi);

                    // now mod switch down our stored value Htemp, and multiply
                    evaluator.mod_switch_to_inplace(Htemp, local_hi.parms_id());
                    evaluator.multiply_inplace(local_hi, Htemp);
                    evaluator.relinearize_inplace(local_hi, relin_keys);
                    evaluator.rescale_to_next_inplace(local_hi);

                    // modify scale of ctemp = 2v, mod switch down, and add
                    ctemp.scale() = local_hi.scale();
                    evaluator.mod_switch_to_inplace(ctemp, local_hi.parms_id());
                    evaluator.add_inplace(local_hi, ctemp);
                }
                
                std::lock_guard<std::mutex> lock(H_mutex);
                H[i] = local_hi;
            });
        }

        tp.wait_work();
        // precompute AllSum: AllSum[i] = sum(zji/2)
        AllSum.clear();
        AllSum.resize(nfeatures);
        std::mutex AllSumMutex;
        for (unsigned i = 0; i < nfeatures; i++) {
            tp.push([i, &evaluator, &dataenc, slot_count, &gal_keys, &AllSumMutex,
                 &AllSum]() {
                auto allsumtemp = dataenc[i];
                Ciphertext ctemp;
                for (int j = 0; j < log2(slot_count); j++) {
                    ctemp = allsumtemp;
                    evaluator.rotate_vector_inplace(ctemp, pow(2, j), gal_keys);
                    evaluator.add_inplace(allsumtemp, ctemp);
                }
                std::lock_guard<std::mutex> lock(AllSumMutex);
                AllSum[i] = allsumtemp;
            });
        }

        tp.wait_work();

        // compute first iteration: beta[i]=H[i]AllSum[i]
        Beta.clear();
        Beta.resize(nfeatures);

        for (unsigned i = 0; i < nfeatures; i++) {
            tp.push([i, &evaluator, &relin_keys, &H, &H_mutex, &AllSum, &Beta]() {
                auto allsumtemp = AllSum[i];
                auto local_hi = H[i];

                evaluator.mod_switch_to_inplace(allsumtemp, local_hi.parms_id());
                evaluator.multiply_inplace(allsumtemp, local_hi);
                evaluator.relinearize_inplace(allsumtemp, relin_keys);
                evaluator.rescale_to_next_inplace(allsumtemp);
                std::lock_guard<std::mutex> lock(H_mutex);
                Beta[i] = allsumtemp;
            });
        }

        tp.wait_work();
        // we will also need (a copy of) each data vector to be multiplied by 5/8: we
        // do this now
        dataencscale = dataenc;
        encoder.encode(0.625, scale, plain);

        for (unsigned i = 0; i < nfeatures; i++) {
            evaluator.multiply_plain_inplace(dataencscale[i], plain);
            evaluator.rescale_to_next_inplace(dataencscale[i]);
        }

        weights.clear();

        for (unsigned i = 0; i < nfeatures; i++) {
            decryptor.decrypt(Beta[i], plain);
            encoder.decode(plain, input);
            weights.emplace_back(input[0]);
        }

        // start of an iteration: we are performing the update beta[i] <- beta[i] +
        // H[i](AllSum[i] -5/8sum Beta.z(j)/2.z(ji)/2
        for (unsigned k = 1; k < parameter_set.number_of_iterations; k++) {
            // calculate the inner product: this is the only calculation where we can't
            // go feature by feature
            evaluator.mod_switch_to_inplace(dataenc[0], Beta[0].parms_id());
            evaluator.multiply(dataenc[0], Beta[0], Htemp);
            evaluator.relinearize_inplace(Htemp, relin_keys);
            evaluator.rescale_to_next_inplace(Htemp);

            for (unsigned i = 1; i < nfeatures; i++) {
                evaluator.mod_switch_to_inplace(dataenc[i], Beta[i].parms_id());
                evaluator.multiply(dataenc[i], Beta[i], ctemp);
                evaluator.relinearize_inplace(ctemp, relin_keys);
                evaluator.rescale_to_next_inplace(ctemp);
                evaluator.add_inplace(Htemp, ctemp);
            }

            // now we have a vector Htemp {beta.zi/2} i=1,...,n. so we can evaluate the
            // update circuit feature by feature
            for (unsigned i = 0; i < nfeatures; i++) {
                tp.push([&evaluator, &AllSum, &H, &H_mutex, slot_count, &gal_keys, i,
                    &dataencscale, Htemp, &relin_keys, &Beta]() {
                    Ciphertext ctemp;
                    Ciphertext allsumtemp;
                    auto d_temp = dataencscale[i]; // Reading is atomic on std::vectors, 
                                                   // but writes aren't: 
                                                   // so take a copy before operating
                    // first modify 5/8zji/2 so that it can be multiplied by the inner
                    // product
                    evaluator.mod_switch_to_inplace(d_temp, Htemp.parms_id());
                    // now multiply, relin, rescale
                    evaluator.multiply(d_temp, Htemp, ctemp);
                    evaluator.relinearize_inplace(ctemp, relin_keys);
                    evaluator.rescale_to_next_inplace(ctemp);

                    // now we need to allsum this vector:
                    for (unsigned j = 0; j < log2(slot_count); j++) {
                        allsumtemp = ctemp;
                        evaluator.rotate_vector_inplace(allsumtemp, pow(2, j), gal_keys);
                        evaluator.add_inplace(ctemp, allsumtemp);
                    }

                    // subtract this from AllSum[i], first switching down & modifying scale
                    allsumtemp = AllSum[i];
                    auto temp_h = H[i];
                    evaluator.mod_switch_to_inplace(allsumtemp, ctemp.parms_id());
                    allsumtemp.scale() = ctemp.scale();
                    evaluator.negate_inplace(ctemp);
                    evaluator.add_inplace(ctemp, allsumtemp);

                    // now multiply by H[i]
                    evaluator.mod_switch_to_inplace(temp_h, ctemp.parms_id());
                    ctemp.scale() = temp_h.scale();
                    evaluator.multiply_inplace(ctemp, temp_h);
                    evaluator.relinearize_inplace(ctemp, relin_keys);
                    evaluator.rescale_to_next_inplace(ctemp);

                    auto temp_beta = Beta[i];
                    // and finally update Beta[i]
                    evaluator.mod_switch_to_inplace(temp_beta, ctemp.parms_id());
                    temp_beta.scale() = ctemp.scale();
                    evaluator.add_inplace(temp_beta, ctemp);

                    std::lock_guard<std::mutex> lock(H_mutex);
                    H[i] = temp_h;
                    Beta[i] = temp_beta;
                    dataencscale[i] = d_temp;
            });
        }

        tp.wait_work();
        weights.clear();

        for (unsigned i = 0; i < nfeatures; i++) {
            decryptor.decrypt(Beta[i], plain);
            encoder.decode(plain, input);
            weights.emplace_back(input[0]);
        }
        cout << "iteration " << k + 1
             << " accuracy is: " << accuracy_LR(weights, cvtrain[l], 2) << "%\n";
        cout << "AUC is: " << 100 * getAUC(weights, cvtrain[l], 2) << "%\n";
        cout << "final level is: "
             << context->get_context_data(Beta[0].parms_id())->chain_index()
             << std::endl;
    }

    cout << "final level is: "
         << context->get_context_data(Beta[0].parms_id())->chain_index();

    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "fold " << l + 1 << " training done. Time = "
         << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
    cout << "CV accuracy is " << accuracy_LR(weights, cvtest[l], 2) << "%\n";
    cout << "CV AUC is " << getAUC(weights, cvtest[l], 2) << "\n";
    accuracy += accuracy_LR(weights, cvtest[l], 2);
    auc += getAUC(weights, cvtest[l], 2);
 }
 cout << "average CV accuracy is: " << accuracy / 5 << "%\n";
 cout << "average CV AUC is: " << auc / 5 << "%\n";
 return 0;
}


/**
 * Fixed Hessian Compact. Uses the Fixed Hessian Method over a compact, rather than
 * feature, encoding. This variant of the function is not really supported and it is not parallelisable.
 */

int Fixed_Hessian_Compact(bool ringdim) {
  cout << "Running Fixed Hessian with compact encoding and Chebyshev sigmoid "
          "approximation:\n";
  dMat Matrix;
  ImportDataLR(Matrix, "edin.txt", false, 2);

  EncryptionParameters parms(scheme_type::CKKS);
  size_t poly_modulus_degree = ringdim ? 32768 : 65536;

  parms.set_poly_modulus_degree(poly_modulus_degree);
  vector<int> mod;
  int x = ringdim ? 14 : 41;
  mod.push_back(50);
  for (int i = 0; i < x; i++)
    mod.push_back(40);
  mod.push_back(50);
  x = ringdim ? 4 : 12;
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
  cout << "KeyGen time = "
       << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
  CKKSEncoder encoder(context);
  size_t slot_count = encoder.slot_count();
  cout << "Number of slots: " << slot_count << "\n";
  dMatMat cvtrain, cvtest;
  CVrandomSampling(cvtrain, cvtest, Matrix);
  Matrix.clear();

  // set 40 bits of precision
  double scale = pow(2.0, 40);
  Plaintext plain, Cp;
  dVec input;
  Ciphertext ctemp, dataenc, dataencscale, H;
  Plaintext P_T1, P_T2;
  double t1, t2;
  Ciphertext Htemp, allsumtemp, allsum, beta;
  dVec weights;
  double accuracy, auc;
  accuracy = 0;
  auc = 0;
  for (int l = 0; l < 5; l++) {
    cout << "starting fold " << l + 1 << "...\n";
    int n = cvtrain[l].size();
    int nfeatures = cvtrain[l][0].size();
    int nrow = slot_count / 16;
    // pack the data into a single vector which corresponds to a 16 x n matrix,
    // padding 6 zeroes
    input.clear();
    for (unsigned i = 0; i < cvtrain[l].size(); i++) {
      for (unsigned j = 0; j < 10; j++)
        input.push_back(cvtrain[l][i][j]);
      for (unsigned j = 10; j < 16; j++)
        input.push_back(0);
    }
    cout << "Encoding...";
    start = chrono::steady_clock::now();
    encoder.encode(input, scale, plain);
    // create a matrix C with 1s in the first column and zeroes elsewhere
    dVec Cm(n * 16.0, 0);
    for (unsigned i = 0; i < n; i++)
      Cm[16.0 * i] += 1;
    // now plaintext version: this is a low precision multiplier
    encoder.encode(Cm, scale, Cp);
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encoding time = "
         << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";

    cout << "Encrypting...";
    start = chrono::steady_clock::now();
    encryptor.encrypt(plain, dataenc);
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encrypting time = "
         << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
    cout << "Creating H...";
    start = chrono::steady_clock::now();
    // creating H: first allsumrows
    H = dataenc;
    for (unsigned i = 0; i < log2(16); i++) {
      ctemp = H;
      evaluator.rotate_vector(ctemp, 16 * pow(2, i), gal_keys, ctemp);
      evaluator.add_inplace(H, ctemp);
    }
    // now we have to multiply by a matrix with 1s in the first column and 0
    // otherwise
    evaluator.multiply_plain_inplace(H, Cp);
    evaluator.rescale_to_next_inplace(H);
    // we need to replicate this value across the rows:
    for (int i = 0; i < log2(16); i++) {
      ctemp = H;
      evaluator.rotate_vector(ctemp, -pow(2, i), gal_keys, ctemp);
      evaluator.add_inplace(H, ctemp);
    }
    // multiply with the data again
    evaluator.mod_switch_to(dataenc, H.parms_id(), ctemp);
    evaluator.multiply_inplace(H, dataenc);
    // allsum H
    for (int j = 0; j < log2(nrow); j++) {
      ctemp = H;
      evaluator.rotate_vector(ctemp, pow(2, j), gal_keys, ctemp);
      evaluator.add_inplace(H, ctemp);
    }
    // time to calculate 1/H -- first we need our starting point, T1 + T2D
    cout << "Calculating 1/H(i)...";
    t1 = T1(0.25 * (nfeatures * n));
    t2 = T2(0.25 * (nfeatures * n));
    encoder.encode(t1, scale, P_T1);
    encoder.encode(t2, scale, P_T2);

    // T2 needs to be multiplied by something 2 levels down, so:
    evaluator.mod_switch_to_next_inplace(P_T2);
    evaluator.mod_switch_to_next_inplace(P_T2);

    // T1 will need to be added to something 3 levels down, so:
    evaluator.mod_switch_to_next_inplace(P_T1);
    evaluator.mod_switch_to_next_inplace(P_T1);
    evaluator.mod_switch_to_next_inplace(P_T1);

    for (int i = 0; i < nfeatures; i++) {
      // negate and store a copy of H(i,i) for later:
      ctemp = H;
      evaluator.negate_inplace(ctemp);

      // find v0 = T1 + T2H(i,i):
      evaluator.multiply_plain_inplace(H, P_T2);
      evaluator.rescale_to_next_inplace(H);
      P_T1.scale() = H.scale();
      evaluator.add_plain_inplace(H, P_T1);

      // now iterate Newton Raphson: each update is 2v - Hv^2
      for (int j = 0; j < 3; j++) {

        // first double and store the result
        evaluator.add(H, H, Htemp);

        // now square the current value, relin and rescale
        evaluator.square_inplace(H);
        evaluator.relinearize_inplace(H, relin_keys);
        evaluator.rescale_to_next_inplace(H);

        // now mod switch down our stored value ctemp, and multiply
        evaluator.mod_switch_to_inplace(Htemp, H.parms_id());
        evaluator.multiply_inplace(H, ctemp);
        evaluator.relinearize_inplace(H, relin_keys);
        evaluator.rescale_to_next_inplace(H);

        // modify scale of Htemp = 2v, mod switch down, and add
        Htemp.scale() = H.scale();
        evaluator.mod_switch_to_inplace(Htemp, H.parms_id());
        evaluator.add_inplace(H, Htemp);
      }
    }
    cout << "computing allsum...\n";
    // precompute AllSum: the ith column of allsum is constant and = sum(zji/2)
    allsum = dataenc;
    for (int j = 0; j < log2(nrow); j++) {
      ctemp = allsum;
      evaluator.rotate_vector_inplace(ctemp, pow(2, j), gal_keys);
      evaluator.add_inplace(allsum, ctemp);
    }

    // compute first iteration: beta=HxAllSum
    beta = allsum;
    evaluator.mod_switch_to_inplace(beta, H.parms_id());
    evaluator.multiply_inplace(beta, H);
    evaluator.relinearize_inplace(beta, relin_keys);
    evaluator.rescale_to_next_inplace(beta);
    // we will also need (a copy of) each data vector to be multiplied by 5/8:
    // we do this now
    dataencscale = dataenc;
    encoder.encode(0.625, scale, plain);
    for (int i = 0; i < nfeatures; i++) {
      evaluator.multiply_plain_inplace(dataencscale, plain);
      evaluator.rescale_to_next_inplace(dataencscale);
    }
    weights.clear();
    decryptor.decrypt(beta, plain);
    encoder.decode(plain, input);
    for (int i = 0; i < nfeatures; i++)
      weights.push_back(input[i]);
    cout << weights.size();
    cout << "," << cvtrain[l][0].size() << "\n";
    cout << "1st iteration accuracy: " << accuracy_LR(weights, cvtrain[l], 2)
         << "%\n";
    cout << "1st iteration AUC: " << 100 * getAUC(weights, cvtrain[l], 2)
         << "%\n";
    // start of an iteration: we are performing the update beta[i] <- beta[i] +
    // H[i](AllSum[i] -5/8sum Beta.z(j)/2.z(ji)/2

    for (int k = 2; k < x; k++) {
      cout << "creating inner product \n";
      // calculate the inner product: first multiply the data by beta
      evaluator.mod_switch_to_inplace(dataenc, beta.parms_id());
      evaluator.multiply(dataenc, beta, Htemp);
      evaluator.relinearize_inplace(Htemp, relin_keys);
      evaluator.rescale_to_next_inplace(Htemp);
      // now all sum rows
      for (int i = 0; i < log2(16); i++) {
        ctemp = Htemp;
        evaluator.rotate_vector(ctemp, pow(2, i), gal_keys, ctemp);
        evaluator.add_inplace(Htemp, ctemp);
      }
      // multiply by C
      evaluator.mod_switch_to_inplace(Cp, Htemp.parms_id());
      evaluator.multiply_plain_inplace(Htemp, Cp);
      evaluator.rescale_to_next_inplace(Htemp);
      // replicate across rows
      for (int i = 0; i < log2(16); i++) {
        ctemp = Htemp;
        evaluator.rotate_vector(ctemp, -pow(2, i), gal_keys, ctemp);
        evaluator.add_inplace(Htemp, ctemp);
      }
      // now we have a vector Htemp {beta.zi/2} i=1,...,n, and so can evaluate
      // polynomial

      // first modify 5/8zji/2 so that it can be multiplied by the inner product
      evaluator.mod_switch_to_inplace(dataencscale, Htemp.parms_id());
      // now multiply, relin, rescale
      evaluator.multiply(dataencscale, Htemp, ctemp);
      evaluator.relinearize_inplace(ctemp, relin_keys);
      evaluator.rescale_to_next_inplace(ctemp);

      // now we need to allsum columns:
      for (int j = 0; j < log2(nrow); j++) {
        allsumtemp = ctemp;
        evaluator.rotate_vector_inplace(allsumtemp, 16 * pow(2, j), gal_keys);
        evaluator.add_inplace(ctemp, allsumtemp);
      }

      // subtract this from AllSum, first switching down & modifying scale
      allsumtemp = allsum;
      evaluator.mod_switch_to_inplace(allsumtemp, ctemp.parms_id());
      allsumtemp.scale() = ctemp.scale();
      evaluator.negate_inplace(ctemp);
      evaluator.add_inplace(ctemp, allsumtemp);

      // now multiply by H[i]
      evaluator.mod_switch_to_inplace(H, ctemp.parms_id());
      evaluator.multiply_inplace(ctemp, H);
      evaluator.relinearize_inplace(ctemp, relin_keys);
      evaluator.rescale_to_next_inplace(ctemp);

      // and finally update Beta[i]
      evaluator.mod_switch_to_inplace(beta, ctemp.parms_id());
      beta.scale() = ctemp.scale();
      evaluator.add_inplace(beta, ctemp);

      weights.clear();
      decryptor.decrypt(beta, plain);
      encoder.decode(plain, input);
      for (int i = 0; i < nfeatures; i++)
        weights.push_back(input[i]);
      cout << "iteration " << k
           << " accuracy is: " << accuracy_LR(weights, cvtrain[l], 2) << "%\n";
      cout << "AUC is: " << 100 * getAUC(weights, cvtrain[l], 2) << "%\n";
    }
    cout << "final level is: "
         << context->get_context_data(beta.parms_id())->chain_index();

    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "fold " << l + 1 << " training done. Time = "
         << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
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
  ImportDataLR(Matrix, "edin.txt", false, 2);
  int n = Matrix.size();
  int nfeatures = Matrix[0].size();

  EncryptionParameters parms(scheme_type::CKKS);
  size_t poly_modulus_degree = 32768;

  parms.set_poly_modulus_degree(poly_modulus_degree);

  parms.set_coeff_modulus(CoeffModulus::Create(
      poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                            40, 40, 40, 40, 40, 40, 40, 40, 40, 50}));
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
  cout << "KeyGen time = "
       << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
  CKKSEncoder encoder(context);
  size_t slot_count = encoder.slot_count();
  cout << "Number of slots: " << slot_count << "\n";
  cout << "Encoding...";
  start = chrono::steady_clock::now();

  // set 40 bits of precision
  double scale = pow(2.0, 40);

  // encode the data, a column at a time. The ith element of data corresponds to
  // the ith feature
  pVec dataplain;
  dataplain.reserve(nfeatures);
  Plaintext plain;
  dVec input;
  for (int i = 0; i < nfeatures; i++) {
    input.clear();
    for (int j = 0; j < Matrix.size(); j++)
      input.push_back(Matrix[j][i]);
    encoder.encode(input, scale, plain);
    dataplain.push_back(plain);
  }

  end = chrono::steady_clock::now();
  diff = end - start;
  cout << "Encoding time = "
       << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";

  cout << "Encrypting...";
  start = chrono::steady_clock::now();
  // encrypt these plaintexts
  cVec dataenc;
  dataenc.reserve(nfeatures);
  Ciphertext ctemp;
  for (int i = 0; i < dataplain.size(); i++) {
    encryptor.encrypt(dataplain[i], ctemp);
    dataenc.push_back(ctemp);
  }
  end = chrono::steady_clock::now();
  diff = end - start;
  cout << "Encrypting time = "
       << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";
  cout << "Creating H...";
  start = chrono::steady_clock::now();
  // creating H: first add all ciphertexts
  cVec H;

  H.reserve(nfeatures);

  for (int i = 1; i < nfeatures; i++)
    evaluator.add_inplace(ctemp, dataenc[1. * nfeatures - 1 - i]);

  // now create H(i,i)
  Ciphertext Htemp, allsumtemp;
  for (int i = 0; i < nfeatures; i++) {
    evaluator.multiply(dataenc[i], ctemp, Htemp);
    evaluator.relinearize_inplace(Htemp, relin_keys);
    evaluator.rescale_to_next_inplace(Htemp);
    // allsum Htemp
    allsumtemp = Htemp;
    for (int j = 0; j < log2(slot_count); j++) {
      Htemp = allsumtemp;
      evaluator.rotate_vector(Htemp, pow(2, j), gal_keys, Htemp);
      evaluator.add_inplace(allsumtemp, Htemp);
    }
    // push back H(i,i) to the vector H
    H.push_back(allsumtemp);
  }
  // time to calculate 1/H(i,i) -- first we need our starting point, T1 + T2D
  cout << "Calculating 1/H(i)...";

  Plaintext P_T1, P_T2;
  double t1, t2;
  t1 = T1(2.5 * n);
  t2 = T2(2.5 * n);
  encoder.encode(t1, scale, P_T1);
  encoder.encode(t2, scale, P_T2);

  // T2 needs to be multiplied by something 1 level down, so:
  evaluator.mod_switch_to_next_inplace(P_T2);

  // T1 will need to be added to something 2 levels down, so:
  evaluator.mod_switch_to_next_inplace(P_T1);
  evaluator.mod_switch_to_next_inplace(P_T1);

  for (int i = 0; i < nfeatures; i++) {

    // negate and store a copy of H(i,i) for later:
    Htemp = H[i];
    evaluator.negate_inplace(Htemp);

    // find v0 = T1 + T2H(i,i):
    evaluator.multiply_plain_inplace(H[i], P_T2);
    evaluator.rescale_to_next_inplace(H[i]);
    P_T1.scale() = H[i].scale();
    evaluator.add_plain_inplace(H[i], P_T1);

    // now iterate Newton Raphson: each update is 2v - H[i,i]v^2
    for (int j = 0; j < 3; j++) {

      // first double and store the result
      evaluator.add(H[i], H[i], ctemp);

      // now square the current value, relin and rescale

      evaluator.square_inplace(H[i]);
      evaluator.relinearize_inplace(H[i], relin_keys);
      evaluator.rescale_to_next_inplace(H[i]);

      // now mod switch down our stored value Htemp, and multiply
      evaluator.mod_switch_to_inplace(Htemp, H[i].parms_id());
      evaluator.multiply_inplace(H[i], Htemp);
      evaluator.relinearize_inplace(H[i], relin_keys);
      evaluator.rescale_to_next_inplace(H[i]);

      // modify scale of ctemp = 2v, mod switch down, and add
      ctemp.scale() = H[i].scale();
      evaluator.mod_switch_to_inplace(ctemp, H[i].parms_id());
      evaluator.add_inplace(H[i], ctemp);
    }
  }
  // precompute AllSum: AllSum[i] = sum(zji/2)
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

  // compute first iteration: beta[i]=H[i]AllSum[i]
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

  // for the Taylor polynomial, each update only uses AllSum multiplied by H(i)
  AllSum = Beta;

  // we will also need (a copy of) each data vector to be multiplied by H(i): we
  // do this now
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
    for (int j = 0; j < input.size(); j++)
      weights[i] += input[j];
    weights[i] /= input.size();
  }
  Matrix.clear();
  ImportDataLR(Matrix, "edin.txt", false);
  cout << "accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
  cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
  // start of an iteration: we are performing the update beta[i] <- beta[i] +
  // H[i]AllSum[i] -sum Beta.z(j)/2.z(ji)/2

  for (int k = 2; k < 5; k++) {

    // calculate the inner product: this is the only calculation where we can't
    // go feature by feature
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

    // now we have a vector Htemp {beta.zi/2} i=1,...,n. so we can evaluate the
    // update circuit feature by feature
    for (int i = 0; i < nfeatures; i++) {

      // first modify H(i)zji/2 so that it can be multiplied by the inner
      // product
      evaluator.mod_switch_to_inplace(dataencscale[i], Htemp.parms_id());
      // now multiply, relin, rescale
      evaluator.multiply(dataencscale[i], Htemp, ctemp);
      evaluator.relinearize_inplace(ctemp, relin_keys);
      evaluator.rescale_to_next_inplace(ctemp);

      // now we need to allsum this vector:
      for (int j = 0; j < log2(slot_count); j++) {
        allsumtemp = ctemp;
        evaluator.rotate_vector_inplace(allsumtemp, pow(2, j), gal_keys);
        evaluator.add_inplace(ctemp, allsumtemp);
      }

      // subtract this from AllSum[i], first switching down & modifying scale
      allsumtemp = AllSum[i];
      evaluator.mod_switch_to_inplace(allsumtemp, ctemp.parms_id());
      allsumtemp.scale() = ctemp.scale();
      evaluator.negate_inplace(ctemp);
      evaluator.add_inplace(ctemp, allsumtemp);

      // and finally update Beta[i]
      evaluator.mod_switch_to_inplace(Beta[i], ctemp.parms_id());
      Beta[i].scale() = ctemp.scale();
      evaluator.add_inplace(Beta[i], ctemp);
    }
    for (int i = 0; i < nfeatures; i++) {
      decryptor.decrypt(Beta[i], plain);
      encoder.decode(plain, input);
      for (int j = 0; j < input.size(); j++)
        weights[i] += input[j];
      weights[i] /= input.size();
    }
    Matrix.clear();
    ImportDataLR(Matrix, "edin.txt", false);
    cout << "iteration " << k
         << " accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
    cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
  }
  cout << "final level is: "
       << context->get_context_data(Beta[0].parms_id())->chain_index();
  end = chrono::steady_clock::now();
  diff = end - start;
  cout << "done. Time = "
       << chrono::duration<double, milli>(diff).count() / 1000.0 << " s \n";

  for (int i = 0; i < nfeatures; i++) {
    decryptor.decrypt(Beta[i], plain);
    encoder.decode(plain, input);
    for (int j = 0; j < input.size(); j++)
      weights[i] += input[j];
    weights[i] /= input.size();
  }
  Matrix.clear();
  ImportDataLR(Matrix, "edin.txt", false);
  cout << "accuracy is: " << accuracy_LR(weights, Matrix) << "%\n";
  cout << "AUC is: " << 100 * getAUC(weights, Matrix) << "%\n";
  return 0;
}
