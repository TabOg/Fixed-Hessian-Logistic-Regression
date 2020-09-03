#include "seal/seal.h"
#include "databasetools.h"
#include "logregtools.h"
#include <iostream>


using namespace std;
using namespace seal;

int GD() {

    
    dMat Matrix;
    ImportDataLR(Matrix, "edin.txt", false, 8);
    int n = Matrix.size();

    int nfeatures = Matrix[0].size();
    EncryptionParameters parms(scheme_type::CKKS);
    size_t poly_modulus_degree = 32768;
    vector<int> mod;
    mod.push_back(38);
    for (int i = 0; i < 26; i++)mod.push_back(28);
    mod.push_back(38);
    parms.set_poly_modulus_degree(poly_modulus_degree);
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
    cout << "Encoding...";
    start = chrono::steady_clock::now();
    double scale = pow(2.0, 28);
    pVec data;
    Plaintext plain;
    dVec input;
    for (int i = 0; i < nfeatures; i++) {
        input.clear();
        for (int j = 0; j < Matrix.size(); j++)input.push_back(Matrix[j][i]);
        encoder.encode(input, scale, plain);
        data.push_back(plain);
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encoding time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    cVec dataenc;
    start = chrono::steady_clock::now();
    cout << "Encrypting...";
    Ciphertext datatemp;
    for (int i = 0; i < nfeatures; i++) {
        encryptor.encrypt(data[i], datatemp);
        dataenc.push_back(datatemp);
    }
    end = chrono::steady_clock::now();
    diff = end - start;
    cout << "Encrypting time = " << chrono::duration <double, milli>(diff).count() / 1000.0 << " s \n";

    double a0 = 0.5;
    double a1 = -1.20096;
    double a3 = 0.81562;
    double a4 = a1 / a3;
    double sc = 4.0 / (1.0 * n);

    cout << "encoding polynomial coefficients...";
    Plaintext coeff0, coeff1, coeff3, coeff4;
    encoder.encode(a0, scale, coeff0);
    encoder.encode(a1, scale, coeff1);
    encoder.encode(2 * a3, scale, coeff3);
    encoder.encode(a4, scale, coeff4);
    cout << "done \n";
    auto begin = chrono::steady_clock::now();
    cout << "calculating ct.sum...\n";
    start = chrono::steady_clock::now();
    Plaintext scaler, scaler1;
    cVec ctsum;
    //creating the vector ctsum, ctsum[j] = 4/nsum(z_{ij}/8)
    encoder.encode(sc, scale, scaler);
    encoder.encode(5 * sc, scale, scaler1);
    Ciphertext allsum, temp;
    cVec Beta;
    for (int i = 0; i < nfeatures; i++) {
        allsum = dataenc[i];
        for (int i = 0; i < log2(slot_count); i++) {
            temp = allsum;
            evaluator.rotate_vector(temp, pow(2, i), gal_keys, temp);
            evaluator.add_inplace(allsum, temp);
        }
        temp = allsum;
        evaluator.multiply_plain_inplace(allsum, scaler);
        evaluator.multiply_plain_inplace(temp, scaler1);
        evaluator.rescale_to_next_inplace(allsum);
        evaluator.rescale_to_next_inplace(temp);
        ctsum.push_back(allsum);
        Beta.push_back(temp);
    }
    dVec weights;
    for (int i = 0; i < nfeatures; i++) {
        decryptor.decrypt(Beta[i], plain);
        encoder.decode(plain, input);
        weights.push_back(accumulate(input.begin(), input.end(), 0.0) / (1.0 * input.size()));
    }
    Matrix.clear();
    ImportDataLR(Matrix, "edin.txt", false, 1.0, '\t');
    cout << "Modulus chain index for Beta " << context->get_context_data(Beta[0].parms_id())->chain_index() << endl;
    cout << "Beta scale" << log2(Beta[0].scale()) << "\n";
    cout << "1st iteration AUC is " << 100 * getAUC(weights, Matrix) << "%\n";
    cout << "1st accuracy is " << accuracy_LR(weights, Matrix) << "%";
    cVec dataencscale = dataenc;
    for (int i = 0; i < dataenc.size(); i++) {
        evaluator.multiply_plain_inplace(dataencscale[i], coeff3);
        evaluator.rescale_to_next_inplace(dataencscale[i]);
    }
    end = chrono::steady_clock::now();
    diff = end - start;


    cout << "done. Time =" << chrono::duration <double, milli>(diff).count() / 1000.0 << "s. \n";
    Plaintext alphap, p;
    Ciphertext mult, innerprod, square;
    double alpha;
    //iterations
    for (int k = 2; k < 8; k++) {
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
        for (int i = 1; i < nfeatures; i++) {
            evaluator.mod_switch_to(dataenc[i], Beta[i].parms_id(), mult);
            evaluator.multiply_inplace(mult, Beta[i]);
            evaluator.relinearize_inplace(mult, relin_keys);
            evaluator.rescale_to_next_inplace(mult);
            evaluator.add_inplace(innerprod, mult);
        }
        //time to evaluate the polynomial!
        square = innerprod;
        evaluator.square_inplace(square);
        evaluator.relinearize_inplace(square, relin_keys);
        evaluator.rescale_to_next_inplace(square);
        evaluator.mod_switch_to_inplace(coeff4, square.parms_id());
        coeff4.scale() = square.scale();
        evaluator.add_plain_inplace(square, coeff4);
        //we now have wi^2 + a1/a3 -- calculate the value of the polynomial for each j
        for (int j = 0; j < nfeatures; j++) {
            evaluator.mod_switch_to_inplace(dataencscale[j], innerprod.parms_id());
            evaluator.multiply(innerprod, dataencscale[j], mult);
            evaluator.relinearize_inplace(mult, relin_keys);
            evaluator.rescale_to_next_inplace(mult);
            /*cout << context->get_context_data(square.parms_id())->chain_index() << "\n";
            cout << context->get_context_data(squarea4.parms_id())->chain_index() << "\n";*/
            evaluator.multiply_inplace(mult, square);
            evaluator.relinearize_inplace(mult, relin_keys);
            evaluator.rescale_to_next_inplace(mult);
            evaluator.mod_switch_to_inplace(scaler, mult.parms_id());
            evaluator.multiply_plain_inplace(mult, scaler);
            evaluator.rescale_to_next_inplace(mult);

            //time to do AllSum on mult
            allsum = mult;
            for (int i = 0; i < log2(slot_count); i++) {
                temp = mult;
                evaluator.rotate_vector_inplace(temp, pow(2, i), gal_keys);
                evaluator.add_inplace(mult, temp);
            }
            //switch down a copy of ctsum[i], add to all sum, update poly
            evaluator.multiply_plain(ctsum[j], alphap, temp);
            evaluator.rescale_to_next_inplace(temp);
            temp.scale() = mult.scale();
            evaluator.mod_switch_to_inplace(temp, mult.parms_id());
            evaluator.add_inplace(mult, temp);
            //update weights vector
            mult.scale() = Beta[j].scale();
            evaluator.mod_switch_to_inplace(Beta[j], mult.parms_id());
            evaluator.add_inplace(Beta[j], mult);
        }
        cout <<"Modulus chain index for Beta "<< context->get_context_data(Beta[0].parms_id())->chain_index() << endl;
        cout << "Beta scale" << log2(Beta[0].scale()) << "\n";
        weights.clear();
        for (int i = 0; i < nfeatures; i++) {
            decryptor.decrypt(Beta[i], plain);
            encoder.decode(plain, input);
            weights.push_back(input[0]);
        }
        Matrix.clear();
        ImportDataLR(Matrix, "edin.txt", false, 1.0, '\t');
        cout << "iteration " << k << " AUC is " << 100 * getAUC(weights, Matrix) << "%\n";
        cout << "iteration " << k << " accuracy is " << accuracy_LR(weights, Matrix) << "%";
    }    
    auto ending = chrono::steady_clock::now();
    auto total = end - begin;
    cout << "done. Time =" << chrono::duration <double, milli>(total).count() / 1000.0 << "s. \n";
    return 0;
}