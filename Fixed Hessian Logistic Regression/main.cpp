#include "seal/seal.h"
#include <iostream>
#include "databasetools.h"
#include "logregtools.h"
#include "main.h"

using namespace std;
using namespace seal;

int main() {
    
    cout << "+---------------------------------------------------------+" << endl;
    cout << "| Select a procedure to run:                              |" << endl;
    cout << "+---------------------------------------------------------+" << endl;
    cout << "| Methods                    | Source Files               |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    cout << "| 1. Plaintext Gradient      | plaintextLR.cpp            |" << endl;
    cout << "|    Descent                 |                            |" << endl;
    cout << "| 2. Plaintext Nesterov      | plaintextLR.cpp            |" << endl;
    cout << "|    Descent                 |                            |" << endl;
    cout << "| 3. Gradient Descent        | GD.cpp                     |" << endl;
    cout << "| 4. Nesterov GD             | NAG.cpp                    |" << endl;
    cout << "| 5. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Chebyshev approximation |                            |" << endl;
    cout << "| 6. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Taylor approximation    |                            |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    int selection = 0;
    string iternum;
    string ringdim;
    bool invalid = true;
    string encoding;
    do
    {
        cout << endl << "> Run (1 ~ 6) or exit (0): ";
        if (!(cin >> selection))
        {
            invalid = false;
        }
        else if (selection < 0 || selection > 7)
        {
            invalid = false;
        }
        else
        {
            invalid = true;
        }
        if (!invalid)
        {
            cout << "  [Beep~~] Invalid option: type 0 ~ 6" << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    } while (!invalid);

    switch (selection)
    {
    case 1:
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        Plaintext_LR("edin.txt", stoi(iternum));
        break;
    case 2:
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        Plaintext_LR_NAG("edin.txt", stoi(iternum));
        break;
    case 3: {
            cout << "Select the precision you would like to run at:";
            std::string precision;
            cin >> precision;
            while (!is_number(precision)) {
                cout << "Please enter a whole number\n";
                cout << "Select the precision you would like to run at:";
                cin >> precision;
            }
        
            cout << "Enter the number of threads you would like to run on (0 for all threads):";
            std::string threads;
            cin >> threads;
            while(!is_number(threads)) {
                cout << "Please enter a whole number\n";
                cout << "Enter the number of threads you would like to run on";
                cin  >> threads;
            }

            const unsigned int nr_precision =  stoi(precision);
            const unsigned int nr_threads = stoi(threads);
   
            GD(nr_precision, nr_threads);
        }
        break;
    case 4: {
                cout << "Select the precision you would like to run at:";
                string precision;
                cin >> precision;
                while (!is_number(precision)) {
                    cout << "Please enter a whole number\n";
                    cout << "Select the precision you would like to run at:";
                    cin >> precision;
                }
            
                const unsigned int nr_precision = stoi(precision);
                Nesterov_GD(nr_precision);
            }

        break;

    case 5: {
        cout << "Select the precision you would like to run at:";
        std::string precision;
        cin >> precision;
        while (!is_number(precision)) {
            cout << "Please enter a whole number\n";
            cout << "Select the precision you would like to run at:";
            cin >> precision;
        }
        
        cout << "Enter the number of threads you would like to run on (0 for all threads):";
        std::string threads;
        cin >> threads;
        while(!is_number(threads)) {
            cout << "Please enter a whole number\n";
            cout << "Enter the number of threads you would like to run on";
            cin  >> threads;
        }

        const unsigned int nr_precision =  stoi(precision);
        const unsigned int nr_threads = stoi(threads);
        Fixed_Hessian_Chebyshev(nr_precision, nr_threads);
        }
    
        break;
    case 6:
        Fixed_Hessian_Taylor();
        break;
    case 7:
        cout << "number of iterations = ";
        cin >> iternum;
        while (!is_number(iternum)) {
            cout << "iteration number must be an integer!\n";
            cout << "number of iterations = ";
            cin >> iternum;
        }
        Plaintext_LR_lowdeg("edin.txt", stoi(iternum));
        break;

    case 0:
        return 0;
    }
    return 0;
}
