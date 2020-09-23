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
    cout << "| 2. Gradient Descent        | GD.cpp                     |" << endl;
    cout << "| 3. Nesterov GD             | NAG.cpp                    |" << endl;
    cout << "| 4. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Chebyshev approximation |                            |" << endl;
    cout << "| 5. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Taylor approximation    |                            |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    int selection = 0;
    string iternum;
    bool invalid = true;
    do
    {
        cout << endl << "> Run example (1 ~ 5) or exit (0): ";
        if (!(cin >> selection))
        {
            invalid = false;
        }
        else if (selection < 0 || selection > 5)
        {
            invalid = false;
        }
        else
        {
            invalid = true;
        }
        if (!invalid)
        {
            cout << "  [Beep~~] Invalid option: type 0 ~ 5" << endl;
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
        GD();
        break;
    case 3:
        Nesterov_GD();
        break;

    case 4:
        Fixed_Hessian_Chebyshev();
        break;

    case 5:
        Fixed_Hessian_Taylor();
        break;

    case 0:
        return 0;
    }
    return 0;
}
