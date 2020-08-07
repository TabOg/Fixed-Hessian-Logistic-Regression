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
    cout << "| 1. Nesterov GD             | NAG.cpp                    |" << endl;
    cout << "| 2. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Chebyshev approximation |                            |" << endl;
    cout << "| 3. Fixed Hessian with      | Fixed Hessian.cpp          |" << endl;
    cout << "|    Taylor approximation    |                            |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    int selection = 0;
    bool invalid = true;
    do
    {
        cout << endl << "> Run example (1 ~ 3) or exit (0): ";
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
            cout << "  [Beep~~] Invalid option: type 0 ~ 3" << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    } while (!invalid);

    switch (selection)
    {
    case 1:
        Nesterov_GD();
        break;

    case 2:
        Fixed_Hessian_Chebyshev();
        break;

    case 3:
        Fixed_Hessian_Taylor();
        break;

    case 0:
        return 0;
    }
    return 0;
}
