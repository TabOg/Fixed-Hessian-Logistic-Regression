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
    bool ringbool{};
    string encoding;
    do
    {
        cout << endl << "> Run (1 ~ 6) or exit (0): ";
        if (!(cin >> selection))
        {
            invalid = false;
        }
        else if (selection < 0 || selection > 6)
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
    case 3:
        cout << "Select N = 2^16 (0) or N = 2^15 (1):";
        cin >> ringdim;
        while (!is_number(ringdim)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        while (!(stoi(ringdim) == 0 || stoi(ringdim) == 1)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        ringbool = (stoi(ringdim) == 1);
        GD(ringbool);
        break;
    case 4:
        cout << "Select N = 2^16 (0) or N = 2^15 (1):";
        cin >> ringdim;
        while (!is_number(ringdim)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        while (!(stoi(ringdim) == 0 || stoi(ringdim) == 1)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        ringbool = (stoi(ringdim) == 1);

        Nesterov_GD(ringbool);
        break;

    case 5:
        cout << "Select N = 2^16 (0) or N = 2^15 (1):";
        cin >> ringdim;
        while (!is_number(ringdim)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        while (!(stoi(ringdim) == 0 || stoi(ringdim) == 1)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select N = 2^16 (0) or N = 2^15 (1):";
            cin >> ringdim;
        }
        
        ringbool = (stoi(ringdim) == 1);
        cout << "Select Encoding Style: Feature (1) or Compact (0):";
        cin >> encoding;
        while (!is_number(encoding)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select Encoding Style: Feature (1) or Compact (0):";
            cin >> encoding;
        }
        /*
        while (!(stoi(encoding) == 0 || stoi(encoding) == 1)) {
            cout << "Please enter either 0 or 1!\n";
            cout << "Select Encoding Style: Feature (1) or Compact (0):";
            cin >> encoding;
        }
        */
        

        //if (stoi(encoding) == 1)Fixed_Hessian_Chebyshev(ringbool);
        //else Fixed_Hessian_Compact(ringbool);
        Fixed_Hessian_Chebyshev(ringbool);
        break;

    case 6:
        Fixed_Hessian_Taylor();
        break;

    case 0:
        return 0;
    }
    return 0;
}
