#include <iostream>
#include "Pyramide.h"
#include "tests.hpp"

using namespace std;

int main() {


    test1();
    test2();
    test3();
    test4();
    cout << "Tests terminés." << endl;

    int ordre;
    std::cout<< "\n\nOrdre pour calculer la matrice de masse et la matrice de raideur : ";
    std::cin>>ordre;

    Pyramide a = Pyramide(ordre);

    auto M = a.calculerMatriceMasse();
    cout << "\n\n\nMatrice de masse : \n" << endl;
    cout << M << endl;

    auto K = a.calculerMatriceRaideur();
    cout << "\n\n\nMatrice de raideur : \n" << endl;
    cout << K << endl;


	return 0;
}