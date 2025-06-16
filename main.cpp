#include <iostream>
#include "Pyramide.h"
#include "tests.hpp"

using namespace std;

int main() {

    // int ordre;
    // std::cout<< "ordre : ";
    // std::cin>>ordre;

    // Pyramide a = Pyramide(ordre);

    // auto M = a.calculerMatriceMasse();
    // cout << "\n\n\nMatrice de masse : \n" << endl;
    // cout << M << endl;

    // cout << "\n\n\nMatrice de masse transposée : \n" << endl;
    // cout << M.transpose() << endl;

    test1();
    test2();
    test3();
    cout << "Tests terminés." << endl;



	return 0;
}