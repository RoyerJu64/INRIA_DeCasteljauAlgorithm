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
    cout<< "\n\nEntrez l'ordre de la pyramide (1 ou 2) pour afficher la matrice de masse et la matrice de raideur : ";
    cin >> ordre;
    Pyramide a = Pyramide(ordre);

    auto start1 = std::chrono::high_resolution_clock::now();
    auto M = a.calculerMatriceMasse();
    cout << "\n\nMatrice de masse : \n" << endl;
    std::cout << M << std::endl;

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start1;
    std::cout << "Temps d'exécution pour le calcul de M : " << elapsed1.count() << " s\n";

    auto start2 = std::chrono::high_resolution_clock::now();
    auto K = a.calculerMatriceRaideur();
    cout << "\n\nMatrice de raideur : \n" << endl;
    std::cout << K << std::endl;
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = end2 - start2;
    std::cout << "Temps d'exécution pour le calcul de K : " << elapsed2.count() << " s\n";

    std::chrono::duration<double> elapsedtot = end2 - start1;
    std::cout << "Temps d'exécution total : " << elapsedtot.count() << " s\n";

	return 0;
}