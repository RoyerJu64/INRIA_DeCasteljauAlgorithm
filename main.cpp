#include <iostream>
#include "Pyramide.h"

using namespace std;

int main() {

    int ordre;
    std::cout<< "ordre : ";
    std::cin>>ordre;

    Pyramide a = Pyramide(ordre);

//     std::cout << "Noeuds pour ordre " << ordre << " :" << std::endl;
//     for (const auto& noeuds : a.getNoeuds()) {
//         double i = std::get<0>(noeuds);
//         double j = std::get<1>(noeuds);
//         double k = std::get<2>(noeuds);
//         std::cout << "(" << i << ", " << j << ", " << k << ")" << std::endl;
//     }

//     std::cout << "\n\nTriplets pour ordre " << ordre << " :" << std::endl;
//     for (const auto& triplets : a.getTriplets()) {
//         int i = std::get<0>(triplets);
//         int j = std::get<1>(triplets);
//         int k = std::get<2>(triplets);
//         std::cout << "(" << i << ", " << j << ", " << k << ")" << std::endl;
//     }

// std::cout << "\n\nTriplets pour ordre " << ordre << " :" << std::endl;
//     for (const auto& paire : a.getTripletToIndex()) {
//         auto& t = paire.first;
//         int i = std::get<0>(t);
//         int j = std::get<1>(t);
//         int k = std::get<2>(t);
//         std::cout << "("<<i<<j<<k<<")"<<"  "<<paire.second << std::endl;
//     }

    // double x = a.compute_Mijklmn(0, 0, 0, 0, 0, 0, 1, 4);
    // cout << "\n\n" << x << endl;

    auto M = a.calculerMatriceMasse();
    cout << "\n\n\nMatrice de masse : \n" << endl;
    cout << M << endl;


	return 0;
}