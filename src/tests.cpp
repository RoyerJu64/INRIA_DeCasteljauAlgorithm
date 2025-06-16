#include "tests.hpp"

void test1(){

    for(int i = 1; i<3; ++i){
        Pyramide a = Pyramide(i);
        auto M = a.calculerMatriceMasse();


    if (!M.isApprox(M.transpose())) {
        std::cout << "Erreur dans le test 1 pour ordre " << i << std::endl;
} else {
        std::cout << "Test 1 réussi pour ordre " << i << std::endl;
}
    }
}

void test2(){
    for(int i = 1; i<3; ++i){
        Pyramide a = Pyramide(i);
        int N = a.calculerMatriceMasse().rows();  // ou .cols()
        Eigen::VectorXd un = Eigen::VectorXd::Ones(N);
        double V = un.transpose() * a.calculerMatriceMasse() * un;

        std::cout << "Ordre " << i << " : V = " << V << ", attendu = " << 1.0/3.0 << std::endl;
        if(std::abs(V - 1.0/3.0) > 1e-10){
            std::cout << "Erreur dans le test 2 pour ordre " << i << std::endl;
        } else {
            std::cout << "Test 2 réussi pour ordre " << i << std::endl;
        }
    }
}


void test3(){
    for(int i = 1; i<3; ++i){
        Pyramide a = Pyramide(i);
        std::random_device rd;  // source d'entropie
        std::mt19937 gen(rd()); // générateur de nombres pseudo-aléatoires
        std::uniform_real_distribution<> dis(0.0, 1.0);

        double u = dis(gen);
        double v_ = dis(gen);
        double w = dis(gen);
        // If your baseBernsteinPyramide expects u+v+w <= 1, normalize:
        if (u + v_ + w > 1.0) {
            double sum = u + v_ + w;
            u /= sum;
            v_ /= sum;
            w /= sum;
        }
        std::vector<double> v = a.baseBernsteinPyramide(u, v_, w);
        double somme = std::accumulate(v.begin(), v.end(), 0.0);
        if(std::abs(somme - 1.0) > 1e-10){
            std::cout << "Erreur dans le test 3 pour ordre " << i << std::endl;
        } else {
            std::cout << "Test 3 réussi pour ordre " << i << std::endl;
        }
    }
}

