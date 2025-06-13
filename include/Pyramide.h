#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <map>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <Eigen/Dense>

class Pyramide
{
public:
    Pyramide(int ordre);


    void genererNoeuds();
    void genererTriplets();
    void genererTripletToIndex();

    double bernstein(int i, int n, double t);
    std::vector<double> baseBernsteinCube(double x, double y, double z);
    std::vector<double> baseBernsteinPyramide(double x, double y, double z);

    double J(double a, double b);

    // Gauss-Legendre and Gauss-Jacobi quadrature (you must fill these)
    void gauss_legendre(int Q, std::vector<double>& points, std::vector<double>& weights);
    void gauss_jacobi(int Q, double alpha, double beta, std::vector<double>& points, std::vector<double>& weights);

    // Evaluation of M_{ijk,lmn}
    double compute_Mijklmn(int i, int j, int k, int l, int m, int n, int N, int Q);

    // Calcul de la matrice de masse
    Eigen::MatrixXd calculerMatriceMasse();

    std::vector<std::tuple<double, double, double>> getNoeuds();
    std::vector<std::tuple<int, int, int>> getTriplets();
    std::map<std::tuple<int, int, int>, int> getTripletToIndex();


private:
    int ordre;
    int nbNoeuds;
    std::vector<std::tuple<double, double, double>> noeuds;
    std::vector<std::tuple<int, int, int>> triplets;
    std::map<std::tuple<int, int, int>, int> tripletToIndex;

};

