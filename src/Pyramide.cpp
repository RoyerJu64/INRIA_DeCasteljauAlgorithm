#include "Pyramide.h"

Pyramide::Pyramide(int ordre) : ordre(ordre), nbNoeuds(0)
{
	genererNoeuds();
    genererTriplets();
    genererTripletToIndex();
	nbNoeuds = noeuds.size();
}

void Pyramide::genererNoeuds()
{
    noeuds.clear();

    if (ordre == 1) {
        // 5 noeuds : 4 coins de la base + sommet
        noeuds.push_back({ 0.0, 0.0, 0.0 });  // (-1, -1, 0)
        noeuds.push_back({ 1.0, 0.0, 0.0 });  // (1, -1, 0)
        noeuds.push_back({ 1.0, 1.0, 0.0 });  // (1, 1, 0)
        noeuds.push_back({ 0.0, 1.0, 0.0 });  // (-1, 1, 0)
        noeuds.push_back({ 0.5, 0.5, 1.0 });  // (0, 0, 1)
    }
    else if (ordre == 2) {
        // coins base
        noeuds.push_back({ 0.0, 0.0, 0.0 });  // (-1, -1, 0)
        noeuds.push_back({ 1.0, 0.0, 0.0 });  // (1, -1, 0)
        noeuds.push_back({ 1.0, 1.0, 0.0 });  // (1, 1, 0)
        noeuds.push_back({ 0.0, 1.0, 0.0 });  // (-1, 1, 0)
        noeuds.push_back({ 0.5, 0.5, 1.0 });  // (0, 0, 1)

        // milieux des côtés de la base
        noeuds.push_back({ 0.5, 0.0, 0.0 });  // (0, -1, 0)
        noeuds.push_back({ 1.0, 0.5, 0.0 });  // (1, 0, 0)
        noeuds.push_back({ 0.5, 1.0, 0.0 });  // (0, 1, 0)
        noeuds.push_back({ 0.0, 0.5, 0.0 });  // (-1, 0, 0)

        // milieux des arêtes verticales
        noeuds.push_back({ 0.25, 0.25, 0.5 });  // (-0.5, -0.5, 0.5)
        noeuds.push_back({ 0.75, 0.25, 0.5 });  // (0.5, -0.5, 0.5)
        noeuds.push_back({ 0.75, 0.75, 0.5 });  // (0.5, 0.5, 0.5)
        noeuds.push_back({ 0.25, 0.75, 0.5 });  // (-0.5, 0.5, 0.5)

        // centre de la pyramide
        noeuds.push_back({ 0.5, 0.5, 0.5 });    // (0, 0, 0.5)
    }
}

void Pyramide::genererTriplets()
{
    triplets.clear();
    for (int k = 0; k <= ordre; ++k) {
        int max_ij = ordre - k;
        for (int i = 0; i <= max_ij; ++i) {
            for (int j = 0; j <= max_ij; ++j) {
                triplets.emplace_back(i, j, k);
            }
        }
    }
}

void Pyramide::genererTripletToIndex()
{
    tripletToIndex.clear();
    int c = 0;
    for (int k = 0; k <= ordre; ++k) {
        int max_ij = ordre - k;
        for (int i = 0; i <= max_ij; ++i) {
            for (int j = 0; j <= max_ij; ++j) {
                tripletToIndex.emplace(std::make_tuple(i, j, k), c);
                ++c;
            }
        }
    }
}


double Pyramide::bernstein(int i, int n, double t) {
    assert(t >= 0.0 && t <= 1.0);
    assert(i >= 0 && i <= n);

    // Initialisation : points de contrôle pour B_i^n
    std::vector<double> ctrl(n + 1, 0.0);
    ctrl[i] = 1.0;  // B_i^n est défini par une valeur 1 en position i

    // Algorithme de De Casteljau
    for (int r = 1; r <= n; ++r) {
        for (int j = 0; j <= n - r; ++j) {
            ctrl[j] = (1.0 - t) * ctrl[j] + t * ctrl[j + 1];
        }
    }

    return ctrl[0];
}

std::vector<double> Pyramide::baseBernsteinCube(double x, double y, double z)
{
    std::vector<double> values;
    for (int k = 0; k <= ordre; ++k) {
        int max_ij = ordre - k;
        for (int i = 0; i <= max_ij; ++i) {
            for (int j = 0; j <= max_ij; ++j) {
                values.push_back(bernstein(i, ordre - k, x) * bernstein(j, ordre - k, y) * bernstein(k, ordre, z));
            }
        }
    }
    return values;
}

std::vector<double> Pyramide::baseBernsteinPyramide(double x, double y, double z)
{
    return baseBernsteinCube(x/(1-z), y/(1-z), z);
}

double Pyramide::J(double a, double b)
{
    return 1.0;
}

void Pyramide::gauss_legendre(int Q, std::vector<double>& points, std::vector<double>& weights)
{
    points.resize(Q);
    weights.resize(Q);
    
    if (Q == 4) {
        points = { -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 };
        weights = { 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 };
        for (int i = 0; i < Q; ++i) {
            points[i] = 0.5 * (1.0 + points[i]);     // Transformation vers [0,1]
            weights[i] = 0.5 * weights[i];            // Ajustement du poids
        }
    }
    else {
        std::cerr << "Gauss-Legendre: ordre non supporté (Q = " << Q << ")\n";
        exit(1);
    }
}

void Pyramide::gauss_jacobi(int Q, double alpha, double beta, std::vector<double>& points, std::vector<double>& weights)
{
    points.resize(Q);
    weights.resize(Q);

    if (Q == 4 && std::abs(alpha - 2.0) < 1e-10 && std::abs(beta - 0.0) < 1e-10) {
        points = {
            -0.6167071458,
            -0.1656632050,
             0.2852315165,
             0.6888754087
        };
        weights = {
            0.1654099274,
            0.3895801273,
            0.3482550502,
            0.0967548951
        };
        for (int i = 0; i < Q; ++i) {
            points[i] = 0.5 * (1.0 + points[i]);     // Transformation vers [0,1]
            weights[i] = 0.5 * weights[i];            // Ajustement du poids
        }
    }
    else {
        std::cerr << "Gauss-Jacobi: ordre non supporté (Q = " << Q << ", alpha = " << alpha << ", beta = " << beta << ")\n";
        exit(1);
    }
}

double Pyramide::compute_Mijklmn(int i, int j, int k, int l, int m, int n, int N, int Q)
{
    std::vector<double> a_pts(Q), a_wts(Q);
    std::vector<double> b_pts(Q), b_wts(Q);
    std::vector<double> c_pts(Q), c_wts(Q);

    gauss_legendre(Q, a_pts, a_wts);
    gauss_legendre(Q, b_pts, b_wts);
    gauss_jacobi(Q, 2.0, 0.0, c_pts, c_wts);

    double M = 0.0;

    for (int qa = 0; qa < Q; ++qa) {
        double a = a_pts[qa];
        double wa = a_wts[qa];
        double Bi_a = bernstein(i, N-k, a); // B_i^{N-k}(a)
        double Bl_a = bernstein(l, N-n, a); // B_l^{N-n}(a)

        for (int qb = 0; qb < Q; ++qb) {
            double b = b_pts[qb];
            double wb = b_wts[qb];
            double Bj_b = bernstein(j, N-k, b); // B_j^{N-k}(b)
            double Bm_b = bernstein(m, N-n, b); // B_m^{N-k}(b)
            double J_val = J(a, b); // Bilinear Jacobian

            for (int qc = 0; qc < Q; ++qc) {
                double c = c_pts[qc];
                double wc = c_wts[qc];
                double Bk_c = bernstein(k, N, c); // B_k^N(c)
                double Bn_c = bernstein(n, N, c); // B_n^N(c)

                M += wa * wb * wc *
                     Bi_a * Bl_a *
                     Bj_b * Bm_b *
                     Bk_c * Bn_c *
                     J_val;
            }
        }
    }

    return M;
}

Eigen::MatrixXd Pyramide::calculerMatriceMasse()
{   
    Eigen::MatrixXd M(nbNoeuds, nbNoeuds);
    
    for (int k = 0; k <= ordre + 1; k++) {
        for (int n = 0; n < ordre + 1; n++) {
            for (int i = 0; i < ordre - k + 1; i++) {
                for (int j = 0; j < ordre - k + 1; j++) {
                    for (int m = 0; m < ordre - n + 1; m++) {
                        for (int l = 0; l < ordre - n + 1; l++) {
                            M(tripletToIndex[std::make_tuple(i,j,k)],tripletToIndex[std::make_tuple(l,m,n)]) = compute_Mijklmn(i,j,k,l,m,n,ordre,4);
                        }
                    }
                }
            }
        }
    }
    return M;
}



std::vector<std::tuple<double, double, double>> Pyramide::getNoeuds()
{
    return noeuds;
}

std::vector<std::tuple<int, int, int>> Pyramide::getTriplets()
{
    return triplets;
}

std::map<std::tuple<int, int, int>, int> Pyramide::getTripletToIndex()
{
    return tripletToIndex;
}
