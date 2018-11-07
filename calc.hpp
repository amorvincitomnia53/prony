#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

using namespace std::complex_literals;
struct Dipole {
    std::complex<double> z;
    std::complex<double> m;
    Dipole conjugate() const
    {
        return {std::conj(1.0 / z), std::conj(m / (z * z))};
    }
};
constexpr double PI = 3.1415926535897932384626;
double threshold = 1e-7;


std::complex<double> calcPotential(const std::vector<Dipole>& dipoles, const std::complex<double> z)
{
    std::complex<double> ret = 0;
    for (auto& d : dipoles) {
        ret += (1.0 / (2 * PI)) * d.m / (z - d.z);
    }
    return ret;
}

std::complex<double> calcField(const std::vector<Dipole>& dipoles, const std::complex<double> z)
{
    std::complex<double> ret = 0;
    for (auto& d : dipoles) {
        if (d.m == 0.0)
            continue;
        std::complex<double> n = d.m / std::abs(d.m);
        auto rel = (z - d.z) / n;
        ret += (1.0 / (2 * PI)) * d.m / std::conj(rel * rel);
    }
    return ret;
}

struct Entry {
    std::complex<double> z;
    std::complex<double> p;
    std::complex<double> dz;
};

std::vector<Entry> potentialTable(const std::vector<Dipole>& dipoles, int div)
{
    std::vector<Entry> ret(div);
    double tick = PI * 2 / div;
    for (int i = 0; i < div; i++) {
        auto z = std::exp(std::complex<double>{0, tick * i});
        auto p = calcPotential(dipoles, z);
        auto dz = 2.0i * PI * z / double(div);
        ret[i] = {z, p, dz};
    }
    return ret;
}


std::complex<double> calcIntegral(const std::vector<Entry>& table, int pow)
{
    std::complex<double> sum = 0;
    for (auto [z, p, dz] : table) {
        sum += dz * p * std::pow(z, pow);
    }
    return -1.0i * sum;
}
std::vector<std::complex<double>> calcIntegralAll(const std::vector<Entry>& table, int pow_max)
{
    std::vector<std::complex<double>> ret(pow_max + 1);
    for (int i = 0; i <= pow_max; i++) {
        ret[i] = calcIntegral(table, i);
    }
    return ret;
}

std::vector<Dipole> estimate(const std::vector<Entry>& data, int pow_max)
{
    auto c = calcIntegralAll(data, pow_max);
    int n = c.size() / 2;
    Eigen::MatrixXcd matrix(n, n);
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
            matrix(y, x) = c[x + y];
        }
    }
    auto svd = matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

    svd.setThreshold(threshold);
    int r = svd.rank();

    Eigen::MatrixXcd matrix2(n, n);
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
            matrix2(y, x) = c[x + y + 1];
        }
    }
    Eigen::MatrixXcd transition = svd.solve(matrix2);
    //    std::cout << matrix << std::endl;
    //    std::cout << matrix2 << std::endl;
    //    std::cout << transition << std::endl;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> eig = transition.eigenvalues();


    std::sort(eig.data(), eig.data() + n, [](const std::complex<double>& a, const std::complex<double>& b) { return std::abs(a) > std::abs(b); });
    //    for (int i = 0; i < r; i++) {
    //        if (std::abs(eig(i)) > 1.0) {
    //            std::swap(eig(i), eig(--r));
    //        }
    //    }

    if (r == 0) {
        return {};
    }

    //        Eigen::MatrixXd example(data.size(), r);
    //        for (int i = 0; i < data.size(); i++) {
    //            for (int j = 0; j < r; j++) {
    //                Dipole d1 = {eig[j], 1.0};
    //                example(i, j) = potentialTable({d1, d1.conjugate()})
    //            }
    //        }
    Eigen::MatrixXcd vandermonde(c.size(), r);
    for (int i = 0; i < r; i++) {
        std::complex<double> z = eig(i);
        std::complex<double> s = 1;
        for (int j = 0; j < c.size(); j++) {
            vandermonde(j, i) = s;
            s *= z;
        }
    }
    //    std::cout << vandermonde << std::endl;


    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vec(c.size());
    for (int i = 0; i < n; i++) {
        vec[i] = c[i];
    }

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> m = vandermonde.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vec);

    //    std::cout << vec << std::endl;
    std::vector<Dipole> ret(r);
    for (int i = 0; i < r; i++) {
        ret[i].m = m(i);
        ret[i].z = eig(i);
    }

    return ret;
}
