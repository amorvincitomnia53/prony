#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

using namespace std::complex_literals;
struct Dipole {  //双極子を表す
    std::complex<double> z;
    std::complex<double> m;
    Dipole conjugate() const
    {
        return {std::conj(1.0 / z), std::conj(m / (z * z))};
    }
};
constexpr double PI = 3.1415926535897932384626;
double log_threshold = -2;  // 閾値Tの常用対数

constexpr double ln10 = std::log(10.0);

std::complex<double> calcPotential(const Dipole& d, const std::complex<double> z)
{
    //一つの双極子dが与えられたとき、それがzの位置に作る複素ポテンシャルを計算する。
    if (d.m == 0.0)
        return 0.0;
    return (1.0 / (2 * PI)) * d.m / (z - d.z);
}
std::complex<double> calcPotential(const std::vector<Dipole>& dipoles, const std::complex<double> z)
{
    //複数の双極子が与えられたとき、それがzの位置に作る複素ポテンシャルを計算する。
    std::complex<double> ret = 0;
    for (auto& d : dipoles) {
        ret += calcPotential(d, z);
    }
    return ret;
}

std::complex<double> calcField(const std::vector<Dipole>& dipoles, const std::complex<double> z)
{
    //複数の双極子が与えられたとき、それがzの位置に作る磁場を計算する。
    std::complex<double> ret = 0;
    for (auto& d : dipoles) {
        if (d.m == 0.0)
            continue;
        std::complex<double> n = d.m / std::abs(d.m);
        auto rel = (z - d.z) / n;
        ret += (1.0i / (2 * PI)) * d.m / std::conj(rel * rel);
    }
    return ret;
}

struct Entry {                // 推定関数の入力となる測定データ。
    std::complex<double> z;   // 測定位置
    std::complex<double> p;   // 複素ポテンシャル
    std::complex<double> dz;  // 測定位置の間隔。周回積分で用いる
};

std::vector<Entry> potentialTable(const std::vector<Dipole>& dipoles, int div)
{
    // 円周上に測定器を配置したときの測定データをシミュレーションする。
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
    // zのpow乗をかけて複素積分することを数値的に実行する。
    std::complex<double> sum = 0;
    for (auto [z, p, dz] : table) {
        sum += dz * p * std::pow(z, pow);
    }
    return -1.0i * sum;
}
std::vector<std::complex<double>> calcIntegralAll(const std::vector<Entry>& table, int pow_max)
{
    // zのk乗をかけて複素積分するということをk=0..powまで一括で実行する。
    std::vector<std::complex<double>> ret(pow_max + 1);
    for (int i = 0; i <= pow_max; i++) {
        ret[i] = calcIntegral(table, i);
    }
    return ret;
}

std::vector<Dipole> estimate(const std::vector<Entry>& data, int pow_max)
{
    // データとパラメタP(=pow)が与えられたときに双極子推定を行う。
    auto c = calcIntegralAll(data, pow_max);  // 複素積分を全て実行する。
    int n = c.size() / 2;                     // ハンケル行列のサイズを求める。

    // ハンケル行列を作る。
    Eigen::MatrixXcd matrix(n, n);
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
            matrix(y, x) = c[x + y];
        }
    }

    // ハンケル行列C_0を特異値分解する。
    auto svd = matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

    svd.setThreshold(std::exp(ln10 * log_threshold));  // 閾値Tを設定
    int r = svd.rank();                                // ハンケル行列のランクを求める。

    // 右辺のハンケル行列C_1を作る。
    Eigen::MatrixXcd matrix2(n, n);
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
            matrix2(y, x) = c[x + y + 1];
        }
    }
    // 疑似逆行列によって行列Fを計算する
    Eigen::MatrixXcd transition = svd.solve(matrix2);

    // 行列Fの固有値を求める。
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> eig = transition.eigenvalues();

    // 行列Fの固有値を絶対値が大きい順にソートし、初めからr個のうちで絶対値が1より大きくないものを採用する。
    std::sort(eig.data(), eig.data() + n, [](const std::complex<double>& a, const std::complex<double>& b) { return std::abs(a) > std::abs(b); });
    for (int i = 0; i < r; i++) {
        if (std::abs(eig(i)) > 1.0) {
            std::swap(eig(i), eig(--r));
        }
    }

    // もし双極子が存在しないと推定されたら、後続の処理がエラーとなるためこの時点で関数から抜ける。
    if (r == 0) {
        return {};
    }

    // m_iの推定
    // 左辺の大きな行列(example)と右辺のデータベクトル(data_vector)を作る
    Eigen::MatrixXd example(data.size(), r * 2);
    Eigen::VectorXd data_vector(data.size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < r; j++) {
            auto try_p = [&](const std::complex<double>& m) {
                Dipole d11 = {eig[j], m};
                Dipole d12 = d11.conjugate();
                return calcPotential(d11, data[i].z).imag() + calcPotential(d12, data[i].z).imag();
            };
            example(i, j * 2) = try_p(1.0);
            example(i, j * 2 + 1) = try_p(1.0i);
        }
        data_vector(i) = data[i].p.imag();
    }

    // 最小二乗法によりm_iを推定する。
    Eigen::VectorXd sol = example.fullPivHouseholderQr().solve(data_vector);

    // 出力形式を揃える。
    std::vector<Dipole> ret(r);
    for (int i = 0; i < r; i++) {
        ret[i].m = {sol(2 * i), sol(2 * i + 1)};
        ret[i].z = eig(i);
    }

    return ret;
}
