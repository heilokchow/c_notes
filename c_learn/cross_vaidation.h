#ifndef CROSS_VALIDATION_H
#define CROSS_VALIDATION_H

#define EIGEN_DONT_PARALLELIZE 
#include <iostream>
#include <Eigen/Dense>
#include <thread>
#include <mutex>
#include <random>
#include <memory>
#include "timer.h"

using namespace Eigen;
#ifndef CROSS_NUM
#define CROSS_NUM 100
#endif

std::mutex iomutex;
int _cross_num = CROSS_NUM;

struct Linear_regression {
public:
    MatrixXd X, XX;
    VectorXd b;
    std::array<MatrixXd, 10> XX_cross;
    int n, p;

    Linear_regression() :
        n(8000), p(1000) {
        int n1 = n / 10, n2 = 9 * n1, n3 = 2 * n;
        X.resize(n3, p);
        XX.resize(p, p);
        XX = MatrixXd::Zero(p, p);
        b.resize(n3);
        std::random_device device;
        std::mt19937 generator(device());
        std::normal_distribution<double> normal(0.0, 1.0);
        {
            Timer t0;
            for (size_t i = 0; i < p; i++) {
                for (size_t j = 0; j < n; j++) {
                    X(j, i) = normal(generator);
                    X(j + n, i) = X(j, i);
                }
            }
        }

        for (size_t i = 0; i < 10; i++) {
            Map<MatrixXd, 0, Stride<Dynamic, 1>> X1(&X(i * n1, 0), n1, p, Stride<Dynamic, 1>(n3, 1));
            XX_cross[i].resize(p, p);
            XX_cross[i] = X1.transpose() * X1;
            XX += XX_cross[i];
        }

        for (size_t i = 0; i < 10; i++)
            XX_cross[i] = XX - XX_cross[i];
        b = X.rowwise().sum();
    }
};

void cross_validation_bench(const std::shared_ptr<Linear_regression> x) {

    int n = x->n, p = x->p, curr = 0;
    int n1 = n / 10, n2 = 9 * n1, n3 = 2 * n;
    std::random_device device;
    std::mt19937 generator(device());
    std::normal_distribution<double> normal(0.0, 1.0);
    VectorXd e(n3), y(n3), xy(p), beta(p);
    do
    {
        double err = 0;
        iomutex.lock();
        curr = _cross_num--;
        iomutex.unlock();

        for (size_t i = 0; i < n; i++) {
            e(i) = normal(generator);
            e(i + n) = e(i);
        }
        y = x->b + e;

        for (size_t i = 0; i < 10; i++) {
            Map<VectorXd, 0, InnerStride<1>> Y1(&y(i * n1), n1);
            Map<VectorXd, 0, InnerStride<1>> Y2(&y((i + 1)* n1), n2);
            Map<VectorXd, 0, InnerStride<1>> E1(&e(i * n1), n1);
            Map<VectorXd, 0, InnerStride<1>> E2(&e((i + 1) * n1), n2);
            Map<MatrixXd, 0, Stride<Dynamic, 1>> X1(&x->X(i * n1, 0), n1, p, Stride<Dynamic, 1>(n3, 1));
            Map<MatrixXd, 0, Stride<Dynamic, 1>> X2(&x->X((i + 1) * n1, 0), n2, p, Stride<Dynamic, 1>(n3, 1));
            xy.noalias() = X2.transpose() * Y2;
            beta.noalias() = x->XX_cross[i].llt().solve(xy);
            err += (Y1 - X1 * beta).sum();
        }
        std::cout << curr << ", Average error: " << err / n << std::endl;
    } while (_cross_num > 0);
}

#endif