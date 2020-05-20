#include <iostream>
#include <vector>
#include <tuple>
#include "Eigen/Dense"

using namespace Eigen;

std::tuple<double, double, double, double, double> solutfornumber(std::vector<int> number, std::vector<double> Ti, std::vector<double> zi) {
    MatrixXd A(Ti.size(), 4);

    for (int i = 0; i < number.size(); i++) {
        A(i, 0) = 1;
        A(i, 1) = -1;
        A(i, 2) = Ti[number[i]];
        A(i, 3) = -1;
    }

    VectorXd b(zi.size());

    for (int i = 0; i < number.size(); i++) {
        b(i) = log10(zi[number[i]]) - log10(1024);
    }

    Vector4d x = (A.transpose() * A).completeOrthogonalDecomposition().pseudoInverse() * A.transpose() * b;
    double r0, rc, k, t0;
    r0 = pow(10, x(0));
    rc = pow(10, x(1)) - r0;
    k = x(2);
    t0 = x(3) / k;


    double averagefault = 0;
    for (int i = 0; i < Ti.size(); i++) {
        double r = r0 * pow(10, k * (Ti[i] - t0));
        averagefault += abs(1024 * r / (rc + r) - zi[i]);
    }
    averagefault /= Ti.size();
    return std::make_tuple(averagefault, r0, rc, k, t0);
}

void print(double a, double b, double c, double d) {
    std::cout << "r0 = " << a << std::endl;
    std::cout << "rc = " << b << std::endl;
    std::cout << "k = " << c << std::endl;
    std::cout << "t0 = " << d << std::endl;
}

int main() {
    std::vector<double> Ti = {71, 64, 52, 41, 33, 23, 17, 12, 2, 0, 87, -5};
    std::vector<double> zi = {27, 31, 43, 58, 69, 86, 102, 111, 122, 137, 18, 87};
    double a, b, c, d;

    double minfault = 999999;
    std::vector<int> numberforminfault;
    for (int x = 0; x < 10000; x++) {
        std::vector<int> number;
        for (int i = 0; i < Ti.size(); i++) {
            if (rand() % 10 > 1) {
                number.push_back(i);
            }
        }
        if (number.size() < 5)
            continue;

        auto[averagefault, r0, rc, k, t0] = solutfornumber(number, Ti, zi);
        a = r0;
        b = rc;
        c = k;
        d = t0;

        if (averagefault < minfault) {
            minfault = averagefault;
            numberforminfault = number;
        }
    }

    print(a, b, c, d);
    std::cout << "fault: " << minfault << std::endl;
    return 0;
}
