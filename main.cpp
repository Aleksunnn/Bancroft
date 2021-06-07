#include <dlib/matrix.h>
#include <boost/optional.hpp>

using namespace dlib;

constexpr double Vc = 0.299792458 / 1.00026 / 1000.0;

boost::optional<dlib::matrix<double>>BancroftAlgorithm(matrix<double> sensorTOA)
{
    auto n = sensorTOA.nr();
    if (n < 4 || sensorTOA.nc() != 4)
        return boost::none;
    auto LorenzProduct = [](matrix<double> &u, matrix<double> &v)
    {
        auto n = u.nc();
        matrix<double> result = sum_cols(pointwise_multiply(colm(u, range(0, n - 2)), colm(v, range(0, n - 2))))
                - pointwise_multiply(colm(u, n - 1), colm(v, n - 1));
        return result;
    };

    matrix<double> &A = sensorTOA;
    set_colm(A, 3) = colm(A, 3)*(-Vc);
    matrix<double> Ainv = pinv(A);
    matrix<double> b = LorenzProduct(A, A);
    matrix<double> ones = ones_matrix<double>(n, 1);
    matrix<double> d = trans(0.5*(Ainv*ones));
    matrix<double> e = trans(0.5*(Ainv*b));
    double alpha = LorenzProduct(d, d);
    double beta = 2 * LorenzProduct(d, e) - 1;
    double gamma = LorenzProduct(e, e);

    double D = beta*beta - 4*alpha*gamma;
    if (D < 0)
        return boost::none;

    double lambdaPlus = (-1*beta + sqrt(D))/(2*alpha);
    double lambdaMinus = (-1*beta - sqrt(D))/(2*alpha);

    matrix<double> matrixResult = zeros_matrix<double>(4, 2);
    set_colm(matrixResult, 0) = trans(lambdaPlus*d + e);
    set_colm(matrixResult, 1) = trans(lambdaMinus*d + e);
    set_rowm(matrixResult, 3) = rowm(matrixResult, 3)/Vc;

    return matrixResult;
}

int main(int argc, char *argv[])
{
    matrix<double> sensorTOA;
    auto Bancroft = BancroftAlgorithm(sensorTOA);
    return 0;
}

