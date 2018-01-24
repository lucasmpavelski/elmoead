#pragma once

#include <cmath>
#include <functional>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

#include "../aux.h"
#include "normalizers.h"
#include "elm_common.h"

class OnSeqELM : public SLFNTrainer {

private:
    MatrixXd H;
    Eigen::LLT<MatrixXd> solver;
    //Eigen::LDLT<MatrixXd> solver;
    //Eigen::FullPivLU<MatrixXd> solver;
    MatrixXd P;

public:
    OnSeqELM(const int no_ins, const int no_hidden, const int no_outs,
             const ELMActFunc act_fun, double* x, double* y, int nsmp,
              const svectord& in_up_norms, const svectord& in_low_norms,
              const bool in_norm=true, const bool out_norm=false)
        : SLFNTrainer(no_ins, no_hidden, no_outs, act_fun, in_up_norms,
                     in_low_norms, in_norm, out_norm)
    {
        elm_net.resetWeights();

        MatrixXd nx = Eigen::Map<MatrixXd>(x, no_ins, nsmp);
        MatrixXd ny = Eigen::Map<MatrixXd>(y, nsmp, no_outs);
        normalizeData(nx, ny);

        H.resize(no_hidden, nsmp);
        for (int i = 0; i < nsmp; ++i)
            activationFunction(nx.col(i), H.col(i));

        P = (H * H.transpose()).inverse();
        elm_net.out_weights = P * H * ny;
    };

    virtual ~OnSeqELM() {}

    //double* predict(const double raw_x[], const int nsmp) const;
    //double* predict(const double raw_x[], const int nsmp, double scr[]) const;

    void updateWeights(const_matrix_ref mX, const_matrix_ref mY)
    {
        const int nsmp = mX.cols();
        const int no_hidden = elm_net.noHidden();

        H.resize(no_hidden, nsmp);
        for (int i = 0; i < nsmp; ++i)
            activationFunction(mX.col(i), H.col(i));

        P = P - P * H * (MatrixXd::Identity(nsmp, nsmp)
                         + H.transpose() * P * H).inverse() * H.transpose() * P;
        elm_net.out_weights = elm_net.out_weights + P * H * (mY - H.transpose() * elm_net.out_weights);

        //        elm_net.resetWeights();

        //        bool bias = 0;
        //        H.resize(no_hidden, nsmp + bias);
        //        for (int i = 0; i < nsmp; ++i)
        //             activationFunction(mX.col(i), H.col(i));

        //        if (bias)
        //            H.col(nsmp) = VectorXd::Ones(no_hidden);

        //        //MatrixXd preH = elm.inW * mX + elm.bias.replicate(1, nsmp);
        //        //const MatrixXd H = (1 + (-preH.array()).exp()).cwiseInverse();
        //        //A = MatrixXd::Identity(no_hidden, no_hidden) * inv_C + H * H.transpose();

        //        const double inv_C = 1.0 / C;

        //        if (no_hidden <= nsmp)
        //        {
        //            solver.compute(H * H.transpose() + MatrixXd::Identity(no_hidden, no_hidden) * inv_C);
        //            elm_net.out_weights = solver.solve(H * mY);
        //        }
        //        else
        //        {
        //            solver.compute(H.transpose() * H + MatrixXd::Identity(nsmp + bias,
        //                                                                  nsmp + bias) * inv_C);
        //            elm_net.out_weights = H * solver.solve(mY);
        //        }
    }
};

/*
 * #pragma once

#include <cmath>
#include <functional>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

#include "../aux.h"
#include "normalizers.h"
#include "elm_common.h"

class RegELM
{
    ELMNetwork elm_net;
    ELMPredictor elm_pred;

public:
    using ActivationFunc = VectorXd(*)(const VectorXd&,const MatrixXd&,const VectorXd&);

    RegELM(const int dim, const int no_hidden, const int no_out,
             const double C, ActivationFunc act_fun=RegELM::gaussianAF,
             const bool in_norm=true, const bool out_norm=false,
             std::pair<double,double> in_norm_limits=std::make_pair(-1.0, 1.0),
             std::pair<double,double> out_norm_limits=std::make_pair(-1.0, 1.0)) :
        dim(dim), no_hidden(no_hidden), C(C), no_out(no_out),
        activation_function(act_fun),
        inW(dim, no_hidden), bias(no_hidden), outW(no_hidden, no_out),
        in_norm(in_norm), out_norm(out_norm),
        input_normalizer(in_norm_limits.first, in_norm_limits.second),
        outputs_normalizers(no_out,
                            Normalizer<double>(out_norm_limits.first,
                                               out_norm_limits.second)),
        //outputs_normalizer0(-1, 1),
        //outputs_normalizer1(-1, 1),
        x(dim),
        hno(no_hidden)
    {};

    virtual ~RegELM() {};

    double getC() const { return C; };
    void setC(const double& new_C) { C = new_C; };

    void train(const double *X, const int nsmp, const double *Y);

    double predict(const double raw_x[], double out[]) const;
    double* predict(const double raw_x[], const int nsmp) const;
    double* predict(const double raw_x[], const int nsmp, double scr[]) const;

    static inline VectorXd sigmoidAF(const VectorXd& x, const MatrixXd& w,
                            const VectorXd& bias) {
        return ((-w.transpose() * x - bias).array().exp() + 1.0).cwiseInverse();
    };

    static inline VectorXd hardlimitAF(const VectorXd& x, const MatrixXd& w,
                              const VectorXd& bias) {
        return (w.transpose() * x - bias).unaryExpr([](const double& xa){
            return (xa >= 0) ? 1.0 : 0.0;
        });
    };

    static inline VectorXd gaussianAF(const VectorXd& x, const MatrixXd& w,
                             const VectorXd& bias) {
        return (- (w.colwise() - x).colwise().squaredNorm().array()
                * bias.transpose().array().abs()).exp();
    };

    static inline VectorXd multiquadricAF(const VectorXd& x, const MatrixXd& w,
                                 const VectorXd& bias) {
        return ((w.colwise() - x).colwise().squaredNorm() +
                bias.cwiseAbs2().transpose()).cwiseSqrt();
    };

    static inline VectorXd tanhAF(const VectorXd& x, const MatrixXd& w,
                                 const VectorXd& bias) {
        return (w.transpose() * x + bias).unaryExpr(std::ptr_fun(tanh));
    };


    static inline VectorXd sinAF(const VectorXd& x, const MatrixXd& w,
                                 const VectorXd& bias) {
        return (w.transpose() * x + bias).unaryExpr(std::ptr_fun(sin));
    };

    static inline VectorXd rbfAF(const VectorXd& x, const MatrixXd& w,
                             const VectorXd& bias) {
        return (- (w.colwise() - x).colwise().squaredNorm().array()
                / bias.transpose().array().abs()).exp();
    };

private:
    const int dim;
    const int no_hidden;
    const int no_out;
    double C;
    mutable VectorXd x;
    mutable VectorXd hno;
    ActivationFunc activation_function;

    MatrixXd inW;
    VectorXd bias;
    MatrixXd outW;

    MatrixXd mX;
    MatrixXd mY;

    MatrixXd H;
    Eigen::LLT<MatrixXd> solver;
    //Eigen::LDLT<MatrixXd> solver;
    //Eigen::FullPivLU<MatrixXd> solver;

    bool in_norm, out_norm;
    Normalizer<double> input_normalizer;
    std::vector<Normalizer<double>> outputs_normalizers;
    //Normalizer<double> outputs_normalizer0;
    //NonLinearNormalizer<double> outputs_normalizer1;
};
*/
