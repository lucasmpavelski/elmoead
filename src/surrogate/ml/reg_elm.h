#pragma once

#include <cmath>
#include <Eigen/Cholesky>

#include "../aux.h"
#include "normalizers.h"
#include "elm_common.h"

class RegELMTrainer : public SLFN::SLFNTrainer
{
    using matrix = SLFN::matrix;
    using const_matrix_ref = SLFN::const_matrix_ref;

private:
    double C;
    matrix H;
    Eigen::LLT<matrix> solver;
    //Eigen::LDLT<MatrixXd> solver;
    //Eigen::FullPivLU<MatrixXd> solver;

public:

/*    RegELM(const int no_ins, const int no_hidden, const int no_outs,
           const double C, const ELMActFunc act_fun,
           const svectord& in_up_norms, const svectord& in_low_norms,
           const bool in_norm=true, const bool out_norm=true)
        : SLFNTrainer(in_up_norms, in_low_norms, in_norm, out_norm)
        , net{no_ins, no_hidden, no_outs}
        , pred{act_fun, no_hidden}
        , C(C)
    {}
*/
    RegELMTrainer(double C)
        : C(C)
    {}

    double getC() const { return C; }
    void setC(double new_C) { C = new_C; }

    //double* predict(const double raw_x[], const int nsmp) const;
    //double* predict(const double raw_x[], const int nsmp, double scr[]) const;

    virtual void updateWeights(const_matrix_ref mX, const_matrix_ref mY,
                               const SLFN::SLFNPredictor& pred,
                               SLFN::SLFNNetwork& net) override;

};

class RegELM : public SLFN::SLFNRegressor<RegELMTrainer> {
public:
    RegELM(size_t no_ins, size_t no_hid, size_t no_out,
           SLFN::ELMActFunc af, double C, bool in_norm, bool out_norm,
           const svectord& lb, const svectord& ub)
        : SLFNRegressor{SLFN::SLFNNetwork{no_ins, no_hid, no_out},
                        SLFN::SLFNPredictor{af, no_hid},
                        RegELMTrainer{C},
                        no_ins, no_out, in_norm, out_norm, lb, ub}
    {}

    double getC() const { return trainer.getC(); }
    void setC(double c) { trainer.setC(c); }

    friend std::ostream& operator<<(std::ostream& os, const RegELM& elm) {
        return os << "type: regularized_elm\n"
                  << "act_func: " << elm.getActFunc() << '\n'
                  << "norm_input: " << elm.isInputNorm() << '\n'
                  << "norm_output: " << elm.isOutputNorm() << '\n'
                  << "C: " << elm.getC() << '\n';
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
