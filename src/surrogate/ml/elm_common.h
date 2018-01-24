#pragma once

#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include "../aux.h"
#include "../plotter.h"
#include "normalizers.h"

namespace SLFN {

// types
using matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
using vector = Eigen::VectorXd;

using rmatrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using rvector = Eigen::RowVectorXd;

using vector_ref = Eigen::Ref<vector>;
using matrix_ref = Eigen::Ref<matrix>;

using rvector_ref = Eigen::Ref<rvector>;
using rmatrix_ref = Eigen::Ref<rmatrix>;

using const_vector_ref = const Eigen::Ref<const vector>;
using const_matrix_ref = const Eigen::Ref<const matrix>;

using const_rvector_ref = Eigen::Ref<const rvector>;
using const_rmatrix_ref = Eigen::Ref<const rmatrix>;

struct SLFNNetwork
{
public:

    matrix in_weights;
    vector bias_vec;
    matrix out_weights;

    SLFNNetwork(size_t no_ins, size_t no_hidden, size_t no_outs) :
        in_weights(no_ins, no_hidden),
        bias_vec(no_hidden),
        out_weights(no_hidden, no_outs)
    {}

    size_t noInputs()  const { return in_weights.rows();  }
    size_t noHidden()  const { return in_weights.cols();  }
    size_t noOutputs() const { return out_weights.cols(); }

    void setNoInputs(size_t ni) { in_weights.resize(ni, noHidden()); }
    void setNoOutputs(size_t no) { out_weights.resize(noHidden(), no); }
    void setNoHidden(size_t nh) {
        in_weights.resize(noInputs(), nh);
        bias_vec.resize(nh);
    }

    const_matrix_ref inputWeightsCRef() const { return in_weights; }
    const_matrix_ref outputWeightsCRef() const { return out_weights; }
    const_vector_ref biasCRef() const { return bias_vec; }

    matrix_ref inputWeights() { return in_weights; }
    matrix_ref outputWeights() { return out_weights; }
    vector_ref bias() { return bias_vec; }

    void setInputWeights(const_matrix_ref iw) { in_weights = iw; }
    void setOutputWeights(const_matrix_ref ow) { out_weights = ow; }
    void setBias(const_vector_ref b) { bias_vec = b; }

    void resetWeights() {
        in_weights.setRandom(noInputs(), noHidden());
        bias_vec.setRandom(noHidden());
    }

    void resetWeights(const double scale) {
        resetWeights();
        in_weights *= scale;
        bias_vec *= scale;
    }

    void resetWeights(const double lb, const double ub) {
        const int dim = noInputs();
        const int no_hidden = noHidden();
        const double hdt = (ub - lb) / 2.0;
        in_weights = hdt * (matrix::Random(dim, no_hidden) +
                            matrix::Ones(dim, no_hidden)) +
                     lb * matrix::Ones(dim, no_hidden);
        bias_vec = hdt * (vector::Random(no_hidden) + vector::Ones(no_hidden)) +
                   lb * vector::Ones(no_hidden);
    }

    void plot(Plotter& plt) const {
        std::ofstream out("/tmp/plt.dat");
        for (int i = 0; i < noInputs(); ++i)
        {
            for (int j = 0; j < noHidden(); ++j)
            {
                out << 0 << " " << double(i+1) / noInputs() << " 0\n";
                out << 1 << " " << double(j+1) / noHidden() << " " << in_weights(i, j) << '\n';
            }
        }
        for (int i = 0; i < noHidden(); ++i)
        {
            for (int j = 0; j < noOutputs(); ++j)
            {
                out << 1 << " " << double(i+1) / noHidden() << " 0\n";
                out << 2 << " " << double(j+1) / noOutputs() << " " << out_weights(i, j) << '\n';
            }
        }
        out.close();
        std::string cmd = "plot '/tmp/plt.dat' u 1:2:3 with lines palette;";
        plt.command(cmd);
    }
};

enum class ELMActFunc {
    SIGMOID,
    GAUSSIAN,
    MULTIQUADRIC,
    HARDLIMIT,
    RBF,
    TANH,
    SIN,
    COS,
    INV_MULTIQUADRIC
};

std::ostream& operator<<(std::ostream& os, const ELMActFunc& af);


inline ELMActFunc str2ELMActFunc(const std::string& str) 
{
    if (str == "SIGMOID"     ) return ELMActFunc::SIGMOID;
    if (str == "GAUSSIAN"    ) return ELMActFunc::GAUSSIAN;
    if (str == "MULTIQUADRIC") return ELMActFunc::MULTIQUADRIC;
    if (str == "HARDLIMIT"   ) return ELMActFunc::HARDLIMIT;
    if (str == "RBF"         ) return ELMActFunc::RBF;
    if (str == "TANH"        ) return ELMActFunc::TANH;
    if (str == "SIN"         ) return ELMActFunc::SIN;
    if (str == "COS"         ) return ELMActFunc::COS;
    if (str == "INV_MULTIQUADRIC") return ELMActFunc::INV_MULTIQUADRIC;

    throw_assert(false,
        "elm_act_func should be: "
        "(SIGMOID | GAUSSIAN | MULTIQUADRIC | HARDLIMIT | RBF | TANH | SIN | COS | INV_MULTIQUADRIC)"
    );
    return ELMActFunc::SIGMOID;
}

struct SLFNPredictor
{
public:
    ELMActFunc af;
    mutable Eigen::RowVectorXd hno;

public:

    SLFNPredictor(const ELMActFunc af = ELMActFunc::SIGMOID, size_t no_hidden=1)
        : af(af)
        , hno(no_hidden)
    {}

    SLFNPredictor(const ELMActFunc af, const SLFNNetwork& elm_net)
        : SLFNPredictor{af, elm_net.noHidden()}
    {}

    ELMActFunc getActFunc() const { return af; }
    void setActFunc(ELMActFunc af) { this->af = af; }
    /*void predict(const SLFNNetwork& net, const double* x, double* y) const {
        predict(net, x, y, 1);
    }*/

    void predict(const SLFNNetwork& net, const_rmatrix_ref x, rmatrix_ref y) const {
        throw_assert(x.rows() == y.rows(), "data inputs/outputs size differ");
        throw_assert(x.cols() == net.noInputs(), "sample with wrong number of inputs");
        throw_assert(y.cols() == net.noOutputs(), "sample with wrong number of outputs");
        hno.resize(net.noHidden());
        for (size_t i = 0; i < x.rows(); ++i)
        {
            //predict(net, x.row(i), y.row(i));
            activationFunction(net, x.row(i), hno);
            y.row(i) = hno * net.out_weights;
        }
    }

    rvector predict(const SLFNNetwork& net, const_rvector_ref x) const {
        rvector y(net.noOutputs());
        predict(net, x, y);
        return y;
    }

    inline void activationFunction(const SLFNNetwork& elm_net, const_vector_ref x,
                                   vector_ref hno) const {
        using tp = typename vector_ref::Scalar;
        using std::ptr_fun;
        const auto& w = elm_net.inputWeightsCRef();
        const auto& b = elm_net.biasCRef();
        switch (getActFunc())
        {
        case ELMActFunc::SIGMOID :
            hno = ((-w.transpose() * x - b).array().exp()
                   + 1.0).cwiseInverse();
            break;
        case ELMActFunc::GAUSSIAN :
            hno = (- (w.colwise() - x).colwise().squaredNorm().array()
                   * b.transpose().array().abs()).exp();
            break;
        case ELMActFunc::MULTIQUADRIC :
            hno = ((w.colwise() - x).colwise().squaredNorm() +
                   b.cwiseAbs2().transpose()).cwiseSqrt();
            break;
        case ELMActFunc::HARDLIMIT :
            hno = (w.transpose() * x - b).unaryExpr([](const double& xa){
                return (xa >= 0) ? 1.0 : 0.0;
            });
            break;
        case ELMActFunc::RBF :
            hno = (- (w.colwise() - x).colwise().squaredNorm().array()
                   / b.transpose().array().abs()).exp();
            break;
        case ELMActFunc::TANH :
            hno = (w.transpose() * x + b).unaryExpr(ptr_fun<tp,tp>(std::tanh));
            break;
        case ELMActFunc::SIN :
            hno = (w.transpose() * x + b).unaryExpr(ptr_fun<tp,tp>(std::sin));
            break;
        case ELMActFunc::COS :
            hno = (w.transpose() * x + b).unaryExpr(ptr_fun<tp,tp>(std::cos));
            break;
        case ELMActFunc::INV_MULTIQUADRIC :
            hno = ((w.colwise() - x).colwise().squaredNorm() +
                   b.cwiseAbs2().transpose()).cwiseSqrt().cwiseInverse();
            break;
        }
    }
};

class SLFNTrainer
{
//    bool in_norm, out_norm;
//    std::vector<Normalizer<double>>    input_normalizers;
//    std::vector<AgvNormalizer<double>> outputs_normalizers;

//    matrix mX;
//    matrix mY;
//    mutable vector mx;

public:

//    SLFNTrainer(size_t no_ins, size_t no_outs,
//                const bool in_norm, const bool out_norm,
//                const svectord& lb, const svectord& ub) :
//        in_norm(in_norm), out_norm(out_norm),
//        input_normalizers(no_ins),
//        outputs_normalizers(no_outs, AgvNormalizer<double>()),
//        mx(no_ins)
//    {
//        for (int i = 0; i < no_ins; ++i)
//        {
//            input_normalizers[i].setOldRange(lb[i], ub[i]);
//            input_normalizers[i].setNewRange(-1.0, 1.0);
//        }
//    }

//    SLFNTrainer() : SLFNTrainer(0, 0, false, false, svectord{}, svectord{}) {}

//    bool normalizedInput() const { return in_norm; }
//    bool normalizedOutput() const { return out_norm; }

//    virtual ~SLFNTrainer() {}

//    void normalizeInput(matrix_ref mX) {
//        for (int i = 0; i < mX.cols(); ++i)
//            for (int j = 0; j < mX.rows(); ++j)
//                mX(j, i) = input_normalizers[j].normalize(mX(j, i));
//    }

//    void normalizeOutput(matrix_ref mY) {
//        for (int i = 0; i < mY.cols(); ++i)
//        {
//            auto beg = mY.col(i).data(), end = mY.col(i).data() + mY.rows();
//            outputs_normalizers[i].findRangeAndNormalize(beg, end);
//        }
//    }



//    void train(const double X[], const double* Y, size_t no_samples,
//               const SLFNPredictor& pred, SLFNNetwork& net) {
//        train(Eigen::Map<const rmatrix>(X, no_samples, net.noInputs()),
//              Eigen::Map<const rmatrix>(Y, no_samples, net.noOutputs()),
//              pred, net);
//    }

//    void train(const_rmatrix_ref X, const_rmatrix_ref Y,
//               const SLFNPredictor& pred, SLFNNetwork& net) {
//        throw_assert(X.cols() == net.noInputs(), "data with incorrenct"
//            "dimesion (" << X.cols() << " instead of " << net.noInputs());
//        throw_assert(Y.cols() == net.noOutputs(), "data with incorrenct"
//            "number of outputs (" << Y.cols() << " instead of " << net.noOutputs());

//        mX = X.transpose();
//        mY = Y;

//        normalizeData(mX, mY);
//        updateWeights(mX, mY, pred, net);
//    }

    virtual void updateWeights(const_matrix_ref mX, const_matrix_ref mY,
                               const SLFNPredictor& pred, SLFNNetwork& net) = 0;

//    ELMActFunc getActFunc() const { return elm_pred.getActFunc(); }
//    void setActFunc(const ELMActFunc af) { elm_pred.setActFunc(af); }


/*    void activationFunction(const_vector_ref x, vector_ref hno) const {
        elm_pred.activationFunction(elm_net, x, hno);
    }
*/
/*    void train(const double *X, const int nsmp, const double *Y) {
        const int dim = elm_net.noInputs();
        const int no_out = elm_net.noOutputs();

//        cout << Y[0] << " " << Y[1] << " " << Y[8] << endl;

        mX = Eigen::Map<const matrix>(X, dim, nsmp);
        mY = Eigen::Map<const matrix>(Y, nsmp, no_out);

//        cout << mY << endl;

//        cout << "mX:\n" << mX << endl;
//        cout << "mY:\n" << mY << endl;

        normalizeData(mX, mY);
        updateWeights(mX, mY);
    }
*/
/*
    double predict(const double raw_x[], double out[]) const {
        const int dim = elm_net.noInputs();
        const int no_out = elm_net.noOutputs();

        if (in_norm)
        {
            for (int i = 0; i < dim; ++i)
                mx(i) = input_normalizers[i].normalize(raw_x[i]);
           elm_pred.predict(elm_net, mx.data(), out);
        }
        else
        {
            elm_pred.predict(elm_net, raw_x, out);
        }


        if (out_norm)
        {
            for (int i = 0; i < no_out; ++i)
                out[i] = outputs_normalizers[i].denormalize(out[i]);
        }

        return out[0];
    }

    const SLFN& network() const { return elm_net; }

    double* predict(const double raw_x[], const int nsmp, double scores[]) const
    {
        for (int i = 0; i < nsmp; ++i)
            predict(raw_x + elm_net.noInputs() * i,
                    scores + i * elm_net.noOutputs());
        return scores;
    }
    */
};

template<class T>
struct SLFNRegressor {
    SLFNNetwork net;
    SLFNPredictor pred;
    T trainer;

    bool in_norm, out_norm;
    std::vector<Normalizer<double>>    input_normalizers;
    std::vector<AgvNormalizer<double>> outputs_normalizers;

    mutable matrix mX;
    mutable matrix mY;
    mutable rvector mx;

    SLFNRegressor(SLFNNetwork&& net,
                  SLFNPredictor&& pred,
                  T&& trainer,
                  size_t no_ins, size_t no_outs,
                  const bool in_norm, const bool out_norm,
                  const svectord& lb, const svectord& ub)
        : net{std::move(net)}
        , pred{std::move(pred)}
        , trainer{std::move(trainer)}
        , in_norm{in_norm}
        , out_norm{out_norm}
        , input_normalizers(no_ins)
        , outputs_normalizers(no_outs, AgvNormalizer<double>())
    {
        if (in_norm)
        {
            for (int i = 0; i < no_ins; ++i)
            {
                input_normalizers[i].setOldRange(lb[i], ub[i]);
                input_normalizers[i].setNewRange(-1.0, 1.0);
            }
        }
    }

    bool isInputNorm() const { return in_norm; }
    bool isOutputNorm() const { return out_norm; }

    void normalizeInputs(matrix_ref mX) const {
        for (int i = 0; i < mX.cols(); ++i)
        {
            for (int j = 0; j < mX.rows(); ++j)
                mX(j, i) = input_normalizers[j].normalize(mX(j, i));
        }
    }

    void normalizeInputs(rmatrix_ref mX) const {
        for (int i = 0; i < mX.rows(); ++i)
        {
            for (int j = 0; j < mX.cols(); ++j)
                mX(i, j) = input_normalizers[j].normalize(mX(i, j));
        }
    }

    void normalizeInput(rvector_ref mX) const {
        for (int i = 0; i < mX.cols(); ++i)
            mX(i) = input_normalizers[i].normalize(mX(i));
    }

    void normalizeOutputs(matrix_ref mY) {
        for (int i = 0; i < mY.cols(); ++i)
        {
            auto beg = mY.col(i).data(), end = mY.col(i).data() + mY.rows();
            outputs_normalizers[i].findRangeAndNormalize(beg, end);
        }
    }

    void normalizeData(matrix_ref mX, matrix_ref mY) {
        if (in_norm) normalizeInputs(mX);
        if (out_norm) normalizeOutputs(mY);
    }

    // Network
    size_t getNoInputs() const { return net.noInputs(); }
    size_t getNoOutputs() const { return net.noOutputs(); }
    size_t getNoHidden() const { return net.noHidden(); }

    void setNoInputs(size_t n) { net.setNoInputs(n); }
    void setNoOutputs(size_t n) { net.setNoOutputs(n); }
    void setNoHidden(size_t n) { net.setNoHidden(n); }

    // Predictor
    void predict(const double* x, double* y) const {
        predict(x, y, 1);
    }

    void predict(const double* x, double* y, size_t no_samples) const {
        using Eigen::Map;
        Map<const rmatrix> mx(x, no_samples, net.noInputs());
        Map<rmatrix> my(y, no_samples, net.noOutputs());
        predictAll(mx, my);
    }

//    void predict(const double x[], double y[], size_t no_samples) const {
//        pred.predict(net, x, y, no_samples);
//    }

    void predictAll(const_rmatrix_ref x, rmatrix_ref y) const {
        for (int i = 0; i < x.rows(); ++i)
            predict(x.row(i), y.row(i));
    }

    void predict(const_rvector_ref x, rvector_ref y) const {
        if (in_norm)
        {
            mx = x;
            normalizeInput(mx);
            pred.predict(net, mx, y);
        }
        else
        {
            pred.predict(net, x, y);
        }

        if (out_norm)
        {
            for (int i = 0; i < y.cols(); ++i)
                y(i) = outputs_normalizers[i].denormalize(y(i));
        }
    }

    ELMActFunc getActFunc() const { return pred.getActFunc(); }
    void setActFunc(ELMActFunc af) { pred.setActFunc(af); }

    // Trainer
    void train(const double X[], const double Y[], size_t no_samples) {
        train(Eigen::Map<const rmatrix>(X, no_samples, net.noInputs()),
              Eigen::Map<const rmatrix>(Y, no_samples, net.noOutputs()));
    }

    void train(const_rmatrix_ref X, const_rmatrix_ref Y) {
        throw_assert(X.cols() == net.noInputs(), "data with incorrenct"
            "dimesion (" << X.cols() << " instead of " << net.noInputs());
        throw_assert(Y.cols() == net.noOutputs(), "data with incorrenct"
            "number of outputs (" << Y.cols() << " instead of " << net.noOutputs());
        mX = X.transpose();
        mY = Y;
        normalizeData(mX, mY);
        trainer.updateWeights(mX, mY, pred, net);
    }
};

}
