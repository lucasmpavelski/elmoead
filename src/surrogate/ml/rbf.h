#pragma once

#include "elm_common.h"
#include "kmeans.h"
#include "../var/de_operators.h"

namespace SLFN {

class RBF {
    KMeans kmeans;
    SLFNNetwork rbf_net,
                mul_net,
                imul_net;
    SLFNPredictor rbf_pred,
                  mul_pred,
                  imul_pred;
    double rbf_weight,
           mul_weight,
           imul_weight;
    matrix mX;
    matrix mY;
    matrix H;
    Eigen::LLT<matrix> solver;
    double max_diff;
    DifferentialEvolution de;

public:

    RBF(size_t no_ins, size_t no_hidden, size_t no_outs,
        const svectord& lb, const svectord& ub)
        : kmeans(no_ins, no_hidden)
        , rbf_net(no_ins, no_hidden, no_outs)
        , mul_net(no_ins, no_hidden, no_outs)
        , imul_net(no_ins, no_hidden, no_outs)
        , rbf_pred (ELMActFunc::RBF, no_ins)
        , mul_pred (ELMActFunc::MULTIQUADRIC, no_ins)
        , imul_pred(ELMActFunc::INV_MULTIQUADRIC, no_ins)
        , max_diff(0)
        , de(DEOperator(DEType::RAND_1_BIN, 1.0, 0.5),
             PMOperator(1.0 / no_hidden, 20.0))
    {
        for (int i = 0; i < no_ins; ++i)
            if (ub[i] - lb[i] > max_diff)
                max_diff = ub[i] - lb[i];
    }

    void predict(const double* x, double* y) const {
        predict(x, y, 1);
    }

    void predict(const double* x, double* y, size_t no_samples) const {
        using Eigen::Map;
        Map<const rmatrix> mx(x, no_samples, noInputs());
        Map<rmatrix> my(y, no_samples, noOutputs());
        predictAll(mx, my);
    }

    void predictAll(const_rmatrix_ref x, rmatrix_ref y) const {
        for (int i = 0; i < x.rows(); ++i)
            predict(x.row(i), y.row(i));
    }

    void predict(const_rvector_ref x, rvector_ref y) const {
        y = rvector::Constant(noOutputs(), 0);
        y += rbf_weight * rbf_pred.predict(rbf_net, x);
        y += mul_weight * mul_pred.predict(mul_net, x);
        y += imul_weight * imul_pred.predict(imul_net, x);
        //if (in_norm)
        //{
        //    mx = x;
        //    normalizeInput(mx);
        //    pred.predict(net, mx, y);
        //}
        //else
        //{
        //    pred.predict(net, x, y);
        //}

        //if (out_norm)
        //{
        //    for (int i = 0; i < y.cols(); ++i)
        //        y(i) = outputs_normalizers[i].denormalize(y(i));
        //}
    }

    ELMActFunc getActFunc() const {return ELMActFunc::SIGMOID; }
    void setActFunc(ELMActFunc af) { }

    double getC() const { return 0; }
    void setC(double C) {}

    size_t noInputs() const { return rbf_net.noInputs(); }
    size_t noHidden() const { return rbf_net.noHidden(); }
    size_t noOutputs() const { return rbf_net.noOutputs(); }

    // Trainer
    void train(const double X[], const double Y[], size_t no_samples) {
        train(Eigen::Map<const rmatrix>(X, no_samples, noInputs()),
              Eigen::Map<const rmatrix>(Y, no_samples, noOutputs()));
    }

    struct OptimizeDevs : public MOProblem {
        const_rmatrix_ref ins;
        const_rmatrix_ref outs;
        rmatrix pred_outs;
        SLFNNetwork& net;
        SLFNPredictor& pred;

        OptimizeDevs(size_t no_hidden, SLFNNetwork& net, SLFNPredictor& pred,
                     double max, const_rmatrix_ref ins, const_rmatrix_ref outs)
            : MOProblem("optimize_devs", no_hidden, 1, svectord(no_hidden, 0.0), svectord(no_hidden, max))
            , ins(ins)
            , outs(outs)
            , net(net)
            , pred(pred)
            , pred_outs(outs.rows(), outs.cols())
        {}

        void setNetwork(SLFNNetwork& net, SLFNPredictor& pred) {
            this->net = net;
            this->pred = pred;
        }

        void evaluate(const double x[], double y[]) {
            net.bias() = Eigen::Map<const vector>(x, net.noHidden());
            pred.predict(net, ins, pred_outs);
            y[0] = (outs - pred_outs).rowwise().squaredNorm().sum();
        }

        void validate(double x[]) {
            truncateToBounds(x);
        }
    };

    void train(const_rmatrix_ref X, const_rmatrix_ref Y) {
        throw_assert(X.cols() == noInputs(), "data with incorrenct"
            "dimesion (" << X.cols() << " instead of " << noInputs());
        throw_assert(Y.cols() == noOutputs(), "data with incorrenct"
            "number of outputs (" << Y.cols() << " instead of " << noOutputs());

        mX = X.transpose();
        mY = Y;

        auto no_samples = mX.cols();
        auto no_hidden = noHidden();
        //normalizeData(mX, mY);
        //trainer.updateWeights(mX, mY, pred, net);

        kmeans.apply(mX);
        rbf_net.inputWeights() = kmeans.centers;
        mul_net.inputWeights() = kmeans.centers;
        imul_net.inputWeights() = kmeans.centers;

        vector_ref bias = rbf_net.bias();
        bias = vector::Zero(no_hidden);

        for (int i = 0; i < no_samples; ++i) {
            auto idx = kmeans.membership[i];
            bias(idx) += (mX.col(idx) - kmeans.centers.col(idx)).squaredNorm();
        }

        bias = (bias / (no_samples - 1)).cwiseSqrt();

        mul_net.bias() = rbf_net.bias();
        imul_net.bias() = rbf_net.bias();

        H.resize(no_hidden, no_samples);

        for (int i = 0; i < no_samples; ++i)
            rbf_pred.activationFunction(rbf_net, mX.col(i), H.col(i));
        solver.compute(H * H.transpose());
        rbf_net.outputWeights() = solver.solve(H * mY);

        for (int i = 0; i < no_samples; ++i)
            mul_pred.activationFunction(mul_net, mX.col(i), H.col(i));
        solver.compute(H * H.transpose());
        mul_net.outputWeights() = solver.solve(H * mY);

        for (int i = 0; i < no_samples; ++i)
            imul_pred.activationFunction(imul_net, mX.col(i), H.col(i));
        solver.compute(H * H.transpose());
        imul_net.outputWeights() = solver.solve(H * mY);

        //const rmatrix ins = mX.transpose();
        //const rmatrix outs = mY;
        //OptimizeDevs prob(no_hidden, rbf_net, rbf_pred, max_diff, ins, outs);

        //const size_t pop_size = 10;
        //const size_t no_evals = pop_size * 50;
        //Population pop(pop_size, prob, NO_SAMPLING);
        //auto min_obj = [](const Solution& a, const Solution& b) { return a.objs[0] < b.objs[0]; };

        //decltype(pop.begin()) it;

        //pop.vars() = rbf_net.bias().transpose().replicate(pop_size, 1);
        //pop.evaluateAll();
        //de.optimize(no_evals, prob, pop, NullCallback);
        //it = std::min_element(pop.begin(), pop.end(), min_obj);
        //rbf_net.bias() = Eigen::Map<vector>(it->vars, no_hidden);

        //prob.setNetwork(mul_net, mul_pred);
        //pop.vars() = mul_net.bias().transpose().replicate(pop_size, 1);
        //pop.evaluateAll();
        //de.optimize(no_evals, prob, pop, NullCallback);
        //it = std::min_element(pop.begin(), pop.end(), min_obj);
        //mul_net.bias() = Eigen::Map<vector>(it->vars, no_hidden);

        //prob.setNetwork(imul_net, imul_pred);
        //pop.vars() = imul_net.bias().transpose().replicate(pop_size, 1);
        //pop.evaluateAll();
        //de.optimize(no_evals, prob, pop, NullCallback);
        //it = std::min_element(pop.begin(), pop.end(), min_obj);
        //imul_net.bias() = Eigen::Map<vector>(it->vars, no_hidden);

        rbf_weight = 0;
        mul_weight = 0;
        imul_weight = 0;

        rvector out(noOutputs());
        rvector in(noInputs());
        double rbf_mse, mul_mse, imul_mse;
        for (int i = 0; i < no_samples; ++i) {
            in = mX.col(i);

            rbf_pred.predict(rbf_net, in, out);
            rbf_mse  = (out - mY.row(i)).squaredNorm();

            mul_pred.predict(mul_net, in, out);
            mul_mse  = (out - mY.row(i)).squaredNorm();

            imul_pred.predict(imul_net, in, out);
            imul_mse = (out - mY.row(i)).squaredNorm();

            if (rbf_mse <= mul_mse && rbf_mse <= imul_mse)
                rbf_weight += 1.0 / no_samples;
            else if (mul_mse <= rbf_mse && mul_mse <= imul_mse)
                mul_weight += 1.0 / no_samples;
            else
                imul_weight += 1.0 / no_samples;
        }
    }
};

}
