#include "reg_elm.h"

/*void RegELM::train(const double* X, const int nsmp, const double* Y)
{
    mX = Map<const Eigen::MatrixXd>(X, dim, nsmp);
    if (in_norm)
        input_normalizer.findRangeAndNormalize(mX.data(),
                                               mX.data() + nsmp * dim);

    mY = Map<const Eigen::MatrixXd>(Y, nsmp, no_out);

    if (out_norm)
    {
        for (int i = 0; i < no_out; ++i)
            outputs_normalizers[i].findRangeAndNormalize(mY.col(i).data(),
                                                         mY.col(i).data() + nsmp);
    }

    inW.setRandom(dim, no_hidden);
    //inW = 0.5 * (MatrixXd::Random(dim, no_hidden) + MatrixXd::Ones(dim, no_hidden));
    bias.setRandom(no_hidden);
    //bias = (0.5 * (VectorXd::Random(no_hidden) + VectorXd::Ones(no_hidden)));

    H.resize(no_hidden, nsmp);
    for (int i = 0; i < nsmp; ++i)
        H.col(i) = activation_function(mX.col(i), inW, bias);

    //MatrixXd preH = elm.inW * mX + elm.bias.replicate(1, nsmp);
    //const MatrixXd H = (1 + (-preH.array()).exp()).cwiseInverse();
    //A = MatrixXd::Identity(no_hidden, no_hidden) * inv_C + H * H.transpose();

    const double inv_C = 1.0 / C;

    if (no_hidden <= nsmp)
    {
        solver.compute(H * H.transpose() + MatrixXd::Identity(no_hidden, no_hidden) * inv_C);
        outW = solver.solve(H * mY);
    }
    else
    {
        solver.compute(H.transpose() * H + MatrixXd::Identity(nsmp, nsmp) * inv_C);
        outW = H * solver.solve(mY);
    }

   //  plt.plot(mX.data(), 300, 30);


}

double RegELM::predict(const double raw_x[], double out[]) const
{
    if (in_norm)
       input_normalizer.normalizeAll(raw_x, raw_x + dim, x.data());
    else
        x = Map<const Eigen::VectorXd>(raw_x, dim);

    hno = activation_function(x, inW, bias);


    for (int i = 0; i < no_out; ++i)
         out[i] = hno.dot(outW.col(i));

    if (out_norm)
    {
        for (int i = 0; i < no_out; ++i)
            out[i] = outputs_normalizers[i].denormalize(out[i]);
    }

    return out[0];
}

double* RegELM::predict(const double raw_x[], const int nsmp,
                          double scores[]) const
{
    for (int i = 0; i < nsmp; ++i)
        predict(raw_x + dim * i, scores + i * no_out);
    return scores;
}

double* RegELM::predict(const double raw_x[], const int nsmp) const
{
    double* scores = new double[nsmp * no_out];
    return predict(raw_x, nsmp, scores);
}
*/
