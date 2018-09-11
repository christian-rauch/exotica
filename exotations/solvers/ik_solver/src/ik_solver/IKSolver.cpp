/*
 *  Created on: 15 Jul 2014
 *      Author: Yiming Yang
 * 
 * Copyright (c) 2016, University Of Edinburgh 
 * All rights reserved. 
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met: 
 * 
 *  * Redistributions of source code must retain the above copyright notice, 
 *    this list of conditions and the following disclaimer. 
 *  * Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in the 
 *    documentation and/or other materials provided with the distribution. 
 *  * Neither the name of  nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific 
 *    prior written permission. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE. 
 *
 */

#include "ik_solver/IKSolver.h"

#include <lapack/cblas.h>
#include "f2c.h"
#undef small
#undef large
#include <lapack/clapack.h>

REGISTER_MOTIONSOLVER_TYPE("IKsolver", exotica::IKsolver)

namespace exotica
{
IKsolver::IKsolver() : iterations_(-1)
{
}

IKsolver::~IKsolver()
{
}

template <typename Scalar>
struct CwiseClampOp
{
    CwiseClampOp(const Scalar& inf, const Scalar& sup)
        : m_inf(inf), m_sup(sup)
    {
    }
    const Scalar operator()(const Scalar& x) const
    {
        return x < m_inf ? m_inf : (x > m_sup ? m_sup : x);
    }
    Scalar m_inf, m_sup;
};

Eigen::MatrixXd inverseSymPosDef(const Eigen::Ref<const Eigen::MatrixXd>& A_)
{
//    std::cout << "A dim: " << A_.rows() << ", " << A_.cols() << std::endl;
    Eigen::MatrixXd Ainv_ = A_;
    double* AA = Ainv_.data();
    integer info;
    integer nn = A_.rows();
//    std::cout << "nn: " << nn << std::endl;
    // Compute Cholesky
    dpotrf_((char*)"L", &nn, AA, &nn, &info);
    if (info != 0)
    {
        throw_pretty("Can't invert matrix. Cholesky decomposition failed: " << info << std::endl
                                                                            << A_);
    }
    // Invert
    dpotri_((char*)"L", &nn, AA, &nn, &info);
    if (info != 0)
    {
        throw_pretty("Can't invert matrix: " << info << std::endl
                                             << A_);
    }
    Ainv_.triangularView<Eigen::Upper>() = Ainv_.transpose();
//    std::cout << "Ain dim: " << Ainv_.rows() << ", " << Ainv_.cols() << std::endl;
    return Ainv_;
}

void IKsolver::Instantiate(IKsolverInitializer& init)
{
    parameters_ = init;
}

void IKsolver::specifyProblem(PlanningProblem_ptr pointer)
{
    if (pointer->type() != "exotica::UnconstrainedEndPoseProblem")
    {
        throw_named("This IKsolver can't solve problem of type '" << pointer->type() << "'!");
    }
    MotionSolver::specifyProblem(pointer);
    prob_ = std::static_pointer_cast<UnconstrainedEndPoseProblem>(pointer);

    C = Eigen::MatrixXd::Identity(prob_->JN, prob_->JN) * parameters_.C;
    if (parameters_.C == 0.0)
        Cinv = Eigen::MatrixXd::Zero(prob_->JN, prob_->JN);
    else
        Cinv = C.inverse();

    W = prob_->W;
    Winv = W.inverse();

    if (parameters_.Alpha.size() != 1 && prob_->N != parameters_.Alpha.size())
        throw_named("Alpha must have length of 1 or N.");
}

UnconstrainedEndPoseProblem_ptr& IKsolver::getProblem()
{
    return prob_;
}

int IKsolver::getLastIteration()
{
    return iterations_;
}

void IKsolver::Solve(Eigen::MatrixXd& solution)
{
    prob_->resetCostEvolution(getNumberOfMaxIterations() + 1);

    Timer timer;

    if (!prob_) throw_named("Solver has not been initialized!");
    const Eigen::VectorXd q0 = (solution.cols() == prob_->N) ? solution.row(0) : prob_->applyStartState();

    if (prob_->N != q0.rows()) throw_named("Wrong size q0 size=" << q0.rows() << ", required size=" << prob_->N);

    const bool UseNullspace = prob_->qNominal.rows() == prob_->N;

    solution.resize(1, prob_->N);

    Eigen::VectorXd q = q0;
    error = INFINITY;
    int i;
    for (i = 0; i < getNumberOfMaxIterations(); i++)
    {
        prob_->Update(q);
//        std::cout << "ydiff: " << "(" << prob_->Cost.ydiff.rows() << ")" << std::endl << (prob_->Cost.ydiff).transpose() << std::endl;
        Eigen::VectorXd yd = prob_->Cost.S * prob_->Cost.ydiff;
//        std::cout << "yd: " << std::endl << yd.transpose() << std::endl;

        // dirty hack for changing the state within a taskmap
        const Eigen::VectorXd q_task = prob_->getState(0);
        if(q_task.size()>0) {
            q = q_task;
        }

        error = prob_->getScalarCost();

        prob_->setCostEvolution(i, error);

        if (error < parameters_.Tolerance)
        {
            if (debug_) HIGHLIGHT_NAMED("IKsolver", "Reached tolerance (" << error << " < " << parameters_.Tolerance << ")");
            break;
        }

        C = Eigen::MatrixXd::Identity(prob_->Cost.J.rows(), prob_->Cost.J.rows()) * parameters_.C;
        Cinv = C.inverse();
//        std::cout << "J dim: " << prob_->Cost.J.rows() << ", " << prob_->Cost.J.cols() << std::endl;
        Eigen::MatrixXd Jinv = PseudoInverse(prob_->Cost.S * prob_->Cost.J);
//        Eigen::MatrixXd Jinv = PseudoInverse(prob_->Cost.J);
//        std::cout << "Jinv: " << std::endl << Jinv << std::endl;
        Eigen::VectorXd qd = Jinv * yd;
//        std::cout << "qd: " << std::endl << qd.transpose() << std::endl;
        if (UseNullspace) qd += (Eigen::MatrixXd::Identity(prob_->N, prob_->N) - Jinv * prob_->Cost.S * prob_->J) * (q - prob_->qNominal);

        ScaleToStepSize(qd);

        if (parameters_.Alpha.size() == 1)
        {
            q -= qd * parameters_.Alpha[0];
        }
        else
        {
            q -= qd.cwiseProduct(parameters_.Alpha);
        }

        if (qd.norm() < parameters_.Convergence)
        {
            if (debug_) HIGHLIGHT_NAMED("IKsolver", "Reached convergence (" << qd.norm() << " < " << parameters_.Convergence << ")");
            break;
        }
    }
    iterations_ = i + 1;

    solution.row(0) = q;

    planning_time_ = timer.getDuration();
}

Eigen::MatrixXd IKsolver::PseudoInverse(Eigen::MatrixXdRefConst J)
{
    Eigen::MatrixXd Jpinv;

    if (J.cols() < J.rows())
    {
        //(Jt*C^-1*J+W)^-1*Jt*C^-1
        if (parameters_.C != 0)
        {
            Jpinv = inverseSymPosDef(J.transpose() * Cinv * J + W) * J.transpose() * Cinv;
//            std::cout << "Jpinv dim: " << Jpinv.rows() << ", " << Jpinv.cols() << std::endl;
        }
        else
        {
            Jpinv = inverseSymPosDef(J.transpose() * J + W) * J.transpose();
        }
    }
    else
    {
        //W^-1*Jt(J*W^-1*Jt+C)
        Jpinv = Winv * J.transpose() * inverseSymPosDef(J * Winv * J.transpose() + C);
    }

    return Jpinv;
}

void IKsolver::ScaleToStepSize(Eigen::VectorXdRef xd)
{
    double max_vel = xd.cwiseAbs().maxCoeff();
    if (max_vel > parameters_.MaxStep)
    {
        xd = xd * parameters_.MaxStep / max_vel;
    }
}
}
