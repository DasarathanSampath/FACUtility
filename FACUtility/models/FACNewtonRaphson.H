#pragma once

#include "scalar.H"
#include "error.H"
#include "mathematicalConstants.H"

namespace Foam
{

struct NewtonControls
{
    label maxIter {50};
    scalar absTol {1e-24};     // absolute residual tolerance
    scalar relTol {1e-10};     // relative step tolerance
    scalar minStep {1e-18};    // avoid stalling
    scalar maxDampIter {12};   // backtracking attempts
    scalar xMin {1e-20};       // physical lower bound for mH
    scalar xMax {20.0};        // physical upper bound (adjust to your chemistry range)
};

template<class FFunc, class DFunc>
inline bool newtonSolve
(
    scalar& x,                 // in: initial guess, out: solution
    const FFunc& F,            // residual functor: F(x)
    const DFunc& dF,           // derivative functor: dF(x)
    const NewtonControls& ctrl,
    const bool verbose = false
)
{
    x = min(max(x, ctrl.xMin), ctrl.xMax);

    scalar Fx = F(x);
    scalar Fx0 = Fx;

    for (label iter=0; iter<ctrl.maxIter; ++iter)
    {
        const scalar dFx = dF(x);

        if (mag(dFx) < SMALL)
        {
            if (verbose) Info<< "Newton: tiny derivative, abort at iter " << iter << nl;
            return false;
        }

        // Newton step
        scalar dx = -Fx/dFx;

        // Step-based convergence check
        // if (mag(dx) <= max(ctrl.relTol*max(mag(x), scalar(1.0)), ctrl.minStep))
        // {
        //     return true;
        // }
        if (mag(dx) <= max(ctrl.relTol*mag(x), ctrl.minStep))
        {
            x = min(max(x + dx, ctrl.xMin), ctrl.xMax);
            return true;
        }

        // Damped Newton / backtracking to ensure residual reduces
        scalar damp = 1.0;
        scalar xTry = x;
        scalar FTry = Fx;

        bool accepted = false;
        for (label k=0; k<ctrl.maxDampIter; ++k)
        {
            xTry = x + damp*dx;
            xTry = min(max(xTry, ctrl.xMin), ctrl.xMax);

            FTry = F(xTry);

            if (mag(FTry) <= mag(Fx))  // accept if residual decreases
            {
                accepted = true;
                break;
            }
            damp *= 0.5;
        }

        if (!accepted)
        {
            if (verbose) Info<< "Newton: no acceptable damped step, abort\n";
            return false;
        }

        x  = xTry;
        Fx = FTry;

        // Residual convergence check
        if (mag(Fx) <= ctrl.absTol || mag(Fx) <= ctrl.absTol + ctrl.relTol*mag(Fx0))
        {
            return true;
        }
    }

    return false;
}

} // namespace Foam
