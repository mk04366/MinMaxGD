/**
 * Example 4.1 from "Modeling and control of TEGs with shared resources".
 */

#include "../include/lminmaxgd.h"
#include "../src/gd.cpp"
#include "../src/poly.cpp"
#include "../src/serie.cpp"
#include "../src/smatrix.cpp"
#include "../src/tools.cpp"

using namespace mmgd;

// Build the ultimately periodic series (g^start_g d^start_d)(g^period_g d^period_d)*.
serie periodicSeries(long start_g, long start_d, long period_g, long period_d)
{
    // In serie::init, epsilon is p, the first starred term is q, and the period is r.
    poly repeatingStart;
    repeatingStart.init(start_g, start_d);
    gd period(period_g, period_d);
    serie result;
    result.init(epsilon, repeatingStart, period);
    return result;
}

// Build a 1-by-1 matrix containing a reference counter as a series.
smatrix referenceMatrix(poly &counter)
{
    smatrix result(1, 1);
    result(0, 0).init(epsilon, counter, e);
    return result;
}

// Return true when every entry of the calculated output is below the reference.
bool meetsReference(smatrix &output, smatrix &reference)
{
    return inf(output, reference) == output;
}

// Convert one monomial into a series so it can participate in the series operators.
serie monomialSeries(long g, long d)
{
    gd monomial(g, d);
    return monomial;
}

// Compute the scalar product of two counters using the library's semiring product.
serie seriesProduct(serie left, serie right)
{
    return otimes(left, right);
}

// Compute the infimum of two counters.
serie seriesInfimum(serie left, serie right)
{
    return inf(left, right);
}

// Compute the semiring sum used when multiplying a row vector by a column vector.
serie seriesSum(serie left, serie right)
{
    return oplus(left, right);
}

// Compute the left residual of two counters: numerator \ denominator.
serie seriesLeftResidual(serie numerator, serie denominator)
{
    return frac(numerator, denominator);
}

// Compute the dual Hadamard residual represented by ]^{-1} in Equation (4.11).
serie seriesHadamardDualResidual(serie left, serie right)
{
    return hadamard_dualres(right, left);
}

// Evaluate P2 u2, where P2 = [e d5, e d0] from Example 4.1.
serie allocationS2(serie &u20, serie &u21)
{
    serie p21 = monomialSeries(0, 5);
    serie p22 = monomialSeries(0, 0);
    return seriesSum(seriesProduct(p21, u20), seriesProduct(p22, u21));
}

// Evaluate P3 u3, where P3 = [e d0] from Example 4.1.
serie allocationS3(serie &u3)
{
    serie p31 = monomialSeries(0, 0);
    return seriesProduct(p31, u3);
}

// Compute F2 from Equation (4.11), using S1 as the higher-priority subsystem.
smatrix fixedPointMapS2(serie &u20, serie &u21, serie &allocationS1, serie &releaseS1,
                        serie &tracking0, serie &tracking1)
{
    serie resource = monomialSeries(2, 1);
    serie h2 = monomialSeries(0, 4);
    serie allocation = allocationS2(u20, u21);
    // The resource-sharing product in (4.11) is Hadamard, not semiring multiplication.
    serie constrainedAllocation = hadamard_prod(allocationS1, allocation);
    serie resourceResidual = hadamard_res(resource, constrainedAllocation);
    // Equation (4.11)'s ]^{-1} is the library's dual Hadamard residual.
    serie releaseResidual = seriesHadamardDualResidual(resourceResidual, releaseS1);
    // Residuation by H2 P2 is component-wise: each denominator is a column of H2 P2.
    serie bound0 = seriesLeftResidual(releaseResidual, seriesProduct(h2, monomialSeries(0, 5)));
    serie bound1 = seriesLeftResidual(releaseResidual, seriesProduct(h2, monomialSeries(0, 0)));
    smatrix result(2, 1);
    result(0, 0) = seriesInfimum(seriesInfimum(tracking0, bound0), u20);
    result(1, 0) = seriesInfimum(seriesInfimum(tracking1, bound1), u21);
    return result;
}

// Compute F3 from Equation (4.11), using S1 and S2 as higher-priority subsystems.
serie fixedPointMapS3(serie &u3, serie &allocationS1, serie &allocationS2Value,
                      serie &releaseS1, serie &releaseS2, serie &tracking)
{
    serie resource = monomialSeries(2, 1);
    serie h3 = monomialSeries(0, 3);
    serie allocation = allocationS3(u3);
    serie higherAllocation = seriesInfimum(allocationS1, allocationS2Value);
    // The resource-sharing product in (4.11) is Hadamard, not semiring multiplication.
    serie constrainedAllocation = hadamard_prod(higherAllocation, allocation);
    serie resourceResidual = hadamard_res(resource, constrainedAllocation);
    serie higherRelease = seriesInfimum(releaseS1, releaseS2);
    // Apply the dual Hadamard residual to the higher-priority release bound.
    serie releaseResidual = seriesHadamardDualResidual(resourceResidual, higherRelease);
    serie bound = seriesLeftResidual(releaseResidual, seriesProduct(h3, monomialSeries(0, 0)));
    return seriesInfimum(seriesInfimum(tracking, bound), u3);
}

// Iterate a monotone map from its tracking bound until the greatest fixed point stabilizes.
smatrix solveS2(serie &allocationS1, serie &releaseS1, smatrix &tracking)
{
    serie tracking0 = tracking(0, 0);
    serie tracking1 = tracking(1, 0);
    serie current0 = tracking0;
    serie current1 = tracking1;
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        smatrix next = fixedPointMapS2(current0, current1, allocationS1, releaseS1, tracking0, tracking1);
        if (next(0, 0) == current0 && next(1, 0) == current1)
            break;
        current0 = next(0, 0);
        current1 = next(1, 0);
    }
    smatrix result(2, 1);
    result(0, 0) = current0;
    result(1, 0) = current1;
    return result;
}

// Iterate a monotone map from its tracking bound until the greatest fixed point stabilizes.
serie solveS3(serie &allocationS1, serie &allocationS2Value, serie &releaseS1,
              serie &releaseS2, serie &tracking)
{
    serie current = tracking;
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        serie next = fixedPointMapS3(current, allocationS1, allocationS2Value, releaseS1, releaseS2, tracking);
        if (next == current)
            break;
        current = next;
    }
    return current;
}

int main()
{
    // G1 maps the two inputs of subsystem S1 to its output.
    smatrix G1(1, 2);
    G1(0, 0) = periodicSeries(0, 17, 1, 10);
    G1(0, 1) = periodicSeries(0, 15, 1, 10);

    // G2 maps the two inputs of subsystem S2 to its output.
    smatrix G2(1, 2);
    G2(0, 0) = periodicSeries(0, 9, 2, 5);
    G2(0, 1) = periodicSeries(0, 4, 2, 5);

    // G3 maps the single input of subsystem S3 to its output.
    smatrix G3(1, 1);
    G3(0, 0) = periodicSeries(0, 3, 2, 4);

    // z1 requires four output firings by time 52.
    poly z1Counter;
    z1Counter.init(0, 52)(4, infinit);
    smatrix z1 = referenceMatrix(z1Counter);

    // z2 requires three output firings by time 27 and two more at time 39.
    poly z2Counter;
    z2Counter.init(0, 27)(3, 39)(5, infinit);
    smatrix z2 = referenceMatrix(z2Counter);

    // z3 requires three output firings by time 9 and two more at time 35.
    poly z3Counter;
    z3Counter.init(0, 9)(3, 35)(5, infinit);
    smatrix z3 = referenceMatrix(z3Counter);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 1 ---------------------------------
    // -----------------------------------------------------------------
    // For the highest-priority subsystem, no resource restriction comes from another subsystem.
    // Left residuation computes the greatest u1 satisfying G1 otimes u1 <= z1.
    smatrix u1opt = lfrac(z1, G1);
    // P1 = [e d2, e d9(1 d10)*] maps u1 to S1's allocation schedule.
    serie p11 = monomialSeries(0, 2);
    serie p12 = periodicSeries(0, 9, 1, 10);

    // TODO:
    serie allocationS1 = seriesSum(seriesProduct(p11, u1opt(0, 0)), seriesProduct(p12, u1opt(1, 0)));

    // H1 = e d6 maps S1's allocation schedule to its resource-release schedule.
    serie releaseS1 = seriesProduct(monomialSeries(0, 6), allocationS1);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 2 ---------------------------------
    // -----------------------------------------------------------------
    // The tracking bounds for S2 are the two entries of G2 \ z2.
    smatrix s2Tracking = lfrac(z2, G2);

    // Compute u2opt as the greatest fixed point of mapping
    smatrix u2opt = solveS2(allocationS1, releaseS1, s2Tracking);
    // Compute S2's allocation and release schedules for use by the S3 fixed-point map.
    serie allocationS2Value = allocationS2(u2opt(0, 0), u2opt(1, 0));
    serie releaseS2 = seriesProduct(monomialSeries(0, 4), allocationS2Value);

    // The tracking bound for S3 is the scalar residual G3 \ z3.
    smatrix s3TrackingMatrix = lfrac(z3, G3);
    serie s3Tracking = s3TrackingMatrix(0, 0);
    // Compute u3opt as the greatest fixed point of F3 with both higher-priority schedules fixed.
    serie u3optValue = solveS3(allocationS1, allocationS2Value, releaseS1, releaseS2, s3Tracking);
    smatrix u3opt(1, 1);
    u3opt(0, 0) = u3optValue;

    // Apply each transfer matrix to its selected input to obtain the predicted output.
    smatrix y1opt = otimes(G1, u1opt);
    smatrix y2opt = otimes(G2, u2opt);
    smatrix y3opt = otimes(G3, u3opt);

    // Print the transfer functions and references so the encoded paper example is visible.
    std::cout << "G1:" << std::endl
              << G1 << std::endl;
    std::cout << "G2:" << std::endl
              << G2 << std::endl;
    std::cout << "G3:" << std::endl
              << G3 << std::endl;
    std::cout << "z1:" << std::endl
              << z1 << std::endl;
    std::cout << "z2:" << std::endl
              << z2 << std::endl;
    std::cout << "z3:" << std::endl
              << z3 << std::endl;

    // Print the three optimal inputs and the outputs they produce.
    std::cout << "u1opt (computed by lfrac):" << std::endl
              << u1opt << std::endl;
    std::cout << "y1opt:" << std::endl
              << y1opt << std::endl;
    std::cout << "u2opt (computed greatest fixed point):" << std::endl
              << u2opt << std::endl;
    std::cout << "y2opt:" << std::endl
              << y2opt << std::endl;
    std::cout << "u3opt (computed greatest fixed point):" << std::endl
              << u3opt << std::endl;
    std::cout << "y3opt:" << std::endl
              << y3opt << std::endl;

    // Check the defining tracking inequalities yk <= zk for all three subsystems.
    std::cout << "Tracking checks: "
              << "S1=" << (meetsReference(y1opt, z1) ? "OK" : "FAILED")
              << ", S2=" << (meetsReference(y2opt, z2) ? "OK" : "FAILED")
              << ", S3=" << (meetsReference(y3opt, z3) ? "OK" : "FAILED")
              << std::endl;

    return 0;
}
