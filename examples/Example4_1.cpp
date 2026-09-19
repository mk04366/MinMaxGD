/**
 * Example 4.1 from "Modeling and control of TEGs with shared resources".
 */

#include "../include/lminmaxgd.h"
#include "../src/gd.cpp"
#include "../src/poly.cpp"
#include "../src/serie.cpp"
#include "../src/smatrix.cpp"
#include "../src/tools.cpp"
#include <vector>

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

// Compute the Hadamard residual represented by ]^{-1} in Equation (4.11).
serie seriesHadamardResidual(serie left, serie right)
{
    return hadamard_res(right, left);
}

serie seriesHadamardProduct(serie left, serie right)
{
    return hadamard_prod(left, right);
} 

// Reduce a non-empty list of higher-priority schedules with the componentwise hadamard product.
serie bigHadamardSeries(const std::vector<serie> &values)
{
    serie result = values[0];
    for (size_t i = 1; i < values.size(); ++i)
        result = seriesHadamardProduct(result, values[i]);
    return result;
}

// Evaluate a row vector P against the current input vector u using semiring matrix multiplication.
serie allocationFrom(const std::vector<serie> &p, const std::vector<serie> &u)
{
    serie result = seriesProduct(p[0], u[0]);
    for (size_t i = 1; i < p.size(); ++i)
        result = seriesSum(result, seriesProduct(p[i], u[i]));
    return result;
}

// Compute the generic mapping Fk from Equation (4.11).
std::vector<serie> fixedPointMap(const std::vector<serie> &current,
                                 const std::vector<serie> &p,
                                 const serie &h,
                                 const std::vector<serie> &higherAllocation,
                                 const std::vector<serie> &higherRelease,
                                 const std::vector<serie> &tracking,
                                 const serie &resource)
{
    std::vector<serie> result;
    serie allocation = allocationFrom(p, current); //x^k_a = P^k u^k
    serie higherAllocationMeet = bigHadamardSeries(higherAllocation); //hadamard prod of higher x^k_a 
    serie constrainedAllocation = hadamard_prod(higherAllocationMeet, allocation); // hadamard b/w higher allocs and current alloc
    serie resourceValue = resource;
    serie resourceResidual = hadamard_res(resourceValue, constrainedAllocation); // resource left-div by constrainted allocs
    serie releaseResidual = seriesHadamardResidual(resourceResidual, bigHadamardSeries(higherRelease)); // gets constraint from alloc-release (star)

    // compare each column of
    for (size_t i = 0; i < current.size(); ++i)
    {
        // Residuation by Hk Pk is component-wise because each input is a column of Hk Pk.
        serie inputBound = seriesLeftResidual(releaseResidual, seriesProduct(h, p[i])); //h^kp^k left-res star
        result.push_back(seriesInfimum(seriesInfimum(tracking[i], inputBound), current[i]));
    }
    return result;
}

// Iterate the generic monotone map from its tracking bound until the greatest fixed point stabilizes.
smatrix solveFixedPoint(const std::vector<serie> &p,
                        const serie &h,
                        const std::vector<serie> &higherAllocation,
                        const std::vector<serie> &higherRelease,
                        smatrix &tracking,
                        const serie &resource)
{
    std::vector<serie> current;
    for (size_t i = 0; i < p.size(); ++i)
        current.push_back(tracking((int)i, 0));

    std::vector<serie> trackingValues = current;
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        std::vector<serie> next = fixedPointMap(current, p, h, higherAllocation,
                                                higherRelease, trackingValues, resource);
        bool stable = true;
        for (size_t i = 0; i < current.size(); ++i)
        {
            if (!(next[i] == current[i]))
                stable = false;
        }
        current = next;
        if (stable)
            break;
    }

    smatrix result((int)current.size(), 1);
    for (size_t i = 0; i < current.size(); ++i)
        result((int)i, 0) = current[i];
    return result;
}

int main()
{

    serie resource = monomialSeries(2, 1);
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
    serie u12_ht = monomialSeries(0, 0);
    serie h1 = monomialSeries(0, 6);

    serie allocationS1 = seriesSum(seriesProduct(p11, u1opt(0, 0)), seriesProduct(u12_ht, u1opt(1, 0)));

    serie releaseS1 = seriesProduct(h1, allocationS1);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 2 ---------------------------------
    // -----------------------------------------------------------------
    // The tracking bounds for S2 are the two entries of G2 \ z2.
    smatrix s2Tracking = lfrac(z2, G2);

    // Pass S2's P and H together with all higher-priority schedules to the generic solver.
    std::vector<serie> p2;
    p2.push_back(monomialSeries(0, 5));
    p2.push_back(monomialSeries(0, 0));
    std::vector<serie> higherAllocationsForS2;
    higherAllocationsForS2.push_back(allocationS1);

    std::vector<serie> higherReleasesForS2;
    higherReleasesForS2.push_back(releaseS1);

    serie h2 = monomialSeries(0, 4);

    smatrix u2opt = solveFixedPoint(p2, h2,
                                    higherAllocationsForS2, higherReleasesForS2,
                                    s2Tracking, resource);
    // Compute S2's allocation and release schedules for use by the S3 fixed-point map.
    serie allocationS2Value = allocationFrom(p2, {u2opt(0, 0), u2opt(1, 0)});
    serie releaseS2 = seriesProduct(monomialSeries(0, 4), allocationS2Value);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 3 ---------------------------------
    // -----------------------------------------------------------------
    // The tracking bound for S3 is the scalar residual G3 \ z3.
    smatrix s3TrackingMatrix = lfrac(z3, G3);
    // Pass both higher-priority schedules to the same generic solver for S3.
    std::vector<serie> p3;
    p3.push_back(monomialSeries(0, 0));
    std::vector<serie> higherAllocationsForS3;
    higherAllocationsForS3.push_back(allocationS1);
    higherAllocationsForS3.push_back(allocationS2Value);
    std::vector<serie> higherReleasesForS3;
    higherReleasesForS3.push_back(releaseS1);
    higherReleasesForS3.push_back(releaseS2);
    smatrix u3opt = solveFixedPoint(p3, monomialSeries(0, 3),
                                    higherAllocationsForS3, higherReleasesForS3,
                                    s3TrackingMatrix, resource);

    // -----------------------------------------------------------------
    // ------------------- Outputs -------------------------------------
    // -----------------------------------------------------------------
    // Apply each transfer matrix to its selected input to obtain the predicted output.
    smatrix y1opt = otimes(G1, u1opt);
    smatrix y2opt = otimes(G2, u2opt);
    smatrix y3opt = otimes(G3, u3opt);

    // Print the three optimal inputs and the outputs they produce.
    std::cout << "u1opt (computed by lfrac):" << std::endl
              << u1opt << std::endl;
    std::cout << "u2opt (computed greatest fixed point):" << std::endl
              << u2opt << std::endl;
    std::cout << "u3opt (computed greatest fixed point):" << std::endl
              << u3opt << std::endl;

    std::cout << "y1opt:" << std::endl
              << y1opt << std::endl;
    std::cout << "z1:" << std::endl
              << z1 << std::endl;

    std::cout << "y2opt:" << std::endl
              << y2opt << std::endl;
    std::cout << "z2:" << std::endl
              << z2 << std::endl;
              
    std::cout << "y3opt:" << std::endl
              << y3opt << std::endl;
    std::cout << "z3:" << std::endl
              << z3 << std::endl;

    // Check the defining tracking inequalities yk <= zk for all three subsystems.
    std::cout << "Tracking checks: "
              << "S1=" << (meetsReference(y1opt, z1) ? "OK" : "FAILED")
              << ", S2=" << (meetsReference(y2opt, z2) ? "OK" : "FAILED")
              << ", S3=" << (meetsReference(y3opt, z3) ? "OK" : "FAILED")
              << std::endl;

    return 0;
}
