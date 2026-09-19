/**
 * Example 4.1 from "Modeling and control of TEGs with shared resources".
 */

#include "../include/lminmaxgd.h"
#include "../src/gd.cpp"
#include "../src/poly.cpp"
#include "../src/serie.cpp"
#include "../src/smatrix.cpp"
#include "../src/tools.cpp"
#include <fstream>
#include <stdexcept>
#include <string>

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

serie seriesSum(serie left, serie right)
{
    return oplus(left, right);
}

// Compute the repository's Hadamard residual used in Equation (4.11).
serie seriesHadamardResidual(serie left, serie right)
{
    return hadamard_res(right, left);
}

serie seriesHadamardProduct(serie left, serie right)
{
    return hadamard_prod(left, right);
} 

// Write finite and periodic counter terms to a CSV file for the plotting script.
void writeCounterCsv(std::ofstream &file, const std::string &subsystem,
                     const std::string &kind, serie &counter, long horizon)
{
    for (unsigned int i = 0; i < counter.getp().getn(); ++i)
    {
        gd term = counter.getp().getpol(i);
        if (term.getg() != infinit && term.getg() != _infinit &&
            (term.getd() <= horizon || term.getd() == infinit))
        {
            file << subsystem << "," << kind << "," << term.getg() << ",";
            if (term.getd() == infinit)
                file << "inf";
            else
                file << term.getd();
            file << "\n";
        }
    }

    for (long repetition = 0; repetition < horizon + 2; ++repetition)
    {
        for (unsigned int i = 0; i < counter.getq().getn(); ++i)
        {
            gd term = counter.getq().getpol(i);
            long gamma = term.getg();
            long delta = term.getd();
            if (gamma != infinit && gamma != _infinit)
                gamma += repetition * counter.getr().getg();
            if (delta != infinit && delta != _infinit)
                delta += repetition * counter.getr().getd();
            if (delta >= 0 && (delta <= horizon || delta == infinit))
            {
                file << subsystem << "," << kind << "," << gamma << ",";
                if (delta == infinit)
                    file << "inf";
                else
                    file << delta;
                file << "\n";
            }
        }
        if (counter.getr().getd() == infinit ||
            repetition * counter.getr().getd() > horizon)
            break;
    }
}

// Find the largest finite time explicitly represented by a counter.
long largestFiniteTime(serie &counter)
{
    long largest = 0;
    for (unsigned int i = 0; i < counter.getp().getn(); ++i)
        if (counter.getp().getpol(i).getd() != infinit)
            largest = std::max(largest, counter.getp().getpol(i).getd());
    for (unsigned int i = 0; i < counter.getq().getn(); ++i)
        if (counter.getq().getpol(i).getd() != infinit)
            largest = std::max(largest, counter.getq().getpol(i).getd());
    return largest;
}

// Export the counters calculated by this executable instead of duplicating them in Python.
void writePlotData(const std::string &path, smatrix &resource,
                   smatrix &z1, smatrix &z2, smatrix &z3,
                   smatrix &y1, smatrix &y2, smatrix &y3,
                   smatrix &allocationS1, smatrix &allocationS2,
                   smatrix &allocationS3, smatrix &releaseS1,
                   smatrix &releaseS2, smatrix &releaseS3)
{
    std::ofstream file(path.c_str());
    if (!file)
        throw std::runtime_error("Could not open plot data file: " + path);
    file << "subsystem,kind,gamma,delta\n";
    file << "RESOURCE,tokens," << resource(0, 0).getq().getpol(0).getg() << ",inf\n";
    const long horizon = std::max(
        std::max(std::max(largestFiniteTime(z1(0, 0)), largestFiniteTime(y1(0, 0))),
                 std::max(largestFiniteTime(z2(0, 0)), largestFiniteTime(y2(0, 0)))),
        std::max(largestFiniteTime(z3(0, 0)), largestFiniteTime(y3(0, 0))));
    writeCounterCsv(file, "S1", "reference", z1(0, 0), horizon);
    writeCounterCsv(file, "S1", "output", y1(0, 0), horizon);
    writeCounterCsv(file, "S2", "reference", z2(0, 0), horizon);
    writeCounterCsv(file, "S2", "output", y2(0, 0), horizon);
    writeCounterCsv(file, "S3", "reference", z3(0, 0), horizon);
    writeCounterCsv(file, "S3", "output", y3(0, 0), horizon);
    writeCounterCsv(file, "S1", "allocation", allocationS1(0, 0), horizon);
    writeCounterCsv(file, "S2", "allocation", allocationS2(0, 0), horizon);
    writeCounterCsv(file, "S3", "allocation", allocationS3(0, 0), horizon);
    writeCounterCsv(file, "S1", "release", releaseS1(0, 0), horizon);
    writeCounterCsv(file, "S2", "release", releaseS2(0, 0), horizon);
    writeCounterCsv(file, "S3", "release", releaseS3(0, 0), horizon);
}

// Reduce the rows of a non-empty schedule matrix with the Hadamard product.
serie bigHadamardProduct(smatrix &values)
{
    serie result;
    result = values(0, 0);
    for (int i = 1; i < values.getrow(); ++i)
        result = seriesHadamardProduct(result, values(i, 0));
    return result;
}

// Evaluate a row matrix P against a column matrix u using semiring matrix multiplication.
smatrix allocationFrom(smatrix &p, smatrix &u)
{
    return otimes(p, u);
}

// Compute the generic mapping Fk from Equation (4.11).
smatrix fixedPointMap(smatrix &current, smatrix &p, smatrix &h,
                      smatrix &higherAllocation, smatrix &higherRelease,
                      smatrix &tracking, smatrix &resource)
{
    serie allocation = allocationFrom(p, current)(0, 0); // x^k_a = P^k u^k // always a counter-series
    std::cout << allocation << " allocation " << std::endl;
    serie higherAllocationProduct = bigHadamardProduct(higherAllocation);
    std::cout << higherAllocationProduct << " higherAllocationProduct " << std::endl;

    smatrix constrainedAllocation = seriesHadamardProduct(higherAllocationProduct, allocation);
    std::cout << constrainedAllocation(0, 0) << " constrainedAllocation " << std::endl;

    smatrix resourceResidual = lfrac(constrainedAllocation, resource);
    std::cout << resourceResidual(0, 0) << " resourceResidual " << std::endl;

    serie releaseProduct = bigHadamardProduct(higherRelease);
    std::cout << releaseProduct << " releaseProduct " << std::endl;
    
    smatrix releaseResidual = seriesHadamardResidual(resourceResidual(0, 0), releaseProduct);

    std::cout << releaseResidual(0, 0) << " releaseResidual " << std::endl;
    smatrix hp = otimes(h, p);
    smatrix inputBound = lfrac(releaseResidual, hp);
    smatrix result = inf(tracking, inputBound);
    return inf(result, current);
}

// Iterate the generic monotone map from its tracking bound until the greatest fixed point stabilizes.
smatrix solveFixedPoint(smatrix &p, smatrix &h,
                        smatrix &higherAllocation, smatrix &higherRelease,
                        smatrix &tracking, smatrix &resource)
{

    // to find the greatest fixed point, we start with the top-element but here, we will choose
    // to start with G\z because that is something already a part of the mapping and the least-constrainted one theoretically
    smatrix current = tracking;
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        smatrix next = fixedPointMap(current, p, h, higherAllocation,
                                     higherRelease, tracking, resource);
        bool stable = next == current;
        current = next;
        if (stable) {
            break;
        }
    }
    return current;
}

int main(int argc, char **argv)
{
    smatrix resource = monomialSeries(2, 1);
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
    // P12 = e d0: this is the input-to-allocation block before the allocation transition.
    serie p12 = monomialSeries(0, 0);
    serie h1 = monomialSeries(0, 6);

    serie allocationS1 = seriesSum(seriesProduct(p11, u1opt(0, 0)), seriesProduct(p12, u1opt(1, 0)));

    serie releaseS1 = seriesProduct(h1, allocationS1);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 2 ---------------------------------
    // -----------------------------------------------------------------
    // The tracking bounds for S2 are the two entries of G2 \ z2.
    smatrix s2Tracking = lfrac(z2, G2);

    // Pass S2's P and H together with all higher-priority schedules to the generic solver.
    smatrix p2(1, 2);
    p2(0, 0) = monomialSeries(0, 5);
    p2(0, 1) = monomialSeries(0, 0);
    smatrix higherAllocationsForS2(1, 1);
    higherAllocationsForS2(0, 0) = allocationS1;
    smatrix higherReleasesForS2(1, 1);
    higherReleasesForS2(0, 0) = releaseS1;
    smatrix h2(monomialSeries(0, 4));

    smatrix u2opt = solveFixedPoint(p2, h2,
                                    higherAllocationsForS2, higherReleasesForS2,
                                    s2Tracking, resource);
    // Compute S2's allocation and release schedules for use by the S3 fixed-point map.
    smatrix allocationS2Value = allocationFrom(p2, u2opt);
    smatrix h2Release(monomialSeries(0, 4));
    smatrix releaseS2 = otimes(h2Release, allocationS2Value);

    // -----------------------------------------------------------------
    // ------------------- SUBSYSTEM 3 ---------------------------------
    // -----------------------------------------------------------------
    // The tracking bound for S3 is the scalar residual G3 \ z3.
    smatrix s3TrackingMatrix = lfrac(z3, G3);
    // Pass both higher-priority schedules to the same generic solver for S3.
    smatrix p3(1, 1);
    p3(0, 0) = monomialSeries(0, 0);
    smatrix higherAllocationsForS3(2, 1);
    higherAllocationsForS3(0, 0) = allocationS1;
    higherAllocationsForS3(1, 0) = allocationS2Value(0, 0);
    smatrix higherReleasesForS3(2, 1);
    higherReleasesForS3(0, 0) = releaseS1;
    higherReleasesForS3(1, 0) = releaseS2(0, 0);
    smatrix h3(monomialSeries(0, 3));
    smatrix u3opt = solveFixedPoint(p3, h3,
                                    higherAllocationsForS3, higherReleasesForS3,
                                    s3TrackingMatrix, resource);
    smatrix allocationS3Value = allocationFrom(p3, u3opt);
    smatrix releaseS3 = otimes(h3, allocationS3Value);

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

    const std::string plotDataPath =
        argc == 3 && std::string(argv[1]) == "--plot-data"
            ? argv[2]
            : "misc/example4_1_plot_data.csv";
    smatrix allocationS1Matrix(allocationS1);
    smatrix releaseS1Matrix(releaseS1);
    writePlotData(plotDataPath, resource, z1, z2, z3, y1opt, y2opt, y3opt,
                  allocationS1Matrix, allocationS2Value, allocationS3Value,
                  releaseS1Matrix, releaseS2, releaseS3);

    return 0;
}
