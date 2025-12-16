// simpix.cpp
// Simpix simulated annealing: rearrange pixels of source to resemble target
// Uses Metropolis acceptance and geometric cooling.
// Output: final image + optional collage (src, tgt, out)
//
// Compile (ROOT required):
//   g++ -O3 -std=c++17 simpix.cpp $(root-config --cflags --libs) -lASImage -lGui -o simpix
//
// Run:
//   ./simpix A.png B.png out.png collage.png 5000000
//
// Notes:
// - Works in batch mode (no GUI pop-up).
// - For the "both directions" requirement: run again swapping src/tgt.

#include "TROOT.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>

// TIMER
#include <chrono>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;

struct RGB {
    uint8_t r, g, b;
};

static inline RGB unpack_rgb(UInt_t argb) {
    RGB c;
    c.b = (uint8_t)((argb >> 0)  & 0xFF);
    c.g = (uint8_t)((argb >> 8)  & 0xFF);
    c.r = (uint8_t)((argb >> 16) & 0xFF);
    return c;
}

static inline UInt_t pack_argb(uint8_t a, const RGB& c) {
    return (UInt_t(a) << 24) | (UInt_t(c.r) << 16) | (UInt_t(c.g) << 8) | UInt_t(c.b);
}

// Squared Euclidean RGB distance (no sqrt) — faster, monotonic, fine for annealing
static inline double dist2(const RGB& a, const RGB& b) {
    double dr = double(a.r) - double(b.r);
    double dg = double(a.g) - double(b.g);
    double db = double(a.b) - double(b.b);
    return dr*dr + dg*dg + db*db;
}

// Total energy = sum over pixels of dist2(current[i], target[i])
static double total_energy(const std::vector<RGB>& cur, const std::vector<RGB>& tgt) {
    assert(cur.size() == tgt.size());
    double E = 0.0;
    for (size_t i = 0; i < cur.size(); i++) {
        E += dist2(cur[i], tgt[i]);
    }
    return E;
}

// ΔE if we swap pixels at i and j in "cur"
static inline double deltaE_swap(const std::vector<RGB>& cur,
                                 const std::vector<RGB>& tgt,
                                 size_t i, size_t j) {
    double E_old = dist2(cur[i], tgt[i]) + dist2(cur[j], tgt[j]);
    double E_new = dist2(cur[j], tgt[i]) + dist2(cur[i], tgt[j]);
    return E_new - E_old;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage:\n"
             << "  " << argv[0] << " <source.png> <target.png> [out.png] [collage.png] [trials]\n\n"
             << "Example:\n"
             << "  " << argv[0] << " A.png B.png out.png collage.png 5000000\n";
        return 0;
    }

    std::string fsrc = argv[1];
    std::string ftgt = argv[2];
    std::string fout = (argc > 3) ? argv[3] : "out.png";
    std::string fcollage = (argc > 4) ? argv[4] : "collage.png";
    long long trials = (argc > 5) ? std::stoll(argv[5]) : 5000000LL;

    // TIMER start (whole program)
    using clock_t = std::chrono::steady_clock;
    auto t_start = clock_t::now();
    auto t_last  = t_start;

    // ROOT app (needed by TASImage on many setups)
    TApplication app("simpix_app", &argc, argv);
    gROOT->SetBatch(kTRUE); // no GUI

    // Load images
    TASImage srcImg(fsrc.c_str());
    TASImage tgtImg(ftgt.c_str());

    int W = srcImg.GetWidth();
    int H = srcImg.GetHeight();

    if (W <= 0 || H <= 0) {
        cerr << "Error: could not read source image: " << fsrc << endl;
        return 1;
    }
    if (tgtImg.GetWidth() != W || tgtImg.GetHeight() != H) {
        cerr << "Error: images must have same dimensions.\n"
             << "Source: " << W << "x" << H << "\n"
             << "Target: " << tgtImg.GetWidth() << "x" << tgtImg.GetHeight() << endl;
        return 1;
    }

    const size_t N = size_t(W) * size_t(H);
    cout << "Pixel Geometry: " << W << " x " << H << "  (N=" << N << ")\n";
    cout << "Trials: " << trials << "\n";

    // Read pixels
    UInt_t* srcPix = srcImg.GetArgbArray();
    UInt_t* tgtPix = tgtImg.GetArgbArray();
    if (!srcPix || !tgtPix) {
        cerr << "Error: could not access pixel arrays.\n";
        return 1;
    }

    // Current state = pixels from source
    std::vector<RGB> cur(N), tgt(N), bestState(N);

    for (size_t i = 0; i < N; i++) {
        cur[i] = unpack_rgb(srcPix[i]);
        tgt[i] = unpack_rgb(tgtPix[i]);
    }

    // RNG
    std::mt19937_64 rng(1234567ULL); // fixed seed for reproducibility
    std::uniform_int_distribution<size_t> pick(0, N - 1);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    // Initial energy
    double Ecur = total_energy(cur, tgt);
    double Ebest = Ecur;
    bestState = cur;

    cout << "Initial energy: " << Ecur << "\n";

    // Choose T0 by sampling positive ΔE values (melting temperature idea)
    double dEmax = 0.0;
    const int warmSamples = 2000;
    for (int k = 0; k < warmSamples; k++) {
        size_t i = pick(rng);
        size_t j = pick(rng);
        if (i == j) continue;
        double dE = deltaE_swap(cur, tgt, i, j);
        if (dE > dEmax) dEmax = dE;
    }

    double T0 = (dEmax > 0.0) ? 1.2 * dEmax : 1.0;
    double T = T0;

    // Cooling
    const double alpha = 0.999999; // closer to 1 = slower cooling
    const long long reportEvery = std::max(100000LL, trials / 50);

    // Annealing loop
    for (long long step = 1; step <= trials; step++) {
        size_t i = pick(rng);
        size_t j = pick(rng);
        if (i == j) continue;

        double dE = deltaE_swap(cur, tgt, i, j);

        bool accept = false;
        if (dE <= 0.0) {
            accept = true;
        } else if (T > 0.0) {
            double p = std::exp(-dE / T);
            if (uni01(rng) < p) accept = true;
        }

        if (accept) {
            std::swap(cur[i], cur[j]);
            Ecur += dE;

            if (Ecur < Ebest) {
                Ebest = Ecur;
                bestState = cur;
            }
        }

        // cool
        T *= alpha;

        // TIMER progress report
        if (step % reportEvery == 0) {
            auto t_now = clock_t::now();
            double sec_total = std::chrono::duration<double>(t_now - t_start).count();
            double sec_chunk = std::chrono::duration<double>(t_now - t_last).count();
            double rate_total = (sec_total > 0.0) ? (double(step) / sec_total) : 0.0;
            double rate_chunk = (sec_chunk > 0.0) ? (double(reportEvery) / sec_chunk) : 0.0;
            t_last = t_now;

            cout << std::fixed << std::setprecision(3)
                 << "step " << step
                 << "  T=" << T
                 << "  Ecur=" << Ecur
                 << "  Ebest=" << Ebest
                 << "  time=" << sec_total << " s"
                 << "  rate=" << rate_total << " steps/s"
                 << "  rate_recent=" << rate_chunk << " steps/s"
                 << "\n";
        }
    }

    cout << "Done.\nFinal best energy: " << Ebest << "\n";

    // Build output image from bestState
    TASImage outImg(srcImg); // same geometry
    UInt_t* outPix = outImg.GetArgbArray();
    const uint8_t A = 255;

    for (size_t i = 0; i < N; i++) {
        outPix[i] = pack_argb(A, bestState[i]);
    }

    // Save output image
    outImg.WriteImage(fout.c_str());
    cout << "Wrote: " << fout << "\n";

    // Make collage (src, tgt, out)
    TCanvas c1("c1", "simpix", 1200, 900);
    c1.Divide(2, 2);

    c1.cd(1); srcImg.Draw("X");
    c1.cd(2); tgtImg.Draw("X");
    c1.cd(3); outImg.Draw("X");

    c1.Print(fcollage.c_str());
    cout << "Wrote: " << fcollage << "\n";

    // TIMER end (whole program)
    auto t_end = clock_t::now();
    double sec_total = std::chrono::duration<double>(t_end - t_start).count();
    cout << std::fixed << std::setprecision(3)
         << "Runtime (s): " << sec_total << "\n";

    return 0;
}
