// ConTimeFromJulia_timeseries.cc
// Your original continuous-time worldline / FK-cluster code
// with MINIMAL additions to dump a per-config time series for autocorrelations.
//
// Build (Windows / MinGW GCC 6.3):
//   g++ -O3 -std=c++14 .\ConTimeFromJulia_timeseries.cc -o ConTimeFromJulia.exe
//
// Build (Linux/HPC):
//   g++ -O3 -std=c++14 -march=native -ffast-math -o ConTimeFromJulia ConTimeFromJulia_timeseries.cc
//
// Run example:
//   ./ConTimeFromJulia --layers 64,96,128 --ratio_min 0.6 --ratio_max 1.4 --ratio_n 17 \
//     --burn_in 500 --sweeps_between 50 --n_configs 200 --seed 1234 --BC 1 \
//     --ts_file timeseries.dat > summary.dat
//
// Notes:
// - Progress output goes to std::cerr.
// - Summary CSV goes to std::cout (redirect to file).
// - Time-series (one line per cfg) goes to --ts_file (CSV).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>   // NEW
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <chrono>

// ------------------------------------------------------------
// Utilities
// ------------------------------------------------------------
static inline std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

static inline std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> r;
    if (n <= 0) return r;
    if (n == 1) { r.push_back(a); return r; }
    r.reserve((size_t)n);
    for (int i = 0; i < n; ++i) {
        double t = double(i) / double(n - 1);
        r.push_back(a + (b - a) * t);
    }
    return r;
}

// Standard error of mean (NOT autocorr-correct)
static inline double sem(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double mean = std::accumulate(x.begin(), x.end(), 0.0) / double(x.size());
    double var = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double d = x[i] - mean;
        var += d * d;
    }
    var /= double(x.size() - 1);
    return std::sqrt(var / double(x.size()));
}

// Progress printing (STDERR so CSV on STDOUT isn't corrupted)
static inline void progress_line(const std::string& msg) {
    std::cerr << msg << std::endl;
}
static inline void progress_tick(int i, int n, int every, const std::string& prefix) {
    if (every <= 0) return;
    if (i == 0 || i == n - 1 || (i % every) == 0) {
        double pct = (n > 0) ? (100.0 * double(i + 1) / double(n)) : 100.0;
        std::cerr << prefix << " " << (i + 1) << "/" << n
            << " (" << std::fixed << std::setprecision(1) << pct << "%)"
            << std::endl;
    }
}

// ------------------------------------------------------------
// Simulation state (0-based indexing)
// ------------------------------------------------------------
struct Sim {
    int L = 0;                 // number of spatial sites
    double DeltaT = 1.0;       // time extent
    double Gamma = 1.0;        // transverse rate
    double K = 1.0;            // spatial coupling
    int BC = 0;                // 0=open time, 1=periodic time

    // base rails
    std::vector<std::vector<double> > time;   // segment start times, per site
    std::vector<std::vector<int> >    spin;   // spin per segment, per site

    // augmented rails (phantom cuts)
    std::vector<std::vector<double> > aug_time;
    std::vector<std::vector<int> >    aug_spin;

    // FK cluster labels for augmented intervals
    std::vector<std::vector<int> > label;

    // adjacency list (0-based)
    std::vector<std::vector<int> > neighbors;

    Sim(int L_,
        const std::vector<std::vector<int> >& neighbors_,
        double DeltaT_ = 1.0, double Gamma_ = 1.0, double K_ = 1.0, int BC_ = 0)
        : L(L_), DeltaT(DeltaT_), Gamma(Gamma_), K(K_), BC(BC_),
        time((size_t)L_), spin((size_t)L_),
        aug_time((size_t)L_), aug_spin((size_t)L_),
        label((size_t)L_), neighbors(neighbors_)
    {
        if ((int)neighbors.size() != L) {
            throw std::runtime_error("neighbors size != L");
        }
    }
};

// ------------------------------------------------------------
// makeRails: sample Poisson(Gamma*DeltaT) flips, uniform times, sort
// ------------------------------------------------------------
static void makeRails(std::mt19937_64& rng, Sim& p) {
    const double lambdaT = p.Gamma * p.DeltaT;

    std::poisson_distribution<int> pois(lambdaT);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    std::bernoulli_distribution coin(0.5);

    for (int x = 0; x < p.L; ++x) {
        std::vector<double>& tvec = p.time[(size_t)x];
        std::vector<int>& svec = p.spin[(size_t)x];
        tvec.clear();
        svec.clear();

        int N = pois(rng);

        // times in (0, DeltaT), sorted; always prepend 0.0
        tvec.push_back(0.0);
        if (N > 0) {
            std::vector<double> times;
            times.reserve((size_t)N);
            for (int i = 0; i < N; ++i) times.push_back(unif01(rng) * p.DeltaT);
            std::sort(times.begin(), times.end());
            tvec.insert(tvec.end(), times.begin(), times.end());
        }

        // spins: random ±1, then alternate at each boundary
        int s = coin(rng) ? 1 : -1;
        for (size_t i = 0; i < tvec.size(); ++i) {
            svec.push_back(s);
            s = -s;
        }

        // Periodic-time BC: enforce odd number of segments (your original behavior)
        if (p.BC == 1 && (tvec.size() % 2 == 0)) {
            tvec.pop_back();
            svec.pop_back();
            if (tvec.empty()) { // keep at least one segment
                tvec.push_back(0.0);
                svec.push_back(1);
            }
        }
    }
}

// ------------------------------------------------------------
// splitRails: for each base segment [t0,t1), sample Poisson(Gamma*Δ) cuts
// and add uniform cut times, sorted.
// ------------------------------------------------------------
static void splitRails(std::mt19937_64& rng, Sim& p) {
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    for (int x = 0; x < p.L; ++x) {
        std::vector<double>& at = p.aug_time[(size_t)x];
        std::vector<int>& as = p.aug_spin[(size_t)x];
        at.clear();
        as.clear();

        const std::vector<double>& tvec = p.time[(size_t)x];
        const std::vector<int>& svec = p.spin[(size_t)x];

        for (size_t i = 0; i < tvec.size(); ++i) {
            double t0 = tvec[i];
            double t1 = (i + 1 < tvec.size()) ? tvec[i + 1] : p.DeltaT;
            int s = svec[i];

            // always include base segment start
            at.push_back(t0);
            as.push_back(s);

            double dT = t1 - t0;
            if (dT <= 0.0) continue;

            std::poisson_distribution<int> pois(p.Gamma * dT);
            int M = pois(rng);
            if (M > 0) {
                std::vector<double> cuts;
                cuts.reserve((size_t)M);
                for (int k = 0; k < M; ++k) {
                    double u = unif01(rng);
                    cuts.push_back(t0 + u * dT);
                }
                std::sort(cuts.begin(), cuts.end());
                at.insert(at.end(), cuts.begin(), cuts.end());
                as.insert(as.end(), (size_t)M, s);
            }
        }
    }
}

// ------------------------------------------------------------
// clearStaleFlips: compress augmented rails back to base rails by removing
// consecutive identical spins.
// ------------------------------------------------------------
static void clearStaleFlips(Sim& p) {
    for (int x = 0; x < p.L; ++x) {
        const std::vector<double>& at = p.aug_time[(size_t)x];
        const std::vector<int>& as = p.aug_spin[(size_t)x];

        std::vector<double>& tvec = p.time[(size_t)x];
        std::vector<int>& svec = p.spin[(size_t)x];
        tvec.clear();
        svec.clear();

        if (at.empty()) continue;

        tvec.push_back(at[0]);
        svec.push_back(as[0]);
        for (size_t i = 1; i < at.size(); ++i) {
            if (as[i] != svec.back()) {
                tvec.push_back(at[i]);
                svec.push_back(as[i]);
            }
        }
    }
}

// ------------------------------------------------------------
// make_cluster_with_graph: SW-like FK clustering by temporal overlap
// ------------------------------------------------------------
static void make_cluster_with_graph(std::mt19937_64& rng, Sim& p) {
    const int L = p.L;
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    // reset labels
    for (int x = 0; x < L; ++x) {
        p.label[(size_t)x].assign(p.aug_time[(size_t)x].size(), 0);
    }

    int cluster_number = 0;

    struct Seed { int x; int t; };
    auto find_seed = [&](void) -> Seed {
        Seed s; s.x = -1; s.t = -1;
        for (int x = 0; x < L; ++x) {
            for (size_t t = 0; t < p.label[(size_t)x].size(); ++t) {
                if (p.label[(size_t)x][t] == 0) {
                    s.x = x; s.t = (int)t;
                    return s;
                }
            }
        }
        return s;
    };

    while (true) {
        Seed seed = find_seed();
        int seed_x = seed.x;
        int seed_t = seed.t;
        if (seed_x < 0) break;

        if (p.aug_time[(size_t)seed_x].empty()) {
            p.label[(size_t)seed_x].clear();
            continue;
        }

        ++cluster_number;
        int root_spin = p.aug_spin[(size_t)seed_x][(size_t)seed_t];

        std::vector<std::pair<int, int> > st;
        st.reserve(256);

        p.label[(size_t)seed_x][(size_t)seed_t] = cluster_number;
        st.push_back(std::make_pair(seed_x, seed_t));

        // periodic wrap: if label[0] then label[last]
        auto maybe_label_wrap = [&](int x, int tindex) {
            if (p.BC == 1 && p.aug_time[(size_t)x].size() > 1 && tindex == 0) {
                int last = (int)p.aug_time[(size_t)x].size() - 1;
                if (p.label[(size_t)x][(size_t)last] == 0) {
                    p.label[(size_t)x][(size_t)last] = cluster_number;
                    st.push_back(std::make_pair(x, last));
                }
            }
        };
        maybe_label_wrap(seed_x, seed_t);

        while (!st.empty()) {
            std::pair<int, int> cur = st.back();
            st.pop_back();
            int cur_x = cur.first;
            int cur_t = cur.second;

            int size_x = (int)p.aug_time[(size_t)cur_x].size();
            double time1x = p.aug_time[(size_t)cur_x][(size_t)cur_t];
            double time2x = (cur_t == size_x - 1) ? p.DeltaT
                : p.aug_time[(size_t)cur_x][(size_t)(cur_t + 1)];

            const std::vector<int>& nbrs = p.neighbors[(size_t)cur_x];
            for (size_t ny = 0; ny < nbrs.size(); ++ny) {
                int rail_y = nbrs[ny];

                int size_y = (int)p.aug_time[(size_t)rail_y].size();
                if (size_y == 0) continue;

                int max_y = (p.BC == 1 && size_y > 1) ? (size_y - 1) : size_y;

                for (int ty = 0; ty < max_y; ++ty) {
                    double time1y = p.aug_time[(size_t)rail_y][(size_t)ty];
                    double time2y = (ty == size_y - 1) ? p.DeltaT
                        : p.aug_time[(size_t)rail_y][(size_t)(ty + 1)];

                    double overlap = 0.0;
                    if (p.BC == 0) {
                        double lo = std::max(time1x, time1y);
                        double hi = std::min(time2x, time2y);
                        if (hi > lo) overlap = hi - lo;
                    }
                    else {
                        // wrap images k=-1..1
                        for (int k = -1; k <= 1; ++k) {
                            double bt0 = time1y + k * p.DeltaT;
                            double bt1 = time2y + k * p.DeltaT;
                            double lo = std::max(time1x, bt0);
                            double hi = std::min(time2x, bt1);
                            if (hi > lo) overlap += (hi - lo);
                        }
                    }

                    if (overlap > 0.0 &&
                        p.label[(size_t)rail_y][(size_t)ty] == 0 &&
                        p.aug_spin[(size_t)cur_x][(size_t)cur_t] == root_spin &&
                        p.aug_spin[(size_t)rail_y][(size_t)ty] == root_spin)
                    {
                        double p_bond = 1.0 - std::exp(-2.0 * p.K * overlap);
                        if (unif01(rng) < p_bond) {
                            p.label[(size_t)rail_y][(size_t)ty] = cluster_number;
                            st.push_back(std::make_pair(rail_y, ty));

                            // wrap join for ty==0
                            if (p.BC == 1 && size_y > 1 && ty == 0) {
                                int last = size_y - 1;
                                if (p.label[(size_t)rail_y][(size_t)last] == 0) {
                                    p.label[(size_t)rail_y][(size_t)last] = cluster_number;
                                    st.push_back(std::make_pair(rail_y, last));
                                }
                            }
                        }
                    }
                }
            }
        }

        // flip whole cluster with prob 1/2
        if (unif01(rng) < 0.5) {
            for (int x = 0; x < L; ++x) {
                for (size_t i = 0; i < p.label[(size_t)x].size(); ++i) {
                    if (p.label[(size_t)x][i] == cluster_number) {
                        p.aug_spin[(size_t)x][i] *= -1;
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------
// Improved estimators from FK weights
// ------------------------------------------------------------
static inline double total_measure(const Sim& p) {
    return double(p.L) * double(p.DeltaT);
}

static std::vector<double> cluster_weights(const Sim& p) {
    int maxlab = 0;
    for (int x = 0; x < p.L; ++x) {
        const std::vector<int>& labs = p.label[(size_t)x];
        for (size_t i = 0; i < labs.size(); ++i) {
            if (labs[i] > maxlab) maxlab = labs[i];
        }
    }
    if (maxlab == 0) return std::vector<double>();

    std::vector<double> w((size_t)maxlab + 1, 0.0); // 1..maxlab used
    const double DT = p.DeltaT;

    for (int x = 0; x < p.L; ++x) {
        const std::vector<double>& at = p.aug_time[(size_t)x];
        const std::vector<int>& labs = p.label[(size_t)x];
        int nseg = (int)at.size();
        for (int i = 0; i < nseg; ++i) {
            int lab = labs[(size_t)i];
            if (lab == 0) continue;
            double t0 = at[(size_t)i];
            double t1 = (i + 1 < nseg) ? at[(size_t)(i + 1)] : DT;
            w[(size_t)lab] += (t1 - t0);
        }
    }

    w.erase(w.begin()); // drop index 0
    return w;
}

static inline std::pair<double, double> improved_moments_from_weights(const std::vector<double>& w, double N) {
    double s2 = 0.0;
    double s4 = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        double wf = w[i];
        double wf2 = wf * wf;
        s2 += wf2;
        s4 += wf2 * wf2;
    }
    double N2 = N * N;
    double N4 = N2 * N2;
    double m2_hat = (N2 > 0.0) ? (s2 / N2) : 0.0;
    double m4_hat = (N4 > 0.0) ? ((3.0 * s2 * s2 - 2.0 * s4) / N4) : 0.0;
    return std::make_pair(m2_hat, m4_hat);
}

static inline std::pair<double, double> jackknife_U(
    const std::vector<double>& m2,
    const std::vector<double>& m4,
    bool scaled = true)
{
    const int n = (int)m2.size();
    if (n <= 1 || (int)m4.size() != n) return std::make_pair(0.0, 0.0);

    double s2 = std::accumulate(m2.begin(), m2.end(), 0.0);
    double s4 = std::accumulate(m4.begin(), m4.end(), 0.0);

    std::vector<double> Uvals((size_t)n, 0.0);
    for (int i = 0; i < n; ++i) {
        double m2i = (s2 - m2[(size_t)i]) / double(n - 1);
        double m4i = (s4 - m4[(size_t)i]) / double(n - 1);
        double Ui = 0.0;
        if (m2i != 0.0) Ui = 1.0 - m4i / (3.0 * m2i * m2i);
        Uvals[(size_t)i] = scaled ? (1.5 * Ui) : Ui;
    }

    double Ubar = std::accumulate(Uvals.begin(), Uvals.end(), 0.0) / double(n);
    double mean_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = Uvals[(size_t)i] - Ubar;
        mean_sq += d * d;
    }
    mean_sq /= double(n);
    double sigma = std::sqrt(double(n - 1) * mean_sq);
    return std::make_pair(Ubar, sigma);
}

// ------------------------------------------------------------
// 1D chain neighbors (periodic spatial BC)
// ------------------------------------------------------------
static std::vector<std::vector<int>> chain_neighbors(int L) {
    std::vector<std::vector<int>> neighbors((size_t)L);
    for (int x = 0; x < L; ++x) {
        int left = (x == 0) ? (L - 1) : (x - 1);
        int right = (x == L - 1) ? 0 : (x + 1);
        neighbors[(size_t)x] = { left, right };
    }
    return neighbors;
}

// ------------------------------------------------------------
// Binder scan on a fixed graph (CSV rows)
// ------------------------------------------------------------
struct Row {
    int L;
    double ratio;
    double Gamma;
    double K;
    double m2, m2_err;
    double m4, m4_err;
    double U, U_err;
    int n;
    int sweeps_between;
    int burn_in;
    int n_layers;
    int N_total;
    int N_boundary;
};

// NEW: add std::ostream* ts_out (nullptr disables)
static std::vector<Row> run_binder_improved_with_graph(
    const std::vector<std::vector<int> >& neighbors,
    const std::vector<double>& ratios,
    uint64_t seed,
    int BC,
    double baseDeltaT,
    int burn_in,
    int sweeps_between,
    int n_configs,
    bool scaled,
    int n_layers,              // for progress prints
    int q,                     // for progress prints (unused for chain)
    int print_every_burn,
    int print_every_cfg,
    std::ostream* ts_out       // NEW
) {
    (void)q;
    std::mt19937_64 rng(seed);
    const int L = (int)neighbors.size();

    std::vector<Row> rows;
    rows.reserve(ratios.size());

    for (size_t job = 0; job < ratios.size(); ++job) {
        double ratio = ratios[job];
        double Gamma = ratio;
        double K = 1.0;

        double DeltaT = baseDeltaT * double(L);
        Sim p(L, neighbors, DeltaT, Gamma, K, BC);

        {
            std::ostringstream oss;
            oss << "[job] L=" << L
                << " (" << (job + 1) << "/" << ratios.size() << ")"
                << " ratio=" << ratio
                << " Gamma=" << Gamma
                << " K=" << K
                << " DeltaT=" << DeltaT
                << " BC=" << BC
                << " burn_in=" << burn_in
                << " sweeps_between=" << sweeps_between
                << " n_configs=" << n_configs;
            progress_line(oss.str());
        }

        auto t_job0 = std::chrono::steady_clock::now();

        // init + burn-in
        makeRails(rng, p);
        splitRails(rng, p);
        clearStaleFlips(p);

        for (int i = 0; i < burn_in; ++i) {
            splitRails(rng, p);
            make_cluster_with_graph(rng, p);
            clearStaleFlips(p);
            progress_tick(i, burn_in, print_every_burn, "  [burn]");
        }

        const double Ntot = total_measure(p);
        std::vector<double> m2_list; m2_list.reserve((size_t)n_configs);
        std::vector<double> m4_list; m4_list.reserve((size_t)n_configs);

        for (int cfg = 0; cfg < n_configs; ++cfg) {
            for (int sw = 0; sw < sweeps_between; ++sw) {
                splitRails(rng, p);
                make_cluster_with_graph(rng, p);
                clearStaleFlips(p);
            }

            if (print_every_cfg > 0 && (cfg == 0 || cfg == n_configs - 1 || (cfg % print_every_cfg) == 0)) {
                long long segs = 0;
                for (int x = 0; x < p.L; ++x) segs += (long long)p.time[(size_t)x].size();
                double avg_segs = (p.L > 0) ? (double(segs) / double(p.L)) : 0.0;

                auto now = std::chrono::steady_clock::now();
                double elapsed_s = std::chrono::duration<double>(now - t_job0).count();

                std::ostringstream oss;
                oss << "  [cfg] " << (cfg + 1) << "/" << n_configs
                    << " avg_segments_per_site=" << std::fixed << std::setprecision(2) << avg_segs
                    << " elapsed_s=" << std::setprecision(2) << elapsed_s;
                progress_line(oss.str());
            }

            // measurement
            std::vector<double> w = cluster_weights(p);
            double m2hat = 0.0, m4hat = 0.0;
            if (!w.empty()) {
                auto mm = improved_moments_from_weights(w, Ntot);
                m2hat = mm.first;
                m4hat = mm.second;
            }
            m2_list.push_back(m2hat);
            m4_list.push_back(m4hat);

            // NEW: time-series write (one line per cfg)
            if (ts_out) {
                double Ucfg = 0.0;
                if (m2hat != 0.0) Ucfg = 1.0 - m4hat / (3.0 * m2hat * m2hat);
                if (scaled) Ucfg *= 1.5;

                long long segs = 0;
                for (int x = 0; x < p.L; ++x) segs += (long long)p.time[(size_t)x].size();
                double avg_segs = (p.L > 0) ? (double(segs) / double(p.L)) : 0.0;

                // Columns: L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site
                // mc_step is measured in "sweeps" (cfg*sweeps_between)
                (*ts_out)
                    << L << "," << ratio << ","
                    << cfg << "," << ((long long)cfg * (long long)sweeps_between) << ","
                    << std::setprecision(17)
                    << m2hat << "," << m4hat << "," << Ucfg << ","
                    << avg_segs
                    << "\n";
            }
        }

        double mean_m2 = std::accumulate(m2_list.begin(), m2_list.end(), 0.0) / double(m2_list.size());
        double mean_m4 = std::accumulate(m4_list.begin(), m4_list.end(), 0.0) / double(m4_list.size());
        double se_m2 = sem(m2_list);
        double se_m4 = sem(m4_list);

        auto UU = jackknife_U(m2_list, m4_list, scaled);

        Row r;
        r.L = L;
        r.ratio = ratio;
        r.Gamma = Gamma;
        r.K = K;
        r.m2 = mean_m2;
        r.m2_err = se_m2;
        r.m4 = mean_m4;
        r.m4_err = se_m4;
        r.U = UU.first;
        r.U_err = UU.second;
        r.n = (int)m2_list.size();
        r.sweeps_between = sweeps_between;
        r.burn_in = burn_in;
        r.n_layers = n_layers;
        r.N_total = L;
        r.N_boundary = 0;

        rows.push_back(r);

        {
            auto now = std::chrono::steady_clock::now();
            double elapsed_s = std::chrono::duration<double>(now - t_job0).count();
            std::ostringstream oss;
            oss << "[job_done] L=" << L
                << " ratio=" << ratio
                << " U=" << std::setprecision(6) << r.U
                << " +/- " << r.U_err
                << " elapsed_s=" << std::setprecision(2) << elapsed_s;
            progress_line(oss.str());
        }
    }

    return rows;
}

// ------------------------------------------------------------
// CLI
// ------------------------------------------------------------
struct Args {
    int q = 7;
    std::vector<int> layers; // interpreted as chain lengths L
    double ratio_min = 1.2;
    double ratio_max = 0.8;
    int ratio_n = 10;
    int burn_in = 100;
    int sweeps_between = 20;
    int n_configs = 1000;
    uint64_t seed = 1234;
    int BC = 1;
    bool scaled = true;

    int print_every_burn = 50;
    int print_every_cfg = 25;

    // NEW: time-series output file (empty => disabled)
    std::string ts_file = "ts.dat";

    Args() { layers = { 32, 28, 24 }; }
};

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                std::exit(2);
            }
            return std::string(argv[++i]);
        };

        if (k == "--q") a.q = std::stoi(need("--q"));
        else if (k == "--layers") a.layers = parse_int_list(need("--layers"));
        else if (k == "--ratio_min") a.ratio_min = std::stod(need("--ratio_min"));
        else if (k == "--ratio_max") a.ratio_max = std::stod(need("--ratio_max"));
        else if (k == "--ratio_n") a.ratio_n = std::stoi(need("--ratio_n"));
        else if (k == "--burn_in") a.burn_in = std::stoi(need("--burn_in"));
        else if (k == "--sweeps_between") a.sweeps_between = std::stoi(need("--sweeps_between"));
        else if (k == "--n_configs") a.n_configs = std::stoi(need("--n_configs"));
        else if (k == "--seed") a.seed = (uint64_t)std::stoull(need("--seed"));
        else if (k == "--BC") a.BC = std::stoi(need("--BC"));
        else if (k == "--scaled") a.scaled = (std::stoi(need("--scaled")) != 0);
        else if (k == "--print_every_burn") a.print_every_burn = std::stoi(need("--print_every_burn"));
        else if (k == "--print_every_cfg")  a.print_every_cfg = std::stoi(need("--print_every_cfg"));

        // NEW:
        else if (k == "--ts_file") a.ts_file = need("--ts_file");

        else if (k == "--help" || k == "-h") {
            std::cout <<
                "ConTimeFromJulia options:\n"
                "  --layers a,b,c   (comma-separated chain lengths L)\n"
                "  --ratio_min DOUBLE\n"
                "  --ratio_max DOUBLE\n"
                "  --ratio_n INT\n"
                "  --burn_in INT\n"
                "  --sweeps_between INT\n"
                "  --n_configs INT\n"
                "  --seed UINT64\n"
                "  --BC 0|1\n"
                "  --scaled 0|1\n"
                "  --ts_file PATH   (write per-config time series CSV; empty disables)\n"
                "  --print_every_burn INT\n"
                "  --print_every_cfg INT\n";
            std::exit(0);
        }
        else {
            std::cerr << "Unknown arg: " << k << "\n";
            std::exit(2);
        }
    }
    return a;
}

// ------------------------------------------------------------
// main: 1D chain in space, continuous time in DeltaT = baseDeltaT*L
// ------------------------------------------------------------
int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);
    std::vector<double> ratios = linspace(args.ratio_min, args.ratio_max, args.ratio_n);

    // Summary CSV header (STDOUT)
    std::cout
        << "L,ratio,Gamma,K,m2,m2_err,m4,m4_err,U,U_err,n,sweeps_between,burn_in\n";
    std::cout << std::setprecision(17);

    // NEW: time-series file
    std::ofstream ts;
    std::ostream* ts_out = nullptr;
    if (!args.ts_file.empty()) {
        ts.open(args.ts_file.c_str());
        if (!ts) {
            std::cerr << "ERROR: could not open --ts_file " << args.ts_file << "\n";
            return 2;
        }
        ts << "L,ratio,cfg,mc_step,m2,m4,U,avg_segments_per_site\n";
        ts_out = &ts;
    }

    for (int L : args.layers) {
        progress_line("============================================================");
        progress_line("[chain] L=" + std::to_string(L) + "  BC=" + std::to_string(args.BC));

        std::vector<std::vector<int> > neighbors = chain_neighbors(L);

        // beta = DeltaT = baseDeltaT * L. With baseDeltaT=1 => beta=L.
        double baseDeltaT = 1.0;

        std::vector<Row> rows = run_binder_improved_with_graph(
            neighbors,
            ratios,
            args.seed,
            args.BC,
            baseDeltaT,
            args.burn_in,
            args.sweeps_between,
            args.n_configs,
            args.scaled,
            /*n_layers=*/L,
            /*q=*/0,
            args.print_every_burn,
            args.print_every_cfg,
            ts_out // NEW
        );

        for (const Row& r : rows) {
            std::cout
                << L << "," << r.ratio << "," << r.Gamma << "," << r.K << ","
                << r.m2 << "," << r.m2_err << ","
                << r.m4 << "," << r.m4_err << ","
                << r.U << "," << r.U_err << ","
                << r.n << "," << r.sweeps_between << "," << r.burn_in
                << "\n";
        }
    }

    return 0;
}
