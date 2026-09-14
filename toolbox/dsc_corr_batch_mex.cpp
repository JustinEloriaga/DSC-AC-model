// Batched Archakov-Hansen inverse map and Gaussian correlation likelihood.
// [ll,P,L,logdet,iterations] = dsc_corr_batch_mex(r,Z,mask,threads)
// threads is optional, defaults to 1, and must be an integer in [1,64].
// r: T-by-m*(m-1)/2, column-major strict-lower-triangle coordinates.
// Z: T-by-m standardized residuals; mask: T-by-m logical observations.
// P and L: m-by-m-by-T. L contains each observed-space lower Cholesky
// factor in its leading nobs-by-nobs block; all other entries are zero.
// ll excludes Gaussian constants and log-volatility contributions.
// Optional std::thread workers process static contiguous date blocks;
// LAPACK/BLAS thread settings and MATLAB global state are never changed.
#include "mex.h"
#include "lapack.h"
#include "blas.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
struct KernelError : std::runtime_error {
    const char* identifier;
    KernelError(const char* id, const std::string& message)
        : std::runtime_error(message), identifier(id) {}
};

struct JoinThreads {
    std::vector<std::thread>& threads;
    explicit JoinThreads(std::vector<std::thread>& value) : threads(value) {}
    ~JoinThreads() {
        for (auto& worker : threads) if (worker.joinable()) worker.join();
    }
};

void require(bool condition, const char* message) {
    if (!condition) throw KernelError("dsc:KernelInput", message);
}

void correlation_error(const char* identifier, const char* message, mwSize t) {
    throw KernelError(identifier, std::string(message) + " (week " +
                      std::to_string(t + 1) + ").");
}

bool real_full_double(const mxArray* value) {
    return mxIsDouble(value) && !mxIsComplex(value) && !mxIsSparse(value) &&
           mxGetNumberOfDimensions(value) == 2;
}

void evaluate(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nrhs == 3 || nrhs == 4,
            "Expected inputs r, Z, a logical observation mask and optional threads.");
    require(nlhs <= 5, "At most five outputs are available.");
    require(real_full_double(prhs[0]) && real_full_double(prhs[1]),
            "r and Z must be real, full, two-dimensional double arrays.");
    require(mxIsLogical(prhs[2]) && !mxIsSparse(prhs[2]) &&
            mxGetNumberOfDimensions(prhs[2]) == 2,
            "mask must be a full, two-dimensional logical array.");
    mwSize requested_threads = 1;
    if (nrhs == 4) {
        require(real_full_double(prhs[3]) && mxGetNumberOfElements(prhs[3]) == 1,
                "threads must be a double scalar integer in [1,64].");
        const double value = mxGetScalar(prhs[3]);
        require(std::isfinite(value) && value >= 1 && value <= 64 && value == std::floor(value),
                "threads must be a double scalar integer in [1,64].");
        requested_threads = static_cast<mwSize>(value);
    }
    const mwSize T = mxGetM(prhs[1]);
    const mwSize m = mxGetN(prhs[1]);
    require(m >= 1, "Z must contain at least one series.");
    require(m <= static_cast<mwSize>(std::numeric_limits<ptrdiff_t>::max()) &&
            m <= std::numeric_limits<mwSize>::max() / m,
            "The matrix dimension exceeds the native indexing limit.");
    const mwSize mm = m * m;
    const mwSize q = m * (m - 1) / 2;
    require(mxGetM(prhs[0]) == T && mxGetN(prhs[0]) == q,
            "r must be T-by-m*(m-1)/2, consistent with Z.");
    require(mxGetM(prhs[2]) == T && mxGetN(prhs[2]) == m,
            "mask must have the same dimensions as Z.");
    require(T == 0 || mm <= std::numeric_limits<mwSize>::max() / T,
            "The requested batch exceeds the native indexing limit.");
    const double* input_r = mxGetDoubles(prhs[0]);
    const double* input_Z = mxGetDoubles(prhs[1]);
    const mxLogical* input_mask = mxGetLogicals(prhs[2]);
    for (mwSize k = 0; k < mxGetNumberOfElements(prhs[0]); ++k) {
        require(std::isfinite(input_r[k]), "All correlation coordinates must be finite.");
    }
    for (mwSize k = 0; k < mxGetNumberOfElements(prhs[1]); ++k) {
        require(!input_mask[k] || std::isfinite(input_Z[k]),
                "Every observed standardized residual must be finite.");
    }

    double* Pout = nullptr;
    double* Lout = nullptr;
    double* logdet_out = nullptr;
    double* iterations_out = nullptr;
    const mwSize dimensions[3] = {m, m, T};
    if (nlhs >= 2) {
        plhs[1] = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS, mxREAL);
        Pout = mxGetDoubles(plhs[1]);
    }
    if (nlhs >= 3) {
        plhs[2] = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS, mxREAL);
        Lout = mxGetDoubles(plhs[2]);
    }
    if (nlhs >= 4) {
        plhs[3] = mxCreateDoubleMatrix(T, 1, mxREAL);
        logdet_out = mxGetDoubles(plhs[3]);
    }
    if (nlhs >= 5) {
        plhs[4] = mxCreateDoubleMatrix(T, 1, mxREAL);
        iterations_out = mxGetDoubles(plhs[4]);
    }
    if (T == 0) {
        if (nlhs >= 1) plhs[0] = mxCreateDoubleScalar(0.0);
        return;
    }

    // Copy all inputs before starting workers. Workers invoke no mx/mex API
    // and write only to disjoint, already allocated output-array segments.
    std::vector<double> r(T*q), Z(T*m);
    std::vector<unsigned char> mask(T*m);
    if (q != 0) std::copy(input_r, input_r + T*q, r.begin());
    std::copy(input_Z, input_Z + T*m, Z.begin());
    std::copy(input_mask, input_mask + T*m, mask.begin());
    std::vector<double> date_likelihood(T, 0.0);

    auto process_dates = [&](mwSize begin_date, mwSize end_date) {
    const ptrdiff_t n = static_cast<ptrdiff_t>(m);
    const char vectors = 'V', lower = 'L', no_transpose = 'N', nonunit = 'N';
    const ptrdiff_t increment = 1;
    ptrdiff_t info = 0;
    std::vector<double> base(mm, 0.0), eigenvectors(mm), eigenvalues(m),
        exp_eigenvalues(m), diagonal(m, 0.0), C(mm), factor(mm), z(m);
    std::vector<mwSize> observed(m);
    double optimal_workspace = 0.0;
    ptrdiff_t query = -1;
    dsyev(&vectors, &lower, &n, eigenvectors.data(), &n, eigenvalues.data(),
          &optimal_workspace, &query, &info);
    if (info != 0 || !std::isfinite(optimal_workspace) || optimal_workspace < 1 ||
        optimal_workspace >= static_cast<double>(std::numeric_limits<ptrdiff_t>::max())) {
        throw KernelError("dsc:CorrelationNumerical", "LAPACK workspace query failed.");
    }
    const ptrdiff_t lwork = static_cast<ptrdiff_t>(optimal_workspace);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    const double tolerance = 1e-8 * std::sqrt(static_cast<double>(m));

    for (mwSize t = begin_date; t < end_date; ++t) {
        mwSize coordinate = 0;
        for (mwSize j = 0; j < m; ++j) {
            for (mwSize i = j + 1; i < m; ++i, ++coordinate) {
                base[i + m*j] = r[t + T*coordinate];
                base[j + m*i] = base[i + m*j];
            }
        }
        // Each block's first date starts at zero; remaining dates use the
        // preceding solved diagonal. One thread exactly retains the original
        // across-date warm starts. Block boundaries only change initialization.
        auto exponential_eigensystem = [&]() {
            std::copy(base.begin(), base.end(), eigenvectors.begin());
            for (mwSize i = 0; i < m; ++i) eigenvectors[i + m*i] = diagonal[i];
            dsyev(&vectors, &lower, &n, eigenvectors.data(), &n,
                  eigenvalues.data(), work.data(), &lwork, &info);
            if (info != 0) {
                correlation_error("dsc:CorrelationNumerical", "Symmetric eigensolver failed", t);
            }
            for (mwSize k = 0; k < m; ++k) {
                exp_eigenvalues[k] = std::exp(eigenvalues[k]);
                if (!std::isfinite(exp_eigenvalues[k])) {
                    correlation_error("dsc:CorrelationNumerical", "Nonfinite correlation transform", t);
                }
            }
        };

        bool converged = false;
        int iteration = 0;
        for (iteration = 1; iteration <= 200; ++iteration) {
            exponential_eigensystem();
            double residual_squared = 0.0;
            // diag(Q*diag(exp(lambda))*Q') = (Q.^2)*exp(lambda).
            // Avoid materializing the full exponential inside this loop.
            for (mwSize i = 0; i < m; ++i) {
                double value = 0.0;
                for (mwSize k = 0; k < m; ++k) {
                    const double Qik = eigenvectors[i + m*k];
                    value += Qik * Qik * exp_eigenvalues[k];
                }
                const double delta = std::log(value);
                if (!std::isfinite(delta)) {
                    correlation_error("dsc:CorrelationNumerical", "Nonfinite correlation transform", t);
                }
                diagonal[i] -= delta;
                residual_squared += delta * delta;
            }
            if (std::sqrt(residual_squared) < tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            correlation_error("dsc:CorrelationConvergence",
                              "Correlation inverse did not converge after 200 iterations", t);
        }
        // Recompute after the final diagonal update exactly as the MATLAB
        // reference does, then enforce symmetry and unit diagonal.
        exponential_eigensystem();
        for (mwSize j = 0; j < m; ++j) {
            C[j + m*j] = 1.0;
            for (mwSize i = j + 1; i < m; ++i) {
                double value = 0.0;
                for (mwSize k = 0; k < m; ++k) {
                    value += eigenvectors[i + m*k] * exp_eigenvalues[k] *
                             eigenvectors[j + m*k];
                }
                if (!std::isfinite(value)) {
                    correlation_error("dsc:CorrelationNumerical", "Nonfinite correlation matrix", t);
                }
                C[i + m*j] = value;
                C[j + m*i] = value;
            }
        }
        std::copy(C.begin(), C.end(), factor.begin());
        dpotrf(&lower, &n, factor.data(), &n, &info);
        if (info != 0) {
            correlation_error("dsc:CorrelationPD",
                              "Correlation matrix is not numerically positive definite", t);
        }
        if (Pout) std::copy(C.begin(), C.end(), Pout + mm*t);
        if (iterations_out) iterations_out[t] = static_cast<double>(iteration);

        mwSize nobs = 0;
        for (mwSize i = 0; i < m; ++i) {
            if (mask[t + T*i]) observed[nobs++] = i;
        }
        double logdet = 0.0;
        if (nobs != 0) {
            const ptrdiff_t no = static_cast<ptrdiff_t>(nobs);
            if (nobs != m) {
                for (mwSize j = 0; j < nobs; ++j) {
                    for (mwSize i = 0; i < nobs; ++i) {
                        factor[i + m*j] = C[observed[i] + m*observed[j]];
                    }
                }
                dpotrf(&lower, &no, factor.data(), &n, &info);
                if (info != 0) {
                    correlation_error("dsc:CorrelationPD",
                                      "Observed correlation matrix is not positive definite", t);
                }
            }
            for (mwSize i = 0; i < nobs; ++i) {
                logdet += 2.0 * std::log(factor[i + m*i]);
                z[i] = Z[t + T*observed[i]];
            }
            dtrsv(&lower, &no_transpose, &nonunit, &no, factor.data(), &n,
                  z.data(), &increment);
            double quadratic = 0.0;
            for (mwSize i = 0; i < nobs; ++i) quadratic += z[i] * z[i];
            date_likelihood[t] = -0.5 * (logdet + quadratic);
            if (Lout) {
                for (mwSize j = 0; j < nobs; ++j) {
                    for (mwSize i = j; i < nobs; ++i) {
                        Lout[mm*t + i + m*j] = factor[i + m*j];
                    }
                }
            }
        }
        if (logdet_out) logdet_out[t] = logdet;
    }
    };

    const mwSize worker_count = std::min(requested_threads, T);
    if (worker_count == 1) {
        process_dates(0, T);
    } else {
        std::vector<std::exception_ptr> failures(worker_count);
        auto run_worker = [&](mwSize worker) {
            try {
                const mwSize base_length = T / worker_count;
                const mwSize remainder = T % worker_count;
                const mwSize first = worker * base_length + std::min(worker, remainder);
                process_dates(first, first + base_length + (worker < remainder ? 1 : 0));
            } catch (...) {
                failures[worker] = std::current_exception();
            }
        };
        {
            std::vector<std::thread> workers;
            workers.reserve(worker_count - 1);
            JoinThreads joiner(workers);
            for (mwSize worker = 1; worker < worker_count; ++worker) {
                workers.emplace_back(run_worker, worker);
            }
            run_worker(0); // The calling thread is included in the requested count.
        }
        for (const auto& failure : failures) if (failure) std::rethrow_exception(failure);
    }
    // Fixed date-order reduction is deterministic for a fixed thread count,
    // regardless of worker scheduling. Repeated calls do not retain state.
    double ll = 0.0;
    for (mwSize t = 0; t < T; ++t) ll += date_likelihood[t];
    if (!std::isfinite(ll)) ll = -mxGetInf();
    if (nlhs >= 1) plhs[0] = mxCreateDoubleScalar(ll);
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Unwind C++ work buffers before entering MATLAB's nonlocal error exit.
    try {
        evaluate(nlhs, plhs, nrhs, prhs);
    } catch (const KernelError& problem) {
        mexErrMsgIdAndTxt(problem.identifier, "%s", problem.what());
    } catch (const std::bad_alloc&) {
        mexErrMsgIdAndTxt("dsc:KernelMemory", "Unable to allocate native kernel workspace.");
    } catch (const std::exception& problem) {
        mexErrMsgIdAndTxt("dsc:KernelInternal", "%s", problem.what());
    }
}
