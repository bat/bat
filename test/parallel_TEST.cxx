/*
 * Copyright (C) 2012, Sankalp Sardana and Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <config.h>

#if THREAD_PARALLELIZATION

#include <GaussModel.h>
#include <test.h>

#include <BAT/BCAux.h>
#include <BAT/BCEmptyModel.h>
#include <BAT/BCLog.h>

#include <TFile.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TTree.h>

#include <omp.h>
#include <cstdio>

using namespace test;

namespace
{
class Output
{

public:

    /**
     * Print speed up/slow down of serial vs parallel execution.
     * Prints only on every second call.
     *
     * @param time The time of the last execution.
     */
    void PrintValues(const double& time)
    {
        static std::vector<double> real_time_array(0);

        real_time_array.push_back(time);

        // only print when two new values have been stored
        if ((real_time_array.size() % 2) == 1)
            return;

        const double& previous_time = real_time_array.at(real_time_array.size() - 2);
        std::cout << "previous time: " << previous_time << std::endl;

        const double diff =  previous_time - time;
        const double ratio = previous_time / time;

        std::cout << "Difference in time taken between serial and parallel: " << diff << std::endl;
        std::cout << "Time improvement factor: " << ratio << std::endl;
    }
};

struct BCCheckModel : BCEmptyModel {
    BCCheckModel(const std::string& file) :
        // need empty string to call ctor that searches for tree in file
        BCEmptyModel(file, "")
    {
    }

    TTree* Tree() const {return fMCMCTree;}

    double prob() const {return fMCMCTree_Prob;}
    unsigned phase() const {return fMCMCPhase;}
    const std::vector<double>& parameters() const {return fMCMCTree_Parameters;}
};

class RunComparison
{
public:
    struct Config {
        /**
         * The configure option to set the delay time */
        long delay;
        /**
         * Multivariate Student's t degree of freedom.
         */
        double dof;
        /**
         * Multivariate proposal */
        bool multivariate;
        /**
         * The configure option to set #markov chains */
        unsigned num_chains;
        /**
         * The configure option to set the #entries to check in parallel test */
        long num_entries;
        /**
         * The configure option to set the #parameters for the model */
        unsigned num_parameters;
        /**
         * The configure option to set the #iterations for each chain in the model */
        unsigned num_iterations;
        /**
         * The configure option to set the plotting of results */
        bool plot;
        /**
         * The configure option which sets the root filename for serial output */
        std::string rootFileNameSerial;
        /**
         * The configure option which sets the root filename for parallel output */
        std::string rootFileNameParallel;

    private:
        Config():
            delay(1e5),
            dof(0),
            multivariate(false),
            num_chains(4),
            num_entries(300),
            num_parameters(1),
            num_iterations(10000),
            plot(false),
            rootFileNameSerial(BAT_TESTDIR "parallel_TEST_GaussModelSerial.root"),
            rootFileNameParallel(BAT_TESTDIR "parallel_TEST_GaussModelParallel.root")
        {

        }

    public:
        static Config Default()
        {
            return Config();
        }
    };
private:
    Config config;

public:
    std::string GaussModel_plots_Serial;
    std::string GaussModel_plots_Parallel;
    unsigned seed;

    RunComparison(const RunComparison::Config& config) :
        config(config),
        GaussModel_plots_Serial(BAT_TESTDIR "parallel_TEST_GaussModel_plots_Serial.pdf"),
        GaussModel_plots_Parallel(BAT_TESTDIR "parallel_TEST_GaussModel_plots_Parallel.pdf"),
        seed(11)
    {
        CreateOutput(false);
        CreateOutput(true);
        Check();
    }

    void CreateOutput(const bool& parallelization) const
    {
        Output op;

        omp_set_dynamic(0);
        omp_set_num_threads(parallelization ? config.num_chains : 1);

        // open log file
        BCLog::OpenLog("log.txt");
        BCLog::SetLogLevel(BCLog::detail);

        // create new GaussModel object
        GaussModel m(parallelization ? "Parallel evaluation" : "Serial evaluation", config.num_parameters, config.delay);

        // set MCMC precision
        m.SetPrecision(BCEngineMCMC::kMedium);
        m.SetNIterationsRun(config.num_iterations);
        m.SetNChains(config.num_chains);

        // switch writing of Markov Chains on
        m.WriteMarkovChain(parallelization ? config.rootFileNameParallel.c_str() : config.rootFileNameSerial.c_str(), "RECREATE", true);

        m.SetProposeMultivariate(config.multivariate);
        m.SetProposalFunctionDof(config.dof);

        m.SetRandomSeed(seed);

        TStopwatch sw;
        sw.Start();

        // run MCMC and marginalize posterior wrt. all parameters
        // and all combinations of two parameters
        m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        m.MarginalizeAll();
        sw.Stop();
        double real_time = sw.RealTime();
        op.PrintValues(real_time);

        if (config.plot) {
            if (parallelization)
                m.PrintAllMarginalized(GaussModel_plots_Parallel.c_str());
            else
                m.PrintAllMarginalized(GaussModel_plots_Serial.c_str());
        }

        // close log file
        BCLog::CloseLog();

    }

    void Check() throw(TestCaseFailedException)
    {
        // check if files were written
        {
            TFile* rfile0 = TFile::Open(config.rootFileNameParallel.c_str());
            TFile* rfile1 = TFile::Open(config.rootFileNameSerial.c_str());

            if (!rfile0)
                TEST_CHECK_FAILED(std::string("Could not open") + config.rootFileNameParallel);

            if (!rfile1)
                TEST_CHECK_FAILED(std::string("Could not open") + config.rootFileNameSerial);

            // clean up
            rfile0->Close();
            rfile1->Close();
            delete rfile0;
            delete rfile1;
        }

        // reuse BAT's deserialization
        static const unsigned M = 2;
        BCCheckModel models[M] = {
            BCCheckModel(config.rootFileNameParallel),
            BCCheckModel(config.rootFileNameSerial),
        };

        // set all branch addresses
        for (unsigned m = 0; m < M; ++m)
            models[m].BCEngineMCMC::Remarginalize();

        // trees should exist
        if (!models[0].Tree())
            TEST_CHECK_FAILED(std::string("Could not open tree from ") + config.rootFileNameParallel);
        if (!models[1].Tree())
            TEST_CHECK_FAILED(std::string("Could not open tree from ") + config.rootFileNameSerial);

        // length has to match
        long N = models[0].Tree()->GetEntries();
        TEST_CHECK_EQUAL(N, models[1].Tree()->GetEntries());

        // number of parameters has to match
        unsigned npar = models[0].GetNParameters();
        TEST_CHECK_EQUAL(npar, models[1].GetNParameters());


        /* compute mean chi2, include lag for independent samples*/

        unsigned includedSamples = 0;
        const unsigned lag = 20;
        double chi2 = 0.0;
        double meanChi2 = 0.0;

        /* compare data for each iteration */

        for (unsigned n = 0; n < N; ++n) {
            // read in data from iteration n
            for (unsigned m = 0; m < M; ++m)
                models[m].Tree()->GetEntry(n);

            // loop over each dimension
            for (unsigned p = 0; p < npar; ++p) {
                TEST_CHECK_EQUAL(models[0].parameters().at(p), models[1].parameters().at(p));
            }

            TEST_CHECK_EQUAL(models[0].prob(), models[1].prob());
            TEST_CHECK_EQUAL(models[0].phase(), models[1].phase());

            // check only in main run
            if ((models[0].phase() == 1) && (n % lag == 0)) {
                ++includedSamples;
                chi2 = 0;
                for (unsigned p = 0; p < npar; ++p)
                    chi2 += models[0].parameters().at(p) * models[0].parameters().at(p);
                meanChi2 += (chi2 - meanChi2) / includedSamples;
            }
        }

        // clean up output files
        remove(config.rootFileNameParallel.c_str());
        remove(config.rootFileNameSerial.c_str());

        /* check that samples are from a Gaussian */

        // chi2 distribution has std. dev sqrt(2 / npar). Allow 2 sigma range
        TEST_CHECK_NEARLY_EQUAL(meanChi2 / npar, 1.0, 2 * std::sqrt(2.0 / npar));

        for (unsigned c = 0; c < models[0].GetNChains(); ++c) {
            const BCEngineMCMC::Statistics& s = models[0].GetStatistics(c);
            TEST_CHECK_EQUAL(s.n_samples, config.num_iterations);

            // unit Gaussian
            for (unsigned i = 0; i < models[0].GetNParameters(); ++i) {
                TEST_CHECK_NEARLY_EQUAL(s.mean.at(i), 0, 2e-1);
                TEST_CHECK_NEARLY_EQUAL(s.variance.at(i), 1, 5e-1);
                for (unsigned j = i + 1; j < models[0].GetNParameters(); ++j) {
                    TEST_CHECK_NEARLY_EQUAL(s.covariance.at(i).at(j), 0, 5e-1);
                }
            }
        }
    }
};

class ParallelTest :
    public TestCase
{
public:
    ParallelTest():
        TestCase("parallelization")
    {
    }

    virtual void run() const
    {
        /*
         * Run MCMC in serial and parallel and compare the output
         */

        RunComparison::Config config = RunComparison::Config::Default();

        config.num_chains = 2;
        config.num_parameters = 5;
        config.delay = 1e2;

        TEST_SECTION("product proposal", {
            config.multivariate = false;
            RunComparison comparison(config);
        });

        TEST_SECTION("multivariate Gaussian proposal", {
            config.multivariate = true;
            RunComparison comparison(config);
        });
        TEST_SECTION("multivariate Cauchy proposal", {
            config.multivariate = true;
            config.dof = 1.0;
            RunComparison comparison(config);
        });
    }
} parallel_test;
}
#endif
