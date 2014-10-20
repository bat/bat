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

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCModelOutput.h>

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
    void PrintValues(const double & time)
    {
        static std::vector<double> real_time_array(0);

        real_time_array.push_back(time);

        // only print when two new values have been stored
        if((real_time_array.size() % 2) == 1)
            return;

        const double & previous_time = real_time_array.at(real_time_array.size() - 2);
        std::cout << "previous time: " << previous_time << std::endl;

        const double diff =  previous_time - time;
        const double ratio = previous_time / time;

        std::cout << "Difference in time taken between serial and parallel: " << diff << std::endl;
        std::cout << "Time improvement factor: " << ratio << std::endl;
    }
};

class RunComparison
{
public:
    struct Config{
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
         * The configure option to set the lag time */
        long lag;
        /**
         * The configure option which sets the root filename for serial output */
        std::string rootFileNameSerial;
        /**
         * The configure option which sets the root filename for parallel output */
        std::string rootFileNameParallel;

    private:
        Config():
            num_chains(4),
            num_entries(300),
            num_parameters(1),
            num_iterations(100),
            plot(false),
            lag(1e5),
            rootFileNameSerial(BAT_TESTDIR "parallel_TEST_GaussModelSerial.root"),
            rootFileNameParallel(BAT_TESTDIR "parallel_TEST_GaussModelParallel.root")
        {

        }

    public:
        static Config Default(){
            return Config();
        }
    };
private:
    Config config;

public:
    std::string GaussModel_plots_Serial;
    std::string GaussModel_plots_Parallel;
    unsigned seed;
    struct DataHolder
    {
        unsigned fMCMCNIterations;

        double fMCMCprob;
        int fMCMCPhase;
        std::vector<double> parameters;

        DataHolder(TTree *tree) {
            TEST_CHECK(tree);

            tree->SetBranchAddress("Iteration",       &fMCMCNIterations);
            tree->SetBranchAddress("LogProbability",  &fMCMCprob);
            tree->SetBranchAddress("Phase",           &fMCMCPhase);

            // need #parameters to initialize the vector
            // todo dirty hack to remove fixed number that might change if cycle is dropped
            tree->GetEntry(0);
            parameters.resize(tree->GetNbranches() - 3);

            for (unsigned k = 0; k < parameters.size(); ++k)
            {
                std::string parName = "Parameter" + stringify(k);

                tree->SetBranchAddress(parName.c_str(),   &parameters[k]);
            }
        }
    };

    RunComparison(const RunComparison::Config & config) :
        config(config),
        GaussModel_plots_Serial(BAT_TESTDIR "parallel_TEST_GaussModel_plots_Serial.pdf"),
        GaussModel_plots_Parallel(BAT_TESTDIR "parallel_TEST_GaussModel_plots_Parallel.pdf"),
        seed(11)
    {
        CreateOutput(false);
        CreateOutput(true);
        Check();
    }

    void CreateOutput(const bool & parallelization) const {
        Output op;

        omp_set_dynamic(0);
        omp_set_num_threads(parallelization ? config.num_chains : 1);

        // open log file
        BCLog::OpenLog("log.txt");
        BCLog::SetLogLevel(BCLog::detail);

        // create new GaussModel object
        GaussModel m(parallelization ? "Parallel evaluation" : "Serial evaluation", config.num_parameters, config.lag);

        // set MCMC precision
        m.MCMCSetPrecision(BCEngineMCMC::kMedium);
        m.MCMCSetNIterationsRun(config.num_iterations);
        m.MCMCSetNChains(config.num_chains);

        // create new output object
        BCModelOutput mout(&m, parallelization ? config.rootFileNameParallel.c_str() : config.rootFileNameSerial.c_str());

        // switch writing of Markov Chains on
        mout.WriteMarkovChain(true);

        m.MCMCSetRandomSeed(seed);

        TStopwatch sw;
        sw.Start();

        // run MCMC and marginalize posterior wrt. all parameters
        // and all combinations of two parameters
        m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        m.MarginalizeAll();
        sw.Stop();
        double real_time=sw.RealTime();
        op.PrintValues(real_time);

        if (config.plot)
        {
            if (parallelization)
                m.PrintAllMarginalized(GaussModel_plots_Parallel.c_str());
            else
                m.PrintAllMarginalized(GaussModel_plots_Serial.c_str());
        }

        // close log file
        BCLog::CloseLog();

        // close output file
        mout.Close();
    }

    void Check() throw(TestCaseFailedException)
    {
        TFile * rfile1 = TFile::Open(config.rootFileNameParallel.c_str());
        TFile * rfile2 = TFile::Open(config.rootFileNameSerial.c_str());

        if (!rfile1)
          TEST_CHECK_FAILED(std::string("Could not open") + config.rootFileNameParallel);

        if (!rfile2)
          TEST_CHECK_FAILED(std::string("Could not open") + config.rootFileNameSerial);

        // find the beginning of the main run
        // assume there is at least one chain and all chains start main run at the same time
        long long main_begin = 0;
        {
            TTree * one = NULL;
            rfile1->GetObject("MarkovChainTree_0", one);

            if (!one)
                TEST_CHECK_FAILED("Could not locate first Markov chain");

            DataHolder oneData(one);
            for (; main_begin < one->GetEntries(); ++main_begin)
            {
                one->GetEntry(main_begin);
                if (oneData.fMCMCPhase == 2)
                    break;
            }
            delete one;
        }

        for (unsigned ichain = 0; ichain < config.num_chains; ++ichain) {
            TTree * one = NULL;
            TTree * two = NULL;

            rfile1->GetObject(TString::Format("MarkovChainTree_%d", ichain), one);
            rfile2->GetObject(TString::Format("MarkovChainTree_%d", ichain), two);

            if (!one)
              TEST_CHECK_FAILED(std::string("Could not locate tree in ") + config.rootFileNameParallel);
            if (!two)
              TEST_CHECK_FAILED(std::string("Could not locate tree in ") + config.rootFileNameSerial);

            DataHolder oneData(one);
            DataHolder twoData(two);
            long nEntries1 = one->GetEntries();
            long nEntries2 = two->GetEntries();

            TEST_CHECK_EQUAL(nEntries1, nEntries2);

            // check prerun(phase 1) and mainrun(phase 2) for identity
            for (unsigned phase = 1; phase < 3 ; ++phase)
            {
                // loop over iterations
                for(long long int entry = (phase - 1) * main_begin;
                        entry <  (phase - 1) * main_begin + config.num_entries ; ++entry) {
                    one->GetEntry(entry);
                    two->GetEntry(entry);

                    // loop over each dimension
                    for (unsigned k = 0; k < oneData.parameters.size(); k++)
                    {
                        TEST_CHECK_EQUAL(oneData.parameters[k], twoData.parameters[k]);
                    }
                    TEST_CHECK_EQUAL(oneData.fMCMCprob, twoData.fMCMCprob);
                    TEST_CHECK_EQUAL(oneData.fMCMCNIterations, twoData.fMCMCNIterations);
                    TEST_CHECK_EQUAL(oneData.fMCMCPhase, twoData.fMCMCPhase);
                }
            }
#if 0
            // check that two consecutive values change only on at most one parameter
            // take only first two parameters
            if (oneData.parameters.size() > 1)
            {
               one->GetEntry(0);
               two->GetEntry(0);

               double previousX = oneData.parameters[0];
               double previousY = oneData.parameters[1];

               for (unsigned phase = 1; phase < 3 ; ++phase)
               {
                  // loop over iterations
                  for(long long int entry = (phase - 1) * main_begin + 1;
                        entry <  (phase - 1) * main_begin + config.num_entries ; ++entry) {
                     one->GetEntry(entry);
                     two->GetEntry(entry);

                     std::cout << print_container(oneData.parameters) << std::endl;
                     if (oneData.parameters[0] != previousX and oneData.parameters[1] != previousY)
                     {
                        std::string msg = "Two (or more) parameters changed at once in iteration " + stringify(entry) + ":\n";
                        msg += stringify(previousX) + "->" + stringify(oneData.parameters[0]) + ", ";
                        msg += stringify(previousY) + "->" + stringify(oneData.parameters[1]);
                        TEST_CHECK_FAILED(msg);
                     }

                     previousX = oneData.parameters[0];
                     previousY = oneData.parameters[1];
                  }
               }
            }
            else
            {
               TEST_CHECK_FAILED("Need at least two parameters to check one-parameter-at-a-time proposal");
            }
#endif
            delete one;
            delete two;
        }

        // clean up
        rfile1->Close();
        rfile2->Close();
        delete rfile1;
        delete rfile2;
        remove(config.rootFileNameParallel.c_str());
        remove(config.rootFileNameSerial.c_str());
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

        config.num_chains = 4;
        config.num_parameters = 1;
        config.lag = 5e4;
        {
            RunComparison comparison(config);
        }
    }
} parallel_test;
}
#endif
