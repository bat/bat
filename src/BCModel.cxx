/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include "BCDataSet.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCParameter.h"
#include "BCPriorModel.h"
#include "BCPrior.h"
#include "BCConstantPrior.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

// ---------------------------------------------------------
BCModel::BCModel(const std::string& name)
    : BCIntegrate(name)
    , fDataSet(0)
    , fPriorModel(0)
    , fDrawPriorFirst(true)
    , fFactorizedPrior(false)
{
    SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDefaultStyle);
}

// ---------------------------------------------------------
BCModel::BCModel(const std::string& filename, const std::string& name, bool loadObservables)
    : BCIntegrate(filename, name, loadObservables)
    , fDataSet(0)
    , fPriorModel(0)
    , fDrawPriorFirst(true)
    , fFactorizedPrior(false)
{
    SetPriorConstantAll();
    SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDefaultStyle);
}

// ---------------------------------------------------------
BCModel::BCModel(const BCModel& other)
    : BCIntegrate(other),
      fDataSet(other.fDataSet),
      fPriorModel(NULL),
      fBCH1DPriorDrawingOptions(other.fBCH1DPriorDrawingOptions),
      fBCH2DPriorDrawingOptions(other.fBCH2DPriorDrawingOptions),
      fBCH1DPosteriorDrawingOptions(other.fBCH1DPosteriorDrawingOptions),
      fBCH2DPosteriorDrawingOptions(other.fBCH2DPosteriorDrawingOptions),
      fDrawPriorFirst(other.fDrawPriorFirst),
      fFactorizedPrior(other.fFactorizedPrior)
{
}

// ---------------------------------------------------------
BCModel& BCModel::operator=(const BCModel& other)
{
    BCIntegrate::operator=(other);
    fDataSet = other.fDataSet;
    fPriorModel = NULL;
    fBCH1DPriorDrawingOptions = other.fBCH1DPriorDrawingOptions;
    fBCH2DPriorDrawingOptions = other.fBCH2DPriorDrawingOptions;
    fBCH1DPosteriorDrawingOptions = other.fBCH1DPosteriorDrawingOptions;
    fBCH2DPosteriorDrawingOptions = other.fBCH2DPosteriorDrawingOptions;
    fDrawPriorFirst = other.fDrawPriorFirst;
    fFactorizedPrior = other.fFactorizedPrior;

    return *this;
}

// ---------------------------------------------------------
BCModel::~BCModel()
{
    delete fPriorModel;
}

// ---------------------------------------------------------
double BCModel::LogProbabilityNN(const std::vector<double>& parameters)
{
    ThreadLocalStorage& s = fMCMCThreadLocalStorage[GetCurrentChain()];

    // first calculate prior (which is usually cheaper than likelihood)
    s.log_prior = LogAPrioriProbability(parameters);
    // then calculate likelihood, or set to -inf, if prior already invalid
    if (std::isfinite(s.log_prior))
        s.log_likelihood = LogLikelihood(parameters);
    else
        s.log_likelihood = -std::numeric_limits<double>::infinity();

    s.log_probability = s.log_prior + s.log_likelihood;
    return s.log_probability;
}

// ---------------------------------------------------------
double BCModel::LogProbability(const std::vector<double>& parameters)
{
    // check if normalized
    if (GetIntegral() <= 0.) {
        BCLog::OutError("BCModel::LogProbability. Normalization not available or zero.");
        return 0.;
    }

    return LogProbabilityNN(parameters) - log(GetIntegral());
}

// ---------------------------------------------------------
void BCModel::InitializeMarkovChainTree(bool replacetree, bool replacefile)
{
    BCEngineMCMC::InitializeMarkovChainTree(replacetree, replacefile);
    if (!fMCMCTree)
        return;
    fMCMCTree->Branch("LogLikelihood", &fMCMCTree_State.log_likelihood, "log_likelihood/D");
    fMCMCTree->Branch("LogPrior",      &fMCMCTree_State.log_prior,      "log_prior/D");
}

// ---------------------------------------------------------
double BCModel::SamplingFunction(const std::vector<double>& /*parameters*/)
{
    double probability = 1;
    for (unsigned i = 0 ; i < GetNParameters() ; ++i)
        probability *= 1. / GetParameter(i).GetRangeWidth();
    return probability;
}

// ---------------------------------------------------------
double BCModel::HessianMatrixElement(unsigned index1, unsigned index2, const std::vector<double>& point)
{
    // check number of entries in vector
    if (point.size() < GetNParameters()) {
        BCLog::OutError("BCModel::HessianMatrixElement : Invalid number of entries in the vector.");
        return -1;
    }

    // define steps
    double nsteps = 1e5;
    double dx1 = GetVariable(index1).GetRangeWidth() / nsteps;
    double dx2 = GetVariable(index2).GetRangeWidth() / nsteps;

    // define points at which to evaluate
    std::vector<double> xpp = point;
    std::vector<double> xpm = point;
    std::vector<double> xmp = point;
    std::vector<double> xmm = point;

    xpp[index1] += dx1;
    xpp[index2] += dx2;

    xpm[index1] += dx1;
    xpm[index2] -= dx2;

    xmp[index1] -= dx1;
    xmp[index2] += dx2;

    xmm[index1] -= dx1;
    xmm[index2] -= dx2;

    // calculate probability at these points
    double ppp = Likelihood(xpp);
    double ppm = Likelihood(xpm);
    double pmp = Likelihood(xmp);
    double pmm = Likelihood(xmm);

    // return derivative
    return (ppp + pmm - ppm - pmp) / (4. * dx1 * dx2);
}

// ---------------------------------------------------------
void BCModel::PrintShortFitSummary()
{
    BCLog::OutSummary("---------------------------------------------------");
    BCLog::OutSummary(Form("Fit summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("   Number of parameters:  Npar  = %i", GetNParameters()));
    if (GetNDataPoints()) {
        BCLog::OutSummary(Form("   Number of data points: Ndata = %i", GetNDataPoints()));
        BCLog::OutSummary(Form("   Number of degrees of freedom = %i", GetNDoF()));
    }
    if (!GetBestFitParameters().empty())
        BCLog::OutSummary("   Best fit parameters (global):");
    PrintParameters(GetBestFitParameters(), BCLog::OutSummary);

    BCLog::OutSummary("---------------------------------------------------");
}

// ---------------------------------------------------------
void BCModel::PrintHessianMatrix(std::vector<double> parameters)
{
    // check number of entries in vector
    if (parameters.size() != GetNParameters()) {
        BCLog::OutError("BCModel::PrintHessianMatrix : Invalid number of entries in the vector");
        return;
    }

    // print to screen
    BCLog::OutSummary("Hessian matrix elements: ");
    BCLog::OutSummary("Parameter values:");

    for (unsigned i = 0; i < parameters.size(); ++i)
        BCLog::OutSummary(Form("Parameter %d : %f", i, parameters.at(i)));

    BCLog::OutSummary("Hessian matrix:");
    // loop over all parameter pairs
    for (unsigned int i = 0; i < GetNParameters(); i++)
        for (unsigned int j = 0; j < i; j++) {
            // calculate Hessian matrix element
            double hessianmatrixelement = HessianMatrixElement(i, j, parameters);

            // print to screen
            BCLog::OutSummary(Form("%d %d : %f", i, j, hessianmatrixelement));
        }
}

// ---------------------------------------------------------
BCPriorModel* BCModel::GetPriorModel(bool prepare, bool call_likelihood)
{
    if (!fPriorModel)
        fPriorModel = new BCPriorModel(*this, call_likelihood);
    else if (prepare)
        fPriorModel->PreparePriorModel();
    fPriorModel->SetCallLikelihood(call_likelihood);
    return fPriorModel;
}

// ---------------------------------------------------------
BCH1D BCModel::GetPrior(unsigned index)
{
    BCH1D prior;

    if (index > GetNVariables())
        return prior;

    if (index < GetNParameters() && GetParameter(index).Fixed())
        return prior;

    // check for factorized prior
    if (fFactorizedPrior && index < GetNParameters() && GetParameter(index).GetPrior() != NULL) {
        TH1* h = GetVariable(index).CreateH1("getprior1d_temp");
        prior = GetParameter(index).GetPrior()->GetBCH1D(h, Form("%s_%d_prior", GetSafeName().data(), index));
        delete h;

        if (prior.Valid()) {
            // correct for flat prior
            bool const_prior = (dynamic_cast<BCConstantPrior*>(GetParameter(index).GetPrior()) != NULL);
            if (const_prior) {
                prior.SetLocalMode((unsigned)0, GetParameter(index).GetRangeCenter());
                prior.SetNBands(0);
                prior.SetDrawLocalMode(0);
            }
        }
    }

    // else use marginalized prior, if it exists
    if (!prior.Valid() && fPriorModel->MarginalizedHistogramExists(index))
        prior = fPriorModel->GetMarginalized(index);

    if (prior.Valid())
        prior.GetHistogram()->SetTitle(Form("prior;%s;P(%s)", GetVariable(index).GetLatexNameWithUnits().data(), GetVariable(index).GetLatexName().data()));

    return prior;
}

// ---------------------------------------------------------
BCH2D BCModel::GetPrior(unsigned index1, unsigned index2)
{
    BCH2D prior;

    if (index1 > GetNVariables() || index2 > GetNVariables())
        return prior;

    if (index1 < GetNParameters() && GetParameter(index1).Fixed())
        return prior;

    if (index2 < GetNParameters() && GetParameter(index2).Fixed())
        return prior;

    std::string title = "prior";

    // check for factorized priors
    if (fFactorizedPrior &&
            index1 < GetNParameters() && GetParameter(index1).GetPrior() != NULL &&
            index2 < GetNParameters() && GetParameter(index2).GetPrior() != NULL) {

        TH2* h2 = fPriorModel->GetVariable(index1).CreateH2("getprior2d_temp", fPriorModel->GetVariable(index2));
        prior = GetParameter(index1).GetPrior()->GetBCH2D(GetParameter(index2).GetPrior(), h2, Form("h2d_prior_%s_%d_%d", GetName().data(), index1, index2));
        delete h2;

        if (prior.Valid()) {
            // correct for flat prior
            bool const_prior1 = (dynamic_cast<BCConstantPrior*>(GetParameter(index1).GetPrior()) != NULL);
            bool const_prior2 = (dynamic_cast<BCConstantPrior*>(GetParameter(index2).GetPrior()) != NULL);
            if (const_prior1)
                prior.SetLocalMode((unsigned)0, GetParameter(index1).GetRangeCenter());
            if (const_prior2)
                prior.SetLocalMode((unsigned)1, GetParameter(index2).GetRangeCenter());
            if (const_prior1 && !const_prior2)
                title += " (flat in " + fPriorModel->GetVariable(index1).GetLatexName() + ")";
            else if (!const_prior1 && const_prior2)
                title += " (flat in " + fPriorModel->GetVariable(index2).GetLatexName() + ")";
            else if (const_prior1 && const_prior2) {
                title += " (flat in both "
                         + fPriorModel->GetVariable(index1).GetLatexName()
                         + " and "
                         + fPriorModel->GetVariable(index2).GetLatexName()
                         + ")";
                prior.SetNBands(0);
                prior.SetDrawLocalMode(false);
            }
        }
    }

    // else use marginalized prior, if it exists
    if (!prior.Valid() && (fPriorModel->MarginalizedHistogramExists(index1, index2) || fPriorModel->MarginalizedHistogramExists(index2, index1)))
        prior = fPriorModel->GetMarginalized(index1, index2);

    if (prior.Valid())
        prior.GetHistogram()->SetTitle(Form("%s;%s;%s;P(%s, %s)",
                                            title.data(),
                                            GetVariable(index1).GetLatexNameWithUnits().data(),
                                            GetVariable(index2).GetLatexNameWithUnits().data(),
                                            GetVariable(index1).GetLatexName().data(), GetVariable(index2).GetLatexName().data()));

    return prior;
}

// ---------------------------------------------------------
unsigned BCModel::PrintKnowledgeUpdatePlots(const std::string& filename, unsigned hdiv, unsigned vdiv, bool call_likelihood)
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCModel::PrintKnowledgeUpdatePlots : No variables defined!");
        return 0;
    }

    // prepare prior
    if (!GetPriorModel(true, call_likelihood) || fPriorModel->GetNParameters() == 0 )
        return 0;

    std::string newFilename(filename);
    BCAux::DefaultToPDF(newFilename);
    if (newFilename.empty())
        return 0;

    // call prior evalution once to set fFactorizedPrior:
    LogAPrioriProbability(std::vector<double>(GetNParameters(), 0));

    if (!fFactorizedPrior || !GetParameters().ArePriorsSet(true) || GetNObservables() > 0) {
        BCLog::OutSummary("Marginalizing prior...");
        fPriorModel->MarginalizeAll();
    }
    fPriorModel->FindMode();

    // Find nonempty 1D prior--posterior pairs
    std::vector<unsigned> H1Indices = GetH1DPrintOrder();
    std::vector<std::pair<BCH1D, BCH1D> > h1;
    h1.reserve(H1Indices.size());
    for (unsigned i = 0; i < H1Indices.size(); ++i) {
        BCH1D prior = GetPrior(H1Indices[i]);
        if (!prior.Valid())
            continue;
        BCH1D posterior = GetMarginalized(H1Indices[i]);
        if (!posterior.Valid())
            continue;
        posterior.GetHistogram()->SetTitle("posterior");
        if (prior.GetNBands() == 0) {
            prior.CopyOptions(fBCH1DPriorDrawingOptions);
            prior.SetNBands(0);
            prior.SetDrawLocalMode(false);
        } else
            prior.CopyOptions(fBCH1DPriorDrawingOptions);
        posterior.CopyOptions(fBCH1DPosteriorDrawingOptions);
        h1.push_back(std::make_pair(prior, posterior));
    }

    // Find nonempty 2D prior--posterior pairs
    std::vector<std::pair<unsigned, unsigned> > H2Coords = GetH2DPrintOrder();
    std::vector<std::pair<BCH2D, BCH2D> > h2;
    h2.reserve(H2Coords.size());
    for (unsigned i = 0; i < H2Coords.size(); ++i) {
        BCH2D prior = GetPrior(H2Coords[i].first, H2Coords[i].second);
        if (!prior.Valid())
            continue;
        BCH2D posterior = GetMarginalized(H2Coords[i].first, H2Coords[i].second);
        if (!posterior.Valid())
            continue;
        posterior.GetHistogram()->SetTitle("posterior");
        if (prior.GetNBands() == 0) {
            prior.CopyOptions(fBCH2DPriorDrawingOptions);
            prior.SetNBands(0);
            prior.SetDrawLocalMode(false);
        } else
            prior.CopyOptions(fBCH2DPriorDrawingOptions);
        posterior.CopyOptions(fBCH2DPosteriorDrawingOptions);
        h2.push_back(std::make_pair(prior, posterior));
    }

    if (h1.empty() && h2.empty()) {
        BCLog::OutWarning("BCModel::PrintKnowledgeUpdatePlots : No update plots to print.");
        return 0;
    }

    const unsigned nplots = h1.size() + h2.size();
    BCLog::OutSummary(Form("Printing knowledge update plots (%lu x 1D + %lu x 2D = %u) into file %s", h1.size(), h2.size(), nplots, newFilename.data()));
    if (nplots > 100)
        BCLog::OutDetail("This can take a while...");

    if (hdiv < 1) hdiv = 1;
    if (vdiv < 1) vdiv = 1;

    int c_width  = 297 * 4;
    int c_height = 210 * 4;
    if (hdiv < vdiv)
        std::swap(c_width, c_height);

    // create canvas and prepare postscript
    TCanvas c("c", "canvas", c_width, c_height);
    c.Divide(hdiv, vdiv);

    unsigned n = 0;

    // open file
    c.Print((newFilename + "[").data());

    // Draw 1D Knowledge Update Plots
    for (unsigned i = 0; i < h1.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(newFilename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        BCAux::DrawKnowledgeUpdate(h1[i].first, h1[i].second, fDrawPriorFirst, fObjectTrash);

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h1.size() > 0) {
        c.Print(newFilename.c_str());
        c.Clear("D");
    }

    // Draw 2D Knowledge Update Plots
    for (unsigned i = 0; i < h2.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(newFilename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        // prior text?
        BCAux::DrawKnowledgeUpdate(h2[i].first, h2[i].second, fDrawPriorFirst, fObjectTrash);

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h2.size() > 0) {
        c.Print(newFilename.c_str());
        c.Clear("D");
    }
    c.Print(std::string(newFilename + "]").c_str());

    if (nplots > 100 && nplots % 100 != 0)
        BCLog::OutDetail(Form(" --> %d plots done", nplots));

    // return total number of drawn histograms
    return nplots;

}

// ---------------------------------------------------------
void BCModel::SetKnowledgeUpdateDrawingStyle(BCAux::BCKnowledgeUpdateDrawingStyle style)
{
    BCAux::SetKnowledgeUpdateDrawingStyle(fBCH1DPriorDrawingOptions, fBCH1DPosteriorDrawingOptions, style);
    BCAux::SetKnowledgeUpdateDrawingStyle(fBCH2DPriorDrawingOptions, fBCH2DPosteriorDrawingOptions, style);

    switch (style) {

        case BCAux::kKnowledgeUpdateDetailedPosterior:
            SetDrawPriorFirst(false);
            break;

        case BCAux::kKnowledgeUpdateDetailedPrior:
            SetDrawPriorFirst(true);
            break;

        case BCAux::kKnowledgeUpdateDefaultStyle:
        default:
            break;
    }
}
