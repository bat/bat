/*
 * Copyright (C) 2007-2014, the BAT core developer team
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
#include "BCPriorConstant.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

// ---------------------------------------------------------
BCModel::BCModel(std::string name)
    : BCIntegrate(name)
    , fDataSet(0)
    , fPriorModel(0)
    , fDrawPriorFirst(true)
    , fFactorizedPrior(false)
{
    SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDefaultStyle);
}

// ---------------------------------------------------------
BCModel::BCModel(std::string filename, std::string name, bool loadObservables)
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
      fPriorModel(0),
      fBCH1DPriorDrawingOptions(other.fBCH1DPriorDrawingOptions),
      fBCH2DPriorDrawingOptions(other.fBCH2DPriorDrawingOptions),
      fBCH1DPosteriorDrawingOptions(other.fBCH1DPosteriorDrawingOptions),
      fBCH2DPosteriorDrawingOptions(other.fBCH2DPosteriorDrawingOptions),
      fDrawPriorFirst(other.fDrawPriorFirst),
      fFactorizedPrior(other.fFactorizedPrior)
{
}

// ---------------------------------------------------------
BCModel::~BCModel()
{
    delete fPriorModel;
}

// ---------------------------------------------------------
void swap(BCModel& A, BCModel& B)
{
    swap(static_cast<BCModel&>(A), static_cast<BCModel&>(B));
    std::swap(A.fDataSet, B.fDataSet);
    std::swap(A.fPriorModel, B.fPriorModel);
    std::swap(A.fBCH1DPriorDrawingOptions, B.fBCH1DPriorDrawingOptions);
    std::swap(A.fBCH2DPriorDrawingOptions, B.fBCH2DPriorDrawingOptions);
    std::swap(A.fBCH1DPosteriorDrawingOptions, B.fBCH1DPosteriorDrawingOptions);
    std::swap(A.fBCH2DPosteriorDrawingOptions, B.fBCH2DPosteriorDrawingOptions);
    std::swap(A.fDrawPriorFirst, B.fDrawPriorFirst);
    std::swap(A.fFactorizedPrior, B.fFactorizedPrior);
}

// ---------------------------------------------------------
double BCModel::LogProbabilityNN(const std::vector<double>& parameters)
{
    double ll = LogLikelihood(parameters);
    double lp = LogAPrioriProbability(parameters);
    if (MCMCGetCurrentChain() >= 0 and MCMCGetCurrentChain() < (int)fMCMCLogLikelihood_Provisional.size() and MCMCGetCurrentChain() < (int)fMCMCLogPrior_Provisional.size()) {
        fMCMCLogLikelihood_Provisional[MCMCGetCurrentChain()] = ll;
        fMCMCLogPrior_Provisional[MCMCGetCurrentChain()] = lp;
    }
    return ll + lp;
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
    fMCMCTree->Branch("LogLikelihood", &fMCMCTree_LogLikelihood, "log(likelihood)/D");
    fMCMCTree->Branch("LogPrior",      &fMCMCTree_LogPrior,      "log(prior)/D");
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
    if (!GetGlobalMode().empty())
        BCLog::OutSummary("   Best fit parameters (global):");
    PrintParameters(GetGlobalMode(), BCLog::OutSummary);

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

    for (int i = 0; i < int(parameters.size()); i++)
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

    if (index < GetNParameters() and GetParameter(index).Fixed())
        return prior;

    // check for factorized prior
    if (fFactorizedPrior and index < GetNParameters() and GetParameter(index).GetPrior() != NULL) {
        TH1* h = GetVariable(index).CreateH1("getprior1d_temp");
        prior = GetParameter(index).GetPrior()->GetBCH1D(h, Form("%s_%d_prior", GetSafeName().data(), index));
        delete h;

        if (prior.Valid()) {
            // correct for flat prior
            bool const_prior = (dynamic_cast<BCPriorConstant*>(GetParameter(index).GetPrior()) != NULL);
            if (const_prior) {
                prior.SetLocalMode((unsigned)0, GetParameter(index).GetRangeCenter());
                prior.SetNBands(0);
            }
        }
    }

    // else use marginalized prior, if it exists
    if (!prior.Valid() and fPriorModel->MarginalizedHistogramExists(index))
        prior = fPriorModel->GetMarginalized(index);

    if (prior.Valid())
        prior.GetHistogram()->SetTitle(Form("prior;%s;P(%s)", GetVariable(index).GetLatexNameWithUnits().data(), GetVariable(index).GetLatexName().data()));

    return prior;
}

// ---------------------------------------------------------
BCH2D BCModel::GetPrior(unsigned index1, unsigned index2)
{
    BCH2D prior;

    if (index1 > GetNVariables() or index2 > GetNVariables())
        return prior;

    if (index1 < GetNParameters() and GetParameter(index1).Fixed())
        return prior;

    if (index2 < GetNParameters() and GetParameter(index2).Fixed())
        return prior;

    std::string title = "prior";

    // check for factorized priors
    if (fFactorizedPrior and
            index1 < GetNParameters() and GetParameter(index1).GetPrior() != NULL and
            index2 < GetNParameters() and GetParameter(index2).GetPrior() != NULL) {

        TH2* h2 = fPriorModel->GetVariable(index1).CreateH2("getprior2d_temp", fPriorModel->GetVariable(index2));
        prior = GetParameter(index1).GetPrior()->GetBCH2D(GetParameter(index2).GetPrior(), h2, Form("h2d_prior_%s_%d_%d", GetName().data(), index1, index2));
        delete h2;

        if (prior.Valid()) {
            // correct for flat prior
            bool const_prior1 = (dynamic_cast<BCPriorConstant*>(GetParameter(index1).GetPrior()) != NULL);
            bool const_prior2 = (dynamic_cast<BCPriorConstant*>(GetParameter(index2).GetPrior()) != NULL);
            if (const_prior1)
                prior.SetLocalMode((unsigned)0, GetParameter(index1).GetRangeCenter());
            if (const_prior2)
                prior.SetLocalMode((unsigned)1, GetParameter(index2).GetRangeCenter());
            if (const_prior1 and !const_prior2)
                title += " (flat in " + fPriorModel->GetVariable(index1).GetLatexName() + ")";
            else if (!const_prior1 and const_prior2)
                title += " (flat in " + fPriorModel->GetVariable(index2).GetLatexName() + ")";
            else if (const_prior1 and const_prior2) {
                title += " (both flat)";
                prior.SetNBands(0);
            }
        }
    }

    // else use marginalized prior, if it exists
    if (!prior.Valid() and (fPriorModel->MarginalizedHistogramExists(index1, index2) or fPriorModel->MarginalizedHistogramExists(index2, index1)))
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
unsigned BCModel::PrintKnowledgeUpdatePlots(std::string filename, unsigned hdiv, unsigned vdiv, bool call_likelihood)
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCModel::PrintKnowledgeUpdatePlots : No variables defined!");
        return 0;
    }

    // prepare prior
    if ( !GetPriorModel(true, call_likelihood) or fPriorModel->GetNParameters() == 0 )
        return 0;

    BCAux::DefaultToPDF(filename);
    if (filename.empty())
        return 0;

    // call prior evalution once to set fFactorizedPrior:
    LogAPrioriProbability(std::vector<double>(GetNParameters(), 0));

    if (!fFactorizedPrior or !GetParameters().ArePriorsSet(true) or GetNObservables() > 0) {
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
        } else
            prior.CopyOptions(fBCH2DPriorDrawingOptions);
        posterior.CopyOptions(fBCH2DPosteriorDrawingOptions);
        h2.push_back(std::make_pair(prior, posterior));
    }

    if (h1.empty() and h2.empty()) {
        BCLog::OutWarning("BCModel::PrintKnowledgeUpdatePlots : No update plots to print.");
        return 0;
    }

    const unsigned nplots = h1.size() + h2.size();
    BCLog::OutSummary(Form("Printing knowledge update plots (%lu x 1D + %lu x 2D = %u) into file %s", h1.size(), h2.size(), nplots, filename.data()));
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
    c.Print((filename + "[").data());

    // Draw 1D Knowledge Update Plots
    for (unsigned i = 0; i < h1.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(filename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        BCAux::DrawKnowledgeUpdate(h1[i].first, h1[i].second, fDrawPriorFirst);

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h1.size() > 0) {
        c.Print(filename.c_str());
        c.Clear("D");
    }

    // Draw 2D Knowledge Update Plots
    for (unsigned i = 0; i < h2.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(filename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        // prior text?
        BCAux::DrawKnowledgeUpdate(h2[i].first, h2[i].second, fDrawPriorFirst);

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h2.size() > 0) {
        c.Print(filename.c_str());
        c.Clear("D");
    }
    c.Print(std::string( filename + "]").c_str());

    if (nplots > 100 && nplots % 100 != 0)
        BCLog::OutDetail(Form(" --> %d plots done", nplots));

    // return total number of drawn histograms
    return nplots;

}

/*
// ---------------------------------------------------------
unsigned BCModel::PrintPriors(std::string filename, unsigned hdiv, unsigned vdiv, bool call_likelihood)
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCModel::PrintPriors : No variables defined!");
        return 0;
    }

    // prepare prior
    if ( !GetPriorModel(true, call_likelihood) or fPriorModel->GetNParameters() == 0 )
        return 0;

    BCAux::DefaultToPDF(filename);
    if (filename.empty())
        return 0;

    // call prior evalution once to set fFactorizedPrior:
    LogAPrioriProbability(std::vector<double>(GetNParameters(), 0));

    if (!fFactorizedPrior or !GetParameters().ArePriorsSet(true) or GetNObservables() > 0) {
        BCLog::OutSummary("Marginalizing prior...");
         fPriorModel->MarginalizeAll();
    }
    fPriorModel->FindMode();

    // Find nonempty H1's
    std::vector<unsigned> H1Indices = GetH1DPrintOrder();
    std::vector<BCH1D> h1;
    h1.reserve(H1Indices.size());
    for (unsigned i = 0; i < H1Indices.size(); ++i) {
        h1.push_back(GetPrior(H1Indices[i]));
        if (h1.back().Valid())
            h1.back().CopyOptions(fBCH1DdrawingOptions);
        else
            h1.pop_back();
    }

    // Find nonempty H2's
    std::vector<std::pair<unsigned, unsigned> > H2Coords = GetH2DPrintOrder();
    std::vector<BCH2D> h2;
    h2.reserve(H2Coords.size());
    for (unsigned k = 0; k < H2Coords.size(); ++k) {
        unsigned i = H2Coords[k].first;
        unsigned j = H2Coords[k].second;
        h2.push_back(GetPrior(i, j));
        if (h2.back().Valid())
            h2.back().CopyOptions(fBCH2DdrawingOptions);
        else
            h2.pop_back();
    }

    return BCAux::PrintPlots(h1, h2, filename, hdiv, vdiv);
}
*/

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
