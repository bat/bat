/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include "BCDataPoint.h"
#include "BCDataSet.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCParameter.h"
#include "BCPriorModel.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TList.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMarker.h>
#include <TArrow.h>

#include <fstream>
#include <iomanip>
#include <set>

// ---------------------------------------------------------
BCModel::BCModel(const char * name)
	: BCIntegrate(name)
  , fModelAPriori(0)
  , fModelAPosteriori(0)
  , fDataSet(0)
	, fPriorModel(0)
	, fBCH1DPriorDrawingOptions(new BCH1D)
	, fBCH2DPriorDrawingOptions(new BCH2D)
	, fBCH1DPosteriorDrawingOptions(new BCH1D)
	, fBCH2DPosteriorDrawingOptions(new BCH2D)
	, fPriorPosteriorNormalOrder(true)
{
	SetKnowledgeUpdateDrawingStyle(kKnowledgeUpdateDefaultStyle);
}

// ---------------------------------------------------------
BCModel::BCModel(std::string filename, std::string name, bool reuseObservables)
	: BCIntegrate(name.data())
  , fModelAPriori(0)
  , fModelAPosteriori(0)
  , fDataSet(0)
	, fPriorModel(0)
	, fBCH1DPriorDrawingOptions(new BCH1D)
	, fBCH2DPriorDrawingOptions(new BCH2D)
	, fBCH1DPosteriorDrawingOptions(new BCH1D)
	, fBCH2DPosteriorDrawingOptions(new BCH2D)
	, fPriorPosteriorNormalOrder(true)
{
	LoadMCMC(filename,"","",reuseObservables);
	SetPriorConstantAll();
	SetKnowledgeUpdateDrawingStyle(kKnowledgeUpdateDefaultStyle);
}

// ---------------------------------------------------------
BCModel::BCModel(const BCModel & bcmodel)
	: BCIntegrate(bcmodel)
	, fPriorModel(0)
	, fBCH1DPriorDrawingOptions(new BCH1D)
	, fBCH2DPriorDrawingOptions(new BCH2D)
	, fBCH1DPosteriorDrawingOptions(new BCH1D)
	, fBCH2DPosteriorDrawingOptions(new BCH2D)
{
	Copy(bcmodel);
}

// ---------------------------------------------------------
BCModel::~BCModel() {
	delete fPriorModel;
	delete fBCH1DPriorDrawingOptions;
	delete fBCH2DPriorDrawingOptions;
	delete fBCH1DPosteriorDrawingOptions;
	delete fBCH2DPosteriorDrawingOptions;
}

// ---------------------------------------------------------
void BCModel::Copy(const BCModel & bcmodel)
{
   //  called for the second time in copy constructor? do copy-and-swap instead
   //   BCIntegrate::Copy(bcmodel);
   fName                            = bcmodel.fName;
   fModelAPriori                    = bcmodel.fModelAPriori;
   fModelAPosteriori                = bcmodel.fModelAPosteriori;
	 fDataSet                         = bcmodel.fDataSet;

	 fBCH1DPriorDrawingOptions -> CopyOptions(*(bcmodel.fBCH1DPriorDrawingOptions));
	 fBCH2DPriorDrawingOptions -> CopyOptions(*(bcmodel.fBCH2DPriorDrawingOptions));
	 fBCH1DPosteriorDrawingOptions -> CopyOptions(*(bcmodel.fBCH1DPosteriorDrawingOptions));
	 fBCH2DPosteriorDrawingOptions -> CopyOptions(*(bcmodel.fBCH2DPosteriorDrawingOptions));
	 fPriorPosteriorNormalOrder = bcmodel.fPriorPosteriorNormalOrder;
}

// ---------------------------------------------------------
double BCModel::LogProbabilityNN(const std::vector<double> &parameters) {
	double ll = LogLikelihood(parameters);
	double lp = LogAPrioriProbability(parameters);
	if (MCMCGetCurrentChain()>=0 and MCMCGetCurrentChain()<(int)fMCMCLogLikelihood_Provisional.size() and MCMCGetCurrentChain()<(int)fMCMCLogPrior_Provisional.size()) {
		fMCMCLogLikelihood_Provisional[MCMCGetCurrentChain()] = ll;
		fMCMCLogPrior_Provisional[MCMCGetCurrentChain()] = lp;
	}
	return ll + lp;
}

// ---------------------------------------------------------
double BCModel::LogProbability(const std::vector<double> &parameters)
{
   // check if normalized
  if (GetIntegral() <= 0.) {
      BCLog::OutError("BCModel::LogProbability. Normalization not available or zero.");
      return 0.;
   }

  return LogProbabilityNN(parameters) - log(GetIntegral());
}

// ---------------------------------------------------------
void BCModel::InitializeMarkovChainTree(bool replacetree, bool replacefile) {
	BCEngineMCMC::InitializeMarkovChainTree(replacetree,replacefile);
	if (!fMCMCTree)
		return;
	fMCMCTree -> Branch("LogLikelihood", &fMCMCTree_LogLikelihood, "log(likelihood)/D");
	fMCMCTree -> Branch("LogPrior",      &fMCMCTree_LogPrior,      "log(prior)/D");
}

// ---------------------------------------------------------
double BCModel::SamplingFunction(const std::vector<double> & /*parameters*/)
{
   double probability = 1;
   for (unsigned i = 0 ; i < GetNParameters() ; ++i)
		 probability *= 1. / GetParameter(i)->GetRangeWidth();
   return probability;
}

// ---------------------------------------------------------
double BCModel::HessianMatrixElement(const BCParameter * par1, const BCParameter * par2, std::vector<double> point)
{
   // check number of entries in vector
   if (point.size() != GetNParameters()) {
      BCLog::OutError("BCModel::HessianMatrixElement : Invalid number of entries in the vector.");
      return -1;
   }

   // define steps
   double nsteps = 1e5;
   double dx1 = par1->GetRangeWidth() / nsteps;
   double dx2 = par2->GetRangeWidth() / nsteps;

   // define points at which to evaluate
   std::vector<double> xpp = point;
   std::vector<double> xpm = point;
   std::vector<double> xmp = point;
   std::vector<double> xmm = point;

   unsigned idx1 = fParameters.Index(par1->GetName());
   unsigned idx2 = fParameters.Index(par2->GetName());

   xpp[idx1] += dx1;
   xpp[idx2] += dx2;

   xpm[idx1] += dx1;
   xpm[idx2] -= dx2;

   xmp[idx1] -= dx1;
   xmp[idx2] += dx2;

   xmm[idx1] -= dx1;
   xmm[idx2] -= dx2;

   // calculate probability at these points
   double ppp = Likelihood(xpp);
   double ppm = Likelihood(xpm);
   double pmp = Likelihood(xmp);
   double pmm = Likelihood(xmm);

   // return derivative
   return (ppp + pmm - ppm - pmp) / (4. * dx1 * dx2);
}

// ---------------------------------------------------------
void BCModel::PrintShortFitSummary(int chi2flag) {
	BCLog::OutSummary("---------------------------------------------------");
	BCLog::OutSummary(Form("Fit summary for model \'%s\':", GetName().data()));
	BCLog::OutSummary(Form("   Number of parameters:  Npar  = %i", GetNParameters()));
	if (GetNDataPoints()) {
		BCLog::OutSummary(Form("   Number of data points: Ndata = %i", GetNDataPoints()));
		BCLog::OutSummary(Form("   Number of degrees of freedom = %i", GetNDoF()));
	}
	if (!GetBestFitParameters().empty())
		BCLog::OutSummary("   Best fit parameters (global):");
	PrintParameters(GetBestFitParameters(),BCLog::OutSummary);
	
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
				double hessianmatrixelement = HessianMatrixElement(GetParameter(i),GetParameter(j),parameters);

				// print to screen
				BCLog::OutSummary(Form("%d %d : %f", i, j, hessianmatrixelement));
      }
}

// ---------------------------------------------------------
BCPriorModel * BCModel::GetPriorModel(bool prepare, bool call_likelihood) {
	if (!fPriorModel)
		fPriorModel = new BCPriorModel(this, call_likelihood);
	else if (prepare)
		fPriorModel -> PreparePriorModel();
	fPriorModel -> SetCallLikelihood(call_likelihood);
	return fPriorModel;
}

// ---------------------------------------------------------
int BCModel::DrawKnowledgeUpdatePlot1D(unsigned index, bool flag_slice_post, bool flag_slice_prior) {
	// Get Prior
	BCH1D * bch1d_prior = 0;
	TLine * const_prior = (index<GetNParameters() and GetParameter(index)->GetPriorType()==BCParameter::kPriorConstant) ? new TLine() : 0;
	TF1   * f1_prior    = (const_prior or index>=GetNParameters()) ? 0 : GetParameter(index)->GetPriorTF1();
	TH1   * h1_prior    = (const_prior or index>=GetNParameters()) ? 0 : GetParameter(index)->GetPriorTH1();
	
	if (const_prior) {
		TH1D * h = new TH1D(TString::Format("%s_prior_%d_const",GetSafeName().data(),index),"",1,GetVariable(index)->GetLowerLimit(),GetVariable(index)->GetUpperLimit());
		h -> SetBinContent(1,1);
		bch1d_prior = new BCH1D(h);

	} else if (f1_prior) {
		TH1D * h = new TH1D(TString::Format("%s_prior_%d_f1",GetSafeName().data(),index),"",f1_prior->GetNpx(),GetVariable(index)->GetLowerLimit(),GetVariable(index)->GetUpperLimit());
		h -> Add(f1_prior,1,"I");
		bch1d_prior = new BCH1D(h);

	} else if (h1_prior)
		bch1d_prior = new BCH1D(h1_prior);

	else if (flag_slice_prior and index<fPriorModel->GetNParameters() and fPriorModel->GetNParameters()<=2) {
		double log_max_val;
		if (fPriorModel->GetNParameters()==2) {
			if (index==0)
				bch1d_prior = new BCH1D(fPriorModel->GetSlice(0,1,log_max_val)->ProjectionX(Form("prior_projx_%s_%i",fPriorModel->GetSafeName().data(),index)));
			else
				bch1d_prior = new BCH1D(fPriorModel->GetSlice(0,1,log_max_val)->ProjectionY(Form("prior_projy_%s_%i",fPriorModel->GetSafeName().data(),index)));
		} else
			bch1d_prior = new BCH1D(fPriorModel->GetSlice(index,log_max_val));

	} else if (fPriorModel->MarginalizedHistogramExists(index))
		bch1d_prior = fPriorModel->GetMarginalized(index);
	
	// if prior doesn't exist, exit
	if (!bch1d_prior)
		return 0;

	// Get Posterior
	BCH1D* bch1d_posterior = 0;
	if (flag_slice_post and index<GetNParameters() and GetNParameters()<=2) {
		double log_max_val;
		if (GetNParameters()==2) {
			if (index==0)
				bch1d_posterior = new BCH1D(GetSlice(0,1,log_max_val)->ProjectionX(Form("posterior_projx_%s_%i",GetSafeName().data(),index)));
			else
				bch1d_posterior = new BCH1D(GetSlice(0,1,log_max_val)->ProjectionY(Form("posterior_projy_%s_%i",GetSafeName().data(),index)));
		} else
			bch1d_posterior = new BCH1D(GetSlice(index,log_max_val));

	} else if (MarginalizedHistogramExists(index))
		bch1d_posterior = GetMarginalized(index);
	
	// if marginal doesn't exist, exit
	if (!bch1d_posterior)
		return 0;

	bch1d_prior -> CopyOptions(*fBCH1DPriorDrawingOptions);
	bch1d_prior -> SetDrawLegend(false);
	bch1d_posterior -> CopyOptions(*fBCH1DPosteriorDrawingOptions);
	bch1d_posterior -> SetDrawLegend(false);

	gPad -> SetLogx(fBCH1DPriorDrawingOptions->GetLogx() or fBCH1DPosteriorDrawingOptions->GetLogx());
	gPad -> SetLogy(fBCH1DPriorDrawingOptions->GetLogy() or fBCH1DPosteriorDrawingOptions->GetLogy());

	// get maximum
	double maxy = 1.1 * std::max<double>(bch1d_prior->GetHistogram()->GetMaximum(), bch1d_posterior->GetHistogram()->GetMaximum());
	double miny = 0.0;
	if (gPad->GetLogy()) {
		maxy *= 2;
		miny = 0.5 * std::min<double>(bch1d_prior->GetHistogram()->GetMinimum(0),bch1d_posterior->GetHistogram()->GetMinimum(0));
	}

	// draw axes
	TH2D * h2_axes = new TH2D(TString::Format("h2_axes_%s_knowledge_update_%d",GetSafeName().data(),index), TString::Format(";%s;P(%s|Data)",GetVariable(index)->GetLatexName().data(),GetVariable(index)->GetLatexName().data()),
														10, GetVariable(index)->GetLowerLimit(), GetVariable(index)->GetUpperLimit(),
														10, miny, maxy);
	h2_axes -> SetStats(false);
	h2_axes -> GetXaxis() -> SetNdivisions(508);
	h2_axes -> Draw();

	// Draw histograms
	if (!fPriorPosteriorNormalOrder) // posterior first
		bch1d_posterior -> Draw();
	bch1d_prior -> Draw();
	if (fPriorPosteriorNormalOrder) // prior first
		bch1d_posterior -> Draw();

	// create / draw legend(s)
	gPad->SetTopMargin(0.02);

	if ( bch1d_prior->GetLegend()->GetNRows() > 0  and  bch1d_posterior->GetLegend()->GetNRows() > 0 ) {
		// both legends have entries, draw both
		bch1d_prior->GetLegend()->SetHeader("prior");
		bch1d_posterior->GetLegend()->SetHeader("posterior");

		// Draw prior legend on top left
		double y1ndc_prior = bch1d_prior -> ResizeLegend();
		bch1d_prior->GetLegend() -> SetX2NDC(bch1d_prior->GetLegend()->GetX1NDC() + 45e-2*(bch1d_prior->GetLegend()->GetX2NDC()-bch1d_prior->GetLegend()->GetX1NDC()));
		bch1d_prior->GetLegend() -> Draw();

		// Draw posterior legend on top right
		double y1ndc_posterior = bch1d_posterior -> ResizeLegend();
		bch1d_posterior->GetLegend() -> SetX1NDC(bch1d_posterior->GetLegend()->GetX1NDC() + 55e-2*(bch1d_posterior->GetLegend()->GetX2NDC()-bch1d_posterior->GetLegend()->GetX1NDC()));
		bch1d_posterior->GetLegend() -> Draw();

		gPad -> SetTopMargin(1-std::min<double>(y1ndc_prior,y1ndc_posterior)+0.01);

	} else {
		// only one legend to draw

		TLegend * legend = 0;

		if ( bch1d_posterior->GetLegend()->GetNRows() > 0 ) {
			// posterior legend alone has entries

			legend = bch1d_posterior->GetLegend();
			for (int i=0; legend->GetListOfPrimitives()->GetEntries(); ++i) {
				TLegendEntry * le = (TLegendEntry*)(legend->GetListOfPrimitives()->At(i));
				if (!le) break;
				if (strlen(le->GetLabel())==0) continue;
				le-> SetLabel(TString::Format("%s of posterior",le->GetLabel()).Data());
			}
			legend -> AddEntry(bch1d_prior->GetHistogram(), "prior", "L");

			bch1d_posterior -> ResizeLegend();

		} else if ( bch1d_prior->GetLegend()->GetNRows() > 0 ) {
			// prior legend alone has entries

			legend = bch1d_prior->GetLegend();
			for (int i=0; legend->GetListOfPrimitives()->GetEntries(); ++i) {
				TLegendEntry * le = (TLegendEntry*)(legend->GetListOfPrimitives()->At(i));
				if (!le) break;
				if (strlen(le->GetLabel())==0) continue;
				le-> SetLabel(TString::Format("%s of prior",le->GetLabel()).Data());
			}
			legend -> AddEntry(bch1d_posterior->GetHistogram(), "posterior", "L");
			bch1d_prior -> ResizeLegend();

		}	else {
			// neither legend has entries

			legend = bch1d_prior->GetLegend();
			legend -> SetNColumns(2);
			legend -> AddEntry(bch1d_prior->GetHistogram(),     "prior", "L");
			legend -> AddEntry(bch1d_posterior->GetHistogram(), "posterior", "L");

			bch1d_prior -> ResizeLegend();
		}

		// Draw legend on top of histogram
		legend -> Draw();

		// rescale top margin
		gPad -> SetTopMargin(1-legend->GetY1NDC()+0.01);		
	}

	gPad -> RedrawAxis();
	gPad -> Update();
	return 1;
}

// ---------------------------------------------------------
int BCModel::DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice) {
	
	if (index1 == index2)
		return 0;
	if (index1 > index2)
		return DrawKnowledgeUpdatePlot2D(index2,index1,flag_slice);

	// Get Posterior
	TH2 * h2d_posterior = 0;
	if (flag_slice and GetNParameters()==2 and index1<GetNParameters() and index2<GetNParameters()) {
		double log_max_val;
		h2d_posterior = GetSlice(index1,index2,log_max_val);
	}	else if (MarginalizedHistogramExists(index1,index2))
		h2d_posterior = GetMarginalizedHistogram(index1,index2);

	if (!h2d_posterior)
		return 0;

	// Get Prior
	bool const_prior1 = index1<GetNParameters() and GetParameter(index1)->GetPriorType()==BCParameter::kPriorConstant;
	TF1 * f1_prior1   = (const_prior1 or index1>=GetNParameters()) ? 0 : GetParameter(index1)->GetPriorTF1();
	TH1 * h1_prior1   = (const_prior1 or index1>=GetNParameters()) ? 0 : GetParameter(index1)->GetPriorTH1();
	bool auto_prior1 = const_prior1 or f1_prior1 or h1_prior1;

	bool const_prior2 = index2<GetNParameters() and GetParameter(index2)->GetPriorType()==BCParameter::kPriorConstant;
	TF1 * f1_prior2   = (const_prior2 or index2>=GetNParameters()) ? 0 : GetParameter(index2)->GetPriorTF1();
	TH1 * h1_prior2   = (const_prior2 or index2>=GetNParameters()) ? 0 : GetParameter(index2)->GetPriorTH1();
	bool auto_prior2 = const_prior2 or f1_prior2 or h1_prior2;

	TH2 * h2d_prior = 0;

	if (!auto_prior1 or !auto_prior2) { // one or both prior pre-defined
		if (flag_slice and GetNParameters()==2 and index1<GetNParameters() and index2<GetNParameters()) {
			double log_max_val;
			h2d_prior = fPriorModel -> GetSlice(index1,index2,log_max_val);
		}	else if (fPriorModel->MarginalizedHistogramExists(index1, index2))
			h2d_prior = fPriorModel -> GetMarginalizedHistogram(index1,index2);
	}

	// if not predefined, use the projection of the marginalization
	if (!auto_prior1 and h2d_prior)
		h1_prior1 = h2d_prior -> ProjectionX(TString::Format("h1_prior1_%s_%d",GetName().data(),index1));
	if (!auto_prior2 and h2d_prior)
		h1_prior2 = h2d_prior -> ProjectionY(TString::Format("h1_prior2_%s_%d",GetName().data(),index2));

	if (!h2d_prior)
		h2d_prior = fPriorModel->GetVariable(index1) -> CreateH2(TString::Format("h2d_prior_%s_%d_%d",GetName().data(),index1,index2).Data(),fPriorModel->GetVariable(index2));

	// Set 2D-prior histogram binning to match 1D binning, for if prior was defined by histogram
	if (h2d_prior) {
		TAxis * xaxis = (h1_prior1) ? h1_prior1->GetXaxis() : h2d_prior->GetXaxis();
		TAxis * yaxis = (h1_prior2) ? h1_prior2->GetXaxis() : h2d_prior->GetYaxis();

		int n_xbins = xaxis->GetNbins();
		double xbins[n_xbins+1];
		xaxis -> GetLowEdge(xbins);
		xbins[n_xbins] = xaxis->GetXmax();

		int n_ybins = yaxis->GetNbins();
		double ybins[n_ybins+1];
		yaxis -> GetLowEdge(ybins);
		ybins[n_ybins] = yaxis->GetXmax();

		h2d_prior -> SetBins(n_xbins,xbins,n_ybins,ybins);
	}
	
	for (int i = 1; i <= h2d_prior->GetNbinsX(); ++i) {
		// x prior
		double x = 1;
		if (f1_prior1)
			x = f1_prior1 -> Eval(h2d_prior->GetXaxis()->GetBinCenter(i));
		else if (h1_prior1)
			x = h1_prior1 -> GetBinContent(h1_prior1->FindFixBin(h2d_prior->GetXaxis()->GetBinCenter(i)));
		
		for (int j = 1; j <= h2d_prior->GetNbinsY(); ++j) {
			// y prior
			double y = 1;
			if (f1_prior2)
				y = f1_prior2 -> Eval(h2d_prior->GetYaxis()->GetBinCenter(j));
			else if (h1_prior2)
				y = h1_prior2 -> GetBinContent(h1_prior2->FindFixBin(h2d_prior->GetYaxis()->GetBinCenter(j)));
			
			h2d_prior -> SetBinContent(i,j,x*y);
		}
	}

	if (!h2d_prior)
		return 0;


	// Create BCH2D's (these normalize the TH2's)
	BCH2D * bch2d_prior     = new BCH2D(h2d_prior);
	BCH2D * bch2d_posterior = new BCH2D(h2d_posterior);

	bch2d_prior     -> CopyOptions(*fBCH2DPriorDrawingOptions);
	bch2d_prior -> SetDrawLegend(false);
	bch2d_posterior -> CopyOptions(*fBCH2DPosteriorDrawingOptions);
	bch2d_posterior -> SetDrawLegend(false);

	if (const_prior1)
		bch2d_prior -> SetLocalMode(0,fPriorModel->GetVariable(index1)->GetRangeCenter());
	if (const_prior2)
		bch2d_prior -> SetLocalMode(1,fPriorModel->GetVariable(index2)->GetRangeCenter());
	std::string prior_text = "";
	if (const_prior1 and !const_prior2)
		prior_text = Form(" (flat in %s)",fPriorModel->GetVariable(index1)->GetLatexName().data());
	else if (!const_prior1 and const_prior2)
		prior_text = Form(" (flat in %s)",fPriorModel->GetVariable(index2)->GetLatexName().data());
	else if (const_prior1 and const_prior2) {
		prior_text = " (both flat)";
		bch2d_prior -> SetNBands(0);
	}

	gPad -> SetLogx(fBCH2DPriorDrawingOptions->GetLogx() or fBCH2DPosteriorDrawingOptions->GetLogx());
	gPad -> SetLogy(fBCH2DPriorDrawingOptions->GetLogy() or fBCH2DPosteriorDrawingOptions->GetLogy());
	gPad -> SetLogz(fBCH2DPriorDrawingOptions->GetLogz() or fBCH2DPosteriorDrawingOptions->GetLogz());

	// draw axes
	TH2D * h2_axes = new TH2D(TString::Format("h2_axes_%s_knowledge_update_%d_%d",GetSafeName().data(),index1,index2), TString::Format(";%s;%s;P(%s %s|Data)",GetVariable(index1)->GetLatexName().data(),GetVariable(index2)->GetLatexName().data(),GetVariable(index1)->GetLatexName().data(),GetVariable(index2)->GetLatexName().data()),
														10, GetVariable(index1)->GetLowerLimit(), GetVariable(index1)->GetUpperLimit(),
														10, GetVariable(index2)->GetLowerLimit(), GetVariable(index2)->GetUpperLimit());
	h2_axes -> SetStats(false);
	h2_axes -> GetXaxis() -> SetNdivisions(508);
	h2_axes -> Draw();

	// ROOT options for both prior and posterior should contain "same" (as they do by default)

	if (!fPriorPosteriorNormalOrder) // posterior first
		bch2d_posterior -> Draw();
	bch2d_prior -> Draw();
	if (fPriorPosteriorNormalOrder) // posterior second
		bch2d_posterior -> Draw();

	// create / draw legend(s)
	if ( bch2d_prior->GetLegend()->GetNRows() > 0  and  bch2d_posterior->GetLegend()->GetNRows() > 0 ) {
		// both legends have entries, draw both

		bch2d_prior->GetLegend()->SetHeader((std::string("prior")+prior_text).data());
		bch2d_posterior->GetLegend()->SetHeader("posterior");

		// Draw prior legend on top left
		double y1ndc_prior = bch2d_prior -> ResizeLegend();
		bch2d_prior->GetLegend() -> SetX2NDC(bch2d_prior->GetLegend()->GetX1NDC() + 45e-2*(bch2d_prior->GetLegend()->GetX2NDC()-bch2d_prior->GetLegend()->GetX1NDC()));
		bch2d_prior->GetLegend() -> Draw();

		// Draw posterior legend on top right
		double y1ndc_posterior = bch2d_posterior -> ResizeLegend();
		bch2d_posterior->GetLegend() -> SetX1NDC(bch2d_posterior->GetLegend()->GetX1NDC() + 55e-2*(bch2d_posterior->GetLegend()->GetX2NDC()-bch2d_posterior->GetLegend()->GetX1NDC()));
		bch2d_posterior->GetLegend() -> Draw();

		gPad -> SetTopMargin(1-std::min<double>(y1ndc_prior,y1ndc_posterior)+0.01);

	} else {
		// only one legend to draw

		TLegend * legend = 0;
		if ( bch2d_posterior->GetLegend()->GetNRows() > 0 ) {
			// posterior legend alone has entries

			legend = bch2d_posterior->GetLegend();
			for (int i=0; legend->GetListOfPrimitives()->GetEntries(); ++i) {
				TLegendEntry * le = (TLegendEntry*)(legend->GetListOfPrimitives()->At(i));
				if (!le) break;
				if (strlen(le->GetLabel())==0) continue;
				le-> SetLabel(TString::Format("%s of posterior",le->GetLabel()).Data());
			}
			legend -> AddEntry(bch2d_prior->GetHistogram(), (std::string("prior")+prior_text).data(), "L");

			bch2d_posterior -> ResizeLegend();

		} else if ( bch2d_prior->GetLegend()->GetNRows() > 0 ) {
			// prior legend alone has entries

			legend = bch2d_prior->GetLegend();
			for (int i=0; legend->GetListOfPrimitives()->GetEntries(); ++i) {
				TLegendEntry * le = (TLegendEntry*)(legend->GetListOfPrimitives()->At(i));
				if (!le) break;
				if (strlen(le->GetLabel())==0) continue;
				le-> SetLabel(TString::Format("%s of prior",le->GetLabel()).Data());
			}
			if (!prior_text.empty())
				legend -> AddEntry((TObject*)0,(std::string("prior: ")+prior_text).data(),"");
			legend -> AddEntry(bch2d_posterior->GetHistogram(), "posterior", "L");

			bch2d_prior -> ResizeLegend();

		}	else {
			// neither legend has entries

			legend = bch2d_posterior->GetLegend();
			legend -> SetNColumns(2);
			legend -> AddEntry(bch2d_prior->GetHistogram(), (std::string("prior")+prior_text).data(), "L");
			legend -> AddEntry(bch2d_posterior->GetHistogram(), "posterior", "L");

			bch2d_posterior -> ResizeLegend();
		}

		// Draw legend on top of histogram
		
		legend -> Draw();
		
		// rescale top margin
		gPad -> SetTopMargin(1-legend->GetY1NDC()+0.01);		
	}

	gPad -> RedrawAxis();
	gPad -> Update();
	return 1;
}

// ---------------------------------------------------------
int BCModel::PrintKnowledgeUpdatePlots(const char * filename, unsigned hdiv, unsigned vdiv, bool flag_slice, bool call_likelihood) {
	// prepare prior
	if ( !GetPriorModel(true,call_likelihood) or fPriorModel->GetNParameters() == 0 )
		// return 0 if failed
		return 0;
	
	fPriorModel -> MarginalizeAll();
	fPriorModel -> FindMode();
	
	std::string file(filename);

	// if file extension is neither .pdf nor .ps, force to .pdf
	if ( file.rfind(".pdf") != file.size()-4 and file.rfind(".ps") != file.size()-3 )
		file += ".pdf";

	// create canvas and prepare postscript
	TCanvas * c = new TCanvas(TString::Format("c_%s_update",fPriorModel->GetName().data()));
	c->cd();
	c->Print(std::string(file + "[").c_str());
	
	if (hdiv<1) hdiv = 1;
	if (vdiv<1) vdiv = 1;
	int npads = hdiv * vdiv;
	
	c -> Divide(hdiv,vdiv);
	
	// loop over all parameters and draw 1D plots
	int ndrawn = 0;
	int nprinted = -1;
	c -> cd(1);
	for (unsigned i = 0; i < GetNVariables(); ++i)
		 if(DrawKnowledgeUpdatePlot1D(i, flag_slice, flag_slice)) {
			 ++ndrawn;
			 if (ndrawn!=0 and ndrawn%npads==0) {
				 c -> Print(file.c_str());
				 nprinted = ndrawn;
				 c -> Clear("D");
			 }
			 c -> cd(ndrawn%npads+1);
		 }
	if (nprinted<ndrawn)
		c -> Print(file.c_str());
	
	c -> Clear("D");
	 
	// loop over all parameter pairs
	ndrawn = 0;
	nprinted = -1;
	c -> cd(1) -> Clear();
	for (unsigned i = 0; i < GetNVariables(); ++i)
		for (unsigned j = i+1; j < GetNVariables(); ++j)
			if (DrawKnowledgeUpdatePlot2D(i,j,flag_slice)) {
				++ndrawn;
				if (ndrawn!=0 and ndrawn%npads==0) {
					c -> Print(file.c_str());
					nprinted = ndrawn;
					c -> Clear();
					c -> Divide(hdiv,vdiv);
				}
				c -> cd(ndrawn%npads+1) -> Clear();
			}
	if (nprinted<ndrawn)
		c -> Print(file.c_str());

	c -> Print(std::string(file + "]").c_str());

	// no error
	return 1;
}

// ---------------------------------------------------------
void BCModel::SetKnowledgeUpdateDrawingStyle(BCModel::BCKnowledgeUpdateDrawingStyle style) {
	switch (style) {
	
	case kKnowledgeUpdateDetailedPosterior:
		// 1D
		SetDrawPriorPosteriorNormalOrder(false);
		fBCH1DPriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH1DPriorDrawingOptions -> SetDrawLocalMode(false);
		fBCH1DPriorDrawingOptions -> SetDrawMean(false);
		fBCH1DPriorDrawingOptions -> SetDrawMedian(false);
		fBCH1DPriorDrawingOptions -> SetDrawLegend(false);
		fBCH1DPriorDrawingOptions -> SetNBands(0);
		fBCH1DPriorDrawingOptions -> SetBandType(BCH1D::kNoBands);
		fBCH1DPriorDrawingOptions -> SetROOToptions("same");
		fBCH1DPriorDrawingOptions -> SetLineColor(13);
		fBCH1DPriorDrawingOptions -> SetMarkerColor(13);
		fBCH1DPriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPosteriorDrawingOptions -> CopyOptions(*fBCH1DPriorDrawingOptions);
		fBCH1DPosteriorDrawingOptions -> SetDrawGlobalMode(true);
		fBCH1DPosteriorDrawingOptions -> SetBandType(BCH1D::kSmallestInterval);
		fBCH1DPosteriorDrawingOptions -> SetNBands(3);
		fBCH1DPosteriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPosteriorDrawingOptions -> SetColorScheme(BCHistogramBase::kGreenYellowRed);
		fBCH1DPosteriorDrawingOptions -> SetLineColor(kBlack);
		fBCH1DPosteriorDrawingOptions -> SetMarkerColor(kBlack);

		// 2D
		fBCH2DPriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH2DPriorDrawingOptions -> SetDrawLocalMode(true,false);
		fBCH2DPriorDrawingOptions -> SetDrawMean(false);
		fBCH2DPriorDrawingOptions -> SetDrawLegend(false);
		fBCH2DPriorDrawingOptions -> SetBandType(BCH2D::kSmallestInterval);
		fBCH2DPriorDrawingOptions -> SetBandFillStyle(-1);
		fBCH2DPriorDrawingOptions -> SetNBands(1);
		fBCH2DPriorDrawingOptions -> SetNSmooth(0);
		fBCH2DPriorDrawingOptions -> SetROOToptions("same");
		fBCH2DPriorDrawingOptions -> SetLineColor(13);
		fBCH2DPriorDrawingOptions -> SetMarkerColor(13);
		fBCH2DPriorDrawingOptions -> SetNLegendColumns(1);
		fBCH2DPosteriorDrawingOptions -> CopyOptions(*fBCH2DPriorDrawingOptions);
		fBCH2DPosteriorDrawingOptions -> SetNBands(3);
		fBCH2DPosteriorDrawingOptions -> SetBandFillStyle(1001);
		fBCH2DPosteriorDrawingOptions -> SetColorScheme(BCHistogramBase::kGreenYellowRed);
		fBCH2DPosteriorDrawingOptions -> SetLineColor(kBlack);
		fBCH2DPosteriorDrawingOptions -> SetMarkerColor(kBlack);
		break;
		
	case kKnowledgeUpdateDetailedPrior:
		// 1D
		SetDrawPriorPosteriorNormalOrder(true);
		fBCH1DPosteriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH1DPosteriorDrawingOptions -> SetDrawLocalMode(false);
		fBCH1DPosteriorDrawingOptions -> SetDrawMean(false);
		fBCH1DPosteriorDrawingOptions -> SetDrawMedian(false);
		fBCH1DPosteriorDrawingOptions -> SetDrawLegend(false);
		fBCH1DPosteriorDrawingOptions -> SetNBands(0);
		fBCH1DPosteriorDrawingOptions -> SetBandType(BCH1D::kNoBands);
		fBCH1DPosteriorDrawingOptions -> SetROOToptions("same");
		fBCH1DPosteriorDrawingOptions -> SetLineColor(13);
		fBCH1DPosteriorDrawingOptions -> SetMarkerColor(13);
		fBCH1DPosteriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPriorDrawingOptions -> CopyOptions(*fBCH1DPosteriorDrawingOptions);
		fBCH1DPriorDrawingOptions -> SetDrawGlobalMode(true);
		fBCH1DPriorDrawingOptions -> SetBandType(BCH1D::kSmallestInterval);
		fBCH1DPriorDrawingOptions -> SetNBands(3);
		fBCH1DPriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPriorDrawingOptions -> SetColorScheme(BCHistogramBase::kGreenYellowRed);
		fBCH1DPriorDrawingOptions -> SetLineColor(kBlack);
		fBCH1DPriorDrawingOptions -> SetMarkerColor(kBlack);

		// 2D
		fBCH2DPosteriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH2DPosteriorDrawingOptions -> SetDrawLocalMode(true,false);
		fBCH2DPosteriorDrawingOptions -> SetDrawMean(false);
		fBCH2DPosteriorDrawingOptions -> SetDrawLegend(false);
		fBCH2DPosteriorDrawingOptions -> SetBandType(BCH2D::kSmallestInterval);
		fBCH2DPosteriorDrawingOptions -> SetBandFillStyle(-1);
		fBCH2DPosteriorDrawingOptions -> SetNBands(1);
		fBCH2DPosteriorDrawingOptions -> SetNSmooth(0);
		fBCH2DPosteriorDrawingOptions -> SetROOToptions("same");
		fBCH2DPosteriorDrawingOptions -> SetLineColor(13);
		fBCH2DPosteriorDrawingOptions -> SetMarkerColor(13);
		fBCH2DPosteriorDrawingOptions -> SetNLegendColumns(1);
		fBCH2DPriorDrawingOptions -> CopyOptions(*fBCH2DPosteriorDrawingOptions);
		fBCH2DPriorDrawingOptions -> SetNBands(3);
		fBCH2DPriorDrawingOptions -> SetColorScheme(BCHistogramBase::kGreenYellowRed);
		fBCH2DPriorDrawingOptions -> SetBandFillStyle(1001);
		fBCH2DPriorDrawingOptions -> SetLineColor(kBlack);
		fBCH2DPriorDrawingOptions -> SetMarkerColor(kBlack);
		break;

	case kKnowledgeUpdateDefaultStyle:
	default:
		// 1D
		fBCH1DPriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH1DPriorDrawingOptions -> SetDrawLocalMode(false);
		fBCH1DPriorDrawingOptions -> SetDrawMean(false);
		fBCH1DPriorDrawingOptions -> SetDrawMedian(false);
		fBCH1DPriorDrawingOptions -> SetDrawLegend(false);
		fBCH1DPriorDrawingOptions -> SetNBands(0);
		fBCH1DPriorDrawingOptions -> SetBandType(BCH1D::kNoBands);
		fBCH1DPriorDrawingOptions -> SetROOToptions("same");
		fBCH1DPriorDrawingOptions -> SetLineColor(kRed);
		fBCH1DPriorDrawingOptions -> SetMarkerColor(kRed);
		fBCH1DPriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPosteriorDrawingOptions -> CopyOptions(*fBCH1DPriorDrawingOptions);
		fBCH1DPosteriorDrawingOptions -> SetNLegendColumns(1);
		fBCH1DPosteriorDrawingOptions -> SetLineColor(kBlue);
		fBCH1DPosteriorDrawingOptions -> SetMarkerColor(kBlue);

		// 2D
		fBCH2DPriorDrawingOptions -> SetDrawGlobalMode(false);
		fBCH2DPriorDrawingOptions -> SetDrawLocalMode(true,false);
		fBCH2DPriorDrawingOptions -> SetDrawMean(false);
		fBCH2DPriorDrawingOptions -> SetDrawLegend(false);
		fBCH2DPriorDrawingOptions -> SetBandType(BCH2D::kSmallestInterval);
		fBCH2DPriorDrawingOptions -> SetBandFillStyle(-1);
		fBCH2DPriorDrawingOptions -> SetNBands(1);
		fBCH2DPriorDrawingOptions -> SetNSmooth(0);
		fBCH2DPriorDrawingOptions -> SetROOToptions("same");
		fBCH2DPriorDrawingOptions -> SetLineColor(kRed);
		fBCH2DPriorDrawingOptions -> SetMarkerColor(kRed);
		fBCH2DPriorDrawingOptions -> SetNLegendColumns(1);
		fBCH2DPosteriorDrawingOptions -> CopyOptions(*fBCH2DPriorDrawingOptions);
		fBCH2DPosteriorDrawingOptions -> SetLineColor(kBlue);
		fBCH2DPosteriorDrawingOptions -> SetMarkerColor(kBlue);
		break;
	}
}
