#ifndef __BCFITTER__H
#define __BCFITTER__H

/*!
 * \class BCFitter
 * \brief A base class for all fitting classes
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2013
 * \detail A base class for all fitting classes
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "../../BAT/BCModel.h"

// ---------------------------------------------------------

class BCFitter : public BCModel
{
 public:

  /** \name Constructors and destructors */
  /* @{ */

  /**
   * The default constructor. */
  BCFitter();

  /**
   * Constructor
   * @param name name of the model */
  BCFitter(const char * name);

  /**
   * The default destructor. */
  ~BCFitter();

  /* @} */

  /** \name Member functions (get) */
  /* @{ */

  /**
   * @return pointer to the error band */
  TGraph * GetErrorBand()
  { return fErrorBand; };

  /**
     const BCParameter * GetParameter(const char * name);
     * @return The 2-d histogram of the error band. */
  TH2D * GetErrorBandXY() const
  { return fErrorBandXY; }

  TH2D * GetErrorBandXY_yellow(double level=.68, int nsmooth=0) const;

  /**
   * Returns a vector of y-values at a certain probability level.
   * @param level The level of probability
   * @return vector of y-values */
  std::vector<double> GetErrorBand(double level) const;

  TGraph * GetErrorBandGraph(double level1, double level2) const;

  TGraph * GetFitFunctionGraph(const std::vector<double> &parameters);

  TGraph * GetFitFunctionGraph()
  { return GetFitFunctionGraph(GetBestFitParameters()); }

  TGraph * GetFitFunctionGraph(const std::vector<double> &parameters, double xmin, double xmax, int n=1000);

  void FixDataAxis(unsigned int index, bool fixed);

  bool GetFixedDataAxis(unsigned int index) const;

  /* @} */

  /** \name Member functions (set) */
  /* @{ */

  /**
   * Turn on or off the filling of the error band during the MCMC run.
   * @param flag set to true for turning on the filling */
  void SetFillErrorBand(bool flag = true)
  { fFlagFillErrorBand=flag; }

  /**
   * Sets errorband histogram */
  void SetErrorBandHisto(TH2D * h)
  { fErrorBandXY = h; }

  /**
   * Turn off filling of the error band during the MCMC run.
   * This method is equivalent to SetFillErrorBand(false) */
  void UnsetFillErrorBand()
  { fFlagFillErrorBand=false; }

  /**
   * Sets index of the x values in function fits.
   * @param index Index of the x values */
  void SetFitFunctionIndexX(int index)
  { fFitFunctionIndexX = index; }

  /**
   * Sets index of the y values in function fits.
   * @param index Index of the y values */
  void SetFitFunctionIndexY(int index)
  { fFitFunctionIndexY = index; }

  /**
   * Sets indices of the x and y values in function fits.
   * @param indexx Index of the x values
   * @param indexy Index of the y values */
  void SetFitFunctionIndices(int indexx, int indexy)
  { SetFitFunctionIndexX(indexx);
    SetFitFunctionIndexY(indexy); }

  /**
   * Sets the error band flag to continuous function */
  void SetErrorBandContinuous(bool flag);

  /**
   * Defines a fit function.
   * @param parameters A set of parameter values
   * @param x A vector of x-values
   * @return The value of the fit function at the x-values given a set of parameters */
  virtual double FitFunction(const std::vector<double> &/*x*/, const std::vector<double> &/*parameters*/)
  { return 0; }

  /* @} */
  /** \name Member functions (miscellaneous methods) */
  /* @{ */

  /**
   * Read */
  int ReadErrorBandFromFile(const char * file);

  /**
   * Performs the fit.
   * @return An error code. */
  virtual int Fit() = 0;

  /**
   * Draw the fit in the current pad. */
  virtual void DrawFit(const char * options, bool flaglegend = false) = 0;

  /**
   * Overloaded from BCEngineMCMC */
  void MCMCIterationInterface();

  /**
   * Overloaded from BCIntegrate. */
  void MarginalizePreprocess();

  /**
   * Overloaded from BCIntegrate. */
  void MarginalizePostprocess()
  {;};

  /**
   * Fill error band histogram for curreent iteration. This method is called from MCMCIterationInterface() */
  void FillErrorBand();

  /* @} */

 private:

  /**
   * Pointer to the error band (for legend) */
  TGraph * fErrorBand;

  /**
   * Flag whether or not to fill the error band */
  bool fFlagFillErrorBand;

  /**
   * The indices for function fits */
  int fFitFunctionIndexX;
  int fFitFunctionIndexY;

  /**
   * A flag for single point evaluation of the error "band" */
  bool fErrorBandContinuous;
  std::vector<double> fErrorBandX;

  /**
   * Number of X bins of the error band histogram */
  unsigned fErrorBandNbinsX;

  /**
   * Number of Y bins of the error band histogram */
  unsigned fErrorBandNbinsY;

 protected:

  /**
   * The error band histogram */
  TH2D * fErrorBandXY;

};

// ---------------------------------------------------------

#endif
