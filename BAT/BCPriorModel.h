#ifndef __BCPRIORMODEL__H
#define __BCPRIORMODEL__H

/*!
 * \class BCPriorModel
 * \brief Class for sampling from prior of a BCModel
 * \author Daniel Greenwald
 * \version 1.0
 * \date 09.2014
 * \detail This class acts as a BCModel using the prior of another BCModel as its posterior
 * for the purpose of knowledge-update plotting.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <limits>

#include "BCModel.h"

// ---------------------------------------------------------

class BCPriorModel : public BCModel
{

public:

	/**
	 * constructor. */
	BCPriorModel(BCModel * model);

	/**
	 * destructor. */
	~BCPriorModel();

	/**
	 * Returns a constant prior. */
	double LogAPrioriProbability(const std::vector<double> &parameters)
	{ return 0; }

   /**
		* Returns prior of model as posterior of PriorModel. */
	double LogLikelihood(const std::vector<double> &parameters)
	{ return (fModel) ? fModel->LogAPrioriProbability(parameters) : -std::numeric_limits<double>::infinity(); }

	/**
	 * Calculates user observables according to the model. */
	void CalculateObservables(const std::vector<double> &parameters)
	{ if (fModel) fModel->CalculateObservables(parameters); }

	/**
	 * Prepare PriorModel from Model. */
	bool PreparePriorModel();

protected:
	BCModel * fModel;

	/* Hide BCModel functions related to BCPriorModel. */
	BCPriorModel * GetPriorModel(bool prepare=true) {return 0;}
	int DrawKnowledgeUpdatePlot1D(unsigned index, std::string options_post="", std::string options_prior="") {return 0;}
	int DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice=false, double interval_content=68e-2) {return 0;}
	int PrintKnowledgeUpdatePlots(const char * filename = "update.pdf", unsigned hdiv=1, unsigned vdiv=1, std::string options="-", double interval_content=68e-2) {return 0;}



};

#endif
