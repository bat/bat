/*
 * Copyright (C) 2013, Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GaussModel.h"
#include "test.h"

#include <BAT/BCParameter.h>
#include <BAT/BCSummaryTool.h>

using namespace test;

class BCSummaryToolTest :
   public TestCase
{
public:
   BCSummaryToolTest():
      TestCase("BCSummaryTool"),
      file_corr("BCSummaryTool_TEST_summary_corr.pdf"),
      file_par("BCSummaryTool_TEST_summary_par.pdf"),
      file_upd("BCSummaryTool_TEST_summary_upd.pdf")
   {
   }

   void print(BCModel & m) const
   {
      m.MarginalizeAll();
      m.FindMode();
			std::cout << "SUMMARY" << std::endl;
      BCSummaryTool s(&m);
			std::cout << "MODEL INSTANTIATED" << std::endl;
      s.PrintCorrelationPlot(file_corr.c_str());
			std::cout << "CORRELATIONS PRINTED" << std::endl;
      s.PrintParameterPlot(file_par.c_str());
			std::cout << "PARAMETERPLOT PRINTED" << std::endl;
      s.PrintKnowledgeUpdatePlots(file_upd.c_str());
			std::cout << "KNOWLEDGEUPDATE PRINTED" << std::endl;
   }

   void mcmc() const
   {
      GaussModel m("mcmc", 3);
      m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
      print(m);
   }

   void fixed() const
   {
      GaussModel m("fixed", 3);
      m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
      m.GetParameter(1)->Fix(6);
      print(m);
   }

   void slice() const
   {
      GaussModel m("slice", 2);
      m.SetMarginalizationMethod(BCIntegrate::kMargGrid);
      print(m);
   }

   virtual void run() const
   {
      mcmc();
			std::cout << "MCMC DONE" << std::endl;
      fixed();
			std::cout << "FIXED DONE" << std::endl;
      slice();
			std::cout << "SLICE DONE" << std::endl;
   }

private:
   std::string file_corr, file_par, file_upd;
} bcsummarytool_test;
