//
// This ROOT macro shows how one can easily analyze the MCMC stored
// in the ROOT file using the default rootOutput example from the BAT
// distribution. The macro is part of BAT and can only be run if BAT
// was installed correctly.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x analyzeMCMCTree.C
//
// or
//
//    root[1] .L analyzeMCMCTree.C
//    root[2] analyzeMCMCTree()
//
// or from the command line
//
//    $ root analyzeMCMCTree.C
//
// The macro accepts two optional arguments:
//   <filename>   the name of the ROOT file containing the MCMC
//   <nchains>    number of chains stored in the file
// The defaults for the two arguments are set to the default settings
// of the rootOutput example.
//
// Macro shows how to produce simple marginalized distributions for
// 1 or 2 parameters. It also shows how to obtain the probability
// distributions for any quantity derived from the parameters of the
// model.

// ---------------------------------------------------------
analyzeMCMCTree(const char *filename="output.root", int nchains=5)
{
   // open root file
   TFile *rfile = TFile::Open(filename);

   // prepare some variables
   // pointers to chains
   TChain **ch = new TChain*[nchains];
   // marginalized distributions
   TH1D *h1d0; // 1D: par1
   TH1D *h1d1; // 1D: par2
   TH1D *h1dx; // 1D: function of par1 and par2
   TH2D *h2d;  // 2D: par1 vs. par2

   // read all chains from the file and fill marginaliyed distributions
   for (int ichain=0; ichain < nchains; ichain++) {

      cout<<TString::Format("MarkovChainTree_%d",ichain)<<endl;

      // read chain
      ch[ichain] = rfile->Get(TString::Format("MarkovChainTree_%d",ichain));

      // extract 1D and 2D marginalized distributions after the chains have converged (Phase==2)
      ch[ichain]->Draw("Parameter0>>h1","Phase==2");
      ch[ichain]->Draw("Parameter1>>h2","Phase==2");
      ch[ichain]->Draw("log(Parameter0*Parameter0)>>h3","Phase==2"); // fill f(par1,par2) = log(par0^2)
      ch[ichain]->Draw("Parameter0:Parameter1>>h4","Phase==2");
      if (ichain==0) {
         h1d0 = (TH1D*)h1;
         h1d0->GetXaxis()->SetTitle("Parameter0");
         h1d1 = (TH1D*)h2;
         h1d1->GetXaxis()->SetTitle("Parameter1");
         h1dx = (TH1D*)h3;
         h1dx->GetXaxis()->SetTitle("f(Parameter0,Parameter1)");
         h2d = (TH2D*)h4;
         h2d->GetXaxis()->SetTitle("Parameter1");
         h2d->GetYaxis()->SetTitle("Parameter0");
      }
      // if there's more than one chain, add the marginaliyed distributions together
      else {
         h1d0->Add(h1);
         h1d1->Add(h2);
         h1dx->Add(h3);
         h2d->Add(h4);
      }
   }

   // remove the stats drawing
   h1d0->SetStats(0);
   h1d1->SetStats(0);
   h1dx->SetStats(0);
   h2d->SetStats(0);

   h2d->SetTitle("");

   // create BAT specific histogram objects
   BCH1D * bh1d0 = new BCH1D(h1d0);
   BCH1D * bh1d1 = new BCH1D(h1d1);
   BCH1D * bh1dx = new BCH1D(h1dx);
   BCH2D * bh2d = new BCH2D(h2d);

   TCanvas * c = new TCanvas();
   c->Divide(2,2);

   // draw marginalized distributions in BAT style
   c->cd(1);
   bh1d0->Draw();
   c->cd(2);
   bh1d1->Draw();
   c->cd(3);
   bh1dx->Draw();
   c->cd(4);
   bh2d->Draw(2);

}

// ---------------------------------------------------------
