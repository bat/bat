#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCEmptyModel.h>
#include <BAT/BCLog.h>

#include <TCanvas.h>

int main()
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new BCEmptyModel that reads in file
    // empty string for second argument tells BAT to search for model in file.
    // If you run the example repeatedly, results are added.
    // To avoid errors, we give the name of the model explicitly.
    BCEmptyModel m("GaussModel.root", "gausMod");

    m.Remarginalize();

    // Draw the 1D distributions and the x:y 2D distribution
    // to a new canvas.
    TCanvas* C = new TCanvas();
    C->Divide(2, 2);

    C->cd(1);
    BCH1D h0 = m.GetMarginalized(0);
    h0.Draw();

    C->cd(2);
    BCH1D h1 = m.GetMarginalized(1);
    h1.Draw();

    C->cd(3);
    BCH1D h2 = m.GetMarginalized(2);
    h2.Draw();

    C->cd(4);
    BCH2D h01 = m.GetMarginalized(0, 1);
    h01.Draw();

    C->Update();
    C->Print("GaussModel_loaded_plots.pdf");

    // close log file
    BCLog::CloseLog();

    return 0;
}
