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
    m.GetMarginalized(0).Draw();
    C->cd(2);
    m.GetMarginalized(1).Draw();
    C->cd(3);
    m.GetMarginalized(2).Draw();
    C->cd(4);
    m.GetMarginalized(0, 1).Draw();

    C->Print("GaussModel_loaded_plots.pdf");

    // close log file
    BCLog::CloseLog();

    return 0;
}
