#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCEmptyModel.h>
#include <BAT/BCLog.h>

#include <TCanvas.h>

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // BCEmptyModel reads in file
    // leaving out the second argument tells BAT to search for a model in the file.
    BCEmptyModel m("GaussModel.root");

    m.Remarginalize();

    // Draw the 1D distributions and the x:y 2D distribution
    // to a new canvas.
    TCanvas C;
    C.Divide(2, 2);

    C.cd(1);
    BCH1D h0 = m.GetMarginalized(0);
    h0.Draw();

    C.cd(2);
    BCH1D h1 = m.GetMarginalized(1);
    h1.Draw();

    C.cd(3);
    BCH1D h2 = m.GetMarginalized(2);
    h2.Draw();

    C.cd(4);
    BCH2D h01 = m.GetMarginalized(0, 1);
    h01.Draw();

    C.Update();
    C.Print((m.GetSafeName() + "_loaded_plots.pdf").data());

    // close log file
    BCLog::CloseLog();

    return 0;
}
