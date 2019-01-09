// A script to display the trajectory as a Z-X plot and X-Y plot
// Recieves the CSV data generated with SIMION and displays all the points on a TGraph

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TRint.h>

using namespace std;


int main(int argc, char** argv){

    TRint rootapp("app",&argc,argv);

    TCanvas *c1 = new TCanvas();
    c1->Divide(2,1);
    
    // SIMION output CSV file
    // in the form of
    // | ion# | posX | posy | posZ |
    const char* zx_file = "./../testplane-bpm.csv";
    string line;
    int ion_number;

    // Specify the origin so that the trajectory is displayed in Inventor coordinates
    double sx = 200.0;
    double sy = 210.0;
    double sz = 213.85;

    // Create graphs to draw
    TGraphErrors *traj_zx = new TGraphErrors();
    traj_zx->SetTitle("Trajectory in ZX View; X_{Inventor} (mm); Z_{Inventor} (mm)");
    TGraphErrors *traj_xy = new TGraphErrors();
    traj_xy->SetTitle("Trajectory in XY View; X_{Inventor} (mm); Y_{Inventor} (mm)");


    // Write out CSV file to graph
    ifstream simion_output( zx_file );
    if (!simion_output){
        cout << "Cannot open SIMION output!" << endl;
        exit(1);
    }else{
        simion_output.ignore(100,'\n');

        int k = 0; // line counter

        while ( getline(simion_output,line) ){ // read loop over the entire CSV file
            istringstream linestream(line);
            string item;
            while ( getline(linestream,item,',') ){ // read loop for each line
                ++k;
                try{
                    int tester = stoi(item);
                    double X_position, Y_position, Z_position;
                    if (k%4 == 1){
                        ion_number = stod(item);
//                        cout << "Reading line " << k << endl;
                    }else if (k%4 == 2){
                        X_position = stod(item);
                    }else if (k%4 == 3){
                        Y_position = stod(item);
                    }else{
                        Z_position = stod(item);
                        traj_zx->SetPoint(k/4 - 1,X_position-sx,Z_position-sz);
                        traj_xy->SetPoint(k/4 - 1,X_position-sx,Y_position-sy);
//                        cout << "Ion #" << ion_number << " is at (" << X_position-sx << ", " << Y_position-sy << ", " << Z_position-sz << ")_(Inventor)" << endl;
                    }
                }
                catch (const invalid_argument& e){
                    cout << "Read " << k << " lines for " << ion_number << " ions" << endl;
                    k = -1;
                }
            }   
        }




        c1->cd(1);
        traj_zx->Draw("AP");
        TF1 *extractionzx = new TF1("extractionzx","x",-10.,300.);
        extractionzx->SetLineStyle(3);
        extractionzx->SetLineColor(kRed);
        extractionzx->Draw("SAME");

        TF1 *extractionfit = new TF1("extractionfit","[0]*x+[1]",50.,150.);
        extractionfit->SetParameters(1.0,0.0);
        extractionfit->SetLineColor(kBlue);
        traj_zx->Fit("extractionfit","M","",50.,150.);
        extractionfit->Draw("SAME");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(12);
        tl->SetTextSize(0.03);
        tl->DrawLatex(200.,50.,"Linear Fit");
        tl->DrawLatex(200.,25.,Form("grad=(%g+-%g)",extractionfit->GetParameter(0),extractionfit->GetParError(0)));
        tl->DrawLatex(200.,0.,Form("ofs=(%g+-%g)",extractionfit->GetParameter(1),extractionfit->GetParError(1)));

        c1->cd(2);
        traj_xy->Draw("AP");
        TF1 *extractionxy = new TF1("extractionxy","0.0",-10.,300.);
        extractionxy->SetLineStyle(3);
        extractionxy->SetLineColor(kRed);
        extractionxy->Draw("SAME");

    
        c1->Update();
        c1->Modified();

        rootapp.Run();
    }


    return 0;
}
