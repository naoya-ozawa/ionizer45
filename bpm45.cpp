// A script to 
// 1. define a plane P(w) such that z = -x + sqrt(2)*w in Inventor coordinates
// 2. show the beam profile at the selected P(w) and calculate the center/RMS
// Developed from the focus_analysis.cpp
// exclusively for the design of the CNS FrEDM thermal ionizer

// The trajectory should fly in the z-x plane in the (1, 0 , 1) direction in SIMION/Inventor coordinates
// The output file (testplane-bpm.csv) should be in the form of 
// | ion # | X position | Y position | Z position |
// and the first line should be the header and the last line should be some character that is not a number


#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TEventList.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>

#include <TH2D.h>
#include <TLatex.h>

#include <TRint.h>

using namespace std;

double position_on_target(double x_simion, double y_simion, double xcent_simion, double ycent_simion, const char* x_or_y){
    // Returns the position of the generated ions on the surface of the target
    // as if it is obsereved by an MCP right above the target
    // Given the center point in SIMION coordinates, it is converted in the following way:
    // 1. revolve -90deg around z-axis, assuming (xcent_simion,ycent_simion,zcent_simion+h)_(simion) = (0,0,0)_(inventor)
    // Returns -9999. for any invalid input
    if ( x_or_y == "x" ){
        return y_simion - ycent_simion;
    }else if ( x_or_y == "y" ){
        return xcent_simion - x_simion;
    }else{
        return -9999.;
    }
}

double dist_testplane(double w, double x_0, double x_d, double z_0, double z_d, double height){
    // Checks whether a point (trkx,trky,trkz)_(SIMION) is past the plane P(w)
    // Returns the (signed) distance from P(w)
    // given that the ion is ejected from (sx,sy,sz-height)_(SIMION)
    // positive: past P(w), 0: on P(w), negative: not past P(w)
    // The target surface is located <height> mm below the origin (in Inventor coordinates)
    double w_trk = ( (x_0 + x_d + z_0 + z_d)/TMath::Sqrt(2.0) ) - w;
//    cout << "w_trk = " << w_trk << " mm around (" << x_0 << ", " << z_0 << ")" << endl;
    return w_trk;
}

double t_trj(double w, double x_0, double x_d, double z_0, double z_d){
    // Returns the parameter to define the intersection between the trajectory and P(w)
    // In the case of x_d + z_d = 0 (the ion is making no progress), returns -9999.
    if ( x_d+z_d == 0){
        cout << "Caught ion making no progress around (" << x_0 << ", " << z_0 << ")_(inventor)" << endl;
        return -9999.;
    }else{
        return ( (TMath::Sqrt(2.0)*w) - x_0 - z_0 ) / ( x_d + z_d );
    }
}

double position_on_bpm(double x_inv, double y_inv, double z_inv, const char* x_or_y){
    // Returns the position of the trajectory point t_trj(w) = (trjptx,trjpty,trjptz)_(inventor) on the MCP
    // Converted in the following way:
    // 1. shift origin vec0 = (0,0,0)_(inventor) to BPM(w) = (w/2sqrt2 , 0 , w/2sqrt2)_(inventor)
    //    (trjptx,trjpty,trjptz)_(inventor) = (trjptx-w/2sqrt2,trjpty,trjptz-w/2sqrt2)_(P)
    // 2. revolve +45deg around y-axis
    // 3. revolve -90deg around z-axis
    // Returns -9999. for any invalid input
    if ( x_or_y == "x" ){
        return y_inv;
    }else if ( x_or_y == "y" ){
        return (z_inv-x_inv)/TMath::Sqrt(2.0);
    }else{
        return -9999.;
    }
}


int main (int argc, char** argv){

    // Check input number here (w must be an input parameter)
    if (argc != 2){
        cout << "usage: ./bpm45 <w>" << endl;
        exit(1);
    }else{	
        double w = strtod(argv[1],NULL);
        cout << "At w = " << w << " mm:" << endl;

        TRint rootapp("app",&argc,argv);

        TCanvas *c1 = new TCanvas();
        c1->Divide(2,2);

        // SIMION output CSV file
        // in the form of
        // | ion# | posX | posY | posZ |
        const char* zx_file = "./../testplane-bpm.csv";
        string line;

        // Data will be stored to a ROOT file as a TTree with step number given to each step
        TFile data_table("trajectory.root","RECREATE","All step data from SIMION");

        int ion_number = 0;
        int step_number = 0;
        double X_position = 0.0;
        double Y_position = 0.0;
        double Z_position = 0.0;

        TTree *zxfocus = new TTree("trajectory","ion_number:step_number:X_position:Y_position:Z_position");
        zxfocus->Branch("ion_number",&ion_number,"ion_number/I");
        zxfocus->Branch("step_number",&step_number,"step_number/I");
        zxfocus->Branch("X_position",&X_position,"X_position/D");
        zxfocus->Branch("Y_position",&Y_position,"Y_position/D");
        zxfocus->Branch("Z_position",&Z_position,"Z_position/D");

        // Beam Profile Monitoring resolution
        int pix = 65;

        // Target lowered from the center
        double h = 30.0; // mm

        // Convert CSV to ROOT
        ifstream simion_output( zx_file );
        if (!simion_output){
            cout << "Cannot open SIMION output!" << endl;
            exit(1);
        }else{
            cout << "SIMION output opened." << endl;
        
            simion_output.ignore(100,'\n'); // skip the header

            int k = 0; // line counter

            while ( getline(simion_output,line) ){ // read loop for the entire CSV file

//                cout << "Read " << line << endl;            
                istringstream linestream(line);
                string item;
                while ( getline(linestream,item,',') ){ // read loop for each line
                    ++k;
                    try{
                        int tester = stoi(item);

                        if (k%4 == 1){        
                            if (stoi(item) != ion_number){
                                step_number = 0;
                            }else{
                                step_number = step_number + 1;
                            }
                            ion_number = stoi(item);   
                        }else if (k%4 == 2){
                            X_position = stod(item);
                        }else if (k%4 == 3){
                            Y_position = stod(item);
                        }else{
                            Z_position = stod(item);
                        }
                    }
                    catch (const invalid_argument& e){
                        // catch the EOF
                        cout << "Read " << k << " lines for " << ion_number << " flown ions." << endl;
                        k = -1;
                    }
                }
                if (k != -1){
                    // Don't fill at the last line 
                    zxfocus->Fill();
                }
            }
        }
    
        data_table.Write();
        // Finished conversion to ROOT TTree






        // Find the number of flown ions    
        int Nions = ion_number;
        cout << Nions << " ions flown" << endl;




        // Find the source center point and draw upper-view 2d histo (x-->-y, y-->x)
        // This will define the origin in Inventor coordinates
        c1->cd(1);
        zxfocus->Draw(">>src_elist","step_number == 0","goff");
        TEventList *src_elist = (TEventList *)gDirectory->Get("src_elist");
        double src_CPx = 0.0;
        double src_CPy = 0.0;
        double src_CPz = 0.0;
        double src_RMSx = 0.0;
        double src_RMSy = 0.0;
        TH2D *src = new TH2D("src","Generated Ions at Target;X (mm);Y (mm)",pix,-15.5,15.5,pix,-15.5,15.5);
        for (int i = 0; i < Nions; ++i){
            zxfocus->GetEntry(src_elist->GetEntry(i));
            src_CPx += X_position;
            src_CPy += Y_position;
            src_CPz += Z_position;
        }
        src_CPx /= Nions;
        src_CPy /= Nions;
        src_CPz /= Nions;

        for (int i = 0; i < Nions; ++i){
            zxfocus->GetEntry(src_elist->GetEntry(i));
//            src->Fill(Y_position-src_CPy,src_CPx-X_position);
            double srcx = position_on_target(X_position,Y_position,src_CPx,src_CPy,"x");
            double srcy = position_on_target(X_position,Y_position,src_CPx,src_CPy,"y");
            src->Fill(srcx,srcy);
        }
        src->Draw("colz");


        c1->cd(3);
        for (int i = 0; i < Nions; ++i){
            zxfocus->GetEntry(src_elist->GetEntry(i));
            double srcx = position_on_target(X_position,Y_position,src_CPx,src_CPy,"x");
            double srcy = position_on_target(X_position,Y_position,src_CPx,src_CPy,"y");
            src_RMSx += srcx*srcx;
            src_RMSy += srcy*srcy;
        }
        src_RMSx = TMath::Sqrt(src_RMSx/Nions);
        src_RMSy = TMath::Sqrt(src_RMSy/Nions);
//    cout << "Ion Source at (" << src_CPx << ", " << src_CPy << ", " << src_CPz << ") (mm in SIMION coordinates)" << endl;  

        TLatex l_src;
        l_src.SetTextAlign(12);
        l_src.SetTextSize(0.05);
        l_src.DrawLatex(0.15,0.9,"Ion Source");
        l_src.DrawLatex(0.15,0.8,"MEANx = 0 [mm]"); // By definition
        l_src.DrawLatex(0.15,0.7,"MEANy = 0 [mm]"); // By definition
        l_src.DrawLatex(0.15,0.6,Form("RMSx = %g [mm]",src_RMSx));
        l_src.DrawLatex(0.15,0.5,Form("RMSy = %g [mm]",src_RMSy));
        l_src.DrawLatex(0.15,0.4,Form("Number of Flown Particles = %d",Nions));





        // Show beam profile at the given w with end view (y-->x, +z-x-->y) 
        c1->cd(2);

        double stepfrontX = 0.0;
        double stepbackX = 0.0;
        double stepfrontY = 0.0;
        double stepbackY = 0.0;
        double stepfrontZ = 0.0;
        double stepbackZ = 0.0;

        double trjptx[Nions];
        double trjpty[Nions];
        double trjptz[Nions];
        for (int i = 0; i < Nions; ++i){
            trjptx[i] = -9999.;
            trjpty[i] = -9999.;
            trjptz[i] = -9999.;
        }


        TH2D *mcp = new TH2D("mcp","Beam at P(w);X (mm);Y (mm)",pix,-15.5,15.5,pix,-15.5,15.5);
        double mcp_cpx = 0.0;
        double mcp_cpy = 0.0;
        double mcp_rmsx = 0.0;
        double mcp_rmsy = 0.0;
        int mcp_hits = 0;

        for (int i = 0; i < Nions; ++i){

            zxfocus->Draw(">>trj",Form("ion_number == %d",i+1),"goff");
            TEventList *trj = (TEventList *)gDirectory->Get("trj");
            for (int j = 0; j < trj->GetN(); ++j){
                zxfocus->GetEntry(trj->GetEntry(j));
                // Store the position for the "front"/current point vec{P_f}
                stepfrontX = X_position;
                stepfrontY = Y_position;
                stepfrontZ = Z_position;

                // vec{P_b}
                double x_0 = stepbackX - src_CPx;
                double y_0 = stepbackY - src_CPy;
                double z_0 = stepbackZ - src_CPz - h;
                // vec{d} = vec{P_f} - vec{P_b}
                double x_d = stepfrontX - stepbackX;
                double y_d = stepfrontY - stepbackY;
                double z_d = stepfrontZ - stepbackZ;


                // Check the position wrt P(w)
                double d_p = dist_testplane(w,x_0,x_d,z_0,z_d,h);
     
                if ( d_p > 0.0 ){// The ion passed P(w)
                    // t_trj
                    trjptx[i] = t_trj(w,x_0,x_d,z_0,z_d)*x_d + x_0;
                    trjpty[i] = t_trj(w,x_0,x_d,z_0,z_d)*y_d + y_0;
                    trjptz[i] = t_trj(w,x_0,x_d,z_0,z_d)*z_d + z_0;
                    ++mcp_hits;
                    cout << "Ion #" << i+1 << " passed P(w) at (" << trjptx[i] << ", " << trjpty[i] << ", " << trjptz[i] << ")_(inventor)" << endl;
                    // connect the two steps with a line and find the intercept with P(w)
                    // use the intercept for RMS and CP calculation
                    break;
                }else if ( d_p == 0.0 ){// The ion is at P(w) <=> t_trj = 1
                    trjptx[i] = x_d + x_0;
                    trjpty[i] = y_d + y_0;
                    trjptz[i] = z_d + z_0;
                    ++mcp_hits;
                    cout << "Ion #" << i+1 << " hit P(w) at (" << trjptx[i] << ", " << trjpty[i] << ", " << trjptz[i] << ")_(inventor)" << endl;
                    break;
                }else{
                    // Store the position for the "back"/previous point vec{P_b}
                    stepbackX = X_position;
                    stepbackY = Y_position;
                    stepbackZ = Z_position;
                }
            }

            if (trjptx[i] != -9999.0){
                // Record only the ions that reached P(w)
                double bpmx = position_on_bpm(trjptx[i],trjpty[i],trjptz[i],"x");
                double bpmy = position_on_bpm(trjptx[i],trjpty[i],trjptz[i],"y");
                mcp->Fill( bpmx , bpmy );
                mcp_cpx += bpmx;
                mcp_cpy += bpmy;
            }else{
                cout << "Ion #" << i+1 << " did not survive!" << endl;
            }

        }
        // Find the beam center
        mcp_cpx /= mcp_hits;
        mcp_cpy /= mcp_hits;

        // Calculate the RMS
        for (int i = 0; i < Nions; ++i){
            if (trjptx[i] != -9999.0){
                double bpmx = position_on_bpm(trjptx[i],trjpty[i],trjptz[i],"x");
                double bpmy = position_on_bpm(trjptx[i],trjpty[i],trjptz[i],"y");
                mcp_rmsx += (bpmx - mcp_cpx)*(bpmx - mcp_cpx);
                mcp_rmsy += (bpmy - mcp_cpy)*(bpmy - mcp_cpy);
            }
        }
        mcp_rmsx = TMath::Sqrt(mcp_rmsx/mcp_hits);
        mcp_rmsy = TMath::Sqrt(mcp_rmsy/mcp_hits);

        mcp->Draw("colz");


        c1->cd(4);
   
        TLatex l_mcp;
        l_mcp.SetTextAlign(12);
        l_mcp.SetTextSize(0.05);
        l_mcp.DrawLatex(0.15,0.9,Form("BPM at w = %g [mm]",w));
        l_mcp.DrawLatex(0.15,0.8,Form("MEANx = %g [mm]",mcp_cpx));
        l_mcp.DrawLatex(0.15,0.7,Form("MEANy = %g [mm]",mcp_cpy));
        l_mcp.DrawLatex(0.15,0.6,Form("RMSx = %g [mm]",mcp_rmsx));
        l_mcp.DrawLatex(0.15,0.5,Form("RMSy = %g [mm]",mcp_rmsy));
        l_mcp.DrawLatex(0.15,0.4,Form("Transmission rate %g%%",100.*double(mcp_hits)/double(Nions)));

        c1->Update();
        c1->Modified();

        rootapp.Run();

        data_table.Close();
    }

    return 0;

}
