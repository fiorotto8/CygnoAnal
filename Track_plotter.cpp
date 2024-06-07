#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <numeric>
#include "Analyzer.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TGraphMultiErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TLine.h"
#include "TROOT.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TTree.h"
#include "TBrowser.h"

using namespace std;
namespace fs = std::experimental::filesystem;

vector <string> trim_name( string full_name, char delimiter); 
// void print_histogram(TH1F *histo, string title, string x_axis, string y_axis);
// void print_graph_simple(TGraph *graph, string title, string x_axis, string y_axis);
// void create_and_print_wf_graph_simple (string filename, vector<int> time, shared_ptr<vector<double>> ampl, string tag);
// void create_and_print_wf_graph_lines (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start, double end, double level);
// void print_graph_lines (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1);
// void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
//  double l_xy, double a_xy, double l_z, int d_z, double a_z, double length);

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);

int main(int argc, char**argv) {
    TApplication *myapp=new TApplication("myapp",0,0); cout << "\n" << endl;


    bool batch_mode = 1;
    if (batch_mode == 1) gROOT->SetBatch(1);              // to avoid prints


    string filename_cam = argv[ 1 ];

    string outputfile = argv[ 2 ]; 
    TFile* file_root = new TFile(outputfile.c_str(),"recreate");

    /* ****************************************  Definition of variables and graphs  ******************************************************************  */

    vector<int> BeginScPix;
    vector<int> EndScPix;
    int nSc_red=0;

    // electron drift velocity in LIME, He:CF4 60:40, 800 V/cm 
    const double drift_vel = 5.8; // cm/µs

    //impact point and directionality
    // For alphas (David)
    Int_t NPIP=2000;
    Float_t wFac=3.5;
  
    // For ER (original from Samuele)
    // Int_t NPIP=80;
    // Float_t wFac=2.5;

    double x_impact, y_impact;
    double xbar, ybar;
    int quadrant_cam;
    double angle_cam;
    const double granularity = 0.0155; // cm/pixel 

    /* ****************************************  Opening root recoed file -- CAMERA ******************************************************************  */

    TFile *reco_data_cam = TFile::Open(filename_cam.c_str());
    vector <string> name_trim1 = trim_name(filename_cam, '/');
    cout << "CAM Reco data file openend: " << name_trim1.back() << endl;

    TTree *tree_cam = (TTree*)reco_data_cam->Get("Events");

    int cam_run;        tree_cam->SetBranchAddress("run",   &cam_run);
    int cam_event;      tree_cam->SetBranchAddress("event", &cam_event);

    vector<float> sc_integral;    sc_integral.reserve(150000);    tree_cam->SetBranchAddress("sc_integral",   sc_integral.data());
    vector<float> sc_rms;         sc_rms.reserve(150000);         tree_cam->SetBranchAddress("sc_rms",        sc_rms.data());
    vector<float> sc_tgausssigma; sc_tgausssigma.reserve(150000); tree_cam->SetBranchAddress("sc_tgausssigma",sc_tgausssigma.data());
    vector<float> sc_length;      sc_length.reserve(150000);      tree_cam->SetBranchAddress("sc_length",     sc_length.data());
    vector<float> sc_width;       sc_width.reserve(150000);       tree_cam->SetBranchAddress("sc_width",      sc_width.data());
    vector<float> sc_xmean;       sc_xmean.reserve(150000);       tree_cam->SetBranchAddress("sc_xmean",      sc_xmean.data());  
    vector<float> sc_ymean;       sc_ymean.reserve(150000);       tree_cam->SetBranchAddress("sc_ymean",      sc_ymean.data());  
    vector<float> sc_nhits;       sc_nhits.reserve(150000);       tree_cam->SetBranchAddress("sc_nhits",      sc_nhits.data());  

    UInt_t nSc;                     tree_cam->SetBranchAddress("nSc",   &nSc);
    UInt_t Nredpix=0;               tree_cam->SetBranchAddress("nRedpix",&Nredpix);

    vector<float> sc_redpixID;      sc_redpixID.reserve(150000);    tree_cam->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
    vector<int> XPix;               XPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iy",       XPix.data());  
    vector<int> YPix;               YPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_ix",       YPix.data());
    vector<float> ZPix;             ZPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iz",       ZPix.data());  

    /* **********************************************  Analysis  CAMERA  ********************************************** */

    file_root->cd();

    Int_t nentries = (Int_t)tree_cam->GetEntries();
    for (Int_t cam_entry = 0; cam_entry < nentries; cam_entry++) {
        
        tree_cam->GetEntry(cam_entry);
        // cout << "Cam run: " << cam_run << "; event(pic): " << cam_event << "; nSc: " << nSc << endl;

        file_root->mkdir(Form("Run_%i_ev_%i", cam_run, cam_event));
        file_root->cd(Form("Run_%i_ev_%i", cam_run, cam_event));

        cout << "Cam run: " << cam_run << "; event(pic): " << cam_event << "; nSc: " << nSc << endl;

        sc_redpixID.clear();
        ScIndicesElem(nSc, Nredpix, sc_redpixID.data(), nSc_red, BeginScPix, EndScPix);

        for ( int sc_i = 0; sc_i < nSc_red; sc_i++ ) {

            cout << "\nCluster ID: " << sc_i << endl;

            // if ( (sc_rms[nc] > 6) && (0.152 * sc_tgausssigma[nc] > 0.3) && (0.152 * sc_length[nc] < 80) && (sc_integral[nc] > 1000) && (sc_width[nc] / sc_length[nc] > 0.8)            //cam cut
            // && ( ( (sc_xmean[nc]-1152)*(sc_xmean[nc]-1152) + (sc_ymean[nc]-1152)*(sc_ymean[nc]-1152) ) <(800*800)) ) {
            if ( sc_integral[sc_i]/sc_nhits[sc_i] > 25 && sc_length[sc_i] > 100 && sc_width[sc_i] > 50 ) {   //Alpha cut fropm Giorgio

                //----------- Build Analyser track  -----------//

                Analyzer Track(Form("Track_run_%i_ev_%i", cam_run, cam_event),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[sc_i],EndScPix[sc_i]);

                Track.SetWScal(wFac);
                Track.SetNPIP(NPIP);
                Track.ApplyThr();
                Track.RemoveNoise(50);
                Track.ImpactPoint(Form("Track_event%i_run%i", cam_run, cam_event));
                Track.ScaledTrack(Form("Track_event%i_run%i", cam_run, cam_event));
                Track.Direction();
                Track.ImprCorrectAngle();
                Track.BuildLineDirection();

                //----------- Get important track parameters  -----------//

                angle_cam = Track.GetDir()/TMath::Pi()*180.;
                xbar = Track.GetXbar();
                ybar = Track.GetYbar();
                x_impact = Track.GetXIP() * granularity;
                y_impact = Track.GetYIP() * granularity;

                // The y needs to be inverted for this calculation from the way root makes the plots.
                // TO BE CROSS CHECKED WITH NEW DATA.
                if      ( xbar < 2305./2. && ybar > 2305./2.) quadrant_cam = 1;
                else if ( xbar > 2305./2. && ybar > 2305./2.) quadrant_cam = 2;
                else if ( xbar > 2305./2. && ybar < 2305./2.) quadrant_cam = 3;
                else if ( xbar < 2305./2. && ybar < 2305./2.) quadrant_cam = 4; 

                //----------- Verbose information  -----------//

                cout << "\n\tTrack information: \n" << endl; 
                cout << "--> Position barycenter: " << "x: " << xbar << "; y: " << ybar << endl;
                cout << "--> Quadrant: " << quadrant_cam << endl;
                cout << "--> Angle: " << angle_cam << " degrees." << endl;
                cout << "--> Length (cm): " << sc_length[sc_i] * granularity << endl;


                // Track.SavetoFile("test1");               // Saves TH2D in root file
                // Track.SavetoFileDavid("test5");          // Saves TH2D, colz, in root file
                // Track.SavePic("test2.png");              // saves png of track in folder with a bunch of variables.
                Track.SavePicDir(Form("run_%i_evt_%i.png", cam_run, cam_event));           // Saves 3 plots with the real direction calculation with Samuele's code
                // Track.SaveRootFile("test4.root");        // Saves a couple plots in a different root file.
                // Track.PlotandSavetoFileCOLZ(Form("Track_ev%i_run%i", cam_run, cam_event));            // Saves TH2D, colz, in root file
                // Track.PlotandSavetoFileDirectionalFull("X-Y Analyser");        // Saves the directionality plots all together, in root file

                // Track profiles
                // TCanvas* c_profile = new TCanvas("c_profile","c_profile",1000,500); 
                // c_profile->Divide(1,2);

                // c_profile->cd(1);
                // TH1D* TrackProfileTrans = Track.FillProfile(false);
                // TrackProfileTrans->SetTitle("Transversal Profile");
                // TrackProfileTrans->Draw();

                // c_profile->cd(2); 
                // TH1D* TrackProfileLongi = Track.FillProfile(true);
                // TrackProfileLongi->SetTitle("Longitudinal Profile"); 
                // TrackProfileLongi->Draw();
                // c_profile->Write("Track profiles",TObject::kWriteDelete);
                // c_profile->DrawClone();
                // delete c_profile;
            }   
        }
        sc_redpixID.resize(nSc); // Resize to avoid problem with vector sizes

        // if (cam_entry == 100) break;
    }
    reco_data_cam->Close();
    


    file_root->Close();
    cout << "**Finished**" << endl;
    if (batch_mode == 0) myapp->Run();

    return 0;
}


void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E){
    nSc_red=0;
    B.clear();
    E.clear();

    vector<float> sc_redpix_start;
    sc_redpix_start.push_back(0); // Initialize the start indices

    int parcount=0; // Counter for iterating through pixels

    // Loop through each supercluster to process its reduced pixels
    for(int i=0; i<nSc; i++){
        if(sc_redpixID[i]>0){
            sc_redpix_start.push_back(sc_redpixID[i]); // Add the start index of the reduced pixels
        }
    }

    nSc_red = sc_redpix_start.size(); // Number of processed superclusters

    sc_redpix_start.push_back(npix); // Ensure the last pixel is included

    // Calculate the begin and end indices for each cluster's reduced pixels
    for(int i=0;i<sc_redpix_start.size()-1;i++){
        B.push_back(sc_redpix_start[i]);
        E.push_back(sc_redpix_start[i+1]);
    }

    sc_redpix_start.clear(); // Clear the temporary storage
}

vector <string> trim_name( string full_name, char delimiter) {

    string trimmed_name;
    vector <string> tokens;
    stringstream check1(full_name);
    string intermediate;
    while(getline(check1, intermediate, delimiter)) tokens.push_back(intermediate);

    return tokens;
}
