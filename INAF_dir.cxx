// To compile this program, use:
//   g++ Analyzer.cxx INAF_dir.cxx -O3 -o INAF `root-config --libs --cflags` -lSpectrum
//
// To run the program, use:
//   ./INAF <path_to_rootfile> <output_directory> <NPIP> <wfactor> <threshold> <remove_noise_value>
//
// Example:
//   ./INAF_dir data.root ./output 10 1.0 20 30

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <chrono>
#include <cmath>  // for M_PI
#include <fstream>
#include <filesystem>

using namespace std;

// Function to calculate indices for clusters based on reduced pixels
// nSc: number of superclusters
// npix: total number of pixels
// sc_redpixID: array of indices for the first reduced pixel in each cluster
// nSc_red: number of processed superclusters (output)
// B: begin indices for reduced pixels (output)
// E: end indices for reduced pixels (output)
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

int main(int argc, char** argv){
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " <path_to_rootfile> <output_directory> <NPIP> <wfactor> <threshold> <remove_noise_value>" << std::endl;
    return 1;
  }

  int NPIP = std::stoi(argv[3]);
  double wfactor = std::stod(argv[4]);
  int threshold = std::stoi(argv[5]);
  int remove_noise_value = std::stoi(argv[6]);

  // Open the input ROOT file and get the "Events" TTree
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");

  // Variable declarations for event data
  int nmax=2000000000;
  int nscmax=5000;
  int npixel=2304;
  int npixelsmall=250;
  float slimnesslimit=0.6;

  unsigned int nSc=0;
  int nSc_red=0;
  UInt_t Nredpix=0;
  Int_t sc_npix=0;
  int run;
  int event, event_out;

  // Reserve space for vectors to store pixel and cluster data
  vector<float> sc_redpixID;
  sc_redpixID.reserve(nmax);
  vector<UInt_t> ScNpixels;
  ScNpixels.reserve(nscmax);
  vector<int> XPix;
  XPix.reserve(nmax);
  vector<int> YPix;
  YPix.reserve(nmax);
  vector<float> ZPix;
  ZPix.reserve(nmax);
  // Additional vectors for analysis
  vector<float> xmean;
  xmean.reserve(nscmax);
  vector<float> ymean;
  ymean.reserve(nscmax);
  vector<float> ymin;
  ymin.reserve(nscmax);
  vector<float> ymax;
  ymax.reserve(nscmax);
  vector<float> xmin;
  xmin.reserve(nscmax);
  vector<float> xmax;
  xmax.reserve(nscmax);
  vector<float> scsize;
  scsize.reserve(nscmax);
  vector<float> scnhits;
  scnhits.reserve(nscmax);
  vector<float> v_sc_rms;
  v_sc_rms.reserve(nscmax);
  vector<float> v_sc_tgausssigma;
  v_sc_tgausssigma.reserve(nscmax);
  vector<float> v_sc_theta;
  v_sc_theta.reserve(nscmax);
  vector<float> width;
  width.reserve(nscmax);
  vector<float> length;
  length.reserve(nscmax);
  vector<float> integral;
  integral.reserve(nscmax);

  // Prepare for the analysis by setting branch addresses for tree
  // This connects the variables defined above with the corresponding data in the ROOT file
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
  tree->SetBranchAddress("nRedpix",&Nredpix);
  
  tree->SetBranchAddress("redpix_ix",XPix.data());
  tree->SetBranchAddress("redpix_iy",YPix.data());

  tree->SetBranchAddress("redpix_iz",ZPix.data());
  tree->SetBranchAddress("sc_width",width.data());
  tree->SetBranchAddress("sc_length",length.data());
  tree->SetBranchAddress("sc_integral",integral.data());
  tree->SetBranchAddress("sc_xmean",xmean.data());
  tree->SetBranchAddress("sc_ymean",ymean.data());
  tree->SetBranchAddress("sc_xmin",xmin.data());
  tree->SetBranchAddress("sc_ymin",ymin.data());
  tree->SetBranchAddress("sc_xmax",xmax.data());
  tree->SetBranchAddress("sc_ymax",ymax.data());
  tree->SetBranchAddress("sc_size",scsize.data());
  tree->SetBranchAddress("sc_nhits",scnhits.data());
  tree->SetBranchAddress("sc_theta",v_sc_theta.data());
  tree->SetBranchAddress("sc_rms",v_sc_rms.data());
  tree->SetBranchAddress("sc_tgausssigma",v_sc_tgausssigma.data());

  // Prepare vectors for storing analysis results
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;

  // Variables for storing intermediate analysis results
  double phi_DIR=0, phi_DIR_deg=0, phi_PCA=0, phi=0;
  double phi_HT_skew=0, phi_HT_maxpx=0, phi_HT_maxpxRebin=0, phi_HT_peakprof=0, phi_HT_maxprof=0, phi_HT_integral=0;
  double scint = 0;
  double skewness_L=0, skewness_T=0, kurtosis_L=0, kurtosis_T=0;
  double xbar=0, ybar=0;
  double xl=0, yl=0, xr=0, yr=0;
  double reco_theta=0, x_mean=0, y_mean=0, x_min=0, x_max=0, y_min=0, y_max=0;
  int procTime=0;
  double profileMean=0, profileSkew=0; 
  bool puCandidate=false;

  //impact point and directionality
  double xIP=0, yIP=0;

  vector<std::pair<double,double>> peakslong;
  vector<std::pair<double,double>> peakstrans;
  int npeaks=0;
  double peak_density=0, tracklength=0, track_width=0, recolength=0, recowidth=0, reco_sc_rms=0, reco_sc_tgausssigma=0;
  vector<double> peak_distance, peak_width;

  //ANALYSIS
  int counter=0;
  int ired=0;

  // Prepare output file path and open output file
  // Assuming 'argv[1]' is the input file and 'argv[2]' is the output directory
  string inputFilePath = argv[1];
  string outputDir = argc > 2 ? argv[2] : "."; // Use current directory if not provided
  string filename = inputFilePath.substr(inputFilePath.find_last_of("/\\") + 1);
  string outputFilePath = outputDir + "/Analysis_" + filename;

  cout<<filename<<endl;
  TFile* fout = new TFile(outputFilePath.c_str(), "recreate");
  cout<<Form("%s/Analysis_%s",argv[2],filename.c_str())<<endl;
  fout->cd();
  fout->mkdir("Tracks");

  TTree* out_tree = new TTree("track_info","track_info");
  out_tree->Branch("run",&run);
  out_tree->Branch("event",&event_out);
  out_tree->Branch("nSc",&nSc);
  out_tree->Branch("nSc_red",&nSc_red);
  out_tree->Branch("Integral",&scint);
  out_tree->Branch("ScSize",&sc_npix);
  out_tree->Branch("AngleDIR",&phi_DIR);
  out_tree->Branch("DegreeDIR",&phi_DIR_deg);
  out_tree->Branch("AnglePCA",&phi_PCA);
  out_tree->Branch("RecoTheta",&reco_theta);
  out_tree->Branch("RecoScRMS",&reco_sc_rms);
  out_tree->Branch("RecoScTGaussSigma",&reco_sc_tgausssigma);
  out_tree->Branch("X_ImpactPoint",&xIP);
  out_tree->Branch("Y_ImpactPoint",&yIP);
  out_tree->Branch("Xmean",&x_mean);
  out_tree->Branch("Ymean",&y_mean);
  out_tree->Branch("Xmin",&x_min);
  out_tree->Branch("Ymin",&y_min);
  out_tree->Branch("Xmax",&x_max);
  out_tree->Branch("Ymax",&y_max);
  out_tree->Branch("XBar",&xbar);
  out_tree->Branch("YBar",&ybar);
  out_tree->Branch("procTime",&procTime);
  out_tree->Branch("ProfileMean",&profileMean);
  out_tree->Branch("ProfileSkew",&profileSkew);
  out_tree->Branch("PileUpCandidate",&puCandidate);

  int pileUpCounter=0;
  cout<<"this run has "<<tree->GetEntries()<<" entries"<<endl;
  for(int k=0;k<tree->GetEntries();k++) //only for FUSION He/CF4 INAF Nov24
  //for(int k=100000;k<110000;k++) //only for FUSION He/CF4 INAF Nov24
  {
    sc_redpixID.clear();
    tree->GetEntry(k);
    //for reduced pixels:
    ScIndicesElem(nSc,Nredpix,sc_redpixID.data(),nSc_red,BeginScPix,EndScPix);

    //Start the cycle on the supercluster of the event
    int pixcounter =0;
    for(int i=0;i<nSc_red;i++)
    {
      scint = integral[i];
      recolength=length[i];
      recowidth=width[i];
      sc_npix = scnhits[i];
      reco_theta=v_sc_theta[i];
      x_mean=xmean[i];
      y_mean=ymean[i];
      x_min=xmin[i];
      x_max=xmax[i];
      y_min=ymin[i];
      y_max=ymax[i];
      reco_sc_rms=v_sc_rms[i];
      reco_sc_tgausssigma=v_sc_tgausssigma[i];

      pixcounter += EndScPix[i] - BeginScPix[i];
      //! Condition to filter out certain events based on physical properties
      //For Polarized 8Kev photon in MANGO
      //if(scint>25000 && scint<50000 && recowidth/recolength>0.7 && recowidth/recolength<1 && x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 && run>22700)
      //if(scint>25000 && scint<50000 && recowidth/recolength>0.7 && recowidth/recolength<1 && x_mean>950 && x_mean<1050 && y_mean<950 && y_mean>1050 && run>22700)
      //For Polarized 17Kev photon in MANGO FUSION He/CF4
      //if( x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 && scint<90000 && scint>60000 && sc_npix<6000  && run>22700)
      //For Polarized 17Kev photon in MANGO QUEST EHD Ar/CF4 0deg
      if( x_mean>1300 && x_mean<2500 && y_mean>400 && y_mean<1500 && scint<45000 && scint>25000)
      //For Polarized 17Kev photon in MANGO QUEST EHD Ar/CF4 90deg
      //if( x_mean>1800 && x_mean<2000 && y_mean<1125 && y_mean>925 && scint<38000 && scint>24000)
      // For LIME 55Fe
      //if (y_max < 1250 && y_min > 1050 && x_max < 1250 && x_min > 1050 && scint>2000 && reco_sc_rms>5 && reco_sc_tgausssigma>2.63 && reco_sc_tgausssigma<4.5 && recowidth/recolength>0.6 )
      {
        event_out = event; // Store the event number for output
        // Start timers for the entire track processing
        auto t0 = std::chrono::steady_clock::now();
        //! Analyzer constructor 
        Analyzer Traccia(Form("Track%i_event%i_run%i_entry%i", counter, event_out, run,k),
            XPix.data(), YPix.data(), ZPix.data(),
            BeginScPix[i], EndScPix[i],false); // true for rebinning

        if (Traccia.Getbinmax()<=0) continue; //skip if the track is empty

        auto t1 = std::chrono::steady_clock::now();
        //! Set directionality parameters
        // (We skip timing SetWScal and SetNPIP now, or just ignore them in the printout)
        Traccia.SetWScal(wfactor);
        Traccia.SetNPIP(NPIP);
        //! ApplyThr
        //Traccia.ApplyThr(10);// He/CF4 fusion
        Traccia.ApplyThr(threshold);
        auto t4 = std::chrono::steady_clock::now();
        //! RemoveNoise
        //Traccia.RemoveNoise(30);// He/CF4 fusion
        Traccia.RemoveNoise(remove_noise_value);
        auto t5 = std::chrono::steady_clock::now();
        //! Check for pile-up (i.e., Traccia.PileUpCandidate())
        auto t5_pileup_start = std::chrono::steady_clock::now();
        puCandidate = Traccia.PileUpCandidate(false, counter, true, 0., 0.); //for 17keV
        //puCandidate = Traccia.PileUpCandidate(false, counter, false, 0.,0.); //for 8keV
        if (puCandidate) {
          pileUpCounter++;
          //Traccia.TrackProfilePlotSave(Form("Track%i_event%i_run%i_entry%i.png", counter, event, run,k));
          continue;
        }
        auto t5_pileup_end   = std::chrono::steady_clock::now();
        //! Get track profile statistics
        Traccia.GetTrackProfileStats(profileMean, profileSkew);
        //! ImpactPoint
        Traccia.ImpactPoint(Form("TrackIPRegion%i_run%i_evt%i", k, run, counter));
        auto t6 = std::chrono::steady_clock::now();
        //! ScaledTrack
        Traccia.ScaledTrack(Form("TrackScaled%i_run%i_evt%i", k, run, counter));
        auto t7 = std::chrono::steady_clock::now();
        //! Direction
        Traccia.Direction();
        auto t8 = std::chrono::steady_clock::now();
        //! ImprCorrectAngle (skipped in timing printout)
        Traccia.ImprCorrectAngle();
        //! BuildLineDirection
        Traccia.BuildLineDirection();
        auto t10 = std::chrono::steady_clock::now();
        // Grab results
        xIP         = Traccia.GetXIP();
        yIP         = Traccia.GetYIP();
        phi_DIR     = Traccia.GetDir();
        phi_DIR_deg = phi_DIR * (180.0 / M_PI);
        phi_PCA     = Traccia.AngleLineMaxRMS();
        xbar        = Traccia.GetXbar();
        ybar        = Traccia.GetYbar();
        // End timing for the entire block
        auto tEnd = std::chrono::steady_clock::now();

        //! Print only every N events
        if (k % 10000000 == 0)
        {
            // Save a diagnostic image
            Traccia.SavePicDir(Form("Track%i_event%i_run%i_entry%i.png", counter, event_out, run,k));
            // Print the results to the console
            std::cout << "Processing entry " << k << " / " << tree->GetEntries() << std::endl;
            std::cout << "counter: " << counter << std::endl;
            std::cout << "XIP: " << xIP << "  YIP: " << yIP << std::endl;
            std::cout << "XIPPrev: " << Traccia.GetXIPPrev() << "  YIPPrev: " << Traccia.GetYIPPrev() << std::endl;
            std::cout << "Degree: " << phi_DIR_deg << "  tan(angle): " << std::tan(phi_DIR) << std::endl;
            // Compute durations in microseconds
            long long dt_ctor = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
            long long dt_applyThr = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t1).count();
            long long dt_removeNoise = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
            long long dt_pileUpCandidate = std::chrono::duration_cast<std::chrono::microseconds>(t5_pileup_end - t5_pileup_start).count();
            long long dt_impactPoint = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5_pileup_end).count();
            long long dt_scaledTrack = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
            long long dt_direction   = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
            long long dt_buildLineDirection = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t8).count();
            long long dt_total = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
            // Print each as "us" and as "% of total"
            auto printTime = [&](const char* label, long long dt) {
                double pct = 100.0 * (double)dt / (double)dt_total;
                std::cout << label << dt << " us (" << pct << " %)\n";
            };
            printTime("Time Analyzer ctor:         ", dt_ctor);
            printTime("Time ApplyThr:              ", dt_applyThr);
            printTime("Time RemoveNoise:           ", dt_removeNoise);
            printTime("Time PileUpCandidate:       ", dt_pileUpCandidate);
            printTime("Time ImpactPoint:           ", dt_impactPoint);
            printTime("Time ScaledTrack:           ", dt_scaledTrack);
            printTime("Time Direction:             ", dt_direction);
            printTime("Time BuildLineDirection:    ", dt_buildLineDirection);
            // Lastly, the total is always 100%
            std::cout << "TOTAL time for all steps: " << dt_total << " us (100%)\n\n";
          }
        // Fill the tree with the total time for this event
        procTime = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
        if (procTime < 0) {
          std::cerr << "Warning: Negative processing time for event " << k << ", resetting to 0." << std::endl;
          procTime = 0; // Reset to zero if negative
          Traccia.SavePicDir(Form("Track%i_event%i_run%i_entry%i.png", counter, event, run,k));
        }
        out_tree->Fill();
      }
    counter++;
    }//!superclusters
  sc_redpixID.resize(nSc);//Madonna fai il resize altriemnti lui si ricorda la dimensione precedente
  }//!ttree entries
  out_tree->Write();
  fout->Close();
  cout<<"pile up percentage: "<<(double)pileUpCounter/(double)counter<<endl;
  return 0;
}
