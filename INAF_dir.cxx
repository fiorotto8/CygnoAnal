//This in particular compile using  g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum
//Then use as ./nameprog path_to_rootfile

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


using namespace std;

/*old
void ScIndicesElem(int nSc,int* ScNelements, vector<int>& B, vector<int>& E){
  B.clear();
  E.clear();

  int parcount=0;

  for(int i=0;i<nSc;i++){
    B.push_back(parcount);
    E.push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}
*/

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

// Function to find the divisor of A closest to B
int closest_divisor(int A, float B){
  vector<float>  div_dist;
  vector<int> div_candidate;

  // Find all divisors of A and their distance to B
  for(int i=1; i<=A; i++){
    if(A%i==0){ // Check if i is a divisor of A
      div_dist.push_back(abs((float)i-B));
      div_candidate.push_back(i);
    }
  }

  // Find the divisor closest to B
  auto result = min_element(div_dist.begin(), div_dist.end());
  return div_candidate[result-div_dist.begin()];
}

int main(int argc, char** argv){

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
  int event;

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
  //tree->SetBranchAddress("sc_size",&sc_npix);
  tree->SetBranchAddress("nRedpix",&Nredpix);
  tree->SetBranchAddress("redpix_ix",YPix.data());
  tree->SetBranchAddress("redpix_iy",XPix.data());
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

  //impact point and directionality
  //Int_t NPIP=300;
  //Float_t wFac=2.;
  double xIP=0, yIP=0;

  vector<std::pair<double,double>> peakslong;
  vector<std::pair<double,double>> peakstrans;
  int npeaks=0;
  double peak_density=0, tracklength=0, track_width=0, recolength=0, recowidth=0, reco_sc_rms=0, reco_sc_tgausssigma=0;
  vector<double> peak_distance, peak_width;

  //ANALYSIS
  int counter=0;
  int ired=0;
  //int counterall=0;

  // Prepare output file path and open output file
  // Assuming 'argv[1]' is the input file and 'argv[2]' is the output directory
  string inputFilePath = argv[1];
  string outputDir = argc > 2 ? argv[2] : "."; // Use current directory if not provided
  string filename = inputFilePath.substr(inputFilePath.find_last_of("/\\") + 1);
  string outputFilePath = outputDir + "/Analysis_" + filename;

  //string filename(argv[1]);
  //filename = filename.substr(filename.find_last_of("/\\")+1);
  cout<<filename<<endl;
  TFile* fout = new TFile(outputFilePath.c_str(), "recreate");
  cout<<Form("%s/Analysis_%s",argv[2],filename.c_str())<<endl;
  fout->cd();
  fout->mkdir("Tracks");

  TTree* out_tree = new TTree("track_info","track_info");
  out_tree->Branch("run",&run);
  out_tree->Branch("event",&event);
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
  /*
  out_tree->Branch("Npeaks",&npeaks);
  out_tree->Branch("PeakDensity",&peak_density);
  out_tree->Branch("TrackLength",&tracklength);
  out_tree->Branch("TrackWidth",&track_width);
  out_tree->Branch("RecoLength",&recolength);
  out_tree->Branch("RecoWidth",&recowidth);
  out_tree->Branch("PeakDistance",&peak_distance);
  out_tree->Branch("Skewness",&skewness_L);
  out_tree->Branch("Kurtosis",&kurtosis_T);
  out_tree->Branch("phi_HT_skew",&phi_HT_skew);
  out_tree->Branch("phi_HT_maxpx",&phi_HT_maxpx);
  out_tree->Branch("phi_HT_maxpxRebin",&phi_HT_maxpxRebin);
  out_tree->Branch("phi_HT_peakprof",&phi_HT_peakprof);
  out_tree->Branch("phi_HT_maxprof",&phi_HT_maxprof);
  out_tree->Branch("phi_HT_integral",&phi_HT_integral);
  */

  cout<<"this run has "<<tree->GetEntries()<<" entries"<<endl;
  for(int k=0;k<tree->GetEntries();k++)
  //for(int k=0;k<10000;k++)
  //for(int k=0;k<1;k++)
  {
    //cout<<"Entry "<<k<<endl;
    sc_redpixID.clear();
    tree->GetEntry(k);
    //cout << "Nev: "<< k << "\nnSc:  " << nSc << " event "<< event <<endl;
    //for reduced pixels:
    ScIndicesElem(nSc,Nredpix,sc_redpixID.data(),nSc_red,BeginScPix,EndScPix);

    //cout<<"nSc "<<nSc<<" nSc_red "<<nSc_red<<" Nredpix "<<Nredpix<<endl;

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
      //cout<<"SC "<<i<<endl;

      pixcounter += EndScPix[i] - BeginScPix[i];
      //cout<<"counted pix: "<<pixcounter<<endl;

      // Condition to filter out certain events based on physical properties
      //For Polarized 8Kev photon in MANGO
      //if(scint>25000 && scint<45000 && recowidth/recolength>0.7 && recowidth/recolength<1 && x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 )
      //For Polarized 8Kev photon in MANGO
      if( x_mean>900 && x_mean<1350 && y_mean<1350 && y_mean>900 && scint<82000 && scint>61000 && sc_npix<6000 )
      // For LIME 55Fe
      //if (y_max < 1250 && y_min > 1050 && x_max < 1250 && x_min > 1050 && scint>2000 && reco_sc_rms>5 && reco_sc_tgausssigma>2.63 && reco_sc_tgausssigma<4.5 && recowidth/recolength>0.6 )
      {
        // Start timers for each step
        auto t0 = std::chrono::high_resolution_clock::now();

        Analyzer Traccia(Form("Track%i_event%i_run%i", counter, k, run),
                        XPix.data(), YPix.data(), ZPix.data(),
                        BeginScPix[i], EndScPix[i]);
        auto t1 = std::chrono::high_resolution_clock::now();

        // (We skip timing SetWScal and SetNPIP now, or just ignore them in the printout)
        Traccia.SetWScal(2.);
        Traccia.SetNPIP(300);

        // ApplyThr
        Traccia.ApplyThr(10);
        auto t4 = std::chrono::high_resolution_clock::now();

        // RemoveNoise
        Traccia.RemoveNoise(30);
        auto t5 = std::chrono::high_resolution_clock::now();

        // ImpactPoint
        Traccia.ImpactPoint(Form("TrackIPRegion%i_run%i_evt%i", k, run, counter));
        auto t6 = std::chrono::high_resolution_clock::now();

        // ScaledTrack
        Traccia.ScaledTrack(Form("TrackScaled%i_run%i_evt%i", k, run, counter));
        auto t7 = std::chrono::high_resolution_clock::now();

        // Direction
        Traccia.Direction();
        auto t8 = std::chrono::high_resolution_clock::now();

        // ImprCorrectAngle (we skip this in printout, so no new timer needed)
        Traccia.ImprCorrectAngle();

        // BuildLineDirection
        Traccia.BuildLineDirection();
        auto t10 = std::chrono::high_resolution_clock::now();

        // Grab results
        xIP         = Traccia.GetXIP();
        yIP         = Traccia.GetYIP();
        phi_DIR     = Traccia.GetDir();
        phi_DIR_deg = phi_DIR * (180.0 / M_PI);
        phi_PCA     = Traccia.AngleLineMaxRMS();
        xbar        = Traccia.GetXbar();
        ybar        = Traccia.GetYbar();

        // End timing for the entire block
        auto tEnd = std::chrono::high_resolution_clock::now();

        // Print times only every 10k events
        if (k % 10000 == 0)
        {
            // Save a diagnostic image
            Traccia.SavePicDir(Form("Track%i_event%i_run%i.png", counter, event, run));

            std::cout << "Processing entry " 
                      << k << " / " << tree->GetEntries() << std::endl;


            cout<<"counter: "<<counter<<endl;
            cout<<"XIP: "<<xIP<<" YIP: "<<yIP<<endl;
            cout<<"XIPPrev: "<<Traccia.GetXIPPrev()<<" YIPPrev: "<<Traccia.GetYIPPrev()<<endl;
            cout<<" Degree: "<<phi_DIR_deg<<" tan angle: "<<tan(phi_DIR)<<endl;


            // ---------------------------------------------------------------------
            // (1) Compute dudrations in microseconds
            // ---------------------------------------------------------------------
            auto dt_ctor           = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
            auto dt_applyThr       = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t1).count(); 
            //  ^ note: "t3" was SetNPIP, which we don't print, so you might prefer (t4 - t0) minus something, 
            //    or specifically measure (t4 - <endOfSetNPIP>), etc. 
            //    But for clarity, let's keep as-is and accept that part includes SetNPIP if there's no t3 timing.
            //    If you'd like a true "ApplyThr only" measurement, put a time point *right before* ApplyThr.
            
            auto dt_removeNoise    = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
            auto dt_impactPoint    = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
            auto dt_scaledTrack    = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
            auto dt_direction      = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
            auto dt_buildDirection = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t8).count(); 
            //  ^ "t9" is ImprCorrectAngle's start, which we aren't showing, 
            //    so you might want to measure (t10 - t8) to *combine* Direction + ImprCorrectAngle + BuildLineDirection, 
            //    or add an extra time-point right after ImprCorrectAngle if you want a pure BuildLineDirection measure.

            auto dt_total = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();

            // ---------------------------------------------------------------------
            // (2) Print each as "us" and as "% of total"
            // ---------------------------------------------------------------------
            // We'll define a helper lambda to keep code clean:
            auto printTime = [&](const char* label, long long dt) {
                double pct = 100.0 * (double)dt / (double)dt_total;
                std::cout << label << dt << " us (" << pct << " %)\n";
            };

            printTime("Time Analyzer ctor:       ", dt_ctor);
            printTime("Time ApplyThr:            ", dt_applyThr);
            printTime("Time RemoveNoise:         ", dt_removeNoise);
            printTime("Time ImpactPoint:         ", dt_impactPoint);
            printTime("Time ScaledTrack:         ", dt_scaledTrack);
            printTime("Time Direction:           ", dt_direction);
            printTime("Time BuildLineDirection:  ", dt_buildDirection);

            // Lastly, the total is always 100%
            std::cout << "TOTAL time for all steps: "
                      << dt_total << " us (100%)\n" << std::endl;
        }

        procTime=std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
        out_tree->Fill();
      }
    counter++;
    }//superclusters

  sc_redpixID.resize(nSc);//Madonna fai il resize altriemnti lui si ricorda la dimensione precedente

  }//ttree entries
  out_tree->Write();
  fout->Close();
  return 0;
}
