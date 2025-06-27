//g++ Analyzer.cxx INAF_dir_MP.cxx -O3 -o INAF_MP `root-config --libs --cflags` -lSpectrum#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <thread>
#include <filesystem>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <chrono>
#include <cmath>
#include <TROOT.h>
#include <TFileMerger.h>
#include "Analyzer.h"
// Include ROOT's TProcessExecutor for multiprocessing
#include <ROOT/TProcessExecutor.hxx>

using namespace std;

// Function to calculate indices for clusters based on reduced pixels (unchanged)
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red,vector<int>& B, vector<int>& E) {
    nSc_red = 0;
    B.clear();
    E.clear();
    vector<float> sc_redpix_start;
    sc_redpix_start.push_back(0);             // start index of first cluster
    // Determine start index of reduced pixels for each supercluster
    for (int i = 0; i < nSc; ++i) {
        if (sc_redpixID[i] > 0) {
            sc_redpix_start.push_back(sc_redpixID[i]);
        }
    }
    nSc_red = sc_redpix_start.size();
    sc_redpix_start.push_back(npix);          // ensure last pixel index is included
    // Compute begin and end indices for each cluster's pixels
    for (size_t i = 0; i < sc_redpix_start.size() - 1; ++i) {
        B.push_back((int)sc_redpix_start[i]);
        E.push_back((int)sc_redpix_start[i+1]);
    }
    sc_redpix_start.clear();
}

/**
 * @brief Checks if a given value passes minimum and maximum cut criteria specified in a configuration map.
 *
 * This function looks for keys named "<name>_min" and "<name>_max" in the provided configuration map.
 * If these keys exist, it compares the given value against the corresponding minimum and maximum thresholds.
 * The value passes the cut if it is not less than the minimum (if specified) and not greater than the maximum (if specified).
 *
 * @param config A map containing configuration parameters as string key-value pairs.
 * @param name The base name of the parameter to check (used to form "<name>_min" and "<name>_max" keys).
 * @param value The value to be tested against the cut criteria.
 * @return true if the value passes the cut; false otherwise.
 */
bool pass_cut(const std::map<std::string, std::string>& config, const std::string& name, double value) {
    auto it_min = config.find(name + "_min");
    auto it_max = config.find(name + "_max");
    if (it_min != config.end() && value < std::stod(it_min->second)) return false;
    if (it_max != config.end() && value > std::stod(it_max->second)) return false;
    return true;
}

/**
 * @brief Parses a configuration file into a map of key-value pairs.
 *
 * This function reads the specified configuration file line by line,
 * ignoring comments (lines starting with '#') and empty lines.
 * Each valid line should contain a key-value pair separated by an '=' character.
 * Leading and trailing whitespace around keys and values is trimmed.
 *
 * @param filename The path to the configuration file to parse.
 * @return std::map<std::string, std::string> A map containing the parsed key-value pairs.
 */
std::map<std::string, std::string> parse_config(const std::string& filename) {
    std::ifstream fin(filename);
    std::map<std::string, std::string> config;
    std::string line;
    while (std::getline(fin, line)) {
        // Remove whitespace and comments
        size_t pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);
        if (line.empty()) continue;
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq+1);
        // Trim
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        val.erase(0, val.find_first_not_of(" \t"));
        val.erase(val.find_last_not_of(" \t") + 1);
        config[key] = val;
    }
    return config;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.txt>\n";
        return 1;
    }
    auto config = parse_config(argv[1]);    

    Long64_t printout = std::stoll(config["print_output"]);
    string inputFilePath = config["input_file"];
    string outputDir = config["output_dir"];
    if (outputDir.empty()) outputDir = ".";  // default to current directory
    int NPIP = std::stoi(config["NPIP"]);
    double wfactor = std::stod(config["wfactor"]);
    int threshold = std::stoi(config["threshold"]);
    int remove_noise_value = std::stoi(config["remove_noise_value"]);
    unsigned int nWorkers = config.count("n_workers") ? std::stoi(config["n_workers"]) : std::thread::hardware_concurrency();
    
    std::filesystem::path inputPath(inputFilePath);
    std::string inputFileNameOnly = inputPath.filename().string();
    std::string finalOutputPath = outputDir + "/AnalysisMP_" + inputFileNameOnly;

    // Enable ROOT thread-safety (good practice for concurrent file access):contentReference[oaicite:3]{index=3}
    ROOT::EnableThreadSafety();

    // Open input file and retrieve TTree to determine total events
    std::cout << "Trying to open file: [" << inputFilePath << "]" << std::endl;
    TFile *fin = TFile::Open(inputFilePath.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file " << inputFilePath << "\n";
        return 1;
    }
    TTree *tree = (TTree*)fin->Get("Events");
    if (!tree) {
        std::cerr << "Error: no TTree named 'Events' in input file\n";
        fin->Close();
        return 1;
    }
    Long64_t totalEntries = tree->GetEntries();
    fin->Close();

    std::cout << "This run has " << totalEntries << " entries" << std::endl;
    if (totalEntries == 0) {
        // No events: create an output file with an empty track_info tree (same branches)
        std::filesystem::path inputPath(inputFilePath);
        std::string inputFileNameOnly = inputPath.filename().string();
        std::string finalPath = outputDir + "/AnalysisMP_" + inputFileNameOnly;
        TFile *fout = TFile::Open(finalPath.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            std::cerr << "Error: cannot create output file " << finalPath << "\n";
            return 1;
        }
        fout->mkdir("Tracks");  // preserve directory structure
        TTree *out_tree = new TTree("track_info", "track_info");
        // Define variables for branches (types match original)
        unsigned int run = 0;
        int event = 0;
        unsigned int nSc = 0;
        int nSc_red = 0;
        double scint = 0;
        Int_t sc_npix = 0;
        double phi_DIR = 0, phi_DIR_deg = 0, phi_PCA = 0;
        double reco_theta = 0, reco_sc_rms = 0, reco_sc_tgausssigma = 0;
        double xIP = 0, yIP = 0;
        double x_mean = 0, y_mean = 0;
        double x_min = 0, x_max = 0, y_min = 0, y_max = 0;
        double xbar = 0, ybar = 0;
        Int_t procTime = 0;
        double profileMean = 0, profileSkew = 0;
        bool puCandidate = false;
        // Create tree branches (same names as original output)
        out_tree->Branch("run", &run);
        out_tree->Branch("event", &event);
        out_tree->Branch("nSc", &nSc);
        out_tree->Branch("nSc_red", &nSc_red);
        out_tree->Branch("Integral", &scint);
        out_tree->Branch("ScSize", &sc_npix);
        out_tree->Branch("AngleDIR", &phi_DIR);
        out_tree->Branch("DegreeDIR", &phi_DIR_deg);
        out_tree->Branch("AnglePCA", &phi_PCA);
        out_tree->Branch("RecoTheta", &reco_theta);
        out_tree->Branch("RecoScRMS", &reco_sc_rms);
        out_tree->Branch("RecoScTGaussSigma", &reco_sc_tgausssigma);
        out_tree->Branch("X_ImpactPoint", &xIP);
        out_tree->Branch("Y_ImpactPoint", &yIP);
        out_tree->Branch("Xmean", &x_mean);
        out_tree->Branch("Ymean", &y_mean);
        out_tree->Branch("Xmin", &x_min);
        out_tree->Branch("Ymin", &y_min);
        out_tree->Branch("Xmax", &x_max);
        out_tree->Branch("Ymax", &y_max);
        out_tree->Branch("XBar", &xbar);
        out_tree->Branch("YBar", &ybar);
        out_tree->Branch("procTime", &procTime);
        out_tree->Branch("ProfileMean", &profileMean);
        out_tree->Branch("ProfileSkew", &profileSkew);
        out_tree->Branch("PileUpCandidate", &puCandidate);
        // Write empty tree and clean up
        out_tree->Write();
        fout->Close();
        std::cout << "pile up percentage: 0" << std::endl;
        return 0;
    }

    // Determine number of worker processes (default to number of CPU cores)
    if (nWorkers == 0) nWorkers = 2;  // use 2 workers if concurrency not detectable
    if (totalEntries < (Long64_t)nWorkers) {
        nWorkers = (unsigned int) totalEntries;
        if (nWorkers < 1) nWorkers = 1;
    }
    // For debugging, force single worker (remove after debugging)
    // nWorkers = 1;

    // Create list of worker IDs for mapping tasks
    std::vector<int> workerIds;
    workerIds.reserve(nWorkers);
    for (unsigned int i = 0; i < nWorkers; ++i) {
        workerIds.push_back(i);
    }

    // Set up parallel processing with TProcessExecutor:contentReference[oaicite:4]{index=4}
    ROOT::TProcessExecutor executor(nWorkers);
    auto results = executor.Map([=](int workerId) {
        // Each worker opens the input file and sets up its own tree and data buffers
        TFile *fIn = TFile::Open(inputFilePath.c_str(), "READ");
        if (!fIn || fIn->IsZombie()) {
            std::cerr << "Worker " << workerId << ": cannot open input file\n";
            if (fIn) fIn->Close();
            return std::make_pair(-1LL, -1LL);
        }
        TTree *tree = (TTree*)fIn->Get("Events");
        if (!tree) {
            std::cerr << "Worker " << workerId << ": 'Events' tree not found in file\n";
            fIn->Close();
            return std::make_pair(-1LL, -1LL);
        }

        // Allocate and reserve memory for event data (reduced buffer sizes)
        int nmax = 1000000;   // 1 million, adjust as needed for your data
        int nscmax = 10000;   // 10 thousand, adjust as needed for your data
        unsigned int nSc = 0;
        int nSc_red = 0;
        UInt_t Nredpix = 0;
        Int_t sc_npix = 0;
        int run = 0;
        int event = 0;
        std::vector<float> sc_redpixID; sc_redpixID.resize(nmax);
        std::vector<UInt_t> ScNpixels;   ScNpixels.resize(nscmax);
        std::vector<int>    XPix;        XPix.resize(nmax);
        std::vector<int>    YPix;        YPix.resize(nmax);
        std::vector<float>  ZPix;        ZPix.resize(nmax);
        std::vector<float>  xmean;       xmean.resize(nscmax);
        std::vector<float>  ymean;       ymean.resize(nscmax);
        std::vector<float>  ymin;        ymin.resize(nscmax);
        std::vector<float>  ymax;        ymax.resize(nscmax);
        std::vector<float>  xmin;        xmin.resize(nscmax);
        std::vector<float>  xmax;        xmax.resize(nscmax);
        std::vector<float>  scsize;      scsize.resize(nscmax);
        std::vector<float>  scnhits;     scnhits.resize(nscmax);
        std::vector<float>  v_sc_rms;    v_sc_rms.resize(nscmax);
        std::vector<float>  v_sc_tgausssigma; v_sc_tgausssigma.resize(nscmax);
        std::vector<float>  v_sc_theta;  v_sc_theta.resize(nscmax);
        std::vector<float>  width;       width.resize(nscmax);
        std::vector<float>  length;      length.resize(nscmax);
        std::vector<float>  integral;    integral.resize(nscmax);
        // Set branch addresses for input tree
        tree->SetBranchAddress("run",         &run);
        tree->SetBranchAddress("event",       &event);
        tree->SetBranchAddress("nSc",         &nSc);
        tree->SetBranchAddress("sc_redpixIdx", sc_redpixID.data());
        tree->SetBranchAddress("nRedpix",     &Nredpix);
        tree->SetBranchAddress("redpix_ix",   XPix.data());   // note: X/Y swapped as in original
        tree->SetBranchAddress("redpix_iy",   YPix.data());
        tree->SetBranchAddress("redpix_iz",   ZPix.data());
        tree->SetBranchAddress("sc_width",    width.data());
        tree->SetBranchAddress("sc_length",   length.data());
        tree->SetBranchAddress("sc_integral", integral.data());
        tree->SetBranchAddress("sc_xmean",    xmean.data());
        tree->SetBranchAddress("sc_ymean",    ymean.data());
        tree->SetBranchAddress("sc_xmin",     xmin.data());
        tree->SetBranchAddress("sc_ymin",     ymin.data());
        tree->SetBranchAddress("sc_xmax",     xmax.data());
        tree->SetBranchAddress("sc_ymax",     ymax.data());
        tree->SetBranchAddress("sc_size",     scsize.data());
        tree->SetBranchAddress("sc_nhits",    scnhits.data());
        tree->SetBranchAddress("sc_theta",    v_sc_theta.data());
        tree->SetBranchAddress("sc_rms",      v_sc_rms.data());
        tree->SetBranchAddress("sc_tgausssigma", v_sc_tgausssigma.data());

        // Prepare analysis result containers
        vector<int> BeginScPix, EndScPix;
        // (BeginScallPix and EndScallPix not used in logic, omitted)
        double phi_DIR = 0, phi_DIR_deg = 0, phi_PCA = 0, phi = 0;
        double phi_HT_skew = 0, phi_HT_maxpx = 0, phi_HT_maxpxRebin = 0;
        double phi_HT_peakprof = 0, phi_HT_maxprof = 0, phi_HT_integral = 0;
        double scint = 0;
        double skewness_L = 0, skewness_T = 0;
        double kurtosis_L = 0, kurtosis_T = 0;
        double xbar = 0, ybar = 0;
        double xl = 0, yl = 0, xr = 0, yr = 0;
        double reco_theta = 0;
        double x_mean = 0, y_mean = 0, x_min = 0, x_max = 0, y_min = 0, y_max = 0;
        double reco_sc_rms = 0, reco_sc_tgausssigma = 0;
        int procTime = 0;
        double profileMean = 0, profileSkew = 0;
        bool puCandidate = false;
        // Impact point variables
        double xIP = 0, yIP = 0;
        //stokes parameters (ik, qk, uk) not used in this worker logic
        double ik = 0, qk = 0, uk = 0;
        // (peakslong, peakstrans, etc. are not explicitly used further in code)
        int npeaks = 0;
        double peak_density = 0, tracklength = 0, track_width = 0;
        double recolength = 0, recowidth = 0;
        // Prepare output file for this worker
        std::string tempFilename = outputDir + "/track_info_temp_worker" + std::to_string(workerId) + ".root";
        TFile *fOut = TFile::Open(tempFilename.c_str(), "RECREATE");
        if (!fOut || fOut->IsZombie()) {
            std::cerr << "Worker " << workerId << ": cannot create output file\n";
            if (fOut) fOut->Close();
            fIn->Close();
            return std::make_pair(-1LL, -1LL);
        }
        TTree *out_tree = new TTree("track_info", "track_info");
        // Create output tree branches (same structure as original output)
        out_tree->Branch("run",          &run);
        out_tree->Branch("event",        &event);
        out_tree->Branch("nSc",          &nSc);
        out_tree->Branch("nSc_red",      &nSc_red);
        out_tree->Branch("Integral",     &scint);
        out_tree->Branch("ScSize",       &sc_npix);
        out_tree->Branch("AngleDIR",     &phi_DIR);
        out_tree->Branch("DegreeDIR",    &phi_DIR_deg);
        out_tree->Branch("AnglePCA",     &phi_PCA);
        out_tree->Branch("RecoTheta",    &reco_theta);
        out_tree->Branch("RecoScRMS",    &reco_sc_rms);
        out_tree->Branch("RecoScTGaussSigma", &reco_sc_tgausssigma);
        out_tree->Branch("X_ImpactPoint",&xIP);
        out_tree->Branch("Y_ImpactPoint",&yIP);
        out_tree->Branch("Xmean",        &x_mean);
        out_tree->Branch("Ymean",        &y_mean);
        out_tree->Branch("Xmin",        &x_min);
        out_tree->Branch("Ymin",        &y_min);
        out_tree->Branch("Xmax",        &x_max);
        out_tree->Branch("Ymax",        &y_max);
        out_tree->Branch("XBar",        &xbar);
        out_tree->Branch("YBar",        &ybar);
        out_tree->Branch("procTime",    &procTime);
        out_tree->Branch("ProfileMean", &profileMean);
        out_tree->Branch("ProfileSkew", &profileSkew);
        out_tree->Branch("PileUpCandidate", &puCandidate);
        out_tree->Branch("ik", &ik);
        out_tree->Branch("qk", &qk);
        out_tree->Branch("uk", &uk);

        // Determine this worker's event range [start, end)
        Long64_t startIndex = (workerId * totalEntries) / nWorkers;
        Long64_t endIndex   = ((Long64_t)(workerId + 1) * totalEntries) / nWorkers;
        if (workerId == (int)nWorkers - 1) {
            endIndex = totalEntries;  // ensure last worker goes to the end
        }

        long long localCounter = 0;
        long long localPileUpCounter = 0;
        for (Long64_t k = startIndex; k < endIndex; ++k) {
            tree->GetEntry(k);
            // Compute cluster pixel index ranges for this event
            ScIndicesElem(nSc, Nredpix, sc_redpixID.data(), nSc_red, BeginScPix, EndScPix);
            // Compute cluster pixel index ranges for this event
            ScIndicesElem(nSc, Nredpix, sc_redpixID.data(), nSc_red, BeginScPix, EndScPix);
            // Loop over superclusters in this event
            int pixcounter = 0;

            for (int i = 0; i < nSc_red; ++i) {
                // Set per-cluster variables from arrays
                scint       = (i < (int)integral.size() ? integral[i] : 0);
                recolength  = (i < (int)length.size()   ? length[i]   : 0);
                recowidth   = (i < (int)width.size()    ? width[i]    : 0);
                sc_npix     = (i < (int)scnhits.size()  ? (Int_t)scnhits[i] : 0);
                reco_theta  = (i < (int)v_sc_theta.size()? v_sc_theta[i]: 0);
                x_mean      = (i < (int)xmean.size()    ? xmean[i]    : 0);
                y_mean      = (i < (int)ymean.size()    ? ymean[i]    : 0);
                x_min       = (i < (int)xmin.size()     ? xmin[i]     : 0);
                x_max       = (i < (int)xmax.size()     ? xmax[i]     : 0);
                y_min       = (i < (int)ymin.size()     ? ymin[i]     : 0);
                y_max       = (i < (int)ymax.size()     ? ymax[i]     : 0);
                reco_sc_rms = (i < (int)v_sc_rms.size() ? v_sc_rms[i] : 0);
                reco_sc_tgausssigma = (i < (int)v_sc_tgausssigma.size()? v_sc_tgausssigma[i] : 0);

                pixcounter += EndScPix[i] - BeginScPix[i];
                //! **Selection condition**: filter events based on properties (same as original)
                std::map<std::string, double> event_vars = {
                    {"x_mean", x_mean},
                    {"y_mean", y_mean},
                    {"scint", scint},
                    // ...add more event variables here as needed
                };      
                bool pass_all = true;
                for (const auto& kv : config) {
                    std::string key = kv.first;
                    if (key.size() > 4 && key.substr(key.size() - 4) == "_min") {
                        std::string var = key.substr(0, key.size() - 4);
                        if (event_vars.count(var) && event_vars[var] < std::stod(kv.second)) pass_all = false;
                    }
                    if (key.size() > 4 && key.substr(key.size() - 4) == "_max") {
                        std::string var = key.substr(0, key.size() - 4);
                        if (event_vars.count(var) && event_vars[var] > std::stod(kv.second)) pass_all = false;
                    }
                }
                if (pass_all)
                {//!
                    // Start timing the track processing
                    auto t0 = std::chrono::high_resolution_clock::now();
                    // Construct Analyzer for this track (use false for no rebinning, as in original)
                    Analyzer Traccia(Form("Track%lli_event%i_run%i_entry%lli", localCounter, event, run, k),
                                    XPix.data(), YPix.data(), ZPix.data(), 
                                    BeginScPix[i], EndScPix[i], false);
                    if (Traccia.Getbinmax() <= 0) {
                        // Skip if track is empty
                        continue;
                    }
                    auto t1 = std::chrono::high_resolution_clock::now();
                    // Set directionality parameters
                    Traccia.SetWScal((float)wfactor);
                    Traccia.SetNPIP(NPIP);
                    // Apply threshold and remove noise
                    Traccia.ApplyThr((double)threshold);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    Traccia.RemoveNoise((double)remove_noise_value);
                    auto t5 = std::chrono::high_resolution_clock::now();
                    // Check for pile-up candidate:contentReference[oaicite:5]{index=5}
                    auto t5_pileup_start = std::chrono::high_resolution_clock::now();
                    puCandidate = Traccia.PileUpCandidate(false, localCounter, true, 0., 0.);
                    if (puCandidate) {
                        localPileUpCounter++;
                        // (Optionally save track image for pile-up, disabled in original)
                        //Traccia.TrackProfilePlotSave(Form("Track%lli_event%i_run%i_entry%lli.png", localCounter, event, run, k));
                        continue;  // skip further analysis and do not fill output for this track
                    }
                    auto t5_pileup_end = std::chrono::high_resolution_clock::now();
                    // Compute track profile statistics
                    Traccia.GetTrackProfileStats(profileMean, profileSkew);
                    // Determine impact point and scaled track (saves internal images, as in original)
                    Traccia.ImpactPoint(Form("TrackIPRegion%lli_run%i_evt%lli", k, run, localCounter));
                    auto t6 = std::chrono::high_resolution_clock::now();
                    Traccia.ScaledTrack(Form("TrackScaled%lli_run%i_evt%lli", k, run, localCounter));
                    auto t7 = std::chrono::high_resolution_clock::now();
                    // Compute track direction
                    Traccia.Direction();
                    auto t8 = std::chrono::high_resolution_clock::now();
                    // Improve angle and build line direction (not timed in output)
                    Traccia.ImprCorrectAngle();
                    Traccia.BuildLineDirection();
                    Traccia.BuildStokesParams();                    auto t10 = std::chrono::high_resolution_clock::now();
                    // Retrieve results from Analyzer
                    xIP         = Traccia.GetXIP();
                    yIP         = Traccia.GetYIP();
                    phi_DIR     = Traccia.GetDir();
                    phi_DIR_deg = phi_DIR * (180.0 / M_PI);
                    phi_PCA     = Traccia.AngleLineMaxRMS();
                    xbar        = Traccia.GetXbar();
                    ybar        = Traccia.GetYbar();
                    ik          = Traccia.GetIk();
                    qk          = Traccia.GetQk();
                    uk          = Traccia.GetUk();
                    // End timing for this track
                    auto tEnd = std::chrono::high_resolution_clock::now();
                    // Print progress for a single random worker only (e.g., worker 0)
                    static int random_worker = rand() % nWorkers;
                    static int last_percent_printed = -1;
                    if (workerId == random_worker) {
                        double frac = (double)(k - startIndex + 1) / (double)(endIndex - startIndex);
                        int percent = static_cast<int>(frac * 100.0);
                        if (percent % 10 == 0 && percent != last_percent_printed && percent != 0) {
                            std::cout << "[Worker " << workerId << "] Progress: " << percent << "% ("
                                    << (k - startIndex + 1) << "/" << (endIndex - startIndex) << ")" << std::endl;
                            last_percent_printed = percent;
                        }
                    }
                    if (k % printout == 0) {
                        Traccia.SavePicDir(Form("Track%lli_event%i_run%i_entry%lli.png", 
                                                localCounter, event, run, k));
                        std::cout << "Processing entry " << k << " / " << totalEntries << std::endl;
                        std::cout << "counter: " << localCounter << std::endl;
                        std::cout << "XIP: " << xIP << "  YIP: " << yIP << std::endl;
                        std::cout << "XIPPrev: " << Traccia.GetXIPPrev() 
                                << "  YIPPrev: " << Traccia.GetYIPPrev() << std::endl;
                        std::cout << "Degree: " << phi_DIR_deg 
                                << "  tan(angle): " << std::tan(phi_DIR) << std::endl;
                        // Compute durations in microseconds for each step
                        long long dt_ctor          = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
                        long long dt_applyThr      = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t1).count();
                        long long dt_removeNoise   = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
                        long long dt_pileUpCheck   = std::chrono::duration_cast<std::chrono::microseconds>(t5_pileup_end - t5_pileup_start).count();
                        long long dt_impactPoint   = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5_pileup_end).count();
                        long long dt_scaledTrack   = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
                        long long dt_direction     = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
                        long long dt_buildLineDir  = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t8).count();
                        long long dt_total         = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
                        // Helper to print timing and percentage
                        auto printTime = [&](const char* label, long long dt) {
                            double pct = 100.0 * (double)dt / (double)dt_total;
                            std::cout << label << dt << " us (" << pct << " %)" << std::endl;
                        };
                        printTime("Time Analyzer ctor:         ", dt_ctor);
                        printTime("Time ApplyThr:              ", dt_applyThr);
                        printTime("Time RemoveNoise:           ", dt_removeNoise);
                        printTime("Time PileUpCandidate:       ", dt_pileUpCheck);
                        printTime("Time ImpactPoint:           ", dt_impactPoint);
                        printTime("Time ScaledTrack:           ", dt_scaledTrack);
                        printTime("Time Direction:             ", dt_direction);
                        printTime("Time BuildLineDirection:    ", dt_buildLineDir);
                        std::cout << "TOTAL time for all steps: " << dt_total << " us (100%)\n" << std::endl;
                    }
                    // Record total processing time for this track and fill output tree
                    procTime = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - t0).count();
                    out_tree->Fill();
                } // end if (selection condition)
                // Increment counter for each cluster processed (skipped clusters via continue do not increment)
                localCounter++;
            } // end supercluster loop
            // Resize sc_redpixID for safety (per original code)
            sc_redpixID.resize(nSc);
        } // end event loop

        // Write this worker's result tree and close files
        out_tree->Write();
        fOut->Close();
        fIn->Close();
        // Return the count of processed tracks and pile-up tracks from this worker
        return std::make_pair(localCounter, localPileUpCounter);
    }, workerIds);  // Map lambda over all worker IDs

    // Check for errors from any worker
    bool workerError = false;
    for (auto& pr : results) {
        if (pr.first < 0 || pr.second < 0) {
            workerError = true;
            break;
        }
    }
    if (workerError) {
        std::cerr << "Error: one or more worker processes failed.\n";
        // Clean up any partial output files
        for (unsigned int i = 0; i < nWorkers; ++i) {
            std::string tempName = outputDir + "/track_info_temp_worker" + std::to_string(i) + ".root";
            std::filesystem::remove(tempName);
        }
        return 1;
    }

    // Merge all worker output ROOT files into the final output using TFileMerger
    TFileMerger merger(kFALSE); // "kFALSE" for not verbose

    // Set the output file
    merger.OutputFile(finalOutputPath.c_str(), "RECREATE");

    // Add all temporary files
    for (unsigned int i = 0; i < nWorkers; ++i) {
        std::string tempName = outputDir + "/track_info_temp_worker" + std::to_string(i) + ".root";
        merger.AddFile(tempName.c_str());
    }

    // Actually merge!
    if (!merger.Merge()) {
        std::cerr << "Error: ROOT TFileMerger failed to merge worker files!" << std::endl;
        return 1;
    }

    // Open final file to create "Tracks" directory if needed (optional, as before)
    TFile *fmerged = TFile::Open(finalOutputPath.c_str(), "UPDATE");
    if (fmerged && !fmerged->IsZombie()) {
        if (!fmerged->GetDirectory("Tracks")) {
            fmerged->mkdir("Tracks");
        }
        fmerged->Close();
    }


    // Remove temporary worker files
    for (unsigned int i = 0; i < nWorkers; ++i) {
        std::string tempName = outputDir + "/track_info_temp_worker" + std::to_string(i) + ".root";
        std::filesystem::remove(tempName);
    }

    // Compute total counters for pile-up percentage from results
    long long totalProcessedTracks = 0;
    long long totalPileUpTracks = 0;
    for (auto& pr : results) {
        totalProcessedTracks += pr.first;
        totalPileUpTracks += pr.second;
    }
    double pileUpRatio = (totalProcessedTracks > 0 ? (double)totalPileUpTracks / (double)totalProcessedTracks : 0.0);
    std::cout << "pile up percentage: " << pileUpRatio << std::endl;

    // Write config parameters as individual branches (one entry) in a new TTree in the output file
    TFile *fout = TFile::Open(finalOutputPath.c_str(), "UPDATE");
    if (fout && !fout->IsZombie()) {
        TTree *config_tree = new TTree("config", "Configuration parameters used for this analysis");
        // Store config values with correct types
        std::map<std::string, int>    int_vars;
        std::map<std::string, double> double_vars;
        std::map<std::string, std::string> string_vars;
        std::map<std::string, void*> branch_ptrs;
        std::vector<std::string> branch_types;
        for (const auto& kv : config) {
            // Try to parse as int
            try {
                size_t idx;
                int val = std::stoi(kv.second, &idx);
                if (idx == kv.second.size()) {
                    int_vars[kv.first] = val;
                    branch_ptrs[kv.first] = &int_vars[kv.first];
                    branch_types.push_back("I");
                    continue;
                }
            } catch (...) {}
            // Try to parse as double
            try {
                size_t idx;
                double val = std::stod(kv.second, &idx);
                if (idx == kv.second.size()) {
                    double_vars[kv.first] = val;
                    branch_ptrs[kv.first] = &double_vars[kv.first];
                    branch_types.push_back("D");
                    continue;
                }
            } catch (...) {}
            // Otherwise, treat as string
            string_vars[kv.first] = kv.second;
            branch_ptrs[kv.first] = (void*)string_vars[kv.first].c_str();
            branch_types.push_back("C");
        }
        // Create branches
        size_t idx = 0;
        for (const auto& kv : config) {
            const std::string& name = kv.first;
            const std::string& type = branch_types[idx++];
            if (type == "I") {
                config_tree->Branch(name.c_str(), (int*)branch_ptrs[name], (name + "/I").c_str());
            } else if (type == "D") {
                config_tree->Branch(name.c_str(), (double*)branch_ptrs[name], (name + "/D").c_str());
            } else {
                config_tree->Branch(name.c_str(), (char*)branch_ptrs[name], (name + "/C").c_str());
            }
        }
        // Fill the tree with one entry
        config_tree->Fill();
        config_tree->Write();
        fout->Close();
    }

    return 0;
}
