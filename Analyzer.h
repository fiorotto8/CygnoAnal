/**
 * @file Analyzer.h
 * @brief Provides the Analyzer class for analyzing two-dimensional histogram data.
 *
 * The Analyzer class offers a comprehensive suite of tools for the analysis of
 * two-dimensional histograms, particularly useful in fields such as particle
 * physics. It enables operations like noise removal, histogram scaling,
 * barycenter calculation, peak finding, and more. The class supports deep
 * copying through a copy constructor, ensuring safe copy semantics for
 * objects managing dynamic memory.
 */

#ifndef ANALYZER_H
#define ANALYZER_H

#include <vector>
#include <cmath>
#include "TMatrixD.h"
#include "TDecompSVD.h"
class TH2F;
class TF1;
class TGraph;
class TH1D;

class Analyzer {
public:
    /**
     * Default constructor.
     * Initializes an empty Analyzer object.
     */
    Analyzer();

    /**
     * Copy constructor.
     * @param source The source Analyzer object to copy.
     */
    Analyzer(const Analyzer& source);

    /**
     * Specialized constructor for creating an Analyzer with a resized track.
     * @param nometh2 Name for the new histogram.
     * @param npixel Number of pixels for the new track.
     * @param Tracklarge Original histogram with the track.
     * @param npixelorig Number of pixels in the original histogram.
     */
    Analyzer(const char* nometh2, int npixel ,TH2F* Tracklarge, int npixelorig );

    /**
     * Constructor for generic track analysis.
     * @param nometh2 Name for the new histogram.
     * @param Tracklarge Histogram of the track to analyze.
     */
    Analyzer(const char* nometh2, TH2F *Tracklarge);

    /**
     * Constructor taking arrays of x, y, z coordinates of the track.
     * @param nometh2 Name for the new histogram.
     * @param X Array of x coordinates.
     * @param Y Array of y coordinates.
     * @param Z Array of z coordinates (intensities).
     * @param B Start index for the track hits.
     * @param E End index for the track hits.
     */
    Analyzer(const char* nometh2,int* X,int* Y,float* Z, int B,int E);

    /**
     * Destructor.
     * Cleans up dynamic resources managed by the Analyzer.
     */
    virtual ~Analyzer();

    //initializer
    void BuildLineMaxRMS();
    void BuildLineDirection();

    //Simple returners
    inline double GetIntegral() const {return fintegral;}
    inline double GetHeight() const {return fheight;}
    inline double GetRadius() const {return fradius;}

    inline TH2F* GetHistoTrack() const {return fTrack;}
    inline TH2F* GetHistoTrackTail() const {return fTrackTail;}
    inline TH2F* GetHistoScaledTrack() const {return fScaledTrack;}

    inline void GetCentre(double &x, double &y) const {x=fxcentr; y=fycentr;};
    inline double GetXbar() const {return fXbar;}
    inline double GetYbar() const  {return fYbar;}
    inline double GetRMSOnMainAxis() const {return fRMSOnMainAxis;}
    inline double GetSkewOnMainAxis() const {return fSkewOnLine;}
    inline double GetPointSkew(double X, double Y) const {return ( (X-fXbar)*cos(fPhiMainAxis) + (Y-fYbar)*sin(fPhiMainAxis) )/fSkewOnLine;}
    inline TF1* GetLineMaxRMS() const  {return fLineMaxRMS;}
    inline double GetXIP() const {return fXIP;}
    inline double GetYIP() const {return fYIP;}
    inline double GetDir() const {return fPhiDir;}
    inline double PDistCm(double X, double Y) const {return sqrt( ( (X-fXbar)*(X-fXbar) + (Y-fYbar)*(Y-fYbar) ));}

    //Setters for directionality
    inline void SetWScal(float a) {fwScal=a;}
    inline void SetNPIP(int a) {fNPIP=a;}

    // Public interface for Analyzer functionality
    void Reset();
    void Integral();
    void SavetoFile(const char* nometh2);
    void SavePic(const char* nometh2);
    void SavePicDir(const char* nomepic);
    void SaveRootFile(const char* nomefile);

    void PlotandSavetoFileCOLZ(const char* nometh2);
    void PlotandSavetoFileDirectionalFull(const char* nomepic);

    void Barycenter();
    void Barycenter(TH2F* Tr,double *X, double *Y);
    double AngleLineMaxRMS();
    double RMSOnLine(double Phi);
    double SkewOnMainAxis();
    void RemoveNoise(double a=0);
    void ApplyThr();
    void ImpactPoint(const char* nometh2);
    void ScaledTrack(const char* nometh2);

    void Direction();
    void ImprCorrectAngle();

    void Edges(double &Xl, double &Yl, double &Xr, double &Yr, double slope);
    TH1D* FillProfile(bool longitudinal, float x1=0, float x2=2304);
    TH1D* CutProfile(TH1D* profile, double height=0.0025);
    TH1D* FillProfileX();
    TH1D* FillProfileY();
    void AnglePCA(double &ang);
    void FindNPeaks(TH1D* h, std::vector<std::pair<double,double>> &foundpeaks);
    void FindPeak(double &xpeak, double &ypeak, double &xpeak_rebin, double &ypeak_rebin);
    // Additional methods and their documentation...
    int Execute_Atul_script(std::string pyvers, std::string inputfile, std::string outfolder, int entries, bool plot=false, bool text=false ) const;

private:
    // Private member variables
    //histogram parameters
    int fminx,fminy,fmaxx,fmaxy;

    double fintegral;
    double fradius;
    double fheight;
    double fxcentr;
    double fycentr;

    TH2F *fTrack;
    TH2F *fTrackTail;
    TH2F *fScaledTrack;
    //TH2 Parameters
    int fnpixelx;
    int fnpixely;

    //directionaliy parameters
    int fNPIP;
    float fwScal;
    double fXbar;
    double fYbar;
    TGraph* fBarPlot;
    double fPhiMainAxis;
    TF1* fLineMaxRMS;
    double fRMSOnMainAxis;
    double fSkewOnLine;

    double fXIPPrev;
    double fYIPPrev;
    double fXIP;
    double fYIP;
    TGraph* fIPPlot;

    double fPhiDir;
    TF1* fLineDirection;
};

#endif // ANALYZER_H