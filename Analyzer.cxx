/**
 * @file Analyzer.cpp
 * @brief Implementation of the Analyzer class for advanced histogram analysis.
 *
 * This file provides the implementation of the Analyzer class, which is designed
 * for the analysis and manipulation of two-dimensional histograms, particularly
 * in the fields of particle physics and similar research areas. It includes
 * functionalities for noise reduction, peak finding, barycenter calculation,
 * and other analytical capabilities to aid in the study and interpretation
 * of experimental or simulation data.
 */

#include <cstdio>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <Riostream.h>
#include <TFile.h>
#include "Analyzer.h"
#include "TAxis.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TSpectrum.h>
#include <TMatrixD.h>
#include "TDecompSVD.h"
#include "TLatex.h"
#include "TLine.h"

using namespace std;

/**
 * A utility function for calculating the Root Mean Square (RMS) on a line
 * defined by its slope and position with respect to a reference point.
 *
 * @param XBar The x-coordinate of the reference point.
 * @param YBar The y-coordinate of the reference point.
 * @param Phi The angle (in radians) defining the slope of the line.
 * @return The calculated RMS value along the line.
 */
double RMSOnLine(double XBar, double YBar, double Phi);

// Constructors and Destructor
/**
 * Default constructor. Initializes an Analyzer object with default values.
 */
Analyzer::Analyzer():
// Initialization list for setting default values
fRange(0.),
fminx(0.),
fminy(0.),
fmaxx(0.),
fmaxy(0.),
fintegral(0.),
fradius(0.),
fheight(0.),
fxcentr(0.),
fycentr(0.),
fTrack(nullptr),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fnpixelx(0),
fnpixely(0),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fPhiMainAxis(0.),
fLineMaxRMS(nullptr),
fRMSOnMainAxis(0.),
fSkewOnLine(0.),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr), 
fPhiDir(0.),
fLineDirection(nullptr)
{
}


/**
 * Copy constructor. Creates a deep copy of an existing Analyzer object.
 *
 * @param source A reference to the source Analyzer object to copy.
 */
Analyzer::Analyzer(const Analyzer& source):
// Initialization list for deep copying member variables
fRange(0.),
fminx(source.fminx),
fminy(source.fminy),
fmaxx(source.fmaxx),
fmaxy(source.fmaxy),
fintegral(source.fintegral),
fradius(source.fradius),
fheight(source.fheight),
fxcentr(source.fxcentr),
fycentr(source.fycentr),
fnpixelx(source.fnpixelx),
fnpixely(source.fnpixely),
fNPIP(source.fNPIP),
fwScal(source.fwScal),
fXbar(source.fXbar),
fYbar(source.fYbar),
fPhiMainAxis(source.fPhiMainAxis),
fRMSOnMainAxis(source.fRMSOnMainAxis),
fSkewOnLine(source.fSkewOnLine),
fXIPPrev(source.fXIPPrev),
fYIPPrev(source.fYIPPrev),
fXIP(source.fXIP),
fYIP(source.fYIP),
fPhiDir(source.fPhiDir)
{
  // Deep copy logic, including conditional cloning of TH2F objects and others
  if(source.fTrack==nullptr)	   fTrack=nullptr;
  else fTrack=(TH2F*)source.fTrack->Clone("Trackcopy");
  if(source.fTrackTail==nullptr)	   fTrackTail=nullptr;
  else fTrackTail=(TH2F*)source.fTrackTail->Clone("TrackTailcopy");
  if(source.fScaledTrack==nullptr)	   fScaledTrack=nullptr;
  else fScaledTrack=(TH2F*)source.fScaledTrack->Clone("TrackScaledcopy");
  if(source.fBarPlot==nullptr)	   fBarPlot=nullptr;
  else fBarPlot= new TGraph(*source.fBarPlot);
  if(source.fIPPlot==nullptr)	   fIPPlot=nullptr;
  else fIPPlot=new TGraph(*source.fIPPlot);
  if(source.fLineMaxRMS==nullptr)	   fLineMaxRMS=nullptr;
  else fLineMaxRMS=(TF1*)source.fLineMaxRMS->Clone("LineMaxRMScopy");
  if(source.fLineDirection==nullptr)	   fLineDirection=nullptr;
  else fLineDirection=(TF1*)source.fLineDirection->Clone("LineDirectioncopy");

}

/**
 * Specialized constructor for creating an Analyzer with a resized track.
 * Does not evaluate the integral in this constructor.
 *
 * @param nometh2 Name for the new histogram.
 * @param npixel Number of pixels for the new track.
 * @param Tracklarge Pointer to the original 2D histogram with the track.
 * @param npixelorig Number of pixels in the original 2D histogram.
 */
Analyzer::Analyzer(const char* nometh2, int npixel, TH2F *Tracklarge, int npixelorig):
// Initialization list and constructor logic
fRange(0.),
fintegral(0.),
fradius(0.),
fheight(0.),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fPhiMainAxis(0.),
fLineMaxRMS(nullptr),
fRMSOnMainAxis(0.),
fSkewOnLine(0.),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr),
fPhiDir(0.),
fLineDirection(nullptr)
{
	fycentr=Tracklarge->GetMaximumBin()/(npixelorig+2);
	fxcentr=Tracklarge->GetMaximumBin()%(npixelorig+2);

	fminx=npixelorig+100;
	fminy=npixelorig+100;
	fmaxx=-1;
	fmaxy=-1;
	for(int i=fxcentr-npixel/2;i<fxcentr+npixel/2;i++) {

		for(int j=fycentr-npixel/2;j<fycentr+npixel/2;j++) {

			double x=Tracklarge->GetBinContent(i,j);
			
      if(i<fminx && x>0)  fminx=i;
			if(i>fmaxx && x>0)  fmaxx=i;
			if(j<fminy && x>0)  fminy=j;
			if(j>fmaxy && x>0)  fmaxy=j;
		}
	}

	fnpixelx=fmaxx-fminx+1+10;
	fnpixely=fmaxy-fminy+1+10;
	fTrack=new TH2F(nometh2,nometh2,fnpixelx,0,fnpixelx,fnpixely,0,fnpixely);
	for(int j=1;j<=fnpixelx;j++) {

		for(int k=1;k<=fnpixely;k++) {

			fTrack->SetBinContent(j,k,Tracklarge->GetBinContent(fminx+j-1-5,fminy+k-1-5));
		}
	}

	TH1D *tax=fTrack->ProjectionX();
	tax->Fit("gaus","Q");
	fxcentr=tax->GetFunction("gaus")->GetParameter(1);
	TH1D *tay=fTrack->ProjectionY();
	tay->Fit("gaus","Q");
	fycentr=tay->GetFunction("gaus")->GetParameter(1);
	delete tax;
	delete tay;

	for(int j=1;j<=fnpixelx;j++) {

		for(int k=1;k<=fnpixely;k++) {

			if(fTrack->GetBinContent(j,k)>0) {

				double x=sqrt((fycentr-k)*(fycentr-k) + (fxcentr-j)*(fxcentr-j) );
				if(x>fradius)  fradius=x;
			}
		}
	}
}

/**
 * Constructor for generic track analysis.
 *
 * @param nometh2 Name for the new histogram.
 * @param Tracklarge Pointer to the histogram of the track to be analyzed.
 */
Analyzer::Analyzer(const char* nometh2, TH2F *Tracklarge):
// Initialization list and constructor logic
fRange(0.),
fradius(0.),
fheight(0.),
fxcentr(0),
fycentr(0),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr),
fPhiDir(0.),
fLineDirection(nullptr)
{
	fintegral=0.;
	fminx=3000;
	fminy=3000;
	fmaxx=0;
	fmaxy=0;

	for(int i=0;i<Tracklarge->GetXaxis()->GetNbins();i++) {
    
    for(int j=0;j<Tracklarge->GetYaxis()->GetNbins();j++) {

      double x=Tracklarge->GetBinContent(i,j);

      if(i<fminx && x>0)  fminx=i;
      if(i>fmaxx && x>0)  fmaxx=i;
      if(j<fminy && x>0)  fminy=j;
      if(j>fmaxy && x>0)  fmaxy=j;
    }
	}

	fmaxx=fmaxx+30;
	fminx=fminx-30;
	fmaxy=fmaxy+30;
	fminy=fminy-30;

	fnpixelx=fmaxx-fminx;
	fnpixely=fmaxy-fminy;

	//std::cout << fnpixelx<<"\t"<<fminx <<"\t"<< fmaxx << "\t" << fnpixely <<"\t"<< fminy <<"\t"<< fmaxy << std::endl;

	fTrack=new TH2F(Form("A%s",nometh2),Form("A%s",nometh2),fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
	int np=0;

	for(int j=0;j<=fnpixelx;j++)
	{
		for(int k=0;k<=fnpixely;k++)
		{
      if(Tracklarge->GetBinContent(fminx+j,fminy+k)>0){
        fTrack->SetBinContent(j,k,Tracklarge->GetBinContent(fminx+j,fminy+k));
        fintegral+=Tracklarge->GetBinContent(fminx+j,fminy+k);
      }
		}
	}

	fPhiMainAxis=AngleLineMaxRMS();
	BuildLineMaxRMS();		//defines fLineMaxRMS
	fRMSOnMainAxis=GetRMSOnMainAxis();
	fSkewOnLine=SkewOnMainAxis();

	TCanvas* OriTrack = new TCanvas();
	Tracklarge->GetXaxis()->SetRangeUser(fminx,fmaxx);
	Tracklarge->GetYaxis()->SetRangeUser(fminy,fmaxy);
	Tracklarge->Draw("COLZ");
	OriTrack->SaveAs(Form("Tracks/%s.png",Tracklarge->GetName()));
	delete OriTrack;

	//fTrack->Rebin2D(2,2);
}

/**
 //! Usually this is the standard Constructor for the Analyzer class.
 * Constructor taking the arrays of x, y, z coordinates of the track.
 *
 * @param nometh2 Name for the new histogram.
 * @param X Array of x coordinates.
 * @param Y Array of y coordinates.
 * @param Z Array of z coordinates (intensities).
 * @param B Start index for the track hits.
 * @param E End index for the track hits.
 */
Analyzer::Analyzer(const char* nometh2, int* X, int* Y, float* Z, int B, int E):
// Initialization list and constructor logicfradius(0.),
fRange(0.),
fheight(0.),
fxcentr(0),
fycentr(0),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr),
fPhiDir(0.),
fLineDirection(nullptr)
{
	fintegral=0.;
  /*
  * Debugging: Check if B and E are within valid range
  cout << "Debug Info:" << endl;
  cout << "B (start index): " << B << endl;
  cout << "E (end index): " << E << endl;
  cout << "Number of elements (E-B): " << E-B << endl;

  * Check if B and E are within the bounds of X and Y arrays.
  * Assuming npix is the size of X and Y arrays. You need to replace npix with the actual size if it's different.
  cout << "Checking if B and E are within valid range..." << endl;
  if (B < 0 || E > 500000) {
      cout << "Error: B is less than 0 or E is greater than the size of the arrays." << endl;
  } else if (E <= B) {
      cout << "Error: E is less than or equal to B, which means the range is invalid." << endl;
  } else {
      cout << "B and E indices appear to be valid." << endl;
  } 
  */

  //  Original version
	fminx = TMath::MinElement(E-B,&X[B]);
	fmaxx = TMath::MaxElement(E-B,&X[B]);
	fminy = TMath::MinElement(E-B,&Y[B]);
	fmaxy = TMath::MaxElement(E-B,&Y[B]);

	fmaxx=fmaxx+10;
	fminx=fminx-10;
	fmaxy=fmaxy+10;
	fminy=fminy-10;

	fnpixelx=fmaxx-fminx;
	fnpixely=fmaxy-fminy;

	//std::cout << fnpixelx<<"\t"<<fminx <<"\t"<< fmaxx << "\t" << fnpixely <<"\t"<< fminy <<"\t"<< fmaxy << std::endl;

	fTrack=new TH2F(Form("A%s",nometh2),Form("A%s",nometh2),fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
	for(int i=B;i<E;i++) {
    
    if (Z[i]>0) {
      
      fTrack->SetBinContent(fTrack->GetXaxis()->FindBin(X[i]),fTrack->GetYaxis()->FindBin(Y[i]),Z[i]);
      fintegral+=Z[i];
    }
	}
  /*
  *   David's version, careful with rotations

	fmaxx = 2305 - TMath::MinElement(E-B,&X[B]);
	fminx = 2305 - TMath::MaxElement(E-B,&X[B]);
	fmaxy = 2305 - TMath::MinElement(E-B,&Y[B]);
	fminy = 2305 - TMath::MaxElement(E-B,&Y[B]);

  fmaxx=fmaxx+30;
	fminx=fminx-30;
	fmaxy=fmaxy+30;
	fminy=fminy-30;

  fnpixelx=fmaxx-fminx;
	fnpixely=fmaxy-fminy;

  fTrack=new TH2F(Form("%s",nometh2),Form("%s",nometh2),fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);

	for(int i=B;i<E;i++){

    if (Z[i]>0) {
      
      fTrack->SetBinContent(fTrack->GetXaxis()->FindBin( 2305 - X[i]),fTrack->GetYaxis()->FindBin( 2305 - Y[i]),Z[i]); // cahnge david (rotation)
      fintegral+=Z[i];
    }
	}
  */

	//fTrack->Rebin2D(2,2);
	fPhiMainAxis=AngleLineMaxRMS();
	BuildLineMaxRMS();
	fRMSOnMainAxis=GetRMSOnMainAxis();
	fSkewOnLine=SkewOnMainAxis();	
}

/**
 * Destructor. Cleans up dynamic resources managed by the Analyzer.
 */
Analyzer::~Analyzer()
{
  // Cleanup logic for dynamic resources
	if(fTrack!=nullptr)	 	 	 delete fTrack;
	if(fTrackTail!=nullptr)	 	 delete fTrackTail;
	if(fScaledTrack!=nullptr)	 delete fScaledTrack;
	if(fBarPlot!=nullptr)	 	 delete fBarPlot;
	if(fLineMaxRMS!=nullptr)	 delete fLineMaxRMS;
	if(fIPPlot!=nullptr)	 	 delete fIPPlot;
	if(fLineDirection!=nullptr)	 delete fLineDirection;
}


// Implementation of Other Member Functions

/**
 * Get the sigma of the track along the main axis in a range (Xbar, Ybar) - (-range, +range)
 * @return The calculated sigma (spread) along the main axis.
 */
double Analyzer::GetSigmaAroundBar() {
  
  double spread = 0;
  double totalWeight = 0;
  double sumDistances = 0;
  double sumDistancesSquared = 0;

  double cosPhi = cos(fPhiMainAxis);
  double sinPhi = sin(fPhiMainAxis);
  double cosPerpPhi = cos(fPhiMainAxis + TMath::Pi() / 2);
  double sinPerpPhi = sin(fPhiMainAxis + TMath::Pi() / 2);

  // Calculate the coordinates for the perpendicular lines
  double x1 = fXbar - fRange * cosPhi;
  double y1 = fYbar - fRange * sinPhi;
  double x2 = fXbar + fRange * cosPhi;
  double y2 = fYbar + fRange * sinPhi;

  for (int i = 1; i <= fnpixelx; ++i) {
    for (int j = 1; j <= fnpixely; ++j) {

      double binContent = fTrack->GetBinContent(i, j);
      if (binContent > 0) {

        double x = fTrack->GetXaxis()->GetBinCenter(i);
        double y = fTrack->GetYaxis()->GetBinCenter(j);

        // Project the distance onto the perpendicular line
        double distance1 = (x - x1) * cosPerpPhi + (y - y1) * sinPerpPhi;
        double distance2 = (x - x2) * cosPerpPhi + (y - y2) * sinPerpPhi;

        // Choose the projection that falls within the range
        double distance = (fabs(distance1) <= fRange) ? distance1 : distance2;

        // Accumulate the weighted sums if within range
        if (fabs(distance) <= fRange) {

          sumDistances += distance * binContent;
          sumDistancesSquared += distance * distance * binContent;
          totalWeight += binContent;
        }
      }
    }
  }

  if (totalWeight > 0) {

    double meanDistance = sumDistances / totalWeight;
    double meanDistanceSquared = sumDistancesSquared / totalWeight;
    spread = sqrt(meanDistanceSquared - meanDistance * meanDistance);
  }

  return spread;
}

/**
 * Initializes fLineMaxRMS based on current state.
 */
void Analyzer::BuildLineMaxRMS(){

  fLineMaxRMS= new TF1("LineMaxRMS","[0]*(x-[1])+[2]",0,2304);

  fLineMaxRMS->SetParameter(1,fXbar);
  fLineMaxRMS->SetParameter(2,fYbar);
  fLineMaxRMS->SetParameter(0,TMath::Tan(fPhiMainAxis));
}

/**
 * Initializes fLineDirection based on current state.
 */
void Analyzer::BuildLineDirection(){

  fLineDirection= new TF1("LineDirection","[0]*(x-[1])+[2]",0,2304);

  fLineDirection->SetParameter(1,fXIP);
  fLineDirection->SetParameter(2,fYIP);
  fLineDirection->SetParameter(0,TMath::Tan(fPhiDir));
}

/**
 * Resets the Analyzer object to its default state, freeing any dynamic resources.
 */
void Analyzer::Reset()
{
	fminx=0;
	fminy=0;
	fmaxx=0;
	fmaxy=0;
	fintegral=0;
	fheight=0;
	fradius=0;
	fxcentr=0;
	fycentr=0;
	fnpixelx=0;
	fnpixely=0;
	fNPIP=0;
	fwScal=0;
	fXbar=0;
  fYbar=0;
  fPhiMainAxis=0;
  fRMSOnMainAxis=0;
  fSkewOnLine=0;
  fXIPPrev=0;
  fYIPPrev=0;
  fXIP=0;
  fYIP=0;
  fPhiDir=0;

	delete fTrack;
	fTrack=nullptr;
  delete fTrackTail;
  fTrackTail=nullptr;
  delete fScaledTrack;
	fScaledTrack=nullptr;
	delete fLineDirection;
	fLineDirection=nullptr;
	delete fIPPlot;
	fIPPlot=nullptr;
	delete fLineMaxRMS;
	fLineMaxRMS=nullptr;
	delete fBarPlot;
	fBarPlot=nullptr;
}

/**
 * Calculates the integral of the track's histogram.
 */
void Analyzer::Integral()
{
	fintegral=0.;
	for(int i=fminx;i<fmaxx;i++) {

		for(int j=fminy;j<fmaxy;j++) {

		  fintegral+=fTrack->GetBinContent(i,j);
		}
	}
}

/**
 * Saves the current histogram to a ROOT file.
 *
 * @param nometh2 The name for the histogram in the file.
 */
void Analyzer::SavetoFile(const char* nometh2)
{
	fTrack->SetName(nometh2);
	fTrack->Write();
	return;
}

/**
 * Plots and saves the current canvas of the Track in COLZ to a ROOT file.
 *
 * @param nometh2 The name for the histogram in the file.
 */
void Analyzer::PlotandSavetoFileCOLZ(const char* nometh2)
{
  TCanvas* canv = new TCanvas("canv","canv",1500,1500);

  fTrack->Draw("COLZ");
  canv->SetName(nometh2);
  canv->Write();
  canv->DrawClone();
  delete canv;

  return;
}

/**
 * Generates and saves a pictorial representation of the current histogram.
 *
 * @param nomepic The filename for the saved picture.
 */
void Analyzer::SavePic(const char* nomepic)
{
  TCanvas* canv = new TCanvas("canv","canv",1500,1500);

  TGraph* g = new TGraph();
  g->SetPoint(0,fXbar,fYbar);
  g->SetMarkerStyle(8);

  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0,Form("NPx=%i",fnpixelx*fnpixely)); 	//assuming the track is a rectangle... doubtful
  l->AddEntry((TObject*)0,Form("Int=%f",fintegral));
  l->AddEntry((TObject*)0,Form("Dens=%f",fintegral/(fnpixelx*fnpixely)));
  l->AddEntry((TObject*)0,Form("Skew=%f",fSkewOnLine));
  l->AddEntry((TObject*)0,Form("SkewNorm=%f",fSkewOnLine/fintegral));

  l->AddEntry((TObject*)0,Form("RMS=%f",fRMSOnMainAxis));
  l->AddEntry((TObject*)0,Form("RMSNorm=%f",fRMSOnMainAxis/fintegral));

  //fTrack->SetName(nomepic);
  fTrack->Draw("COLZ");
  g->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");
  l->Draw("SAME");

  canv->SaveAs(Form("Tracks/%s",nomepic));

  delete canv;

  return;
}

/**
 * Saves a pictorial representation showing the directionality of the track, with lines indicating a specified range.
 *
 * @param nomepic The filename for the saved picture.
 */
void Analyzer::SavePicDirWithRange(const char* nomepic) {
  TCanvas* canv = new TCanvas("canv", "canv", 1500, 1500); // Adjusted canvas size to 1500x1500

  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0, Form("%f", fPhiDir / TMath::Pi() * 180));

  fTrack->Draw("COLZ");
  fBarPlot->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");

  // Calculate the coordinates for the range lines
  double cosPhi = cos(fPhiMainAxis);
  double sinPhi = sin(fPhiMainAxis);
  double cosPerpPhi = cos(fPhiMainAxis + TMath::Pi() / 2);
  double sinPerpPhi = sin(fPhiMainAxis + TMath::Pi() / 2);
  
  double x1 = fXbar - fRange * cosPhi;
  double y1 = fYbar - fRange * sinPhi;
  double x2 = fXbar + fRange * cosPhi;
  double y2 = fYbar + fRange * sinPhi;

  // Extend the lines to cover the visible range of the canvas
  double xMin = fTrack->GetXaxis()->GetXmin();
  double xMax = fTrack->GetXaxis()->GetXmax();
  double yMin = fTrack->GetYaxis()->GetXmin();
  double yMax = fTrack->GetYaxis()->GetXmax();

  TLine* line1 = new TLine(x1 + (xMin - fXbar) * cosPerpPhi, y1 + (yMin - fYbar) * sinPerpPhi,
                            x1 + (xMax - fXbar) * cosPerpPhi, y1 + (yMax - fYbar) * sinPerpPhi);
  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->Draw("SAME");

  TLine* line2 = new TLine(x2 + (xMin - fXbar) * cosPerpPhi, y2 + (yMin - fYbar) * sinPerpPhi,
                            x2 + (xMax - fXbar) * cosPerpPhi, y2 + (yMax - fYbar) * sinPerpPhi);
  line2->SetLineColor(kRed);
  line2->SetLineStyle(2);
  line2->Draw("SAME");

  l->Draw("SAME");

  canv->SaveAs(Form("Tracks/%s", nomepic));

  delete canv;
  delete line1;
  delete line2;
}

/**
 * Saves a pictorial representation showing the directionality of the track.
 *
 * @param nomepic The filename for the saved picture.
 */
void Analyzer::SavePicDir(const char* nomepic){

  TCanvas* canv = new TCanvas("canv","canv",3500,1500);
  canv->Divide(3,1);

  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0, Form("%f",fPhiDir/TMath::Pi()*180));

  canv->cd(1);
  fTrack->Draw("COLZ");
  fBarPlot->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");

  canv->cd(2);
  fTrackTail->Draw("COLZ");
  fIPPlot->Draw("SAMEP");

  canv->cd(3);
  fScaledTrack->Draw("COLZ");
  fIPPlot->Draw("SAMEP");
  fLineDirection->Draw("SAME");
  canv->cd(3)->SetLogz();
  l->Draw("SAME");

  canv->SaveAs(Form("Tracks/%s",nomepic));

  delete canv;
}

/**
 * Plots and saves a representation of the directionality process.
 *
 * @param nomepic The filename for the saved picture.
 */
void Analyzer::PlotandSavetoFileDirectionalFull(const char* nomepic){

  TCanvas* canv = new TCanvas("canv","canv",800,800);
  canv->Divide(2,2);

  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0, Form("%f",fPhiDir/TMath::Pi()*180));

  TLegend* l2 = new TLegend();
  l2->AddEntry((TObject*)0,Form("NPx=%i",fnpixelx*fnpixely)); 	//assuming the track is a rectangle... doubtful
  l2->AddEntry((TObject*)0,Form("Int=%f",fintegral));
  l2->AddEntry((TObject*)0,Form("Dens=%f",fintegral/(fnpixelx*fnpixely)));
  l2->AddEntry((TObject*)0,Form("Skew=%f",fSkewOnLine));
  l2->AddEntry((TObject*)0,Form("SkewNorm=%f",fSkewOnLine/fintegral));

  l2->AddEntry((TObject*)0,Form("RMS=%f",fRMSOnMainAxis));
  l2->AddEntry((TObject*)0,Form("RMSNorm=%f",fRMSOnMainAxis/fintegral));

  canv->cd(1);
  fTrack->Draw("COLZ");
  fBarPlot->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");
  l2->Draw("SAME");

  // Have do to this workaround to name the first pad with a different name.
  TPad *padtitle = new TPad("padtitle", "padtitle",0.2,0.90,0.8,0.99);
  padtitle->Draw("SAME");
  padtitle->cd();
  // padtitle->SetFillStyle(1);
  padtitle->SetFillColor(kWhite);
  auto tex = new TLatex(0.5,0.5,"Original");
  tex->SetTextAlign(22);
  tex->SetTextSize(0.5);
  tex->Draw();

  canv->cd(2);
  fTrackTail->SetTitle("Track tail + IP");
  fTrackTail->Draw("COLZ");
  fIPPlot->Draw("SAMEP");

  canv->cd(3);
  fScaledTrack->SetTitle("Impact Point");
  fScaledTrack->Draw("COLZ");
  fIPPlot->Draw("SAMEP");
  fLineDirection->Draw("SAME");
  canv->cd(3)->SetLogz();
  l->Draw("SAME");

  canv->cd(4);
  fTrack->SetTitle("Final Directionality Angle");
  fTrack->Draw("COLZ");
  fIPPlot->Draw("SAMEP");
  fLineDirection->SetLineColor(kBlack);
  fLineDirection->SetLineWidth(2);
  fLineDirection->SetLineStyle(9);
  fLineDirection->Draw("SAME");

  // canv->SaveAs(Form("Tracks/%s.png",nomepic));
  canv->SetName(nomepic);
  canv->Write();
  canv->DrawClone();
  delete canv;
}

/**
 * Saves the current state to a ROOT file.
 *
 * @param nomefile The filename for the saved ROOT file.
 */
void Analyzer::SaveRootFile(const char* nomefile)
{  
  TFile* f = new TFile(Form("Tracks/%s",nomefile),"recreate");
  f->cd();

  fTrack->Write();
  fTrackTail->Write();
  fScaledTrack->Write();
  fBarPlot->Write();
  fIPPlot->Write();
  fLineDirection->Write();

  f->Save();
  f->Close();
}

/**
 * Calculates and updates the barycenter of the histogram associated with the track.
 * The barycenter is calculated using the weighted average positions of the histogram bins,
 * where the weights are the bin contents (assumed to represent charge or hit density).
 */
void Analyzer::Barycenter()
{
  double Xb=0;
  double Yb=0;
  double Z=0;
  double ChargeTot=0;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fTrack->GetBinContent(i,j);
      
      if(Z>0) {

        Xb+=(Z*fTrack->GetXaxis()->GetBinCenter(i));
        Yb+=(Z*fTrack->GetYaxis()->GetBinCenter(j));
        ChargeTot+=Z;
      }
    }
  }

  Xb/=ChargeTot;
  Yb/=ChargeTot;

  fXbar=Xb;
  fYbar=Yb;

  return;
}

/**
 * Calculates and updates the barycenter of the histogram associated with the track.
 * The barycenter is calculated using the weighted average positions of the histogram bins,
 * where the weights are the bin contents (assumed to represent charge or hit density).
 */
void Analyzer::Barycenter(TH2F* Tr,double *X,double *Y)
{
  double Xb=0;
  double Yb=0;
  double Z=0;
  double ChargeTot=0;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=Tr->GetBinContent(i,j);
      
      if(Z>0) {

        Xb+=(Z*Tr->GetXaxis()->GetBinCenter(i));
        Yb+=(Z*Tr->GetYaxis()->GetBinCenter(j));
        ChargeTot+=Z;
      }
    }
  }

  Xb/=ChargeTot;
  Yb/=ChargeTot;

  *X=Xb;
  *Y=Yb;

  return;
}

/**
 * Computes the angle of the line passing through the barycenter of a 2D histogram
 * that maximizes the root mean square (RMS) deviation of the data points from the line.
 *
 * The function first calculates the barycenter of the track (fXbar, fYbar) and initializes
 * a plot for it. It then iterates over all bins of the histogram `fTrack`, calculating
 * two sums, `Sum1` and `Sum2`, which are used to determine the angle `Phi` of the line
 * passing through the barycenter. The RMS deviation of the data points is calculated for
 * this line and its perpendicular line. The function returns the angle of the line with
 * the larger RMS deviation.
 *
 * @return The angle (in radians) of the line that maximizes the RMS deviation of the data points.
 */
double Analyzer::AngleLineMaxRMS()
{
  double Sum1=0;
  double Sum2=0;
  double Z=0.;
  double Phi;
  double RmsAng;
  double RmsAngPerp;
  Barycenter(fTrack,&fXbar,&fYbar);

  fBarPlot=new TGraph();
  fBarPlot->SetName("BarycenterPlot");
  fBarPlot->SetPoint(0,fXbar,fYbar);
  fBarPlot->SetMarkerStyle(8);

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fTrack->GetBinContent(i,j);
      if(Z>0) {

        Sum1+= Z*(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*(fTrack->GetYaxis()->GetBinCenter(j)-fYbar);
        Sum2+= Z*( (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*(fTrack->GetYaxis()->GetBinCenter(j)-fYbar) - (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)  );
      }
    }
  }

  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);

  RmsAng=RMSOnLine(Phi);
  RmsAngPerp=RMSOnLine(Phi+TMath::Pi()/2);

  if( RmsAng > RmsAngPerp ) {

    fRMSOnMainAxis=RmsAng;
    return Phi;
  }
  else {

    fRMSOnMainAxis=RmsAngPerp;
    if(Phi+TMath::Pi()/2>TMath::Pi()/2)     return Phi+TMath::Pi()/2-TMath::Pi();
    else     return Phi+TMath::Pi()/2;
  }
}


/**
 * Helper method called by AngleLineMaxRMS to calculate the RMS value along a line at a given angle.
 *
 * @param Phi The angle (in radians) at which to calculate the RMS value.
 * @return The RMS value calculated along the line at the specified angle Phi.
 */
double Analyzer::RMSOnLine(double Phi)
{
  double RMS=0;
  double ChargeTot=0;
  double Z=0.;


  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fTrack->GetBinContent(i,j);
      
      if(Z!=0) {

        RMS+= Z*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(Phi) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(Phi) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(Phi) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(Phi) );
        ChargeTot+=Z;
      }
    }
  }

  return RMS;
}

/**
 * Calculates the skewness along the main axis of the track histogram.
 *
 * @return The calculated skewness value, providing insight into the asymmetry of the hit distribution
 * along the main axis of the track.
 */
double Analyzer::SkewOnMainAxis(){

  Float_t Skew=0;
  Float_t ChargeTot=0;
  Float_t Z;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fTrack->GetBinContent(i,j);
      
      if(Z>0) {

        Skew+= Z*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) );
        ChargeTot+=Z;
      }
    }
  }

  fSkewOnLine=Skew;

  return Skew;
}

/**
 * Removes isolated noise hits from the track histogram.
 * A hit is considered isolated and is removed if it has no neighboring hits within a one-bin radius.
 *
 * @param a The threshold for neighboring hit energy sum. Hits with neighboring energy sum less than or equal to 'a' are removed.
 */
void Analyzer::RemoveNoise(double EnSum){
  std::vector<Int_t> ToRemoveXBin;
  std::vector<Int_t> ToRemoveYBin;

  Float_t ETemp;

  for(int i = 1; i < fnpixelx - 1; i++){
    
    for(int j = 1; j < fnpixely - 1; j++){
      
      if(fTrack->GetBinContent(i, j) > 0) {

        ETemp = 0;
        if(fTrack->GetBinContent(i + 1, j) >= 0) ETemp += fTrack->GetBinContent(i + 1, j);
        if(fTrack->GetBinContent(i - 1, j) >= 0) ETemp += fTrack->GetBinContent(i - 1, j);
        if(fTrack->GetBinContent(i, j + 1) >= 0) ETemp += fTrack->GetBinContent(i, j + 1);
        if(fTrack->GetBinContent(i, j - 1) >= 0) ETemp += fTrack->GetBinContent(i, j - 1);

        if(ETemp <= EnSum){
          ToRemoveXBin.push_back(i);
          ToRemoveYBin.push_back(j);
        }
      }
    }
  }

  for(int i = 0; i < ToRemoveXBin.size(); i++){
    fTrack->SetBinContent(ToRemoveXBin[i], ToRemoveYBin[i], 0);
  }
}

/**
 * Applies a threshold to the track histogram, setting bin contents below a specified value to zero.
 * This method is useful for further noise reduction and data cleaning.
 */
void Analyzer::ApplyThr(double EnThr){
  Int_t XBinMin=fTrack->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+fTrack->GetXaxis()->GetNbins();

  Int_t YBinMin=fTrack->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+fTrack->GetYaxis()->GetNbins();

  Float_t z;

  for(int i=XBinMin; i<XBinMax;i++){
    
    for(int j=YBinMin;j<YBinMax;j++){
      
      z=fTrack->GetBinContent(i,j);
      
      if(z>0 && z<=EnThr){
	      fTrack->SetBinContent(i,j,0);
      }
    }
  }
}

/**
 * Identifies the impact point of the track by analyzing the tail of the distribution.
 *
 * @param nometh2 Name for the histogram representing the track's tail.
 * This method updates the internal state related to the impact point coordinates.
 */
void Analyzer::ImpactPoint(const char* nometh2){

  fTrackTail = new TH2F(nometh2,nometh2,fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);

  //fTrackTail->Rebin2D(2,2);

  Float_t rminN=0.;
  //Float_t rminN=0.7;
  Double_t X,Y,Z;
  Int_t NSelPoints=0;
  Float_t PointSkew,PointDistCm;

  std::vector<float> XIntPointPrev;
  std::vector<float> YIntPointPrev;

// Original version
  // do{
  //   rminN+=0.03;
  //   NSelPoints=0;

  // first point is the Barycenter from PCA only
  Barycenter(fTrack,&fXIPPrev,&fYIPPrev);
  XIntPointPrev.push_back(fXIPPrev);
  YIntPointPrev.push_back(fYIPPrev);

  // Updated version for speed
  do{
    if(NSelPoints-fNPIP>6*fNPIP){
      rminN+=2.0;
    } else {
      //rminN+=2.0;
      rminN+=.5;
    }

    //Original version 
    //rminN+=0.03;

    NSelPoints=0;

    fTrackTail->Reset();

    for(int j=0;j<fTrack->GetXaxis()->GetNbins();j++){
      for(int l=0;l<fTrack->GetYaxis()->GetNbins();l++){

        X=fTrack->GetXaxis()->GetBinCenter(j);
        Y=fTrack->GetYaxis()->GetBinCenter(l);
        Z=fTrack->GetBinContent(j,l);

        PointSkew=GetPointSkew(X,Y);
        PointDistCm=PDistCm(X,Y);

        if(PointSkew>0 && PointDistCm> rminN && Z>0){
          fTrackTail->SetBinContent(fTrackTail->GetXaxis()->FindBin(X),fTrackTail->GetYaxis()->FindBin(Y),Z);
          NSelPoints++;
        }//chiudo if selection

      }//chiuso for l
    }//chiudo for j (fill histos)


    Barycenter(fTrackTail,&fXIPPrev,&fYIPPrev);

    XIntPointPrev.push_back(fXIPPrev);
    YIntPointPrev.push_back(fYIPPrev);


  }while(NSelPoints>fNPIP);

  Barycenter(fTrackTail,&fXIP,&fYIP);

  int PrevIndex=(int)(XIntPointPrev.size()/2);
  if ((int)XIntPointPrev.size()<=2){
    PrevIndex=0;
  }

  fXIPPrev=XIntPointPrev[PrevIndex];
  fYIPPrev=YIntPointPrev[PrevIndex];

  fIPPlot = new TGraph();
  fIPPlot->SetName("IPPLot");
  fIPPlot->SetPoint(0,fXIP,fYIP);
  fIPPlot->SetMarkerStyle(8);

}

/**
 * Generates a scaled version of the track histogram, emphasizing the region around the impact point.
 *
 * @param nometh2 Name for the scaled histogram.
 * The scaling is performed based on the distance from the impact point, with further points being
 * de-emphasized.
 */
void Analyzer::ScaledTrack(const char* nometh2){

  fScaledTrack = new TH2F(nometh2,nometh2,fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
  //fScaledTrack->Rebin2D(2,2);

  double X,Y,Z;

  for(int j=0;j<fnpixelx;j++){
    for(int l=0;l<fnpixely;l++){
      
      if(fTrack->GetBinContent(j,l)>0){
      
        X=fTrack->GetXaxis()->GetBinCenter(j);
        Y=fTrack->GetYaxis()->GetBinCenter(l);
        Z=fTrack->GetBinContent(j,l);

        fScaledTrack->SetBinContent(fScaledTrack->GetXaxis()->FindBin(X),fScaledTrack->GetYaxis()->FindBin(Y),Z*exp(-( sqrt( (X-fXIP)*(X-fXIP)+(Y-fYIP)*(Y-fYIP) ) )/fwScal ) );
      }
    }//chiudo l
  }//chiudo for j
}

/**
 * Calculates the direction of the track based on the scaled histogram.
 * This method updates the internal state related to the track's direction.
 */
void Analyzer::Direction()
{
  double Sum1=0;
  double Sum2=0;
  double Z=0.;
  double Phi;
  double RmsAng;
  double RmsAngPerp;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fScaledTrack->GetBinContent(i,j);
      if(Z>0) {

        Sum1+= Z*(fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)*(fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP);
        Sum2+= Z*( (fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP)*(fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP) - (fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)*(fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)  );
      }
    }
  }

  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);

  RmsAng=RMSOnLine(Phi);
  RmsAngPerp=RMSOnLine(Phi+TMath::Pi()/2);

  if( RmsAng > RmsAngPerp ) {

    fPhiDir=Phi;
  }
  else {

    fRMSOnMainAxis=RmsAngPerp;
    if(Phi+TMath::Pi()/2>TMath::Pi()/2)     fPhiDir = Phi+TMath::Pi()/2-TMath::Pi();
    else     fPhiDir= Phi+TMath::Pi()/2;
  }

}

/**
 * Corrects the track's direction angle to ensure consistency in the representation.
 * This method adjusts the direction angle based on geometric considerations.
 */
void Analyzer::ImprCorrectAngle(){

  Float_t AnglePerp = fPhiDir+TMath::Pi()/2;

  Float_t qIP= fYIP-TMath::Tan(AnglePerp)*fXIP;
  Float_t qPIP= fYIPPrev-TMath::Tan(AnglePerp)*fXIPPrev;

  if(fPhiDir>0){
    if(qIP<qPIP) {
      return;
    } else {
      fPhiDir= fPhiDir-TMath::Pi();
    }//chiudo else qIP
  } 
  else {
    if(qIP<qPIP){
      fPhiDir= fPhiDir+TMath::Pi();
    } else {
      return;
    }//chiudo else qIP
  }//chiudo esle angolo

}



/**
 * Identifies the edges of the track along the main axis, useful for profiling and analysis.
 *
 * @param Xl, Yl References to store the coordinates of the leftmost edge point.
 * @param Xr, Yr References to store the coordinates of the rightmost edge point.
 * @param slope The slope of the main axis, used to determine the direction for edge detection.
 */
void Analyzer::Edges(double &Xl, double &Yl, double &Xr, double &Yr, double slope) 
{
  double Xp, Yp, Zp;
  double dist, tempdist_r=0, tempdist_l=0;
  double ii=0, jj=0;

  Barycenter(fTrack, &fXbar, &fYbar);

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Zp=fTrack->GetBinContent(i,j);
      if(Zp!=0) {

        ii=fTrack->GetXaxis()->GetBinCenter(i);
        jj=fTrack->GetYaxis()->GetBinCenter(j);
        Xp = (1./(1+pow(slope,2)))*(ii+fXbar*pow(slope,2)+slope*(jj-fYbar));
        Yp = fYbar+(slope/(1+pow(slope,2)))*(ii-fXbar+slope*(jj-fYbar));
        dist = sqrt(pow((Xp-fXbar),2)+pow((Yp-fYbar),2));
        
        if(dist>tempdist_r && Xp>fXbar) {

          tempdist_r = dist;
          Xr = Xp; Yr = Yp;
        }
        else if(dist>tempdist_l && Xp<fXbar) {

          tempdist_l = dist;
          Xl = Xp; Yl = Yp;
        }
      }
    }
  }

  return;
}


/**
 * Generates a profile histogram along the main axis or perpendicular to it.
 *
 * @param longitudinal If true, profiles along the main axis; otherwise, profiles perpendicularly.
 * @param x1, x2 Range limits for profiling. Defaults cover the entire range.
 * @return A pointer to the generated TH1D profile histogram.
 */
TH1D* Analyzer::FillProfile(bool longitudinal, float x1, float x2)
{

  double xl,yl,xr,yr;
  double slope;
  double ii=0, jj=0;

  if(longitudinal) slope =  tan(AngleLineMaxRMS());
  else slope = tan(TMath::Pi()/2.+AngleLineMaxRMS());

  Edges(xl,yl,xr,yr,slope);
  int binmax = (int)sqrt(pow((xl-xr),2)+pow((yl-yr),2));
  TH1D* TrackProfile=new TH1D("TrackProf","TrackProf",binmax+2,0,binmax+2);

  double Xp, Yp, Zp;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Zp=fTrack->GetBinContent(i,j);
      
      if(Zp!=0 ) {

        ii=fTrack->GetXaxis()->GetBinCenter(i);
        jj=fTrack->GetYaxis()->GetBinCenter(j);
        
        if(ii>x1 && ii<x2) {

          Xp = (1/(1+pow(slope,2)))*(ii+fXbar*pow(slope,2)+slope*(jj-fYbar));
          Yp = fYbar+(slope/(1+pow(slope,2)))*(ii-fXbar+slope*(jj-fYbar));
          TrackProfile->Fill(sqrt(pow((Xp-xl),2)+pow((Yp-yl),2)),Zp);
        }
      }
    }
  }

  return TrackProfile;
}

/**
 * Trims a profile histogram to focus on the main region of interest around the peak.
 *
 * @param profile The original profile histogram to trim.
 * @param height Fraction of the peak height to use as a cutoff for trimming.
 * @return A pointer to the trimmed TH1D profile histogram.
 */
TH1D* Analyzer::CutProfile(TH1D* profile, double height)
{
//height: cut bins with content lower than height*maximum intensity in the profile  (default is 0.25%)

  int binmin = profile->FindFirstBinAbove(height*profile->GetMaximum(),1);
  int binmax = profile->FindLastBinAbove(height*profile->GetMaximum(),1);
  int bins_cut = 0;

  std::string title = profile->GetTitle();
  title += "_cut";
  TH1D* profile_cut = new TH1D(title.c_str(),title.c_str(),binmax-binmin,0,binmax-binmin);

  for(int bins=binmin; bins<binmax; bins++) {

    profile_cut->SetBinContent(bins_cut,profile->GetBinContent(bins));
    bins_cut++;
  }

  return profile_cut;

}

/**
 * Generates a profile histogram along the X-axis of the track histogram.
 *
 * @return A pointer to the generated TH1D X-axis profile histogram.
 */
TH1D* Analyzer::FillProfileX()
{
  double Zp;

  TH1D* TrackProfileX=new TH1D("TrackProfX","TrackProfX",fnpixelx,fminx,fmaxx);

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Zp=fTrack->GetBinContent(i,j);
      
      if(Zp!=0) {

        TrackProfileX->Fill(fTrack->GetXaxis()->GetBinCenter(i),Zp);
      }
    }
  }

  return TrackProfileX;

}

/**
 * Generates a profile histogram along the Y-axis of the track histogram.
 *
 * @return A pointer to the generated TH1D Y-axis profile histogram.
 */
TH1D* Analyzer::FillProfileY()
{
  double Zp;

  TH1D* TrackProfileY=new TH1D("TrackProfY","TrackProfY",fnpixely,fminy,fmaxy);

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Zp=fTrack->GetBinContent(i,j);
      if(Zp!=0) {

        TrackProfileY->Fill(fTrack->GetYaxis()->GetBinCenter(j),Zp);
      }
    }
  }

  return TrackProfileY;

}

/**
 * Calculates the principal component analysis (PCA) angle of the track histogram.
 *
 * @param ang Reference to a double to store the calculated angle.
 * The PCA angle provides another measure of the track's main direction.
 */
void Analyzer::AnglePCA(double &ang)
{
  TMatrixD M(2,2);
  TMatrixD V(2,2);
  double x,y,x2,y2,xy;
  double Z;

  for(int i=0; i<fnpixelx; i++) {

    for(int j=0; j<fnpixely; j++) {

      Z=fTrack->GetBinContent(i,j);
      
      if(Z!=0){
      
        x=(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*Z;
        y=(fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*Z;
        x2 += x*x;
        y2 += y*y;
        xy += x*y;
      }
    }
  }

  M(0,0) = x2;
  M(0,1) = M(1,0) = xy;
  M(1,1) = y2;

  TDecompSVD SVDmatrix(M);
  if(SVDmatrix.Decompose()) V = SVDmatrix.GetV();

  ang = TMath::ATan(V(1,0)/V(0,0));

}

/**
 * Finds peaks within a one-dimensional histogram using the TSpectrum class.
 *
 * @param h Pointer to the TH1D histogram to analyze.
 * @param foundpeaks Reference to a vector to store the found peaks (position, sigma).
 */
void Analyzer::FindNPeaks(TH1D* h, std::vector<std::pair<double,double>> &foundpeaks)
{

  TSpectrum* s = new TSpectrum();

  int npeaks;
  double* peaksPos(nullptr);

  std::vector< std::pair<double,double> > peaks;
  std::vector<double> peaks_tocompare;

  for(int i=2; i<16; i=i+1){ //scan for different sigma

    npeaks = s->Search(h,i,"nobackground",0.1); //number of peaks with current sigma (npeaks2=number of peaks with previous sigma)
    peaksPos = s->GetPositionX(); //positions of peaks with current sigma (peaks2=positions of peaks with previous sigma)
    if (npeaks == 0){continue;} //if no peaks are found, go to next sigma
    else{ //if peaks are found
      //if any of the peaks found in the previous iteration is equal (within one sigma) to the current one, ignore it; otherwise, save the position of the new peak and increase the total number of peaks
      for(int j=0; j<npeaks; j++){ //loop over number of new peaks

        if(!peaks_tocompare.empty()){

          auto it = std::find_if(peaks_tocompare.begin(), peaks_tocompare.end(), [&](double p){ return (p>(peaksPos[j]-(double)i) && p<(peaksPos[j]+(double)i)); });
          if(it != peaks_tocompare.end()){ //if it's already stored
            continue;
          }
	        else{ //it's a new peak, save it
            peaks.push_back( std::make_pair(peaksPos[j],(double)i) );
            continue;
          }
        }
        peaks.push_back(std::make_pair(peaksPos[j],(double)i));
      } //end loop sui picchi trovati
    } //end if peaks were found

    peaks_tocompare.clear();
    for(int k=0; k<npeaks; k++){peaks_tocompare.push_back(peaksPos[k]);}

  } //end scan on different sigma

foundpeaks = peaks;

delete s;

return;

}

/**
 * Identifies the coordinates of the maximum intensity pixel in the track histogram.
 *
 * @param xpeak, ypeak References to store the coordinates of the peak in the original histogram.
 * @param xpeak_rebin, ypeak_rebin References to store the coordinates of the peak in a rebinned version of the histogram.
 */
void Analyzer::FindPeak(double &xpeak, double &ypeak, double &xpeak_rebin, double &ypeak_rebin) 
{

  int maxbin = fTrack->GetMaximumBin();
  int x,y,z;
  fTrack->GetBinXYZ(maxbin, x, y, z);

  xpeak=fTrack->GetXaxis()->GetBinCenter(x);
  ypeak=fTrack->GetYaxis()->GetBinCenter(y);

  TH2F* TrackRebin = (TH2F*)fTrack->Clone();
  //TrackRebin->Rebin2D(2,2);
  maxbin = TrackRebin->GetMaximumBin();
  TrackRebin->GetBinXYZ(maxbin, x, y, z);
  xpeak_rebin=TrackRebin->GetXaxis()->GetBinCenter(x);
  ypeak_rebin=TrackRebin->GetYaxis()->GetBinCenter(y);

}

//
/*TF1* Analyzer::LeastSquareLine()
{
  double sum1=0, sum2=0, sum3=0, sum4=0, sum5=0, a=0 , b=0;
  double Z=0;

  for(int i=1;i<fnpixelx;i++) {

    for(int j=1;j<fnpixely;j++) {

      Z=fTrack->GetBinContent(i,j);
      
      if(Z!=0) {

        sum1+= Z*(fTrack->GetXaxis()->GetBinCenter(i))*(fTrack->GetXaxis()->GetBinCenter(i));
        sum2+= Z*(fTrack->GetYaxis()->GetBinCenter(j));
        sum3+= Z*(fTrack->GetXaxis()->GetBinCenter(i));
        sum4+= Z*(fTrack->GetXaxis()->GetBinCenter(i))*(fTrack->GetYaxis()->GetBinCenter(j));
        sum5+= Z;
      }
    }
  }

  a = (sum1*sum2-sum3*sum4)/(sum5*sum1-(sum3*sum3));
  b = (sum5*sum4-sum3*sum2)/(sum5*sum1-(sum3*sum3));

  TF1* line = new TF1("leastsquare","[0]*x+[1]",1000,2000);
  line->SetParameters(b,a);

  return line;

}*/


/**
 * Executes an external Python script for additional data analysis, potentially including
 * discrimination variables calculation, plotting, and textual output.
 *
 * @param pyvers String indicating the Python version to use.
 * @param inputfile The path to the input file for the script.
 * @param outfolder The directory where the script's output should be saved.
 * @param entries The number of entries from the input file that the script should process.
 * @param plot Boolean indicating whether the script should generate plots.
 * @param text Boolean indicating whether the script should generate textual output.
 * @return An integer indicating the success (0) or failure (non-zero) of the operation.
 */
int Analyzer::Execute_Atul_script(std::string pyvers, std::string inputfile, std::string outfolder, int entries, bool plot, bool text ) const
{
	int i=system("python"+pyvers+" discriminating_vars_BaData.py -I "+inputfile+ " -E "+entries+ " -P "+plot+" -T "+text+"  -O "+ outfolder);
	return i;		//0 if the command was successful, 1 if not
}
