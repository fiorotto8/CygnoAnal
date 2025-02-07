import uproot
import numpy as np
import ROOT
import argparse

def hist(data, x_name, channels=100, linecolor=4, linewidth=4,write=True,save=False):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, 0.99*np.min(array), 1.01*np.max(array))
    # Use numpy vectorization to fill histogram
    for x in array:
        hist.Fill(x)
    # Set visual attributes and axis titles
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    # Set maximum digits on axes to manage display
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetMaxDigits(3)
    if write: hist.Write()
    return hist
def grapherr(x,y,ex,ey,x_string, y_string,name=None, color=4, markerstyle=22, markersize=2,write=True):
    plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") , np.array(   ex   ,dtype="d"),np.array( ey   ,dtype="d"))
    if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
    else: plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)#blue
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write==True: plot.Write()
    return plot

parser = argparse.ArgumentParser(description="Fit modulation data from a ROOT file.")
parser.add_argument("--energy", type=float, default=17.0, help="Energy value for the photoelectron direction distribution.")
parser.add_argument("--file", type=str, default="Analysis_reco_Pol_Fusion.root", help="Path to the ROOT file.")
args = parser.parse_args()

name = f"{args.energy}keV photoelectron"

file = uproot.open(args.file)
# Access the TTree named "track_info"
track_info_tree = file["track_info"]
# Extract arrays from the TTree
directions = track_info_tree["DegreeDIR"].array()
energies = track_info_tree["Integral"].array()
pu_candidate = track_info_tree["PileUpCandidate"].array()

# Filter out entries where puCandidate == 1
mask = pu_candidate != 1
directions = directions[mask]
energies = energies[mask]

#! Define histo for energy
histoEnergy=hist(energies,name+" Energy spectrum",channels=500,write=False)
# Define the fitting function for energy
fit_func_energy = ROOT.TF1("fit_func_energy", "gaus", histoEnergy.GetXaxis().GetXmin(), histoEnergy.GetXaxis().GetXmax())
fit_func_energy.SetParameters(histoEnergy.GetBinContent(histoEnergy.GetMaximumBin()), histoEnergy.GetMean(), histoEnergy.GetRMS())
histoEnergy.Fit(fit_func_energy, "R")

fit_func_energy = ROOT.TF1("fit_func_energy", "gaus(0) + pol1(3)", histoEnergy.GetXaxis().GetXmin(), histoEnergy.GetXaxis().GetXmax())
fit_func_energy.SetParameters(histoEnergy.GetBinContent(histoEnergy.GetMaximumBin())-histoEnergy.GetBinContent(5), histoEnergy.GetMean(), histoEnergy.GetRMS(),histoEnergy.GetBinContent(5),0) 
print(histoEnergy.GetBinContent(histoEnergy.GetMaximumBin())-histoEnergy.GetBinContent(5), histoEnergy.GetMean(), histoEnergy.GetRMS(),histoEnergy.GetBinContent(5),0)
histoEnergy.Fit(fit_func_energy, "R")

EnResolution=fit_func_energy.GetParameter(2)/fit_func_energy.GetParameter(1)
# Calculate the area under the linear background (pol1) to estimate pile-up
xmin = histoEnergy.GetXaxis().GetXmin()
xmax = histoEnergy.GetXaxis().GetXmax()
slope = fit_func_energy.GetParameter(4)
intercept = fit_func_energy.GetParameter(3)
pileUp = 0.5 * (intercept + (intercept + slope * (xmax - xmin))) * (xmax - xmin)
goodEvents = fit_func_energy.Integral(xmin, xmax)
allEvents = pileUp + goodEvents
pileUpFraction=pileUp/allEvents

# Draw the fit function on the histogram
fit_func_energy.Draw("same")
canvas = ROOT.TCanvas(name+" Energy spectrum", name+" Energy spectrum", 1000, 1000)
histoEnergy.Draw()
# Add statistics box to the canvas
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0000)

# Create TPaveText to display the pile up fraction
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.03)
latex.SetTextColor(4)
#latex.DrawLatex(0.15, 0.35, f"Good Events: {goodEvents:.0f}")
#latex.DrawLatex(0.15, 0.3, f"Pile Up Events: {pileUp:.0f}")
#latex.DrawLatex(0.15, 0.25, f"All events: {histoEnergy.GetEntries():.0f}")
latex.DrawLatex(0.15, 0.25, f"Resolution(%): {100*EnResolution:.1f}")
latex.DrawLatex(0.15, 0.2, f"Background Fraction: {pileUpFraction*100:.2f}%")
latex.DrawLatex(0.15, 0.15, f"Expected Background Fraction: 25%")

canvas.Update()
canvas.SaveAs(f"Energy_{args.energy}.png")

#! define histogram for directions
histo=hist(directions,name,write=False)
# Define the fitting function
fit_func = ROOT.TF1("fit_func", "[0] + [1]*cos((x - [2])*3.14159/180)**2", -180, 180)
fit_func.SetParameters(250, 150, 90)
histo.Fit(fit_func, "RQ")

PolDegree=fit_func.GetParameter(1)/(2*fit_func.GetParameter(0)+fit_func.GetParameter(1))

# Draw the fit function on the histogram
fit_func.Draw("same")
canvas1 = ROOT.TCanvas(name+" direction distribution", name+" direction distribution", 1000, 1000)
histo.Draw()
# Add statistics box to the canvas
#ROOT.gStyle.SetOptFit(1111)

# Create TPaveText to display the modulation factor
latex.DrawLatex(0.15, 0.22, f"Modulation Factor: {PolDegree/0.75:.3f}")
latex.DrawLatex(0.15, 0.17, f"Purity of the line: 75%")
latex.DrawLatex(0.15, 0.12, "Fitting Function: [0] + [1]*cos^2((x - [2]))")

canvas1.Update()

canvas1.SaveAs(f"ModualtionCurve_{args.energy}.png")


