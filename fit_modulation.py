import uproot
import numpy as np
import ROOT
import argparse
import os

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
parser.add_argument("--draw", action="store_true", help="Flag to draw the histograms.",default=False)
parser.add_argument("--run", help="append to ouput name",default=None)
parser.add_argument("--out_folder", type=str, default="output/", help="Output folder for the plots.")
args = parser.parse_args()

line_purity=0.75
line_purity_err=0.15

name = f"{args.energy}keV photoelectron"

ROOT.gROOT.SetBatch(True)

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

# Get the bin with the largest number of entries
max_bin = histoEnergy.GetMaximumBin()
ymax = histoEnergy.GetBinContent(max_bin)
histoEnergy.GetYaxis().SetRangeUser(0, ymax*1.5)

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
latex.DrawLatex(0.15, 0.85, f"Resolution(%): {100*EnResolution:.1f}")
latex.DrawLatex(0.15, 0.8, f"Background Fraction: {pileUpFraction*100:.2f}%")
latex.DrawLatex(0.15, 0.75, f"Expected Background Fraction: 25%")
latex.DrawLatex(0.15, 0.7, f"chi2/ndf: {fit_func_energy.GetChisquare()/fit_func_energy.GetNDF():.1f}")

canvas.Update()
if args.draw is True: 
    if args.run is None: canvas.SaveAs(f"Energy_{args.energy}.png")
    else: canvas.SaveAs(f"{args.out_folder}/Energy_{args.energy}_run{args.run}.png")

#! define histogram for directions
histo=hist(directions,name,channels=200,write=False)
# Define the fitting function
fit_func = ROOT.TF1("fit_func", "[0] + [1]*cos((x - [2])*3.14159/180)**2", -180, 180)
fit_func.SetParameters(250, 150, 90)
histo.Fit(fit_func, "RQ")

PolDegree=fit_func.GetParameter(1)/(2*fit_func.GetParameter(0)+fit_func.GetParameter(1))/line_purity
PolDegree_err = PolDegree * np.sqrt(
    (fit_func.GetParError(1) / fit_func.GetParameter(1))**2 +
    (fit_func.GetParError(0) / fit_func.GetParameter(0))**2 +
    (line_purity_err / line_purity)**2
)

# Draw the fit function on the histogram
fit_func.Draw("same")
canvas1 = ROOT.TCanvas(name+" direction distribution", name+" direction distribution", 1000, 1000)
histo.Draw()
# Add statistics box to the canvas
#ROOT.gStyle.SetOptFit(1111)

# Get the bin with the largest number of entries
max_bin = histo.GetMaximumBin()
ymax = histo.GetBinContent(max_bin)
histo.GetYaxis().SetRangeUser(0, ymax*1.8)


# Create TPaveText to display the modulation factor
latex.DrawLatex(0.15, 0.85, f"Modulation Factor: {PolDegree:.3f} +/- {PolDegree_err:.3f}")
latex.DrawLatex(0.15, 0.8, f"Purity of the line: 75%")
latex.DrawLatex(0.15, 0.75, "Fitting Function: [0] + [1]*cos^2((x - [2]))")
if fit_func.GetNDF()!=0: latex.DrawLatex(0.15, 0.7, f"chi2/ndf: {fit_func.GetChisquare()/fit_func.GetNDF():.1f}")
latex.DrawLatex(0.15, 0.65, f"Phase (deg): {fit_func.GetParameter(2):.1f} +/- {fit_func.GetParError(2):.1f}")


canvas1.Update()

if args.draw is True: 
    if args.run is None: canvas1.SaveAs(f"ModualtionCurve_{args.energy}.png")
    else: canvas1.SaveAs(f"{args.out_folder}/ModualtionCurve_{args.energy}_run{args.run}.png")

if args.run is None: output_file = f"{args.out_folder}/output.txt"
else: output_file = f"{args.out_folder}/output_{args.run}.txt"
modulation_factor = PolDegree
if fit_func.GetNDF()!=0: reduced_chi2 = fit_func.GetChisquare() / fit_func.GetNDF()
else: reduced_chi2=0

# Check if the file exists, if not create it
if not os.path.exists(output_file):
    with open(output_file, "w") as f:
        pass

# Append the modulation factor and reduced chi2 to the file
with open(output_file, "a") as f:
    f.write(f"{modulation_factor:.3f};{reduced_chi2:.3f}\n")

