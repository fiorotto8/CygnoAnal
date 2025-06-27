import uproot
import numpy as np
import ROOT
import argparse
import os

def hist(data, x_name, channels=100, linecolor=4, linewidth=4,write=True,save=False):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, np.min(array), np.max(array))
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
n_sigma = 2

def fit_gaus_plus_linear(
    data, 
    histo_name="Spectrum", 
    nbins=200, 
    xmin=None, 
    xmax=None, 
    draw=True, 
    output_name=None,
    out_folder=".", 
    run=None
):
    # Create histogram
    if xmin is None: xmin = min(data)
    if xmax is None: xmax = max(data)
    histo = ROOT.TH1F(histo_name, histo_name, nbins, xmin, xmax)
    for value in data:
        histo.Fill(value)

    # Initial Gaussian fit
    fit_init = ROOT.TF1("fit_init", "gaus", xmin, xmax)
    fit_init.SetParameters(histo.GetBinContent(histo.GetMaximumBin()), histo.GetMean(), histo.GetRMS())
    histo.Fit(fit_init, "RQ")

    # Final fit: Gaussian + linear background
    fit_func = ROOT.TF1("fit_func", "gaus(0) + pol1(3)", xmin, xmax)

    # Estimate slope from bin difference
    nbins = histo.GetNbinsX()
    bin_low = 0
    bin_high = nbins -0
    x_low = histo.GetXaxis().GetBinCenter(bin_low)
    x_high = histo.GetXaxis().GetBinCenter(bin_high)
    y_low = histo.GetBinContent(bin_low)
    y_high = histo.GetBinContent(bin_high)
    slope_estimate = (y_high - y_low) / (x_high - x_low)

    fit_func.SetParameters(
        histo.GetBinContent(histo.GetMaximumBin()) - histo.GetBinContent(bin_low),
        histo.GetMean(),
        histo.GetRMS(),
        histo.GetBinContent(bin_low),
        slope_estimate
    )

    """
    print("Fit Parameter Estimates:")
    print("  Amplitude:", histo.GetBinContent(histo.GetMaximumBin()) - histo.GetBinContent(bin_low))
    print("  Mean:", histo.GetMean())
    print("  Sigma:", histo.GetRMS())
    print("  Intercept:", histo.GetBinContent(bin_low))
    print("  Slope:", slope_estimate)
    """
    
    histo.Fit(fit_func, "RQ")

    # Energy resolution
    resolution = fit_func.GetParameter(2) / fit_func.GetParameter(1)

    # Background estimation
    xlow_bg = histo.GetXaxis().GetXmin() + bin_low
    xhigh_bg = histo.GetXaxis().GetXmax() - (nbins - bin_high)
    slope = fit_func.GetParameter(4)
    intercept = fit_func.GetParameter(3)
    background_events = 0.5 * (intercept + intercept + slope * (xhigh_bg - xlow_bg)) * (xhigh_bg - xlow_bg)

    # Signal estimation
    signal_func = ROOT.TF1("signal_func", "gaus", xlow_bg, xhigh_bg)
    for i in range(3):
        signal_func.SetParameter(i, fit_func.GetParameter(i))
    signal_events = signal_func.Integral(xlow_bg, xhigh_bg)
    all_events = fit_func.Integral(xlow_bg, xhigh_bg)
    background_fraction = background_events / all_events if all_events > 0 else 0

    # Plotting
    max_bin = histo.GetMaximumBin()
    ymax = histo.GetBinContent(max_bin)
    histo.GetYaxis().SetRangeUser(0, ymax * 1.5)

    canvas = ROOT.TCanvas(histo_name, histo_name, 1000, 1000)
    histo.Draw()
    fit_func.Draw("same")

    # Background line
    background_func = ROOT.TF1("background_func", "pol1", xlow_bg, xhigh_bg)
    background_func.SetParameters(intercept, slope)
    background_func.SetLineColor(ROOT.kRed)
    background_func.SetLineStyle(2)
    background_func.Draw("same")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    # Display results
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextColor(4)
    param_names = ["Amplitude", "Mean", "Sigma", "Background Intercept", "Background Slope"]
    for i in range(fit_func.GetNpar()):
        param = fit_func.GetParameter(i)
        error = fit_func.GetParError(i)
        name = param_names[i] if i < len(param_names) else f"Param[{i}]"
        latex.DrawLatex(0.15, 0.85 - i*0.05, f"{name}: {param:.3f} \u00B1 {error:.3f}")
    latex.DrawLatex(0.15, 0.85 - fit_func.GetNpar()*0.05, f"chi2/ndf: {fit_func.GetChisquare() / fit_func.GetNDF():.1f}")

    canvas.Update()

    # Print the range of n sigma from the Gaussian mean
    print(f"Fitting results for {histo_name}:")
    mean = fit_func.GetParameter(1)
    sigma = fit_func.GetParameter(2)
    lower = mean - n_sigma * abs(sigma)
    upper = mean + n_sigma * abs(sigma)
    print(f"{n_sigma} sigma range: [{lower:.3f}, {upper:.3f}]")

    if draw:
        if output_name is None:
            output_name = f"{histo_name.replace(' ', '_')}.png"
        if run is not None:
            output_name = f"{out_folder}/{histo_name.replace(' ', '_')}_run{run}.png"
        canvas.SaveAs(output_name)

    return {
        "resolution": resolution,
        "background_fraction": background_fraction,
        "fit": fit_func,
        "signal_func": signal_func,
        "background_func": background_func,
        "canvas": canvas,
        "histogram": histo
    }

parser = argparse.ArgumentParser(description="Fit modulation data from a ROOT file.")
parser.add_argument("--energy", type=float, default=17.0, help="Energy value for the photoelectron direction distribution.")
parser.add_argument("--file", type=str, default="Analysis_reco_Pol_Fusion.root", help="Path to the ROOT file.")
parser.add_argument("--draw", action="store_true", help="Flag to draw the histograms.",default=False)
parser.add_argument("--run", help="append to ouput name",default=None)
parser.add_argument("--out_folder", type=str, default="output/", help="Output folder for the plots.")
parser.add_argument("--pileup_cut", action="store_true", default=False, help="Enable the cut on pileup candidates.")
args = parser.parse_args()

# Create the output folder if it does not exist
if not os.path.exists(args.out_folder):
    os.makedirs(args.out_folder)

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

xmeans = track_info_tree["Xmean"].array()
ymeans = track_info_tree["Ymean"].array()
xIP = track_info_tree["X_ImpactPoint"].array()
yIP = track_info_tree["Y_ImpactPoint"].array()

# Filter out entries where puCandidate == 1
if args.pileup_cut:
    mask = pu_candidate != 1
    directions = directions[mask]
    energies = energies[mask]
    xmeans = xmeans[mask]
    ymeans = ymeans[mask]
    xIP = xIP[mask]
    yIP = yIP[mask]


#! Define histo for energy
histoEnergy=hist(energies,name+" Energy spectrum",channels=200,write=False,)
# Define the fitting function for energy
fit_func_energy = ROOT.TF1("fit_func_energy", "gaus", histoEnergy.GetXaxis().GetXmin(), histoEnergy.GetXaxis().GetXmax())
fit_func_energy.SetParameters(histoEnergy.GetBinContent(histoEnergy.GetMaximumBin()), histoEnergy.GetMean(), histoEnergy.GetRMS())
histoEnergy.Fit(fit_func_energy, "RQ")

fit_func_energy = ROOT.TF1("fit_func_energy", "gaus(0) + pol1(3)", histoEnergy.GetXaxis().GetXmin(), histoEnergy.GetXaxis().GetXmax())
# Estimate the slope as the difference between the content of the last and first bins divided by the range
# Estimate the slope using the range minus the first and last 10 bins
nbins = histoEnergy.GetNbinsX()
xmin = histoEnergy.GetXaxis().GetXmin()
xmax = histoEnergy.GetXaxis().GetXmax()
bin_low = 1
bin_high = nbins - 1
x_low = histoEnergy.GetXaxis().GetBinCenter(bin_low)
x_high = histoEnergy.GetXaxis().GetBinCenter(bin_high)
y_low = histoEnergy.GetBinContent(bin_low)
y_high = histoEnergy.GetBinContent(bin_high)
slope_estimate = (y_high - y_low) / (x_high - x_low)
fit_func_energy.SetParameters(
    histoEnergy.GetBinContent(histoEnergy.GetMaximumBin()) - histoEnergy.GetBinContent(bin_low),
    histoEnergy.GetMean(),
    histoEnergy.GetRMS(),
    histoEnergy.GetBinContent(bin_low),
    slope_estimate
)
# 1. Fit parameters (amplitude, mean, sigma, background intercept, background slope):
print("1. Amplitude estimate:", histoEnergy.GetBinContent(histoEnergy.GetMaximumBin())-histoEnergy.GetBinContent(5))
print("2. Mean value:", histoEnergy.GetMean())
print("3. Sigma (RMS):", histoEnergy.GetRMS())
print("4. Background intercept estimate:", histoEnergy.GetBinContent(bin_low))
print("5. Background slope:", slope_estimate)
histoEnergy.Fit(fit_func_energy, "RQ")

EnResolution=fit_func_energy.GetParameter(2)/fit_func_energy.GetParameter(1)
# Calculate the area under the linear background (pol1) to estimate background events
xmin = histoEnergy.GetXaxis().GetXmin()+bin_low
xmax = histoEnergy.GetXaxis().GetXmax()- bin_high
slope = fit_func_energy.GetParameter(4)
intercept = fit_func_energy.GetParameter(3)
# Integral of background (pol1) over the range
backgroundEvents = 0.5 * (intercept + (intercept + slope * (xmax - xmin))) * (xmax - xmin)
# Integral of signal (gaus) over the range
signal_func = ROOT.TF1("signal_func", "gaus", xmin, xmax)
for i in range(3):
    signal_func.SetParameter(i, fit_func_energy.GetParameter(i))
signalEvents = signal_func.Integral(xmin, xmax)
# All events from the fit (gaus + pol1)
allEvents = fit_func_energy.Integral(xmin, xmax)
# Fraction of background events
pileUpFraction = backgroundEvents / allEvents if allEvents > 0 else 0

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

# Draw background and signal components separately
# Background: pol1 (parameters 3 and 4)
background_func = ROOT.TF1("background_func", "pol1", xmin, xmax)
background_func.SetParameters(fit_func_energy.GetParameter(3), fit_func_energy.GetParameter(4))
background_func.SetLineColor(ROOT.kRed)
background_func.SetLineStyle(2)
background_func.Draw("same")


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

print(f"Fitting results for {name} Energy:")
mean = signal_func.GetParameter(1)
sigma = signal_func.GetParameter(2)
lower = mean - n_sigma * abs(sigma)
upper = mean + n_sigma * abs(sigma)
print(f"{n_sigma} sigma range: [{lower:.3f}, {upper:.3f}]")


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

#! crazy test
""" # This is a test to fit the modulation curve with an additional pol1 background

# Take the previous fit function and parameters, and fit again with an additional pol1 background
fit_func = ROOT.TF1("fit_func", "[0] + [1]*cos((x - [2])*3.14159/180)**2 + pol1(3)", -180, 180)
# Set initial parameters: previous fit + pol1 (set to 0,0)
for i in range(3):
    fit_func.SetParameter(i, fit_func.GetParameter(i))
fit_func.SetParameter(3, 200)  # pol1 intercept
fit_func.SetParameter(4, -100)  # pol1 slope
histo.Fit(fit_func, "RQ")
"""
#!

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

if args.run is None:
    output_file = None
else:
    output_file = f"{args.out_folder}/output_{args.run}.txt"
    if os.path.exists(output_file):
        os.remove(output_file)
    # Create the file so it exists for later appending
    with open(output_file, "w") as f:
        pass


modulation_factor = PolDegree
if fit_func.GetNDF()!=0: reduced_chi2 = fit_func.GetChisquare() / fit_func.GetNDF()
else: reduced_chi2=0

# Fit Xmean and Ymean distributions using fit_gaus_plus_linear
fit_result_xmean = fit_gaus_plus_linear(
    xmeans,
    histo_name=f"{name} Xmean",
    nbins=200,
    draw=args.draw,
    output_name=f"{args.out_folder}/Xmean_{args.energy}_run{args.run}.png" if args.run else f"Xmean_{args.energy}.png",
    out_folder=args.out_folder,
    run=args.run
)

fit_result_ymean = fit_gaus_plus_linear(
    ymeans,
    histo_name=f"{name} Ymean",
    nbins=200,
    draw=args.draw,
    output_name=f"{args.out_folder}/Ymean_{args.energy}_run{args.run}.png" if args.run else f"Ymean_{args.energy}.png",
    out_folder=args.out_folder,
    run=args.run
)

fit_result_xip = fit_gaus_plus_linear(
    xIP,
    histo_name=f"{name} XIP",
    nbins=200,
    draw=args.draw,
    output_name=f"{args.out_folder}/XIP_{args.energy}_run{args.run}.png" if args.run else f"XIP_{args.energy}.png",
    out_folder=args.out_folder,
    run=args.run
)

fit_result_yip = fit_gaus_plus_linear(
    yIP,
    histo_name=f"{name} YIP",
    nbins=200,
    draw=args.draw,
    output_name=f"{args.out_folder}/YIP_{args.energy}_run{args.run}.png" if args.run else f"YIP_{args.energy}.png",
    out_folder=args.out_folder,
    run=args.run
)

#! Using stokes parameters to create a graph with error bars


#compute modulation with stokes parameters
iks = track_info_tree["ik"].array()
qks = track_info_tree["qk"].array()
uks = track_info_tree["uk"].array()

mu=1#to measure the modulation factor we need to put mu=1
q=np.mean(qks)/mu
u=np.mean(uks)/mu
I=sum(iks)
err_q=np.sqrt((1/(I-1))*( (2/mu**2)-q**2 )  )
err_u=np.sqrt((1/(I-1))*( (2/mu**2)-u**2 )  )

Pol_deg_stokes=(1/mu)*(q**2 + u**2)**0.5
Pol_deg_stokes_err = Pol_deg_stokes * np.sqrt(
    (err_q / q)**2 +
    (err_u / u)**2
)
Pol_angle_stokes = 0.5 * np.arctan2(u, q) * 180 / np.pi
Pol_angle_stokes_err = 0.5 * np.sqrt(
    (err_u / u)**2 +
    (err_q / q)**2
) * 180 / np.pi

# Print the results
print(f"Modulation factor from Stokes parameters: {Pol_deg_stokes:.3f} +/- {Pol_deg_stokes_err:.3f}")
print(f"Polarization angle: {Pol_angle_stokes:.3f} +/- {Pol_angle_stokes_err:.3f} degrees")   


if output_file is not None:
    # Check if the file exists, if not create it
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            # Write header
            f.write("modulation_factor_curve,error_curve,chi2_curve,modulation_factor_stokes,error_stokes,polarization_angle,polarization_angle_err,other_params\n")

    # Append the results to the file
    with open(output_file, "a") as f:
        f.write(
            f"{args.run},"
            f"{PolDegree:.3f},"
            f"{PolDegree_err:.3f},"
            f"{reduced_chi2:.3f},"
            f"{Pol_deg_stokes:.3f},"
            f"{Pol_deg_stokes_err:.3f},"
            f"{Pol_angle_stokes:.3f},"
            f"{Pol_angle_stokes_err:.3f},"
            f"\n"
        )
