import Analyzer
import argparse
from glob import glob
import os
import uproot
import ROOT
from tqdm import tqdm
import numpy as np
import re

# Hard coded value for pixel size in micrometers
pixel_size = 64.5  # um
# Function to calculate distance based on hole number
def distance(hole_num):  # in mm
    return 5.4 + 10.6 * hole_num

err_V=0.01
err_LoverE=0.1

def grapherr(x, y, ex, ey, x_string, y_string, name=None, color=4, markerstyle=22, markersize=2, write=True):
    plot = ROOT.TGraphErrors(len(x), np.array(x, dtype="d"), np.array(y, dtype="d"), np.array(ex, dtype="d"), np.array(ey, dtype="d"))
    if name is None:
        plot.SetNameTitle(y_string + " vs " + x_string, y_string + " vs " + x_string)
    else:
        plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)  # Set marker color (default is blue)
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write:
        plot.Write()  # Write the plot to a ROOT file if write is True
    return plot

def extract_numbers(tree_name):
    values=tree_name.split('-')
    return float(values[0]), float(values[1])

def calculate_mean_and_error(data):
    mean = np.mean(data)
    error = np.std(data) / np.sqrt(len(data))
    return mean, error

def process_ttree(tree):
    branch_means = {}
    branch_errors = {}
    branch_names = tree.keys()  # Get all branch names
    arrays = tree.arrays(branch_names)  # Retrieve the data for all branches
    
    for branch_name in branch_names:
        data = arrays[branch_name].to_numpy()  # Convert to numpy array
        mean, error = calculate_mean_and_error(data)
        branch_means[branch_name] = mean
        branch_errors[branch_name] = error
        
    return branch_means, branch_errors

def create_multigraph_plots_hole(mean_data):
    variables = [key.replace("_mean", "") for key in mean_data if key.endswith("_mean") and key not in ["hole_mean", "DriftV_mean"]]
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kViolet, ROOT.kPink]
    marker_styles = [20, 21, 22, 23, 24, 25, 26, 27]

    for var in variables:
        multigraph = ROOT.TMultiGraph()
        canvas = ROOT.TCanvas(f"c_{var}_DriftScan", f"Canvas for {var}",1000,1000)
        legend = ROOT.TLegend(0.1, 0.7, 0.3, 0.9)

        for idx, hole in enumerate(set(mean_data["hole_mean"])):
            # Filter data for the current hole
            indices = [i for i, h in enumerate(mean_data["hole_mean"]) if h == hole]
            drift_vs = [mean_data["DriftV_mean"][i] for i in indices]
            drift_vs_error = [0 for _ in indices]  # Assuming no error for DriftV

            y_key = f"{var}_mean"
            error_key = f"{var}_error"
            y_values = [mean_data[y_key][i] for i in indices]
            y_errors = [mean_data[error_key][i] for i in indices]

            # Generate the plot
            name = f"{var}_vs_DriftV_mean_hole_{hole}"
            plot = grapherr(drift_vs, y_values, drift_vs_error, y_errors, "DriftV_mean", var, name, color=colors[idx % len(colors)], markerstyle=marker_styles[idx % len(marker_styles)], write=False)
            multigraph.Add(plot)
            legend.AddEntry(plot, f"Hole {hole}", "p")
        
        multigraph.Draw("AP")
        multigraph.GetXaxis().SetTitle("DriftV_mean")
        multigraph.GetYaxis().SetTitle(var)
        legend.Draw()
        canvas.BuildLegend(0.1, 0.7, 0.3, 0.9)  # Legend position: top-left corner
        canvas.Write()
        del canvas  # Ensure the canvas is properly deleted to avoid memory issues

def create_multigraph_plots_driftv(mean_data):
    variables = [key.replace("_mean", "") for key in mean_data if key.endswith("_mean") and key not in ["hole_mean", "DriftV_mean"]]
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kViolet, ROOT.kPink]
    marker_styles = [20, 21, 22, 23, 24, 25, 26, 27]

    for var in variables:
        multigraph = ROOT.TMultiGraph()
        canvas = ROOT.TCanvas(f"c_{var}_HoleScan", f"Canvas for {var}", 1000, 1000)
        legend = ROOT.TLegend(0.1, 0.7, 0.3, 0.9)

        for idx, driftv in enumerate(set(mean_data["DriftV_mean"])):
            # Filter data for the current DriftV
            indices = [i for i, dv in enumerate(mean_data["DriftV_mean"]) if dv == driftv]
            holes = [mean_data["hole_mean"][i] for i in indices]
            holes_error = [0 for _ in indices]  # Assuming no error for hole

            y_key = f"{var}_mean"
            error_key = f"{var}_error"
            y_values = [mean_data[y_key][i] for i in indices]
            y_errors = [mean_data[error_key][i] for i in indices]

            # Generate the plot
            name = f"{var}_vs_Hole_mean_driftv_{driftv}"
            plot = grapherr(holes, y_values, holes_error, y_errors, "Hole_mean", var, name, color=colors[idx % len(colors)], markerstyle=marker_styles[idx % len(marker_styles)], write=False)
            multigraph.Add(plot)
            legend.AddEntry(plot, f"DriftV {driftv}", "p")
        
        multigraph.Draw("AP")
        multigraph.GetXaxis().SetTitle("Hole_mean")
        multigraph.GetYaxis().SetTitle(var)
        legend.Draw()
        canvas.BuildLegend(0.1, 0.7, 0.3, 0.9)  # Legend position: top-left corner
        canvas.Write()
        del canvas  # Ensure the canvas is properly deleted to avoid memory issues

def calculate_le(distance_mm, drift_v_kv_per_cm):
    """
    Output the L/E value in um^2/V.
    """
    # Convert distance from mm to um
    distance_um = distance_mm * 1000  # 1 mm = 1000 um
    # Convert DriftV from kV/cm to V/um
    drift_v_v_per_um = drift_v_kv_per_cm * 0.1  # 1 kV/cm = 0.1 V/um
    #print(distance_mm,distance_um,drift_v_kv_per_cm,drift_v_v_per_um,distance_um / drift_v_v_per_um)
    return distance_um / drift_v_v_per_um

def extract_pressure_from_filename(filename):
    match = re.search(r'(\d+\.?\d*)mbar', filename)
    if match:
        return float(match.group(1))
    else:
        raise ValueError("No pressure value followed by 'mbar' found in the filename.")

def fit_and_calculate_diffT(L_over_E, diffusionCUT_squared_mean, L_over_E_error, diffusionCUT_squared_error, mean_data, threshold_value, err_V, target_drift_v_kv_per_cm=0.5):
    # Filter out points based on sc_integral_mean threshold
    valid_indices = [i for i in range(len(mean_data["sc_integral_mean"])) if mean_data["sc_integral_mean"][i] >= threshold_value]
    
    L_over_E = np.array([L_over_E[i] for i in valid_indices], dtype="d")
    diffusionCUT_squared_mean = np.array([diffusionCUT_squared_mean[i] for i in valid_indices], dtype="d")
    L_over_E_error = np.array([L_over_E_error[i] for i in valid_indices], dtype="d")
    diffusionCUT_squared_error = np.array([diffusionCUT_squared_error[i] for i in valid_indices], dtype="d")
    
    # Create and fit the graph
    test_graph = grapherr(L_over_E, diffusionCUT_squared_mean, L_over_E_error, diffusionCUT_squared_error, "L/E (um^{2}/V)", "diff^{2} (um^{2})", write=False)
    fit_function = ROOT.TF1("fit_function", "[0] + [1]*x", np.min(L_over_E) * 0.9, np.max(L_over_E) * 1.1)
    test_graph.Fit("fit_function", "RQ")
    test_graph.Write()
    
    slope = fit_function.GetParameter(1)  # V
    err_slope = fit_function.GetParError(1)  # V

    # Calculate DriftV and convert to V/μm
    DriftV = np.array([mean_data["DriftV_mean"][i] for i in valid_indices], dtype="d")
    DriftV_err = err_V * DriftV
    drift_v_v_per_um = 0.1 * DriftV  # 1 kV/cm = 0.1 V/μm
    drift_v_error_v_per_um = err_V * drift_v_v_per_um

    # Compute DiffT in μm²/V using NumPy
    DiffT = np.sqrt(slope / drift_v_v_per_um)
    
    # Propagate the error for DiffT using NumPy
    term1 = (0.5 * np.sqrt(1 / (drift_v_v_per_um * slope)) * err_slope) ** 2
    term2 = (0.5 * slope / (drift_v_v_per_um ** (3/2)) * drift_v_error_v_per_um) ** 2
    DiffT_error = np.sqrt(term1 + term2)

    # Calculate DiffT and its error at target_drift_v_kv_per_cm (converted to V/μm)
    target_drift_v_v_per_um = 0.1 * target_drift_v_kv_per_cm  # Convert kV/cm to V/μm
    target_DiffT = np.sqrt(slope / target_drift_v_v_per_um)
    target_DiffT_error = np.sqrt((0.5 * np.sqrt(1 / (target_drift_v_v_per_um * slope)) * err_slope) ** 2 +
                                 (0.5 * slope / (target_drift_v_v_per_um ** (3/2)) * (0.05 * target_drift_v_v_per_um)) ** 2)

    sigma_0 = np.sqrt(fit_function.GetParameter(0))  # um
    sigma_0_error = (0.5 * fit_function.GetParError(0)) / sigma_0  # um

    return target_DiffT, target_DiffT_error, sigma_0, sigma_0_error

def create_diffusion_plot_with_cut(integral_cuts, DT_scan, DT_scan_err, y_line_value,range=None):
    # Create the graph
    graph = grapherr(integral_cuts, DT_scan, [0] * len(DT_scan), DT_scan_err, "Integral cut", "Diffusion [um^{1/2}]", write=False)

    # Create a canvas
    canvas = ROOT.TCanvas("c_diffusion_plot", "Diffusion Plot with Integral Cut", 800, 600)

    # Draw the graph
    graph.Draw("AP")

    # Draw a horizontal line at y = y_line_value
    line = ROOT.TLine(min(integral_cuts), y_line_value, max(integral_cuts), y_line_value)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.SetLineStyle(2)  # Dashed line
    line.Draw("same")

    # Set graph title and axis labels
    graph.SetTitle("Diffusion vs. Integral Cut with Horizontal Line")
    graph.GetXaxis().SetTitle("Integral cut")
    graph.GetYaxis().SetTitle("Diffusion [um^{1/2}]")
    if range is not None: graph.GetYaxis().SetRangeUser(range[0],range[1])  # Set the y-axis range

    # Update and draw the canvas
    canvas.Update()
    canvas.Write()

    return canvas, graph, line

# Create the parser for command-line arguments
parser = argparse.ArgumentParser(description="Plot analyzed Diffusion data.")
parser.add_argument('inFile', type=str, help='Input analyzed root file')
parser.add_argument('-v', '--verbose', help='print more info', action='store_true')
args = parser.parse_args()

# Open the ROOT file
root_file = uproot.open(args.inFile)

# Create a list to store TTree names
tree_names = []

mean_data = {
"chi2_mean": [],
"chi2_error": [],
"GausAmp_mean": [],
"GausAmp_error": [],
"Gausmean_mean": [],
"Gausmean_error": [],
"Gaussigma_mean": [],
"Gaussigma_error": [],
"offset_mean": [],
"offset_error": [],
"nSc_mean": [],
"nSc_error": [],
"sc_integral_mean": [],
"sc_integral_error": [],
"sc_length_mean": [],
"sc_length_error": [],
"sc_width_mean": [],
"sc_width_error": [],
"sc_nhistsc_size_mean": [],
"sc_nhistsc_size_error": [],
"sc_integralsc_nhits_mean": [],
"sc_integralsc_nhits_error": [],
"sc_tgausssigma_mean": [],
"sc_tgausssigma_error": [],
"diffusionCUT_mean": [],
"diffusionCUT_error": [],
"sc_integralsc_length_mean": [],
"sc_integralsc_length_error": [],
"hole_mean": [],
"DriftV_mean": []
}
# Iterate over all items in the ROOT file
for key in root_file.keys():
    tree_name = key[:-2]
    
    # Process the TTree to get branch means and errors
    tree = root_file[tree_name]
    branch_means, branch_errors = process_ttree(tree)
    
    # Fill the dictionary with mean and error values for each branch
    for branch, mean in branch_means.items():
        mean_key = f"{branch}_mean"
        error_key = f"{branch}_error"
        if mean_key in mean_data:
            mean_data[mean_key].append(mean)
        if error_key in mean_data:
            mean_data[error_key].append(branch_errors[branch])
# Open a ROOT file for writing plots
root_file = ROOT.TFile.Open(f"Plots_{args.inFile}", "RECREATE")

if args.verbose:
    print("Creating multigraph plots for DriftScan and HoleScan")
    root_file.mkdir("DriftScan")
    root_file.cd("DriftScan")
    create_multigraph_plots_driftv(mean_data)
    root_file.mkdir("HoleScan")
    root_file.cd("HoleScan")
    create_multigraph_plots_hole(mean_data)

root_file.cd()
distances_mm = [distance(hole) for hole in mean_data["hole_mean"]]
l_over_e = [calculate_le(dist, drift_v) for dist, drift_v in zip(distances_mm, mean_data["DriftV_mean"])]
l_over_e_error = err_LoverE*np.array(l_over_e,dtype="d")  # Assume 10% error for L/E

# Calculate the square of diffusionCUT_mean in μm and propagate the error
diffusionCUT_squared_mean = [(x * pixel_size) ** 2 for x in mean_data["diffusionCUT_mean"]]
diffusionCUT_squared_error = [2 * x * pixel_size * ex * pixel_size for x, ex in zip(mean_data["diffusionCUT_mean"], mean_data["diffusionCUT_error"])]

# Create and fit a graph for L/E vs. diffusion squared
test_graph = grapherr(l_over_e, diffusionCUT_squared_mean, l_over_e_error, diffusionCUT_squared_error, "L/E (um^{2}/V)", "diff^{2} (um^{2})", write=False)
fit_function = ROOT.TF1("fit_function", "[0] +[1]*x", np.min(l_over_e) * 0.9, np.max(l_over_e) * 1.1)
test_graph.Fit("fit_function", "RQ")
test_graph.Write()
slope = fit_function.GetParameter(1)  # V
err_slope = fit_function.GetParError(1)  # V
sigma_0 = np.sqrt(fit_function.GetParameter(0))  # um
sigma_0_error = (0.5 * fit_function.GetParError(0)) / sigma_0  # um
print("sigma_0 [um]:", sigma_0, "+/-", sigma_0_error) 

DriftV=np.array(mean_data["DriftV_mean"],dtype="d")
DriftV_err=err_V*DriftV
# Convert DriftV_mean from kV/cm to V/μm
drift_v_v_per_um = 0.1*DriftV  # 1 kV/cm = 0.1 V/μm
drift_v_error_v_per_um = err_V*drift_v_v_per_um

# Compute DiffT in um^{1/2} using NumPy
DiffT = np.sqrt(slope / drift_v_v_per_um)
# Propagate the error for DiffT using NumPy
term1 = (0.5 * np.sqrt(1 / (drift_v_v_per_um * slope)) * err_slope) ** 2
term2 = (0.5 * slope / (drift_v_v_per_um ** (3/2)) * drift_v_error_v_per_um) ** 2
DiffT_error = np.sqrt(term1 + term2)
#print(DiffT,DriftV)

# Create the final graph for diffusion coefficient vs. E
test_graph = grapherr(DriftV, DiffT, DriftV_err, DiffT_error, "E(V/cm)", "Diff coef. [um^{2}/V]")

pressure = extract_pressure_from_filename(args.inFile)
if pressure==650: GarDiff=1.306
if pressure==575: GarDiff=1.332
if pressure==500: GarDiff=1.279

root_file.mkdir("LYscan")
root_file.cd("LYscan")

integral_cuts = np.linspace(min(mean_data["sc_integral_mean"]), max(mean_data["sc_integral_mean"]) * 0.8, 20)
#print(len(integral_cuts))

DT_scan,DT_scan_err,sigma_0_scan,sigma_0_scan_err = [],[],[],[]
for cut in tqdm(integral_cuts, desc="Scanning Cut"):
    temp=fit_and_calculate_diffT(l_over_e, diffusionCUT_squared_mean, l_over_e_error, diffusionCUT_squared_error, mean_data, cut, err_V)
    DT_scan.append(temp[0])
    DT_scan_err.append(temp[1])
    sigma_0_scan.append(temp[2])
    sigma_0_scan_err.append(temp[3])

root_file.cd()
create_diffusion_plot_with_cut(integral_cuts, DT_scan, DT_scan_err, GarDiff,range=[0.9*min(DT_scan), max([1.35,1.1*max(DT_scan)])])
create_diffusion_plot_with_cut(integral_cuts, sigma_0_scan, sigma_0_scan_err, 1000000)