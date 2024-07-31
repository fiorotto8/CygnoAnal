import Analyzer
import argparse
from glob import glob
import os
import uproot
import ROOT
from tqdm import tqdm
import numpy as np
import pandas as pd
import uproot3
import awkward as ak
import sys
import re
import time

def hist(data, x_name, channels=100, linecolor=4, linewidth=4, write=True):
    array = np.array(data, dtype="d")
    hist = ROOT.TH1D(x_name, x_name, channels, 0.9 * np.min(array), 1.1 * np.max(array))
    hist.FillN(len(array), array, np.ones(len(array)))  # Optimized histogram filling
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetMaxDigits(3)
    if write:
        hist.Write()
    return hist

def convert_data_dict_to_tree_format(data_dict):
    tree_dict = {key: np.array(value, dtype=object) if isinstance(value, list) else value for key, value in data_dict.items()}
    return tree_dict

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory {directory} created.")
    else:
        print(f"Directory {directory} already exists.")

def append_event_data(event_data, events_arrays, j, k, stdCUT, trk, gaus_pars, offset, chi2,phi):
    # Append the data to the event_data dictionary
    event_data["nSc"].append(events_arrays["nSc"][j])
    event_data["sc_integral"].append(events_arrays["sc_integral"][j][k])
    event_data["sc_length"].append(events_arrays["sc_length"][j][k])
    event_data["sc_width"].append(events_arrays["sc_width"][j][k])
    event_data["sc_tgausssigma"].append(events_arrays["sc_tgausssigma"][j][k])
    event_data["diffusionCUT"].append(stdCUT)
    event_data["PhiMainAxis"].append(phi)
    event_data["sc_integralsc_nhits"].append(events_arrays["sc_integral"][j][k] / events_arrays["sc_nhits"][j][k])
    event_data["sc_integralsc_length"].append(events_arrays["sc_integral"][j][k] / events_arrays["sc_length"][j][k])
    event_data["sc_nhistsc_size"].append(events_arrays["sc_nhits"][j][k] / events_arrays["sc_size"][j][k])
    event_data["GausAmp"].append(gaus_pars[0] if gaus_pars else None)
    event_data["Gausmean"].append(gaus_pars[1] if gaus_pars else None)
    event_data["Gaussigma"].append(gaus_pars[2] if gaus_pars else None)
    event_data["offset"].append(offset)
    event_data["chi2"].append(chi2)

def process_otherparam_list(param_list, data_dict, data_key):
    if param_list:
        mean_value = np.mean(param_list)
        data_dict[data_key].append(mean_value)

def extract_number_from_filename(filename):
    """
    Extract the number from a given filename string.

    Parameters:
    filename (str): The input filename string.

    Returns:
    float: The extracted number as a float.
    """
    # Use regular expression to find the number in the filename
    match = re.search(r"(\d+\.?\d*)mbar", filename)
    if match:
        # Convert the extracted number to float and return
        return float(match.group(1))
    else:
        raise ValueError("No number found in the given filename.")

# Initialize an empty dictionary to store aggregated data
data_dict = {
    "temperature_mean": [],
    "pressure_mean": [],
    "humidity_mean": [],
    "hole_mean": [],
    "DriftV_mean": []
}

# Create the parser
parser = argparse.ArgumentParser(description="Analyse one diffusion run.")
parser.add_argument('inDir', type=str, help='Input folder with the reco file with DriftV and numHole')
parser.add_argument('-v', '--verbose', help='print more info', action='store_true')
parser.add_argument('-d', '--draw', help='draw and save the tracks', action='store_true')
parser.add_argument('-n', '--num_events', type=int, default=None, help='Number of events to process before stopping. Process all if not specified.')
args = parser.parse_args()

start_time = time.time()

file_list = glob(os.path.join(args.inDir, '**', '*'), recursive=True)

ROOT.gROOT.SetBatch(True)
int_cut =1000
pressure_folder=extract_number_from_filename(args.inDir)
if pressure_folder == 500: length_cut,dEdx_cut = 300,70 
elif pressure_folder == 575: length_cut,dEdx_cut = 300,120
elif pressure_folder==650: length_cut,dEdx_cut = 200,150
else: 
    print("Pressure not found")
    sys.exit(1)

output_dir = "Tracks/"
create_directory(output_dir)

output_file_name = f'output_{args.inDir.replace("/", "")}.root'

dictionaries=[]
# Loop over the file list and open each ROOT file
#for i, file in enumerate(file_list):
for i, file in enumerate(file_list):
    print(f"File: {file}")
    #! open the file TTree
    with uproot.open(file) as root_file:
        #! Events parameters
        if "Events" in root_file:
            #Empty structure to store the data
            event_data = {"nSc": [],
            "sc_integral": [],
            "sc_length": [],
            "sc_width": [],
            "sc_tgausssigma": [],
            "diffusionCUT": [],
            "PhiMainAxis": [],
            "sc_integralsc_nhits": [],
            "sc_integralsc_length": [],
            "sc_nhistsc_size": [],
            "GausAmp": [],
            "Gausmean": [],
            "Gaussigma": [],
            "offset": [],
            "chi2": []}
            events_tree = root_file["Events"]
            events_arrays = events_tree.arrays(["nSc", "nRedpix", "sc_redpixIdx", "sc_integral", "sc_length", "sc_width", "sc_tgausssigma", "redpix_ix", "redpix_iy", "redpix_iz","sc_nhits","sc_size"],library="pd")        
            #! cycle over the events
            for j, sc in tqdm(enumerate(events_arrays["nSc"]), total=len(events_arrays["nSc"]), desc="Processing Events"):
                if args.num_events is not None and j >= args.num_events: break
                if args.verbose: print(f"Event {j}")
                nSc_red, B, E = Analyzer.ScIndicesElem(events_arrays["nSc"][j], events_arrays["nRedpix"][j], events_arrays["sc_redpixIdx"][j])
                lengths = [E[l] - B[l] for l in range(len(B))]
                #! cycle over the superclusers
                for k in range(nSc_red):
                    if args.verbose: print(f"nSc_red {k}")
                    if nSc_red > 0 and lengths[k] > 0 and events_arrays["sc_integral"][j][k] > int_cut and events_arrays["sc_length"][j][k] > length_cut and events_arrays["sc_integral"][j][k]/events_arrays["sc_length"][j][k]<dEdx_cut:
                        trk = Analyzer.Track(f"Event_file{i}_image{j}_SC{k}", events_arrays["redpix_ix"][j], events_arrays["redpix_iy"][j], events_arrays["redpix_iz"][j], B[k], E[k])
                        stdCUT,gaus_pars,offset,chi2 = trk.GetSigmaAroundBar()
                        if args.draw: trk.save_histogram(output_dir)
                        # Use the function to append data
                        if chi2 < 2: 
                            append_event_data(event_data, events_arrays, j, k, stdCUT, trk, gaus_pars, offset, chi2,trk.fPhiMainAxis)
                        #else: trk.save_histogram(output_dir)
                        del trk

            #! File parameters
            if "OtherParam" in root_file:
                otherparam_tree = root_file["OtherParam"]
                otherparam_arrays = otherparam_tree.arrays(library="pd")

                try:
                    temperature_list = otherparam_arrays["KEG_temp"].tolist()
                    pressure_list = otherparam_arrays["KEG_pressure"].tolist()
                    humidity_list = otherparam_arrays["KEG_humidity"].tolist()
                except KeyError:
                    temperature_list = otherparam_arrays["KEG_t"].tolist()
                    pressure_list = otherparam_arrays["KEG_p"].tolist()
                    humidity_list = otherparam_arrays["KEG_h"].tolist()
                
                hole_list = otherparam_arrays["HOLE_number"].tolist()
                DriftV_list = otherparam_arrays["DRIFT_V"].tolist()

                process_otherparam_list(temperature_list, data_dict, "temperature_mean")
                process_otherparam_list(pressure_list, data_dict, "pressure_mean")
                process_otherparam_list(humidity_list, data_dict, "humidity_mean")
                process_otherparam_list(hole_list, data_dict, "hole_mean")
                process_otherparam_list(DriftV_list, data_dict, "DriftV_mean")
            else:
                print(f"'OtherParam' TTree not found in {file}")

        dictionaries.append(event_data)

# Create a new ROOT file
with uproot.recreate(output_file_name) as output_file:
    # Loop over dictionaries and write each one to a TTree
    for m, event_data in enumerate(dictionaries):
        tree_name = f'{data_dict["hole_mean"][m]}-{data_dict["DriftV_mean"][m]:.3f}'
        output_file[tree_name] = {key: np.array(value) for key, value in event_data.items()}
        
stop_time = time.time()
print("Code Execution Time: ", stop_time-start_time)