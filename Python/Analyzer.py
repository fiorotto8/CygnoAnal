import numpy as np
import ROOT

def ScIndicesElem(nSc, npix, sc_redpixID):
    sc_redpixID = np.array(sc_redpixID)
    sc_redpix_start = np.concatenate(([0], sc_redpixID[sc_redpixID > 0], [npix]))
    nSc_red = len(sc_redpix_start) - 1

    B = sc_redpix_start[:-1].tolist()
    E = sc_redpix_start[1:].tolist()

    return nSc_red, B, E

def Theta1Calc(x, y, xbar, ybar, theta):
    """
    Calculate the angle theta1 such that tan(theta1 + theta) = (y - ybar) / (x - xbar).

    Parameters:
    x (float): x-coordinate of the first point.
    y (float): y-coordinate of the first point.
    xbar (float): x-coordinate of the second point.
    ybar (float): y-coordinate of the second point.
    theta (float): base angle.

    Returns:
    float: angle theta1.
    """
    # Calculate the angle of the line through (x, y) and (xbar, ybar) relative to the x-axis
    angle = np.arctan2(y - ybar, x - xbar)
    
    # Subtract the base angle theta to get theta1
    theta1 = angle - theta
    
    return theta1

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
    return plot

class Track:
    def __init__(self, nometh2, X, Y, Z, B, E):
        self.fradius = 0.
        self.fRange = None
        self.fheight = 0.
        self.fxcentr = 0
        self.fycentr = 0
        self.fTrackTail = None
        self.fScaledTrack = None
        self.fNPIP = 0
        self.fwScal = 0.
        self.fXbar = 0.
        self.fYbar = 0.
        self.fBarPlot = None
        self.fXIPPrev = 0.
        self.fYIPPrev = 0.
        self.fXIP = 0.
        self.fYIP = 0.
        self.fIPPlot = None
        self.fPhiMainAxis = 0.
        self.fLineDirection = None
        self.fintegral = 0.
        self.fname = nometh2
        self.SigmaDistrTH= None
        self.chi2 = None

        X, Y, Z = np.array(X), np.array(Y), np.array(Z)
        
        B = int(B)
        E = int(E)

        # Efficiently calculate min and max
        self.fminx, self.fmaxx = np.min(X[B:E]), np.max(X[B:E])
        self.fminy, self.fmaxy = np.min(Y[B:E]), np.max(Y[B:E])

        self.fmaxx += 100
        self.fminx -= 100
        self.fmaxy += 200
        self.fminy -= 50

        self.fnpixelx = int(self.fmaxx - self.fminx)
        self.fnpixely = int(self.fmaxy - self.fminy)
        
        self.fTrack = ROOT.TH2F(f"A{nometh2}", f"A{nometh2}", self.fnpixelx, self.fminx, self.fmaxx, self.fnpixely, self.fminy, self.fmaxy)

        # Use boolean indexing to filter out Z <= 0 before the loop
        valid_indices = Z[B:E] > 0
        X_valid = X[B:E][valid_indices]
        Y_valid = Y[B:E][valid_indices]
        Z_valid = Z[B:E][valid_indices]

        for x, y, z in zip(X_valid, Y_valid, Z_valid):
            self.fTrack.Fill(x, y, z)
            self.fintegral += z
        #self.ApplyThr(1)
        #self.RemoveNoise(1)

        self.fXbar, self.fYbar = self.Barycenter()
        self.fPhiMainAxis = self.AngleLineMaxRMS()

    def __del__(self):
        #print(f"Track object {self.fname} is being deleted.")
        if self.fTrack:
            del self.fTrack
        if self.SigmaDistrTH:
            del self.SigmaDistrTH

    def Barycenter(self):
        Xb = 0
        Yb = 0
        ChargeTot = 0

        x_centers = np.array([self.fTrack.GetXaxis().GetBinCenter(i) for i in range(1, self.fnpixelx + 1)])
        y_centers = np.array([self.fTrack.GetYaxis().GetBinCenter(j) for j in range(1, self.fnpixely + 1)])

        for i in range(1, self.fnpixelx + 1):
            for j in range(1, self.fnpixely + 1):
                Z = self.fTrack.GetBinContent(i, j)
                if Z > 0:
                    Xb += Z * x_centers[i-1]
                    Yb += Z * y_centers[j-1]
                    ChargeTot += Z

        if ChargeTot > 0:
            Xb /= ChargeTot
            Yb /= ChargeTot

        return Xb, Yb

    def AngleLineMaxRMS(self):
        """
        This function calculates the angle of the line that maximizes the Root Mean Square (RMS) value 
        along the primary axis of a 2D histogram. The histogram is represented by `self.fTrack`, and 
        the function identifies the principal axis of the data distribution by evaluating the RMS 
        along different lines passing through the centroid of the histogram data. The function 
        returns the angle (in radians) of the main axis that has the maximum RMS value.
        """
        # Initialize summation variables
        Sum1 = 0
        Sum2 = 0

        # Generate arrays of bin centers for the x and y axes
        x_centers = np.array([self.fTrack.GetXaxis().GetBinCenter(i) for i in range(1, self.fnpixelx + 1)])
        y_centers = np.array([self.fTrack.GetYaxis().GetBinCenter(j) for j in range(1, self.fnpixely + 1)])

        # Iterate over each pixel
        for i in range(1, self.fnpixelx + 1):
            for j in range(1, self.fnpixely + 1):
                # Get the content (Z value) of the current bin
                Z = self.fTrack.GetBinContent(i, j)
                if Z > 0:
                    # Calculate the coordinates of the current bin center
                    x = x_centers[i-1]
                    y = y_centers[j-1]
                    # Update the summations
                    Sum1 += Z * (x - self.fXbar) * (y - self.fYbar)
                    Sum2 += Z * ((y - self.fYbar) ** 2 - (x - self.fXbar) ** 2)

        # Calculate the angle Phi using the arctan2 function
        Phi = -0.5 * np.arctan2(2 * Sum1, Sum2)

        # Calculate the RMS on the line at angle Phi and its perpendicular
        RmsAng = self.RMSOnLine(Phi)
        RmsAngPerp = self.RMSOnLine(Phi + np.pi / 2)

        # Determine which RMS is larger and set the main axis accordingly
        if RmsAng > RmsAngPerp:
            self.fRMSOnMainAxis = RmsAng
            return Phi
        else:
            self.fRMSOnMainAxis = RmsAngPerp
            # Ensure the angle is within the range [-pi/2, pi/2]
            return Phi + np.pi / 2 if Phi + np.pi / 2 <= np.pi / 2 else Phi + np.pi / 2 - np.pi

    def RMSOnLine(self, phi):
        rms = 0
        charge_tot = 0

        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        for i in range(1, self.fnpixelx + 1):
            for j in range(1, self.fnpixely + 1):
                z = self.fTrack.GetBinContent(i, j)
                if z != 0:
                    x_centered = self.fTrack.GetXaxis().GetBinCenter(i) - self.fXbar
                    y_centered = self.fTrack.GetYaxis().GetBinCenter(j) - self.fYbar
                    projection = x_centered * cos_phi + y_centered * sin_phi
                    rms += z * projection * projection
                    charge_tot += z

        if charge_tot > 0:
            rms = np.sqrt(rms / charge_tot)

        return rms

    def plot_histogram(self):
        canvas = ROOT.TCanvas("canvas", "Track Histogram", 1000, 1000)
        self.fTrack.Draw("COLZ")

        # Get the axis ranges of the histogram
        x_min = self.fTrack.GetXaxis().GetXmin()
        x_max = self.fTrack.GetXaxis().GetXmax()
        y_min = self.fTrack.GetYaxis().GetXmin()
        y_max = self.fTrack.GetYaxis().GetXmax()
        
        # Calculate the length of the line to cover the entire histogram
        line_length = max(x_max - x_min, y_max - y_min)

        # Calculate the endpoints of the main axis line within the histogram bounds
        x1 = self.fXbar - 0.5 * line_length * np.cos(self.fPhiMainAxis)
        y1 = self.fYbar - 0.5 * line_length * np.sin(self.fPhiMainAxis)
        x2 = self.fXbar + 0.5 * line_length * np.cos(self.fPhiMainAxis)
        y2 = self.fYbar + 0.5 * line_length * np.sin(self.fPhiMainAxis)
        
        # Draw the main axis line
        main_line = ROOT.TLine(x1, y1, x2, y2)
        main_line.SetLineColor(ROOT.kRed)
        main_line.SetLineWidth(2)
        main_line.Draw("SAME")

        if self.fRange is not None:
            # Calculate the coordinates for the range lines
            cosPhi = np.cos(self.fPhiMainAxis)
            sinPhi = np.sin(self.fPhiMainAxis)
            cosPerpPhi = np.cos(self.fPhiMainAxis + np.pi / 2)
            sinPerpPhi = np.sin(self.fPhiMainAxis + np.pi / 2)

            x1 = self.fXbar - self.fRange * cosPhi
            y1 = self.fYbar - self.fRange * sinPhi
            x2 = self.fXbar + self.fRange * cosPhi
            y2 = self.fYbar + self.fRange * sinPhi

            # Draw the first perpendicular delimiting line
            line1 = ROOT.TLine(x1 + (x_min - self.fXbar) * cosPerpPhi, y1 + (y_min - self.fYbar) * sinPerpPhi,
                               x1 + (x_max - self.fXbar) * cosPerpPhi, y1 + (y_max - self.fYbar) * sinPerpPhi)
            line1.SetLineColor(ROOT.kRed)
            line1.SetLineStyle(2)  # Dashed line
            line1.Draw("SAME")

            # Draw the second perpendicular delimiting line
            line2 = ROOT.TLine(x2 + (x_min - self.fXbar) * cosPerpPhi, y2 + (y_min - self.fYbar) * sinPerpPhi,
                               x2 + (x_max - self.fXbar) * cosPerpPhi, y2 + (y_max - self.fYbar) * sinPerpPhi)
            line2.SetLineColor(ROOT.kRed)
            line2.SetLineStyle(2)  # Dashed line
            line2.Draw("SAME")

        if self.SigmaDistrTH is not None:
            # Create a new pad for the histogram
            pad = ROOT.TPad("pad", "pad", 0.11, 0.5, 0.5, 0.89)
            pad.Draw()
            pad.cd()
            ROOT.gStyle.SetOptFit(1111)
            self.SigmaDistrTH.Draw()
            # Add TPaveText to display chi2 value
            canvas.cd()  # Go back to the main canvas
            pave_text = ROOT.TPaveText(0.4, 0.93, 0.6, 0.98, "NDC")  # Normalized device coordinates
            pave_text.AddText(f"Chi2: {self.chi2:.2f}")
            pave_text.SetFillColor(0)  # Transparent fill
            pave_text.SetTextAlign(22)  # Center alignment
            pave_text.Draw()

        canvas.Draw()      

        input("Press Enter to continue...")
        del canvas, main_line,line1,line2

    def save_histogram(self, directory):
        """
        Saves a pictorial representation showing the directionality of the track, with lines indicating a specified range.
        
        Args:
            directory (str): The directory to save the picture in.
        """
        # Create a canvas to draw the histogram
        canvas = ROOT.TCanvas("canvas", "canvas", 1500, 1500)
        self.fTrack.Draw("COLZ")

        # Get the axis ranges of the histogram
        x_min = self.fTrack.GetXaxis().GetXmin()
        x_max = self.fTrack.GetXaxis().GetXmax()
        y_min = self.fTrack.GetYaxis().GetXmin()
        y_max = self.fTrack.GetYaxis().GetXmax()
        
        # Calculate the length of the line to cover the entire histogram
        line_length = max(x_max - x_min, y_max - y_min)

        # Calculate the endpoints of the main axis line within the histogram bounds
        x1 = self.fXbar - 0.5 * line_length * np.cos(self.fPhiMainAxis)
        y1 = self.fYbar - 0.5 * line_length * np.sin(self.fPhiMainAxis)
        x2 = self.fXbar + 0.5 * line_length * np.cos(self.fPhiMainAxis)
        y2 = self.fYbar + 0.5 * line_length * np.sin(self.fPhiMainAxis)
        
        # Draw the main axis line
        main_line = ROOT.TLine(x1, y1, x2, y2)
        main_line.SetLineColor(ROOT.kRed)
        main_line.SetLineWidth(2)
        main_line.Draw("SAME")

        if self.fRange is not None:
            # Calculate the coordinates for the range lines
            cosPhi = np.cos(self.fPhiMainAxis)
            sinPhi = np.sin(self.fPhiMainAxis)
            cosPerpPhi = np.cos(self.fPhiMainAxis + np.pi / 2)
            sinPerpPhi = np.sin(self.fPhiMainAxis + np.pi / 2)

            x1 = self.fXbar - self.fRange * cosPhi
            y1 = self.fYbar - self.fRange * sinPhi
            x2 = self.fXbar + self.fRange * cosPhi
            y2 = self.fYbar + self.fRange * sinPhi

            # Draw the first perpendicular delimiting line
            line1 = ROOT.TLine(x1 + (x_min - self.fXbar) * cosPerpPhi, y1 + (y_min - self.fYbar) * sinPerpPhi,
                               x1 + (x_max - self.fXbar) * cosPerpPhi, y1 + (y_max - self.fYbar) * sinPerpPhi)
            line1.SetLineColor(ROOT.kRed)
            line1.SetLineStyle(2)  # Dashed line
            line1.Draw("SAME")

            # Draw the second perpendicular delimiting line
            line2 = ROOT.TLine(x2 + (x_min - self.fXbar) * cosPerpPhi, y2 + (y_min - self.fYbar) * sinPerpPhi,
                               x2 + (x_max - self.fXbar) * cosPerpPhi, y2 + (y_max - self.fYbar) * sinPerpPhi)
            line2.SetLineColor(ROOT.kRed)
            line2.SetLineStyle(2)  # Dashed line
            line2.Draw("SAME")

        if self.SigmaDistrTH is not None:
            # Create a new pad for the histogram
            pad = ROOT.TPad("pad", "pad", 0.11, 0.5, 0.5, 0.89)
            pad.Draw()
            pad.cd()
            ROOT.gStyle.SetOptFit(1111)
            self.SigmaDistrTH.Draw()
            # Add TPaveText to display chi2 value
            canvas.cd()  # Go back to the main canvas
            pave_text = ROOT.TPaveText(0.4, 0.93, 0.6, 0.98, "NDC")  # Normalized device coordinates
            pave_text.AddText(f"Chi2: {self.chi2:.2f}")
            pave_text.AddText(f"Theta: {self.fPhiMainAxis:.2f}")
            pave_text.SetFillColor(0)  # Transparent fill
            pave_text.SetTextAlign(22)  # Center alignment
            pave_text.Draw()

        # Save the canvas as a PNG file
        canvas.SaveAs(f"{directory}/{self.fname}.png")

        # Clean up
        del canvas
        if 'line1' in locals():
            del line1
        if 'line2' in locals():
            del line2
        del main_line


    def GetSigmaAroundBar(self, Srange=100, Fit=True):
        """
        Calculate the sigma of the track along the main axis with a cut in the vicinity of the barycenter.

        Parameters:
        - Srange (float): The range within which to consider the distances on the main axis.
        - Fit (bool): Whether to fit the distribution of distances with a Gaussian plus a constant.

        Returns:
        - weighted_std_dev (float): The calculated standard deviation (spread) of the distances perpendicular to the main axis.
        - gaus_pars (list): Parameters from the Gaussian fit (amplitude, mean, sigma).
        - offset (float): The constant offset from the fit.
        - chi2 (float): The reduced chi-squared value from the fit.
        """

        # Set the range within which the sigma is calculated
        self.fRange = Srange

        # Calculate cosine and sine of the main axis angle
        cosPhi = np.cos(self.fPhiMainAxis)
        sinPhi = np.sin(self.fPhiMainAxis)

        # Initialize lists to store distances and corresponding charges
        distances, charges = [], []

        # Precompute bin centers for x and y axes
        x_centers = np.array([self.fTrack.GetXaxis().GetBinCenter(i) for i in range(1, self.fnpixelx + 1)])
        y_centers = np.array([self.fTrack.GetYaxis().GetBinCenter(j) for j in range(1, self.fnpixely + 1)])

        # Loop over each bin in the 2D histogram
        for i in range(1, self.fnpixelx + 1):
            for j in range(1, self.fnpixely + 1):
                binContent = self.fTrack.GetBinContent(i, j)

                if binContent > 0:
                    # Coordinates of the current bin center
                    x = x_centers[i - 1]
                    y = y_centers[j - 1]

                    # Projection of the distance from the bin center to the barycenter onto the main and perpendicular axes
                    distance_onaxis = (x - self.fXbar) * cosPhi + (y - self.fYbar) * sinPhi
                    distance_perpAxis = -(x - self.fXbar) * sinPhi + (y - self.fYbar) * cosPhi

                    # Accumulate the distance perpendicular to the main axis if it is within the specified range
                    if abs(distance_onaxis) <= self.fRange:
                        distances.append(distance_perpAxis)
                        charges.append(binContent)

        # Ensure there is at least one distance value
        if len(distances) == 0:
            return 0, [], 0, 0

        # Convert lists to NumPy arrays for efficient computation
        distances = np.array(distances)
        charges = np.array(charges)

        # Calculate weighted mean and standard deviation
        weighted_mean = np.average(distances, weights=charges)
        weighted_std_dev = np.sqrt(np.average((distances - weighted_mean) ** 2, weights=charges))

        # Histogram binning setup
        n_bins = max(10, len(distances) // 10)  # Ensure a reasonable number of bins
        hist_range = (distances.min(), distances.max())
        hist, bin_edges = np.histogram(distances, bins=n_bins, weights=charges, range=hist_range)

        # Calculate bin centers
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Error on the binned charges (using Poisson statistics as an example, adjust if needed)
        err_binned_charges = np.sqrt(hist)

        # Handle bins with zero total charge to avoid division by zero
        valid_bins = hist > 0
        bin_centers = bin_centers[valid_bins]
        hist = hist[valid_bins]
        err_binned_charges = err_binned_charges[valid_bins]

        # Create a graph for fitting
        self.SigmaDistrTH = grapherr(bin_centers, hist, np.zeros_like(bin_centers), err_binned_charges, "Distance(px)", "Charge(ADC)")

        # Gaussian fit
        gaus_pars = [0, 0, 0]
        chi2 = 0

        if Fit:
            fit_function = ROOT.TF1("fit_function", "[0] * TMath::Gaus(x, [1], [2])", hist_range[0], hist_range[1])
            fit_function.SetParameters(hist.max(), weighted_mean, weighted_std_dev)

            self.SigmaDistrTH.Fit(fit_function, "RQ")

            gaus_pars = [fit_function.GetParameter(i) for i in range(3)]
            if fit_function.GetNDF() != 0:
                chi2 = fit_function.GetChisquare() / fit_function.GetNDF()

        self.chi2 = chi2

        return weighted_std_dev, gaus_pars, 0, chi2





    def GetSigmaAroundBar_old(self, Srange=100, Fit=True):
        """
        Calculate the sigma of the track along the main axis with a cut in the vicinity of the barycenter.

        Parameters:
        - Srange (float): The range within which to consider the distances on the main axis.
        - Fit (bool): Whether to fit the distribution of distances with a Gaussian plus a constant.

        Returns:
        - weighted_std_dev (float): The calculated standard deviation (spread) of the distances perpendicular to the main axis.
        - gaus_pars (list): Parameters from the Gaussian fit (amplitude, mean, sigma).
        - offset (float): The constant offset from the fit.
        - chi2 (float): The reduced chi-squared value from the fit.
        """
        
        # Set the range within which the sigma is calculated
        self.fRange = Srange

        # Calculate the cosine and sine of the main axis angle (fPhiMainAxis)
        # These are used for projecting onto the main and perpendicular axes
        cosPhi = np.cos(self.fPhiMainAxis)
        sinPhi = np.sin(self.fPhiMainAxis)

        # Initialize lists to store distances and corresponding charges
        distances, charges, err_percent = [], [], 0.1

        # Precompute bin centers for the histogram along x and y axes
        x_centers = np.array([self.fTrack.GetXaxis().GetBinCenter(i) for i in range(1, self.fnpixelx + 1)])
        y_centers = np.array([self.fTrack.GetYaxis().GetBinCenter(j) for j in range(1, self.fnpixely + 1)])

        # Loop over each bin in the 2D histogram
        for i in range(1, self.fnpixelx + 1):
            for j in range(1, self.fnpixely + 1):
                # Get the content (weight or charge) of the current bin
                binContent = self.fTrack.GetBinContent(i, j)

                # If the bin has content, proceed with the calculation
                if binContent > 0:
                    # Get the center coordinates of the current bin
                    x = x_centers[i - 1]
                    y = y_centers[j - 1]

                    # Calculate the projection of the distance from the bin center to the barycenter onto the main axis
                    distance_onaxis = (x - self.fXbar) * cosPhi + (y - self.fYbar) * sinPhi
                    
                    # Calculate the projection of the distance from the bin center to the barycenter onto the perpendicular axis
                    distance_perpAxis = -(x - self.fXbar) * sinPhi + (y - self.fYbar) * cosPhi
                    
                    # Accumulate the distance perpendicular to the main axis if it is within the specified range
                    if abs(distance_onaxis) <= self.fRange:
                        distances.append(distance_perpAxis)
                        charges.append(binContent)

        # Ensure there is at least one distance value
        if len(distances) == 0:
            distances.append(0)
            distances.append(0)
            charges.append(1)  # Add a default charge to avoid division by zero
            charges.append(1)  # Add a default charge to avoid division by zero

        # Calculate the weighted mean of the distances
        weighted_mean = np.average(distances, weights=charges)
        
        # Calculate the weighted variance of the distances
        weighted_variance = np.average((distances - weighted_mean) ** 2, weights=charges)
        
        # Calculate the weighted standard deviation (spread) from the weighted variance
        weighted_std_dev = np.sqrt(weighted_variance)

        # Define the range and number of bins for the histogram of distances
        min_distance, max_distance, n_bins = min(distances), max(distances), int(len(distances)/50)
        
        # Create the bin edges for the histogram
        bins = np.linspace(min_distance, max_distance, n_bins + 1)

        # Determine which bin each distance falls into
        bin_indices = np.digitize(distances, bins)

        # Initialize arrays to store the total charge and the count of charges per bin
        binned_charges, charge_counts, binned_charges_sq = np.zeros(n_bins), np.zeros(n_bins), np.zeros(n_bins)

        # Accumulate charges and charge counts for each bin
        for charge, bin_idx in zip(charges, bin_indices):
            if bin_idx > 0 and bin_idx <= n_bins:
                binned_charges[bin_idx - 1] += charge
                charge_counts[bin_idx - 1] += 1
                binned_charges_sq[bin_idx - 1] += charge ** 2
        
        # Calculate total charges and relative error for each bin
        with np.errstate(divide='ignore', invalid='ignore'):
            # Total charge per bin
            total_charges = binned_charges

            # Relative error of total charge per bin
            # Here, err_binned_charges is calculated as the standard deviation of the charge in each bin
            std_devs = np.sqrt(np.true_divide(binned_charges_sq, charge_counts) - np.square(total_charges / charge_counts))
            std_devs[charge_counts == 0] = 0

            err_binned_charges = std_devs / total_charges
            err_binned_charges[charge_counts == 0] = 0

        # Calculate the centers of each bin
        bin_centers = (bins[:-1] + bins[1:]) / 2
        # Calculate the error on bin centers, with safe indexing
        if len(bins) > 1:
            err_bin_centers = ((bins[1] - bins[0]) / 2) * np.ones(n_bins)
        else:
            err_bin_centers = np.zeros(n_bins)
        # Create a mask for bins where the total charge is greater than zero
        positive_charge_mask = total_charges > 0

        # Apply the mask to all relevant arrays to filter out bins with charge <= 0
        bin_centers = bin_centers[positive_charge_mask]
        total_charges = total_charges[positive_charge_mask]
        err_bin_centers = err_bin_centers[positive_charge_mask]
        err_binned_charges = err_binned_charges[positive_charge_mask]
        #forget this and use percent of the nominal value
        err_binned_charges = err_percent * total_charges

        # Create or reset the TH1 histogram of distances
        if not hasattr(self, 'SigmaDistr'):
            # If SigmaDistrTH doesn't exist, create a new graph with error bars
            self.SigmaDistrTH = grapherr(
                bin_centers, total_charges, err_bin_centers, err_binned_charges, "Distance(px)", "Charge(ADC)"
            )
        else:
            # Reset the histogram if it already exists
            self.SigmaDistrTH.Reset()

        # Define the Gaussian plus constant function for fitting
        fit_function = ROOT.TF1("fit_function", "[0] * TMath::Gaus(x, [1], [2])", min_distance, max_distance)

        # Set initial parameters for the fit: amplitude, mean, sigma, and constant
        fit_function.SetParameters(total_charges[np.argmin(np.abs(bin_centers))], np.mean(distances), np.std(distances))
        
        # Set names for the parameters
        fit_function.SetParNames("Amplitude", "Mean", "Sigma")
        """
        # Set parameter limits for the fit
        fit_function.SetParLimits(3, -2, 10)   # Constant offset
        fit_function.SetParLimits(2, 0, 30)    # Sigma
        fit_function.SetParLimits(1, -10, 10)  # Mean
        fit_function.SetParLimits(0, 0, 20)    # Amplitude
        """
        # Initialize variables for the fit parameters and chi-squared value
        gaus_pars = []
        offset = 0
        chi2 = 0

        # Perform the fit if the Fit flag is set to True
        if Fit:
            self.SigmaDistrTH.Fit("fit_function", "RQ")  # Perform the fit with quiet mode
            # Extract fit parameters: amplitude, mean, sigma
            gaus_pars = [fit_function.GetParameter(i) for i in range(3)]
            # Extract constant offset from the fit
            offset = 0
            # Calculate reduced chi-squared value if degrees of freedom are not zero
            if fit_function.GetNDF() != 0:
                chi2 = fit_function.GetChisquare() / fit_function.GetNDF()

        # Store the chi-squared value
        self.chi2 = chi2

        # Return the calculated spread, fit parameters, offset, and chi-squared value
        return weighted_std_dev, gaus_pars, offset, chi2

    def RemoveNoise(self, EnSum):
        """
        Removes isolated noise hits from the track histogram.
        A hit is considered isolated and is removed if it has no neighboring hits within a one-bin radius.

        Parameters:
        EnSum (float): The threshold for neighboring hit energy sum. Hits with neighboring energy sum less than or equal to 'EnSum' are removed.
        """
        ToRemove = []

        # Pre-fetch bin contents for efficiency
        bin_contents = np.array([[self.fTrack.GetBinContent(i, j) for j in range(self.fnpixely)] for i in range(self.fnpixelx)])

        for i in range(1, self.fnpixelx - 1):
            for j in range(1, self.fnpixely - 1):
                if bin_contents[i, j] > 0:
                    ETemp = (
                        (bin_contents[i + 1, j] if bin_contents[i + 1, j] >= 0 else 0) +
                        (bin_contents[i - 1, j] if bin_contents[i - 1, j] >= 0 else 0) +
                        (bin_contents[i, j + 1] if bin_contents[i, j + 1] >= 0 else 0) +
                        (bin_contents[i, j - 1] if bin_contents[i, j - 1] >= 0 else 0)
                    )

                    if ETemp <= EnSum:
                        ToRemove.append((i, j))

        for i, j in ToRemove:
            self.fTrack.SetBinContent(i, j, 0)

    def ApplyThr(self, EnThr):
        """
        Applies a threshold to the track histogram, setting bin contents below a specified value to zero.
        This method is useful for further noise reduction and data cleaning.

        Parameters:
        EnThr (float): The energy threshold. Bin contents below this value are set to zero.
        """
        XBinMin = self.fTrack.GetXaxis().GetFirst()
        XBinMax = XBinMin + self.fTrack.GetXaxis().GetNbins()

        YBinMin = self.fTrack.GetYaxis().GetFirst()
        YBinMax = YBinMin + self.fTrack.GetYaxis().GetNbins()

        # Pre-fetch bin contents for efficiency
        bin_contents = np.array([[self.fTrack.GetBinContent(i, j) for j in range(YBinMin, YBinMax)] for i in range(XBinMin, XBinMax)])

        for i in range(XBinMin, XBinMax):
            for j in range(YBinMin, YBinMax):
                z = bin_contents[i - XBinMin, j - YBinMin]
                if z > 0 and z <= EnThr:
                    self.fTrack.SetBinContent(i, j, 0)
        """
        Applies a threshold to the track histogram, setting bin contents below a specified value to zero.
        This method is useful for further noise reduction and data cleaning.

        Parameters:
        EnThr (float): The energy threshold. Bin contents below this value are set to zero.
        """
        XBinMin = self.fTrack.GetXaxis().GetFirst()
        XBinMax = XBinMin + self.fTrack.GetXaxis().GetNbins()

        YBinMin = self.fTrack.GetYaxis().GetFirst()
        YBinMax = YBinMin + self.fTrack.GetYaxis().GetNbins()

        bin_contents = np.array([[self.fTrack.GetBinContent(i, j) for j in range(YBinMin, YBinMax)] for i in range(XBinMin, XBinMax)])
        bin_contents[(bin_contents > 0) & (bin_contents <= EnThr)] = 0

        for i in range(XBinMin, XBinMax):
            for j in range(YBinMin, YBinMax):
                self.fTrack.SetBinContent(i, j, bin_contents[i - XBinMin, j - YBinMin])