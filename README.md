[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/logistic5)

📘 Overview
logistic5 is a MATLAB toolbox implementing the Five-Parameter Logistic (5PL) regression model, an extension of the widely used 4PL model in immunoassays, dose-response analysis, and bioassay standard curve fitting. The additional parameter E introduces asymmetry into the sigmoid curve, enabling more flexible modeling of real-world biological data where the curve is not perfectly symmetrical around its inflection point.

The 5PL model is defined as:
    F(x) = D + (A - D) / (1 + (x / C)^B)^E

where:
- A : Minimum asymptote (response at zero analyte concentration)
- B : Hill slope (defines curve steepness and direction)
- C : Inflection point or EC50
- D : Maximum asymptote (response at infinite concentration)
- E : Asymmetry factor (E = 1 reduces to the symmetric 4PL model)

This repository provides two key functions:
- L5P    : Performs 5PL curve fitting using nonlinear least squares
- L5Pinv : Computes the inverse of the 5PL model, estimating concentrations from response values

✨ Features
- Robust nonlinear regression of 5PL curves
- Automatic detection of initial parameter estimates when not supplied
- Bounds for parameters automatically chosen based on data if not provided
- Support for replicate measurements: weights derived from row-wise standard deviations
- Returns a cfit MATLAB object plus detailed goodness-of-fit statistics
- Includes full inverse evaluation for interpolation/extrapolation
- Suitable for assays requiring asymmetric sigmoid fitting, including biochemical, pharmacological, and immunological workflows

📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/logistic5

2. Add the folder to your MATLAB path:
      addpath('path_to_logistic5')

3. Verify successful installation:
      which L5P
      which L5Pinv

⚙️ Requirements
- MATLAB (any recent version)
- Curve Fitting Toolbox (required for fit, fittype, cfit models)

📈 Usage
Fitting a 5PL curve:

    x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
    y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];

    [cf, G] = L5P(x, y);

Plotting the results:

    plot(x, y, 'ro');
    hold on;
    plot(cf, 'r');
    hold off;

Interpolating unknown samples:

    yq = 1.8;
    x_est = L5Pinv(cf, yq);

Using explicit parameters instead of a cfit object:

    params = [A, B, C, D, E];
    x_est = L5Pinv(params, yq);

🔢 Inputs
L5P(x, y)
L5P(x, y, st)
L5P(x, y, st, L, U)

- x : Column vector of concentrations (N×1)
- y : Column vector OR matrix (N×M) of responses
- st: Optional starting values [A0 B0 C0 D0 E0]
- L : Optional lower bounds for [A B C D E]
- U : Optional upper bounds for [A B C D E]

L5Pinv(cf, y)
- cf : cfit object from L5P OR numeric vector [A B C D E]
- y  : Query response values (any shape)

📤 Outputs
L5P returns:
- cf : cfit object for the fitted 5PL curve
- G  : Structure containing sse, rsquare, adjrsquare, rmse, dfe

L5Pinv returns:
- x : Interpolated concentration values corresponding to y

🧠 Interpretation
- E introduces curve asymmetry; E > 1 bends the curve in one direction, E < 1 in the opposite direction
- B controls the steepness of the curve, positive or negative
- C represents the EC50 and defines where the curve transitions most rapidly
- Strong fits show:
  • High R² and adjusted R²  
  • Low RMSE and SSE  
  • Visually coherent sigmoidal shape following experimental points  

📌 Notes
- When providing replicates (matrix y), L5P automatically computes row means and uses standard deviations as weights.
- Good starting values and reasonable parameter bounds can improve convergence.
- Values of y equal to A or D cannot be inverted directly; L5Pinv will warn about extrapolation.

🧾 Citation
If you use logistic5 in research, analysis, or publications, please cite:

Cardillo G. (2025). logistic5: Five-parameter logistic regression tools in MATLAB (L5P and L5Pinv).  
Available at: https://github.com/dnafinder/logistic5

👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

📄 License
logistic5 is distributed under the terms specified in the LICENSE file:  
https://github.com/dnafinder/logistic5
