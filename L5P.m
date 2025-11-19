function [cf, G] = L5P(x, y, varargin)
%L5P Five-parameter logistic regression (5PL).
%
%   [cf, G] = L5P(x, y)
%   [cf, G] = L5P(x, y, st)
%   [cf, G] = L5P(x, y, st, L, U)
%
%   Description:
%       L5P fits a five-parameter logistic (5PL) model to experimental
%       data, commonly used in bioassays and immunoassays (e.g., ELISA),
%       dose-response curves, and standard curves. The 5PL model extends
%       the classic 4PL by adding an asymmetry parameter E to handle
%       non-symmetrical sigmoidal curves.
%
%       The 5PL equation is:
%
%           F(x) = D + (A - D) / (1 + (x / C)^B)^E
%
%       where:
%           A : Minimum asymptote (response at zero concentration)
%           B : Hill slope (steepness and direction of the curve)
%           C : Inflection point (EC50)
%           D : Maximum asymptote (response at infinite concentration)
%           E : Asymmetry factor (E = 1 gives a symmetric 4PL curve)
%
%   Inputs:
%       x   - Column vector (N x 1) of x-values (e.g., concentrations).
%       y   - Column vector (N x 1) of responses, or matrix (N x M) of
%             replicates. If y is a matrix, L5P uses the row-wise mean as
%             the response and the row-wise standard deviation as weights.
%
%       st  - (Optional) 1x5 row vector of starting values:
%                 [A0, B0, C0, D0, E0]
%             If empty or omitted, starting points are estimated from data.
%
%       L   - (Optional) 1x5 row vector of lower bounds for [A, B, C, D, E].
%             If empty or omitted, bounds are inferred from data.
%
%       U   - (Optional) 1x5 row vector of upper bounds for [A, B, C, D, E].
%             If empty or omitted, bounds are inferred from data.
%
%   Outputs:
%       cf  - cfit object representing the fitted 5PL curve.
%       G   - Structure with goodness-of-fit measures:
%                 G.sse, G.rsquare, G.adjrsquare, G.dfe, G.rmse
%
%   Example:
%       x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L5P(x, y);
%
%       plot(x, y, 'ro'); hold on; plot(cf, 'r'); hold off;
%
%   GitHub repository:
%       https://github.com/dnafinder/logistic5
%
%   Citation:
%       Cardillo G. (2025). logistic5: Five-parameter logistic regression
%       tools in MATLAB (L5P and L5Pinv). Available at:
%       https://github.com/dnafinder/logistic5
%
%   License:
%       This function is distributed under the terms specified in the
%       LICENSE file of the logistic5 repository.
%
%   Author:
%       Giuseppe Cardillo
%       giuseppe.cardillo.75@gmail.com
%
%   Created:
%       2012-01-01 (original concept)
%
%   Updated:
%       2025-11-19 (refactored and documented version)
%
%   Version:
%       1.1.0
%
%   See also:
%       L5Pinv, L4P, L4Pinv, L3P, L3Pinv

% -----------------------------
% Input parsing and validation
% -----------------------------
p = inputParser;

% x: numeric, real, finite, non-empty column vector
addRequired(p, 'x', @(v) validateattributes(v, ...
    {'numeric'}, {'column', 'real', 'finite', 'nonnan', 'nonempty'}));

% y: numeric, real, finite, 2D, non-empty (vector or matrix)
addRequired(p, 'y', @(v) validateattributes(v, ...
    {'numeric'}, {'2d', 'real', 'finite', 'nonnan', 'nonempty'}));

% st: optional starting values [A0, B0, C0, D0, E0]
addOptional(p, 'st', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'finite', 'nonnan', 'numel', 5}));

% L: optional lower bounds
addOptional(p, 'L', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'nonnan', 'numel', 5}));

% U: optional upper bounds
addOptional(p, 'U', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'nonnan', 'numel', 5}));

parse(p, x, y, varargin{:});

x  = p.Results.x;
y  = p.Results.y;
st = p.Results.st;
L  = p.Results.L;
U  = p.Results.U;

clear p;

% Ensure compatible sizes and minimum number of data points
assert(size(x, 1) == size(y, 1), ...
    'L5P:SizeMismatch', 'x and y must have the same number of rows.');
assert(size(x, 1) >= 5, ...
    'L5P:NotEnoughData', 'At least 5 data points are required to fit a 5PL model.');

% -----------------------------
% Handle replicate measurements
% -----------------------------
% If y is a matrix, use row-wise means as responses and row-wise standard
% deviations as weights. Otherwise, y is treated as a single response vector.
if ~isvector(y)
    we = std(y, 0, 2);   % weights: standard deviation per row
    y  = mean(y, 2);     % response: mean per row
else
    y  = y(:);           % ensure column vector
    we = zeros(size(x)); % no weights
end

% -----------------------------
% Estimate starting values
% -----------------------------
% Heuristic estimates:
%   A0: min(y)
%   B0: sign of the slope between first and last data point
%   C0: x at which y is closest to the midpoint between min(y) and max(y)
%   D0: max(y)
%   E0: 1 (no asymmetry -> reduces to 4PL)

slope = (y(end) - y(1)) / (x(end) - x(1));

if isempty(st)
    yMin   = min(y);
    yMax   = max(y);
    yMid   = yMin + (yMax - yMin) / 2;  % midpoint in y-range
    [~, Idx] = min(abs(y - yMid));
    
    A0 = yMin;
    B0 = sign(slope);
    C0 = x(Idx);
    D0 = yMax;
    E0 = 1;          % initial guess: symmetric curve
    
    st = [A0, B0, C0, D0, E0];
end

% -----------------------------
% Parameter bounds
% -----------------------------
% Default logic:
%   A, C, D, E are constrained to be >= 0 (unless user overrides).
%   B is constrained based on the observed slope:
%       - if slope >= 0: B >= 0
%       - if slope < 0 : B can be negative, but capped at 0 from above.

if isempty(L)
    L = zeros(1, 5);  % [A, B, C, D, E]
    if slope < 0
        L(2) = -Inf;  % allow negative slopes
    end
end

if isempty(U)
    U = Inf(1, 5);    % [A, B, C, D, E]
    if slope < 0
        U(2) = 0;     % prevent positive slopes in a decreasing curve
    end
end

% -----------------------------
% Define model and fit options
% -----------------------------
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', st, ...
    'Lower',      L, ...
    'Upper',      U);

% If all weights are positive, use them for weighted regression
if all(we)
    set(fo, 'Weights', we);
end

% Define the 5PL model as a fittype object
ft = fittype('D + (A - D) / ((1 + (x / C)^B)^E)', ...
    'dependent',   {'y'}, ...
    'independent', {'x'}, ...
    'coefficients',{'A', 'B', 'C', 'D', 'E'});

% -----------------------------
% Fit the 5PL model
% -----------------------------
[cf, G] = fit(x, y, ft, fo);

end
