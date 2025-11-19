function x = L5Pinv(cf, y)
%L5PINV Inverse of the five-parameter logistic (5PL) equation.
%
%   x = L5Pinv(cf, y)
%
%   Description:
%       L5Pinv computes the inverse of the five-parameter logistic (5PL)
%       model used in L5P. Given a fitted 5PL curve (or an explicit 5PL
%       parameter vector) and a set of response values y, this function
%       returns the corresponding x values that satisfy:
%
%           y = D + (A - D) / (1 + (x / C)^B)^E
%
%       where:
%           A : Minimum asymptote
%           B : Hill slope (steepness and direction)
%           C : Inflection point (EC50)
%           D : Maximum asymptote
%           E : Asymmetry factor (E = 1 → symmetric 4PL)
%
%   Syntax:
%       x = L5Pinv(cf, y)
%
%   Inputs:
%       cf  - Either:
%              • A cfit object returned by L5P, containing parameters
%                A, B, C, D, and E of the 5PL model, or
%              • A numeric vector [A, B, C, D, E] (1x5 or 5x1).
%
%       y   - Numeric array of response values for which you want to
%             compute the corresponding x values. y must be real, finite,
%             non-empty. It may be a scalar, vector, or matrix; x will have
%             the same size as y.
%
%   Outputs:
%       x   - Numeric array of x values such that the 5PL model evaluated
%             at x returns (approximately) the values in y. The size of x
%             matches the size of y.
%
%   Example:
%       % Fit a 5PL curve using L5P:
%       xdata = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       ydata = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L5P(xdata, ydata);
%
%       % Invert the fitted curve at a specific response:
%       yq = 1.782315;
%       xq = L5Pinv(cf, yq);
%
%       % Alternatively, using an explicit parameter vector:
%       params = coeffvalues(cf);       % [A, B, C, D, E]
%       xq2   = L5Pinv(params, yq);     % should match xq
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
%       L5P

% -----------------------------
% Input parsing and validation
% -----------------------------

p = inputParser;

% cf must be either a fit-like object or a numeric vector of length 5.
addRequired(p, 'cf', @(v) ...
    isobject(v) || ...
    (isnumeric(v) && isvector(v) && numel(v) == 5 && ...
     all(isfinite(v(:)) & isreal(v(:)) & ~isnan(v(:)))));

% y must be numeric, real, finite, non-empty; any shape allowed.
addRequired(p, 'y', @(v) ...
    validateattributes(v, {'numeric'}, ...
                       {'real', 'finite', 'nonnan', 'nonempty'}));

parse(p, cf, y);
cf = p.Results.cf;
y  = p.Results.y;

clear p;

% -----------------------------
% Extract 5PL parameters (A, B, C, D, E)
% -----------------------------
% If cf is a fit object, extract coefficients; otherwise use numeric vector.

if isobject(cf)
    params = coeffvalues(cf); % expects [A, B, C, D, E]
else
    params = cf(:).';         % ensure row vector [A, B, C, D, E]
end

A = params(1);
B = params(2);
C = params(3);
D = params(4);
E = params(5);

% -----------------------------
% Basic sanity checks on y
% -----------------------------
% For a standard 5PL curve, the response typically lies between A and D.
% Values outside (min(A, D), max(A, D)) may still be computed, but can
% correspond to extrapolation or non-physical x values.

yMinModel = min(A, D);
yMaxModel = max(A, D);

if any(y(:) <= yMinModel)
    warning('L5Pinv:BelowModelRange', ...
        ['Some response values are <= the lower asymptote. ' ...
         'Inversion may be extrapolative or unreliable.']);
end

if any(y(:) >= yMaxModel)
    warning('L5Pinv:AboveModelRange', ...
        ['Some response values are >= the upper asymptote. ' ...
         'Inversion may be extrapolative or unreliable.']);
end

% -----------------------------
% Invert the 5PL equation
% -----------------------------
% From:
%   y = D + (A - D) / (1 + (x / C)^B)^E
%
% Rearranging:
%   y - D = (A - D) / (1 + (x / C)^B)^E
%   (A - D) / (y - D) = (1 + (x / C)^B)^E
%   ((A - D) / (y - D))^(1 / E) = 1 + (x / C)^B
%   (x / C)^B = ((A - D) / (y - D))^(1 / E) - 1
%   x / C     = ( ((A - D) / (y - D))^(1 / E) - 1 )^(1 / B)
%   x         = C * ( ((A - D) / (y - D))^(1 / E) - 1 ).^(1 / B)
%
% The computation is performed element-wise over y.

ratio = (A - D) ./ (y - D);
inner = ratio .^ (1 ./ E) - 1;
x     = C .* (inner .^ (1 ./ B));

end
