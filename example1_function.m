clear; clc;

% Example 1: corrected scalar objectives and their derivatives.
% ff{k}: scalar objective; df{k}: derivative of ff{k}.
% f{i}: two-variable local objective; grad_f{i}: local gradient.

ff = cell(200,1);
df = cell(200,1);

ff{1}  = @(x) sqrt(x.^4+3) ./ 5 - 7 .* sin(x).^2 ./ 10;
ff{2}  = @(x) 2 .* sin(x) - (x.^2+2).^(1/3) ./ 10;
ff{3}  = @(x) 3 .* x.^2 ./ (10 .* sqrt(x.^2+1));
ff{4}  = @(x) -sin(x) - sqrt(x.^4+3) ./ 20;
ff{5}  = @(x) 2 .* sin(x).^2 - x.^2 ./ (5 .* sqrt(x.^2+1));
ff{6}  = @(x) -x.^2 ./ (10 .* sqrt(x.^2+1)) - sqrt(x.^4+3) ./ 10;
ff{7}  = @(x) -sin(x);
ff{8}  = @(x) x.^2 - 3 .* sin(x).^2 ./ 10;
ff{9}  = @(x) 2 .* sin(x).^2 + (x.^2+2).^(1/3) ./ 5;
ff{10} = @(x) -(x.^2+2).^(1/3) ./ 10;
ff{11} = @(x) sin(x).^2;
ff{12} = @(x) (x.^4-2 .* x.^2) ./ (1+x.^2);
ff{13} = @(x) (x.^2 ./ (1+x.^2)) .* sin(x).^2;
ff{14} = @(x) -1 ./ (1+x.^2);
ff{15} = @(x) (log(1+x.^2)-x.^2).^2 ./ (1+x.^2);
ff{16} = @(x) log(1+x.^2).^2;
ff{17} = @(x) (1-exp(-x.^2)).^2;
ff{18} = @(x) (x-sin(x)).^2 ./ (1+x.^2);
ff{19} = @(x) sin(x).^2 ./ (1+x.^2);
ff{20} = @(x) (cos(x)-1).^2;

df{1}  = @(x) 2 .* x.^3 ./ (5 .* sqrt(x.^4+3)) - 7 .* cos(x) .* sin(x) ./ 5;
df{2}  = @(x) 2 .* cos(x) - x ./ (15 .* (x.^2+2).^(2/3));
df{3}  = @(x) 3 .* x ./ (5 .* sqrt(x.^2+1)) - 3 .* x.^3 ./ (10 .* (x.^2+1).^(3/2));
df{4}  = @(x) -cos(x) - x.^3 ./ (10 .* sqrt(x.^4+3));
df{5}  = @(x) 4 .* cos(x) .* sin(x) - 2 .* x ./ (5 .* sqrt(x.^2+1)) + x.^3 ./ (5 .* (x.^2+1).^(3/2));
df{6}  = @(x) x.^3 ./ (10 .* (x.^2+1).^(3/2)) - x ./ (5 .* sqrt(x.^2+1)) - x.^3 ./ (5 .* sqrt(x.^4+3));
df{7}  = @(x) -cos(x);
df{8}  = @(x) 2 .* x - 3 .* cos(x) .* sin(x) ./ 5;
df{9}  = @(x) 4 .* cos(x) .* sin(x) + 2 .* x ./ (15 .* (x.^2+2).^(2/3));
df{10} = @(x) -x ./ (15 .* (x.^2+2).^(2/3));
df{11} = @(x) 2 .* sin(x) .* cos(x);
df{12} = @(x) 2 .* x - 6 .* x ./ (1+x.^2).^2;
df{13} = @(x) 2 .* x .* sin(x).^2 ./ (1+x.^2).^2 ...
             + 2 .* (x.^2 ./ (1+x.^2)) .* sin(x) .* cos(x);
df{14} = @(x) 2 .* x ./ (1+x.^2).^2;
df{15} = @(x) 2 .* (log(1+x.^2)-x.^2) ...
             .* (2 .* x ./ (1+x.^2)-2 .* x) ./ (1+x.^2) ...
             - 2 .* x .* ((log(1+x.^2)-x.^2) ./ (1+x.^2)).^2;
df{16} = @(x) 4 .* x .* log(1+x.^2) ./ (1+x.^2);
df{17} = @(x) 4 .* x .* (1-exp(-x.^2)) .* exp(-x.^2);
df{18} = @(x) 2 .* (x-sin(x)) .* (1-cos(x)) ./ (1+x.^2) ...
             - 2 .* x .* ((x-sin(x)) ./ (1+x.^2)).^2;
df{19} = @(x) (2 .* sin(x) .* cos(x) .* (1+x.^2)-2 .* x .* sin(x).^2) ./ (1+x.^2).^2;
df{20} = @(x) -2 .* (cos(x)-1) .* sin(x);

% Keep the 20 base terms fixed when expanding to 200 scalar terms.
base_ff = ff(1:20);
base_df = df(1:20);
for k = 1:200
    base_idx = mod(k-1, 20) + 1;
    old_f    = base_ff{base_idx};
    old_df   = base_df{base_idx};
    E        = @(x) exp(-(sin(x).^k) ./ k);
    dE       = @(x) -E(x) .* cos(x) .* sin(x).^(k-1);
    ff{k} = @(x) old_f(x) + E(x);
    df{k} = @(x) old_df(x) + dE(x);
end

f = cell(100,1);
grad_f = cell(100,1);
for i = 1:100
    fL = ff{2*i-1};
    fR = ff{2*i};
    dL = df{2*i-1};
    dR = df{2*i};
    f{i} = @(x1,x2) fL(x1) + fR(x2);
    grad_f{i} = @(x1,x2) [dL(x1); dR(x2)];
end

% Complete objective and gradient for x = [x1; x2].
objective = @(x) sum(cellfun(@(fun) fun(x(1)), ff(1:2:end))) ...
              + sum(cellfun(@(fun) fun(x(2)), ff(2:2:end)));
gradient = @(x) [sum(cellfun(@(fun) fun(x(1)), df(1:2:end))); ...
                 sum(cellfun(@(fun) fun(x(2)), df(2:2:end)))];

% Derivative alias for callers using the original corrected library.
df_funcs = df;
