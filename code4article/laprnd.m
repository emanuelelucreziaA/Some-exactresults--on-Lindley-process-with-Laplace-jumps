function y  = laprnd(m, n, mu, sigma)
%LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution

%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07

%Check inputs
if nargin < 2
    error('At least two inputs are required');
end

if nargin == 2
    mu = 0; sigma = 1;
end

if nargin == 3
    sigma = 1;
end

% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2*abs(u));

%% Commento di Emanuele
%{
y = F_L^{-1}F_U(u)
dove F_L, F_U sono le distribuzioni di una Laplace e uniforme,
rispettivamente, con:
c = F_U(u) = u+0.5
e
y = F_L^{-1}(c) = 
  = mu + b*log(2*c) se c < 0.5
  = mu - b*log(2-2*c) se c > 0.5

%}
