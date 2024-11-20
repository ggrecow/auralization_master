function [g_maxrE] = ita_hoa_rEmaxWeights(N, varargin)
%ITA_HOA_REMAXWEIGHTS Get max rE weighths
% References:
%   [1] Daniel; Représentation de champs acoustiques, application à la
%       transmission et à la reproduction de scènes sonores complexes dans un
%       contexte multimédia; 2001
%   [2] Zotter, Frank; All-Round Ambisonic Panning and Decoding; 2012
%
% Arguments:
%   N: int, maximum spherical harmonics order
% Options: (default)
%   {'definition': ('Daniel'), use rE max defintion from [1]
%                  'Zotter'  , or from [2]}
%
% Returns:
%   g_maxrE: [], max rE weights

%% Parse input
params.definition = 'Daniel';
[params, ~] = ita_parse_arguments(params, varargin);

%% Calculate Weights
g_maxrE = zeros((N+1)^2, 1);
switch lower(params.definition)
    case {'daniel', 'd'}
        % see [1] p.312, A.67/A.68
        % max_rE is largest root of P_(N+1)
        syms re
        % Legendre Polynom n = N+1, m = 0
        P = symfun(1/2^(N+1)/factorial(N+1)*diff(((re^2-1)^(N+1)),(N+1)), re);
        % find maximum root(Nullstelle)
        re = double(max(vpasolve(P)));
        
    case {'zotter', 'z', 'allrad'}
        % see [2] eq.(10)
        re = cosd(137.9 / (N+1.51));
        
    otherwise
        error('I do not know that max rE definiton!');
end

% g_n=P_n(max_rE)
ndx = ita_sph_degreeorder2linear(0:N);
for n = 0:N
    P = legendre(n, abs(re));
    if n == 0
        g_maxrE(ndx(n+1)) = P(1);
    else
        g_maxrE(ndx(n)+1:ndx(n+1)) = P(1);
    end
end

end

