function [g_inPhase] = ita_hoa_inPhaseWeights(N, varargin)
%ITA_HOA_INPHASEWEIGHTS Calculate in-phase weights
% References:
%   [1] Daniel; Représentation de champs acoustiques, application à la
%       transmission et à la reproduction de scènes sonores complexes dans un
%       contexte multimédia; 2001
%
% Arguments:
%   N: int, maximum spherical harmonics order
% Options: (default)
%   {}
%
% Returns:
%   g_inPhase: [], in-phase weights

%% Calculate Weights
g_inPhase = zeros((N+1)^2, 1);
% see [1] p.314
% g_0 = sqrt(N*(2*M+1)/(M+1)^2)
g_inPhase(1)=sqrt(N*(2*N+1)/(N+1)^2);
ndx = ita_sph_degreeorder2linear(0:N);
for n = 1:N
    g_n = factorial(N) * factorial(N+1) / factorial(N+n+1) / factorial(N-n);
    g_inPhase(ndx(n)+1:ndx(n+1)) = g_n;
end

end

