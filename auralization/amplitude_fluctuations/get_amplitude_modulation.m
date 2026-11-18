function u = get_amplitude_modulation(nSamples, dt, fc, rngStream)
% GET_AMPLITUDE_MODULATION  Seeded, unit-variance random modulation signal for
% the turbulence-induced amplitude fluctuations.
%
%   u = get_amplitude_modulation(nSamples, dt, fc, rngStream)
%
% White Gaussian noise filtered by a first-order low-pass with cutoff
% frequency fc, as in Lincke, Kawai & Pieren (2025), Sec. 3.3 (fc = 0.5 Hz).
% Its power spectral density is S(f) ~ 1/(1 + (f/fc)^2), i.e. a
% Gauss-Markov (Ornstein-Uhlenbeck) process with correlation time 1/(2 pi fc).
% It is generated with the exact discrete-time recursion
%
%     u(n) = r u(n-1) + sqrt(1 - r^2) e(n),   r = exp(-2 pi fc dt),  e ~ N(0,1)
%
% started from the stationary distribution, u(1) ~ N(0,1). Hence u is
% stationary with zero mean and unit variance from the first sample on, for
% any sampling interval dt. The same rngStream state always gives the same u.
%
% INPUTS
%   nSamples  : number of samples
%   dt        : sampling interval [s] (e.g. the STFT hop size / fs)
%   fc        : cutoff frequency of the first-order low-pass [Hz] (LKP25: 0.5 Hz)
%   rngStream : (optional) RandStream; default: global stream
%
% OUTPUT
%   u         : [nSamples x 1] modulation signal, zero mean, unit variance
%
% Reference
%   D. Lincke, C. Kawai, R. Pieren, "Atmospheric turbulence affects noise
%   annoyance from aircraft flyovers", Acta Acustica 9, 69 (2025).
%   doi:10.1051/aacus/2025053
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(3, 4);
if nargin < 4 || isempty(rngStream), rngStream = RandStream.getGlobalStream; end

e = randn(rngStream, nSamples, 1);
r = exp(-2*pi*fc*dt);

u = zeros(nSamples, 1);
if nSamples == 0, return; end
u(1) = e(1);
u(2:end) = filter(sqrt(1 - r^2), [1 -r], e(2:end), r*e(1));   % same as the recursion above

end
