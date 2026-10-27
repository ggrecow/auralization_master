function [B, bhat, Ccoh] = get_decorrelation_filter(w, nfft, fs, rngStream, opts)
% GET_DECORRELATION_FILTER  Random-phase partial decorrelation filter for the
% ground-reflected sound, following Pieren & Lincke (2022), Sec. III.
%
%   [B, bhat, Ccoh] = get_decorrelation_filter(w, nfft, fs, rngStream, opts)
%
% Designs one realization of a random-phase all-pass filter whose
% expectation reproduces the target coherence factor
%
%     Ccoh(f) = exp(-w * f^2),     i.e.  E[B(f)] = Ccoh(f)
%
% so that direct + filtered reflected sound gives, in the mean, the
% partially coherent ground effect (Pieren & Lincke Eq. 3, Ostashev Eq. 22).
%
% Design (Kendall 1995, extended by Pieren & Lincke):
%     B(f) = exp( j * beta(f) * phi ),  phi ~ U(-pi, pi),  0 <= beta <= 1
% with Ccoh = sinc(beta) = sin(pi*beta)/(pi*beta).
%
% Each call draws one random phase vector phi from rngStream; the function
% itself has no notion of time. The update rate is set by the caller:
% get_FIR assigns one RandStream substream per filter epoch (0.1 s, epochs
% defined in apply_propagation, following Pieren & Lincke 2022), so all
% time steps within an epoch share the same phi, while beta still follows
% the w of each time step.
%
% INPUTS
%   w         : coherence parameter, Ccoh = exp(-w*f^2). w <= 0 -> B = 1.
%   nfft      : double-sided FFT length of the transfer functions
%               (numFreqBins_double_sided in get_FIR, 2^14 by default)
%   fs        : sampling frequency [Hz]
%   rngStream : (optional) RandStream the random phase is drawn from. Set its
%               substream before each call to control which filters share a
%               phase (see get_FIR). Default: global stream.
%   opts      : (optional) struct
%       .N                  filter taps, default 2^11 (43 ms @ 48 kHz)
%       .tukey              Tukey cosine fraction, default 0.25
%       .betaApprox         1 = power law (Pieren Eq. 21)
%                           2 = raised cosine (Pieren Eq. 22)
%                           3 = exact inverse of sinc (DEFAULT)
%                           Use 1 to reproduce Pieren & Lincke exactly.
%       .normalizeMagnitude false (default) keeps the paper's design;
%                           true forces |B| = 1 after windowing.
%
% OUTPUTS
%   B    : [nfft/2+1 x 1] single-sided spectrum, zero-delay (zero-centred
%          filter), ready to multiply with the reflected-path TF
%   bhat : [N x 1] windowed, centred impulse response (for plotting)
%   Ccoh : [nfft/2+1 x 1] target coherence factor on the TF grid
%
% NOTES
%   * Zero-centred placement in the nfft buffer compensates the N/2 filter
%     delay; otherwise the comb would shift towards lower frequencies.
%   * The Tukey window lowers |B| to ~0.85 (-1.4 dB) where Ccoh -> 0. On the
%     total (direct + reflected, |Q| = 1) this is at most ~ -0.6 dB.
%   * The central peak of bhat is NOT reliably the maximum of |bhat| when
%     the decorrelation is strong. Do not locate the reflected-path delay
%     on the decorrelated IR; see get_FIR.
%
% References
%   R. Pieren, D. Lincke, JASA 151(4), 2453-2460 (2022).
%   G. Kendall, Computer Music Journal 19(4), 71-87 (1995).
%
% Author: Gil Felix Greco, Braunschweig (27.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(3, 5);
if nargin < 4 || isempty(rngStream), rngStream = RandStream.getGlobalStream; end
if nargin < 5, opts = struct; end

if ~isfield(opts, 'N'),                  opts.N = 2^11;               end
if ~isfield(opts, 'tukey'),              opts.tukey = 0.25;           end
if ~isfield(opts, 'betaApprox'),         opts.betaApprox = 3;         end
if ~isfield(opts, 'normalizeMagnitude'), opts.normalizeMagnitude = false; end

N     = opts.N;
nBins = nfft/2 + 1;
f     = (0:nBins-1).' * fs/nfft;
Ccoh  = exp(-w * f.^2);

assert(mod(N,2) == 0 && N < nfft, 'opts.N must be even and smaller than nfft.')

% ---- trivial case: no coherence loss -> unit filter ----
if w <= 0
    B = ones(nBins, 1);
    bhat = zeros(N,1); bhat(N/2+1) = 1;
    return
end

% ---- design grid of the filter ----
fd    = (0:N/2).' * fs/N;
Cd    = exp(-w * fd.^2);
beta  = il_beta_from_coherence(Cd, opts.betaApprox);

% ---- random-phase half spectrum, Eq. (17) ----
phi        = (rand(rngStream, N/2+1, 1) - 0.5) * 2*pi;   % U(-pi, pi)
H          = exp(1j * beta .* phi);
H(1)       = 1;      % DC must be real
H(end)     = 1;      % Nyquist must be real

% ---- Hermitian spectrum and IDFT, Eqs. (23)-(25) ----
Bfull = [H; conj(flipud(H(2:end-1)))];
b     = ifft(Bfull, 'symmetric');

% ---- centre and window, Eq. (26) ----
bhat = circshift(b, N/2) .* tukeywin(N, opts.tukey);
c    = N/2 + 1;                          % index of the central sample

% ---- place zero-centred in the nfft buffer (delay compensation) ----
bpad              = zeros(nfft, 1);
bpad(1:N-c+1)     = bhat(c:end);         % causal half
bpad(end-c+2:end) = bhat(1:c-1);         % anti-causal half, wrapped

Bd = fft(bpad);
B  = Bd(1:nBins);

if opts.normalizeMagnitude
    B = B ./ max(abs(B), eps);
end

end

% =====================================================================
function beta = il_beta_from_coherence(Ccoh, approx)
% Inverse of Ccoh = sinc(beta), first lobe only (0 <= beta <= 1).

Ccoh = min(max(Ccoh, 0), 1);

switch approx
    case 1   % power law, Pieren & Lincke Eq. (21)
        beta = (1 - Ccoh).^(2/3);
    case 2   % raised cosine, Pieren & Lincke Eq. (22)
        beta = (2/pi) * acos(Ccoh.^(5/6));
    case 3   % exact, via lookup table of the first sinc lobe
        bt = linspace(0, 1, 2001).';
        Ct = sin(pi*bt) ./ (pi*bt);  Ct(1) = 1;  Ct(end) = 0;
        beta = interp1(flipud(Ct), flipud(bt), Ccoh, 'pchip');
    otherwise
        error('opts.betaApprox must be 1, 2 or 3.')
end

beta = min(max(beta, 0), 1);

end
