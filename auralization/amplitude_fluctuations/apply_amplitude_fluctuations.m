function [y, info] = apply_amplitude_fluctuations(x, fs, t_ctrl, f_ctrl, sigma_ctrl, opts)
% APPLY_AMPLITUDE_FLUCTUATIONS  Applies turbulence-induced amplitude
% (level) fluctuations to a receiver signal as a smooth, time-variant,
% zero-phase gain in the short-time Fourier transform (STFT) domain.
%
%   [y, info] = apply_amplitude_fluctuations(x, fs, t_ctrl, f_ctrl, sigma_ctrl, opts)
%
% For each STFT frame m and frequency bin f, the level is changed by
%
%     dL(f, m) = sigma(f, t_m) * u(t_m) - (ln(10)/20) * sigma(f, t_m)^2      [dB]
%
% where sigma(f, t) is the standard deviation of the level fluctuations
% (e.g. from get_log_amplitude_fluctuation) and u(t) a unit-variance random
% process (get_amplitude_modulation). The second term is the energy
% correction of Lincke, Kawai & Pieren (2025), Eq. (16), ln(10)/20 = 0.115:
% for Gaussian dL it gives E[10^(dL/10)] = 1, i.e. the AVERAGE ENERGY IS
% UNCHANGED at every frequency.
%
% A single process u(t) is used for all frequencies (fully correlated level
% fluctuations across frequency, as in all approaches in the literature) and
% for all channels (binaural: same gain at both ears). Applied to the summed
% direct + ground-reflected signal, this is the same as applying identical
% fluctuations to both paths, which leaves the expected cross term, and thus
% the coherence factor of the coherence-loss filters, unchanged (no double
% counting of decorrelation).
%
% The gain changes from frame to frame (hop = nfft/4, 10.7 ms at 48 kHz and
% nfft = 2048) and is cross-faded by the overlapping windows, independent of
% the 0.1 s FIR blocks of the propagation filters. Because u(t) has a
% correlation time of ~0.3 s, there are no switching artifacts.
%
% INPUTS
%   x          : [N x nCh] signal (e.g. mono or binaural receiver signal)
%   fs         : sampling frequency [Hz]
%   t_ctrl     : [nT x 1] times [s] of the rows of sigma_ctrl, on the time
%                axis of x (t = 0 at the first sample)
%   f_ctrl     : [1 x nF] frequencies [Hz] of the columns of sigma_ctrl
%   sigma_ctrl : [nT x nF] standard deviation of the level fluctuations [dB]
%                (interpolated linearly in time and in log-frequency, held
%                constant outside the given ranges)
%   opts       : (optional) struct
%       .nfft  STFT size, default 2048 (power of 2)
%       .fc    cutoff frequency of the modulation low-pass [Hz], default 0.5 (LKP25)
%       .seed  seed of the random modulation, default 1
%
% OUTPUTS
%   y    : [N x nCh] signal with amplitude fluctuations
%   info : struct with
%       .t_frame  frame centre times [s]
%       .u        modulation process at the frames (unit variance)
%       .f_bin    bin frequencies [Hz]
%       .dL       applied level change [nBins x nFrames] (dB), only if
%                 requested with opts.keep_dL = true (memory!)
%
% NOTES
%   * Analysis/synthesis with a square-root periodic Hann window and 75 %
%     overlap (weighted overlap-add). With sigma = 0 the signal is
%     reconstructed exactly (up to rounding).
%   * The STFT parameters are a framework choice, not taken from the
%     literature. Lincke et al. (2025) used a time-variant second-order IIR
%     shelving filter fitted to sigma(f) instead.
%
% Reference
%   D. Lincke, C. Kawai, R. Pieren, Acta Acustica 9, 69 (2025). doi:10.1051/aacus/2025053
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(5, 6);
if nargin < 6, opts = struct; end
if ~isfield(opts, 'nfft'),    opts.nfft = 2048;   end
if ~isfield(opts, 'fc'),      opts.fc = 0.5;      end
if ~isfield(opts, 'seed'),    opts.seed = 1;      end
if ~isfield(opts, 'keep_dL'), opts.keep_dL = false; end

nfft = opts.nfft;
hop  = nfft/4;
nBins = nfft/2 + 1;
[N, nCh] = size(x);

% square-root periodic Hann window: sum of w^2 over the 4 overlapping frames = 2
n   = (0:nfft-1).';
win = sqrt(0.5 - 0.5*cos(2*pi*n/nfft));
norm_ola = 2;

% zero padding at both ends, so that every sample is covered by 4 frames
pad   = nfft;
nPad  = N + 2*pad;
nFrames = floor((nPad - nfft)/hop) + 1;
nPad  = (nFrames - 1)*hop + nfft;
xp    = zeros(nPad, nCh);
xp(pad+1:pad+N, :) = x;

t_frame = ((0:nFrames-1).'*hop + nfft/2 - pad)/fs;     % frame centre on the time axis of x
f_bin   = (0:nBins-1).' * fs/nfft;

% random modulation (seeded, reproducible)
rs = RandStream('mrg32k3a', 'Seed', opts.seed);
u  = get_amplitude_modulation(nFrames, hop/fs, opts.fc, rs);

% sigma on the frame/bin grid: time (linear, held at the ends), then log-frequency
t_ctrl = t_ctrl(:);
f_ctrl = f_ctrl(:).';
if numel(t_ctrl) == 1
    sig_t = repmat(sigma_ctrl(1,:), nFrames, 1);
else
    tq    = min(max(t_frame, t_ctrl(1)), t_ctrl(end));
    sig_t = interp1(t_ctrl, sigma_ctrl, tq, 'linear');            % [nFrames x nF]
end
if numel(f_ctrl) == 1
    sig_tf = repmat(sig_t, 1, nBins);
else
    fq     = log(min(max(f_bin, f_ctrl(1)), f_ctrl(end)));
    sig_tf = interp1(log(f_ctrl).', sig_t.', fq, 'linear').';     % [nFrames x nBins]
end

% level change per frame and bin, Lincke et al. (2025) Eqs. (15)-(16)
dL = sig_tf .* u - (log(10)/20) * sig_tf.^2;                       % [nFrames x nBins]
G  = 10.^(dL/20);

% weighted overlap-add
yp = zeros(nPad, nCh);
for m = 1:nFrames
    idx = (m-1)*hop + (1:nfft);
    X   = fft(xp(idx, :) .* win);
    Xh  = X(1:nBins, :) .* G(m, :).';
    seg = real(ifft([Xh; conj(Xh(end-1:-1:2, :))]));
    yp(idx, :) = yp(idx, :) + seg .* win;
end
y = yp(pad+1:pad+N, :) / norm_ola;

info.t_frame = t_frame;
info.u       = u;
info.f_bin   = f_bin;
if opts.keep_dL
    info.dL = dL.';
end

end
