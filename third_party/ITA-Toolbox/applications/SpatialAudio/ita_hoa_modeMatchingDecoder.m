function [D_mmad] = ita_hoa_modeMatchingDecoder(LS_pos, varargin)
%ITA_HOA_SAMPLINGDECODER Mode Matching Ambisonic Decoding (MMAD)
% The MMAD uses the Moore-Penrose pseudoinverse to calculate the decoder
% matrix.
% Taken from:
% [1] Poletti; Three-Dimensional Surround Sound Systems Based on Spherical Harmonics; 2005
% [2] Zotter et al.; Comparison of energy-preserving and all-round Ambisonic decoders; 2013
%
% Arguments:
%   LS_pos: loudspeaker setup / SH sampling positions
%   Options: (default)
%   {N           : SH truncation order, should not be greater than floor(sqrt(L)-1)
%    'definition': 'Poletti' , (Y.'*Y + lambda*I)^(-1) * Y.'
%                  'Zotter') , Y.' * (Y*Y.' + lambda*I)^(-1)
%                  ('svd')   , singular value decomposition
%    'lambda'    : (L/(4*pi)), regularization factor, L = #loudspeaker
%    'real'      : (true)    , real valued SH
%                  false     , complex valued SH
%
% Returns:
%   D_mmad: MMAD matrix

%% Parse input
if mod(numel(varargin), 2)
    N = varargin{1};
    varargin(1) = [];
else
    N = floor(sqrt(LS_pos.nPoints)-1);
end
% Default values
params.definition = 'svd';
params.lambda = LS_pos.nPoints / (4*pi);
params.real = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Get SH base functions sampled at LS_pos
Y = ita_sph_base(LS_pos, N, 'real', params.real, varargin{:}).';

%% Calculate Decoder Matrix
switch lower(params.definition)
    case {'poletti', 'p'}
        D_mmad = ((Y.'*Y + params.lambda*eye(LS_pos.nPoints)) \ Y.');
    case {'zotter', 'z'}
        D_mmad = (Y.' / (Y*Y.' + params.lambda*eye((N+1)^2)));
    case 'svd'
        D_mmad = pinv(Y);
    otherwise
        error('I do not know that MMAD definition');
end

end

