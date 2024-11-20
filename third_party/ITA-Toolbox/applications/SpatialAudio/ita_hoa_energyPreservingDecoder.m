function [D_epad] = ita_hoa_energyPreservingDecoder(LS_pos, varargin)
%ITA_HOA_SAMPLINGDECODER Energy Preserving Ambisonics Decoder (EPAD)
% The EPAD uses singular value decomposition to calculate the decoder
% matrix by dropping the singular values from the decomposition.
% Taken from:
% [1] Zotter et al.; Comparison of energy-preserving and all-round Ambisonic decoders; 2013
%
% Arguments:
%   LS_pos: loudspeaker setup / SH sampling positions
%   Options: (default)
%   {N           : SH truncation order, should not be greater than floor(sqrt(L)-1)
%    'real'      : (true)            , real valued SH
%                  false             , complex valued SH}
%
% Returns:
%   D_epad: EPAD matrix

%% Parse input
if mod(numel(varargin), 2)
    N = varargin{1};
    varargin(1) = [];
else
    N = floor(sqrt(LS_pos.nPoints)-1);
end
% Default values
params.real = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Get SH base functions sampled at LS_pos
Y = ita_sph_base(LS_pos, N, 'real', params.real, varargin{:}).';

%% Calculate Decoder Matrix
[U, S, V] = svd(Y, 'eco');
D_epad = V*U.' / max(S, [],'all');

end

