function [D_sad] = ita_hoa_samplingDecoder(LS_pos, varargin)
%ITA_HOA_SAMPLINGDECODER Sampling Ambisonic Decoding (SAD)
% The SAD uses an approximation of the Moore-Penrose pseudoinverse to 
% calculate the decoder matrix. It is the optimal decoder for t-designs.
% Taken from:
%   [1] Zotter, Frank; All-Round Ambisonic Panning and Decoding; 2012
%   [2] Zotter et al.; Comparison of energy-preserving and all-round Ambisonic decoders; 2013
%
% Arguments:
%   LS_pos: loudspeaker setup / SH sampling positions
% Options: (default)
%   {N        : SH truncation order, should not be greater than floor(sqrt(L)-1)
%    'real'   : (true)             , real valued SH
%               false              , complex valued SH}
%
% Returns:
%   D_sad: SAD matrix

%% Parse input
if mod(numel(varargin), 2)
    N = varargin{1};
    varargin(1) = [];
else
    N = floor(sqrt(LS_pos.nPoints)-1);
end
% Default values
params.norm = 'N3D';
params.real = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Get SH base functions sampled at LS_pos
Y = ita_sph_base(LS_pos, N, 'real', params.real, varargin{:}).';

%% Calculate Decoder Matrix
D_sad = Y.'/LS_pos.nPoints*4*pi;
% D_sad = Y.'/(N+1);

end

