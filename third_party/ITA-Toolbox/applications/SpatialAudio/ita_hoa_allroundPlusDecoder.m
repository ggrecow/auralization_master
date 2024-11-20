function [D_allradplus] = ita_hoa_allroundPlusDecoder(LS_pos, tDesign, varargin)
%ITA_HOA_ALLROUNDPLUSDECODER This decoder combines the decoder strategies 
%   of SAD and AllRAD into one decoder matrix by weighted superposition.
% Taken from:
%   [1] Zotter et al.; Comparison of energy-preserving and all-round Ambisonic decoders; 2013
%
% Arguments:
%   LS_pos : itaCoordinates, loudspeaker setup
%   tDesign: itaCoordinates, Ambisonics Reproduction setup
% Options: (default)
%   {N          : SH truncation order , should not be greater than floor(sqrt(L)-1)}
% 
% Returns:
%   D_allradplus: AllRAD+ matrix

%% Parse input
if mod(numel(varargin), 2)
    N = varargin{1};
    varargin(1) = [];
else
    N = floor(sqrt(LS_pos.nPoints)-1);
end
% Default values
params.decoder = 'SAD';
[params, ~] = ita_parse_arguments(params, varargin);
% Enforce SAD decoding
if ~strcmpi(params.decoder, 'sad')
    idx = find(lower(varargin) == 'sad');
    varargin{idx+1} = 'SAD';
end

%% Calculate full Decoder Matrix
% AllRAD Matrix
D_allrad = ita_hoa_allroundDecoder(LS_pos, tDesign, N, varargin{:});

% SAD Matrix
D_sad = ita_hoa_samplingDecoder(LS_pos, N, varargin{:});

% calculate AllRAD+ Decoder
D_allradplus = D_allrad/sqrt(2) + D_sad/sqrt(8);

end
