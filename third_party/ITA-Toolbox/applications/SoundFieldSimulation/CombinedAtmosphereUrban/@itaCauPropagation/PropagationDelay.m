function delay = PropagationDelay(obj, atmoDelay, urbanPaths, constSpeedOfSound)
%PROPAGATIONDELAY Calculates corrected delay of a set of sound paths of a
%combined atmospheric urban simulation.
% All input variables are optional.

if nargin < 2; [~, atmoDelay] =  obj.calcAtmosphericParameters(); end
if nargin < 3; urbanPaths = obj.urbanPaths; end
if nargin < 4; constSpeedOfSound = obj.constSpeed; end

%% Derive parameters from paths
virtSourceDistance = norm( urbanPaths(1).propagation_anchors{end}.interaction_point -...
                           urbanPaths(1).propagation_anchors{1}.interaction_point );
urbanPathLengths = ita_propagation_path_length( urbanPaths );

%% Calculation of delay
delay = zeros(1, numel(urbanPaths));
for idPath = 1:numel(urbanPaths)
    delay(idPath) = adjust_delay(atmoDelay, virtSourceDistance, urbanPathLengths(idPath), constSpeedOfSound);
end

%% Equivalent complex linlear phase response
% if nargout > 1
%     phaseFactor = zeros(numel(freqVector), numel(urbanPaths));
%     phaseFactor(2:end, :) = exp( 1j*2*pi*delay.*freqVector(2:end) );
% end

%% Adjust delay
function delay = adjust_delay(atmoDelay, virtSourceDistance, pathLength, c)
%ITA_CAU_DELAY adjusts urban propagation delay based on atmospheric path

deltaDelay = atmoDelay - virtSourceDistance/c;
delay = deltaDelay + pathLength/c;
