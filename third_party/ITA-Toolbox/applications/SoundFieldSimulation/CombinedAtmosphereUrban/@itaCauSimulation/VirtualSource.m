function [virtualSource] = VirtualSource(obj, atmoPath, receiverPosition, delay)
%Calculates the virtual source position based on an atmospheric path and a
%receiver Position.
%   Optionally, the propagation delay can be passed as third parameter
%   (e.g. when smoothing the delay in a dynamic scenario).

if nargin < 4
    delay = atmoPath.t(end);
end

distance = delay * obj.atmosphere.c(receiverPosition(3));
switch (obj.virtualSourceModus)
    case 'delay'
    case 'spreading_loss'
        try
            distance = 1 / atmoRay.spreadingLoss;
        catch
            disp('Warning: Ray Area at receiver could not be determined! Switched to propagation delay mode.')
        end
    otherwise
        disp('Warning: Unknown Virtual Source mode! Switched to propagation delay mode.')
end

wfNormalAtReceiver = atmoPath.n.cart(end,:);
% Make sure virtual source doesn't end up below ground
if wfNormalAtReceiver(3) > 0
    wfNormalAtReceiver(3) = 0;
end
wfNormalAtReceiver = wfNormalAtReceiver / norm(wfNormalAtReceiver);


virtualSource = receiverPosition - wfNormalAtReceiver * distance;
%Set minimum height to 0
if virtualSource(3) < 0
    warning('Height of virtual source corrected to minimum of 0 m')
    virtualSource(3) = 0;
end

