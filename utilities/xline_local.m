function h = xline_local(varargin)
%XLINE_COMPAT  Version-safe wrapper for xline.
%
% Behavior:
%   - MATLAB >= R2018b: calls native xline
%   - MATLAB <= R2018a: does nothing (no error)
%
% Usage:
%   xline_compat(...)
%   a = xline_compat(...)

if verLessThan('matlab','9.5')   % R2018b = 9.5
    % Old MATLAB: do nothing
    h = [];
    return;
end

% New MATLAB: call native xline
h = xline(varargin{:});
end
