function h = yline_local(varargin)
%YLINE_COMPAT  Version-safe wrapper for yline.
%
% Behavior:
%   - MATLAB >= R2018b: calls native yline
%   - MATLAB <= R2018a: does nothing (no error)
%
% Usage:
%   yline_compat(...)
%   a = yline_compat(...)

if verLessThan('matlab','9.5')   % R2018b = 9.5
    % Old MATLAB: do nothing
    h = [];
    return;
end

% New MATLAB: call native yline
h = yline(varargin{:});
end
