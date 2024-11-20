function [gitCommitID, currentBranch] = ita_git_commit_id(varargin)
%ITA_GIT_COMMIT_ID - Returns the ID of your current local git commit of the
%ITA Toolbox.
%  Optionally, the current brancht is returned as second output argument.
%  Returns an empty string (''), if ID / branch cannot be determined.
%
%  Syntax:
%   [gitCommitID, currentBranch] = ita_git_commit_id()
%

%   Reference page in Help browser
%      <a href="matlab:doc ita_git_commit_id">doc ita_git_commit_id</a>
%
% Autor: Philipp Schaefer -- Email: psc@akustik.rwth-aachen.de
% Created:  15 Apr 2021

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved.
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder.
% </ITA-Toolbox>

%% Defaults
gitCommitID = '';
currentBranch = '';

%% Find git folder
gitFolder = fullfile(ita_toolbox_path, '.git');
if ~exist(gitFolder, 'dir')
    ita_verbose_info('No git folder found. Did you download the ITA Toolbox instead of cloning it?')
    return;
end

%% Check HEAD
headFile = fullfile(gitFolder, 'HEAD');
if ~exist(headFile, 'file')
    ita_verbose_info('Cannot find git head!')
    return;
end

txt = txt_fileread(headFile);
numLines = numel(txt);
if numLines == 0; return; end

bDetachedHead = ~strcmp(txt{1}, 'ref:');
bValidBranch = ~bDetachedHead && numLines >= 2;

%% Get commit ID of detached HEAD
if(bDetachedHead)
    gitCommitID = txt{1};
end

%% Get current commit ID from branch
if ~bValidBranch; return; end

branchFile = fullfile(gitFolder, txt{2});
currentBranch = strrep(txt{2}, 'refs/heads/', '');
if ~exist(branchFile, 'file')
    ita_verbose_info('Cannot find commit of current branch!')
    return;
end
txt = txt_fileread(branchFile);
if isempty(txt); return; end

gitCommitID = txt{1};


function txt = txt_fileread(filename)
fileID = fopen(filename);
try
    txt = textscan(fileID, '%s');
    txt = txt{1};
catch err
    fclose(fileID);
    rethrow(err);
end
fclose(fileID);