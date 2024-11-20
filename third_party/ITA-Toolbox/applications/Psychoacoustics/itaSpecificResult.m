classdef itaSpecificResult

    %   ITASPECIFICRESULT - class for specific results from psychoacoustic
    %   parameters with 2-dimesional data (typ. bark-bin over time)
    %
    % This class should be used for the output of ita-toolbox functions
    % that provide data with time and frequency information.
    %
    %   %   Reference page in Help browser
    %       <a href="matlab:doc itaResult">doc itaResult</a>

    % <ITA-Toolbox>
    % This file is part of the ITA-Toolbox. Some rights reserved.
    % You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder.
    % </ITA-Toolbox>

    % Author: HBR & JSE
    % Created:  22-Sept-2022

    properties(Access = private)
        % Internal fields, no access from outside the class
        mTimeVector %values of x-Axis
        mFreqVector %values of y-Axis
        mData       %[X,Y,Channel]
        mChannelUnits = {};
        mChannelDataType = {}; %e.g. {'Specific Loudness','Specific Roughness'}
    end

    properties
        xUnit = 's';
        yUnit = 'Bark'; %or 'Bark'
    end


    properties(Dependent = true, Hidden = false)
        data
        timeVector %time Vector (x-Axis)
        freqVector %frequency Vector (y-Axis)
        nChannels   % how many data Channels
        nSamples    % numebr of time-Axis samples
        nBins       % number of frequency/bark axis samples
        channelUnits %Unit per Channel
        channelDataType %e.g. {'Specific Loudness','Specific Roughness'}
    end

    methods
        function this = itaSpecificResult(data,timeVector,freqVector)
            % Constructor
            %   itaSpecificResult(dataMatrix,timeVector,freqVector)
            this.data = data;
            this.timeVector = timeVector;
            this.freqVector = freqVector;
        end

        % to split channels for plot
        function result = ch(this, iChannel)
                result = this.mData(:,:,iChannel);
        end

        %% Get/set Stuff
        function result = get.data(this)
            result = this.mData;
        end

        function this = set.data(this, Data)
            this.mData = Data;
        end

        function result = get.timeVector(this)
            result = this.mTimeVector;
        end

        function this = set.timeVector(this,timeVector)
            %ToDo: add check for matching dimensions
            this.mTimeVector = timeVector(:);
        end

        function result = get.freqVector(this)
            result = this.mFreqVector;
        end

        function this = set.freqVector(this,freqVector)
            %ToDo: add check for matching dimensions
            this.mFreqVector = freqVector(:);
        end

        function result = get.channelUnits(this)
            result = this.mChannelUnits;
        end 

        function this = set.channelUnits(this, channelUnits)
            this.mChannelUnits = channelUnits;
        end 

         function result = get.channelDataType(this)
            result = this.mChannelDataType;
         end

         function this = set.channelDataType(this, channelDataType)
             this.mChannelDataType = channelDataType;
         end

        function result = get.nSamples(this)
            result = size(this.mData, 1);
        end

        function result = get.nBins(this)
            result = size(this.mData, 2);
        end

        function result = get.nChannels(this)
            result = size(this.mData, 3);
        end

     
        %% Plot
        function plot_specific(this)
            %figure out best layout
            subPLotLayout = this.bestSubplotLayout(this.nChannels);

            figure
            for iChannel = 1:this.nChannels
                subplot(subPLotLayout(1),subPLotLayout(2),iChannel)
                [X,Y] = meshgrid(this.timeVector, this.freqVector);
                plotData = this.ch(iChannel);
                p1 = pcolor(X,Y,plotData.');
                p1.FaceColor = 'interp';
                p1.EdgeColor = 'none';
                xlabel(sprintf('Time in [%s]', this.xUnit));
                ylabel(sprintf('Frequency in [%s]',this.yUnit));
                p2 = colorbar;
                
               if ~isempty(this.channelDataType)
                    title(sprintf('%s in %s',this.channelDataType{iChannel},this.channelUnits{iChannel}));
                     ylabel(p2,sprintf('%s in [%s]', this.channelDataType{iChannel}, this.channelUnits{iChannel}),'FontSize',10,'Rotation',270);
                    p2.Label.Position(1) = 4;
               else 
                   ita_verbose_info('Please define the channelDataType and the channelUnits to create the plot`s title and legend.')
               end 
                if strcmp(this.xUnit,'Hz')
                    [fTick, fTickLabel] = ita_plottools_ticks('log');
                    set(gca, 'YScale','log', 'YTick', fTick, 'YTickLabel', fTickLabel);
                end
            end
        end


    end

    %% internal
    methods(Access=private,Hidden=true)
        function subPlotLayout = bestSubplotLayout(this, nPlots)
            %figure out the best way to arrange n plots
            subPlotLayout = [floor(sqrt(nPlots)),ceil(sqrt(nPlots))];
            if subPlotLayout(1)*subPlotLayout(2) < this.nChannels
                subPlotLayout = [ceil(sqrt(nPlots)),ceil(sqrt(nPlots))];
            end
        end
    end
    %% static methods
    methods(Static, Hidden = true)
        function this = loadobj(sObj)
            % Called when an object is loaded
            superclass = false;
            % change mpo: not relying on svn properties anymore
            try
                if sObj.classrevision > 2600 ...
                        || isnan(sObj.classrevision) % for non-SVN clients
                    superclass = true;
                end
            catch
                % all right: it was NO superclass
            end

            if ~superclass
                if isfield(sObj,'header') %has header
                    this = ita_import_old(sObj);
                    return
                else    % is headerless
                    sObj.dimensions = sObj.dims;
                    sObj.channelCoordinates = itaCoordinates(sObj.channelcoordinates);
                    sObj.channelOrientation = itaCoordinates(sObj.channelorientation);
                    sObj = rmfield(sObj,{'dims','channelcoordinates','channelorientation'});
                end
            end

            % change mpo: not relying on svn properties anymore
            try
                sObj = rmfield(sObj,{'classrevision', 'classname'}); % Remove these fields cause no more needed, but we could use them for special case handling (e.g. old versions)
            catch
                % fields were not there obviously
            end

            % change mpo: not relying on svn properties anymore
            if isstruct(sObj)
                this = itaSpecificResult(sObj); % Just call constructor, he will take care
            else
                this = sObj;
            end
        end

        function revision = classrevision
            % Return last revision on which the class definition has been changed (will be set automatic by svn)
            rev_str = '$Revision: 12902 $'; % Please dont change this, will be set by svn
            revision = str2double(rev_str(isstrprop(rev_str,'digit')));
        end

        function result = propertiesSaved
            result = {'abscissa', 'resultType'};
        end

    end

end