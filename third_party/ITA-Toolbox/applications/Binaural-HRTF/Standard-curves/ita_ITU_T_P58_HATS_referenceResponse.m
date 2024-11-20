function res = ita_ITU_T_P58_HATS_referenceResponse(fieldtype,plot)
% get the ITU P.58-2021 definition for artificial head transfer functions
%
% Arguments:
%   fieldtype: ('free') 'diffuse' Normalization definition 
%   plot:      (true)

% based on Recommendation P.58 (06/21) https://www.itu.int/rec/T-REC-P.58-202106-I/en
arguments
    fieldtype = 'free'
    plot = false
end

if nargout < 1
    plot = true;
end

switch fieldtype
    case {'free','freeField','freefield','free-field'}
        responseName = 'free-field response';
        freq = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000]';
        data = [0 0 0 0 .5 1 1.5 2 2.5 3.5 3.5 3.5 5 12.5 18.5 15.5 13 11 5 2 7 9 5.5 -3]';
        tol  = [-1.5 -1.5 -1.5 -1.5 -2 -2 -2 -1.5 -1.5 -1.5 -2 -3 -3.5 -4 -4.5 -2.5 -2.5 -3.5 -3 -5 -7.5 -6 -5 -5.5;...
                1.5 1.5 1.5 1.5 1.5 1.5 1.5 2 2 3 3.5 3 2.5 1.5 1.5 5.5 5 5 8.5 9 3.5 5 8 12.5]';
    case {'diffuse'}
        responseName = 'diffuse-field response';
        freq = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000]';
        data = [0 0 0 0 .5 .5 1 1.5 2 4 5 6.5 8 10.5 14 12 11.5 11 8 6.5 10.5 4 1 -3]';
        tol  = [-1.5 -1.5 -1.5 -1.5 -1.5 -1 -1 -1.5 -1.5 -2 -2.5 -3 -2.5 -2.5 -3 -1 -2 -2 -4 -4 -10 -3 -4 -5;...
                1 1 1 1 1 1.5 1.5 1.5 2 2 2 1.5 1.5 2 2 6 6 5 6.5 8 2 6.5 8.5 9.5]';
    otherwise
            error('Unknown field type argument');
end


response = itaResult;
response.freqVector = freq;
response.freqData = 10.^(data/20);
response.channelNames = {responseName};
response.plotLineProperties = {'LineStyle','-','Color',[1 1 1]*0};

limits = itaResult;
limits.freqVector = freq;
limits.freqData = 10.^((data+tol)/20);
limits.channelNames = {'lower limit','upper limit'};
limits.plotLineProperties = {'LineStyle','--','Color',[1 1 1]*0.5};


if plot
    f = ita_plot_freq(response);
    hold on
    ita_plot_freq(limits,'figure_handle',f,'hold',true)
end

res = merge(response,limits);

end