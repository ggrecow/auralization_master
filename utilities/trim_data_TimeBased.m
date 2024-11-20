function [idx_lower, idx_upper] = trim_data_TimeBased(input, TrimTime, input_4, tag_auralization)
%
% function [idx_lower, idx_upper] = trim_data_TimeBased(input, TrimTime, input_4, tag_auralization)
%
% this functions help to trim the input data to a specified time duration,
% which is desired for comaprison purposes.
%
% Therefore, this function identifies the max SPL over time on the input
% data and then provide lower and upper indexes that can be used to trim the data.
% The lower and upper idxs are obtained based on a (equal) time before and after
% the Max spl on the input data
%
% input - input data to be trimmed. Preferebly, an <OASPL> vector as
%         outputed from the <PANAM_SQAT_data_conversion>, because this function was
%         designed to work with this data
%
% TrimTime - time span to trim the data before and after the max spl
%
% Gil Felix Greco - Braunschweig 21.09.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global save_mat_fig

for i = 1:size(input,1)
    SPL_vs_time(i) = input{i}.OASPL_overall; % get SPL over time
    time(i) = input{i}.source_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
    %     time(i) = input{i}.retarded_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
end

dt = time(2) - time(1); % get dt from time vector. Assumes a constant dt
TrimTimebins = TrimTime/dt; % number of bins required to have <TrimTime>

[max_SPL,idx_max] = max(SPL_vs_time); % find idx of max value

idx_lower = idx_max - TrimTimebins; % idx of TrimTime before the Max. SPL
idx_upper = idx_max + TrimTimebins; % idx of TrimTime after the Max. SPL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_factor = 1000; % conversion factor (m-->km)

if idx_lower < 1 || idx_upper > length(time) % if TrimTime is larger than the input signal's time
    
    idx_lower = 1; % keep original signal
    idx_upper = length(time);
    
    fprintf( '\n*--------------------------------------------------------------------------*' );
    fprintf( '\nLog of trim function (time based) \n' );
    fprintf( '\n- Max. SPL on the data: %.4g dB SPL', max_SPL ) ;
    fprintf( '\n- Desired time before/after max. SPL: %.4g (s)', TrimTime ) ;
    fprintf( '\n- WARNING: desired trim is larger than signal''s time length.' );
    fprintf( '\n- WARNING: Original signal''s length will be maintained!!!' );
    fprintf( '\n- The auralised signal will be %.4g (s) long', (dt*length(time)) - dt ) ;
    fprintf( '\n*--------------------------------------------------------------------------*\n' );
    
    %%% PLOT SPL vs time
    
    figure('name','PROCESSING - data trimming function (time-based)');
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    xmin = time(idx_lower);
    xmax = time(idx_upper) - time(idx_lower);
    ymin = min(SPL_vs_time);
    ymax = (max_SPL - min(SPL_vs_time) ) + 5;
    
    rectangle('position',[xmin ymin xmax ymax],'FaceColor',[0.95,0.95,0.95]); hold on;
    
    plot( time, SPL_vs_time ); hold on;
    plot( time(idx_max), max_SPL,'o' );
    line(NaN,NaN,'LineWidth',10,'Color',[0.95,0.95,0.95]); % dummy plot. work around for rectangle legend
    
    xlabel('PANAM receiver time, $t_{\mathrm{P,i}}$ (s)','Interpreter','Latex');
    ylabel('SPL, $L_{\mathrm{p,Z}}$ (dB re 20 $\mu$Pa)','Interpreter','Latex');
    
    legend('PANAM input data',...
        sprintf('Max. SPL = %.4g dB', max_SPL),...
        sprintf('Total auralization time is %.4g (s)', (dt*length(time)) - dt ),...
        'Location','Best');
    
    ylim([ymin max_SPL + 5]);
    
    set(gcf,'color','w');

    %%% save plot

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_TrimData_SPL');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

    %%% PLOT flight profile (altitude vs. x distance)
    
    figure('name','PROCESSING - data trimming function (time-based)');
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    xmin = (input_4.x(idx_lower))./conv_factor;
    xmax = (input_4.x(idx_upper))./conv_factor - (input_4.x(idx_lower))./conv_factor;
    ymin = 0;
    ymax = (max( (input_4.z)./conv_factor ) ) + 5./conv_factor;
    
    rectangle('position',[xmin ymin xmax ymax],'FaceColor',[0.95,0.95,0.95]); hold on;
    
    plot( (input_4.x)./conv_factor, (input_4.z)./conv_factor ); hold on;
    line(NaN,NaN,'LineWidth',10,'Color',[0.95,0.95,0.95]); % dummy plot. work around for rectangle legend
    
    legend('PANAM input: altitude','Location','Best');
    
    xlabel('Distance, $x$~(km)','Interpreter','Latex');
    ylabel('Altitude AGL, $h_{\mathrm{AGL}}$ (km)','Interpreter','Latex');
    
    set(gcf,'color','w');
    
    %%% save plot

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_TrimData_flight_profile');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

else % trim is possible
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf( '\n*--------------------------------------------------------------------------*' );
    fprintf( '\nLog of trim function (time based) \n' );
    fprintf( '\n- Max. SPL on the data: %.4g dB SPL', max_SPL ) ;
    fprintf( '\n- Desired time before/after max. SPL: %.4g (s)', TrimTime ) ;
    fprintf( '\n- The auralised signal will be %.4g (s) long', 2*TrimTime ) ;
    fprintf( '\n*--------------------------------------------------------------------------*\n\n' );
    
    %%% PLOT SPL vs time
    
    figure('name','PROCESSING - data trimming function (time-based)');
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    xmin = time(idx_lower);
    xmax = time(idx_upper) - time(idx_lower);
    ymin = min(SPL_vs_time);
    ymax = (max_SPL - min(SPL_vs_time) ) + 5;
    
    rectangle('position',[xmin ymin xmax ymax],'FaceColor',[0.95,0.95,0.95]); hold on;
    
    plot( time, SPL_vs_time ); hold on;
    plot( time(idx_max), max_SPL,'o' );
    line(NaN,NaN,'LineWidth',10,'Color',[0.95,0.95,0.95]); % dummy plot. work around for rectangle legend
    
    xlabel('PANAM receiver time, $t_{\mathrm{P,i}}$ (s)','Interpreter','Latex');
    ylabel('SPL (dB re 20 $\mu$Pa)','Interpreter','Latex');
    
    legend('PANAM input data',...
        sprintf('Max. SPL = %.4g dB', max_SPL),...
        sprintf('Total auralization time is %.4g (s)', 2*TrimTime),...
        'Location','S');
    
    ylim([ymin max_SPL + 5]);
    
    set(gcf,'color','w');
    
    %%% save plot


    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_TrimData_SPL');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

    %%% PLOT flight profile (altitude vs. x distance)
    
    figure('name','PROCESSING - data trimming function (time-based)');
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    xmin = (input_4.x(idx_lower))./conv_factor;
    xmax = (input_4.x(idx_upper))./conv_factor - (input_4.x(idx_lower))./conv_factor;
    ymin = 0;
    ymax = (max( (input_4.z)./conv_factor ) ) + 5./conv_factor;
    
    rectangle('position',[xmin ymin xmax ymax],'FaceColor',[0.95,0.95,0.95]); hold on;
    
    plot( (input_4.x)./conv_factor, (input_4.z)./conv_factor ); hold on;
    line(NaN,NaN,'LineWidth',10,'Color',[0.95,0.95,0.95]); % dummy plot. work around for rectangle legend
    
    legend('PANAM input: altitude','Location','Best');
    
    xlabel('Distance, $x$~(km)','Interpreter','Latex');
    ylabel('Altitude AGL, $h_{\mathrm{AGL}}$ (km)','Interpreter','Latex');
    
    set(gcf,'color','w');
    
    %%% save plot

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_TrimData_flight_profile');
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

end

end


