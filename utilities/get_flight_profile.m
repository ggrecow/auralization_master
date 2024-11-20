function OUT = get_flight_profile(PATH, show, proc, save_figs, tag, figuresdir)
% function get_flight_profile(PATH, show, proc, save_figs, tag, figuresdir)
%
% plot a flight profile from PANAM
%
% INPUTS:
%       * PATH : string - with the path to data (usually called: geschw_hoehe_verlauf.dat')
%       * show : boolean - plot flight profile  0 (no); 1(yes),
%       * proc : boolean - 0 (approach); 1 (departure); 2 (flyover) - only changes the x-axis label
%       * save_figs : boolean - 0 (no); 1(yes),
%       * tag : string - name tag to save figure
%       * figuresdir : string - dir path to save figs
%
% OUTPUTS:
%       OUT : struct
%       OUT.x - position in x-axis, in meters
%       OUT.y - position in y-axis, in meters
%       OUT.z - altitude in meters
%       OUT.TAS -  true airspeed in meters/sec
%       OUT.thrust - thrust in kN
%
% Author: Gil Felix Greco, 09.06.2023 (updated 07.12.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global save_mat_fig

%% load data

[HEADER,DATA] = mhdrload(PATH);

% for plotting
conv_factor = 1000; % conversion factor (m-->km)
x = DATA(:,1)./conv_factor;
altitude = DATA(:,3)./conv_factor;
TAS = DATA(:,4);
thrust = DATA(:,5);

% outputs
OUT.x = DATA(:,1); % position in x-axis, in meters
OUT.y = DATA(:,2); % position in y-axis, in meters
OUT.z = DATA(:,3); % altitude in meters
OUT.TAS = DATA(:,4); % true airspeed in meters/sec
OUT.thrust =  DATA(:,5); % thrust in kN

%% plot

if show == true
    
    figure('name','INPUT - flight profile from PANAM');
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(h,'defaultAxesColorOrder',[left_color; right_color]);
    hold all;
    
    % altitude
    yyaxis left
    a=plot(x,altitude);
    ylabel('Altitude AGL, $h_{\mathrm{AGL}}$ (km)','Interpreter','Latex');
    
    if proc == 0 % (approach)
        xlabel('Distance from runway threshold, $x_{\mathrm{app}}$~(km)','Interpreter','Latex');
    elseif proc == 1 % (departure)
        xlabel('Distance from brake release, $x_{\mathrm{dep}}$~(km)','Interpreter','Latex');
    elseif proc == 2 % (flyover)
        %         xlabel('Distance from receiver, $x_{\mathrm{rec}}$~(km)','Interpreter','Latex');
        xlabel('Distance, $x$~(km)','Interpreter','Latex');
    end
    
    % ylim([0 2500])
    
    % TAS
    yyaxis right
    b=plot(x,TAS,'b-');
    
    % thrust
    yyaxis right
    c=plot(x,thrust,'r-');
    ylabel('True airspeed (m/s), Thrust (kN)','Interpreter','latex');
    
    % xlim([-18026 -3026])
    % ylim([0 250])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gcf,'color','w');
    
    legend([a b c],{'Altitude','True airspeed','Thrust'},'Location','Best');
    
    legend boxoff
    box on
    
    % ax = gca;
    
    % ax.XTick = [-3026-15000 -3026-10000 -3026-5000 -3026];
    % ax.XTickLabels = {'-15','-10','-5','0'};
    %
    % ax.XAxis.MinorTick = 'on';
    % ax.XAxis.MinorTickValues =  -3026-15000:500:-3026;
    %
    % ax.YAxis(2).MinorTick = 'on';
    % ax.YAxis(2).MinorTickValues =  0:10:250;
    %
    % ax.YAxis(1).MinorTick = 'on';
    % ax.YAxis(1).MinorTickValues =  0:100:2500;
    %
    % x0=0;
    % y0=0;
    % width=20;
    % height=10;
    % set(gcf,'position',[x0,y0,width,height])
    

    if save_figs==1

        filename = strcat(figuresdir, 'flight_profile_', tag);
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );

    else
    end
    
else
end

end