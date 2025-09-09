function PLOT_eigenrays(receiver, source, eigenrays)

figure

% plot receiver
r = plot3(receiver(1), receiver(2), receiver(3), 'o', 'MarkerEdgeColor', [0, 0 ,0], 'MarkerFaceColor', [0.39,0.83,0.07], 'MarkerSize', 8 ); hold all;

% --- Stretch axes in X and Y ---
stretchY = 0.9;
stretchX = 1.5;
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3)*stretchX, pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)*stretchX, pos(4)*stretchY])
movegui(h, 'center');


% plot limits
b = 250;
% b = 1;
xmin = source(1,1)-b; xmax = source(end,1)+b;
ymin = -b; ymax = b;
zmin = 0; zmax = source(end,3)+b;
set( gca,'XLim',[xmin xmax],'YLim',[ymin ymax],'ZLim',[zmin zmax] );
% set( gca,'XLim',[receiver(1)-5 receiver(1)+5],'YLim',[receiver(2)-5 receiver(2)+5],'ZLim',[0 receiver(3)+5] ); % close view on the receiver

% plot trajectory
t = plot3(source(:,1), source(:,2), source(:,3), 'k-');

view(0,0);
xlabel('$x$ (m)','Interpreter','Latex'); 
ylabel('Lateral Distance, $y$ (m)','Interpreter','Latex'); 
zlabel('Altitude AGL, $h_\mathrm{AGL}$ (m)','Interpreter','Latex'); 
set(gcf,'color','w');

% --- cosmetics ---
FontSize = 18; 
ax = gca;
ax.FontSize = FontSize;   % all font sizes
ax.LineWidth = 1;   % box line width
% ax.XTick = [];      % remove x ticks
% ax.XTickLabel = []; % remove x label
ax.XTick = linspace(ax.XLim(1), ax.XLim(2), 5);

%%% define gif filename
save_gif = 0; % boolean flag
gifname = 'eigenrays.gif';

for i = 1:size(source,1)

    head = plot3( source(i,1), source(i,2), source(i,3), 'o', 'MarkerEdgeColor', [0, 0 ,0], 'MarkerFaceColor', [0.30,0.75,0.93], 'MarkerSize', 8 ); % plot aircraft/source
    rays = plot3( eigenrays(i,1).r.x, eigenrays(i,1).r.y, eigenrays(i,1).r.z, 'b-' ); % plot (0 order) rays
    % rays_1 = plot3( eigenrays(i,2).r.x, eigenrays(i,2).r.y, eigenrays(i,2).r.z, 'b--' ); % plot (1st order) rays

    % legend([r, head, t, rays, rays_1], 'Receiver', 'Source', 'Trajectory','Direct path', '1st  order reflection','Location','NW','Interpreter','Latex')
    legend([r, head, t, rays], 'Receiver', 'Aircraft', 'Trajectory','Direct sound path', 'Location','northoutside', 'numcolumns', 4, 'Interpreter','Latex')
    set(findobj(gcf,'Type','Legend'),'FontSize',FontSize);
    drawnow;

    if save_gif == 1
        %%% capture frame and write to GIF
        frame = getframe(gcf);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256,'dither');
        if i == 1
            imwrite(A,map,gifname,"gif","LoopCount",Inf,"DelayTime",0.05);
        else
            imwrite(A,map,gifname,"gif","WriteMode","append","DelayTime",0.05);
        end
    else
    end

    pause(0.01);
    delete(head); % delete source position of current frame
    delete(rays); % delete direct path from current frame
    % delete(rays_1); % delete 1st order reflection from current frame

end

end