function PLOT_eigenrays(receiver, source, eigenrays)

figure

% plot receiver
r = plot3(receiver(1), receiver(2), receiver(3), ' ok', 'MarkerFaceColor', [0 0.6 0]); hold all;

% plot limits
b = 250;
xmin = source(1,1)-b; xmax = source(end,1)+b;
ymin = -b; ymax = b;
zmin = 0; zmax = source(end,3)+b;
set( gca,'XLim',[xmin xmax],'YLim',[ymin ymax],'ZLim',[zmin zmax] );
% set( gca,'XLim',[receiver(1)-5 receiver(1)+5],'YLim',[receiver(2)-5 receiver(2)+5],'ZLim',[0 receiver(3)+5] ); % close view on the receiver

% plot trajectory
t = plot3(source(:,1), source(:,2), source(:,3), 'k-');

view(0,0);
xlabel('Distance, $x$ (m)','Interpreter','Latex'); 
ylabel('Lateral Distance, $y$ (m)','Interpreter','Latex'); 
zlabel('Altitude AGL, $h_\mathrm{AGL}$ (m)','Interpreter','Latex'); 
set(gcf,'color','w');
    
 for i = 1:size(source,1)
     
    head = plot3( source(i,1), source(i,2), source(i,3), ' ok', 'MarkerFaceColor', [0.6 0 0] ); % plot receiver
    rays = plot3( eigenrays(i,1).r.x, eigenrays(i,1).r.y, eigenrays(i,1).r.z, 'b-' ); % plot (0 order) rays
    rays_1 = plot3( eigenrays(i,2).r.x, eigenrays(i,2).r.y, eigenrays(i,2).r.z, 'b--' ); % plot (1st order) rays
      
    legend([r, head, t, rays, rays_1], 'Receiver', 'Source', 'Trajectory','Direct path', '1st  order reflection','Location','NW','Interpreter','Latex')
    
    drawnow;
    pause(0.01);
    delete(head); % delete source position of current frame
    delete(rays); % delete direct path from current frame
    delete(rays_1); % delete 1st order reflection from current frame    
    
end

end