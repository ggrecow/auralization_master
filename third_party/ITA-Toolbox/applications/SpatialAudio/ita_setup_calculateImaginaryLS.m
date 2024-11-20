function [LS_imag] = ita_setup_calculateImaginaryLS(LS_pos, contours)
%ITA_SETUP_CALCULATIMAGINARYLS Fill holes in LS setups with imaginary
%loudspeakers.

%% Sort contour points & calculate imaginary speaker positions
% X = @(phi) [[1 0 0];[0 cos(phi) -sin(phi)];[0 sin(phi) cos(phi)]];
Y = @(phi) [[cos(phi) 0 -sin(phi)];[0 1 0];[sin(phi) 0 cos(phi)]];
Z = @(phi) [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]];
if ~isempty(contours)
    for idx = 1:numel(contours)
        % get contour points
        contour = LS_pos.n(contours{idx});
        % mean vector r_mean
        r_mean = itaCoordinates([mean(contour.x), mean(contour.y), mean(contour.z)]);
        contour_mod = contour;
        
%         % scale point vectors so that they lay on a plane perpendicular to r_mean
%         t = r_mean.r^2 ./ contour.dot(repmat(r_mean, contour.nPoints));
%         contour_mod.r = t .* contour.r;
%         contour_mod.r(t<0, :) = -contour_mod.r(t<0, :);
        
%         % project points to plane by scaled addition of r_mean
%         t = 1 - contour.dot(repmat(r_mean, contour.nPoints)) / r_mean.r^2;
%         r_project = repmat(r_mean, contour.nPoints);
%         r_project.r = t .* r_project.r;
%         contour_mod = contour + r_project;

        % Rotate projected plane
        R = Z(r_mean.phi) * Y(r_mean.theta) * Z(-r_mean.phi);
        contour_mod = itaCoordinates(contour_mod.cart*R');
%         contour_mod.phi(t<0) = contour_mod.phi(t<0) + pi;
        
        [~, jdx] = sort(contour_mod.phi);
        contours{idx} = contour.n([jdx; jdx(1)]);
    end
    
    % calculate imaginary speaker positions
    LS_imag = itaCoordinates(zeros(numel(contours), 3));
    for idx = 1:numel(contours)
        for jdx = 1:contours{idx}.nPoints-2
            v1 = contours{idx}.n(jdx+1) - contours{idx}.n(jdx);
            v2 = contours{idx}.n(jdx+2) - contours{idx}.n(jdx+1);
            LS_imag.cart(idx, :) = LS_imag.cart(idx, :) + v1.cross(v2).cart;
        end
        LS_imag.r(idx) = LS_pos.r(1);
    end
else
    LS_imag = itaCoordinates();
end

end
