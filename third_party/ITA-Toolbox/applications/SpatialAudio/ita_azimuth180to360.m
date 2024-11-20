function azimuth_deg_360 = ita_azimuth180to360(azimuth_deg_180)
%ITA_AZIMUTH360TO180 Converts the azimuth angle range to [-180 180]
%   Detailed explanation goes here
for k=1:numel(azimuth_deg_180)
    if (abs(azimuth_deg_180(k))>180)
        error('Invalid input range');
    elseif azimuth_deg_180(k)>=0
        azimuth_deg_360(k)=azimuth_deg_180(k);
    else
        azimuth_deg_360(k)=360+azimuth_deg_180(k);
    end
end