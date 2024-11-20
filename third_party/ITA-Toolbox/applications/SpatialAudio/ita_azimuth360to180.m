function azimuth_deg_180 = ita_azimuth360to180(azimuth_deg_360)
%ITA_AZIMUTH360TO180 Converts the azimuth angle range to [-180 180]
%   Detailed explanation goes here
for k=1:numel(azimuth_deg_360)
    if (azimuth_deg_360(k)<0)
        error('Negative input values not allowed');
        sign=-1;
    else
        sign=1;
    end
    azimuth_deg_360(k)=rem(azimuth_deg_360(k),360);
    if azimuth_deg_360(k)<180
        azimuth_deg_180(k)=sign*azimuth_deg_360(k);
    else
        azimuth_deg_180(k) = sign*rem(azimuth_deg_360(k),180)-180;
    end
end