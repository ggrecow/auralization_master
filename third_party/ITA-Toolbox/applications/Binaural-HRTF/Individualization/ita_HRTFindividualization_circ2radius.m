function dS = ita_HRTFindividualization_circ2radius(wS,cS)
% function to calculate the (half) head depth [mm] with given (half) head width [mm]
% and (half) head circumference [mm]
%
% INPUT:
% wS ... (half) head width of subject, measured fromt tragus to tragus with caliper [mm] -> wS = (head_width - 12)/2
% cS ... (half) head circumference, measured from tragus to tragus over
%               bend of nose (shortest distance) [mm]
%
% OUTPUT:
% dS ... (half) head depth, measured horizontally from bend of nose to ear canal entrance [mm]
%
% Author: Ramona Bomhardt, rbo@akustik.rwth-aachen.de

U2 = 2*cS;
sigma10 = 0.177*sqrt(wS.*(3*U2-12.567*wS));
sigma6 = 0.2122*U2-0.7778*wS+sigma10;
dS = sigma6;


%(bis wS = a / 2 recht genau, bei sehr flachen Ellipsen ungenauer)
end