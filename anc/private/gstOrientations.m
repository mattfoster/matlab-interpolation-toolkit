function phi_2 = gstOrientations(lam_1, lam_2, gxy, gx2, gy2) 
% function [phi_1, phi_2] = gstOrientations(lam_1, lam_2, gxy, gx2, gy2) 
% 	Calculate gst orientations from eigenvalues and orientation prod.
%----------------------------------------------------
% 
% phi_1 = atan((lam_1 - gx2) ./ gxy);
% phi_1(gxy == 0) = 0;
% 
phi_2 = atan(gxy ./ (lam_2 - gy2));
phi_2(isnan(phi_2)) = 0;
