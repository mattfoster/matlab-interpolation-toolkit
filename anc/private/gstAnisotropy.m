function [aniso, energy] = gstAnisotropy(lam_1, lam_2)
% function [aniso, energy] = gstAnisotropy(lam_1, lam_2)
%	Compute anisotropy and energy from given eigenvaulues.

energy = lam_1 + lam_2;
diff = lam_1 - lam_2;
aniso = diff ./ energy;

aniso(energy == 0) = 0; % necessary?
aniso(diff == 0) = 0;
