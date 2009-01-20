function [lam_1, lam_2] = gstEigenvalues(gxy, gx2, gy2)
% function [lam_1, lam_2] = gstEigenValues(gxy, gx2, gy2)
% Create gst eigenvalues, and sort them into lam_1 and lam_2.

[rr, cc] = size(gxy);

% Preallocation for speed:
lam_1 = zeros(size(gxy));
lam_2 = zeros(size(gxy));

% Loop, calculating the eigenvalues.
for i = 1:rr
	for j = 1:cc
		eigenValues = eig([gx2(i,j), gxy(i,j); gxy(i,j), gy2(i,j)]);
		lam_1(i,j) = max(eigenValues);
		lam_2(i,j) = min(eigenValues);
	end
end
