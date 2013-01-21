function chi2 = fdbk_fun(X,phi1,phi2,r);

%   chi = fdbk_fun(X,phi1,phi2,r);
%
%   Function to be minimized by fminsearch to find feedback phase for
%   LiTrack 13-card.
%
%   INPUTS:     X:      phase [rad]
%               phi1:   phase of 1st section [rad]
%               phi2:   phase of 2nd section [rad]
%               r:      ratio of (E-E0)/eV0 [ ]
%   OUTPUTS:    chi2:   gets minimized to zero
%======================================================================

chi2 = ((cos(phi1) + cos(phi2))*cos(X) + (sin(phi2) - sin(phi1))*sin(X) - r)^2;