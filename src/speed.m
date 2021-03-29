% SPEED estimates the speed of a propagating front from arrival times (x,t)
%       by centered finite difference
%
% v = speed(x,t)
%
% INPUTS	x(:)	positions
%		t(:)	arrival times at positions x
% 
% OUTPUTS	v(:)	front speed
%
function v = speed(x,t)

v = zeros(size(t));
v(1) = (x(2)-x(1))./(t(2)-t(1));
v(2:end-1) = (x(3:end)-x(1:end-2))./(t(3:end)-t(1:end-2));
v(end) = (x(end)-x(end-1))./(t(end)-t(end-1));
