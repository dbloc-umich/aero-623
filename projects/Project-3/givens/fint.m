function fint
% This code integrates a function over a triangular domain using
% quadrature.

% load mesh
V = load('V.txt');
E = load('E.txt');
figure(1); clf; plotmesh(V,E); axis equal; axis off;

% obtain quadrature rule in reference triangle
[xiq, etaq, wq] = getquad;

% loop over elements (rows of E)
I = 0; % this will be the desired integral
for k = 1:size(E,1),
  N = E(k,:); % 3 nodes for element k
  x = V(N,1); % x-locations of these 3 nodes
  y = V(N,2); % y-locations of these 3 nodes
  J = [x(2)-x(1), x(3)-x(1);
       y(2)-y(1), y(3)-y(1)]; % Jacobian matrix
  detJ = det(J);
  
  % map quadrature points to global space
  xq = zeros(length(wq),1); yq = xq;
  for q = 1:length(wq),
    xq(q) = x(1) + (x(2)-x(1))*xiq(q) + (x(3)-x(1))*etaq(q);
    yq(q) = y(1) + (y(2)-y(1))*xiq(q) + (y(3)-y(1))*etaq(q);
  end
  
  % evaluate f at the quadrature points (global space)
  fq = getf(xq, yq);
  I = I + wq'*fq*detJ  % add to total integral
  
  plot(xq, yq, 'rx', 'linewidth', 2); pause;
end
I

%-----------------------------------
function plotmesh(V,E)
for k = 1:size(E,1);
  N = E(k,:); % 3 nodes for element k
  x = V(N,1); % x-locations of these 3 nodes
  y = V(N,2); % y-locations of these 3 nodes
  x = [x; x(1)]; y = [y; y(1)];
  plot(x,y, 'k-'); hold on;
end

%-----------------------------------
% function that we want to integrate
function f = getf(x,y)
%f = 5.0*ones(size(x));
f = x;


%-----------------------------------
% quadrature rule for a triangle
function [xq, yq, wq] = getquad
% Order 3 Dunavant Points  (4 points total)
xq = [1/3, 0.6, 0.2, 0.2]';
yq = [1/3, 0.2, 0.2, 0.6]';
wq = [-0.28125, 0.260416666666667, 0.260416666666667, 0.260416666666667]';

