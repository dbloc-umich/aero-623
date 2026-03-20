function curved

N = 20;  % number of intervals for plotting

xi = linspace(0,1,N+1);
eta = linspace(0,1,N+1);

% V = coordinates of plotting points in ref space
V = zeros((N+1)*(N+2)/2, 2);
M = zeros(N+1,N+1); % mapping matrix
k = 0;
for j = 0:N
  for i = 0:(N-j)
    k = k + 1;
    V(k,1) = xi(i+1);
    V(k,2) = eta(j+1);
    M(i+1,j+1) = k;  % node number in spot i,j
  end
end

% E = subelement matrix
E = zeros(N*N,3);
k = 0;
for j = 1:N
  for i = 1:(N+1-j)
    k = k+1;
    E(k,:) = [M(i,j), M(i+1,j), M(i,j+1)];
    if (i < N+1-j)
      k = k+1;
      E(k,:) = [M(i+1,j), M(i+1,j+1), M(i,j+1)];
    end
  end
end

% add a V mapping for curved elements
Q = 2; %  geometry order

xyQ = [2,1; 4,2; 6,2; 2.5,3; 4.5,3.5; 2,5]; % 6 geometry nodes

%xyQ = [2,1; 5,2.5; 6,2; 2.5,3; 4.5,3.5; 2,5]; % 6 geometry nodes

%xyQ = [0,0; 1,0; 2,0; 0,1; 1,1; 0,2;]; % simple linear triangle
V = Ref2Glob(V, Q, xyQ);

% plot the submesh
figure(1); clf;
for k = 1:size(E,1)
  X = V(E(k,:), :);
  I = [1,2,3,1];
  plot(X(I,1), X(I,2), 'k-'); hold on;
end
plot(xyQ(:,1), xyQ(:,2), 'go', 'linewidth', 3, 'markersize', 10);


% get quadrature points in 2D, plot them in global space 
% note, you do not need them in global space for integration
[xq, yq, wq] = quad2d;
xyg = Ref2Glob([xq', yq'], Q, xyQ);
plot(xyg(:,1), xyg(:,2), 'rx', 'linewidth', 2);

%return

% integrate over the curved element: find the area
Area = 0.;
for iq = 1:length(xq)
  Jq = ElemJacobian([xq(iq), yq(iq)], Q, xyQ);
  Area = Area + det(Jq)*wq(iq);
end
fprintf(1, 'Area = %.5f\n', Area);

%return

% plot some normals on one of the edges
%edge = 2; xi_sigma = 0; eta_sigma = -1;
edge = 1; xi_sigma = -1; eta_sigma = 1;
[xq1, wq1] = quad1d;
xref = RefEdge2RefElem(edge, xq1);
xglob = Ref2Glob(xref, Q, xyQ);
plot(xglob(:,1), xglob(:,2), 'ms', 'linewidth', 2, 'markersize', 12);
edgelen = 0;
for i=1:size(xref,1)
  J = ElemJacobian(xref(i,:), Q, xyQ);
  stan = J(:,1)*xi_sigma + J(:,2)*eta_sigma; % tangent * ds/dsigma
  nvec = [stan(2), -stan(1)];
  edgelen = edgelen + norm(nvec)*wq1(i);
  plot([xglob(i,1), xglob(i,1)+nvec(1)*wq1(i)], ...
       [xglob(i,2), xglob(i,2)+nvec(2)*wq1(i)], 'm-', 'linewidth', 3);
end
edgelen
axis equal;


%--------------------------------
function xglob = Ref2Glob(xref, Q, xyQ)
xglob = xref;
for n = 1:size(xref,1)
  phi = TriLagrange(xref(n,1), xref(n,2), Q);
  xy = [0,0];
  for i = 1:size(xyQ,1)
    xy = xy + phi(i)*xyQ(i,:);
  end
  xglob(n,:) = xy;
end

%--------------------------------
function xref = RefEdge2RefElem(edge, xedge)
sigma = zeros(length(xedge), 1);
sigma(:,1) = xedge; Z = zeros(size(sigma));
if (edge == 1)
  xref = [1-sigma, sigma];
elseif (edge == 2)
  xref = [Z, 1-sigma];
elseif (edge == 3)
  xref = [sigma, Z];
else
  error('edge out of bounds');
end

%--------------------------------
function J = ElemJacobian(xref, Q, xyQ)
[phix, phiy] = GradTriLagrange(xref(1),xref(2),Q);
J = zeros(2,2);
for i = 1:size(xyQ,1)
  J = J + [xyQ(i,:)'*phix(i), xyQ(i,:)'*phiy(i)];
end


%--------------------------------
function phi = TriLagrange(x,y,order)

if (order == 2)
  phi = zeros(6,1);
  phi(1) = 1.0-3.0*x-3.0*y+2.0*x*x+4.0*x*y+2.0*y*y;
  phi(2) = 4.0*x-4.0*x*x-4.0*x*y;
  phi(3) = -x+2.0*x*x;
  phi(4) = 4.0*y-4.0*x*y-4.0*y*y;
  phi(5) = 4.0*x*y;
  phi(6) = -y+2.0*y*y;
else
  error('basis order not supported')
end


%-------------------------------
function [phix, phiy] = GradTriLagrange(x,y,order)

if (order == 2)
  phix = zeros(6,1);
  phix(1) =  -3.0+4.0*x+4.0*y;
  phix(2) =  4.0-8.0*x-4.0*y;
  phix(3) =  -1.0+4.0*x;
  phix(4) =  -4.0*y;
  phix(5) =  4.0*y;
  phix(6) =  0.0;

  phiy = zeros(6,1);
  phiy(1) =  -3.0+4.0*x+4.0*y;
  phiy(2) =  -4.0*x;
  phiy(3) =  0.0;
  phiy(4) =  4.0-4.0*x-8.0*y;
  phiy(5) =  4.0*x;
  phiy(6) =  -1.0+4.0*y;

else
  error('basis gradient order not supported')
end
  
%-------------------------------
function [xq, wq] = quad1d
n = 7; % Order 13 Gauss-Legendre points
xq = [
  0.025446043828621, 0.129234407200303, 0.297077424311301, 0.500000000000000, ...
  0.702922575688699, 0.870765592799697, 0.974553956171379
     ];
wq = [
  0.064742483084435, 0.139852695744638, 0.190915025252560, 0.208979591836735, ...
  0.190915025252560, 0.139852695744638, 0.064742483084435
     ];


%-------------------------------
function [xq, yq, wq] = quad2d
n = 16; % order 8 Dunavant points
xyq = [
  0.333333333333333, 0.333333333333333, 0.081414823414554, 0.459292588292723, ...
  0.459292588292723, 0.459292588292723, 0.459292588292723, 0.081414823414554, ...
  0.658861384496480, 0.170569307751760, 0.170569307751760, 0.170569307751760, ...
  0.170569307751760, 0.658861384496480, 0.898905543365938, 0.050547228317031, ...
  0.050547228317031, 0.050547228317031, 0.050547228317031, 0.898905543365938, ...
  0.008394777409958, 0.263112829634638, 0.263112829634638, 0.728492392955404, ...
  0.728492392955404, 0.008394777409958, 0.263112829634638, 0.008394777409958, ...
  0.728492392955404, 0.263112829634638, 0.008394777409958, 0.728492392955404
  ];
xq = xyq(1:2:end);
yq = xyq(2:2:end);  
wq = [
  0.072157803838894, 0.047545817133642, 0.047545817133642, 0.047545817133642, ...
  0.051608685267359, 0.051608685267359, 0.051608685267359, 0.016229248811599, ...
  0.016229248811599, 0.016229248811599, 0.013615157087217, 0.013615157087217, ...
  0.013615157087217, 0.013615157087217, 0.013615157087217, 0.013615157087217
     ];

