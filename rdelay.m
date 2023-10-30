function [alpha,beta,nbisec,nQZ,Lmax,Lopt] = rdelay(A,B,tau,tol_eig)
% RDELAY computes the distance to the imaginary axis
% from the delayed system
%   dx/dt(t) = Ax(t) + Bx(t-tau)
% where tau > 0 is the delay, by bisection.
% The output interval [alpha beta] contains the distance.
% nbisec is the number of bisections; nQZ is the number of QZ calls.
% Algorithm uses the subintervals [2L-1,2L+1], where -Lmax <= L <= Lmax
% Lopt is the interval of the answer [alpha,beta]
%
% tol_eig is the admitted maximum deviation of "the imaginary eigenvalues"
% from the imaginary axis: e.g. tol_eig = 1e-10
%
% Choose a suitable relgap > 1 inside the function: relgap = 1.0001
% if closer bounds are preferred, or relgap = 2 if a larger gap between
% alpha and beta is satisfactory.

% Alexander N. Malyshev, October 30, 2023

n = size(A,1);
n2 = n*2;

rs = min(svd(A+B));
mu = rs + norm(A) + norm(B);
theta = tau*mu;
Lmax = ceil((theta-1)/2);
if isreal(A) && isreal(B)
  Lmin = 0;
else
  Lmin = -Lmax;
end

tol0 = tol_eig*mu; % threshold for the zero

% Pade denominator: q = q_1 + q_2*z + ... + q_{N+1}*z^N
N = 6;
q = zeros(1,N+1);
q(1) = 1;
for k = 1:N
  if theta >= 1
    q(k+1) = (N-k+1)/(2*N-k+1)/k * q(k); % q ~ exp(-z)
  else
    q(k+1) = (N-k+1)/(2*N-k+1)/k * q(k)*theta; % q ~ exp(-theta*z)
  end
end

BB = eye(n2*(N+1));
AA = zeros(n2*(N+1));
AA((n2+1):(n2*(N+1)+1):end) = -1;

nbisec = 0;
nQZ = 0;
relgap = 1.0001; % choose relgap > 1, e.g. 2 or 1.0001
alpha = 0;
beta = rs;
Lopt = 0;

while beta > relgap*max(alpha,tol0)
  nbisec = nbisec+1;
  s = sqrt(beta*max(alpha,tol0));
  for L = Lmin:Lmax
    % form the block Frobenius pencil AA+lambda*BB
    if theta >= 1
      C = eye(n)*(-q(N+1)/tau);
    else
      C = eye(n)*(-q(N+1)*mu);
    end
    BB(1:n2,1:n2) = [zeros(n) C; C*(-1)^(N+1) zeros(n)]; % P_{N+2};
    if theta >= 1
      C = (A-eye(n)*(2i*L/tau))*q(1) + B*(exp(-2i*L)*q(1));
    else
      C = (A+B)*q(1);
    end
    AA(1:n2,1+n2*N:n2*(N+1)) = ...
      [-s*q(1)*eye(n) C; C' -s*q(1)*eye(n)]; % P_{1}
    for k = 2:N+1
      if theta >= 1
        C = (A-eye(n)*(2i*L/tau))*q(k) + B*(exp(-2i*L)*(-1)^(k-1)*q(k)) ...
          - eye(n)*(q(k-1)/tau);
      else
        C = A*q(k) + B*((-1)^(k-1)*q(k)) - eye(n)*(q(k-1)*mu);
      end
      AA(1:n2,1+n2*(N+1-k):n2*(N+2-k)) = ... % P_{k}
        [-s*q(k)*eye(n) C; C'*(-1)^(k-1) -s*(-1)^(k-1)*q(k)*eye(n)];
    end
    % scale the 1st block row of the pencil AA+lambda*BB
    rho = norm([BB(1:n2,1:n2) AA(1:n2,:)]);
    BB(1:n2,1:n2) = BB(1:n2,1:n2)/rho;
    AA(1:n2,:) = AA(1:n2,:)/rho;
    % is there an eigenvalue of AA+lambda*BB near 1i*[-1,1]?
    nQZ = nQZ + 1;
    E = eig(AA,-BB);
    imag_eig = any(abs(real(E)) < tol_eig & abs(imag(E)) <= 1+tol_eig);
    if imag_eig, break, end
  end
  % update
  if imag_eig % there is an imaginary eigenvalue
    beta = s;
    Lopt = L;
  else
    alpha = s;
  end
end