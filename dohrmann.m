function [x,v,a,j,V] = dohrmann(data, f, B)

% BRIEF DESCRIPTION
% This function filters data using smoothing splines. The smoothing spline
% fit is then differentiated analytically, so that the outputs include the 
% filtered data with its first, second and third derivatives.
% The amount of smoothing and the particular implementation is taken from
% the following paper:
% Dohrmann CR, Busby HR, Trujillo DM (1988) Smoothing Noisy Data Using Dynamic
% Programming and Generalized Cross-Validation. J Biomech Eng 110:37-41.
%
% Steven Charles, February 2006

% INPUTS
% data: an N-by-1 vector
% f:    frequency at which data was sampled (in Hz)
% B:    smoothing parameter (this can be determined using Generalized
%       Cross-Validation with a fellow matlab script, wahba

% OUTPUTS
% x:    filtered data
% v:    first derivative of x
% a:    second derivative of x
% j:    third derivative of x
% V:    parameter involved in estimation of optimal smoothing parameter (is
%       used as input in matlab script wahba)

% OVERVIEW
% 1. Initialize parameters
% 2. Do backward sweep
% 3. Do forward sweep
% 4. Isolate filtered data and its derivatives
% 5. Compute V(B)


% 1. INITIALIZATION
% kinematics
N = length(data);
h = 1/f;
 M = [1 h (h^2)/2; 0 1 h; 0 0 1];
Theta = zeros(3,1,N);
P = [(h^3)/6; (h^2)/2; h];
% general
e = zeros(3,1,N);
e(1,1,:) = data;
U = [1 0 0; 0 0 0; 0 0 0];
I = eye(3);
% Backward sweep
D = zeros(N,1);
H = zeros(1,3,N);
R = zeros(3,3,N);
s = zeros(3,1,N);
% Forward sweep
g = zeros(N,1);
% Computing V(B)
M_hat = zeros(3,3,N);
E = zeros(3,3,N);
Q = zeros(3,3,N);   % I'm only looking at Qkj, where k=j


% 2. BACKWARD SWEEP
% Final condition
R(:,:,N) = U;
s(:,:,N) = -2*U*e(:,:,N);

% Loop from N-1 to 1
for k = N-1:-1:1
    D(k) = 1/(2*B + 2*P'*R(:,:,k+1)*P);
    H(:,:,k) = (2*R(:,:,k+1)*P)';
    R(:,:,k) = U + M'*(R(:,:,k+1) - 0.5*H(:,:,k)'*D(k)*H(:,:,k))*M;
    s(:,:,k) = -2*U*e(:,:,k) + M'*(I - H(:,:,k)'*D(k)*P')*s(:,:,k+1);
end


% 3. FORWARD SWEEP
% Initial Condition
Theta(:,:,1) = -0.5 * inv(R(:,:,1)) * s(:,:,1);
Q(:,:,1) = inv(R(:,:,1));

% Loop from 2 to N-1
for k = 1:N-1
    g(k) = -inv(2*B + 2*P'*R(:,:,k+1)*P)*P'*(s(:,:,k+1) + 2*R(:,:,k+1)*M*Theta(:,:,k));
    Theta(:,:,k+1) = M*Theta(:,:,k) + P*g(k);
    M_hat(:,:,k+1) = (M' * (I - H(:,:,k)'  * D(k) * P'))';
    E(:,:,k+1) = -P * D(k) * P';
    Q(:,:,k+1) = M_hat(:,:,k+1) * Q(:,:,k) * (M_hat(:,:,k+1)') - 2*E(:,:,k+1);
end


% 4. ISOLATE FILTERED DATA AND ITS DERIVATIVES
x = Theta(1,1,:);
x = x(:);
v = Theta(2,1,:);
v = v(:);
a = Theta(3,1,:);
a = a(:);
j = g;       % the g(k) computed above is jerk

% For comparison: unfiltered data and its derivatives
data_prime = 100*diff(data);
data_prime(end+1) = data_prime(end);
data_double_prime = 100*diff(data_prime);
data_double_prime(end+1) = data_double_prime(end);


% 5. COMPUTE MINIMIZING FUNCTION V(B)
Tr_A = sum(Q(1,1,:));
% If the variance of the noise is known, an optimal value of B can be
% chosen as the value of B which minimizes R_sigma:
%sigma = sqrt(1/12);
%R_sigma = (1/N) * dot((x-data),(x-data)) - (2*(sigma^2)/N)*(N - Tr_A) + sigma^2;
V = (1/N) * dot((x-data),(x-data)) / ((1/N) * (N - Tr_A))^2;