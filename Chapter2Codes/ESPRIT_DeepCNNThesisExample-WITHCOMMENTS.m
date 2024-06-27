clear all; close all; clc;
% setup file for Matlab figures
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultTextFontSize',18)

%% Problem Creation and Eigenvalues/Eigenvectors created from Correlation Matrix
randseed=11; rand('seed',randseed);randn('state',0)
k=2*pi;
cs=2;

N=5; d=0.5;   % number of array elements = N and spacing d 1/2 wavelength
w=ones(1,N); % element weights originally uniform all set to 1 
x=[-N/2*d+d/2:d:N/2*d-d/2];

Ns=2;  % number of signals
sa=[1 1]; ts=[-50 30];
tsr=ts*pi/180;
K=1000; % number of signal samples
s=sign(randn(Ns,K)).*(sa'*ones(1,K));

% Gaussian noise
sig0=.01;
nois=sqrt(sig0)*randn(N,K);

A=exp(1j*k*x'*sin(tsr)); % --> Signal Array Manifold Vector
% In this example eigenbeams no interferer signals 

Rss=s*s'/K;      % signal correlation matrix
Rnn=(nois*nois')/K;         %noise correlation matrix
Rns=(nois*s')/K;           % noise/signal correlation matrix
Rsn=(s*nois')/K;           % signal/noise correlation matrix
Rrr=1/K*A*Rss*A'+A*Rsn+Rns*A'+Rnn;        % array correlation matrix
%[eigvec, eigval] = eig(Rrr);

% Try using overall signal instead of addition of individual components 
total_signal = A*s+nois; % X(t)
Rx = 1/K*(total_signal*total_signal');
[eigvec,eigval] = eig(Rx); % same result as Rrr

%% ESPRIT Algorithm from Van Trees
% Compare LS, TLS
% LS Method: 

% sort eigenvalues 
[~,indices] = sort(diag(eigval), 'descend');
% D_sorted = D(indices, indices);
eigvec_sorted = eigvec(:,indices);
% extract first K eigenvectors (columns):
% K = Ns the number of signal sources 2 from the 5 eigenvectors selected
output_eigvec = eigvec_sorted(:, 1:Ns);

% create the subarray based on those defined with the corresponding ASV
%subarray definition 1: 1 2 3 4 --> 2: 2 3 4 5 
subarr1 = output_eigvec(1:N-1,:); % same as V = ASV
subarr2 = output_eigvec(2:end,:); % ASV for subarray2

% Use Least Squares Method to find psi values (LS):
% N-1xN-1 matrix take inverse then multiply by 
% [181x181] x [181x4] = [181x4] x [4x181] = [181 x 181]
psi = inv(subarr1'*subarr1)*subarr1'*subarr2;

% eigen decomposition of psi:
[eigvec_psi,eigval_psi] = eig(psi);

% extract using psi = k*d*cos(theta)
% use angle for tan-1 of resulting eig_val_psi
% diag because other values in matrix should be near 0
doa_positions = angle(diag(eigval_psi));
doa_theta_LS = asind(doa_positions/(k*d))

%% For plot: just output with ESPRIT not entire range
% this is because ESPRIT is not an exhaustive search - just uses invariance of the array itself
fig_ESPRIT = figure(1);

legend_text = [];
for i = 1:Ns
    xline(doa_theta_LS(i));
    %plot(doa_theta_LS(i),1,'o','MarkerSize',14);
    hold on;
    legend_text = [legend_text "ESPRIT LS DOA Estimate for Source "+string(i)+": "+string(doa_theta_LS(i))];
end

title("ESPRIT DOA: LS & TLS")
xlim([-75 75])
ylim([0 1.25])
xlabel("\theta \circ")
legend(legend_text,FontSize=12,Location='best');

%% TLS ESPRIT
C = [subarr1';subarr2']*[subarr1 subarr2];

[U,~,~] = svd(C);
[eigvec_subarr,eigval_subarr] = eig(C);

% D = number of signals in psi-space so Ns
% U = eigenvectors associated with decomposition of signal space x = A*s + nois
V12 = U(1:Ns,Ns+1:2*Ns);         % D x D = Ns x Ns
V22 = U(Ns+1:2*Ns,Ns+1:2*Ns);     % D x D = Ns x Ns

psi = -V12/V22;               % Eq. (9.122) in [1]
[psieigvec,psieigval] = eig(psi); % Eq. (9.123) in [1]
%   Extract angle information estimated from two subarrays based on the
%   distance of the phase center between subarrays.
doas_TLS = angle(diag(psieigval));
doas_TLS = asind(doas_TLS/(k*d))

for i=1:Ns
    plot(doas_TLS(i),1,'squre','MarkerSize',20); hold on;
    legend_text = [legend_text "ESPRIT TLS DOA Estimate for Source "+string(i)+": "+string(doas_TLS(i))];
end
legend(legend_text,FontSize=12,Location='best');
setTightMargins(0.05)
saveas(fig_ESPRIT,'ESPRIT_LS_TLS.jpg');

%% Moved Comments:

% https://github.com/shantistewart/MUSIC-and-ESPRIT-Algorithms/blob/master/MATLAB/ESPRIT.m
% for theirs, X = signal, K = # sources, T = snapshots taken, N =  array size (determined from X)
% NOT SURE WHY YOU HAVE TO SORT EIGEN VECTORS 
% use eigen values and vectors and sort them to use for subarray
% DON'T USE THE d value for the output sorted output eigenvectors, do by hand line below***** 
% Try with EigenDecomp EigenVectors:
%V12 = eigvec_subarr(1:Ns,Ns+1:2*Ns);         % D x D = Ns x Ns
%V22 = eigvec_subarr(Ns+1:2*Ns,Ns+1:2*Ns);     % D x D = Ns x Ns
% No Unitary Used 

%% MATLAB: espritdoa script using TLS from Van Trees OAP
%{ 
% Row weighting
Ns = M-saSpacing; %number of elements in a subarray
ms = rweight;
w = min(ms,Ns-ms+1);                             % Eq 9.133 in [1]
weights = diag(sqrt([1:w-1 w*ones(1,Ns-2*(w-1)) w-1:-1:1])); % Eq 9.132 in [1]
O = zeros(Ns,saSpacing);

% Selection Matrices
Js1 = [weights O]; % Eq 9.134 in [1]
Js2 = [O weights]; % Eq 9.135 in [1]

% Selecting subarray signal subspaces
Us1 = Js1*eigenvects(:,1:D);
Us2 = Js2*eigenvects(:,1:D);
% TLS-ESPRIT
C = [Us1';Us2']*[Us1 Us2];    % Eq. (9.123) in [1]
[U,~,~] = svd(C);             % C is 2*D x 2*D
V12 = U(1:D,D+1:2*D);         % D x D
V22 = U(D+1:2*D,D+1:2*D);     % D x D
psi = -V12/V22;               % Eq. (9.122) in [1]
psieig = eig(psi);
%   Extract angle information estimated from two subarrays based on the
%   distance of the phase center between subarrays.
doas = 1/saSpacing*angle(psieig);
%}

%% Unitary ESPRIT

%% ESPRIT with changes in antenna positions


%% Github Try to Used:
% https://github.com/shantistewart/MUSIC-and-ESPRIT-Algorithms/blob/master/MATLAB/ESPRIT.m

%{

K = 2;
% [5x2] x [5x1000] = [5x1000] + [5x1000] = [5x1000] = [NxK] K = T for them = # of snapshots taken
X = A*s+nois;

thetas = ESPRIT(X,K)
% Function Description: implements ESPRIT algorithm for direction of
% arrival estimation.
function thetas = ESPRIT(X, K)
% wavelength:
lambda = 2;
% sensor separation:
dist = 1;

% get dimensions of X:
dim = size(X);
N = dim(1);
T = dim(2);

% estimate autocorrelation matrix:
Rx = (1/T)*(X*X');

% perform eigendecomposition of autocorrelation matrix:
[V, D] = eig(Rx);
% sort eigenvectors in decreasing order of their eigenvalues:
[~,indices] = sort(diag(D), 'descend');
% D_sorted = D(indices, indices);
V_sorted = V(:,indices);
% extract first K eigenvectors (columns):
V1 = V_sorted(:, 1:K);

% W1 and W2 matrices:
W1 = V1(1:N-1, :);
W2 = V1(2:end, :);

% solve optomization by LS:
psi = inv(W1'*W1) * W1' * W2;

% perform eigendecomposition of psi:
[V_psi, D_psi] = eig(psi);

% extract theta_k's:
angles = angle(diag(D_psi));
sines = angles / (-2*pi*(dist/lambda));
thetas = asin(sines);
% convert to degrees:
thetas = (180/pi) * thetas;

end

%}

%% Notes about ESPRIT from "Changing Array Positions" Script:

% Estimation of Signal Parameters via Rotational Invariance Techniques 
% subarray from original array separated by distance d 
% typically with even number array elements
% ------------------------------------------------------------ %
% Least Squares (LS): Toeplitz Approximation Method (TAM) - page 1200 pdf
% 2-non-overlapping subarrays - example 10 element ULA
% where Ψ= kz*d= (2*π/λ)*cosΘ*d = 2π/λ*uz*d and D is the number of signals
% Ns is the # elements in the subarray 
% --> Description of Required Matrices <--
% Subdivided into 1,3,5,7,9 and 2,4,6,8,10 arrays Js1 and Js2
% Js1 = 5x9 matrix with 1's in positions [1,1],[2,3],[3,5],[4,7],[5,9]
% Js2 = 5x9 matrix with 1's in positions [1,2],[2,4],[3,6],[4,8],[5,10]
% Create Subarray manifold vectors from Array Manifold/Steering Vector (ASV): V
% V1 = Js1*V and V2 = Js2*V 
% Exploit shift invariance property of array: V2 = V1*Φ 
% where Φ = diag [exp(j*ds*Ψ1),exp(j*ds*Ψ2),...,exp(j*ds*ΨD)] 

% --> Description of Obtaining Θ from LS ESPRIT Method 

% ------------------------------------------------------------ %
% Total Least Squares (TLS): 

% ------------------------------------------------------------ %
% Unitary ESPRIT: 