function angle_out = ESPRIT_DOA_Estimate(signal_in,varargin)
% ------------------------------------------ % 
% Obtain the output ESPRIT DOA from inputs:
% INPUTS:
%
% signal_in: Signal matrix with K snapshots signal_in (signal_in = s for BPSK K snapshots of [1,-1])
% rand_seed: Random Seed to use: Default = 11 (for repeatability) [varargin]
% theta_in: input signal(s) array in degrees [used for ASV]
% element_pos: Element Position Array: Default 0.5λ spacing 5 element 
% SNR_in: Input SNR in dB ex: 20dB -> sig0 = 0.01 (noise power): Default 20dB
% SNR EX2: dB = 10*log10(x) x = 0.01 -> dB = -20 for noise
% 
% OUTPUTS: 
% 
% angle_out: output array of angles calculated by ESPRIT
% Note: No Search Space or Spectrum Output by ESPRIT compared to MUSIC
% 
% Note: varargin is a cell array - but requires consistent input order 
% Note: ASV = Array Steering Vector (exponential term describing phase
% shift between elements in Uniform Linear Array (ULA))
% ------------------------------------------ % 
%% Code for parsing and carrying out ESPRIT 

% Parse Number of input arguments 
% error if none or > 7 arguments provided

if nargin == 0 || nargin > 5
    error("ERROR: Incorrect number of parameters input, view the comments in the function");
end
% default rand_seed & theta_in & element_pos & SNR
if nargin == 1
    rand_seed = 11;
    theta_in = 30; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default theta_in & element_pos & SNR
if nargin == 2
    rand_seed = varargin{2};
    theta_in = 30; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default element_pos & SNR 
if nargin == 3
    rand_seed = varargin{2};
    theta_in = varargin{3}; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default SNR 
if nargin == 4
    rand_seed = varargin{2};
    theta_in = varargin{3};
    element_pos = varargin{4};
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% all defined for nargin == 6 - also output ASV for user 
if nargin == 5
    rand_seed = varargin{1};
    theta_in = varargin{2};
    element_pos = varargin{3};
    sig0 = 10^(-varargin{4}/10);
end

% define error if theta provided outside of range -90° to 90°
if ~isempty(find(abs(theta_in) > 90))
    error("ERROR: theta(Θ) provided is outside range [-90°,90°]"); 
end

% define theta in terms of radians and ASV 
theta_in_rad = theta_in*pi/180;
% define expected ULA element positions - create ASV with this
N = 5; d = 0.5;
element_pos_expected = [-N/2*d+d/2:d:N/2*d-d/2];

% create expected Array Steering Vector (ASV) 
%A=exp(1j*(2*pi)*element_pos_expected'*sin(theta_in_rad));
% for ESPRIT just use perturbed ASV
A=exp(1j*(2*pi)*element_pos'*sin(theta_in_rad));

% ------------------------------------------ % 

% grab number of snapshots, K for signal_in
K = size(signal_in,2);
% create overall signal
nois=sqrt(sig0)*randn(N,K);
% Create covariance of signal and noise (near 0-valued)
% first signal * hermitian second signal (' notation) 
Rsn = signal_in*nois'; 
Rns = nois*signal_in';
% Create covariance of noise and signal then multiply all but noise by ASV
Rss = signal_in*signal_in'; % cov(signal)
Rnn = nois*nois'; % cov(nois)
% Multiply signal by expected ASV (A) - Normalize by # snapshots: K
% create the cov matrix (Rrr) for input signal [x*x^H / K]
% Don't create Cov from ASV because ASV not applicable unless theta known
%signal A*A' | sig/nois A | nois/sig A' | all normalized by 1/K
cov_matrix = 1/K*((A*Rss*A' + A*Rsn + Rns*A') + Rnn); 

% determine the eigenvalue decomposition
[eigvec, eigval] = eig(cov_matrix);
num_signals = length(theta_in); 

% ------------------------------------------ %
% Compare LS, TLS, and Unitary at minimum
% LS Method: 
% https://github.com/shantistewart/MUSIC-and-ESPRIT-Algorithms/blob/master/MATLAB/ESPRIT.m
% for theirs, X = signal, K = # sources, T = snapshots taken, N =  array size (determined from X)
% ------------------------------------------ %

% --- ESPRIT Algorithm --- %

% NOT SURE WHY YOU HAVE TO SORT EIGEN VECTORS 
% use eigen values and vectors and sort them to use for subarray
% DON'T USE THE d value for the output sorted output eigenvectors, do by hand line below***** 
[~,indices] = sort(diag(eigval), 'descend');
% D_sorted = D(indices, indices);
eigvec_sorted = eigvec(:,indices);
% extract first K eigenvectors (columns):
% K = num_signals the number of signal sources 2 from the 5 eigenvectors selected
output_eigvec = eigvec_sorted(:, 1:num_signals);

% create the subarray based on those defined with the corresponding ASV
%subarray definition 1: 1 2 3 4 --> 2: 2 3 4 5 
% Define Subarray from eigenvectors of perturbed array
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
% Output Angle via Least Squares Method 
angle_out = asind(doa_positions/((2*pi)*d));


% ------------------------------------------ %

% end of function
end 