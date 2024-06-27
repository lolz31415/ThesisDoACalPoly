function [music_doa_estimation,angle_out] = MUSIC_DOA_Estimate(signal_in,varargin)
% ------------------------------------------ % 
% Obtain the output MUSIC Spectrum from inputs:
% INPUTS:
%
% signal_in: Signal matrix with K snapshots signal_in (signal_in = s + n)
% rand_seed: Random Seed to use: Default = 11 (for repeatability) [varargin]
% search_space Input Search Space: Default -90° to 90° - 181 samples (1° resoltuion)
% theta_in: input signal(s) array in degrees [used for ASV]
% element_pos: Element Position Array: Default 0.5λ spacing 5 element 
% SNR_in: Input SNR in dB ex: 20dB -> sig0 = 0.01 (noise power): Default 20dB
% SNR EX2: dB = 10*log10(x) x = 0.01 -> dB = -20 for noise
% 
% OUTPUTS: 
% 
% music_doa_estimation: 
% angle_out: output array of angles calculated by MUSIC
% search_space: angular search used by MUSIC - default 
%
% apparently "Note:" at the front makes it look different
% Note: varargin is a cell array - but requires consistent input order 
% Note: ASV = Array Steering Vector (exponential term describing phase
% shift between elements in Uniform Linear Array (ULA))
% ------------------------------------------ % 

% Parse Number of input arguments 
% error if none or > 7 arguments provided

if nargin == 0 || nargin > 6
    error("ERROR: Incorrect number of parameters input, view the comments in the function");
end
% default rand_seed & theta_in & element_pos & SNR
if nargin == 1
    rand_seed = 11;
    % L is number of points for search space
    L = 181; % 181 = 1° resolution ... 
    search_space = linspace(-90,90,L);
    theta_in = 30; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default search_space & theta_in & element_pos & SNR
if nargin == 2
    rand_seed = varargin{1};
    % L is number of points for search space
    L = 181; % 181 = 1° resolution ... 
    search_space = linspace(-90,90,L);
    theta_in = 30; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default theta_in & element_pos & SNR
if nargin == 3
    rand_seed = varargin{1};
    search_space = varargin{2};
    theta_in = 30; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default element_pos & SNR 
if nargin == 4
    rand_seed = varargin{1};
    search_space = varargin{2};
    theta_in = varargin{3}; 
    % # elements and spacing btwn elements
    N = 5; d = 0.5;
    element_pos = [-N/2*d+d/2:d:N/2*d-d/2];
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% default SNR 
if nargin == 5
    rand_seed = varargin{1};
    search_space = varargin{2};
    theta_in = varargin{3};
    element_pos = varargin{4};
    % noise power [dB] = 10^(-SNR[dB]/10)
    sig0 = 10^(-20/10);
end

% all defined for nargin == 6 - also output ASV for user 
if nargin == 6
    rand_seed = varargin{1};
    search_space = varargin{2};
    theta_in = varargin{3};
    element_pos = varargin{4};
    sig0 = 10^(-varargin{5}/10);
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
A=exp(1j*(2*pi)*element_pos_expected'*sin(theta_in_rad));

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

% create new ASV from element position - depends on imperfections
theta_search_rad = search_space*pi/180;
% output from search_space: N x L (ang resolution size) 
A_new = exp(1j*(2*pi)*element_pos'*sin(theta_search_rad));

% MUSIC Algorithm for determining DoA
num_signals = length(theta_in); 
% want only minimum eigenvalue to determine eigenvectors to use as element
% weights in the array 
noise_eigenvalue = min(diag(eigval));


% find noise eigenvalue column position in original eigenvalue matrix
[~,min_eigenval_idx] = find(noise_eigenvalue==eigval);
% using known column of minimum (noise) eigenvalue determine noise eigenvector
noise_eigenvector = eigvec(:,min_eigenval_idx);

% size N x L (search_space size 181 for 1° resolution) 
% find music spatial spectrum 

%music_doa_estimation =  1./abs(A_new'*noise_eigenvector).^2;
% sum across dimension 2 of the array 
% music_doa_estimation = sum(abs(A_new'*noise_eigenvector).^2,2)+eps(1); %ep 9.43 Van Trees
music_doa_estimation = 1./abs(A_new'*noise_eigenvector).^2; %+eps(1); %ep 9.43 Van Trees
music_doa_estimation = 10*log10((abs(music_doa_estimation)));% in dB

% element-wise transpose of sqrt of MUSIC output spectrum linear units 
% music_doa_estimation = sqrt(1./music_doa_estimation).';

% spectrum_angle = search_space; % MATLAB calls these spec and scanAng

% ------------------------------------------ % 
% Another Example:
% https://github.com/tanweer-mahdi/MUSIC-and-root-MUSIC/blob/master/MUSIC.m
% ZC = Parameter search space. For DOA estimation, each columns are potential array steering vectors
% gains(i) = (ZC(:,i)'*ZC(:,i))/(ZC(:,i)'*(eign*eign')*ZC(:,i));
% ------------------------------------------ % 

% find num_signals # of peaks in the inverse of the output values
% Use largest peaks after since NP/num_signals might yield a peak that is
% not the largest (local peak/maximum)
[peaks,locations] = findpeaks(music_doa_estimation,'SortStr','descend');%,num_signals); 
MUSIC_pk_num_signals_vals = maxk(peaks,num_signals);
% Determine if the number of found DOA "locations" is less than the number
% of signals which is when
D = min(num_signals,length(locations));
% disp(locations(1:D)); % debug: determine where in array answer is
% disp(search_space(locations(1:D))); % debug: Output the determined DOA

% assert throws error if condition is false: somehow more DOAs than
% num_signals (impossible)
assert(D <= num_signals);

% If D is valid then determine the output DOA locations 
% iterate to num_signals not D 
if D>0
    pk_val_indx = zeros(1,num_signals); actual_sig_pk_loc = zeros(1,num_signals);
    for i = 1:num_signals
        pk_val_indx(i) = find(MUSIC_pk_num_signals_vals(i)==peaks);
        actual_sig_pk_loc(i) = locations(pk_val_indx(i));
    end
    % plug the locations into the search space 
    angle_out = search_space(actual_sig_pk_loc);
else
    % outputs 1x0 empty double row vector
    angle_out = zeros(1,0);
end

% end of function
end 