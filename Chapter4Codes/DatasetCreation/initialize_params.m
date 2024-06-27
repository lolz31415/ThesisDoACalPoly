%% Function: Initialize Parameters 

function [signal,noise,x_expected] = initialize_params(theta_signals,varargin)
% ------------------------------------------ % 
% Function Purpose: Initialize Array - Output Array Element Locations & Signal for Covariance Matrix Creation
%
% INPUTS:
% theta_signals: determines the number of signals 
% * N: number of elements in array
% * d_expected: spcaing between elements in the ULA
% * SNR_dB: converted to linear for use in creaing AWGN noise variance
% * signal_amplitude: 1 for all signals unless otherwise specified
% * randseed: random seed number used for repeatable numbers
% * num_snapshots: number of time snapshots to capture signal and noise 
% OUTPUTS: 
% signal: output signal for covariance matrix use later
% noise: output noise vector based on input SNR
% x_expected: Random values within specified range to be added to
%
% Note: Using Uniform Linear Array - set distance | Warn user for array spacing > 0.5 λ
% maybe in the future have controllable weights (then not ULA?) 
% Note: '*' means varargin input values to the function in the description
% ------------------------------------------ 
 
    % realize SNR by setting noise variance linear power     
    narginchk(1,7); % ensure 1-6 input arguments
    % Check Number of Input Arguments
    
    % define default values 
    default_N=5;  default_SNR_dB=20; default_random_seed=11;default_d_expected=0.5;default_snapshots=1000;
    % use -> var_name = varargin{1}; % assigning all varargin vars
    switch nargin
        case 1 % N, d_expected, SNR_dB, signal_amplitude, random_seed - not set
            N = default_N; % number of elements
            d_expected = default_d_expected; % spacing between elements
            SNR_dB = default_SNR_dB; 
            signal_amplitude = ones(length(theta_signal));
            randseed=default_random_seed;
            num_snapshots = default_snapshots;
        case 2 % d_expected, SNR_dB, signal_amplitude, random_seed - not set
            N = varargin{1};
            d_expected = default_d_expected; % spacing between elements
            SNR_dB = default_SNR_dB; 
            signal_amplitude = ones(length(theta_signal));
            randseed=default_random_seed;
            num_snapshots = default_snapshots;
        case 3 % SNR_dB, signal_amplitude, random_seed - not set
            N = varargin{1}; 
            d_expected = varargin{2}; 
            SNR_dB = 20; 
            signal_amplitude = ones(length(theta_signal));
            randseed=default_random_seed; 
            num_snapshots = default_snapshots;
        case 4 % signal_amplitude, random_seed - not set
            N = varargin{1}; 
            d_expected = varargin{2}; 
            SNR_dB = varargin{3};           
            signal_amplitude = ones(length(theta_signal));
            randseed=default_random_seed; 
            num_snapshots = default_snapshots;
        case 5 % random_seed - not set
            N = varargin{1}; 
            d_expected = varargin{2}; 
            SNR_dB = varargin{3};           
            signal_amplitude = varargin{4}; 
            randseed=default_random_seed;
            num_snapshots = default_snapshots;
        case 6 % num_snapshots - not set
            N = varargin{1}; 
            d_expected = varargin{2}; 
            SNR_dB = varargin{3};           
            signal_amplitude = varargin{4}; 
            randseed=varargin{5}; 
            num_snapshots = default_snapshots;
        case 7 % set variables to varargin values
            N = varargin{1}; 
            d_expected = varargin{2}; 
            SNR_dB = varargin{3};           
            signal_amplitude = varargin{4}; 
            randseed=varargin{5}; 
            num_snapshots = varargin{6};
        otherwise 
            error('Invalid number of input arguments.');
    end
    
    if d_expected > 0.5
        warning('Array Spacing'+char(8805)+'0.5'+char(955)); % char(8805) = ≥ char(955) = λ
    end
    
    %{
    k=2*pi; % Spatial Wavelength 2*pi/lambda -> for ASV: exp(j*k*x_perturbe*sin(theta_signal_rad)
    w=ones(1,N); % Element weights originally uniform all set to 1 
    %}

    % user applied "random_seed"
    rand('seed',randseed);randn('state',0); % create a repeatable random seed for signal creation

    x_expected = -N/2*d_expected+d_expected/2:d_expected:N/2*d_expected-d_expected/2; % Create matrix of element position in array
    signal=sign(randn(length(theta_signals),num_snapshots)).*(signal_amplitude'*ones(1,num_snapshots)); % Creation of either +/-1 1000 sample signal using "sign" and "randn"
    % Gaussian noise - no correlated noise or phase/gain offset
    % noise power [dB] = 10^(-SNR[dB]/10)
    noise_variance = 10^(-SNR_dB/10); % noise variance levels σ^2 (sig0)
    noise=sqrt(noise_variance)*randn(N,num_snapshots); % creating a Gaussian RV for noise scaled by σ

end