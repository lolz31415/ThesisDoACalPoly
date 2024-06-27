%% TRY HIGH LEVEL H5 COMMANDS: 
% h5create, h5write, h5disp, h5read, h5info 
% Note: Still use low leve H5T (Datatype) Commands to make compound datatypes
%
% Just use Datasets for all variables to store, no groups 
% add the perturbation multiplier as a dataset so 7 total 
% For either use of HDF2Struct later or by hand convert to workspace
% 4D array? 
% Separate files for each perturbation multipliers -> Datastore? 

clear all; close all; clc; 

% file_id = H5F.create(name,flags,fcpl_id,fapl_id)
%filename = 'training_data_array_perturbations.h5';
filename = 'validation_data_array_perturbations.h5';
%filename = 'testing_data_array_perturbations.h5';

if exist(filename,"file")
    delete(filename);
end
%%%%% Line 52 change number of multipliers & Line 40 Change # Samples %%%%%%%%%
%% initialization 
% Automatically overwrites 
% How to close H5 file? | Script doesn't terminate

theta_signals = 30; N = 5; d_exp = 0.5; SNR_dB = 20; sig_amp = ones(1,length(theta_signals)); randseed=11; num_snapshots = 1000;
[signal,noise,x_expected] = initialize_params(theta_signals,N,d_exp,SNR_dB,sig_amp,randseed,num_snapshots);

% to help define dataset sizes for creation in HDF5
input_signals = length(theta_signals); 

% datatype to allow storing complex numbers 
%%% Create separate Validation and Testing Files by Using this Script and %%% changing a few values %%%
training_data = 1; % 100e3 samples per multiplier (set below) 
validation_data = 0; % 50e3 sampler per multiplier (set below) 
testing_data = 0; % 25e3 sampler per multiplier (set below) 

% samples_per_perturbation_mult = 100e3; % after test on small number
% 5e6 samples might be a lot to produce: Try making 100 for now. 
base_sample_number = 500; 
% didn't want to change variable 'training_sample_number' = samples per perturbation multipliers
if contains(filename,'training') 
    training_sample_number = base_sample_number;
elseif contains(filename,'validat')
    training_sample_number = base_sample_number/2; 
elseif contains(filename,'test')
    training_sample_number = base_sample_number/4;
end

% for example: 
% 25 pts: 0.02λ | 500pts: 0.001λ | 1000pts: 0.0005λ 
% originally 25
perturbation_multiplier_amount = 50; 
% Don't start with 0 -> trying to learn the mapping might mess it up 
% start with a small epsilon value of 0.005
epsilon_start = 0.005; 
perturbation_multipliers = linspace(epsilon_start,0.5,perturbation_multiplier_amount);

% Just do 1 dataset for each variable (7) in each dataset
% not doing groups for perturb mult so include into dataset 
theta_dataset_name = strcat("dataset","_theta");
arr_perturb_mult_dataset_name = strcat("dataset","_arr_perturb_mult"); 
arr_perturb_dataset_name = strcat("dataset","_arr_perturb_elem");
cov_real_dataset_name = strcat("dataset","_cov_real");
cov_imag_dataset_name = strcat("dataset","_cov_imag");
cov_abs_dataset_name = strcat("dataset","_cov_abs");
cov_phase_dataset_name = strcat("dataset","_cov_angle");
groups_num_sample_dataset_name = strcat("dataset","_groups_",string(perturbation_multiplier_amount)...
    ,"_num_samples_",string(training_sample_number));

%% Create File
% Open the HDF5 file (replace 'your_file.h5' with your desired filename)

percent_done = 0; 

% not required anymore using h5info can obtain this info: 
% dataset_names = [theta_dataset_name, arr_perturb_mult_dataset_name, arr_perturb_dataset_name, cov_real_dataset_name,...
%    cov_imag_dataset_name,cov_abs_dataset_name,cov_phase_dataset_name];

% write files to storage for datset_names
% If "string" was specified as the datatype in the corresponding call to h5create, data is a MATLAB string array. 
% The string array dimensions must match those specified in the call to h5create.

% Can create the H5 files outside of the loop to have 3rd dimension: 
% num_samples * perturbation_multiplier_amount

% store the perturbation multiplier and just store it
% should we shuffle training data or will the training process
% automatically grab data randomly? 
% Separate HDF5 Files for each perturbation multiplier? 

% set names for create/write here since group value changes
theta_name = strcat('/',theta_dataset_name); 
perturb_mult_name = strcat('/', arr_perturb_mult_dataset_name); 
perturb_name = strcat('/', arr_perturb_dataset_name); 
cov_real_name = strcat('/', cov_real_dataset_name);
cov_imag_name = strcat('/', cov_imag_dataset_name);
cov_abs_name = strcat('/', cov_abs_dataset_name);
cov_phase_name = strcat('/', cov_phase_dataset_name); 
store_multipliers_num_samples = strcat('/',groups_num_sample_dataset_name);

%h5create(filename,group_name, [training_sample_number,5,5]);
% Create manually each dataset for each group here so can be written to
total_num_samples = perturbation_multiplier_amount * training_sample_number;
h5create(filename, theta_name, [1, input_signals, total_num_samples]); 
h5create(filename, perturb_mult_name, [1, 1, total_num_samples]); 
h5create(filename, perturb_name, [1, N, total_num_samples]); 
h5create(filename, cov_real_name, [N, N, total_num_samples]); 
h5create(filename, cov_imag_name, [N, N, total_num_samples]); 
h5create(filename, cov_abs_name, [N, N, total_num_samples]); 
h5create(filename, cov_phase_name, [N, N, total_num_samples]);
% store & write num groups and num samples per group 
h5create(filename, store_multipliers_num_samples, [1, 2]);
h5write(filename, store_multipliers_num_samples,...
    [perturbation_multiplier_amount training_sample_number]);

start_write_time = tic;

for i = 1:perturbation_multiplier_amount
    % define current iteration of perturbation multipliers 
    current_perturbation_multiplier = perturbation_multipliers(i);
    % Define Percent Completion:
    percent_done = i/perturbation_multiplier_amount * 100; 
    fprintf("Iteration Percentage Complete: "+string(percent_done)+"%% \n");
    
    % loop through trainins samples for each perturbation multiplier
    for j = 1:training_sample_number
        
        % create output data
        [perturbations, covariance_matrix] = ...
            create_perturbed_cov_matrix(N,theta_signals,signal,...
            noise,d_exp,num_snapshots,current_perturbation_multiplier,x_expected);
       % store output data in struct format for each sample 
       output_data.theta = theta_signals; % 1 x length(theta_signals)
       output_data.perturbationMultiplier = current_perturbation_multiplier; % 1 x 1
       output_data.elementPositionPerturbations = perturbations; % 1 x N
       output_data.realCovariance = real(covariance_matrix); % N x N
       output_data.imaginaryCovariance = imag(covariance_matrix); % N x N
       output_data.magnitudeCovariance = abs(covariance_matrix); % N x N
       output_data.phaseCovariance = 180/pi * angle(covariance_matrix); % N x N
       % write the output data to the appropriate dataset_name
       
       % h5write(FILENAME,DATASETNAME,DATA,START,COUNT)
       % note use (i-1)*training_sample_number+j for starting bounds 
       start_bounds = (i - 1) * training_sample_number + j; 
       h5write(filename, theta_name, output_data.theta, [1,input_signals,start_bounds], [1,input_signals,1]);
       h5write(filename, perturb_mult_name, output_data.perturbationMultiplier,[1,1,start_bounds], [1,1,1]);      
       h5write(filename, perturb_name, output_data.elementPositionPerturbations,[1,1,start_bounds], [1,N,1]);
       h5write(filename, cov_real_name, output_data.realCovariance,[1,1,start_bounds], [N,N,1]);
       h5write(filename, cov_imag_name, output_data.imaginaryCovariance,[1,1,start_bounds], [N,N,1]);
       h5write(filename, cov_abs_name, output_data.magnitudeCovariance,[1,1,start_bounds], [N,N,1]);
       h5write(filename, cov_phase_name, output_data.phaseCovariance,[1,1,start_bounds], [N,N,1]);       
    end
    fprintf(string(training_sample_number)+ ...
        ' Samples have been written to dataset! New Perturbation Multiplier \n');
end

stop_write_time = toc(start_write_time);

current_directory = pwd; 
file_location = strcat(current_directory,filename);
fprintf("File Creation Complete! Your File can be found at: \n %s \n",file_location);
disp("Time to Write Dataset: "+string(stop_write_time));

%% Read Outputs to File

%read_file = h5info('training_data_array_perturbations.h5'); 
%read_file_groups = read_file.Groups; 
%read_file_groups_dataset = read_file_groups.Datasets

% obtain dataset names to use for h5read 
% h5disp(filename); % provides information about the structure of the h5 file
data_info = h5info(filename);
% normal matrix doesn't work 
dataset_names = strcat('/',{data_info.Datasets.Name}); 

% obtain data from datasets within files 
% collect data in struct format
% use for loop and switch statement for 'contains' since dataset names are
% alphabetical not in order created - 7 case switch statement 

read_back_data = struct('theta', [], 'perturb_mult', [], ...
                       'perturb_elem', [], 'cov_real', [], ...
                       'cov_imag', [], 'cov_abs', [], ...
                       'cov_angle', []);

% switch / case requires EXACT match not == 1 to run code 

for l = 1:length(dataset_names)
    current_name = dataset_names{l};
    if contains(current_name,'theta')
            read_back_data.theta = h5read(filename,dataset_names{l});
    elseif contains(current_name,'perturb_mult')
            read_back_data.perturb_mult = h5read(filename,dataset_names{l});
    elseif contains(current_name,'perturb_elem')
            read_back_data.perturb_elem = h5read(filename,dataset_names{l});
    elseif contains(current_name,'cov_real')
            read_back_data.cov_real = h5read(filename,dataset_names{l});
    elseif contains(current_name,'cov_imag')
            read_back_data.cov_imag = h5read(filename,dataset_names{l});
    elseif contains(current_name,'cov_abs')
            read_back_data.cov_abs = h5read(filename,dataset_names{l});
    elseif contains(current_name,'cov_angle')
            read_back_data.cov_angle = h5read(filename,dataset_names{l});
    end
end
% why not writing to struct to store data 
% use HD5Datastore -> multiple H5 Files for each perturbation multiplier 

%% Functions Sections at bottom of script

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

%% Create 1 Sample Perturbed Array Data 

function [perturbations, covariance_matrix] = ...
    create_perturbed_cov_matrix(N,theta_signals,signal,noise,d_expected,num_snapshots,perturbation_multiplier,varargin)
% ------------------------------------------ % 
% Function Purpose: Output Array Perturbations & Covariance Matrix  
% *** Kind of clunky function, but accomplishes what we need it to do:
% Can call this function in a for loop to create dataset for CNN to utilize
%
% INPUTS:
% N: number of elements in array same as input to "initialize_params" function
% theta_signals: input signal locations that go into "initialize_params" function
% signal: output signal from "initialize_params" function
% noise: output noise from "initialize_params" function
% d_expected: distance between elements from "initialize_params" function
% num_snapshots(K): snapshots received used for processing (averaged signal over)
% perturbation_multiplier: maximum element perturbation 
% * x_expected: array element positions (assumed to be 0.5 if no input)
%
% Note: Perturbation Multipliers will already be known for output / labeling data 
% Note: '*' are varargin not required inputs to the function   
%
% OUTPUTS: 
% perturbations: Random values within specified range to be added to
% DO WE NEED THIS? [element_positions]
% covariance_matrix: produced covariance matrix from perturbations 
%   
% Note: can use narginchk(min,max) and nargin to use switch/case statement
% for # of inputs and what to set as default values
% ------------------------------------------ % 
    
    narginchk(7,8); % ensure 7-8 input arguments
    % Check Number of Input Arguments
    switch nargin
        case 7 % input_arr not defined use 0.5λ spacing 
            x_expected = -N/2*d_expected+d_expected/2:d_expected:N/2*d_expected-d_expected/2; % Create matrix of element position in array
        case 8 % assign to variable argument if 8 total inputs 
            x_expected = varargin{1}; 
        otherwise 
            error('Invalid number of input arguments.');
    end

    % When we make the covariance matrix, have it based on if the array elements were 0.5 spacing 
    d = zeros(1,N); x_perturbed = zeros(1,N);
    % We will create the Correlation Matrix Used for EVD then use expected values to plot the result and see MUSIC output in comparison
    for i=1:N
        % ensure values +/-0.1 and round to 3 decimal places
        %d(i) = round(0.1*randn,4);
        d(i) = perturbation_multiplier*randn;        % random distance offset scaled by perturbation_multiplier
        x_perturbed(i) = x_expected(i)+d(i); % random element location perturbed from original location x_expected
    end
    
    k = 2*pi/d_expected; % phase constant where d_expected is in λ 
    theta_signal_rad = theta_signals*pi/180; % Angle of incoming signal relative to array elements
    A=exp(1j*k*x_perturbed'*sin(theta_signal_rad)); % Signal Array Manifold Vector w/ perturbed element locations
    
    % R = a a^H / K take the number of snapshots K, to scale the correlation
    % matrix by, and multiply the hermitian matrix (complex conjugate) by the original value
    
    Rss=signal*signal'/num_snapshots;                 % signal correlation matrix
    Rnn=(noise*noise')/num_snapshots;                 %noise correlation matrix
    Rns=(noise*signal')/num_snapshots;                % noise/signal correlation matrix - should be near 0
    Rsn=(signal*noise')/num_snapshots;                % signal/noise correlation matrix - should be near 0 
    Rrr=A*Rss*A'+A*Rsn+Rns*A'+Rnn;                    % array correlation matrix

    % assign outputs of matrix 
    covariance_matrix = Rrr; 
    perturbations = d; 

end