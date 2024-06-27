%% Using 1 input source: Determine the absolute theta error for a single set of array perturbations
% Define incoming signals & Create element array positions
close all; clear all; clc; 
% setup file for Matlab figures
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',18)
legend_font_size = 12;

% number of snapshots / signals acquired for each element in array  
K = 1000; 
theta_in = -15; % Original is 50° but -75° showed worse absolute error
num_signals = length(theta_in);
% uniform signal input
signal_amplitude=ones(1,num_signals);
% BPSK input signal for K snapshots (N x K)
signal_in = sign(randn(num_signals,K)).*(signal_amplitude'*ones(1,K));  
% initialize constant random number generator (rng) for this file
rand_seed =11; rand('seed',rand_seed);randn('state',0)
L = 181; search_space = linspace(-90,90,L);
N = 5; d =0.5; 
% create expected element location array
element_pos_expected = -N/2*d+d/2:d:N/2*d-d/2;
% define size of element perturbation maaaaaatrix and positions after perturbation
element_perturbation = zeros(1,N); element_position = zeros(1,N);
% define input SNR in dB: 0.01 for gaussian random/white noise power level
SNR_in = 20;

% Samples Per perturbation_multiplier
samples_per_multiplier = 50; 
% define the multipliers for the random element perturbations (mult.)
perturbation_multipliers = 0.1:0.1:0.5; 
num_perturb_mults = length(perturbation_multipliers);

% pre-allocate arrays to size based on (mult., samples per mult.)
element_perturbation_avg = zeros(length(perturbation_multipliers),samples_per_multiplier);
MUSIC_angle_out = zeros(length(perturbation_multipliers),samples_per_multiplier);
MUSIC_theta_error = zeros(length(perturbation_multipliers),samples_per_multiplier);
ESPRIT_angle_out = zeros(length(perturbation_multipliers),samples_per_multiplier);
ESPRIT_theta_error = zeros(length(perturbation_multipliers),samples_per_multiplier);
% iterate through each perturbation multiplier 
for j = 1:num_perturb_mults
    % iterate through # samples per multiplier
    for k = 1:samples_per_multiplier
        % add perturbation to each element in array
        for i=1:N
            % random distance offset scaled by current perturbation multiplier
            element_perturbation(i) = perturbation_multipliers(i)*randn;  
            % random element location perturbed from original location x_expected
            element_position(i) = element_pos_expected(i)+element_perturbation(i); 
        end
    % take MUSIC spectrum
    [music_doa_estimation,MUSIC_angle_out_current] = MUSIC_DOA_Estimate(signal_in,rand_seed,search_space,...
        theta_in,element_position,SNR_in);
    ESPRIT_angle_out_current = ESPRIT_DOA_Estimate(signal_in,rand_seed,theta_in,element_position,SNR_in);

    % if no peakfinder output, then the determined angle is NaN
    if isempty(MUSIC_angle_out_current)
        MUSIC_angle_out_current = NaN;
    end
    % save output angle, error, and avg. perturbation 
    element_perturbation_avg(j,k) = 1/length(element_perturbation)*sum(abs(element_perturbation)); 
    ESPRIT_angle_out(j,k) = ESPRIT_angle_out_current;
    ESPRIT_theta_error(j,k) =  abs(ESPRIT_angle_out(j,k) - theta_in);
    MUSIC_angle_out(j,k) = MUSIC_angle_out_current; 
    MUSIC_theta_error(j,k) = abs(MUSIC_angle_out(j,k) - theta_in);
    end
end

avgPerturb_MUSIC_ESPRIT_fig = figure(100);
for j = 1:length(perturbation_multipliers)
    % plot samples_per_multiplier in each series
    scatter(element_perturbation_avg(j,:),MUSIC_theta_error(j,:),'b');
    hold on;
    scatter(element_perturbation_avg(j,:),ESPRIT_theta_error(j,:),'r');
end
% theta didn't fit into margins for title
%title("|\theta_{error}| vs. Avg Perturbation for \theta = "+string(theta_in)+"^\circ");
title("|\theta_{error}| vs. Avg Perturbation");

xlabel("Average Perturbation (\lambda)"); ylabel("| \theta_{error} |");
% legend(string(perturbation_multipliers))
legend(["MUSIC","ESPRIT"]);

setTightMargins(0.05); % Normally 0.05
% Plotting average error of these samples in: 
% "PlotAvgPerturbationAvgTrialsVSThetaError.m" 
num_samples = length(perturbation_multipliers)*samples_per_multiplier;
saveas(avgPerturb_MUSIC_ESPRIT_fig,...
    "MUSIC_and_ESPRIT_performanceVSAvgPerturbation_"+string(num_perturb_mults)+"SampleNumber_"+string(num_samples)+"DOASingleSource_"+string(theta_in)+".png");

%% Run algorithm 
% Input perturbed elements array into MUSIC: "element_position"
% "music_doa_estimation" are linear power outputs of MUSIC noise subspectrum
% even with perturbed array

%[music_doa_estimation,angle_out] = MUSIC_DOA_Estimate(signal_in,rand_seed,search_space,...
%    theta_in,element_position,SNR_in);

% With IDEAL Positions for Testing Purposes
%[music_doa_estimation,angle_out] = MUSIC_DOA_Estimate(signal_in,rand_seed,search_space,...
%    theta_in,element_pos_expected,SNR_in);
%disp("Output Angle: "+string(angle_out));

% Base Error: (Measured - Expected) / Expected
% Take Euclidean Error for > 1 Source - when get there
%theta_error = angle_out-theta_in;
%disp("Error: "+string(theta_error));

% Output Avg. Element perturbation
%disp("Avg Element perturbation: "+string(element_perturbation_avg));

%% Average Element Error for MUSIC and ESPRIT - 1 source input
% Create x (start with K) samples 
% Take the average element error for magnitude |0.1λ| to |0.25λ|
% use the array element locations to produce and estimated angle of arrival for each algorithms 
% Take Error and plot vs. average element position error
% Add-in Root-MUSIC 

% Previous File Name: "PlotAvgperturbationVSAngleEstimationError" 

% ----------------------------------------------------- % 
% Now take "pertbutation_multiplier" from 0.1λ-0.5λ do 100 samples 
% store "element_perturbation_avg" ,"angle_out" measured, "theta_error"