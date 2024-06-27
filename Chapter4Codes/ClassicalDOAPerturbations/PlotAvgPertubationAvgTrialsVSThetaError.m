%% Average Element Error for MUSIC and ESPRIT - 1 source input
% Create x (start with K) samples 
% Take the average element error for magnitude |0.1λ| to |0.25λ|
% use the array element locations to produce and estimated angle of arrival for each algorithms 
% Take Error and plot vs. average element position error
% Add-in Root-MUSIC 

% Previous File Name: "DifferentPertubationMultipliersVSThetaError" 

% ----------------------------------------------------- % 
% Now take "pertbutation_multiplier" from 0.1λ-0.5λ do 100 samples 
% store "element_perturbation_avg" ,"angle_out" measured, "theta_error"


%% Define incoming signals & Create element array positions
% Definine perturbations, multipliers, input signals, and MUSIC Spectrum
% Plot MUSIC Spectrum outputs wrt Array Perturbation Value Averages

close all; clear all; clc; 
% number of snapshots / signals acquired for each element in array  
K = 1000; 
theta_in = 10;
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
% define size of element perturbation matrix and positions after perturbation
element_perturbation = zeros(1,N); element_position = zeros(1,N);
% define input SNR in dB: 0.01 for gaussian random/white noise power level
SNR_in = 20;

% Samples Per perturbation_multiplier
samples_per_multiplier = 5000; 

% ----------------------------------------- % 
% define the multipliers for the random element perturbations (mult.)
% Change based on desired resolution of plots: 
% Can't view perturbation multiplier values on labeled plots past about 15

%perturbation_multipliers = 0.1:0.1:0.5; % 5 multipliers
perturbation_multipliers = 0.05:0.05:0.5; % 10 multipliers
%perturbation_multipliers = 0.01:0.01:0.5; % 50 multipliers
%perturbation_multipliers = 0.005:0.005:0.5; % 100 multipliers

% ----------------------------------------- % 

% pre-allocate arrays to size based on (mult., samples per mult.)
element_perturbation_avg = zeros(length(perturbation_multipliers),samples_per_multiplier);
angle_out = zeros(length(perturbation_multipliers),samples_per_multiplier);
theta_error = zeros(length(perturbation_multipliers),samples_per_multiplier);

% iterate through each perturbation multiplier 
for j = 1:length(perturbation_multipliers)
    % iterate through samples_per_multiplier variable
    for k = 1:samples_per_multiplier
        % add perturbation to each element in array
        for i=1:N
            % random distance offset scaled by 0.1
            element_perturbation(i) = perturbation_multipliers(j)*randn;              
            % random element location perturbed from original location x_expected
            %element_position(i) = element_pos_expected(i)+element_perturbation(i); 
            
            % Try just adding absolute value of perturbations? - this can also be done outside loop 
            %element_position(i) = element_pos_expected(i)+abs(element_perturbation(i)); 
            element_position(i) = element_pos_expected(i)+(element_perturbation(i)); 

        end
    
        % Troubleshooting Printing 
        %element_perturbation % apply absolute value of perturbation to position 

        % Troubleshooting Printing 
        %element_position
        
        % -------------------------------------------------- % 
        % Find Avg distance between elements and use for plot 
        % pre-allocate size then use this to determine avg element spacing
        %{
        element_distance = zeros(1,N-1); 
        for m=1:(N-1)
            % difference between elements for new spacing
            element_distance(m)= element_position(m+1)-element_position(m);    
        end

        % Troubleshooting Printing 
        element_distance
        %}
        % -------------------------------------------------- % 

        % take MUSIC spectrum for each sample_per_multiplier
        [music_doa_estimation,angle_out_current] = MUSIC_DOA_Estimate(signal_in,rand_seed,search_space,...
        theta_in,element_position,SNR_in);
        
        % if no peakfinder output, then the determined angle is NaN
        if isempty(angle_out_current)
            angle_out_current = NaN;
        end
        % save output angle, error, and avg. perturbation 
        % avg perturbations: perturbation_multipliers (row) x samples_per_multiplier (col)
        
        % Try just taking a sum of the absolute value of the element_perturbations 
        %element_perturbation_avg(j,k) = 1/length(element_perturbation)*sum((element_perturbation));
        element_perturbation_avg(j,k) = 1/length(element_perturbation)*sum((abs(element_perturbation)));
        

        % element_perturbation_avg(j,k) = 1/length(element_perturbation)*sum(element_perturbation); 
        % each column angle_out is 1 sample per multiplier each row is 1 perturbation multiplier
        angle_out(j,k) = angle_out_current; 
        theta_error(j,k) = angle_out(j,k) - theta_in;
    end
end

%% Average Theta Error 5000 Samples MUSIC - Variable Number Average Perturbation Multipliers
% Take Absolute Value of DOA Estimate Error and use for averaging:


% pre-allocate arrays
avg_element_perturbation_multiplier_MUSIC = zeros(1,size(element_perturbation_avg,1));
avg_theta_error_multiplier_MUSIC = zeros(1,size(theta_error,1));

% Take average over samples_per_multiplier samples for each multiplier
for j = 1:length(perturbation_multipliers)
    % convert 2D matrix to 1D
    angle_out_current = angle_out(j,:);
    element_perturbation_avg_current = element_perturbation_avg(j,:);
    % try absolue value of each theta_error for samples_per_multiplier
    theta_error_current = abs(theta_error(j,:));

    % check and remove NaN results for each multiplier
    if ~isempty(find(isnan(angle_out_current))) 
        nan_index = find(isnan(angle_out_current));
        % can apply [] since now 1D Matrix 
        angle_out_current(nan_index) = [];
        % do for element_perturbation_avg_current and theta_error_current
        theta_error_current(nan_index) = [];
        element_perturbation_avg_current(nan_index) = [];
    end
    % don't have to try absolute value for perturbation here since previously done 
    %avg_element_perturbation_multiplier_MUSIC(j) = sum(element_perturbation_avg_current)/length(element_perturbation_avg_current);
    
    avg_element_perturbation_multiplier_MUSIC(j) = sum(abs(element_perturbation_avg_current))/length(element_perturbation_avg_current);

    avg_theta_error_multiplier_MUSIC(j) = sum(theta_error_current)/length(theta_error_current);
end

% Plot them together if greater than 10 multipliers used 
if length(perturbation_multipliers) <= 10
    
    % Create Plot
    avg_MUSIC_error_fig = figure(100);
    %scatter(avg_element_perturbation_multiplier_ESPRIT,avg_theta_error_multiplier_ESPRIT);
    % abs
    scatter(avg_element_perturbation_multiplier_MUSIC,abs(avg_theta_error_multiplier_MUSIC));
    %scatter(perturbation_multipliers,abs(avg_theta_error_multiplier_MUSIC));
    
    title("|\theta_{error}| vs. Avg. Perturb. MUSIC")
    xlabel('Average Perturbation (\lambda)'); ylabel("|\theta_{error}|");
    
    % Add space to not cut off text in plot 
    extraSpace = 0.01; 
    xlim([min(avg_element_perturbation_multiplier_MUSIC) - extraSpace, ...
        max(avg_element_perturbation_multiplier_MUSIC) + extraSpace]);
    
    % Labeling each point -> Need to shift points to left by dx to ensure all fit on plot
    % And aren't being block by axis labels (in text section)
    xpos = avg_element_perturbation_multiplier_MUSIC; ypos = abs(avg_theta_error_multiplier_MUSIC);    
    dx = 0.0025; dy = 1;
    text(xpos-dx,ypos-dy,sprintfc('%.2f',perturbation_multipliers),FontSize=12);    
    % Saving Figures
    setTightMargins(0.275); % usually 0.05
    saveas(avg_MUSIC_error_fig,"MUSICAvgErrorForAvgPerturbMults_"+string(length(perturbation_multipliers))+"SamplesEachPerturbMult_"+string(samples_per_multiplier)+".jpg");
end

% Finding values that are NaN - didn't find a DoA:
% K = find(isnan(angle_out)); % outputs locations in angle_out to find NaN

%% Using ESPRIT -> Note: Array Perturbations Immediately in ASV for Eigenvalue Decomposition 
% ESPRIT should have a linear trend for error since direct dependence on ASV
% Determining ESPRIT Capability to deal with array perturbations 
% Use the same: "samples_per_multiplier" variable as in MUSIC

% Same # multipliers avg perturbation as for MUSIC
perturbation_multipliers_ESPRIT = perturbation_multipliers;

% reminder: theta_in = 50°
% element_pos_expected -> 50.0978°

% --------------------------------------------- %

% USING SAME VARIABLE NAMES / OVERALL FUNCTIONALITY 
% Turn this part into a function? 

% pre-allocate arrays to size based on (mult., samples per mult.)
element_perturbation_avg_ESPRIT = zeros(length(perturbation_multipliers_ESPRIT),samples_per_multiplier);
angle_out_ESPRIT = zeros(length(perturbation_multipliers_ESPRIT),samples_per_multiplier);
theta_error_ESPRIT = zeros(length(perturbation_multipliers_ESPRIT),samples_per_multiplier);

% iterate through each perturbation multiplier 
for j = 1:length(perturbation_multipliers_ESPRIT)
    % iterate through # samples per multiplier
    for k = 1:samples_per_multiplier
        % add perturbation to each element in array
        for i=1:N
            % random distance offset scaled by 0.1
            element_perturbation(i) = perturbation_multipliers_ESPRIT(j)*randn;  
            % random element location perturbed from original location x_expected
            element_position(i) = element_pos_expected(i)+element_perturbation(i); 
        end
        % take ESPRIT spectrum - angle_out_current only possible output of ESPRIT function 
        angle_out_current_ESPRIT = ESPRIT_DOA_Estimate(signal_in,rand_seed,...
            theta_in,element_position,SNR_in);
    % if no peakfinder output, then the determined angle is NaN
    if isempty(angle_out_current_ESPRIT)
        angle_out_current_ESPRIT = NaN;
    end
    % save output angle, error, and avg. perturbation 
    % Without Absolute Value: Lower Value
    %element_perturbation_avg(j,k) = 1/length(element_perturbation)*sum(element_perturbation); 
    % With Absolute Value
    element_perturbation_avg_ESPRIT(j,k) = sum(abs(element_perturbation)) / length(element_perturbation); 
    angle_out_ESPRIT(j,k) = angle_out_current_ESPRIT; 
    theta_error_ESPRIT(j,k) = angle_out_ESPRIT(j,k) - theta_in;
    end
end

% Average 100 samples over each multiplier (no ABS?!) 
% If this works do more multipliers 0.05 (10) then 0.01 (50) change
% starting to the given multiplier number


% pre-allocate arrays
avg_element_perturbation_multiplier_ESPRIT = zeros(1,size(element_perturbation_avg_ESPRIT,1));
avg_theta_error_multiplier_ESPRIT = zeros(1,size(theta_error_ESPRIT,1));

% Take average over 100 samples for each multiplier
for j = 1:length(perturbation_multipliers_ESPRIT)
    % convert 2D matrix to 1D
    angle_out_current_ESPRIT = angle_out_ESPRIT(j,:);
    element_perturbation_avg_current_ESPRIT = element_perturbation_avg_ESPRIT(j,:);
    theta_error_current_ESPRIT = abs(theta_error_ESPRIT(j,:));

    % check and remove NaN results for each multiplier
    if ~isempty(find(isnan(angle_out_current_ESPRIT))) 
        nan_index = find(isnan(angle_out_current_ESPRIT));
        % can apply [] since now 1D Matrix 
        angle_out_current_ESPRIT(nan_index) = [];
        % do for element_perturbation_avg_current and theta_error_current_ESPRIT
        theta_error_current_ESPRIT(nan_index) = [];
        element_perturbation_avg_current_ESPRIT(nan_index) = [];
    end
    avg_element_perturbation_multiplier_ESPRIT(j) = sum(element_perturbation_avg_current_ESPRIT)/length(element_perturbation_avg_current_ESPRIT);
    avg_theta_error_multiplier_ESPRIT(j) = sum(theta_error_current_ESPRIT)/length(theta_error_current_ESPRIT);
end



% ---------------------------------------------- %
% Y Vs. X = Plot Titles

if length(perturbation_multipliers) <= 10
    
    % Create Plot
    avg_ESPRIT_error_fig = figure(200);
    % scatter(avg_element_perturbation_multiplier_ESPRIT,avg_theta_error_multiplier_ESPRIT);
    % take absolute value of error
    scatter(avg_element_perturbation_multiplier_ESPRIT,abs(avg_theta_error_multiplier_ESPRIT));
    title("|\theta_{error}| vs. Avg. Perturb. ESPRIT")
    xlabel('Average Perturbation (\lambda)'); ylabel("|\theta_{error}|");
    
    % Add space to not cut off text in plot - wider than 0.002 for MUSIC  
    extraSpace = 0.01; 
    xlim([min(avg_element_perturbation_multiplier_ESPRIT) - extraSpace, ...
        max(avg_element_perturbation_multiplier_ESPRIT) + extraSpace]);
    
    % Labeling each point 
    % Move x-positions since 0.05 being cut-off by setting dx = 0.0025
    xpos = avg_element_perturbation_multiplier_ESPRIT; ypos = abs(avg_theta_error_multiplier_ESPRIT); 
    dx = 0.0025; dy = 0.15;
    text(xpos-dx,ypos-dy,sprintfc('%.2f',perturbation_multipliers_ESPRIT),FontSize=12);
    setTightMargins(0.275); % normally 0.05
    saveas(avg_ESPRIT_error_fig,"ESPRITAvgErrorForAvgPerturbMults_"+string(length(perturbation_multipliers_ESPRIT))+"SamplesEachPerturbMult_"+string(samples_per_multiplier)+".jpg");

end

%% If > 10 perturbation multipliers: plot MUSIC with ESPRIT

    overall_fig = figure(300);
    %scatter(avg_element_perturbation_multiplier_ESPRIT,avg_theta_error_multiplier_ESPRIT);
    % abs
    scatter(avg_element_perturbation_multiplier_MUSIC,abs(avg_theta_error_multiplier_MUSIC),'b'); hold on;
    scatter(avg_element_perturbation_multiplier_ESPRIT,abs(avg_theta_error_multiplier_ESPRIT),'r');


    title("|\theta_{error}| vs. Average Perturbation")
    xlabel('Average Perturbation (\lambda)'); ylabel("|\theta_{error}|");
    legend(["MUSIC", "ESPRIT"]);

    % Saving Figures
    setTightMargins(0.05);
    saveas(overall_fig,"OverallAvgErrorForAvgPerturbMults_"+string(length(perturbation_multipliers))+"SamplesEachPerturbMult_"+string(samples_per_multiplier)+".jpg");