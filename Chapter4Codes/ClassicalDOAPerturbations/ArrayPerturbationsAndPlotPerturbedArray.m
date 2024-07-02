% Script Description/Purpose: 
% Used for determining effect of snapshots, different angle for source
% location, different average perturbation multipliers for a single case,
% and applying perturbations or viewing effects on expected array.

%% Clear the workspace between successive runs
clear all; close all; clc;

% setup file for Matlab figures
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',18)
legend_font_size = 12;
%% Problem Creation and Eigenvalues/Eigenvectors created from Correlation Matrix
randseed=11; rand('seed',randseed);randn('state',0)
k=2*pi;

% Originally d=0.5; % inter-element spacing 
N=5;   % number of array elements
% w gets set later to uniform linear array
% w=ones(1,N); % element weights originally uniform all set to 1 

% Original Element Position 
% x=[-N/2*d+d/2:d:N/2*d-d/2];
% set up element position in for loop to reset spacing at each element
d = zeros(1,N); x = zeros(1,N); 

% Expected values prior to random perturbation application
d_exp = 0.5;
x_exp= -N/2*d_exp+d_exp/2:d_exp:N/2*d_exp-d_exp/2;

%% Configuration Values
% Creating Random values and changing element array positions 
% perturbations set if using perturbations or just wide angle -80° source

perturbations = 0; % <-- Value to set if desired view perturbations
average_perturbation_multiplier = 0.1;

% Original Value 5000
K=1000; % number of signal samples

% Original Value [-80]
ts=[-60]; % set angular direction of signal

% Configure SNR where: SNR = 10*log(1/sig0) - power ratio
sig0=1e-2; %sig0=1e-2;% initial value =1e-2 sig0 is noise power level (linear)


%% Configuring the RX signal snapshots, sample covariance, and EVD for use in Eigenbeams and MUSIC

for i=1:N
    % ensure values +/-average_perturbation_multiplier and round to 3 decimal places
    % For just using -80 out of range [-60,60] images
    if perturbations == 1
        d(i) = average_perturbation_multiplier*randn;
    else
        % no perturbations desired case: 
        d(i)=0;
    end

    % add the normal distributed value to the original array
    x(i) = x_exp(i)+d(i); 
end

fprintf('%s %f ','random distances: ',d); 
fprintf('\n'); % \n need semicolon so no output display 
fprintf('%s %f ','array locations: ',x)
fprintf('\n'); % \n need semicolon so no output display 



Ns=length(ts);  % number of signals
sa=ones(1,Ns); % set sig ampl and dir for num sig: Ns
% number of signals: cs
num_signals = length(ts);
if (((num_signals~=Ns) | (length(sa))~=Ns)) % if either is not equivalent to num input signals
    msg = "Check the size of the signal amplitude (sa) or direction (ts)";
    error(msg)
end
tsr=ts*pi/180; % convert input to rad
s=sign(randn(Ns,K)).*(sa'*ones(1,K)); % produce interfering signal 

% Gaussian noise
nois=sqrt(sig0)*randn(N,K);

% Correction: Compute the Covariance Matrix using the perturbed array and
% apply this to MUSIC, not adjust element positions after computing the ASV

% BELOW WRONG: 
% --- Use x_exp the expected array values instead of the real one and see how it affects the output ---
% Array Steering Vector (ASV) = Array Manifold Vector
% A=exp(1j*k*x_exp'*sin(tsr)); % --> Signal Array Manifold Vector
A=exp(1j*k*x'*sin(tsr)); % Create 
% In this example eigenbeams no interferer signals 

Rss=s*s'/K;      % signal correlation matrix
Rnn=(nois*nois')/K;         %noise correlation matrix
Rns=(nois*s')/K;           % noise/signal correlation matrix
Rsn=(s*nois')/K;           % signal/noise correlation matrix
Rrr=A*Rss*A'+A*Rsn+Rns*A'+Rnn;        % array correlation matrix

[eigvec, eigval] = eig(Rrr);

%% Plotting Values to create array steering vector for all angles and more angular resolution

% Find the output where A is multiplied by w and plot over theta
L = 181; % # of theta points to take 1 degree for this case --> Samples
theta = linspace(-90,90,L);
theta_rad = theta*pi/180;
% A_new allows plotting by splitting the Array Manifold Vector based on all
% angles not just maximum points 
A_new = exp(1j*k*x_exp'*sin(theta_rad)); % 1 x L 
% using the actual array element positions
%A_new = exp(1j*k*x'*sin(theta_rad)); % 1 x L 

%% Creating plots of the same size every time
% Plot Elements on the Aperture
xpos = 5; ypos = 1; % N = y-dir M = x-dir 
% See Above:
% x=[-N/2*d+d/2:d:N/2*d-d/2];
elementPos = x;
%fg=figure;
%fg.Position(3:4) = [435 255];

% Run an if statement to determine if the number of elements are enough or
% just make them centered on 0 on other axis
% M elements x N elements y row 11 elementPos for x
array_position_fig = figure(1);
array_position_fig.Name = "Array Positions";
sizeInput = size(elementPos);
if sizeInput(1) == 2
    scatter(elementPos(1,:),elementPos(2,:),'o','filled')
else
    elementPos(2,1:N) = zeros(1,N);
    scatter(elementPos(1,:),elementPos(2,:),'o','filled')
end

axis equal
title("Array Element Positions in terms of \lambda")
xlim([-(xpos.*d_exp)/2 (xpos.*d_exp)/2]); 
ylim([-(ypos*d(1)+1)/2 (ypos*d(N)+1)/2]);
xlabel('x(\lambda)')
ylabel('y(\lambda)')
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])

% Show Element Positions on Plot (for x-pos, y is just for visualization) 
element_pos_string = "Element Positions: "+join(string(elementPos(1,:)));
text(-1.15,-0.25,element_pos_string,FontSize=12);

% box on
setTightMargins(0.05);

% Define Path for Images to be sent to (so figures are not being saved to current directory) 
if perturbations == 1
    saveas(array_position_fig,'ArrayGeometry_ChangingArrayPosition'+join(string(ts))+'.png');
else
    saveas(array_position_fig,'ArrayGeometry_NO_Perturbations'+join(string(ts))+'.png');
end

%% Eigenbeam Plotting 

% Produce the weights for each set of signals by taking 1 column from the
% eigenvectors where the 2 largest are the signal weights
% Pre-allocate Arrray Sizes
w = ones(N); 
y = ones(N,L);
y_dB = y;
eigenbeam_fig = figure(10);
eigenbeam_fig.Name = "Eigenbeams";
% eigval from perturbed array's covariance matrix
eigen_values_plot = diag(real(eigval))';

% subtitle_eigval = join(eigen_values_plot); % gives all eigen values in
% one string delminated by a space
% NOTE: Ensure elements are of the desired data type in MATLAB when
% preallocating (only strings for something like a title below)
subtitle_eigval = string(ones(1,N)); % pre-allocate array

% subtitle_eigval is just eigen_values_plot with some extra text to explain corresponding eigenbeam

% y = w * A_new multiplies eigenvectors by A_new where A_new subdivides in
% 1° sections so plotting can occur and samples for each of the 5 sensors can be viewed
% so A_new here is the MUSIC spectrum obtained - use the A from covariance
% matrix produced by perturbed element locations
legend_vals = [];

for i = 1:N   
    subtitle_eigval(i) = string(i) + ": " + string(eigen_values_plot(i));
    w(:,i) = eigvec(:,i);
    y(i,:) = w(:,i)'*A_new;
    y_dB(i,:) = 20*log10(abs(y(i,:)));
    plot(theta,y_dB(i,:)); hold on;
    legend_vals = [legend_vals "EigenBeam "+string(i)];
end

grid on; 
% Since we know N signals Ns-sa noise eigenbeams: 
P = N - length(sa); % number of noise eigenbeams
[smallest_val, indx] = mink(eigen_values_plot,P)
% use the two smallest eigenvalues to plot line at their minimum point

for i = 1:P
    %peak_vals(P,:)
    [~,locations(P,:)] = findpeaks(-y_dB(indx(P),:)) % find peaks for - of data for minimum values  
    % Plot Minimum Values for each Noise EigenBeam
    for k = 1:length(locations(P,:))
        theta_min = theta(locations(P,k)); 
        xline_var = xline(theta_min,':k',{"\theta = "+string(theta_min)+"^\circ"},LineWidth=2); % second argument is line style
        xline_var.FontSize = 18;
        hold on;
    end
end
hold off; 

% Determine best orientation and Location  of legend
legend(legend_vals,Location="best",FontSize=legend_font_size,Orientation="vertical");
title('Eigenbeam Plots');

% Max Eigenbeam at desired DOAs is the desired eigenvector = use as weights 
subtitle_txt = ["Eigenbeam: "+join(subtitle_eigval(1)+" | "+subtitle_eigval(2)),...
                "Eigenbeam: "+join(subtitle_eigval(3)+" | "+subtitle_eigval(4)),...
                "Eigenbeam: "+subtitle_eigval(5)];
subtitle(subtitle_txt);
xlabel("\theta^\circ")
ylabel("AF(dB)")
setTightMargins(0.05)

if perturbations == 0
    % For large angle outside of [-60°,60°] no perturbations:
    saveas(eigenbeam_fig,'Eigenbeam_NO_Perturbations'+join(string(ts))+'.png');
else
    saveas(eigenbeam_fig,'Eigenbeam_ChangingArrayPosition'+join(string(ts))+'.png');
end

%% MUSIC Algorithm DoA

%%% New Corrected Noise Eigenvector Determination %%%

% There are P = N-M for N = # array elements, M = # signals noise
% eigenvalues/vectors
% Pick from the minimum of these to create the MUSIC spectrum for best
% results
noise_eigval = mink(diag(eigval),P); % cs = # of signals
num_DOAs = num_signals;

% Don't need to iterate, just use minimum eigenvalue of the P noise
% eigenvalues and use corresponding single minimum eigenvector to create 
% MUSIC spectrum
min_eigenval = min(noise_eigval);
% find noise eigenvalue column position in original eigenvalue matrix
[~,min_eigenval_idx] = find(min_eigenval==eigval);
% using known column of minimum (noise) eigenvalue determine noise eigenvector
noise_eigvec = eigvec(:,min_eigenval_idx);
%%% End Revised/New Section %%%


P_theta = 1./abs(A_new'*noise_eigvec).^2;
P_dB = 10*log10(abs(P_theta));

MUSIC_fig = figure(100);
MUSIC_fig.Name = "MUSIC Array Perturbations Single Source";
plot(theta,P_dB,'b');
title("MUSIC Algorithm");
xlabel('\theta \circ');
ylabel("P(\theta) [dB]");

%%% New Peak Finder: Allow just known number of maximum peaks to use %%%
% Find Peaks / DOA on MUSIC Sub-Spectra 
% Use maxk to only grab the num_signals (cs) (1) largest signals 
% DOA determined 4 locations: -83°,-42°,1°,24°
[MUSIC_pk,MUSIC_loc] = findpeaks(P_dB);
MUSIC_pk_num_signals_vals = maxk(MUSIC_pk,num_signals);
% Pre-allocate:
pk_vals_indx = zeros(1,num_signals); MUSIC_pk_num_signals_loc = zeros(1,num_signals);
% iterate through num_signals (num_signals) to find the maximum num_signals values from peakfinder
for i = 1:num_signals
    pk_vals_indx(i) = find(MUSIC_pk_num_signals_vals(i)==MUSIC_pk);
    MUSIC_pk_num_signals_loc(i) = MUSIC_loc(pk_vals_indx(i));
end

%%% New/Revised Peak Finder End %%%

MUSIC_DOAs = theta(MUSIC_pk_num_signals_loc); % for -80° produces at -79,0,60 results
for i = 1:Ns
    xline_var=xline(MUSIC_DOAs(i),':k',["\theta = "+string(MUSIC_DOAs(1))+"^\circ"],LineWidth=2); hold on;
    xline_var.FontSize = 24; 
end
disp('Determined MUSIC Doas: '+join(string(MUSIC_DOAs)));
hold off;

% axis equal % unsupported for multiple coordinate system axes 
% box on
setTightMargins(0.05);

% Define Path for Images to be sent to (so figures are not being saved to current directory) 
% This is used when no perturbations but DOA > |60°| since can cause
% interesting occurrences to MUSIC subspectra: Add on the signals [ts] to saved figure name
if perturbations == 0
    saveas(MUSIC_fig,'MUSIC_NO_Perturbations'+join(string(ts))+"NumSnapshotsK"+string(K)+'.png');
else
    % Used normally when array position imperfections are applied 
    saveas(MUSIC_fig,'MUSIC_ArrayPerturbations_Source'+join(string(ts))+...
        "PerturbMult"+string(average_perturbation_multiplier)+"Lambda"+"NumSnapshotsK"+string(K)+'.png');
end
