% Figure outputs 2-source MUSIC and Eigenbeamformer plots
clear all; close all; clc;

% setup file for Matlab figures
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',18)
legend_font_size = 12;
%% Problem Creation and Eigenvalues/Eigenvectors created from Correlation Matrix
randseed=11; rand('seed',randseed);randn('state',0)
k=2*pi;
cs=2;

% with N = 5 produces peaks also at -3(N=5) and 63(N=4)
N=4; d=0.5;   % number of array elements
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

Rss=s*s'/K;                           % signal correlation matrix
Rnn=(nois*nois')/K;                   %noise correlation matrix
Rns=(nois*s')/K;                      % noise/signal correlation matrix
Rsn=(s*nois')/K;                      % signal/noise correlation matrix
Rrr=A*Rss*A'+A*Rsn+Rns*A'+Rnn;        % array correlation matrix

[eigvec, eigval] = eig(Rrr);

%% Plotting Values to create array steering vector for all angles and more angular resolution

% Find the output where A is multiplied by w and plot over theta
L = 181; % # of theta points to take 1 degree for this case --> Samples
theta = linspace(-90,90,L);
theta_rad = theta*pi/180;
% A_new allows plotting by splitting the Array Manifold Vector based on all
% angles not just maximum points 
A_new = exp(1j*k*x'*sin(theta_rad)); % 1 x L 

%% Eigenbeam Plotting 
% Produce the weights for each set of signals by taking 1 column from the
% eigenvectors where the 2 largest are the signal weights
% Pre-allocate Arrray Sizes
w = ones(N);
y = ones(N,L);
y_dB = y;
figure(1)
eigen_values_plot = diag(real(eigval))';
% pre-allocate array - could also use strings
subtitle_eigval = ["1" "1" "1" "1"]; legend_vals = []; 

for i = 1:N   
    subtitle_eigval(i) = string(i) + ": " + string(eigen_values_plot(i));
    w(:,i) = eigvec(:,i);
    y(i,:) = w(:,i)'*A_new;
    y_dB(i,:) = 20*log10(abs(y(i,:)));
    plot(theta,y_dB(i,:)); hold on;
    legend_vals = [legend_vals "EigenBeam "+string(i)];
end

grid on; 
% Since we know 2 signals 2 noise eigenbeams: 
P = N - length(sa); % number of noise eigenbeams
[smallest_val, indx] = mink(eigen_values_plot,P)

% use the two smallest eigenvalues to plot line at their minimum point
for i = 1:P
    [peak_vals(P,:),locations(P,:)] = findpeaks(-y_dB(indx(P),:)) % find peaks for - of data for minimum values
    % Plot Minimum Values for each Noise EigenBeam
    for k = 1:length(locations(P,:))
        theta_min = theta(locations(P,k)); 
        xline(theta_min,':k',["\theta = "+string(theta_min)+"^\circ"],'LineWidth',2); % second argument is line style
        hold on;
    end
end
hold off; 

% Determine best orientation and Location  of legend
legend(legend_vals,Location="best",FontSize=legend_font_size,Orientation="vertical");
title('Eigenbeam Plots');

% Max Eigenbeam at desired DOAs is the desired eigenvector = use as weights 
subtitle_txt = ["Eigenbeams: "+join(subtitle_eigval(1)+" | "+subtitle_eigval(2)),...
                "Eigenbeams: "+join(subtitle_eigval(3)+" | "+subtitle_eigval(4))];

subtitle(subtitle_txt);
xlabel("\theta^\circ")
ylabel("AF(dB)")
setTightMargins(0.05)
saveas(gcf,'Eigenbeamformer.png')

%% MUSIC Algorithm DoA
% Take minimum columns from eigval to produce noise_eigvectors
% Idea is signal eigenvalues are a lot larger in correlation values than
% noise eigenvalues
% MUSIC fails at low SNR

% There are P = N-M for N = # array elements, M = # signals noise
% eigenvalues/vectors
% Pick from the minimum of these to create the MUSIC spectrum for best
% results
noise_eigval = mink(diag(eigval),P); % cs = # of signals
num_DOAs = cs;

% Don't need to iterate, just use minimum eigenvalue of the P noise
% eigenvalues and use corresponding single minimum eigenvector to create 
% MUSIC spectrum
min_eigenval = min(noise_eigval);
% find noise eigenvalue column position in original eigenvalue matrix
[~,min_eigenval_idx] = find(min_eigenval==eigval);
% using known column of minimum (noise) eigenvalue determine noise eigenvector
noise_eigvec = eigvec(:,min_eigenval_idx);

% 2 sub-spectra for each noise eigenvalue from the overall eigenvalues
P_theta = 1./abs(A_new'*noise_eigvec).^2;
P_dB = 10*log10(abs(P_theta));

figure(10)
%plot(theta,P_theta); %hold on; yyaxis right;
plot(theta',P_dB); hold on; % plot peaks
title("MUSIC Algorithm")
xlabel("\theta^\circ");
%ylabel("P(\theta)"); yyaxis right;
ylabel("P(\theta) [dB]")

% Find Peaks / DOA on MUSIC Sub-Spectra - Annotate ignored if called with
% output arguments
% For this case also finds peak at 0° with value 5.8 juse use known maximum 2 values 
% Use maxk to only grab the cs (2) largest signals
[MUSIC_pk,MUSIC_loc] = findpeaks(P_dB);
MUSIC_pk_cs_vals = maxk(MUSIC_pk,cs);
% Pre-allocate:
pk_vals_indx = zeros(1,cs); MUSIC_pk_cs_loc = zeros(1,cs);
% iterate through cs to find the maximum cs values from peakfinder
for i = 1:cs
    pk_vals_indx(i) = find(MUSIC_pk_cs_vals(i)==MUSIC_pk);
    MUSIC_pk_cs_loc(i) = MUSIC_loc(pk_vals_indx(i));
end

MUSIC_DOAs = theta(MUSIC_pk_cs_loc);
xline(MUSIC_DOAs(1),':k',["\theta = "+string(MUSIC_DOAs(1))+"^\circ"],'LineWidth',2); hold on;
xline(MUSIC_DOAs(2),':k',["\theta = "+string(MUSIC_DOAs(2))+"^\circ"],'LineWidth',2); 
setTightMargins(0.05);
saveas(gcf,'MUSIC2Sources.png');