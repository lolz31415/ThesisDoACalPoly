% Simple script to compare the number of antenna elements to theoretical
% angular resolution: half-power beamwidth 
clc; clear all; close all; 
% setup file for Matlab figures
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',18)
legend_font_size = 12;

d = 0.5; %half wavelength spacing
k = 2*pi;
lambda = 1; % since d is relative to lambda
N = [2 4 5 7 10 12];
% Plot for N = 2,4,5,7,10,12
fig1 = figure(100);
fig1.Name = "VaryArrayElements";

%theta = linspace(-90,90,361);
theta = linspace(0,180,361);

% use sind for degrees
psi = k*d*cosd(theta);

for i = 1:length(N)
    current_N = N(i);
    % use N(i) for each value of N 
    % Take Absolute Value of E_theta? -> some parts are negative
    %E_theta = sin(N(i)/2*psi)./sin(1/2*psi);
     E_theta = 1/N(i)*abs( sin(N(i)/2*psi)./sin(1/2*psi) );
     % E_theta = 1/N(i)*abs( cos(N(i)/2*psi)./cos(1/2*psi) );

    E_theta_dB = 20*log10(E_theta);
    plot(theta,E_theta_dB); hold on; 
end
hold off; 
% ylim to enforce scaling to show scale of nulls for larger arrays 
ylim([-75,25]);
legend(string(N),Location="southeast",FontSize=legend_font_size,Orientation="vertical");
xlabel('\theta \circ');
% E-Field Pattern 
ylabel('Radiation Pattern (dB)');
% Previous Title didn't fit in chart/figure area
% title('Varying Number of Array Elements determine angular resolution ');
title('Array Azimuth E-Field Pattern');

% axis equal % unsupported for multiple coordinate system axes 
% box on
setTightMargins(0.05)
saveas(fig1,'VaryArrayElements.png')
