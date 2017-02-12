% Ahmed Elgabaly %
% 861054818 %
% EE 231 %
% Optimal Downlink Beam Forming using Semidefinite Optimization %

clc; close all; clear all; % clear all variables in the command window, workspace and close all plots

% Variables
Delta = 5:25; % Variation in the angle 
Theta_1 = deg2rad(10); % Angle of first mobile device
Theta_2 = deg2rad(10 + Delta); % Angle of second mobile device
Theta_3 = deg2rad(10 - Delta); % Angle of Third mobile device
SpreadAngle = deg2rad(2); % Spread angle
t = length(Delta); % Length of Delta

R1 = zeros(3,3);% Channel Covariance Matrix for the first mobile device
R2 = zeros(3,3,t); % t dimensional Channel Covariance Matrix for the second mobile device
R3 = zeros(3,3,t); % t dimensional Channel Covariance Matrix for the third mobile device

% R1 Calculation 
for k = 1:3
    for l = 1:3
        R1(k,l) = exp(1i*pi*(k-l)*sin(Theta_1)) *...
            exp( ((pi * (k-l) * SpreadAngle * cos(Theta_1))^2)/(-2) );
    end
end

% R2 Calculation
for c = 1:t
    for k = 1:3
        for l = 1:3
            R2(k,l,c) = exp(1i*pi*(k-l)*sin(Theta_2(c))) *...
            exp( ((pi * (k-l) * SpreadAngle * cos(Theta_2(c)))^2)/(-2) );
        end
    end
end

% R3 Calculation
for c = 1:t
    for k = 1:3
        for l = 1:3
            R3(k,l,c) = exp(1i*pi*(k-l)*sin(Theta_3(c))) *...
            exp( ((pi * (k-l) * SpreadAngle * cos(Theta_3(c)))^2)/(-2) );
        end
    end
end

% SemiDefinite Program %

SINR = 2; % Signal to interference + noise ratio
VARi = 1; % Variance of the noise
Transmitted_Power = zeros(1,t); % Minimum Transmitted power vecotr for all three mobile devices at their current angle configuration

for i = 1:t % loop over every channel covariance matrix for every user at a different angle
    cvx_begin sdp quiet % Define a CVX semidefinite program 
    
    variable W1(3,3) complex semidefinite % define W1 as a complex symmertric positive semidefinite matrix 
    variable W2(3,3) complex semidefinite % define W2 as a complex symmertric positive semidefinite matrix
    variable W3(3,3) complex semidefinite % define W3 as a complex symmertric positive semidefinite matrix
    % Which automatically means that     
    % W1 == conjugate transpose(W1), W2 == conjugate transpose(W2), W3 == conjugate transpose(W3)
    % W1 >= 0, W2 >= 0, W3 >= 0
    minimize ( trace(W1) + trace(W2) + trace(W3) ) % Define the objective function 
    
    subject to % Define the Equality constraints for every mobile device
    trace(R1 * W1) - ( SINR * (trace(R1 * W2) + trace(R1 * W3)) ) == SINR * VARi; % Equality constraint on mobile device one
    trace(R2(:,:,i) * W2) - ( SINR * (trace(R2(:,:,i) * W1) + trace(R2(:,:,i) * W3)) ) == SINR * VARi; % Equality constraint on mobile device two
    trace(R3(:,:,i) * W3) - ( SINR * (trace(R3(:,:,i) * W1) + trace(R3(:,:,i) * W2)) ) == SINR * VARi; % Equality constraint on mobile device three
    
    cvx_end % End CVX program
    
    Transmitted_Power(i) = cvx_optval; % The optimal value is the minimum transmitted power from the base station 
    % to all three mobile devices simultaneously at their current angle
    % with respect to the base station
    cvx_clear; % clear all CVX data for the current model before starting a new one
end % End for loop

plot(Delta,Transmitted_Power); % Plot the angle variation of the mobile devices with respect to the base station VS
% the minimum transmitted power from the base station to all three mobile
% devices simultaneously at their current direction with respect to the
% base station 
xlabel('Theta'); % The X axis of the curve 
ylabel(['Transmitted Power SINR = ' num2str(SINR)]); % The Y axis of the curve