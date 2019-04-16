clear; close all; clc;
addpath(genpath('.'));

%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length
time_axis = (1:N)/Fs;

%% Low Pass
B_low = conv([1 zeros(1,5) -1], [1 zeros(1,5) -1]);
A_low = conv([1 -1], [1 -1]);
%fvtool(B_low,A_low);

y1 = filter(B_low, A_low, data);


%% High Pass
B_high = [-1 zeros(1,15) 32 -32 zeros(1,14) 1];
A_high = [1 -1];

%fvtool(B_high,A_high);

y2 = filter(B_high, A_high, y1);
%% Band Pass
%fvtool(conv(B_low,B_high),conv(A_low,A_high));

y_filtered = filter(conv(B_low,B_high), conv(A_low,A_high), data);

%% Differentitation

% Matlab can't implement non causal filters
% one Solution is to shift our sample two times in time
% We can only implement H(Z)/Z²
% the results comming out of this filer are shifter +2

% we need to define Ts = 1/Fs


B_diff = [1 2 0 -2 -1];
A_diff = [8/Fs];

% this is H(Z)/Z²
% fvtool(B_diff,A_diff);

y_derived = filter(B_diff,A_diff, y_filtered);
plot((1:N)/Fs,y_filtered);

%% Abs²

y_squared = abs(y_derived).^2;

%% Moving Window Integration
M = 5000;
Ts = 0.08;
s_MWI = zeros(1,N - M + 1);

vect = M:N-M+1;

for n=M:N-M+1
    s_MWI(n-M+1) = sum(y_squared(1, n-M+1:n));
end


%%

k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
k(5<k & k<10)



