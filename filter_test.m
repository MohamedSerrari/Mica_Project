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

%% Abs²

y_squared = abs(y_derived).^2;

%plot(time_axis, y_squared);

%% Moving Window Integration
M = floor(Fs*0.09);
h = ones(1,M)/M;
s_MWI = filter(h,1,y_squared);

%% Thresholding
% Manual thresholding
findpeaks(s_MWI);
%th = 6*10^14;
%th = max(s_MWI)/sqrt(3);

%th = input("Enter the threshold:",'s')
answer = inputdlg('Enter the threshold: ','Defining the threshold');
th = str2num(answer{1});

%ind = zeros(1,length(s_MWI));
ind = s_MWI > th;

%% Finding Peaks

figure
hold on

delay = 8;
y_delayed = [zeros(1,delay) y_filtered];

QRS_only = y_delayed.*[ind zeros(1,delay)]; % This signal contains the QRS complexes only 

[pks, loc] = findpeaks(QRS_only(1,1:last)); 
pks(pks < 10^5) = 0; % Another thresholding is needed because findpeaks gives all the local minimas
loc_th = find(pks);
R_locs = loc(loc_th);

last = 1000;

data_plot = y_delayed(1,1:last); 
plot(data_plot);

for i=1:length(R_locs)
    % inter_pks is the distance between two peaks, it is updated everytime
    % for better precision in the code. In the case of the last PQRST
    % complex, inter_pks is not updated and the last value is kept.
    
    if (i < length(R_locs)) % checking is it's the last complex in the signal
        inter_pks = floor((R_locs(i+1) - R_locs(i))/2);
    end

    start = R_locs(i) - inter_pks;  % the interval [start:finish] has exactly one PQRST complex 
    finish = R_locs(i) + inter_pks; % This interval is centered arround R_locs(i)

    [Q_val, Q_loc] = min(data_plot(1,start:R_locs(i))); 
    [S_val, S_loc] = min(data_plot(R_locs(i):finish));

    Q_loc = Q_loc + start - 1;     % Correcting Q_loc 
    S_loc = S_loc + R_locs(i) - 1; % Correcting S_loc 

    [P_val, P_loc] = max(data_plot(1,start:Q_loc));
    [T_val, T_loc] = max(data_plot(1,S_loc: finish));

    P_loc = P_loc + start - 1; % Correcting P_loc 
    T_loc = T_loc + S_loc - 1; % Correcting T_loc 
    
    plot(R_locs(i),data_plot(R_locs(i)),'r*'); text(R_locs(i),data_plot(R_locs(i)),' R ','Color','red','FontSize',14);
    plot(Q_loc, Q_val,'r*'); text(Q_loc,Q_val,' Q ','Color','red','FontSize',14);
    plot(S_loc, S_val,'r*'); text(S_loc,S_val,' S ','Color','red','FontSize',14);
    plot(P_loc, P_val,'r*'); text(P_loc,P_val,' P ','Color','red','FontSize',14);
    plot(T_loc, T_val,'r*'); text(T_loc,T_val,' T ','Color','red','FontSize',14);
end

%% Spectrogram
ecg_normal_1 = load('ecg_normal_1.mat');
ecg_VF = load('ecg_VF.mat');

window = 4*60;
figure; spectrogram(ecg_normal_1.ecg, window);
figure; spectrogram(ecg_VF.ecg, window)

%% Tachycardia/Bradycardia

rythym = 1/sum(R_locs(2:end) - R_locs(1:end-1))*(length(R_locs) - 1)*60*Fs
