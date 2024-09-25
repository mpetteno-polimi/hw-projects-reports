clear all
close all
clc

%% Import Input Audio Signal
[Vin, ~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice, ~] = audioread('outlowsweep.wav');
[OutMidSpice, ~] = audioread('outmidsweep.wav');
[OutHighSpice, FsLTSpice] = audioread('outhighsweep.wav');
TsLTSpice = 1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact = 4;
Fs = FsLTSpice / downSampFact; 

%% Downsample Input Signal
Vin = Vin([1:downSampFact:end]);

%% Sampling Period
Ts = 1/Fs;

%% Number of Samples
Nsamp = length(Vin);
%% Simulated time
tstop = Nsamp*Ts;
%% Parameters of Dynamic Element
L1 = 0.35*10^(-3);
L2 = 0.35*10^(-3);
L3 = 3.5*10^(-3);
L4 = 3.5*10^(-3);
C1 = 2.8*10^(-6);
C2 = 2.8*10^(-6);
C3 = 28*10^(-6);
C4 = 4.7*10^(-6);
C5 = 28*10^(-6);
C6 = 47*10^(-6);
%% Resistive Parameters
R1 = 10;
RspkLow = 8;
R2 = 10;
RspkMid = 8;
RspkHigh = 8;
%% WDF setting of free parameters (adaptation conditions)
% HPF
syms Z4_hpf real;
Z1_hpf = RspkHigh;
Z2_hpf = (2*L1)/Ts;
Z3_hpf = Ts/(2*C1);
% BPF
syms Z8_bpf real;
Z1_bpf = RspkMid;
Z2_bpf = Ts/(2*C4);
Z3_bpf = (2*L3)/Ts;
Z4_bpf = Ts/(2*C2);
Z5_bpf = R1;
Z6_bpf = Ts/(2*C3);
Z7_bpf = (2*L2)/Ts;
% LPF
syms Z6_lpf real;
Z1_lpf = RspkLow;
Z2_lpf = Ts/(2*C6);
Z3_lpf = Ts/(2*C5);
Z4_lpf = R2;
Z5_lpf = (2*L4)/Ts;

%% Computation of Scattering Matrices
% HPF
N_port_hpf = 4;
N_ind_curr_hpf = 2;
Z_hpf = diag([Z1_hpf, Z2_hpf, Z3_hpf, Z4_hpf]);
B_hpf = [eye(N_ind_curr_hpf); [1, 1; -1, -1]]';
S_hpf = eye(N_port_hpf) - 2*Z_hpf*B_hpf'*inv(B_hpf*Z_hpf*B_hpf')*B_hpf;
S_hpf = subs(S_hpf, Z4_hpf, solve(S_hpf(4,4)==0, Z4_hpf)); % adapt junction to non linear element port
S_hpf = double(S_hpf); % convert back to double

% BPF
N_port_bpf = 8;
N_ind_curr_bpf = 4;
Z_bpf = diag([Z1_bpf, Z2_bpf, Z3_bpf, Z4_bpf, Z5_bpf, Z6_bpf, Z7_bpf, Z8_bpf]);
B_bpf = [eye(N_ind_curr_bpf); [0, 1, 0, 0; 1, 1, 1, 0; 1, 1, 1, 1; -1, -1, -1, -1]]';
S_bpf = eye(N_port_bpf) - 2*Z_bpf*B_bpf'*inv(B_bpf*Z_bpf*B_bpf')*B_bpf;
S_bpf = subs(S_bpf, Z8_bpf, solve(S_bpf(8,8)==0, Z8_bpf)); % adapt junction to non linear element port
S_bpf = double(S_bpf); % convert back to double

% LPF
N_port_lpf = 6;
N_ind_curr_lpf = 3;
Z_lpf = diag([Z1_lpf, Z2_lpf, Z3_lpf, Z4_lpf, Z5_lpf, Z6_lpf]);
B_lpf = [eye(N_ind_curr_lpf); [0, 1, 0; 1, 1, 1; -1, -1, -1]]';
S_lpf = eye(N_port_lpf) - 2*Z_lpf*B_lpf'*inv(B_lpf*Z_lpf*B_lpf')*B_lpf;
S_lpf = subs(S_lpf, Z6_lpf, solve(S_lpf(6,6)==0, Z6_lpf)); % adapt junction to non linear element port
S_lpf = double(S_lpf); % convert back to double

%% Initialization of Waves
% (a) waves incident to the junction, that are reflected waves by the elements
% (b) waves reflected by the junction, that are incident waves to the elements
% HPF
a_hpf = zeros(N_port_hpf, 1);
b_hpf = zeros(N_port_hpf, 1);
% BPF
a_bpf = zeros(N_port_bpf, 1);
b_bpf = zeros(N_port_bpf, 1);
% LPF
a_lpf = zeros(N_port_lpf, 1);
b_lpf = zeros(N_port_lpf, 1);

%% Initialize Output Signals
% Low
VoutLow = zeros(size(Vin));
% Mid
VoutMid = zeros(size(Vin));
% High
VoutHigh = zeros(size(Vin));

ii = 0;
while (ii < Nsamp)
    ii = ii + 1;

    %% Manage Dynamic Elements
    % For dynamic elements, current reflected wave (a) depends on previous
    % incident wave (b), then we use (a) to compute the current reflected wave
    % HPF
    a_hpf(2) = -b_hpf(2); % L1
    a_hpf(3) = b_hpf(3); % C1
    % BPF
    a_bpf(2) = b_bpf(2); % C4
    a_bpf(3) = -b_bpf(3); % L3
    a_bpf(4) = b_bpf(4); % C2
    a_bpf(6) = b_bpf(6); % C3
    a_bpf(7) = -b_bpf(7); % L2
    % LPF
    a_lpf(2) = b_lpf(2); % C6
    a_lpf(3) = b_lpf(3); % C5
    a_lpf(5) = -b_lpf(5); % L4
    
    %% Forward Scan
    % Compute waves reflected by the adaptor towards the nonlinear root element
    b_hpf(4) = S_hpf(4,:)*a_hpf;
    b_bpf(8) = S_bpf(8,:)*a_bpf;
    b_lpf(6) = S_lpf(6,:)*a_lpf;

    %% Local Root Scattering
    % Compute wave reflected by the nonlinear root element
    a_hpf(4) = 2*Vin(ii) - b_hpf(4);
    a_bpf(8) = 2*Vin(ii) - b_bpf(8);
    a_lpf(6) = 2*Vin(ii) - b_lpf(6);

    %% Backward Scan
    % Compute waves reflected by the adaptor towards the linear elements
    % (for brevity also recomputes the waves reflected towards the 
    % nonlinear root element as in Forward Scan)
    b_hpf = S_hpf*a_hpf;
    b_bpf = S_bpf*a_bpf;
    b_lpf = S_lpf*a_lpf;

    %% Read Output
    % We can omit a(1) cause reflected waves from resistor are always 0
    VoutHigh(ii) = b_hpf(1)/2;
    VoutMid(ii) = b_bpf(1)/2;
    VoutLow(ii) = b_lpf(1)/2;
    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

