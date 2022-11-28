% =========================================================================
% --                    Passband Digital Transmission                    --
% --                 QPSK Modulator - Coherent detector                  --
% =========================================================================

clear; clc;
load('message')

% Simulation parameters
b=2; M=2^b;                                             % Number of bits per symbol and the modulation order
symbol = [1 0; 0 0; 0 1; 1 1];                          % Set of symbols
Tb=1; Tsym=b*Tb;                                        % Bit duration and Symbol duration
Nb=32; Ns=b*Nb;                                         % Number of sample times in Tb and Tsym
Ts=Tsym/Ns;                                             % Sample time
fc=4/Tsym;                                              % Carrier frequency
t=(0:Ns-1)*Ts;                                          % Time slot
Es=2; sqEs=sqrt(Es);                                    % Energy of signal waveform
SNRdBs=1:10; MaxIter=10000;                             % Range of SNRbdB and number of iterations

% Plot parameters (buffer)
% Note that we are not going to plot every single message. Rather, we are going to plot the very last four simulations
LB=4*Ns;                                                % buffer size
tt=(0:LB-1)*Ts;
signal_waveforms=zeros(1,LB); 
received_signals=zeros(1,LB);
correlator_outputs=zeros(2,LB);

% Generate QPSK signal waveforms
% ===============================To Do=====================================

phase = [pi/4; 3*pi/4; 5*pi/4; -pi/4];

% =========================================================================

% Implement Tx-channel-Rx
% Oscillator
su = sqrt(2/Tsym)*[cos(2*pi*fc*t); -sin(2*pi*fc*t)]; suT=su*Ts;

for iter=1:length(SNRdBs)                               
    SNRbdB = SNRdBs(iter); SNR = 10^(SNRbdB/10);
    N0=2*(Es/b)/SNR; sgmT=sqrt(N0/2/Ts);
    errs=0;                                           
    for k=1:MaxIter
        i=message(k); s=symbol(i,:);
        % ===========================To Do=================================
        % Implement Transmitter, Channel and Correlator
        yr=zeros(2,Ns);
        for n=1:Ns       
            %%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transmitter
            x = (2 * s(1) - 1) * su(1, n) - (2 * s(2) - 1) * su(2, n);
            
            %%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%
            % Bandpass noise (Channel) and noisy received signal
            noise = sgmT * randn;

            %%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%
            r = x + noise;

            %%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%
            % Correlator
            yr_t = [r * suT(1, n); r * suT(2, n)];
            yr = [yr(:, 2:Ns) yr_t];

            % Save last four Tx & Rx signals in buffer (Don't change)
            if iter==length(SNRdBs) && k > MaxIter-5
                signal_waveforms=[signal_waveforms(2:LB) x];
                received_signals=[received_signals(2:LB) r];                 
                correlator_outputs=[correlator_outputs(:,2:LB) yr(:,end)];
            end
        end
        % Implement Coherent Detector (Don't change the variable name d)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        sum_yr = sum(yr, 2);

        arctan = atan(sum_yr(2) / sum_yr(1));
        if sum_yr(1) < 0
            arctan = arctan + pi;
        end
        diffs = abs(phase - arctan);
        argmin = find(diffs == min(diffs));
        
        d = symbol(argmin, :);
        % ================================================================
        errs=errs+sum(s~=d);                         
        if errs>100, break; end
    end
    BER(iter)=errs/(k*b);
end

% Plot
SNRbdBt=0:0.1:10; SNRbt=10.^(SNRbdBt/10);
BER_theoretical=(1+(b>1))/b*qfunc(sqrt(b*SNRbt)*sin(pi/M));       
subplot(221), plot(tt,signal_waveforms,'k',tt,received_signals,'b:'), title('Received signal')
subplot(222), plot(tt,correlator_outputs(1,:),'b:'), title('In-phase correlator output')
subplot(223), plot(tt,correlator_outputs(2,:),'b:'), title('Quadrature correlator output')
subplot(224), semilogy(SNRbdBt,BER_theoretical,'k-',SNRdBs,BER,'b*'), title('BER for QPSK Signaling')