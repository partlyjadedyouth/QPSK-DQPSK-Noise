% =========================================================================
% --                    Passband Digital Transmission                    --
% --                DQPSK Modulator - Non-Coherent detector              --
% =========================================================================

clear; clc;
load('message')

% Simulation parameters
b=2; M=2^b;                                             % Number of bits per symbol and the modulation order
symbol = [1 0; 0 0; 0 1; 1 1];                          % Set of symbols
Tb=1; Tsym=b*Tb;                                        % Bit duration and Symbol duration
Nb=32; Ns=b*Nb;                                         % Number of sample times in Tb and Ts
Ts=Tsym/Ns;                                             % Sample time
fc=4/Tsym;                                              % Carrier frequency
t=(0:Ns-1)*Ts;                                          % Time slot
Es=2; sqEs=sqrt(Es);                                    % Energy of signal waveform
SNRdBs=1:10; MaxIter=10000;                             % Range of SNRbdB and number of iterations

% Plot parameters
% Note that we are not going to plot every single message. Rather, we are going to plot the very last four simulations
LB=4*Ns;                                                % buffer size
tt=(0:LB-1)*Ts;
signal_waveforms=zeros(1,LB); 
received_signals=zeros(1,LB);
correlator_outputs=zeros(2,LB);

% Generate DQPSK Signal Waveforms
% ===============================To Do=====================================
 
phase = [0; pi/2; pi; 3*pi/2];

% =========================================================================

% Implement Tx-channel-Rx
% oscillator
% Modified su b.c. the constellation was changed
su=sqrt(2/Tsym)*[cos(2*pi*fc*t - pi/4); -sin(2*pi*fc*t - pi/4)]; suT=su*Ts;

for iter=1:length(SNRdBs)
    SNRbdB=SNRdBs(iter); SNR=10^(SNRbdB/10);
    N0=2*(Es/b)/SNR; sgmT=sqrt(N0/2/Ts);
    errs=0;
    is0=1;                                              % Initial signal
    th0=0;                                              % Initial guess of signal phase, which is possibly wrong
    for k=1:MaxIter
        i=message(k); s=symbol(i,:);
        % ===========================To Do=================================
        % Implement Transmitter, Channel and Correlator
        yr=zeros(2,Ns);

        %%%%%%%%%%%%%%%%%% To Do : Differnetial encoding %%%%%%%%%%%%%%%%%%
        i_prime = mod(is0 + i - 1, 4) + 1;
        symbol_prime = symbol(i_prime, :);
        for n=1:Ns
            %%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transmitter
            x = (2 * symbol_prime(1) - 1) * su(1, n) - (2 * symbol_prime(2) - 1) * su(2, n);

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
        
        is0 = i_prime - 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Implement Non-Coherent Detector
        sum_yr = sum(yr, 2);

        th1 = atan(sum_yr(2) / sum_yr(1)) - pi/4;
        if sum_yr(1) < 0
            th1 = th1 + pi;
        elseif sum_yr(2) < 0 % (1, 1) -> 3*pi/2  
            th1 = th1 + 2*pi;
        end
        
        % calibrate 'phase difference < 0' case
        delta_th = th1 - th0;
        if delta_th < -pi/6
            delta_th = delta_th + 2*pi;
        end

        diffs = abs(phase - delta_th);
        argmin = find(diffs == min(diffs));
        d = symbol(argmin, :);
        
        th0 = th1;
        % ================================================================
        errs=errs+sum(s~=d); 
        if errs>100 && iter<10, break; end                                 
    end
    BER(iter)=errs/(k*b);
end

% Plot
SNRbdBt=0:0.1:10; SNRbt=10.^(SNRbdBt/10);
BER_theoretical=(1+(b>1))/b*qfunc(sqrt(b*SNRbt/2)*sin(pi/M)); 
subplot(221), plot(tt,signal_waveforms,'k',tt,received_signals,'b:'), title('Received signal')
subplot(222), plot(tt,correlator_outputs(1,:),'b:'), title('In-phase correlator output')
subplot(223), plot(tt,correlator_outputs(2,:),'b:'), title('Quadrature correlator output')
subplot(224), semilogy(SNRbdBt,BER_theoretical,'k-',SNRdBs,BER,'b*'), title('BER for DQPSK Signaling')