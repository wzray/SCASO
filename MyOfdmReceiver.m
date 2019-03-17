
function [decoded_data]= MyOfdmReceiver(data);
 %% run transmitter code to load sts and lts and other parameters 
%  OFDM_TX;
% OTFS_TX
N_OFDM_SYMS             = 500;         % Number of OFDM symbols
MOD_ORDER               =  1;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])
pilots = [1 1 -1 1].';
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
trel = poly2trellis(7, [171 133]);              % Define trellis
channel_coding = .5; % coding rate 
trellis_end_length = 8; % bits for trellis to end 
%% Rx processing params

rx_data = data;          % run OFDM tx code to get raw_rx_dec
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation

% Repeat the following code for each packet

%% Packet Detection
% auto correlation 16
length_samples= length(rx_data) - 200;
a_sample=16;
while( a_sample < length_samples)
    
    output_a(a_sample)= rx_data(a_sample-length(sts_t) + (1:length(sts_t))) * rx_data(a_sample + (1:length(sts_t)))' ./norm(rx_data(a_sample+(1:length(sts_t)))); 
    
    a_sample= a_sample+1;
       
end
output_a= output_a./max(abs(output_a));
% a=1;
% while(abs(output_a(a))<LTS_CORR_THRESH*max(abs(output_a)) ||isnan(output_a(a)))
%     a=a+1;
% end
index=find(abs(output_a)>0.9);
start_index_a=index(1)-15;
start_index = start_index_a;

%% lts cross correlation
lts_sample=64;
length_samples=length(rx_data) - 20000;
while( lts_sample < length_samples)
    
    output_lts(lts_sample)= rx_data(lts_sample-length(lts_t) + (1:length(lts_t))) * lts_t' ./norm(rx_data(lts_sample-length(lts_t) + (1:length(lts_t))));
    
    lts_sample= lts_sample+1;
end
output_lts= output_lts./max(abs(output_lts));
%a=1;
% while(abs(output_lts(a))<LTS_CORR_THRESH*max(abs(output_lts)) ||isnan(output_lts(a)))
%     a=a+1;
% end
a=find(output_lts>LTS_CORR_THRESH);
start_index_lts=a-63-512;
start_index=start_index_lts;


% disp("start index_lts is");
% disp(start_index_lts);
% disp("start index_c is ");
% disp(start_index_c);
disp("start index_a is ");
disp(start_index_a);
figure(1);
title("stf auto");
plot((abs(output_a))); 
% % hold on 
% % plot(abs(rx_data)); 
% % hold off
% figure(2);
% title("stf cross");
% plot((abs(output_a)));
% figure(3);
% plot((abs(output_lts)));
% title("ltf cross correlation");

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices
%% CFO estimation and correction
LTS_packet_1=rx_data(start_index+480+32:start_index+480+32+63);
LTS_packet_2=rx_data(start_index+480+32+64:start_index+480+32+64+63);
avg_diff=mean(angle(LTS_packet_2./LTS_packet_1));  
freq_offset=(avg_diff)/(2*pi*64);
disp("freq_offset");
disp(freq_offset);
ext_packet=rx_data(start_index:end);
cfo_corr_packet=ext_packet .* exp(-1i*2*pi*freq_offset*[0:length(ext_packet)-1]);
LTS_packet_1_corr=cfo_corr_packet(1+480+32:1+480+32+63);
LTS_packet_2_corr=cfo_corr_packet(1+480+32+64:1+480+32+64+63);
payload=cfo_corr_packet(1+640:640+500*80);

% Output: Packet with each value multiplied by CFO correction factor
%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help
payload_mat_cp=(vec2mat(payload,80)).';
payload_mat=payload_mat_cp(17:end,:);

% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)


% FFT
% Refer to IFFT perfomed at TX
payload_mat_f=fft(payload_mat,64);
% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix
channel1=(fft(LTS_packet_1_corr,64))./lts_f;
channel2=(fft(LTS_packet_2_corr,64))./lts_f;
channel_avg=(channel1+channel2)/2;
payload_eq=payload_mat_f(:,1:end)./(channel_avg.');
rx_syms1_mat=payload_eq(SC_IND_DATA,:);                                   
rx_syms1=rx_syms1_mat(:).';
% Output : Symbol equalized matrix in frequency domain of same size 
 % Advanced topics: 
%% SFO estimation and correction using pilots
payload_eq1=payload_eq;
iter_lim=2;
for iter=1:iter_lim
    pilot_data=payload_eq1(SC_IND_PILOTS,:);
    pilot_data=pilot_data.*pilots_mat;
    pilot_x_axis=[14 22 14];
    sfo_ph=(diff(unwrap(angle(pilot_data)),[],1))./(pilot_x_axis.');
    sfo_ph_avg=mean(sfo_ph);                                                  %find mean slope of sfo vs sub carrier index for each symbol
    for i=1:length(sfo_ph_avg)
        payload_no_sfo(:,i)=payload_eq(:,i).*exp(-1i*sfo_ph_avg(i)*(0:63)');  %correct sfo for each symbol
    end
    % Output: Symbol equalized matrix with pilot phase correction applied
    % Phase Error Correction using pilots
%     Extract the pilots and calculate per-symbol phase error
    pilot_data2=payload_no_sfo(SC_IND_PILOTS,:);
    pilot_data2=pilot_data2.*pilots_mat;
    ph_noise=angle(mean(pilot_data2));                                            % find remaining phase error for each symbol.
    payload_no_ph_noise=payload_no_sfo.*(exp(-1i*(ph_noise)).*ones(N_SC,1));
    payload_eq1=payload_no_ph_noise;
end
payload_data=payload_no_ph_noise(SC_IND_DATA, :);
rx_syms=payload_data(:).';
% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms 

% Demodulation

figure(4);
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

figure(5);
scatter(real(payload_data(48,:)), imag(payload_data(48,:)),'filled');
title(' Signal Space of received bits subcarr 48');
xlabel('I'); ylabel('Q');

figure(6);
scatter(real(payload_data(30,:)), imag(payload_data(30,:)),'filled');
title(' Signal Space of received bits subcarr 30');
xlabel('I'); ylabel('Q');

figure(7);
scatter(real(rx_syms1), imag(rx_syms1),'filled');
title(' Signal Space of received bits before sfo');
xlabel('I'); ylabel('Q');

% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
decoded_data = vitdec(Demap_out,trel,7,'trunc','hard');

% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER
