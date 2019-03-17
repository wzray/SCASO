% close all;
% clear
% 
% %% test simulation 
% OFDM_TX; % run base transmitter which you can use for developing. 
% 
% [decoded_data]= MyOfdmReceiver(raw_rx_data,MOD_ORDER);
%  
% % rx_data is the final output corresponding to tx_data, which can be used
% % to calculate BER
% [number,ber_16QAM_1] = biterr(tx_data,decoded_data);
% 
% ber_16QAM_1


%% test real world dat 
% load("tx_data_all.mat"); 
[samples n] = load_samples('Received data/otfs_rx.dat','float32',0);
[decoded_data]= MyOTFSReceiver(samples(1.5e6:2e6).',64,64,16,50,otfs_pre_copy,tx_data_all);
[number,ber_16QAM_2] = biterr(tx_data_all,decoded_data);
ber_16QAM_2

