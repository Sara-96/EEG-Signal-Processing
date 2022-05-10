X;  % 32-channel data
load('Electrodes1') ;

% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = Electrodes1.labels ;
disp_eeg(X,offset,feq,ElecName);
