function [A,recA]=dec_rec(I, N)
im=double(I);

%---------- This is main part -----------------%
wname = 'rbio4.4';
%N=2;
[C,S] = wavedec2(im,N,wname); %% N is the number of decomposition
A = appcoef2(C,S,wname)/2^N; 
na = prod(S(1,:));
recA = waverec2(C.*[ones(1,na) zeros(1,length(C)-na)],S,wname); 