function [dI0_ds,dI1_ds,dI2_ds,dI0,dI1,dI2]=DWT_Decomp3(I, N0, N1, N2)

I=double(I);
[dI0_ds,recI0]=dec_rec(I, N0);
dI0=I;

dI1=I-recI0;
dI1=255*(dI1-min(min(dI1)))/(max(max(dI1-min(min(dI1)))));
[dI1_ds,recI1]=dec_rec(dI1, N1);

dI2=dI1-recI1;
dI2=255*(dI2-min(min(dI2)))/(max(max(dI2-min(min(dI2)))));
[dI2_ds,recI2]=dec_rec(dI2, N2);




