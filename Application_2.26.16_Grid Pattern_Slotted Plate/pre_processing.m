function [I1,I2] = pre_processing(Im1,Im2,level_wavelet,size_filter)


Im1=double(Im1);
Im2=double(Im2);

I1=Im1;
I2=Im2;


% decoposition based on wavelet transform
N0=level_wavelet;
%N0=1;
N1=0;
N2=0;
[dI0_ds_1,dI1_ds_1,dI2_ds_1,dI0_1,dI1_1,dI2_1]=DWT_Decomp3(I1, N0, N1, N2);
[dI0_ds_2,dI1_ds_2,dI2_ds_2,dI0_2,dI1_2,dI2_2]=DWT_Decomp3(I2, N0, N1, N2);

I1=dI0_ds_1;
I2=dI0_ds_2;


% applying a Gaussian filter to images
mask_size=size_filter;
%mask_size=4;
std=mask_size*0.62;
H1=fspecial('gaussian',mask_size,std);
H2=fspecial('gaussian',mask_size,std);
I1=(imfilter(I1,H1)+imfilter(I1,H2))/2;
I2=(imfilter(I2,H1)+imfilter(I2,H2))/2;


% floor_value=5;
% for i=1:length(I1(:,1))
%     for j=1:length(I1(1,:))
%         if I1(i,j)<floor_value
%             I1(i,j)=10;
%         end
%     end
% end
% 
% for i=1:length(I2(:,1))
%     for j=1:length(I2(1,:))
%         if I2(i,j)<floor_value
%             I2(i,j)=10;
%         end
%     end
% end








