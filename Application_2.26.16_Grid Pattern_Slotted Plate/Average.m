clear
clc

%BROWSE FOR DATA FILES
    disp('SELECT DATA FILES')
    [filename, pathname, ~] = uigetfile(...
        {'*.tif','OFI Images (*.tif)';'*.*','All Files (*.*)'},...
        'Pick a file','MultiSelect', 'on');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    if iscell(filename)==0
        NF=1;
    else
        NF=numel(filename);
    end

    
%CROP IMAGE
    clc
    CP=1;
    if NF==1
       file=filename; 
    else
       file=filename{1}; 
    end
    I=imread(file);    
    while CP==1 
        disp('SELECT REGION OF INTEREST AND DOUBLE CLICK')
        I2=figure; imshow(I,[min(I(:)) max(I(:))]);
        [I3,ROI] = imcrop(I2);
        close(I2)
        I4=figure; imshow(I3,[min(I3(:)) max(I3(:))]);
        pause(2)
        clc
        S=input('SELECT A NEW REGION OF INTEREST?\n','s');
        SS=sum(double(S));
        if SS==78||SS==110||SS==157||SS==189||SS==221 %no
            CP=0;
        elseif SS==89||SS==121||SS==241||SS==273||SS==305||SS==337 %yes
            CP=1;
        end
        close(I4)
    end    
    
    
for j=1:1:NF
    if j==1
        Isum=imread(filename{1});
    else
        Isig=imread(filename{j});
        Isum=Isum+Isig;
    end
end
Isum=Isum./NF;
II=imcrop(Isum,ROI);
mean(II(:))
figure
imshow(II,[min(II(:)) max(II(:))])
% file=sprintf('Average_10.tif');
% imwrite(Isum,file)




