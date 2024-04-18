function [out] = Point(in)

%Calibration Coeficients:
    %Tfit = [-67.0043 6.9917 79.3194];
    %Tfit = [-19.5070 3.7725 79.3194];
    Tfit = [-64.3789 6.8534 79.3194]; % ref = 22.5 C
    ROI = [1 1 1608 1208];
    h = fspecial('average',[9 9]);
    a = double(imread(in{1}));
    b = imcrop(a,ROI);
    bg = imfilter(b,h);
    %b = imcrop(a,ROI);
    %sig(1) = mean2(b);

    a = double(imread(in{2}));
    b = imcrop(a,ROI);
    off = imfilter(b,h);
    a = double(imread(in{3}));
    b = imcrop(a,ROI);
    on = imfilter(b,h);

    figure(11)
    imagesc(bg);
    impixelinfo;

    figure(12)
    imagesc(off);
    impixelinfo;

    figure(13)
    imagesc(on);
    impixelinfo;
    pause(1)

    a = (on-bg)./(off-bg);
    b = imfilter(a,h);

    out = polyval(Tfit,b);
    %out = b;

    figure(4)
    imagesc(out,[30 80]);
    impixelinfo;
    colorbar
    impixelinfo
    pause(0.5)
