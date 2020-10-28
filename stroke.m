%% stroke.m
% input strokes on lit are and shadow area from the user

close all;
ar=[88];
impath = 'imgs\original\';% shadow image
path = 'imgs\stroke\'; % output path

len = length(ar);
wz = 6;

for i = 1:len
    im = imread([impath,num2str(ar(i)),'.png']);% shadow image
    imsz = size(im); imhw = imsz(1:2);
    
    imshow(im);           % display image
    msk = zeros(imhw);
    title('Select Lit Pixels');
    while true % get lit pixels
        [~,x,y] = freehanddraw(gca,'color','r','linewidth',wz); % get xy
        if length(x)>2
            tmsk = xy2msk([x,y]',imhw,wz); % convert xy to mask
            msk = msk + tmsk;
        else break;
        end
    end
    title('Select Shadow Pixels');
    while true % get shadow pixels
        [~,x,y] = freehanddraw(gca,'color','y','linewidth',wz); % get xy
        if length(x)>2
            tmsk = xy2msk([x,y]',imhw,wz); % convert xy to mask
            msk = msk + 0.5*tmsk;
        else break;
        end
    end
    imwrite(msk,[path,num2str(ar(i)),'.jpg']);
end
close gcf;
datapath(false);
