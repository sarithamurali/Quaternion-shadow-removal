close all;
clear all;
set(0,'DefaultFigureVisible','off')

ar=[88];
impath = 'imgs\original\';% shadow image
ipath = 'imgs\stroke\'; % user input path
outpath='imgs\removal\'; % removal output
dpath='imgs\detection\'; % shadow detection input

slen=length(ar);
for index = 1:slen
    clearvars -except slen index impath outpath ar ipath dpath
    
    %% input
    img = imread([impath,num2str(ar(index)),'.png']);
    imsk = im2double(imread([ipath,num2str(ar(index)),'.png'])); 
    
    if exist([dpath,num2str(ar(index)),'.png'])
    det=imread([dpath,num2str(ar(index)),'.png']);
    end
    
    I=rgb2ycbcr(img);
    doubI=im2double(I);
    
    %lit pixels mask
    limsk = imsk>0.9;
    
    %shadow pixels mask
    smsk = imsk<0.6&imsk>0.4;
    
    %selected lit pixels values
    [x(:,1),x(:,2)]=find(limsk~=0);
    for i=1:length(x)
        lit(i,1:3)=doubI(x(i,1),x(i,2),:);
    end
    %selected shadow pixels values
    [y(:,1),y(:,2)]=find(smsk~=0);
    for i=1:length(y)
        shad(i,1:3)=doubI(y(i,1),y(i,2),:);
    end
    
    mL=mean(lit,1);
    mS=mean(shad,1);
    
    I1=reshape(I,[],3);
    
    b1=mL(1);b2=mL(2);b3=mL(3);
    q1=mS(1);q2=mS(2);q3=mS(3);
    
    %% evaluating theta1 for Y
    A=(2/3.0)*(b2+b3-q2-q3)-(1/3.0)*(b1-q1);
    B=(1/sqrt(3))*(b3-b2-q3+q2);
    C=(b1-q1);
    
    b=((C^2)-(A^2))/2;
    a=(A^2+C^2-2*A*C+4*B^2)/4;
    c=(A^2+C^2+2*A*C-4*B^2)/4;
    
    D=b^2-4*a*c;
    
    cost1a=(-b+(sqrt(D)))/(2*a);
    theta1a=acosd(cost1a);
    cost1b=(-b-(sqrt(D)))/(2*a);
    theta1b=acosd(cost1b);
    %
    A=qrot3d(mL,[1,1,1],deg2rad(real(theta1a)))-qrot3d(mS,[1,1,1],deg2rad(real(theta1a)));
    B=qrot3d(mL,[1,1,1],deg2rad(real(theta1b)))-qrot3d(mS,[1,1,1],deg2rad(real(theta1b)));
    if abs(A(1))<=abs(B(1))
        theta1=theta1a;
    else
        theta1=theta1b;
    end
    
    ang=deg2rad(theta1);
    rotdata = qrot3d(im2double(I1),[1,1,1],ang);
    Ians1=reshape(rotdata,size(I,1),size(I,2),3);
    Output(:,:,1)=im2uint8(Ians1(:,:,1));
   
    %% evaluating theta2 for Cb
    A=(2/3.0)*(b1+b3-q1-q3)-(1/3.0)*(b2-q2);
    B=(1/sqrt(3))*(b1-b3-q1+q3);
    C=(b2-q2);
    
    b=((C^2)-(A^2))/2;
    a=(A^2+C^2-2*A*C+4*B^2)/4;
    c=(A^2+C^2+2*A*C-4*B^2)/4;
    
    D=b^2-4*a*c;
    cost2a=(-b+(sqrt(D)))/(2*a);
    theta2a=acosd(cost2a);
    cost2b=(-b-(sqrt(D)))/(2*a);
    theta2b=acosd(cost2b);
    
    A=qrot3d(mL,[1,1,1],deg2rad(real(theta2a)))-qrot3d(mS,[1,1,1],deg2rad(real(theta2a)));
    B=qrot3d(mL,[1,1,1],deg2rad(real(theta2b)))-qrot3d(mS,[1,1,1],deg2rad(real(theta2b)));
    if abs(A(2))<=abs(B(2))
        theta2=theta2a;
    else
        theta2=theta2b;
    end
    
    ang=deg2rad(theta2);
    rotdata = qrot3d(im2double(I1),[1,1,1],ang);
    Ians1=reshape(rotdata,size(I,1),size(I,2),3);
    Output(:,:,2)=im2uint8(Ians1(:,:,2));
    
    %% evaluating theta3 for Cr
    A=(2/3.0)*(b1+b2-q1-q2)-(1/3.0)*(b3-q3);
    B=(1/sqrt(3))*(b2-b1-q2+q1);
    C=(b3-q3);
    
    b=((C^2)-(A^2))/2;
    a=(A^2+C^2-2*A*C+4*B^2)/4;
    c=(A^2+C^2+2*A*C-4*B^2)/4;
    
    D=b^2-4*a*c;
    cost3a=(-b+(sqrt(D)))/(2*a);
    theta3a=acosd(cost3a);
    cost3b=(-b-(sqrt(D)))/(2*a);
    theta3b=acosd(cost3b);
    
    A=qrot3d(mL,[1,1,1],deg2rad(real(theta3a)))-qrot3d(mS,[1,1,1],deg2rad(real(theta3a)));
    B=qrot3d(mL,[1,1,1],deg2rad(real(theta3b)))-qrot3d(mS,[1,1,1],deg2rad(real(theta3b)));
    if abs(A(3))<abs(B(3))
        theta3=theta3a;
    else
        theta3=theta3b;
    end
    
    ang=deg2rad(theta3);
    rotdata = qrot3d(im2double(I1),[1,1,1],ang);
    Ians1=reshape(rotdata,size(I,1),size(I,2),3);
    Output(:,:,3)=im2uint8(Ians1(:,:,3));
    
    clearvars a b c A B C cost1a cost2a cost3a cost1b cost2b cost3b b1 b2 b3 q1 q2 q3 theta1a theta1b theta1c theta2a theta2b theta2c theta3a theta3b theta3c;
    
    %% Regain the colors
    
    if exist('det','var')
    d=im2bw(det);d=double(~(d));d(:,:,2)=d;d(:,:,3)=d(:,:,2);
    
    dip=uint8(d).*(img+1);   
    R = dip(:,:,1);
    G = dip(:,:,2);
    B = dip(:,:,3);
    newpixs = [R(:) G(:) B(:)];
    newpixs( all(~newpixs,2), : ) = [];
    
    p=divisors(nnz(d)/3);
    a=p(length(p)/2);
    b=nnz(d)/(3*a);
    
    newpixs=reshape(newpixs,a,b,3);
   
    Id=double(I).*d;
    Io=double(Output).*d;
    Iod=(255-double(Output)).*d;
    if mean2(abs(Id(:,:,1)-Io(:,:,1)))>mean2(abs(Id(:,:,1)-Iod(:,:,1)))
        Output(:,:,1)=255-Output(:,:,1);
    end
    if mean2(abs(Id(:,:,2)-Io(:,:,2)))>mean2(abs(Id(:,:,2)-Iod(:,:,2)))
        Output(:,:,2)=255-Output(:,:,2);
    end
    if mean2(abs(Id(:,:,3)-Io(:,:,3)))>mean2(abs(Id(:,:,3)-Iod(:,:,3)))
        Output(:,:,3)=255-Output(:,:,3);
    end
    end
    if (exist('det','var')==0)
        newpixs=img;
    end
        figure;
        subplot(2,4,1),imshow(I(:,:,1),[]),title('Y');
        subplot(2,4,2),imshow(I(:,:,2),[]),title('Cb');
        subplot(2,4,3),imshow(I(:,:,3),[]),title('Cr');
        subplot(2,4,4),imshow(I),title('YCbCr');
        subplot(2,4,5),imshow(Output(:,:,1),[]),title('new Y');
        subplot(2,4,6),imshow(Output(:,:,2),[]),title('new Cb');
        subplot(2,4,7),imshow(Output(:,:,3),[]),title('new Cr');
        subplot(2,4,8),imshow(Output),title('new YCbCr');

    %% Normalise Output
    minY=min(min(Output(:,:,1)));
    maxY=max(max(Output(:,:,1)));
    minI=min(min(I(:,:,1)));
    maxI=max(max(I(:,:,1)));
    
    X(:,:,1)=minI+((Output(:,:,1)-minY) .*((maxI-minI)/(maxY-minY)));
    
    minY=min(min(Output(:,:,2)));
    maxY=max(max(Output(:,:,2)));
    minI=min(min(I(:,:,2)));
    maxI=max(max(I(:,:,2)));
    
    
    X(:,:,2)=minI+((Output(:,:,2)-minY) .*((maxI-minI)/(maxY-minY)));
    
    minY=min(min(Output(:,:,3)));
    maxY=max(max(Output(:,:,3)));
    minI=min(min(I(:,:,3)));
    maxI=max(max(I(:,:,3)));
    
    X(:,:,3)=minI+((Output(:,:,3)-minY) .*((maxI-minI)/(maxY-minY)));
    
    Output=X;
    fi=ycbcr2rgb(Output);
    final=fi;
   
    minY=min(min(final(:,:,1)));
    maxY=max(max(final(:,:,1)));
    minI=min(min(newpixs(:,:,1)));
    maxI=max(max(newpixs(:,:,1)));
    
    X(:,:,1)=minI+((final(:,:,1)-minY) .*((maxI-minI)/(maxY-minY)));
    
    minY=min(min(final(:,:,2)));
    maxY=max(max(final(:,:,2)));
    minI=min(min(newpixs(:,:,2)));
    maxI=max(max(newpixs(:,:,2)));
    
    X(:,:,2)=minI+((final(:,:,2)-minY) .*((maxI-minI)/(maxY-minY)));
    
    minY=min(min(final(:,:,3)));
    maxY=max(max(final(:,:,3)));
    minI=min(min(newpixs(:,:,3)));
    maxI=max(max(newpixs(:,:,3)));
    
    X(:,:,3)=minI+((final(:,:,3)-minY) .*((maxI-minI)/(maxY-minY)));
    
    %% Color transfer and Output
    
    fil=im2uint8(color_transfer2(newpixs,X));
    figure,imshow(fil)
    imwrite(fil,[outpath,[num2str(ar(index)),'.png']]);
end