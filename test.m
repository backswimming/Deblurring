clear;
clc;
close all;

%% to choose different data and shots
num_shot=3;
number=20;

%% basic file reading and processing
path=[num2str(num_shot),'sh/'];
filename_tra='k_points.mat';
filename_info='parout1.txt';
load([path,filename_tra]);
samples=textread([path,filename_info],'%s',1,'headerlines',16);
samples=samples{1};
index=strfind(samples,'=');
samples=samples(index+1:end);
samples=str2double(samples);
interval=textread([path,filename_info],'%s',1,'headerlines',11);
interval=interval{1};
index=strfind(interval,'=');
interval=interval(index+1:end);
interval=str2double(interval);
during=samples*interval;
load([path,'im_b0_CG']);
img=im_b0_CG;
m_b=img(:,:,:,number);
m_b=sum(m_b,3);
[N,~]=size(m_b);

%% field map
path1='Pre/';
path2='Post/';
fold1='B0_map_2.0mm/';
fold2='mDixon_3.5mm/';
filename='field_map_dicom';
load([path2,fold2,filename]);
b0=field_map(:,:,number);
% clear field_map;
% load([path2,fold2,filename]);
% b0=(field_map(:,:,number)+b0)/2;
delta_f=imresize(b0,[N,N],'bilinear');
delta_f=imrotate(fliplr(delta_f),90);

%% defining const&info
interval = interval/1000;
read_out = during/1000;
% delta_f=rot90(b0_q,3);
% delta_f=flip(delta_f,2);
% delta_f(abs(delta_f)<10)=0;
% figure;
% imshow(abs(m_b),[]);
figure;
imagesc(delta_f); 
caxis([-50,50]);
colorbar;
title('delta_f');
axis image;
% figure;
% xx=1:168;
% [X,Y] = meshgrid(xx,xx);
% surf(X,Y,delta_f);

%% calculate basic MTF & PSF
t1 = time_map1(N, k_points, interval);
t2 = time_map2(N, read_out);
t3 = time_map3(N, read_out);

%% using H to do blurring and deblurring
guard=1;
for i=1:N
    for j=1:N
        f_map=delta_f(i,j)*ones(N,N);
        H_k = exp(1i*2*pi*(f_map.*t1));
        H_k = fftshift(H_k);
        h = fftshift(ifft2(H_k));
        max_h=max(max(abs(h)));
        threshold=0.0015*max_h;
        h(abs(h)<threshold) = 0;
        if(mod(N,2)==1)
            np=(N-1)/2;
            h_p=padarray(h,[np,np],'both');
        else
            np=(N-2)/2;
            h_p=padarray(h,[1,1],'post');
            h_p=padarray(h_p,[np,np],'both');
        end
        tmp=h_p(i:i+N-1,j:j+N-1);
        h_t=reshape(tmp,1,[]);
        len=sum(h_t~=0);
        index_h = find(h_t~=0);
        x(guard:guard+len-1)=(j-1)*N+i;
        y(guard:guard+len-1)=index_h;
        z(guard:guard+len-1)=h_t(index_h);
        guard=guard+len;
    end
end
H=sparse(x,y,z);

%% deblurring
clear x y z
flag=0;
m_b_1=reshape(m_b,[],1);
figure;
subplot(2,2,1);
m_blur_CG = squeeze(sos(img(:,:,:,number),3));
imshow(abs(m_blur_CG),[]);
title('blurring-image');
load('T2W/imim.mat');
fig=double(imim(:,:,number));
fig=imresize(fig,[N,N],'bilinear');
fig=imrotate(fliplr(fig),90);
subplot(2,2,2);
imshow(fig,[]);
title('reference');
fig=rot90(fig,2);
fig_1=reshape(fig,[],1);
m_fig=H*fig_1;
m_fig=reshape(m_fig,N,N);
subplot(2,2,3);
imshow(abs(m_fig),[]);
title('blurring-in-theory');


m_b=squeeze(sosfit(img(:,:,:,number),3));
m_b_1=reshape(m_b,[],1);
xi=H'*m_b_1;
ri=xi-H'*(H*xi);
pi=ri;
N_i=4;
for i=1:N_i
%     m_tmp=reshape(xi,N,N);
%     m_tmp=rot90(m_tmp,2);
%     subplot(3,4,i+1);
%     imshow(abs(m_tmp),[]);
%     title(['iteration-',num2str(i)]);
    ai=ri.'*ri/(pi.'*H'*H*pi);
    xi=xi+ai*pi;
%     er=double(abs(xi-m_1));
%     er=er.^2; 
%     err(i)=sum(er);
%     if(i>1)
%         if((err(i)>err(i-1))&&(~flag))
%             m_d_1=xi;
%             flag=1;
%         end
%     end
    rj=ri-ai*H'*(H*pi);
    pi=rj+rj.'*rj*pi/(ri.'*ri);
    ri=rj;
end
m_deblur=rot90(reshape(xi,N,N),2);

figure;
imshow(abs(m_deblur),[]);
title('deblurring-invivo_r-3.5mm');
