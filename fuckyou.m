clear;
clc;
close all;

%% to choose different data and shots
num_shot=3;
number=20;

%% basic file reading and processing
path=[num2str(num_shot),'sh/'];
load([path,'im_b0_CG']);
img=im_b0_CG;
m_blur_CG = squeeze(sos(img(:,:,:,number),3));
haha=real2jpg(m_blur_CG);
haha_p=haha(30:138,25:138);
haha(30:138,25:138)=histeq(haha_p);
imshow(haha);
title('blurring-image');