function b = real2jpg(a)

b=(a-min(a(:)))./(max(a(:))-min(a(:)))*255;%%aÎªdoubleĞÍ
b=uint8(b);