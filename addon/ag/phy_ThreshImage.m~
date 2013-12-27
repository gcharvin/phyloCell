function [imdist imbw]=phy_ThreshImage(imdata,thr,C)


imdata=phy_scale(imdata);

maxe=max(max(imdata));
meane=mean2(imdata);

imbw=im2bw(imdata,double(meane+thr*(maxe-meane)));
    
imdist=bwdist(imbw|~C);

imdist(~C)=0;
