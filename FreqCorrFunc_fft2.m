function [r_bin,cf,f_res]= FreqCorrFunc_fft2(freq,mask,bsz)

if nargin < 2 || isempty(mask)
   mask= ones(size(freq));
   mask=mask(:);
end
mask=mask(:);
 
Rres=4; 
jj = repmat((1:size(freq,1)),size(freq,2),1);
ii=jj';
F= nanmean(freq(mask));
Normf= (nanstd(freq(mask)))^2;
CF= zeros([(  ceil(sqrt(size(freq,1)^2 + size(freq,2)^2)./Rres)+1 ),1]);
DC= zeros([(  ceil(sqrt(size(freq,1)^2 + size(freq,2)^2)./Rres)+1 ),1]);

for f=1:numel(freq(:))
    clc; disp(100*f/numel(freq(:)))
    if mask(f)~=0; 
        [ci,cj] = ind2sub(size(freq),f);
        distance_map =ceil((round(sqrt((ii-ci).*(ii-ci)+(jj-cj).*(jj-cj)))+1)./(Rres)); %if we fftshit here is only 1 time instead of Nframes/2
        distance_map= distance_map(:);

    
        A = accumarray(distance_map(mask), (freq(ci,cj)-F).*(freq(mask)-F));
        CF(1:numel(A))= CF(1:numel(A))+ A; 
    
        dist_counts = accumarray(distance_map(mask),ones(sum(mask(:)),1));
        DC(1:numel(dist_counts))= DC(1:numel(dist_counts))+dist_counts;
    end   
end

cf= CF./DC; cf=cf/Normf;
r_bin = bsz*(1: ceil(sqrt(size(freq,1)^2 + size(freq,2)^2)./Rres)+1)*Rres;

%[m,im]=max(DC(:));
%r_bin= r_bin(1:im);
%cf = cf(1:im);
 %%% fit with exponential

l = find(cf < 0, 1);
if l==0; l=floor(numel(cf)/2);end;
y=cf(1:l);
x=r_bin(1:l);

 ft = fittype('a*exp(-x/b)+c','independent','x');
 fo = fitoptions('Method','NonlinearLeastSquares',...
 'Lower',[0,min(r_bin(:)),-10],...    %%% [a,b,c]
 'Upper',[Inf,max(r_bin(:))*2,1],...
 'StartPoint',[1,max(r_bin(:))/2, 0]);
%          
 f_res = fit(x',y,ft,fo);
end