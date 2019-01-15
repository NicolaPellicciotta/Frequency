function [r_bin,cf,f_res]= FreqCorrFunc_fft(freq,f_filter)

        %%% calculate the correlation function from a frequency map
        %%% obtained with fft 
        BoxSize=4;

        if nargin < 2 || isempty(f_filter)
            f_filter= ones(size(freq));
        end
         r= repmat(1:size(freq,1)*BoxSize,[size(freq,2),1]);
         c=r';
         r=r(:);c=c(:)
         r=r(f_filter);c=c(f_filter);
         
         R = sqrt(((repmat(r,[1,numel(r)])) - (repmat(r,[1,numel(r)]))').^2 ...
           + ((repmat(c,[1,numel(c)])) - (repmat(c,[1,numel(c)]))').^2 );      
         
         F= nanmean( freq(:)); 
         MatF=repmat(freq,[1,numel(freq)]);
         
         CF= ( (MatF-F).*(MatF'-F));
        
         R_bin= (min(R)+1):(BoxSize):(max(R(:))*2/3);
         
         for i=1:(numel(R_bin)-1)
             
 %           slot=(R>=R_bin(i) & R < (R_bin(i+1)) );
             slot=(R < (R_bin(i+1)) );
             cf(i)= mean(CF(slot)); 
             N_f(i)=  var(MatF(slot)-F);   
             cf(i)= cf(i)/N_f(i);
             N(i)=sum(slot(:));
  %           ecf(i)=  nanstd(CF (R>=R_bin(i) & R < (R_bin(i+1)) ))/sqrt(numel(CF (R>=R_bin(i) & R < (R_bin(i+1)) )));
             r_bin(i)= (R_bin(i)+R_bin(i+1))/2; 

         end

          
         %%%% fit with exponential
         
         ft = fittype('a*exp(-x/b)+c','independent','x');
         fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,BoxSize,-1],...
               'Upper',[1,Inf,1],...
               'StartPoint',[0.9, BoxSize*5, 0]);
         
         f_res = fit(r_bin',cf',ft,fo);
         
         
end



