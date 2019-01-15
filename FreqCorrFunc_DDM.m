function [r_bin,cf,f_res]= FreqCorrFunc(cilia,freq,f_filter)
        
        BoxSize=64;
        N_Boxes_row = floor(cilia.height / BoxSize);   %number of boxes that fit (vertically) in the field of view given BoxSize
        N_Boxes_col = floor(cilia.width  / BoxSize);    %number of boxes that fit (horizontally) in the field of view given BoxSize
        N_Boxes = N_Boxes_row * N_Boxes_col;      
        col_offset = 1;
        row_offset = 1;
        
        clear r; clear c;
        %%% freq - pos
         for i = 1:N_Boxes_row
            for j = 1:N_Boxes_col
                    
                    ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
                    r(ii) =  (row_offset + (i-1)*BoxSize + row_offset + (i)*BoxSize -1)/2;
                    c(ii) =  (col_offset + (j-1)*BoxSize + col_offset + (j)*BoxSize -1)/2;
             end
         end
          
         r=r(cilia.Results.ind_good_boxes); c=c(cilia.Results.ind_good_boxes); r=r(:);c=c(:);
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