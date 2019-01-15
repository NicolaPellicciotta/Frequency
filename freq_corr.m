
%% Load Data from Plate 3

path_dir = '/u/homes/np451/Desktop/ependymal data/November_2017/3.11/Plate3/analysis/'; cd(path_dir); 

subdir=1:1;
mu= [0.75, 0.25, 0.5, 1 , 0 , 1.5 ...
    , 1 , 1.5 , 0.75 ,.5 ,.25 , 0 ...
    ,.25 , 0, 1 , .75 , 1.5 , 0.5 ];
clear c;

BoxSize=64;
px2mu=0.146  %%% 40X

%F=zeros([numel(subdir),12]);
for jj=1:numel(subdir)
    exp=num2str(subdir(jj));
    cd(path_dir); cd(exp);

    
    %%%% load frequency from each exp
    d=dir('*.movie');

    for nf=1:size(d,1)
        filename= d(nf).name;
        Res_filename= strcat(filename(1:end-5),'mat');
        cd(path_dir);cd(exp); cd(filename);   
        load(Res_filename);
        load('freq.mat');
%        freq=ones(size(freq))%+(rand(size(freq))-0.5); 

        F(jj,nf)=nanmean( freq(:)); 
%        c{jj,nf}=freq(:);

        
        
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
         
         R = sqrt(((repmat(r,[1,numel(r)])) - (repmat(r,[1,numel(r)]))').^2 ...
           + ((repmat(c,[1,numel(c)])) - (repmat(c,[1,numel(c)]))').^2 );      
         
         F= nanmean( freq(:)); 
         MatF=repmat(freq,[1,numel(freq)]);
         
         CF= ( (MatF-F).*(MatF'-F));
        
         R_bin= (min(R)+1):(BoxSize):(max(R(:))*2/3);
         
         for i=1:(numel(R_bin)-1)
             
 %            slot=(R>=R_bin(i) & R < (R_bin(i+1)) );
             slot=(R < (R_bin(i+1)) );
             cfs(i)= mean(CF(slot)); 
             N_f(i)=  var(MatF(slot)-F);   %nansum(sqrt(MatF(slot)-F).^2)*nansum(sqrt(MatF(slot)'-F).^2);
             cf(i)= cfs(i)/N_f(i);
             N(i)=sum(slot(:));
  %           ecf(i)=  nanstd(CF (R>=R_bin(i) & R < (R_bin(i+1)) ))/sqrt(numel(CF (R>=R_bin(i) & R < (R_bin(i+1)) )));
             r_bin(i)= (R_bin(i)+R_bin(i+1))/2; 
         end
         
         %errorbar(r_bin*px2mu,cf,ecf,'-o');xlabel('distance [um]'); ylabel('Freq coor');
         plot(r_bin*px2mu,cf,'-o');xlabel('distance [um]'); ylabel('Freq coor');
         
         
     end
    
end
