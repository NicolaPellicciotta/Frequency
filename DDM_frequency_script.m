%%%% script DDM for frequency (on boxes) %%%
%%% to use this script you should load only Documents/MATLAB/DDM_cilia and Documents/MATLAB/DDM/MatLab_Common_Function 
%% 
path = '/home/np451/Desktop/Mouse/synchronisation/23.5.18/normal/';
cd(path);
 mkdir('analysis_multiDDM'); 
ana_dir= strcat(path,'analysis_multiDDM/');
%store_dir= '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/6.4.18/beads/';
store_dir= '/home/np451/Desktop/Mouse/synchronisation/23.5.18/normal/data'
cd(store_dir);
d=dir('*.movie');
box_size=[32];  %%% box size for ddm


for nf=1:size(d,1)
    filename= d(nf).name;
    cd(ana_dir);
    if exist(filename)==0;    
    cd(store_dir);    %%% moving to store to load frames

    cilia = DDM_Analysis_nico(filename);

    cd(ana_dir); mkdir(filename); cd(filename);  %%%% moving to the analysis dit

   
    cilia.VariableBoxSize_Analysis(box_size);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts

    % plot good boxes with std %
    end
end
%% extract frequency and plots from DDM
ana_dir = '/home/np451/Desktop/Mouse/MAR2018/analysis_multiDDM/'
cd(ana_dir);
d=dir('*.movie');
N_exp= numel(dir);


for nf=1:size(d,1)
    filename= d(nf).name;
    cd(ana_dir);cd(filename);
    load(strcat(filename(1:end-6),'.mat'));
    
    minbox= numel(cilia.Results)
    box_size= cilia.Results(minbox).BoxSize;
    
    new_std = cilia.std_fs;
    f_y=cilia.width;
    f_x=cilia.height;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ); 
    new_mask = oversample(cilia.Results(minbox).ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results(minbox).MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    freq=cilia.Results(minbox).MedianFrequencyVec(:);
    save('freq.mat','freq');
    

end

%% Load Data

ana_dir = '/home/np451/Desktop/Mouse/MAR2018/analysis_multiDDM/'
cd(ana_dir);
d=dir('*.movie');
N_exp= numel(dir);


for nf=1:size(d,1)
        filename= d(nf).name;
        cd(ana_dir);cd(filename);
        s = str2num(filename(end-7:end-6));
        m = str2num(filename(end-10:end-9));
        h = str2num(filename(end-13:end-12));
        t(nf)= h*3600 + m*60 + s; 
        
        if isempty(strfind(filename,'Pa'));
            shear(nf)=0;
        else
            shear(nf)=str2num(filename(strfind(filename,'Pa')-1)); 
        end
            
        
%        load(strcat(filename(1:end-6),'.mat')); cilia.gather_results;
%        freq=cilia.Results.MedianFrequencyVec(:);
%        save('freq.mat','freq');
%        F(jj,nf)=nanmean( cilia.Results.MedianFrequencyVec(:));        
        load('freq.mat');
        freq=freq(freq>4 & freq<16);
        F(nf)=nanmean( freq(:));
        N(nf)= numel(~isnan(freq(:)));
        std_F(nf)= nanstd(freq(:))/sqrt(N(nf));
        c{nf}=freq(:);
end    
shear= 10.^(-shear);
shear(shear==1)=0;
[t,I] =sort(t);
t=(t-min(t(:)))/60;
shear=shear(I);
F=F(I);
%%
plot(t,F,'o-');
hold on;
plot(t,shear*1000*max(F(:)));
ylim([min(F(:)),max(F(:))]);
%%

cd(path_dir);

%%% mean without weights with error:  F_t ; eF_t ; rv  (CBF, errCBF, viscosity)
F_t=[mean(mean(F(mu==0,:),1),2),mean(mean(F(mu==0.5,:),1),2),mean(mean(F(mu==1,:),1),2),mean(mean(F(mu==1.5,:),1),2),mean(mean(F(mu==2,:),1),2)];
std_t = [std(std(F(mu==0,:),[],1),[],2),std(std(F(mu==0.5,:),[],1),[],2),std(std(F(mu==1,:),[],1),[],2),std(std(F(mu==1.5,:),[],1),[],2),std(std(F(mu==2,:),[],1),[],2)];
std_N = [numel(F(mu==0,:)),numel(F(mu==0.5,:)),numel(F(mu==1,:)),numel(F(mu==1.5,:)),numel(F(mu==2,:))];
eF_t = std_t./sqrt(std_N);    
load('MC_mu.mat')
rv= M(2,:); %%% real viscosity Pa*s


%% mean with weights with error:  F_t ; eF_t ; rv  (CBF, errCBF, viscosity)
tCBF= zeros(size(mu));
MC=[0,0.5,1,1.5,2];
for k=1:numel(mu); tCBF(k)=sum( (F(k,:).*(std_F(k,:).^2))/(sum(std_F(k,:).^2)) ); teCBF(k)=sum( ((std_F(k,:).^2).*(std_F(k,:).^2))/(sqrt(N_exp)*(sum(std_F(k,:).^2))) ); end;

MC=[0,0.5,1,1.5,2];
clear CBF;clear eCBF;
for k=1:numel(MC); temp_CBF=tCBF(mu==MC(k)); temp_eCBF= teCBF(mu==MC(k))'; CBF(k)=[sum( temp_CBF.*temp_eCBF)/(sum(temp_eCBF))]; eCBF(k)=[std(temp_CBF)/sqrt(numel(temp_CBF))];end; %eCBF(k)=[sum( temp_eCBF.*temp_eCBF)/(sum(temp_eCBF))];end;

figure(); errorbar(rv,CBF,eCBF,'-o'); hold on;
title('CBF vs viscosity weighted average'); legend('exp point');
xlabel('{$\eta[Pa*s]$}','Interpreter','latex','FontSize',15);
ylabel('{$CBF [Hz]$}','Interpreter','latex','FontSize',15);
%saveas(gcf,'CBF_viscosity_lin_weighted','pdf');



%% linear 

figure(); plot(mu,F,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_allfield_viscosity.png');


figure(); plot(rv,F_t,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');


figure(); errorbar(rv,F_t,eF_t,'o');
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');

%% loglog


figure(); errorbar(rv,F_t,eF_t,'o'); hold on;
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_viscosity_lin','jpg');
p= polyfit(log(rv),log(F_t),1);
fit_f=polyval(p, log(rv));

plot(rv,exp(fit_f),'-');

legend('exp',['y=x^{',num2str(p(1),2),'}']);
set(gca,'xscale','linear','yscale','linear');
saveas(gcf,'CBF_viscosity_lin_fit','jpg');

set(gca,'xscale','log','yscale','log');
saveas(gcf,'CBF_viscosity_log_fit','jpg');



%% semilogx
figure(); errorbar(rv,F_t,eF_t,'d'); hold on;

%%%%% fit semilogx CBF = log(visc)

px= polyfit(log(rv),(F_t),1);
fit_fx=polyval(px, log(rv));
plot(rv,(fit_fx),'-');

set(gca,'xscale','log','yscale','lin');
legend('exp',['y=',num2str(px(1),2),'*log(x)']);
title('CBF vs viscosity'); xlabel('viscosity [Pa*s]');ylabel('CBF [Hz]');
saveas(gcf,'CBF_viscosity_logx','jpg');


%% correlate number of cells with frequency

for ii=1:size(c,1)
for i=1:size(c,2);
    plot(numel(c{ii,i}), nanmean(c{ii,i}),'o');hold on;
end;
end

%%
cc=zeros(size(c,1));
for ii=1:size(c,1)
for i=1:size(c,2);
    cc(ii)= cc(ii)+numel(c{ii,i}); 
end;
end

cc=cc(:,1);cc=cc(:);

cc=reshape(cc,[numel(cc)/2,2]);
c_m= mean(cc,2);
ec_m= std(cc,[],2)/2;















%% previuos
path_dir = '/home/np451/Desktop/ependymal data/6.6/NOFLOW/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


%%%% also for the other

path_dir = '/home/np451/Desktop/ependymal data/6.6/0.1mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end



%%%%%%%

path_dir = '/home/np451/Desktop/ependymal data/6.6/0.5mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


%%%%%

path_dir = '/home/np451/Desktop/ependymal data/6.6/1mlmin/'
cd(path_dir);
mkdir('analysis');
d=dir('*.movie');
R=[]

for nf=10:20%size(d,1);
    
    filename= d(nf).name;
    cd(path_dir);
    cilia=DDM_Analysis(filename);
    cd('analysis');

    cilia.VariableBoxSize_Analysis([box_size]);
    save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save sruslts



    % plot good boxes with std %


    new_std = cilia.std_fs;
    new_std=   new_std(1:floor(f_x/box_size)*box_size , 1:floor(f_y/box_size)*box_size ) 
    new_mask = oversample(cilia.Results.ind_good_boxes,box_size);
    mask_overlay(new_std, new_mask, [1 0 0], 0.3);
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-5),'jpg']);
    close 'all'

    % plot frequencies
    cilia.gather_results;
    histogram(cilia.Results.MedianFrequencyVec); %% plot histogram 
    fig = get(groot,'CurrentFigure');
    saveas(fig,[cilia.Filename(1:end-6),'_freq_.jpg']);
    
    
    
    
    % add frequencies to R
    R=cat(1,cilia.Results.MedianFrequencyVec(:));
end


