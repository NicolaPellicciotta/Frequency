ana_dir = '/home/np451/Desktop/Mouse/MAR2018/ATP/bau/';
plot_dir='/home/np451/Desktop/Mouse/MAR2018/ATP/plots';
cd(ana_dir);
d=dir('40X_pos1*.movie');
N_exp= numel(dir);
clear R;
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
            
        
        load(strcat(filename(1:end-6),'v2.mat'));
        if isempty(cilia.Results(5).row_offset)
            cilia.Results(5).row_offset=cilia.Results(4).row_offset;
            cilia.Results(5).col_offset=cilia.Results(4).col_offset;
        end
        cilia.gather_results;

        [bs,D,eD,dD]=Sigmoid(cilia);
        [bs,bs_ind] =sort(bs);
        R(nf).bs=bs; R(nf).D=D(bs_ind); R(nf).eD= eD(bs_ind); R(nf).dD=dD(bs_ind);

        
        %%% measure frequency using the smallest box 
        bs_vec=[];
        for jj=1:numel(cilia.Results); bs_vec(jj)= cilia.Results(jj).BoxSize;end
        [bsz,smaller_box]=min(bs_vec);
        smaller_box=5;
        
        freq= cilia.Results(smaller_box).MedianFrequencyVec(:);
        mask=freq>4 & freq<25;
        freq=freq(mask);
        R(nf).F=nanmean( freq(:));
        R(nf).N = numel(~isnan(freq(:)));
        R(nf).std_F= nanstd(freq(:))/sqrt(R(nf).N);
        R(nf).freq=freq(:);
 
        
        %%% measure dependence of Damping over q using the largest boxsize
       [max_bsz,larger_box]=max(bs_vec);
       R(nf).D_q= cilia.Results(larger_box).MedianDamping_vs_q;
        
        

        %% correlation spatial frequency
        figure();
        plot(R(nf).bs,1./R(nf).D,'o','LineWidth',2); hold on;
        xlabel('distance[px]');
        %xlim([0,512]);
        ylabel('1/damping');
        saveas(gca,strcat(filename(1:end-5),'Sigmoid.png'));
        close all


end    
[t,I] =sort(t);
t=(t-min(t(:)))/60;
shear=shear(I);
R=R(I);
cd(plot_dir)
save('variables_pos2DDM.mat');
%% nice plots with shear stress
figure();
yyaxis right
area(t,shear);
alpha(.2)
set(gca, 'YScale', 'Log','Xscale','Lin');
ylabel('shear stress [Pa.s]')
yyaxis left
plot(t,F,'o-','LineWidth',2);
hold on;
xlabel('time [s]');
ylabel('CBF [Hz]');

%% plot DDM_correlation function of the main frequencies

figure();hold on;
a1=1;a2=2;a3=3;
plot(R(a1).bs,R(a1).D,'-o','LineWidth',2); hold on;
plot(R(a2).bs,R(a2).D,'-o','LineWidth',2); hold on;
plot(R(a3).bs,R(a3).D,'-o','LineWidth',2); hold on;

% errorbar(R(a1).bs,R(a1).D,R(a1).eD,'-o','LineWidth',2); hold on;
% errorbar(R(a2).bs,R(a2).D,R(a2).eD,'-o','LineWidth',2); hold on;
% errorbar(R(a3).bs,R(a3).D,R(a3).eD,'-o','LineWidth',2); hold on;
% 
% 

xlabel('distance[px]');
%xlim([0,512]);
ylabel('freq correlation');

legend(num2str(R(a1).F),num2str(R(a2).F),num2str(R(a3).F));
set(gca, 'YScale', 'Lin','Xscale','Log');
%ylim([min(F(:)),max(F(:))]);
%saveas(gca,'pos1_sigmoids_2peak.png');

%% plot damping average with same freuency (for the cilia with more boxes)


figure();
Db=mean([R(1).D',R(3).D'],2);
Fb=mean([R(1).F,R(3).F]);
Da=mean([R(2).D',R(4).D'],2);
Fa=mean([R(2).F,R(4).F]);


plot(R(a1).bs,Db,'-o','LineWidth',2); hold on;
plot(R(a2).bs,Da,'-o','LineWidth',2); hold on;

xlabel('distance[px]');
%xlim([0,512]);
ylabel('freq correlation');

legend(num2str(Fb),num2str(Fa));
set(gca, 'YScale', 'Lin','Xscale','Lin');
saveas(gca,'pos1_sigmoids_averagePeak_moreboxes_32-950.png');

%%
%% plot DDM Damping vs q for different frequecies
px2mu=0.146; 
figure();hold on;
a1=1;a2=2;a3=3;a4=4
q= 2*pi*(1:numel(R(a1).D_q))/ (max_bsz*px2mu) ;
plot(q,R(a1).D_q,'-o','LineWidth',2); hold on;
plot(q,R(a2).D_q,'-o','LineWidth',2); hold on;
plot(q,R(a3).D_q,'-o','LineWidth',2); hold on;
plot(q,R(a4).D_q,'-o','LineWidth',2); hold on;

xlabel('mode[px]');
%xlim([0,512]);
ylabel('DDM damping [1/s]');

legend(num2str(R(a1).F),num2str(R(a2).F),num2str(R(a3).F),num2str(R(a4).F));
set(gca, 'YScale', 'Lin','Xscale','Log');
%ylim([min(F(:)),max(F(:))]);
saveas(gca,'pos1_D_q.png');



%%


a1=10;
filename='40X_middlebaseline.07Mar2018_14.19.25.mat';
q_max= size(cilia.Results(1).Box(1).Iqtau,1); bsz= q_max*2./(1:q_max);
Damp=cilia.Results(1).Box(1).Damping;
plot(bsz(1:numel(Damp)),1./Damp,'-o');
hold on;
plot(R(a1).bs,1./R(a1).D,'-o','LineWidth',2); hold on;

%% compare frequency correlation function with multi-DDM

ana_DDM = '/home/np451/Desktop/Mouse/MAR2018/analysis_multiDDM/';
ana_fft = '/home/np451/Desktop/Mouse/MAR2018/analysis_fft/';
cd(ana_DDM);
load('variables_DDM.mat');
Rddm=R;
cd(ana_fft);
load('variables.mat');


figure();hold on;
a1=10;
tc= 0.73*(1./Rddm(a1).D - min(1./Rddm(a1).D(:)))/abs(max(1./Rddm(a1).D(:))-min(1./Rddm(a1).D(:)))
%tc= 1./Rddm(a1).D ;
plot(Rddm(a1).bs,tc,'-o','LineWidth',2); hold on;
plot(r_bin, R(a1).cf,'-o','LineWidth',2);
xlabel('distance[px]');
xlim([0,512]);
ylabel('freq correlation');
legend('multiDDM','corrFunc');
