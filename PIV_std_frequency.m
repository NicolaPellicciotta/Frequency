long= {'0/BF','100/BF','200/BF'};
FL_long={'0/FL','100/FL','200/FL'};

aaa=0;
for insert=long;
insert=insert{1};
aaa=aaa+1;

data_dir = strcat('/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/18.6.18/',insert);
a_folder=strcat('/home/np451/Desktop/ependymal/June2018/18.6.18/PIV/',insert); 
piv_folder= strcat('/home/np451/Desktop/ependymal/June2018/18.6.18/PIV/',FL_long{aaa});

cd(data_dir); 
suffix='.movie';
direc = dir('*BF*.movie');
N_files= size(direc,1);

%%

for i=1:N_files
    cd(data_dir)
    filename = direc(i).name;
    exp_name=filename(1:end-6);
    cd(a_folder);
    if exist(strcat(a_folder,'/',filename)) == 0
        cd(data_dir)    
        mo=moviereader(filename);
        FR=mo.FrameRate;       
        fs=mo.read();
        s=std(double(fs(:,:,1:200)),[],3);
        p_val=1.5; 
        ss=s;ss(ss<p_val)=0;ss=logical(ss);
        N_frames=size(fs,3);
        fs_roi= fs(repmat(ss,[1,1,N_frames]));
        roi=reshape(fs_roi,[sum(ss(:)),N_frames]);
        roi=double(roi)- mean(roi,2);

        roi=roi';  %%% for periodogram, he likes in this way
        N_seg=3;
        window = hann(floor(N_frames)/N_seg);
        fft_a= zeros([size(roi,1),numel(window)]);
        Npx=1000;
        N_rep= floor(size(roi,2)/Npx);
        
        clear Pxx;
        for ff=1:N_rep
            ind_l=(ff-1)*Npx +1;ind_r=(ff)*Npx; 
            temp_roi= roi(:,ind_l:ind_r);
            [pxx, fq]=pwelch(double(temp_roi),window,50,size(temp_roi,1),FR);
            Pxx(:,ind_l:ind_r)=pxx;
        end
        
        temp_roi= roi(:,ind_r+1:end);
        [pxx, fq]=pwelch(double(temp_roi),window,50,size(temp_roi,1),FR);
        Pxx(:,ind_r+1:ind_r+size(pxx,2))=pxx;
        
        
        Pxx=Pxx';
        [col,row]=meshgrid(1:size(ss,2),1:size(ss,1));
        col=col(ss==1);
        row=row(ss==1);
        
        %%%% average in boxes PIV-like
        bbs=32;
        bs=bbs/2;
        
        indxL=strfind(filename,'_X')+2;
        indyL=strfind(filename,'Y')+1;
        indxR=strfind(filename,'Y')-1;
        indyR=strfind(filename,'18Jun')-2;
        posx= str2num(filename(indxL :indxR));
        posy= str2num(filename(indyL:indyR));
    
        %%%%%load the std from BF videos
        cd(piv_folder);
        d_FL=dir(strcat('*X',num2str(posx),'Y',num2str(posy),'*.mat') );
        FL_res= load(d_FL.name);
                
 %       gridx= 1:bbs:(floor(size(fs,2)/bbs)*bbs);
 %       gridy= 1:bbs:(floor(size(fs,1)/bbs)*bbs);
 %       bXX= repmat(gridx,[numel(gridy),1])+bs;%bXX=bXX';
 %       bYY= repmat(gridy,[1,numel(gridx)])+bs;%bYY=bYY';
        
        x=FL_res.x;
        y=FL_res.y;
        bXX=x;
        bYY=y;
        %%%%% average over a box and find peak
         for b=1:numel(bXX);
           ind= (col > bXX(b)-bs) & (col <= bXX(b)+bs) & (row> bYY(b)-bs) & (row <= bYY(b)+bs); 
           fft_m= mean(Pxx(ind(:),:),1);
           if sum(ind)/(32^2)>0.08
               
                ind_fq= fq>10 & fq<40;
                [pks,locs] = findpeaks(fft_m(ind_fq),fq(ind_fq));                
                [~,where] = max(pks);
                frequency_strip(b) = locs(where);
           else
               frequency_strip(b) = nan;
           end
        
         end
        
        freq_strip=reshape(frequency_strip,size(bXX));
        cd(a_folder);
        save(strcat('freq_',exp_name,'.mat'),'freq_strip','ss');
        
    end
end
end

%% frequency map

%%%%%%%%%% Load data to find a complete map of the frequencies in the
%%%%%%%%%% channel

%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date
ppm= 0.14*2;  %%% 20X


%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/18.6.18/PIV/';

cd(path_dir); 

%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
%BFdir={'0/BF','30A/BF','30B/BF','100A/BF','100B/BF','1000/BF'};
subdir= {'0/FL','100/FL','200/FL'};
BFdir={'0/BF','100/BF','200/BF'};


%BFdir = {'100B/BF'};

for jj=1:numel(subdir)
    
start=0;
cd(path_dir);cd(subdir{jj}); d=dir('*FL*.mat');

clear xx;clear yy;clear uu;clear vv;clear uu1;clear vv1;clear ff;
clear XX;clear YY;clear UU;clear VV;clear UU1;clear VV1;


for ii=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    filename= d(ii).name;

    cd(path_dir);cd((subdir{jj}));
    load(filename);
    subdir= {'0/FL','100/FL','200/FL'};
    BFdir={'0/BF','100/BF','200/BF'};
    cd(path_dir);cd(subdir{jj});
    d=dir(strcat('*FL*.mat'));
    
    indxL=strfind(filename,'_X')+2;
    indyL=strfind(filename,'Y')+1;
    indxR=strfind(filename,'Y')-1;
    indyR=strfind(filename,'18Jun')-2;
    posx= str2num(filename(indxL :indxR));
    posy= str2num(filename(indyL:indyR));
    
    %%%%%load the std from BF videos
    cd(path_dir);cd(BFdir{jj});
    d_BF=dir(strcat('freq*X',num2str(posx),'Y',num2str(posy),'*.mat') );
    BF_res= load(d_BF.name);
    ss=BF_res.ss;
    freq_strip=BF_res.freq_strip;
    cd(path_dir);cd((subdir{jj}));
    
    ss(ss>0)=1;
    
    bs=16;
    for b=1:numel(x);
        y_box=floor(y(b)-bs);
        x_box=floor(x(b)-bs);
        temp_ss =ss(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs),:);
        box_ss(b)=mean(temp_ss(:));
    end
    
    ind= box_ss>0.08;
        
    xx(start+1:start+numel(x(ind))) = x(ind) +posy/ppm;
    yy(start+1:start+numel(x(ind))) = y(ind) -posx/ppm;
    um=nanmean(u,3);
    vm=nanmean(v,3);
    u1m=nanmean(u1,3);
    v1m=nanmean(v1,3);
    
    
    uu(start+1:start+numel(x(ind))) = um(ind);
    vv(start+1:start+numel(x(ind))) = vm(ind);
    uu1(start+1:start+numel(x(ind))) = u1m(ind);
    vv1(start+1:start+numel(x(ind))) = v1m(ind);
    ff(start+1:start+numel(x(ind))) = freq_strip(ind);
%%%%%%%copiato
       
    start= start+ numel(x(ind));
end


 XX= xx;
 YY= yy;
 UU= uu;
 VV= vv;
 UU1= uu1;
 VV1= vv1;

M= sqrt(UU1.^2 +VV1.^2);
nu= UU1./M;
nv= VV1./M;

% %%%%%% end correlation
% 
% %%%%% polarisation
% 
M= sqrt(UU1.^2 +VV1.^2);
nu= UU1./M;
nv= VV1./M;
Mm=median(M);
pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2);
UM= nanmean(nu);
VM= nanmean(nv);

save('freq_post.mat')






cd(path_dir);
end


%% plots
%% this script find the ciliated cell density in the channel and polarisation at different distance from the channel

ppm= 0.14*2;  %%% 20X
%path_dir = '/home/np451/Desktop/ependymal/June2018/13.6.18/PIV/';
path_dir = '/home/np451/Desktop/ependymal/June2018/18.6.18/PIV/';
cd(path_dir); 
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
subdir= {'0/FL','100/FL','200/FL'};

xres=11;
x_slices= (1:xres)*3800/xres;
for sub=1:numel(subdir)

cd(path_dir);cd(subdir{sub});
load('freq_post.mat');
%subdir= {'0/FL','30A/FL','30B/FL','100A/FL','100B/FL','1000/FL'};
subdir= {'0/FL','100/FL','200/FL'};
for g=1:(numel(x_slices)-1)
    ind_x= (XX>x_slices(g) & XX<=x_slices(g+1));
    Pcilia(g)= numel(XX(XX>x_slices(g) & XX<=x_slices(g+1)));
    pol_xslice(g)= sqrt(nanmean(UU1(ind_x))^2 + nanmean(VV1(ind_x))^2); 
    FF(g)= nanmean(ff(ind_x)); 
end


R{sub}.Pcilia=Pcilia/(sum(Pcilia(:)));
R{sub}.x_bin = (x_slices(1:end-1) + x_slices(2:end))./2;
R{sub}.density=numel(XX);
R{sub}.pol= pol_xslice;
R{sub}.FF= FF;
clear Pcilia
end


%%%%%% plot polarisation at different distances from the channel
%leg={'0','30-A','30-B','100-A','100-B','200'};
%flow=[0,30,30,100,100,200];
leg={'0','100','200'};
flow=[0,100,200];


figure()
for b=1:numel(R)
 plot(R{b}.x_bin,R{b}.FF,'-o'); hold on;

end
xlabel('X position [px]'); ylabel('Cilia Probaility ditribution');

legend(leg)


%%%%%% plot polarisation at different distances from the channel

figure()
for b=1:numel(R)
 plot(R{b}.x_bin,R{b}.pol,'-o','LineWidth',2); hold on;

end
xlabel('X position [px]'); ylabel('polarisation');

legend(leg)

