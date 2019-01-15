%% scrip to enlarge the range of boxsize in a cilia variable class

path = '/home/np451/Desktop/Mouse/MAR2018/ATP/';
cd(path);
 mkdir('analysis_multiDDM'); 
ana_dir= strcat(path,'analysis_multiDDM/');
store_dir= '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/15.3.18/';

cd(ana_dir);
d=dir('40X_*.movie');
%box_size=[32,64,128,256,700];  %%% box size for ddm
box_size=[8,16,512,950];  %%% box size for ddm


for nf=1:size(d,1)
    filename= d(nf).name;   
 

   cd(ana_dir); cd(filename);  %%%% moving to the analysis dit
   load(strcat(filename(1:end-6),'.mat'));
   cilia.load_movie();
   
    cilia.VariableBoxSize_Analysis(box_size);
    save([cilia.Filename(1:end-6),'2.mat'],'cilia');  %%% save sruslts

    % plot good boxes with std %
end