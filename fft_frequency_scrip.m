
path = '/home/np451/Desktop/Mouse/synchronisation/23.5.18/psico/';
cd(path);
mkdir('analysis_fft');
ana_dir= strcat(path,'analysis_fft/');
%store_dir= '/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/6.4.18/';
store_dir= '/home/np451/Desktop/Mouse/synchronisation/23.5.18/psico/data';


cd(store_dir);
d=dir('40X_*.movie');

for nf=1:size(d,1)
    filename= d(nf).name;
    cd(ana_dir);  
if exist(filename)~=7
    cd(store_dir); 
    mo=moviereader(filename);
    fs=mo.read(); 

    cd(ana_dir); mkdir(filename); cd(filename);  %%%% moving to the analysis dit

    bsz = 8;
    
% parameters
    height = floor(size(fs,1)/4)*4;
    width = floor(size(fs,2)/4)*4;
    Nframes = size(fs,3);
    FR = mo.FrameRate%  mo.FrameRate; % frame rate
    fs=fs(1:height,1:width,:);
% fft window
    window = hann(floor(Nframes/2));
% wwindow = hann(floor(Nframes/2));

% find how long will the frequency vector be
    [~,dummyf] = periodogram(AutoCorr(double(squeeze(fs(1,1,:))),floor(Nframes/2)),window,floor(Nframes/2),FR);
    Nfreqs = numel(dummyf);
    clear dummyf;

% create binning map

    binningmap = create_binning_map([height, width],bsz);

% N_frames, height*box matrix
    freq_ind = repmat((1:Nfreqs)',1,height*bsz);

% row vector with the index within a binning coloumn (always start from 1,
% then proper placement in frequency map will be done by frequency_map(:,cc) = frequency_strip;
    temp_ind = reshape(binningmap(:,1:bsz),height * bsz,1)';
% put into matrix form for accumarray
    temp_ind_mat = repmat(temp_ind,Nfreqs,1);
    
    
% initialise frequency map
    frequency_map = nan(height/bsz,width/bsz);

% for loop on columns of blocks

for cc = 1:ceil(width/bsz)
    fprintf('%.2d/%.2d', cc, ceil(width/bsz) );
    ccleft = (cc-1) * bsz + 1;
    ccright = ccleft + bsz-1;
    
    % take the first block, reshape it to be a N_frames, height*box matrix
    % (good for periodogram)
    % Each column is a time signal
    temp_fs = reshape(fs(:,ccleft:ccright,:),height * bsz,Nframes)'; 
    
    % subtract the mean
    temp_fs = double(temp_fs) - mean(temp_fs,'double');
    
    % take autocorrelation over time
    temp_fs = AutoCorr(temp_fs,floor(Nframes/2));
    
    % make Fourier transforms
    % pxx is a matrix with the spectrum of each pixel as a column
    [pxx, frequencies] = periodogram(temp_fs,window,size(temp_fs,1),FR);
%     [pwxx, frequencies] = pwelch(temp_fs,wwindow,[],numel(wwindow),FR);

    % average points of the spectra with same box index and same frequency
    boxavg_pxx = accumarray({freq_ind(:),temp_ind_mat(:)},pxx(:)) ./ (bsz*bsz);
    
    % times a downsampled spectrum
    % resample by averaging two close frequencies
%     resampled_boxavg_pxx = bin_matrix(boxavg_pxx,[2 1]); % maybe easier to do boxavg_pxx(1:2:end,:)
%     resampled_frequencies = bin_matrix(frequencies,[2 1]);
%     prod_boxavg_pxx = boxavg_pxx(1:size(resampled_boxavg_pxx,1),:) .* resampled_boxavg_pxx;
%     prod_frequencies = frequencies(1:size(resampled_boxavg_pxx,1));
    
    % peak detection
    frequency_strip = nan(size(boxavg_pxx,2),1);
    parfor ii = 1:size(boxavg_pxx,2)
        if any(isnan(abs(boxavg_pxx(:,ii))))==0;
        % find peaks
        [pks,locs] = findpeaks(boxavg_pxx(:,ii),frequencies);
        
        % index higher peak
        [~,where] = max(pks);
        
        % write
        frequency_strip(ii) = locs(where);
        end
    end %parfor
    
    frequency_map(:,cc) = frequency_strip;
    fprintf(repmat('\b',1,5))
end
save('freq.mat','frequency_map');
figure(); imagesc(frequency_map);colorbar;
saveas(gca,'frequency_map.pdf');
close all;
end

end

%%


%% Load Data

ana_dir = '/home/np451/Desktop/Mouse/MAR2018/analysis_fft/'
cd(ana_dir);
d=dir('*.movie');
N_exp= numel(dir);
clear R;
bsz=8;
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
        freq=frequency_map(frequency_map>4 & frequency_map<25);
        F(nf)=nanmean( freq(:));
        N(nf)= numel(~isnan(freq(:)));
        std_F(nf)= nanstd(freq(:))/sqrt(N(nf));
        R(nf).freq=freq(:);
        mask=frequency_map>4 & frequency_map<25;

        %% correlation spatial frequency
        [r_bin,cf,f_res]= FreqCorrFunc_fft2(frequency_map,mask,bsz );
        
        L(nf)= f_res.b;
        R(nf).cf=cf;
        R(nf).f_res=f_res;

        figure();
        plot(r_bin,cf,'o'); hold on;
        plot(f_res);
        xlabel('distance[px]');
        xlim([0,512]);
        ylabel('freq correlation');
        saveas(gca,strcat(filename(1:end-5),'_PostPro.png'));
        close all


end    
shear= 10.^(-shear);
shear(shear==1)=1e-4;
[t,I] =sort(t);
t=(t-min(t(:)))/60;
shear=shear(I);
F=F(I);
L=L(I);
R=R(I);
%%
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

%ylim([min(F(:)),max(F(:))]);

%%
figure();
hold on;
for i=1:numel(R);
    plot(R(i).cf,'o');
end