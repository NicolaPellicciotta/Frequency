function [Noise] = noise_fft(F_rest,temp_BW)
%%%%% noise measuring, idea taken from paper 
%Goldstain and Polin "Noise and Synchronization in Pairs of Beating
%Eukaryotic Flagella Raymond"
% it reads a movie file called noise that usually is a longer movie (20 sec) 
d=dir('*noise*movie');
Noise.movie='noise';

if isempty(d); d=dir('*P0*movie');Noise.movie='P0';end
mo=moviereader(d(1).name);
FR=mo.FrameRate;
fs=mo.read();
N_frames= size(fs,3);
Noise.FR=FR;Noise.N_frames=N_frames;

frequency_T = F_rest;
fq_min= F_rest-7;fq_max= F_rest+7;

%%%% set parameters for 
beat_rep=25;  %%%% for how many beat you want to use to calculate the frequency,should be large enough to have decent frequency resolution; 
N_beat_t_array= floor(beat_rep/5);  %%%% each fft time window will start after how many beats?
frames_beat=floor(FR/ frequency_T);  %%%%% how many frames last one beat
dt= frames_beat*beat_rep;  %%%%% number of frames total of the fft time window


%% find the mean frequency over a time dt, starting from time t_array(tt)
t_array= 1:frames_beat*N_beat_t_array:(N_frames-dt);%%% this are the times where the fft time window starts
frequency=zeros([1,numel(t_array)]);

for jj=1:numel(t_array)
tt=t_array(jj);
BW_tt= repmat(temp_BW,[1,1,dt]);  %%% repeating the mask over a number of frames dt
fs_tt=fs(:,:,tt:dt+tt);  %%% cut the frames to the time window 
fs_roi=fs_tt(BW_tt);    %%%% taking data only in the mask
fs_roi=reshape(fs_roi,[sum(temp_BW(:)),dt]);  %%%% reshaping to have the right size for the fft
roi=double(fs_roi)- mean(fs_roi,2); %%%% removing background
%%% find mean freq
window = hann(floor(dt));   %%%%% window for fft
window= repmat(window,[1,size(roi,1)])';
n= floor(N_frames);   %%%%% this I think it is an average, help to have smoother data
%if mod(n,2)==0; n= n-1;end
pxx= abs(fft(double(roi).*window,n,2)).^2;
m_pxx= mean(pxx(:,1:floor(n/2)),1);
fq= (0:(FR./n):(FR./2-FR./n));    

f_range=fq>fq_min & fq< fq_max;
[pks,locs,w,p] = findpeaks(m_pxx(f_range),fq(f_range));%%%% find the peaks frequency in the selected freq range
baseline= mean((m_pxx(f_range)));
[~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index     %%%%% load results in Rc
[~,where] = max(pks);
frequency(jj)=locs(end);
%%%%%%
    
end

%%%%% noise calculation
Noise.N_beat_t_array= N_beat_t_array;
Noise.t_array=t_array;
Noise.frames_beat=frames_beat;
Noise.beat_rep=beat_rep;
Noise.dt=dt;
Noise.frequencies= frequency;
Noise.value= beat_rep*var(frequency)/(mean(frequency)^2); 
Noise.f1_min=fq_min;
Noise.f1_min=fq_max;
end