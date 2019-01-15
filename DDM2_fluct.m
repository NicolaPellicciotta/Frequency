function [ Iqtau, Iq2tau ] = DDM2_fluct(frame_stack, N_couple_frames_to_average, tau,flag )
%DDM2_core_for_anisotropy Computational core for DDM calculation


frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
N_px_frame = frame_size*frame_size;
fft2_norm_factor = 1/N_px_frame;
N_frames = size(frame_stack,3);

frame_stack= frame_stack(1:frame_size,1:frame_size,:);


if nargin < 2 || isempty(N_couple_frames_to_average)
    N_couple_frames_to_average = 100;
end


%% general-purpose variables

max_q = floor(frame_size/2);


%% window for fft2, following giavazzi's 2017 paper
% there is no (-1).^j term because I want the window to be ==1 at x==0, this
% is achieved by shifting the x of he funciont in the paper by
% frame_size/2, this is equivalent to multiplying by another (-1)^j that cancels out the first one

jj = repmat(1:frame_size,frame_size,1);
ii = jj';
cc = max_q + 1;
distance_map = fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);
dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

if flag
a = [0.3635819; 0.4891775; 0.1365995; 0.0106411];
j = (0:3)';

win2 = sum(a .* cos(2*pi*j.*distance_map(:)'./frame_size),1);
win2 = reshape(win2,frame_size .* [1 1]);
win2 = fftshift(win2);
else
    win2=1;
end


%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)


Iq2tau = zeros(frame_size, frame_size, N_frames-tau, 'double');


Iq2tau = abs(fft2( win2 .* ...
            (single(frame_stack(:,:,1:end-tau))-single(frame_stack(:,:,(tau+1):end))) ) * fft2_norm_factor).^2;

        
Iqtau = zeros(max_q, N_frames-tau, 'double');        
for i=1:(N_frames-tau)
    I_temp= Iq2tau(:,:,i);
    oneD_power_spectrum = accumarray(distance_map,I_temp(:))./dist_counts;	%radial average
    Iqtau(:,i) = oneD_power_spectrum(2:max_q+1); 
end


Iq2tau = fftshift(fftshift(Iq2tau,1),2); %put central mode in the middle

end