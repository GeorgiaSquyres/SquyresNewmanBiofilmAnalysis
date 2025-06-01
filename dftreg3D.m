function [output,Greg] = dftreg3D(buf1ft,buf2ft)
% function [output Greg] = dft3D(buf1ft,buf2ft,usfac);

% GRS edit to extent dftregistration code into 3D; does not support
% upsampling. Original dftreg docs copied below: 
%
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]

% Outputs
% output =  [net_x_shift,net_y_shift,net_z_shift]
% net_x_shift net_y_shift net_z_shift Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
[x,y,z]=size(buf1ft);
% compute cross-correlation
CC = ifftn(buf1ft.*conj(buf2ft));
clear buf1ft buf2ft
% find peak location
CCmax = max(max(max(CC)));
i = find(CC==CCmax);
[xloc,yloc,zloc] = ind2sub([x,y,z],i);
clear CC

% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
%{
% compute error
rfzero = sum(abs(buf1ft(:)).^2)/(x*y*z);
rgzero = sum(abs(buf2ft(:)).^2)/(x*y*z); 
error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
error = sqrt(abs(error));
% compute phase difference
diffphase=atan2(imag(CCmax),real(CCmax));
%}

% compute xyz offsets
xd2 = fix(x/2); 
yd2 = fix(y/2);
zd2 = fix(z/2);
if xloc > xd2
    x_shift = xloc - x - 1;
else
    x_shift = xloc - 1;
end

if yloc > yd2
    y_shift = yloc - y - 1;
else
    y_shift = yloc - 1;
end

if zloc > zd2
    z_shift = zloc - z - 1;
else
    z_shift = zloc - 1;
end

output=[x_shift,y_shift,z_shift];

% Compute registered version of buf2ft
%{
if (nargout > 1)
    [nx,ny,nz]=size(buf2ft);
    Nx = ifftshift([-fix(nx/2):ceil(nx/2)-1]);
    Ny = ifftshift([-fix(ny/2):ceil(ny/2)-1]);
    Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
    [Nx,Ny,Nz] = ndgrid(Nx,Ny,Nz);
    Greg = buf2ft.*exp(i*2*pi*(-x_shift*Nx/nx-y_shift*Ny/ny-z_shift*Nz/nz));
    Greg = Greg*exp(i*diffphase);
end
%}
return
