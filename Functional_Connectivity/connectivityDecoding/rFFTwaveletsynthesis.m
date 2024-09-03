function x=rFFTwaveletsynthesis(w,FFTsynthesisfilters,J);

% rFFTWAVELETANALYSIS FFT-based implementation of the redundant inverse wavelet transform.
% 	
% 	Author: Dimitri Van De Ville
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	
% 	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, Vol. 42, No. 1, pp. 43--67, January 2000.
% 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
% 	Proc. SPIE, Wavelet Applications in Signal and Image Processing VII,
%     Denver, CO, USA, 19-23 July, 1999, vol. 3813, pp. 422-431. 
% 	[3] T. Blu and M. Unser, "The fractional spline wavelet transform: definition and 
%	implementation," Proc. IEEE International Conference on Acoustics, Speech, and 
%	Signal Processing (ICASSP'2000), Istanbul, Turkey, 5-9 June 2000, vol. I, pp. 512-515 .

M=length(w);
M=M/(J+1);
if length(FFTsynthesisfilters)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

%
% Reconstruction of the signal from its
% bandpass components
%

nu=0:1/M:(1-1/M);

G=conj(FFTsynthesisfilters(1,:));
H=conj(FFTsynthesisfilters(2,:));

H=H.*exp(-2*i*pi*nu);

y=w(length(w)+((-M+1):0));
w=w(1:(length(w)-M));
Y=fft(y,M);
for j=J:-1:1
	z=w(length(w)+((-M+1):0));
	w=w(1:(length(w)-M));
	Z=fft(z,M);
	%M=2*M;
	
	H1=H(1:2^(j-1):length(H)); H1=repmat(H1,[1 2^(j-1)]);
	G1=G(1:2^(j-1):length(G)); G1=repmat(G1,[1 2^(j-1)]);

	Y=(G1.*Y+H1.*Z)/2;
%	Y0=G1(1:M/2).*Y+H1(1:M/2).*Z;
%	Y1=G1(M/2+(1:M/2)).*Y+H1(M/2+(1:M/2)).*Z;
%	Y=[Y0 Y1];
end
%case for j=1 ... redundant part
%z=w(length(w)+((-M+1):0));
%w=w(1:(length(w)-M));
%Z=fft(z,M);
%M=2*M;

%H1=H;
%G1=G;

%Y0=G1(1:M/2).*Y+H1(1:M/2).*Z;
%Y1=G1(M/2+(1:M/2)).*Y+H1(M/2+(1:M/2)).*Z;
%Y=[Y0 Y1];

%x=real(ifft(Y,M));
x=ifft(Y,M);
