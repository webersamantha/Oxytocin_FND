function c=FFTpostfilter(x,FFTinterpfilter);
% 	
% 	Author: Dimitri Van De Ville, October 2008
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	

M=length(x);

if length(FFTinterpfilter)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	c=[];
	return
end

% Fourier transform of the signal
X=fft(x,M);

% Prefilter
X=X.*FFTinterpfilter;

% Inverse Fourier transform
c=ifft(X);
%c=real(ifft(X));
