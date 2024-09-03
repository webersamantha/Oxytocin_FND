function w=rFFTwaveletanalysis(x,FFTanalysisfilters,J);
% rFFTWAVELETANALYSIS FFT-based implementation of the redundant wavelet transform.
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

M=length(x);
%if M~=2^round(log(M)/log(2))
%	disp(' ')
%	disp('The size of the input signal must be a power of two!')
%	disp(' ')
%	w=[];
%	return
%end

if length(FFTanalysisfilters)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

% Fourier transform of the signal
X=fft(x,M);

nu=0:1/M:(1-1/M);

G=FFTanalysisfilters(1,:);
H=FFTanalysisfilters(2,:);

H=H.*exp(i*2*pi*nu);

w=[];

Y=G.*X;
Z=H.*X;
%         Y=1/2*(Y(1:M/2)+Y(M/2+(1:M/2)));
%         Z=1/2*(Z(1:M/2)+Z(M/2+(1:M/2)));
z=ifft(Z,M);
w=[w z];

%         M=M/2;
X=Y;

%         G=G(1:2:length(G));
%         H=H(1:2:length(H));

% if J >= 2 %if J<2 the for-loop is not used anyway
    for j=2:J
        %
        % Computation of the outputs y and z
        %
	G=[G(1:2:length(G)) G(1:2:length(G))];
	H=[H(1:2:length(H)) H(1:2:length(H))];

        Y=G.*X;
        Z=H.*X;
        %Y=1/2*(Y(1:M/2)+Y(M/2+(1:M/2)));
        %Z=1/2*(Z(1:M/2)+Z(M/2+(1:M/2)));
        z=ifft(Z);
        w=[w z];
        
        X=Y;
        
%        G=G(1:2:length(G));
%        H=H(1:2:length(H));
    end
% end
w=real([w ifft(X,M)]);

