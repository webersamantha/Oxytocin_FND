function TCS=my_wavelet(TCS,alpha,tau,nSubbands,subbandN)
% 
% Compute wavelet decomposition of a time series
% 
% alpha:    b-spline param: degree.
% tau:      b-spline param: shift.
%               0: b-spline always symmetrical
%               (alpha+1) / 2: causal
% nSubbands - number of wavelet subbands
% subbandN - subband to return
% 
% typical usage: alpha=3, tau=0, 4th subband for default mode with TR=1.1
%
% v1.0 Sep 2009 Dimitri van de Ville
% length of temporal signal
sz=size(TCS,2);

% number of wavelet subbands
J=nSubbands;

% sanity check: signal length must be at least 2^J long
if sz<2^J
    error(['Signal length must be a minimum of ' num2str(2^J) '.']);
end

% adjust size
%sz=floor(sz0/2^J)*2^J;
%A=A0(1:sz);

% parameters of the wavelet transform
% alpha: degree of wavelet (alpha+1 = order)
type='ortho';     % flavor of wavelet: ortho, bspline, dual

% compute filters
FP=FFTfractinterpfilter(sz,alpha,tau,type);
[FA,FS]=FFTfractsplinefilters(sz,alpha,tau,type);

PA=zeros(size(TCS));

for idx=1:size(TCS,1),

% prefiltering
PA(idx,:)=FFTprefilter(TCS(idx,:),FP);

% perform non-redundant wavelet analysis
%w=FFTwaveletanalysis(PA,FA,J,0);

% perform redundant wavelet analysis
w(idx,:)=rFFTwaveletanalysis(PA(idx,:),FA,J);

% perform non-redundant wavelet synthesis
%PR=FFTwaveletsynthesis(w,FS,J,0);

% perform non-redundant wavelet synthesis
%PR=rFFTwaveletsynthesis(w,FS,J);

% postfiltering
%R=FFTpostfilter(PR,FP);
end;

TCS=w(:,1+(subbandN-1)*sz:subbandN*sz);
