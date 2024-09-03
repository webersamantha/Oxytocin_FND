function PC=jPartialCorrWhittaker(CM)
% computes partial correlation using Whittaker's technique
% na?ve implementation

CMi=inv(CM);
PC=zeros(size(CMi));
for i=1:size(CMi,1)
    for j=1:size(CMi,2)
        PC(i,j)=-CMi(i,j)/sqrt(CMi(i,i)*CMi(j,j));
    end
end