function [pdc] =test_muthu(N,D,p,a2)
s1=size(a2);
q=s1(2)/p;

Phat= size(a2, 2)/size(a2, 1);
results.Ahat= a2;
results.Phat= Phat;

Apoly= reshape([eye(D) -results.Ahat ], D, D, Phat+1);
Aspec= fftshift(fft(Apoly, N, 3), 3);
Aspec= Aspec(:, :, N/2+1:end);

PDCnorma= sum(Aspec.*conj(Aspec), 1);

d1= 1;
for d1= 1:D
    results.PDC(d1, :, :)= abs(Aspec(d1, :, :))./sqrt(PDCnorma(1, :, :));

end
pdc=results.PDC;
end


