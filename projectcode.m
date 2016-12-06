M = 1e6;
k = 1e3;
str=['T12:';'T21:';'T13:';'T43:';'T2 :';'T3 :';'T4 :';'G1 :';'G2 :';'M1 :';...
    'Z12:';'Z13:';'Z24:';'Z34:'];
% old S  T12    T21   T13   T43   T2     T3     T4      G1   G2    M1
Sb_old = [200*M, 60*M, 50*M, 50*M, 100*M, 100*M, 100*M, 60*M, 50*M, 50*M];
Sb = 100*M;
%   old V  T12     T21    T13   T43   T2     T3      T4     G1       G2     M1
Vb_old = [13.8*k, 200*k, 69*k, 66*k, 6.9*k, 13.8*k, 34*k, 12.8*k, 13.8*k, 2.4*k];

Vb = [13.8*k, 200*k, 69*k, 69*k, 6.9*k, 13.8*k, 34.5*k, 13.8*k, 6.9*k, 2.4353*k];

Zpu_old = [0.4, 0.06, 0.05, 0.082, 0.1, 0.2, 0.2, 0.35, 0.2, 0.1];
Zb_old =  Vb_old.^2 ./ Sb_old;
Z = Zpu_old.*Zb_old;

Zb = Vb.^2 ./ Sb;
Zpu = Zpu_old.*Zb_old./Zb;
Zpu_line = round([120, 4.76, 2.38, 7.14]./Zb([2,3,7,4]),2);%[200*k, 69*k, 34.5*k, 69*k];
Zpu = [Zpu, Zpu_line];
for i = 1:14
    fprintf('%s %f\n',str(i,:),Zpu(i));
end

