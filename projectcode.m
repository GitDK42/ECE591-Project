%% Setup

M = 1e6;
k = 1e3;
TijStr = ['T12:';'T21:';'T13:';'T43:'];
TkStr = ['T2 :';'T3 :';'T4 :'];
XkStr = ['G1 :';'G2 :';'M1 :'];
ZijStr =['Z12:';'Z13:';'Z24:';'Z34:'];
str=[TijStr; TkStr; XkStr; ZijStr];

%% a.) Change of base 
%T12 = 200*M; T21=60*M; T13=60*M; T43*M;
% Tij         T12    T21   T13   T43
TijSb_old = [200*M, 60*M, 50*M, 50*M];
TijVb_old = [13.8*k, 200*k, 69*k, 66*k];
TijVb = [13.8*k, 200*k, 69*k, 69*k]; % new
Tijpu_old = [0.4, 0.06, 0.05, 0.082];
Tij = (1:4);

% Tk         T2     T3     T4
TkSb_old = [100*M, 100*M, 100*M];
TkVb_old = [6.9*k, 13.8*k, 34*k];
TkVb = [6.9*k, 13.8*k, 34.5*k];      % new
Tkpu_old  = [0.1, 0.2, 0.2];
Tk = (5:7);
% Xk         G1    G2    M1
XkSb_old = [60*M, 50*M, 50*M];
XkVb_old = [12.8*k, 13.8*k, 2.4*k];
XkVb = [13.8*k, 6.9*k, 2.4353*k];    % new
Xkpu_old  = [0.35, 0.2, 0.1];
Xk = (8:10);
% Old S vector
Sb_old = [TijSb_old, TkSb_old, XkSb_old];
Sb = 100*M; % new S bases
% Old V vector 
Vb_old = [TijVb_old, TkVb_old, XkVb_old];
Vb = [TijVb, TkVb, XkVb]; % new V bases
% Old Zpu vector
Zpu_old = [Tijpu_old, Tkpu_old, Xkpu_old];

Zb_old =  Vb_old.^2 ./ Sb_old;
Z = Zpu_old.*Zb_old;

Zb = Vb.^2 ./ Sb;
Zpu = Zpu_old.*Zb_old./Zb;
%                 Z12   Z13   Z24   Z34
Zpu_line = round([120, 4.76, 2.38, 7.14]./Zb([2,3,7,4]),2);%[200*k, 69*k, 34.5*k, 69*k];
Zpu = [Zpu, Zpu_line];
Zij = (11:14);
for i = 1:14
    fprintf('%s %f\n',str(i,:),Zpu(i));
end

%% b.) Ybus

% [ 1/(T12+Z12+T21+T13+Z13)  -1/(T12+Z12+T21)  -1/(T13+Z13)          0     ]    
% [  -1/(T12+Z12+T21)     1/(T12+Z12+T21+Z24)        0            -1/(Z24) ] 
% [     -1/(T13+Z13)              0       1/(Z31+T31+Z34+T34)  -1/(Z34+T34)]
% [            0               -1/(Z24)      -1/(Z34+T34)   1/(Z24+Z34+T43)]
Ztop = sum(Zpu([Tij([1,2]),Zij(1)]));
Zleft = sum(Zpu([Tij(3),Zij(2)]));
Zright = sum(Zpu(Zij(3)));
Zbottom = sum(Zpu([Tij(4),Zij(4)]));
Y1 = [-1/(Ztop)-1/(Zleft),...
      1/(Ztop), ...
      1/(Zleft),...
     0];
Y2 = [1/(Ztop),...
      -1/(Ztop)-1/(Zright), ...
       0,...
      1/(Zright)];
Y3 = [1/(Zleft),...
       0, ...
      -1/(Zleft)-1/(Zbottom),...
      1/(Zbottom)];
Y4 = [ 0,...
      1/(Zright), ...
      1/(Zbottom),...
      -1/(Zbottom)-1/(Zright)];
Y = 1i*[Y1;Y2;Y3;Y4] % Ybus is all Bs as Z values are purely imaginary.
G = real(Y);
B = imag(Y);
%% c.) Newton-Raphson Method
V1 = 1; V2 = 0.95; V3 = 1.02;
d1 = 0; 
P2 = 0.25; P3 = 0.35; P4 = -0.72;
Q4 = -0.19;

f1 = @(d2,d3,d4,V4) ...
    ( -P2 + V2*B(2,1)*sin(d2-d1)+V2*V4*B(2,4)*sin(d2-d4) );  % P2 eq
f2 = @(d2,d3,d4,V4) ...
    ( -P3 + V3*B(3,1)*sin(d3-d1)+V3*V4*B(3,4)*sin(d3-d4) );  % P3 eq
f3 = @(d2,d3,d4,V4) ...
    ( -P4 + V4*V2*B(4,2)*sin(d4-d2)+V4*V3*B(4,3)*sin(d4-d3) );  % P4 eq
f4 = @(d2,d3,d4,V4) ...
    ( -Q4 - B(4,4)*V4^2 - V4*V2*B(4,2)*cos(d4-d2) - V4*V3*B(4,3)*cos(d4-d3) ); % Q4 eq

J = @(d2,d3,d4,V4)...
 [V2*B(2,1)*cos(d2-d1)+V2*V4*B(2,4)*cos(d2-d4), 0, -V2*V4*B(2,4)*cos(d2-d4), V2*B(2,4)*sin(d2-d4); ...
 0, V3*B(3,1)*cos(d3-d1)+V3*V4*B(3,4)*cos(d3-d4), -V3*V4*B(3,4)*cos(d3-d4), V3*B(3,4)*sin(d3-d4);...
 -V4*V2*B(4,2)*cos(d4-d2), -V4*V3*B(4,3)*cos(d4-d3), V4*V2*B(4,2)*cos(d4-d2)+V4*V3*B(4,3)*cos(d4-d3), V2*B(4,2)*sin(d4-d2)+V3*B(4,3)*sin(d4-d3);...
 -V4*V2*B(4,2)*sin(d4-d2), -V4*V3*B(4,3)*sin(d4-d3), V4*V2*B(4,2)*sin(d4-d2)+V4*V3*B(4,3)*sin(d4-d3), -2*V4*B(4,4)-V2*B(4,2)*cos(d4-d2)-V3*B(4,3)*cos(d4-d3)];

% new values = [d2;d3;d4;V4] - J^-1 * [f1;f2;f3;f4]
xNew = @(d2,d3,d4,V4) ...
 ( [d2;d3;d4;V4] - J(d2,d3,d4,V4)\[f1(d2,d3,d4,V4);...
                                   f2(d2,d3,d4,V4);...
                                   f3(d2,d3,d4,V4);...
                                   f4(d2,d3,d4,V4)] );
                                 
% initial guess
d2_0 = 0; d3_0 = 0; d4_0 = 0;
V4_0 = 1;
d2=d2_0;d3=d3_0;d4=d4_0;
V4=V4_0;

dx = [1,1,1,1];

threshold = 0.001;
iterations=0;
allowedIter = 100000;

while(abs(dx(1)) >= threshold ...
        || abs(dx(2)) >= threshold ...
        || abs(dx(3)) >= threshold ...
        || abs(dx(4)) >= threshold)
    % calculate new theta
    tnew = xNew(d2,d3,d4,V4);
    d2_1 = tnew(1);
    d3_1 = tnew(2);
    d4_1 = tnew(3);
    V4_1 = tnew(4);
    % Calculate error
    dx = [d2_1-d2; d3_1-d3; d4_1-d4; V4_1-V4];
    % update theta1 and theta2
    d2=d2_1;
    d3=d3_1;
    d4=d4_1;
    V4=V4_1;
    % increment iterations catch
    iterations=iterations+1;
    if(iterations > allowedIter)
        warning(['Divergent case suspected. Pick new initial guesses,'...
            'change threshold, or increase number of allowed iterations.']);
%         fprintf(['Theta1 Error = %f\nV1 Error = %f\nTheta3 Error = %f\n'...
%             'Num Iterations = %i\n\n'],...
%             dx(1), dx(3), dx(2), iterations);
        break;
    end
end
rad2deg=@(x) 180*x/pi;
disp('Error Output:');
fprintf(['  Delta_2 Error = %f radians\n  Delta_3 Error = %f radians\n  Delta_4 Error = %f radians\n'...
    '  V4 Error = %f Volts\n'...
    'Num Iterations = %i\n\n'],...
    dx(1), dx(2), dx(3), dx(4),iterations);

% Bus Powers
P1 =  V2*B(1,2)*sin(d1-d2) + V3*B(1,3)*sin(d1-d3);
Q1 = -V1^2*B(1,1) - V2*B(1,2)*cos(d1-d2) - V3*B(1,3)*cos(d1-d3);
P2 =  V2*B(2,1)*sin(d2-d1)+V2*V4*B(2,4)*sin(d2-d4);
Q2 = -V2^2*B(2,2) - V2*B(2,1)*cos(d2-d1) - V2*V4*B(2,4)*cos(d2-d4);
P3 =  V3*B(3,1)*sin(d3-d1)+V3*V4*B(3,4)*sin(d3-d4)
Q3 = -V3^2*B(3,3) - V3*B(3,1)*cos(d3-d1)+V3*V4*B(3,4)*cos(d3-d4);
P4 =  V4*V2*B(4,2)*sin(d4-d2)+V4*V3*B(4,3)*sin(d4-d3);
Q4 = -B(4,4)*V4^2 - V4*V2*B(4,2)*cos(d4-d2) - V4*V3*B(4,3)*cos(d4-d3);

% Line Currents:
Itop = (V1*(cos(d1)+1i*sin(d1))-V2*(cos(d2)+1i*sin(d2)))/Ztop;       %(V1-V2)/Ztop
Ileft = (V1*(cos(d1)+1i*sin(d1))-V3*(cos(d3)+1i*sin(d3)))/Zleft;     %(V1-V3)/Zleft
Iright = (V2*(cos(d2)+1i*sin(d2))-V4*(cos(d4)+1i*sin(d4)))/Zright;   %(V2-V4)/Zright
Ibottom = (V3*(cos(d3)+1i*sin(d3))-V4*(cos(d4)+1i*sin(d4)))/Zbottom; %(V3-V4)/Zbottom

% Line Powers:
Qtop = abs(Itop)^2 * Ztop;
Qleft = abs(Ileft)^2 * Zleft;
Qright = abs(Iright)^2 * Zright;
Qbottom = abs(Ibottom)^2 * Zbottom;

disp('Calculated Values (Per Unit):');
disp('Bus 1:');
fprintf(['  V1 = %+.3f Volts\n  d1 = %+.3f Degrees\n  P1 = %+.3f Watts\n  Q1 = %+.3f VAr\n\n'],...
    V1, rad2deg(d1), P1, Q1);
disp('Bus 1:');
fprintf(['  V2 = %+.3f Volts\n  d2 = %+.3f Degrees\n  P2 = %+.3f Watts\n  Q2 = %+.3f VAr\n\n'],...
    V2, rad2deg(d2), P2, Q2);
disp('Bus 1:');
fprintf(['  V3 = %+.3f Volts\n  d3 = %+.3f Degrees\n  P3 = %+.3f Watts\n  Q3 = %+.3f VAr\n\n'],...
    V3, rad2deg(d3), P3, Q3);
disp('Bus 4:');
fprintf(['  V4 = %+.3f Volts\n  d4 = %+.3f Degrees\n  P4 = %+.3f Watts\n  Q4 = %+.3f VAr\n\n'],...
    V4, rad2deg(d4), P4, Q4);
disp('Total Real Grid Power:');
fprintf('  P_total = %+.10f\n',P1+P2+P3+P4);
disp('Total Reactive Grid Power:');
fprintf('  Q_total = %+.10f\n\n',Q1+Q2+Q3+Q4+Qtop+Qleft+Qright+Qbottom);
% fprintf(['  Delta_2 = %+.3f Degrees\n  Delta_3 = %+.3f Degrees\n  Delta_4 = %+.3f Degrees\n'...
%     '  V4 = %+.3f Volts\n  P1 = %+.3f Watts\n  Q1 = %+.3f VAr\n  Q2 = %+.3f VAr\n  Q3 = %+.3f VAr\n'...
%     '  P4 = %+.3f Watts\n  Q4 = %+.3f VAr\n\n'],...
%     rad2deg(d2), rad2deg(d3), rad2deg(d4), V4, P1, Q1, Q2, Q3, P4, Q4);

%% d.) Eaf_G2 and delta_G2
X2 = Zpu(Xk(2)) + Zpu(Tk(1));  % TG2 + XG2 = 1.7
dG2 = atan( (P2*X2/V2) / (Q2*X2/V2 + V2) ); % relative to Bus 2
e_afG2 = X2*P2/(V2*sin(dG2+d2));
disp('G2 V and Delta relative to slack bus');
fprintf('  dG2 = %+.3f Degrees\n  e_afG2 = %+.3f\n\n',dG2+d2, e_afG2);


