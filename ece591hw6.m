clear; clc;

%% Problem 1 
% define functions 
f1 = @(theta1, theta2) ( 1-2*cos(3*theta1)+2*cos(3*theta2) );
f2 = @(theta1, theta2) ( 1-2*cos(5*theta1)+2*cos(5*theta2) );

% J = [df1/dt1, df1/dt2; df2/dt1, df2/dt2]
J = @(theta1,theta2) [6*sin(3*theta1), -6*sin(3*theta2); ...
                    10*sin(5*theta1), -10*sin(5*theta2)];

% new thetas = [t1;t2] - J^-1 * [f1;f2]
xNew = @(t1,t2) ( [t1; t2] - J(t1,t2)\[f1(t1,t2);f2(t1,t2)] );

threshold = 0.01;
% guess initial
t1_0 = 0.1; 
t1=t1_0;

t2_0 = 0.5; 
t2=t2_0;
dx = 1;
iterations=0;
allowedIter = 100000;
while(abs(dx(1)) >= threshold || abs(dx(2)) >= threshold)
    % calculate new theta
    tnew = xNew(t1, t2);
    t1_1 = tnew(1);
    V2_1 = tnew(2);
    % Calculate error
    dx = [t1_1-t1; V2_1-t2];
    % update theta1 and theta2
    t1=t1_1;
    t2=V2_1;
    % increment iterations catch
    iterations=iterations+1;
    if(iterations > allowedIter)
        warning(['Divergent case suspected. Pick new initial guesses,'...
            'change threshold, or increase number of allowed iterations.']);
        fprintf('Theta1 Error = %f\nTheta2 Error = %f\nNum Iterations = %i\n',...
            dx(1), dx(2), iterations);
        break;
    end
end
disp('Problem 1 Output:');
fprintf('  Theta 1 Error = %f\n  Theta 2 Error = %f\n  Num Iterations = %i\n',...
            dx(1), dx(2), iterations);
disp('Problem 1, part b:');
fprintf('  Theta 1 = %f degrees\n  Theta 2 = %f degrees\n\n',...
            t1*180/pi, t2*180/pi);
%% Problem 1 Code Output:
% Problem 1 Output:
%   Theta 1 Error = 0.001008
%   Theta 2 Error = 0.001228
%   Num Iterations = 5
% Problem 1, part b:
%   Theta 1 = 23.645141 degrees
%   Theta 2 = 33.327816 degrees
%% Problem 2%% define functions and variables
P2 = -0.1358;
Q2 = -0.2291;
Z12 = 0.01+1i*0.2; 
Y12 = Z12^-1;
% admittance matrix:
% [G11+B11  -G12-B12]
% [-G21-B21  G22+B22]
% G11=G12=G21=G22
% B11=B12=B21=B22
G12 = real(Y12);
B12 = imag(Y12);
f1 = @(theta2, V2) ( -P2-V2*G12*cos(theta2)-V2*B12*sin(theta2)+V2^2*G12 );
f2 = @(theta2, V2) ( -Q2-V2*G12*sin(theta2)+V2*B12*cos(theta2)-V2^2*B12 );

% J = [df1/dt1, df1/dV2; df2/dt1, df2/dV2]
J = @(t2, V2) [V2*G12*sin(t2)-V2*B12*cos(t2), -G12*cos(t2)-B12*sin(t2)+2*V2*G12;...
               +V2*G12*cos(t2)-V2*B12*sin(t2), -G12*sin(t2)+B12*cos(t2)-2*V2*B12];

% new thetas = [t1;t2] - J^-1 * [f1;f2]
xNew = @(t2, V2) ( [t2; V2] - J(t2,V2)\[f1(t2,V2);f2(t2,V2)] );

V1=1; t1=0;
threshold = 0.01;
% guess initial
t2_0 = 0; 
t2=t2_0;

V2_0 = 1; 
V2=V2_0;
dx = 1;
iterations=0;
allowedIter = 100000;
while(abs(dx(1)) >= threshold || abs(dx(2)) >= threshold)
    % calculate new theta
    tnew = xNew(t2, V2);
    t2_1 = tnew(1);
    V2_1 = tnew(2);
    % Calculate error
    dx = [t2_1-t2; V2_1-V2];
    % update theta1 and theta2
    t2=t2_1;
    V2=V2_1;
    % increment iterations catch
    iterations=iterations+1;
    if(iterations > allowedIter)
        warning(['Divergent case suspected. Pick new initial guesses,'...
            'change threshold, or increase number of allowed iterations.']);
        fprintf('Theta2 Error = %f\nV2 Error = %f\nNum Iterations = %i\n\n',...
            dx(1), dx(2), iterations);
        break;
    end
end
disp('Problem 2 Output:')
fprintf('  Theta_2 Error = %f\n  V_2 Error = %f\n  Num Iterations = %i\n',...
            dx(1), dx(2), iterations);
disp('Problem 2, part d:');
fprintf('  Theta_2 = %f degrees\n  V_2 = %f volts\n',...
            t2*180/pi, V2);
disp('Problem 2, part e:');
fprintf('  Experimentally derived values:\n  P2 = %f W\n  Q2 = %f VA\n\n',...
            f1(t2,V2)+P2,f2(t2,V2)+Q2);
% I12=((V1-V2)/Z12);
% Sline=I12*conj(I12*Z12);
% S1 = V1*conj(-I12);
% fprintf('Total Power = %f W\nTotal Reactive Power = %f VA\n\n',...
%     real(S1)+real(Sline)+P2,imag(S1)+imag(Sline)+Q2);

%% Problem 2 Code Output:
% Problem 2 Output:
%   Theta_2 Error = -0.000904
%   V_2 Error = -0.006884
%   Num Iterations = 2
% Problem 2, part d:
%   Theta_2 = -1.507062 degrees
%   V_2 = 0.956535 volts
% Problem 2, part e:
%   Experimentally derived values:
%   P2 = -0.135757 W
%   Q2 = -0.199439 VA
%% Problem 3
%define functions and variables
P1 = -0.8913;
Q1 = -0.1434;

V2 = 1;
theta2 = 0;

V3 = 1.02;
P3 = -0.8184;
% admittance matrix
B11 = 15; B12 = -5; B13 = -10;
B21 = -5; B22 = 15; B23 = -10;
B31 = -10; B32 = -10; B33 = 20;
f1 = @(t1,t3,V1) ...
    ( -P1 + V1*B12*sin(t1) + V1*V3*B13*sin(t1-t3) );  % P1 eq
f2 = @(t1,t3,V1) ...
    ( -P3 + V3*V1*B31*sin(t3-t1) + V3*B32*sin(t3));  % P3 eq
f3 = @(t1,t3,V1) ...
    ( -Q1 - V1^2*B11 - V1*B12*cos(t1) - V1*V3*B13*cos(t1-t3) ); % Q1 eq

% J = [df1/dt1, df1/dt3 df1/dV1 
%      df2/dt1, df2/dt3 df2/dV1
%      df3/dt1, df3/dt3 df3/dV1]
J = @(t1,t3,V1)...
 [V1*B12*cos(t1)+V1*V3*B13*cos(t1-t3), -V1*V3*B13*cos(t1-t3), B12*sin(t1)+V3*B13*sin(t1-t3); ...
 -V3*V1*B31*cos(t3-t1), V3*V1*B31*cos(t3-t1)+V3*B32*cos(t3), V3*B31*sin(t3-t1);...
 V1*B12*sin(t1)+V1*V3*B13*sin(t1-t3), -V1*V3*B13*sin(t1-t3), -2*V1*B11-B12*cos(t1)-V3*B13*cos(t1-t3)];

% new thetas = [t1;t2] - J^-1 * [f1;f2]
xNew = @(t1, t3, V1) ...
 ( [t1; t3; V1] - J(t1, t3, V1)\[f1(t1, t3, V1);...
                                 f2(t1, t3, V1);...
                                 f3(t1, t3, V1)] );

threshold = 0.01;
% guess initial
t1_0 = 0; 
t1=t1_0;

t3_0 = 0; 
t3=t3_0;

V1_0 = .75; 
V1=V1_0;

dx = [1,1,1];

iterations=0;
allowedIter = 100000;

while(abs(dx(1)) >= threshold ...
        || abs(dx(2)) >= threshold ...
        || abs(dx(3)) >= threshold )
    % calculate new theta
    tnew = xNew(t1, t3, V1);
    t1_1 = tnew(1);
    t3_1 = tnew(2);
    V1_1 = tnew(3);
    % Calculate error
    dx = [t1_1-t1; t3_1-t3;V1_1-V1];
    % update theta1 and theta2
    t1=t1_1;
    t3=t3_1;
    V1=V1_1;
    % increment iterations catch
    iterations=iterations+1;
    if(iterations > allowedIter)
        warning(['Divergent case suspected. Pick new initial guesses,'...
            'change threshold, or increase number of allowed iterations.']);
        fprintf(['Theta1 Error = %f\nV1 Error = %f\nTheta3 Error = %f\n'...
            'Num Iterations = %i\n\n'],...
            dx(1), dx(3), dx(2), iterations);
        break;
    end
end
disp('Problem 3 Output:');
fprintf(['  Theta_1 Error = %f\n  V_1 Error = %f\n  Theta_3 Error = %f\n'...
    'Num Iterations = %i\n\n'],...
            dx(1), dx(3), dx(2),iterations);
P2 = V2*V1*B21*sin(-t1) + V2*V3*B23*sin(-t3);
Q2 = -V2*V1*B21*cos(-t1) - V2^2*B22 - V2*V3*B23*cos(-t3);
Q3 = -V3*V1*B31*cos(t3-t1) - V3*V2*B32*cos(t3) - V3^2*B33;
disp('Problem 3, part e:');
fprintf(['  Theta_1 = %f degrees\n  V_1 = %f volts\nP_2 = %f W\n  '...
        'Q_2 = %f VA\n  Theta_3 = %f degrees\n  Q_3 = %f VA\n\n'],...
            t1*180/pi, V1, P2, Q2, t3*180/pi, Q3);
 % Q = I^2*X      
%  I12 = (V1-V2)/(1i*0.2); I13 = (V1-V3)/(1i*0.1); I23 = (V2-V3)/(1i*0.1);
% %  Q_line12=((V1-V2)/0.2)^2*0.2;
% %  Q_line13=((V1-V3)/0.1)^2*0.1;
% %  Q_line23=((V2-V3)/0.1)^2*0.1;
% S_line12 = I12*conj(I12*(1i*0.2));
% S_line13 = I13*conj(I13*(1i*0.1));
% S_line23 = I23*conj(I23*(1i*0.1));
% fprintf('Total Real Power = %f W\nTotal Reactive Power = %f VA\n\n',...
%     P1+real(S_line12)+real(S_line13)+real(S_line23)+P2+P3,...
%     Q1+imag(S_line12)+imag(S_line13)+imag(S_line23)+Q2+Q3);
%% Problem 3 Code Output:
% Problem 3 Output:
%   Theta_1 Error = -0.000002
%   V_1 Error = -0.000317
%   Theta_3 Error = 0.000009
%   Num Iterations = 4
% 
% Problem 3, part e:
%   Theta_1 = 7.307992 degrees
%   V_1 = 1.019815 volts
% P_2 = 1.709700 W
%   Q_2 = 0.202312 VA
%   Theta_3 = 5.971152 degrees
%   Q_3 = -0.264061 VA