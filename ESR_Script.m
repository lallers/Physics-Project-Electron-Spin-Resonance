clc
small.dc =  [ 1.337 1.353 1.334 1.341; 
             1.507 1.515 1.511 1.511;
             1.678 1.695 1.682 1.681;
             1.848 1.874 1.866 1.818;
             2.028 2.045 2.038 2.030;
             2.198 2.207 2.213 2.196];
small.rf =  [ 75.46 75.34 75.16 75.41;
             85.08 85.03 85.29 85.03;
             95.15 95.29 95.01 95.18;
             105.08 105.05 105.28 105.01;
             115.42 115.14 115.07 115.41;
             125.03 125.05 125.12 125.01]*10^6;
medium.dc = [ 0.558 0.738 0.909 1.089 1.268;
             0.567 0.721 0.910 1.081 1.259;
             0.561 0.729 0.903 1.088 1.257;
             0.562 0.731 0.901 1.081 1.251]';
medium.rf = [30.09 39.99 50.54 59.91 70.02;
             30.57 39.92 50.00 59.98 69.96;
             30.41 40.01 49.92 59.92 69.97;
             30.25 39.96 50.25 59.90 70.01]'*10^6;
big.dc    = [ .250 .301 .344 .403 .464 .507;
              .258 .300 .353 .420 .463 .515;
              .243 .310 .349 .419 .478 .518;
              .249 .300 .341 .420 .463 .511]';
big.rf    = [13.05 16.03 19.18 22.02 25.02 28.07;
             13.06 16.01 19.01 22.02 25.07 28.04;
             13.13 16.09 19.03 22.08 25.21 28.04;
             13.01 16.02 19.10 22.04 25.01 28.01]'*10^6;
         

         
%% SMALL
small.meandc = mean(small.dc,2);
small.stddc  = std(small.dc,0,2);
small.meanrf = mean(small.rf,2);
small.stdrf  = std(small.rf,0,2);


%% Medium
medium.meandc = mean(medium.dc,2);
medium.stddc  = std(medium.dc,0,2);
medium.meanrf = mean(medium.rf,2);
medium.stdrf  = std(medium.rf,0,2);


%% Big
big.meandc = mean(big.dc,2);
big.stddc  = std(big.dc,0,2);
big.meanrf = mean(big.rf,2);
big.stdrf  = std(big.rf,0,2);


%% B Calc
p.R   = 13.5 / 2 * 10^(-2);
p.mu  = 1.256*10^(-6);
p.N   = 320 ;
p.con = (4/5)^(3/2);
p.h   = 6.2582*10^(-16);
p.mub = 5.788 * 10^(-9);


small.B    = (small.meandc * p.mu * p.con *p.N * 1/p.R)/2;
small.Bf    = (small.dc * p.mu * p.con *p.N * 1/p.R)/2;
medium.B   = (medium.meandc * p.mu * p.con *p.N * 1/p.R)/2;
medium.Bf   = (medium.dc * p.mu * p.con *p.N * 1/p.R)/2;
big.B      = (big.meandc * p.mu * p.con *p.N * 1/p.R)/2;
big.Bf      = (big.dc * p.mu * p.con *p.N * 1/p.R)/2;

small.fit  = polyfit(small.B,small.meanrf,1);
small.val  = polyval(small.fit, small.B);

medium.fit  = polyfit(medium.B,medium.meanrf,1);
medium.val  = polyval(medium.fit, medium.B);

big.fit  = polyfit(big.B,big.meanrf,1);
big.val  = polyval(big.fit, big.B);

small.g  = p.h .* small.meanrf./(p.mub .* small.B) ;
medium.g = p.h .* medium.meanrf./(p.mub .* medium.B) ;
big.g    = p.h .* big.meanrf./(p.mub .* big.B) ;

small.g_fit = (small.fit * p.h/p.mub)*10^(-3) - 1 ;
medium.g_fit = (medium.fit * p.h/p.mub)*10^(-3) - 1;
big.g_fit = (big.fit * p.h/p.mub)*10^(-3) - 1;

%%
figure(1)
hold on

plot(small.Bf,small.rf,'x')
plot(medium.Bf,medium.rf,'x')
plot(big.Bf,big.rf,'x')

plot(small.B,small.val)
plot(medium.B,medium.val)
plot(big.B,big.val)

hold off
legend('Small','Medium','Large')
xlabel('Magentic Field Average ( T )');ylabel('Resonant Frequency (Hz)')

%%
B  = [big.B; medium.B; small.B;];
Freq = [big.meanrf; medium.meanrf;  small.meanrf];
Fit = polyfit(B',Freq',1);
Val = polyval(Fit,B);

figure(2)
hold on
scatter(small.B,small.meanrf,'o')
scatter(medium.B,medium.meanrf,'x')
scatter(big.B, big.meanrf,'d')
plot(B,Val)
legend('Small','Medium','Large','fit')
xlabel('Magentic Field Average ( T )');ylabel('Resonant Frequency (Hz)')
%% Resonance for each point
B =@(I) I * p.mu * p.con *p.N * 1/p.R;
small.g = (p.h*small.rf./(p.mub*B(small.dc)))*10^(-3);
smallgmean = mean(small.g,2)
medium.g = (p.h*medium.rf./(p.mub*B(medium.dc)))*10^(-3);
big.g = (p.h*big.rf./(p.mub*B(big.dc)))*10^(-3);

smallgmean = mean(small.g,2);
smallgstd = std(std(small.g,0,2))
medgmean = mean(medium.g,2)
medstd = std(std(medium.g,0,2))
biggmean = mean(big.g,1)
biggstd = std(std(big.g,0,2))


