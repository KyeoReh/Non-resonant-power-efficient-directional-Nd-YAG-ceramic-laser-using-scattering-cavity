% For given theta_fiber, the mean transport length of light in the spherical caivity is caculated.
% Notice the theta_fiber is  
%        asin(cavity diameter / aperture diameter)      for   0  < theta_fiber <= pi/2
%  and   pi - asin(cavity diameter / aperture diameter) for pi/2 < theta_fiber <= pi

%% set path
cd(fileparts(matlab.desktop.editor.getActiveFilename))

%% main
xN = 10;      % angle sampling number. higher xN will give presice result, but the calculation takes long.
theta_fiber = linspace(0,pi,xN+1);
theta_fiber = theta_fiber(1:end-1).';
meanTransporLength = zeros(xN,1);

for kk = 1:1:xN
    meanTransporLength(kk) = calcMeanTransportLength(theta_fiber(kk),500);
    plot(theta_fiber(1:kk),meanTransporLength(1:kk))
    pause(0.01);    
end

%% save
save('meanTransportLength.mat','theta_fiber','meanTransporLength')

