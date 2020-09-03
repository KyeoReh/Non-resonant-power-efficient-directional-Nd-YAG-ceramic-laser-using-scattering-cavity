function meanTransporLength = calcMeanTransportLength(theta_fiber,N)
    % unit of meanTransporLength = radius of cavity
    
    phiVec = linspace(0,2*pi,2*N);
    % thetaVec = acos(linspace(cos(theta_fiber),-1,N));
    thetaVec = linspace(theta_fiber,pi,N);
    theta0Vec = thetaVec;
    [theta0Grid, thetaGrid, phiGrid] = meshgrid(theta0Vec,thetaVec,phiVec);    
    
    expr = sin(thetaGrid).*sin(theta0Grid)...
        .*sqrt(1-sin(thetaGrid).*sin(theta0Grid).*cos(phiGrid)-cos(thetaGrid).*cos(theta0Grid));
    clearvars phiGrid
    % Fry, E.S., et al., Integrating cavities: temporal response. Applied Optics, 2006. 45(36): p. 9053-9065.
    
    normFactor = sin(thetaGrid).*sin(theta0Grid);
    clearvars thetaGrid theta0Grid    
    meanTransporLength = real(sqrt(2)*sum(expr(:))/sum(normFactor(:)));    
end