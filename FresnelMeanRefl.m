function R_mean = FresnelMeanRefl(n12,thetaN)
% calculates mean reflectance based on Fresnel coefficients.

if nargin == 1
    thetaN = 1000;
end

theta = linspace(0,pi/2,thetaN);    % theta domain[0,pi/2]
Rs_theta = abs(...
    (n12*cos(theta) - sqrt(1-(n12*sin(theta)).^2))...
    ./(n12*cos(theta) + sqrt(1-(n12*sin(theta)).^2))...
    ).^2;   % reflectance of s-pol (power ratio)

Rp_theta = abs(...
    (n12*sqrt(1-(n12*sin(theta)).^2) - cos(theta))...
    ./(n12*sqrt(1-(n12*sin(theta)).^2) + cos(theta))...
    ).^2;   % reflectance of p-pol (power ratio)

Rnon_theta = 1/2 * (Rs_theta + Rp_theta);    % non-polarization refelctance

%% weight factor calculatoin considering sphere shape.

% Now we are considering the mean reflectance in the spherical scattering
% cavity. According to the Lambertian cosine law, the incident beam is
% identically  distributed on the cavity surface. Therfore, the weight
% factor will proportional to the area of sphere shell that provides the
% same incident angle theta, which is sin(2*theta).

R_mean= sum(Rnon_theta.*sin(2*theta)) / sum(sin(2*theta));    

