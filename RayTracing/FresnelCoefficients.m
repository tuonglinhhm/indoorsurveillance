% Fresnel Transmission and Reflection Coefficients
% Created by: Salaheddin Hosseinzadeh
% Created on: 03 /  May / 2016
% Last Revision: 16.05.2019
% Notes: 
% Please read notes blow before asking question
% 1- Fresnel equation in the form blow (in this code) does not require the frequency of
% the transmission. This is due to a simple approximation I believe.
% 2- Second simplification is that the magnetization of the materials in the transmission
% environment as result of the magnetic field is negligible. (This is
% usually true in case of radio propagation as most propagation takes place
% in the air). So with all these simplifications the refractive index (n) simplifies to the permittivity.it's only the permittivity that is important.
% er is the epsilon (E) of the second material (1st is always air in radio propagation) and therefore
% discarded from input. 
% theta is the incident angle
% sin(theta1)/sin(theta2) = n2/n1, this is how transmission angle is
% calculated. This angle is only important for finding the energy that
% reflects or transmits through, and does not impact the trajectory of the
% beam as this angle is compensated on the exit from second media to the first (from n2 to n1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% er0 & er are the epsilon (E) of the first and second material (1st is air)
% theta1 & theta2 are the incident and transmission angles in degree
% sin(theta1)/sin(theta2) = n2/n1, this is how transmission angle is calculated
% n1 & n2 are the reflection index of the materials and if permuability of
% the two medium are equal u1 = u2 then n1 = sqrt(e1)
% [TEReflFac,TETransFac,TMReflFac,TMTransFac] = FresnelCoefficients (e1,e2,theta1)
function [TEReflFac,TETransFac,TMReflFac,TMTransFac] = FresnelCoefficients (er0,er,theta1,demoMode)

% clear all
er0 = 1; % for air or vacuum
% e2 = 2;


f = 1; % frequency in GHz

c = 0; % set to zero to cancel conductivity and metal support
d = 1.60;
erC = c.* f.^d;  % imaginary part of relative perm, function of freq
eta = er - i.*erC; 
% None of the above matters cuz c=0 % unless you want otherwise

er = sqrt(er); % converting permittivity to refractive index

theta1 = 0:90;
theta = theta1;
theta1 = pi*theta1./180; % theta1 is in radian for ease of use in MATLAB

thetat = asin(sin(theta1).*sqrt(1/er.^2));


TEReflFac = abs((er0.*cos(theta1) - er.*cos(thetat))./(er0.*cos(theta1)+ er.*cos(thetat))).^2;
TETransFac = 1- TEReflFac;

TMReflFac = abs((er.*cos(theta1) - er0.*cos(thetat))./(er.*cos(theta1)+ er0.*cos(thetat))).^2;

TMTransFac = 1 - TMReflFac;


% Converting to power


if demoMode == 1
    figure('Name',['Fresnel Coeffs for Er = ',num2str(e2)])
    plot(theta,TEReflFac);
    hold on
    plot(theta,TETransFac,'r');

    % figure
    plot(theta,TMReflFac,'LineStyle','--')
    hold on
    plot(theta,TMTransFac,'LineStyle','--','Color','red')


    legend('TE Ref','TE Trans','TM Ref','TM Trans')
end




%{
n1 = sqrt(e1);
n2 = sqrt(e2);

% demo = 0;

 
% theta1 = 0:90;
theta = theta1;
theta1 = pi*theta1./180;

theta2 = asin(n1./n2.*sin(theta1));


TEReflFac = abs((n1.*cos(theta1) - n2.*sqrt(1 - (n1./n2 .* sin(theta1)).^2))./...
    (n1.*cos(theta1) + n2.*sqrt(1 - (n1./n2 .* sin(theta1)).^2))); % S Polarization

TETransFac = 1 - TEReflFac;


TMReflFac = abs((n1.*sqrt(1 - (n1/n2 .* sin(theta1)).^2) - n2.*cos(theta1))./...
    (n1.*(sqrt(1 - (n1/n2.*sin(theta1)).^2)) + n2.*cos(theta1))); % P Polarization

TMTransFac = (2.*n1.*cos(theta1))./(n1.*cos(theta2) + n2.*cos(theta1));

if demoMode == 1
    figure('Name',['Fresnel Coeffs for Er = ',num2str(e2)])
    plot(theta,TEReflFac);
    hold on
    plot(theta,TETransFac,'r');

    % figure
    plot(theta,TMReflFac,'LineStyle','--')
    hold on
    plot(theta,TMTransFac,'LineStyle','--','Color','red')


    legend('TE Ref','TE Trans','TM Ref','TM Trans')
end
%}

































