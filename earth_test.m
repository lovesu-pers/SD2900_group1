function [earth, axes] = earth_test(J2000,RE)

load topo topo;    % load data 
% whos('topo','topomap1');

omega_earth = 7.292115855377074e-005; % (rad/sec)  
Go = 0; %1.727564365843028; % (rad) http://www.amsat.org/amsat/articles/g3ruh/106.html 
% GMST = Go + omega_earth*86400*(J2000 + 0.5);
GMST = pi+omega_earth*J2000;



[x,y,z] = sphere(50);          % create a sphere 
earth = surface(RE*x,RE*y,RE*z);            % plot spherical surface

earth.FaceColor = 'texturemap';    % use texture mapping
earth.CData = topo;                % set color data to topographic data
earth.EdgeColor = 'none';          % remove edges
earth.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces
earth.SpecularStrength = 0.4;      % change the strength of the reflected light
colormap('winter')


lim = 2*RE;
Xaxis_ECI= line([0 lim],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis_ECI = line([0 0],[0 lim],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');

Xaxis_ECEF = line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
Yaxis_ECEF = line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-');
rotate(Xaxis_ECI, [0 0 1], 0, [0,0,0]);
rotate(Yaxis_ECI, [0 0 1], 0, [0,0,0]);

clim([-1,1])
% light('Position',[-1 0 1])     % add a light


%%
rotate(earth, [0 0 1], GMST*180/pi, [0,0,0]);
rotate(Xaxis_ECEF, [0 0 1], 180+GMST*180/pi, [0,0,0]);
rotate(Yaxis_ECEF, [0 0 1], 180+GMST*180/pi, [0,0,0]);

axes = [Xaxis_ECI, Yaxis_ECI, Xaxis_ECEF, Yaxis_ECEF];
end