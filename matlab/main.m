% The following script tests the areaPolygon function. Points A-E are
% collected into an area called 'Pts', which is then sent to areaPolygon.

A = [0 0 0];
B = [6 0 0];
C = [6 3 0];
D = [3 5 0];
E = [0 3 0];

Pts = [A;B;C;D;E];

area = areaPolygon(Pts);