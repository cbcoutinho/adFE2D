function [ area ] = areaPolygon( Pts )
%AREAPOLYGON Determines the area of a polygon in 2D
%   The area of a polygon is calculated by decomposing the polygon into a
%   triangles. The area of each triangle is calculated and combined into
%   the total area.

[m,~] = size(Pts);

area = 0;

for ii = 1:m-2
    x1 = Pts(1,:);
    x2 = Pts(ii+1,:);
    x3 = Pts(ii+2,:);
    
    area = area + areaTriangle(x1,x2,x3);

end

end

