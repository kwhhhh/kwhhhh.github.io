close all
clc
clear

set_degree = 5 / 111.0;  %transform the distance to degree.
load('../download/sample.mat');
center=nmfms(1,1:2);
ard=nmfms(:,3:4);
ctlon = center(:,1); ctlat = center(:,2);
ardlon = ard(:,1); ardlat = ard(:,2);
count=0;
for i = 1:size(ard,1)
    latdum = (ardlat(i)>(ctlat-set_degree)&ardlat(i)<(ctlat+set_degree));
    londum = (ardlon(i)>(ctlon-set_degree/cos(ardlat(i)*pi/180))&ardlon(i)<(ctlon+set_degree/cos(ardlat(i)*pi/180)));
    index1 = find((londum)==1); index2 = find((latdum)==1); [index,i1,i2]  = intersect(index1,index2);
    if ~isempty(index)
        count=count+1;
    end
end