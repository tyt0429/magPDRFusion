function [mag_az]=magnetheading(north_mag,sublocsfin,nnindex,location,coordinate,orientationmap)

clear PQ NNPOINT
PQ=[location(1,size(location,2)) location(2,size(location,2))];
NNPOINT=dsearchn(coordinate',PQ);
mag_az=north_mag(sublocsfin(size(sublocsfin,1)),1)-(orientationmap(NNPOINT))+90;