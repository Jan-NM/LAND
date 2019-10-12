# Some functionality for LAND
```
Orte = simCluster('default', 2);
cell01 = ClusterAnalysis(Orte);
cell01.scatterPlot
cell01.DBSCAN(30,6,0,2000,1); % radius, minPoints, isRandom, maxDiameter, showImage
cell01.distanceAnalysis(200,0,1); % maxDistance, isRandom, showPlot



```
