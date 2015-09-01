% This demo does not work yet.  Please don't use

for i=1:4;
[data(i),truth] = nirs.testing.simData;
end
j = nirs.modules.OpticalDensity();
j = nirs.modules.Resample(j);
j = nirs.modules.BeerLambertLaw(j);
data2=j.run(data);

j = nirs.modules.Connectivity();
ConnStats=j.run(data2);

j=nirs.modules.MixedEffectsConnectivity();
 