function APE = sectionAPE(x,z,b)


bMean = mean(b,2,'omitnan');
bPrime = b - bMean;
[N2] = gradient(bMean,z);
L = range(x);
dz = mean(diff(z));


APE = 1/2 * dz * sum( mean(bPrime.^2,2,'omitnan')./N2, 'omitnan');





end