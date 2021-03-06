clear;
N=1000;
k=35;
p=4;

xx = linspace(-1,1,100);
yy = VMFMeanDirDensity(xx,k,p);
yy = yy/sum(yy(:));
[t] = randVMFMeanDir(N, k, p);
h=hist(t, xx);
h = h/sum(h(:));
figure; 
plot(xx, h, 'b');
hold on;
plot(xx,yy, 'r');


N=1000;
k=5;
p=4;

xx = linspace(-1,1,100);
yy = VMFMeanDirDensity(xx,k,p);
yy = yy/sum(yy(:));
[t] = randVMFMeanDir(N, k, p);
h=hist(t, xx);
h = h/sum(h(:));
figure; 
plot(xx, h, 'b');
hold on;
plot(xx,yy, 'r');


N=1000;
k=300;
p=4;

xx = linspace(0.9,1,100);
yy = VMFMeanDirDensity(xx,k,p);
yy = yy/sum(yy(:));
[t] = randVMFMeanDir(N, k, p);
h=hist(t, xx);
h = h/sum(h(:));
figure; 
plot(xx, h, 'b');
hold on;
plot(xx,yy, 'r');