load data.mat

%plot Busse balloon in wave number
figure(1)
hold on
plot(a1,2*pi./WN1,'b')
plot(a2,2*pi./WN2,'b')
plot(a3,2*pi./WN3,'b')
plot(a4,2*pi./WN4,'b')
plot(a5,2*pi./WN5,'b')
plot(a6,2*pi./WN6,'b')
hold off

%plot Busse balloon in pulse number
figure(2)
hold on
plot(a1,500./WN1,'b')
plot(a2,500./WN2,'b')
plot(a3,500./WN3,'b')
plot(a4,500./WN4,'b')
plot(a5,500./WN5,'b')
plot(a6,500./WN6,'b')
hold off
