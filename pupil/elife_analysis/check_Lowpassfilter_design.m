clear
samplerate=500;
t = 0:1/samplerate:10-1/samplerate;
x = sin(2*pi*0.1*t)+sin(2*pi*0.5*t)+sin(2*pi*1*t)+sin(2*pi*2*t)+sin(2*pi*3*t)+sin(2*pi*4*t)+sin(2*pi*5*t)+sin(2*pi*6*t);


y = fft(x);     
f = (0:length(y)-1)*samplerate/length(y);

plot(f,abs(y))
title('Magnitude')

[z,p,k] = butter(1,3.75/(samplerate/2),'low');  % first order filter? why?
[sos,g] = zp2sos(z,p,k);
Hd = dfilt.df2tsos(sos,g);

xx=filter(Hd,x);

yy = fft(xx);     
f = (0:length(yy)-1)*samplerate/length(yy);


subplot(2,2,1)
plot(t,x)
subplot(2,2,2)
plot(f,abs(y))
xlim([0 7])
title('Magnitude')
subplot(2,2,3)
plot(t,xx);
subplot(2,2,4)
plot(f,abs(yy))
xlim([0 7])
title('Magnitude')