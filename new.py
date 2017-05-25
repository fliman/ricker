import sys, ctypes
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

#_lib = ctypes.CDLL("/RD_work/ywang/pycapi/libG2.so", mode=1)
sys.path.append('/RD_work/ywang/pycapi/')

#import pycapi
import salmon.create
import salmon.plot
import matplotlib.pyplot as plt

from scipy import signal




sr = 2.0

n = 4001

t = (n-1)*sr / 1000

print "t ",t

t0 = 1.0

fmax = 500/sr

fcut = fmax * 0.8

fmin = 0


wmin = 0 * 2.0 * np.pi
#wmin =  -np.pi

dw = 2.0*np.pi/t

df = 1.0 / t

nf= int(fmax /df) 
frq = np.linspace(0, fmax, int(nf))

tim = np.linspace(0,t, n)

fpeak = fmax / 7.0 *2.0*np.pi

fpeak = 100

wpeak = fpeak *  2.0 * np.pi
wcut = 200 * 2.0 * np.pi

fq = []
for i in range(n/2+1):
  w = wmin + i*dw
  f = fmin + i*df
  print w 
#  but = 1.0/(1+(w/wcut) * (w/wcut) * (w/wcut) * (w/wcut)*(w/wcut) * (w/wcut) * (w/wcut) * (w/wcut))
  ex = 2.0  * np.exp(w*-1j*t0) * np.exp(-w*w/wpeak/wpeak) * w * w / np.sqrt(np.pi) / wpeak / wpeak / wpeak 
#  if(f >= fcut): 
#    ang = np.pi * 0.5 * (f - fcut) / (fmax - fcut)
#    ex = ex * np.cos(ang) * np.cos(ang)
#    if(f >= fmax * 0.98):
#      ex = 0.0
  fq.append(ex)
  #fq.append(ex*np.sqrt(1.0))


tt=salmon.create.ricker(n, sr, 100, t0, 1.0)

ww = fft.fft(tt)

print ww


pad=[0 for i in range(n)]

fhalf = np.array(fq[0:n/2+1])
shalf =  np.conj(fq[n/2+1:0:-1])

newfq=np.concatenate((fhalf, shalf))
print newfq
newfq = np.array(newfq)
print len(newfq)
#print len(fq[nf:0:-1]), len(fq[0:nf]), nf
kk=fft.ifft(newfq) *13*n/4

print kk

plt.figure()
plt.plot(frq, 20 * np.log10(np.abs(newfq[range(nf)])))
#plt.plot(frq, 20 * np.log10(np.abs(ww[range(nf)])))
#plt.semilogx(frq, 20 * np.log10(np.abs(newfq[range(nf)])))
#plt.plot(frq, fq[nf:0:-1])
plt.figure()
plt.plot(tim[0:len(kk)], kk)



 
tim1=[]
for i in range(n):
	
  dt = i * sr / 1000.0
  bb =  np.pi * np.pi * fpeak * fpeak * (dt-t0) * (dt - t0)
  ex = (1.0 - 2 * bb) * np.exp(-bb) + 3.0
#  if(f >= fcut): 
#    ang = np.pi * 0.5 * (f - fcut) / (fmax - fcut)
#    ex = ex * np.cos(ang) * np.cos(ang)
#    if(f >= fmax * 0.98):
#      ex = 0.0
  tim1.append(ex)

#plt.plot(tim[0:len(kk)], tim1[0:len(kk)])

plt.show()
