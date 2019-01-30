#!/usr/bin/python

##################################################################
#  Yigit Dallilar  30.01.2013
#
#  FOR CHAOTIC PENDULUM PROBLEM, PROGRAM CALCULATES RES OUTPUT
#  BY USING SOLVE FUNCTIONS. THEN ONE CAN DRAW DESIRED PLOTS 
#  BY OTHER FUNCTIONS.
#
#  **TO DO :
#   (1) ANIMATION OF THE MOVEMENT OF THE PENDULUM
#   (2) STROBE START ANGLE SHOULD BE ADDED (NOW REGARDED AS 0)
#   (3) BASINS OF ATTRACTION GRAPHS
#   (4) FFT
##################################################################
# variables :
# w  : angular frequency
# th : angle
# ph : driven angle
# wd : driven angular frequency
# q  : damping factor
# g  : driven force amplitude
##################################################################
# ************** HOW TO USE **************
#
# VARIABLES TO USE IN FUNCTIONS : 
# - var = numpy.array([w,th,ph])
# - cnst = numpy.array([q,g,wd])
# - res is the output of solve functions
# - dt is time interval (optional)
# - plotstep defines step between plotting points (optional)
# - steady defines after how many periods to draw. (optional)
# - strobe is wd/ws ratio (optional)
#
# FUNCTIONS WITH DESCRIPTIONS :
# -solve by runge-kutta method :
# solve(var,cnst,'dt','steps','plotstep')
# -solve by euler method : 
# solve2(var,cnst,'dt','steps','plotstep')
# - plot phase space :
# phase_space(res,cnst,'plotstep','steady')
# - plot poincare section :
# poincare_sec(res,cnst,'strobe','plotstep','steady')
# - plot t versus * graphs :
# - calculates fourier transform :
# - plots trajectory :
# - plot basins of attraction :
##################################################################

from numpy import sin,cos,pi,array,zeros,rint
#import matplotlib 
#matplotlib.use("gtkcairo")
import pylab as plt

##################################################################
# dw/dt
def fw (t,var,cnst) :
  return -var[0]/cnst[0]-sin(var[1])+cnst[1]*cos(var[2])

# dth/dt
def fth (t,var,cnst) :
  return var[0]

# dph/dt  
def fph (t,var,cnst) :
  return cnst[2]

##################################################################
# just runge-kutta for this problem  
def solve (var,cnst,dt=0.01,steps=100000,plotstep=100) :
  f = (fw,fth,fph)
  t = 0
  tmp = var
  res = ([var[0]],[var[1]],[var[2]],[t])
  k = zeros((4,3))
  
  i = 0
  while (i < steps) :
    k[0] = dt * array([f[0](t,tmp,cnst),f[1](t,tmp,cnst),f[2](t,tmp,cnst)])
    t = t + 0.5*dt
    tmp1 = tmp+k[0]*0.5
    k[1] = dt * array([f[0](t,tmp1,cnst),f[1](t,tmp1,cnst),f[2](t,tmp1,cnst)])
    tmp1 = tmp+k[1]*0.5
    k[2] = dt * array([f[0](t,tmp1,cnst),f[1](t,tmp1,cnst),f[2](t,tmp1,cnst)])
    t = t + 0.5*dt
    tmp1 = tmp+k[2]
    k[3] = dt * array([f[0](t,tmp1,cnst),f[1](t,tmp1,cnst),f[2](t,tmp1,cnst)])
    tmp = array([res[0][i],res[1][i],res[2][i]]) + 1./6. * (k[0] + k[3] + 2*(k[1] + k[2]))
    res[0].append(tmp[0])
    res[1].append(tmp[1]) 
    res[2].append(tmp[2])
    res[3].append(t)
    
    i = i + 1  
  
  return res
  
##################################################################
# euler method
def solve2(var,cnst,dt=0.0001,steps=100000,plotstep=100) :
  f = (fw,fth,fph)
  t = 0
  tmp = var
  res = ([var[0]],[var[1]],[var[2]],[t])
  k = zeros(3)
  
  i = 0
  while (i < steps) :
    k = dt * array([f[0](t,tmp,cnst),f[1](t,tmp,cnst),f[2](t,tmp,cnst)])
    tmp = array([res[0][i],res[1][i],res[2][i]]) + k
    res[0].append(tmp[0])
    res[1].append(tmp[1]) 
    res[2].append(tmp[2])
    res[3].append(t+dt)
    t = t + dt
    i = i + 1
  
  return res
  
##################################################################  
# plots phase space
def phase_space(res,cnst,plotstep=100,steady=30) :
    
  arr = ([],[])
  time = 2*pi/cnst[2]*steady
  i = 0
  while (i < len(res[0])) :
    if (time <= res[3][i]) :
      arr[0].append(res[0][i])
      arr[1].append(res[1][i])
    i = i + 1
    
  i = 0
  while (i < len(arr[0])) :
    if ((arr[1][i] % (2*pi)) > (arr[1][i] % pi)) : 
      arr[1][i] = - (pi - (arr[1][i] % pi))
    else : 
      arr[1][i] = arr[1][i] % pi
    i = i + 1 
  

  plt.plot(arr[1][0::plotstep],arr[0][0::plotstep],'k.',markersize=1)
  plt.xlim(-pi,pi)
  plt.ylim(-pi,pi)
  plt.show()  

##################################################################
# plots poincare sections
def poincare_sec(res,cnst,strobe=1.,plotstep=100,steady=30) :
    
  sect_res = ([],[],[])
  time = 2*pi*strobe/cnst[2]
  
  i = 0
  j = 0
  while (i < len(res[0])) :
    if (time*j <= res[3][i]) :
      if (j >= steady) :
        sect_res[0].append(res[0][i])
        sect_res[1].append(res[1][i])
        sect_res[2].append(res[2][i])
      j = j + 1
    i = i + 1
  
  i = 0
  while (i < len(sect_res[0])) :
    if ((sect_res[1][i] % (2*pi)) > (sect_res[1][i] % pi)) : 
      sect_res[1][i] = - (pi - (sect_res[1][i] % pi))
    else :
      sect_res[1][i] = sect_res[1][i] % pi
    i = i + 1 
  
  plt.plot(sect_res[1][0::plotstep],sect_res[0][0::plotstep],'k*',markersize=5)
  plt.xlim(-pi,pi)
  plt.ylim(-pi,pi)
  plt.show()  
  
##################################################################
# plots t versus th or w
def plott(res,opt="w",timeint=0.1) :
  
  i = 0
  while (i < len(res[0])) :
    if ((res[1][i] % (2*pi)) > (res[1][i] % pi)) : 
      res[1][i] = - (pi - (res[1][i] % pi))
    else :
      res[1][i] = res[1][i] % pi
    i = i + 1 
  
  arr = ([],[],[],[])
  i = 0
  j = 0
  while (i < len(res[0])) :
    if (timeint*j <= res[3][i]) :
      arr[0].append(res[0][i])
      arr[1].append(res[1][i])
      arr[2].append(res[2][i])
      arr[3].append(res[3][i])
      j = j + 1
    i = i + 1
  
  if (opt == "w") : plt.plot(arr[3],arr[0])
  if (opt == "th") : plt.plot(arr[3],arr[1])
  if (opt == "ph") : plt.plot(arr[3],arr[2])
  plt.show()  
    
##################################################################      
# fourier transform
def fourier(res) :

  dt = res[3][1] - res[3][0]
  term1 = lambda t,w,dt : res[1][int(rint(t/dt))]*cos(w*t)*dt
  term2 = lambda t,w,dt : -res[1][int(rint(t/dt))]*sin(w*t)*dt
  
  power = ([],[])
  i = 0
  while (i < 100) :
    w = i/50.
    sum1 = 0
    sum2 = 0
    for t in res[3] :
      sum1 = sum1 + term1(t,w,dt)
      sum2 = sum2 + term1(t,w,dt)
    power[0].append((sum1 + sum2)**2)
    power[1].append(w)
    i = i + 1
  return power
  
##################################################################
#trajectory of the pendulum
def trajectory(res,cnst,steady=[30,34]) :

  place = ([],[])
  radius = 1.
  i = 1
  dir = "0"
  time = 2*pi/cnst[2]*array(steady)
  while(i < len(res[0])) :
    if (time[0] <= res[3][i]) :
      if (dir == "0") :
        place[0].append(radius * cos(res[1][i]))
        place[1].append(radius * sin(res[1][i])) 
        if (res[1][i+1] > res[1][i]) : 
          dir = "+"
          point = (res[1][i],radius)
        else : 
          dir = "-"
          point = (res[1][i],radius)
      elif (res[1][i] > res[1][i-1] and dir == "-") :
        dir = "+"
        radius = radius + 0.2 
        point = (res[1][i],radius)
        place[0].append(radius * cos(res[1][i-1]))
        place[1].append(radius * sin(res[1][i-1])) 
        place[0].append(radius * cos(res[1][i]))
        place[1].append(radius * sin(res[1][i])) 
      elif (res[1][i] < res[1][i-1] and dir == "+") : 
        dir = "-"
        radius = radius + 0.2
        point = (res[1][i],radius)
        place[0].append(radius * cos(res[1][i-1]))
        place[1].append(radius * sin(res[1][i-1])) 
        place[0].append(radius * cos(res[1][i]))
        place[1].append(radius * sin(res[1][i])) 
      else : 
        if ((radius == point[1]) and ((dir == '+' and point[0]+2*pi <= res[1][i+5]) or (dir == '-' and point[0]-2*pi >= res[1][i+5]))) :
          place[0].append(radius * cos(res[1][i]))
          place[1].append(radius * sin(res[1][i])) 
          radius = radius + 0.2
          point = (res[1][i],radius)
          place[0].append(radius * cos(res[1][i]))
          place[1].append(radius * sin(res[1][i])) 
        else :  
          place[0].append(radius * cos(res[1][i]))
          place[1].append(radius * sin(res[1][i]))
    if (time[1] <= res[3][i]) :
      break
    i = i + 1 
    
  #return place
  plt.plot([place[0][0],place[0][-1]],[place[1][0],place[1][-1]],'r*')
  #plt.plot(place[0][20::60],place[1][20::60],'b*')
  #plt.plot(place[0][40::60],place[1][40::60],'g*')
  #plt.plot(place[0][60::60],place[1][60::60],'y*')
  plt.plot(place[0],place[1],'k-',linewidth = 1.5)
  plt.show()  