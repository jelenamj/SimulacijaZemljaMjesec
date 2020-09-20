import numpy as np
import matplotlib.pyplot as plt 
from  decimal import Decimal as D
import decimal 

decimal.getcontext().prec = 30

#Runge-Kutta method for Earth
def RK4E(t,Z,V):
    Zn=0.0
    if  Z==Re[0]:
        Zn=Rm[0]        
    elif Z==Re[1]:
         Zn=Rm[1]        
    else:
        print('error')
        
    k1=V
    l1=F2(Z,Zn)
    
    k2=V+dt*l1/2
    l2=F2(Z+dt*k1/2,Zn)
    
    k3=V+dt*l2/2
    l3=F2(Z+dt*k2/2,Zn)
    
    k4=V+dt*l3
    l4=F2(Z+dt*k3,Zn)
    
    return [dt*(l1 + D(2)*l2 + D(2)*l3+l4)/D(6),dt*(k1+ D(2)*k2+ D(2)*k3+k4)/D(6)]

def F2(Z,Zn): return -G*MassSun*Z/Res**3-G*MassMoon*(Z-Zn)/Rme**3

#Runge-Kutta method for Moon
def RK4M(t,Z,V):
    
    Zn=0.0
    if  Z==Rm[0]:
        Zn=Re[0]        
    elif Z==Rm[1]:
         Zn=Re[1]        
    else:
        print('error')
    
    k1=V
    l1=G2(Z,Zn)
    
    k2=V+dt*l1/2
    l2=G2(Z+dt*k1/2,Zn)
    
    k3=V+dt*l2/2
    l3=G2(Z+dt*k2/2,Zn)
    
    k4=V+dt*l3
    l4=G2(Z+dt*k3,Zn)
    
    return [dt*(l1 + D(2)*l2 + D(2)*l3+l4)/D(6),dt*(k1+ D(2)*k2+ D(2)*k3+k4)/D(6)]

def G2(Z,Zn): return -G*MassSun*Z/Rms**3-G*MassEarth*(Z-Zn)/Rme**3


#Masses of Sun, Earth and Moon
MassSun=D(1.989e30)
MassEarth=D(5.9722e24)
MassMoon=D(7.34767e22)

#Universal gravitational constant
G=D(6.67e-11)


#Distances Earth-Sun (Perihelion), Earth-Moon , Sun-Moon
Res=D(147099760000.0)
Rme=D(405696000.0)
Rms=Rme+Res


RmsPrevioust=Rms
ResPrevioust=Res

#Position and velocity of Earth and Moon
Re=np.array([Res,D(0.0)],np.dtype(decimal.Decimal))
Rm=np.array([Rms,D(0.0)],np.dtype(decimal.Decimal))
Vz=np.array([D(0.0),D(30280.0)],np.dtype(decimal.Decimal))
Vm=np.array([D(0.0),D(30280.0-1024.0)],np.dtype(decimal.Decimal))

RKEarth=[D(0.0),D(0.0)]
RKMoon=[D(0.0),D(0.0)]
Rz1=[D(0.0),D(0.0)]
Rm1=[D(0.0),D(0.0)]
Vz1=[D(0.0),D(0.0)]
Vm1=[D(0.0),D(0.0)]

#All X and Y components 
Xe=[]
Ye=[]
Xm=[]
Ym=[]

MonthSplit=[]

dt=D(60.0)
time=D(0.0)

while time<D(60*60*24*366) :
    Rms=np.linalg.norm([Rm[0],Rm[1]])
    Res=np.linalg.norm([Re[0],Re[1]])
    Rme=np.linalg.norm([(Rm[0]-Re[0]),(Rm[1]-Re[1])])
      
    #Calculation of position and velocity for Earth in t+dt
    #X component
    RKEarth=RK4E(time, Re[0], Vz[0])
    Rz1[0]=Re[0]+RKEarth[1]
    Vz1[0]=Vz[0]+RKEarth[0]
    #Y component
    RKEarth=RK4E(time, Re[1], Vz[1])
    Rz1[1]=Re[1]+RKEarth[1]
    Vz1[1]=Vz[1]+RKEarth[0]
    
    #Calculation of position and velocity for Moon in t+dt
    #X component
    RKMoon=RK4M(time,Rm[0],Vm[0])
    Rm1[0]=Rm[0]+RKMoon[1]
    Vm1[0]=Vm[0]+RKMoon[0]
    #Y component
    RKMoon=RK4M(time,Rm[1],Vm[1])
    Rm1[1]=Rm[1]+RKMoon[1]
    Vm1[1]=Vm[1]+RKMoon[0]
    
    #All X and Y components 
    #Earth
    Xe.append(Rz1[0])
    Ye.append(Rz1[1])
    #Moon
    Xm.append(Rm1[0])
    Ym.append(Rm1[1])
    
    Re=Rz1
    Rm=Rm1
    Vm=Vm1
    Vz=Vz1
    time=time+dt
    
    Day360=60*60*24*360
    AproximChangeDistance=2000000
    
    
    if (time>D(Day360) and (Ye[-1])>D(0.0) and (Ye[-2])<D(0.0)):
        print('time before new year (s) ',time-dt,',  in days',(time-dt)/D(60*60*24))
        print(' Coordinates of Earth before new year     X= ',Xe[-2],'     Y=  ',Ye[-2])
        
        print('time after new year (s) ',time,',  in days',time/D(60*60*24))
        print(' Coordinates of Earth after new year     X= ',Xe[-1],'     Y=  ',Ye[-1])
        
    if (RmsPrevioust-ResPrevioust)>0 and (Rms-Res)<0:
        MonthSplit.append(time)
        
    RmsPrevioust=Rms
    ResPrevioust=Res


LenghtMonths=[]
for i in range(len(MonthSplit)):
    if i>0:
        LenghtMonths.append(MonthSplit[i]-MonthSplit[i-1])
        
Month=sum(LenghtMonths)/len(LenghtMonths)
print('Average lenght of synodic month is (s):   ',Month)
print('Average lenght of synodic month is (days):   ',Month/D(60*60*24))
plt.plot(Xm, Ym,'b')
plt.plot(Xe, Ye,'y')
plt.title('h=60s')
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
plt.savefig('h60_oneyear.png')
plt.show()
