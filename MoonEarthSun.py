import numpy as np
import matplotlib.pyplot as plt 
from  decimal import Decimal as D
import decimal 

decimal.getcontext().prec = 30

class Body:
    def __init__(self, mass, position, velocity):
        self.Mass=mass
        self.Velocity=velocity
        self.Position=position
        self.PositionNextdt=[D(0.0),D(0.0)]
        self.VelocityNextdt=[D(0.0),D(0.0)]
        self.X=[]
        self.Y=[]
    
#Rs-distance object from Sun
    @staticmethod
    def F2(Z,Zn,Rs,Mass): return -G*Sun.Mass*Z/Rs**3-G*Mass*(Z-Zn)/Rme**3
    
    @staticmethod
    def RungeKutta4(t,Z,V,Rs,Mass):
            Zn=0.0
            if  Z==Earth.Position[0]:
                Zn=Moon.Position[0]        
            elif Z==Earth.Position[1]:
                 Zn=Moon.Position[1] 
            elif  Z==Moon.Position[0]:
                Zn=Earth.Position[0]        
            elif Z==Moon.Position[1]:
                 Zn=Earth.Position[1]     
                
            k1=V
            l1=Body.F2(Z,Zn,Rs,Mass)
            k2=V+dt*l1/2
            l2=Body.F2(Z+dt*k1/2,Zn,Rs,Mass)
            k3=V+dt*l2/2
            l3=Body.F2(Z+dt*k2/2,Zn,Rs,Mass)
            k4=V+dt*l3
            l4=Body.F2(Z+dt*k3,Zn,Rs,Mass)
            
            return [dt*(l1 + D(2)*l2 + D(2)*l3+l4)/D(6),dt*(k1+ D(2)*k2+ D(2)*k3+k4)/D(6)]

#Universal gravitational constant
G=D(6.67e-11)

#Distances (Perihelion) Earth-Sun, Earth-Moon , Sun-Moon
Res=D(147099760000.0)
Rme=D(405696000.0)
Rms=Rme+Res
#Position and velocity of Earth and Moon
Re=np.array([D(Res),D(0.0)],np.dtype(decimal.Decimal))
Rm=np.array([Rms,D(0.0)],np.dtype(decimal.Decimal))
Ve=np.array([D(0.0),D(30280.0)],np.dtype(decimal.Decimal))
Vm=np.array([D(0.0),D(30280.0-1024.0)],np.dtype(decimal.Decimal))

#Instances
Moon=Body(D(7.34767e22),Rm,Vm)
Earth=Body(D(5.9722e24),Re,Ve)
Sun=Body(D(1.989e30),[D(0.0),D(0.0)],[D(0.0),D(0.0)])

RmsPrevioust=Rms
ResPrevioust=Res
MonthSplit=[]

dt=D(60.0)
time=D(0.0)
Day1=D(60*60*24)

while time<D(60*60*24*366) :
    Rms=np.linalg.norm([Moon.Position[0],Moon.Position[1]])
    Res=np.linalg.norm([Earth.Position[0],Earth.Position[1]])
    Rme=np.linalg.norm([(Moon.Position[0]-Earth.Position[0]),(Moon.Position[1]-Earth.Position[1])])
      
    #Calculation of position and velocity for Earth in t+dt
    #X component
    RKEarth=Body.RungeKutta4(time, Earth.Position[0], Earth.Velocity[0],Res,Moon.Mass)
    Earth.PositionNextdt[0]=Earth.Position[0]+RKEarth[1]
    Earth.VelocityNextdt[0]=Earth.Velocity[0]+RKEarth[0]
    #Y component
    RKEarth=Body.RungeKutta4(time, Earth.Position[1], Earth.Velocity[1],Res,Moon.Mass)
    Earth.PositionNextdt[1]=Earth.Position[1]+RKEarth[1]
    Earth.VelocityNextdt[1]=Earth.Velocity[1]+RKEarth[0]
    
    #Calculation of position and velocity for Moon in t+dt
    #X component
    RKMoon=Body.RungeKutta4(time,Moon.Position[0],Moon.Velocity[0],Rms,Earth.Mass)
    Moon.PositionNextdt[0]=Moon.Position[0]+RKMoon[1]
    Moon.VelocityNextdt[0]=Moon.Velocity[0]+RKMoon[0]
    #Y component
    RKMoon=Body.RungeKutta4(time,Moon.Position[1],Moon.Velocity[1],Rms,Earth.Mass)
    Moon.PositionNextdt[1]=Moon.Position[1]+RKMoon[1]
    Moon.VelocityNextdt[1]=Moon.Velocity[1]+RKMoon[0]
    
    #All X and Y components 
    #Earth
    Earth.X.append(Earth.PositionNextdt[0])
    Earth.Y.append(Earth.PositionNextdt[1])
    #Moon
    Moon.X.append(Moon.PositionNextdt[0])
    Moon.Y.append(Moon.PositionNextdt[1])
    
    Earth.Position=Earth.PositionNextdt
    Moon.Position=Moon.PositionNextdt
    Moon.Velocity=Moon.VelocityNextdt
    Earth.Velocity=Earth.VelocityNextdt
    time=time+dt
    
    Day363=60*60*24*363
    
    if (time>D(Day363) and (Earth.Y[-1])>D(0.0) and (Earth.Y[-2])<D(0.0)):
        print('time before new year (s) ',time-dt,',  in days',(time-dt)/Day1)
        print(' Coordinates of Earth before new year     X= ',Earth.X[-2],'     Y=  ',Earth.Y[-2])
        print('time after new year (s) ',time,',  in days',time/Day1)
        print(' Coordinates of Earth after new year     X= ',Earth.X[-1],'     Y=  ',Earth.Y[-1])
        
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
print('Average lenght of synodic month is (days):   ',Month/Day1)
plt.plot(Moon.X, Moon.Y,'b')
plt.plot(Earth.X, Earth.Y,'y')
plt.title('h=60s')
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
plt.savefig('h60_oneyear.png')
plt.show()
