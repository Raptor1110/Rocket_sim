import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as int1d
import numpy.linalg as la

class Rocket: #Objeto Foguete
    def __init__(self,dry_mass,prop_mass,diammeter,body_length,thrust_curve):
        self.dry_mass = dry_mass
        self.propellant_mass = prop_mass
        self.diammeter = diammeter
        self.body_length =body_length
        self.thrust_curve = thrust_curve

        # Motor

        Initial_Thrust = np.array([[0,0]])                                          # Inicio da queima com t = 0 e empuxo = 0 
        Burn = np.loadtxt(self.thrust_curve,skiprows=11)                            # Curva de empuxo do arquivo
        self.Thrust_curve = np.append(Initial_Thrust,Burn,axis=0)                   # Curva completa de empuxo

        shape = self.Thrust_curve.shape
        burn_time = self.Thrust_curve[shape[0]-1,0]
            
        self.Total_burn_time = burn_time                                            # Tempo total de queima do motor

    def Thrust(self,t):
        if t<= self.Total_burn_time:
            Engine_Thrust = int1d(self.Thrust_curve[:,0],self.Thrust_curve[:,1])    # Interpolação do empuxo
            return Engine_Thrust(t)
        else:
            return 0

    def mass(self,t):                                                               # Massa em tempo real
        if t < self.Total_burn_time:
            mass_data = np.array([ self.dry_mass + self.propellant_mass ,self.dry_mass])
            mass_time = np.array([0,self.Total_burn_time])
            Mass = int1d(mass_time,mass_data)
            return Mass(t)
        else:
            return self.dry_mass

    def get_thrust_curve(self):                 # Plota a curva de empuxo utilizada
        fig = plt.figure(2)
        plt.plot(self.Thrust_curve[:,0],self.Thrust_curve[:,1])
        plt.xlabel('Time')
        plt.ylabel('Thrust')
        plt.grid()
        plt.show()

def Wind(Wind_speed,Wind_direction):
    if Wind_speed >0:
        Cd = 1.1

        Wx = Wind_speed * math.cos(math.radians(Wind_direction))
        Wy = -Wind_speed * math.sin(math.radians(Wind_direction))
        Wz = 0  

        Wind_velocity = np.array([Wx,Wy,Wz])
    
        Wind_drift = 0.5 * 1.225 * 0.01 * np.square(Wind_velocity) * Cd * (Wind_velocity/Wind_speed)
    else:
        Wind_drift = 0
    return Wind_drift

def Rocket_attitude(v,t):

    v_norm = la.norm(v)
    Attitude_vector = ( v/v_norm ) 


    return Attitude_vector

def Drag(v,t):

    Cd = Vesper.Cd
    S = math.pi * ((Vesper.diammeter * 0.5)**2)
    rho = 1.225
    # Versor do arrasto
    Attitude = Rocket_attitude(v,t)
    
    # Aerodynamic Drag
    D = Attitude * -0.5 * rho * S * Cd * np.square(v)
    return D

def Thrust(v,t):
    
    Attitude = Rocket_attitude(v,t)
    Thrust = Attitude * Vesper.Thrust(t)
    
    return Thrust

def acceleration(t,v):
    Mass = Vesper.mass(t)
    if t>=0.05: # Foguete no solo 
        g = np.array([0,0,-9.81])
    else:       # Motor gera empuxo
        g = np.array([0,0,0])

    Wind_speed = 0
    Wind_direction = 0


    total_burn = Vesper.Total_burn_time
    Total_Drag= Drag(v,t)
    Wind_Force = Wind(Wind_speed,Wind_direction) # Wind Force 
    
    if t < total_burn:
        Engine_Thrust = Thrust(v,t) 
        a = (Engine_Thrust/Mass) + (Total_Drag/Mass) + g + (Wind_Force/Mass)
    
    #elif t>5 and t < 5.4:
        #Engine_Thrust = Thrust(v,t-5)
        #a = (Engine_Thrust/Mass) + (Total_Drag/Mass) + g + (Wind_Force/Mass)
    else:
        a = (Total_Drag/Mass) + g + (Wind_Force/Mass)
    return a

def RK4(f,v0,t0,tf,dt): # Metodo RungeKutta

    t = np.arange(t0,tf,dt) # vetor tempo de t0 até tf com intervalos de tempo dt
    

    nt = t.size # numero de tempos
    n_row,n_col = v0.shape # numero de variaveis em x
    v = np.zeros((n_row,nt))
    x = np.zeros((n_row,nt)) # Trajetoria 
    


    v[:,0] = v0[:,0]
    

    for k in range(0,nt-1):

        k1 = dt * f(t[k],v[:,k])
        k2 = dt * f(t[k] + dt/2 ,v[:,k] + k1/2)
        k3 = dt * f(t[k] + dt/2 ,v[:,k] + k2/2)
        k4 = dt * f(t[k] + dt ,v[:,k] + k3)


        dv = (k1 + (2*k2) + (2*k3) + (k4))/6 # diferencial de x
        

        v[:,k+1] = v[:,k] + dv # vetor de estados no momento desejado 

        dx = v[:,k]*dt

        x[:,k+1] = x[:,k] + dx # vetor de espaços
        
        if x[2,k+1]<0: #Chegou no chão
            break


    return x,v,t

##################### DATA #############################

Trhust_Curve = 'meteor-RASP_teste_1_.eng'    #  'Vesper.eng'# Nome do arquivo de curva de empuxo

Dry_Mass = 2.25                                                             # Massa sem propelente 
Prop_mass = 0.25                                                            # Massa do propelente 
diammeter = 0.085                                                           # Diametro do bodytube
body_length = 1074                                                          # Comprimento do bodytube
Vesper = Rocket(Dry_Mass ,Prop_mass, diammeter, body_length, Trhust_Curve)
Vesper.Cd = 0.41

#################### SIM CONFIG #########################

Launch_angle = 0
Launch_Heading = 0
Launch_speed = 0.01

save_data =  False

###################### DATA ARAYS #########################

vx = Launch_speed * math.sin(math.radians(Launch_angle)) * math.cos(math.radians(Launch_Heading))
vy = Launch_speed * math.sin(math.radians(Launch_angle)) * -math.sin(math.radians(Launch_Heading))
vz = Launch_speed * math.cos(math.radians(Launch_angle))

v0 = np.array([[vx],
               [vy],
               [vz]])

x,v,t = RK4(acceleration,v0,0,22,0.001)

if save_data:

    np.savetxt('Vertical_speed.txt',v[2,:])
    np.savetxt('time.txt',t)
    np.savetxt('Altitude.txt',x[2,:])

########################## PLOTS ############################

plt.subplot(2,2,1)
plt.plot(t,x[2,:])
plt.xlabel('Time')
plt.ylabel('Altitude')
plt.grid()

plt.subplot(2,2,2)
plt.plot(t,v[2,:])
plt.xlabel('Time')
plt.ylabel('Vertical Speed')

ax = plt.subplot(2,2,3,projection='3d')
ax.plot(x[0,:],x[1,:],x[2,:])
ax.set_xlabel('Downrange X')
ax.set_ylabel('Downrange Y')
ax.set_zlabel('Altitude Z')
plt.show()

Vesper.get_thrust_curve()

