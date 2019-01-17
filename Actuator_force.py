import numpy as np

#stall torque from motor (116 rpm servocity premium gearmotor)
tmax =1.13 #N-m

#efficiency lead screw
zeta_s = 0.15

#efficiency of pulley system
zeta_p = 0.97

#lead screw lead
lead = 1/1000.0 #m/rotation

# Max Actuator Force
F_act = tmax*zeta_p*zeta_s*2*np.pi/(2*lead)

#-----------------------------------------------------

#applied Loading
F_applied = 22.24*3.0/5.0 #newtons #3 lbs

#geometry of Scissor Linkage
l_bent = 4.0
d = 8.0
l_r =42
l_total = 5.0*l_r

#-----------------------------------------------------

#Critical Parameters
w_min = 24.0 #mm
w_max = 82.47 #mm


def L_theta(theta,l_bent, l_total):
    return l_total*np.sin(theta) + l_bent*np.cos(theta)

def w_theta(theta, l_r, l_total):
    return 2*l_r*np.cos(theta) - l_bent*np.sin(theta)

def theta_min(l_r):
    return np.arcsin(d/(2*l_r))

theta_max = 68/180*np.pi

def E_ratio(l_b,l_o, l_total):
    return L_theta(theta_max, l_b, l_total)/(L_theta(theta_min(l_o), l_b, l_total) + l_b)

def F_nut(theta, F_applied,l_bent, l_r, l_total ):
    return F_applied/2.0*( l_total*np.cos(theta) - l_bent*np.sin(theta) ) 			\
					/(l_r*np.sin(theta) + l_bent*np.cos(theta) )

# F_req = lambda theta:  F_applied/2.0*( l_total*np.cos(theta) - l_bent*np.sin(theta) ) 			\
# 					/(l_r*np.sin(theta) + l_bent*np.cos(theta) )

#-----------------------------------------------------


L_max = L_theta(theta_max, l_bent, l_total)
L_min = L_theta(theta_min(l_r), l_bent, l_total)

#inertia
I_com = 0.85 #kg/m^2
I_stage1 = I_com + 339.5/1000.0*(L_max/1000.0)**2
I_stage2 = I_com + 339.5/1000.0*(2*L_max/1000.0)**2
I_base = I_com + I_stage1 + I_stage2

omega = (np.pi/2.0)/5.0
alpha = omega/2.0
g = 9.8
m_stage = 339.5/1000.0


F_app = 1/(168.24/1000.0)*(I_base*alpha + m_stage*g*L_max/1000.0 + 2*m_stage*g*L_max/1000.0)

def w_l(L, l_bent, l_o, l_e):
    a1 = 2.0/(l_bent**2+ l_e**2)
    a2 = (l_bent**2 + l_o*l_e)*np.sqrt(l_bent**2+ l_e**2 - L**2)
    a3 = (L*l_bent)*(l_e - l_o)
    return a1*(a2-a3)


def print_stats():
    print("F_nut: " + str( F_nut(theta_min(l_r), F_app, l_bent, l_r, l_total) ) )
    print("F_act: " + str(F_act) )
    print("w_max: " + str(w_theta(theta_min(l_r), l_r, l_total)) )
    print("theta_min: " + str( theta_min(l_r)*180.0/(np.pi) ) )
    print("E_ratio: " + str( E_ratio(l_bent, l_r, l_total) ) )
    print("Min Length: " + str(L_min) )
    print("Max Length: "+ str(L_max) )

    print("I_stage1: " + str(I_stage1) )
    print("I_stage2: " + str(I_stage2) )
    print("I_base: " + str(I_base) )
    
    print("alpha: " + str(alpha) )
    print("F_app Grav Approximation: " + str(F_app) )