import numpy as np
from math import erf


bounds = [1000,1000]


# Parameters for user distribution, For the paper simulation they used a truncated guassian distribution, however any arbitrary distribution can be used
mean = [250,300]
std_dev = [1000,1000]
scale_factor = 2*np.pi*std_dev[0]*std_dev[1]*erf((bounds[0]-mean[0])/(np.sqrt(2)*std_dev[0]))*erf((bounds[1]-mean[1])/(np.sqrt(2)*std_dev[1]))
def f(x,y):
    return 1/scale_factor*np.exp(-((x-mean[0])**2/(2*std_dev[0]**2) + (y-mean[1])**2/(2*std_dev[1]**2)))
N = 300 #number of users
u = 10e6 #load per user, 10Mbps

#UAV parameters
num_UAVs = 5
x_locations = np.random.uniform(0,bounds[0],num_UAVs)
y_locations = np.random.uniform(0,bounds[1],num_UAVs)
h_UAV = 1000*np.ones(num_UAVs) #meters, altitude of UAV
hover_time = 30*60*np.ones(num_UAVs) #30 minutes, hover time of UAV

#communication parameters
f_c = 2e9 #2GHz, carrier frequency
P = .5*np.ones(num_UAVs) #W, transmit power
B = 1 #1MHz, bandwidth
N_O = -170 #dBm/Hz noise power
mu_los = 10**(3/10) #dB, Additional path loss for line of sight
mu_nlos = 10**(23/10) #dB, Additional path loss for non-line of sight
alpha = np.ones(num_UAVs)*0.01 #control time factor
b1 = 0.36
b2 = 0.21
beta_interference_factor = 1



#parameters for numerical integration
num_x_points = 40
num_y_points = 40

x_int = np.linspace(0,bounds[0],num_x_points)
y_int = np.linspace(0,bounds[1],num_y_points)
dx = x_int[1]-x_int[0]
dy = y_int[1]-y_int[0]






