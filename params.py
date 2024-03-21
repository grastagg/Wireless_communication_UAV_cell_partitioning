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
h_UAV = 200*np.ones(num_UAVs) #meters, altitude of UAV

#communication parameters
f_x = 2e9 #2GHz, carrier frequency
P = .5*np.ones(num_UAVs) #W, transmit power
B = 1e6 #1MHz, bandwidth
N_O = -170 #dBm/Hz noise power
mu_los = 3 #dB, Additional path loss for line of sight
mu_nlos = 23 #dB, Additional path loss for non-line of sight
alpha = np.ones(num_UAVs)*0.01 #control time factor
b1 = 0.36
b2 = 0.21




