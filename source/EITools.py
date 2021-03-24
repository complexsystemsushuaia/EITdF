#############################################################################################
#
# Tools module 
# (part of EITdF)
#
# Daniel Badagnani
# Instituto de Ciencias Polares, Ambiente y Recursos Naturales
# Universidad Nacional de Tierra del Fuego - Ushuaia - Argentina
#
# 2020
#
#############################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import io

"""
# Esta funcionalidad se implementó en la clase EITdata
def add_noise(signal,SNR,noise_covariance):
	# Adds autocorrelated noise to the signal
	# signal is a Nmeas x Nframes array.
	# SNR: Signal to Noise Ratio
	# noise_covariance is a Nmeas x Nmeas array.
	#
	# Get eigenValues and eigenVectors of the covariance matrix
	eVa, eVe = np.linalg.eig(noise_covariance)
	# Get the rotation from uncorrelated to correlated noise
	T = eVe.dot(np.diag(np.sqrt(eVa)))
	#
	[Nmeas, Nframes] = np.shape(signal)
	# Generate frames of noise
	uncorrelated_noise = [np.normal(0,1,Nframes) for k in range(Nmeas)]
	uncorrelated_noise_frames = np.vstack(uncorrelated_noise)
	noise_frames = T.dot(uncorrelated_noise_frames)
	# Get renormalization factor for noise to the specified SNR
	noise_level = np.linalg.norm(noise_frames)/(SNR*np.linalg.norm(signal))
	# Add noise frames to data frames as returned data
	return signal + noise_frames * noise_level
"""
	
def show_image_frame(frame, nodes, triangles, boundary, electrode_positions): #Add palette options in params
	color=frame[:]             #color=delta_normalizado[:,delta]
	color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
	thisplot = plt
	thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color,cmap = cm.PuOr)
	thisplot.plot(boundary[:,0],boundary[:,1])
	thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])
	return

def show_image_frame_dif(frame1, frame2, nodes, triangles, boundary, electrode_positions):
	color=frame2[:]-frame1[:]             #color=delta_normalizado[:,delta]
	color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
	thisplot = plt
	thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color,cmap = cm.PuOr)
	thisplot.plot(boundary[:,0],boundary[:,1])
	thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])
	return
	
def get_frame_pixel_from_coordinates(frame, x, y):
	pass
	
def get_frame_pixel_from_pseudopolar(frame, r, theta):
	# r is a normalized radius (being 0 the origin and 1 the boundary for the given theta
	# theta: angle in degrees respect to x axis
	pass
	
 
