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


def show_gridimage_frame(frame, recmodel, vmin, vmax, hyper_fraction, invert=False):
  gridframe = np.empty((recmodel.Ygridsize,recmodel.Xgridsize))
  for i in range(0,recmodel.Ygridsize):
    for j in range(0,recmodel.Xgridsize):
      indx = round(recmodel.condmap[i,j])
      fillout = vmax
      if invert: fillout = vmin
      if (indx == -1):
        gridframe[i,j] = fillout
      else:
        gridframe[i,j] = frame[indx]
  thisplot = plt
  thisplot.imshow(gridframe, cmap = cm.Blues, interpolation='bicubic', extent=(recmodel.x0, recmodel.xf, recmodel.y0, recmodel.yf), vmin=vmin, vmax=vmax)
  thisplot.colorbar()
  thisplot.title(("Hyperparameter = ", hyper_fraction ," * w0"))
  thisplot.plot(recmodel.boundary[:,0], recmodel.boundary[:,1])
  thisplot.scatter(recmodel.electrode_positions[:,0], recmodel.electrode_positions[:,1])
	
	
def show_image_frame(frame, nodes, triangles, boundary, electrode_positions, min = 1, max = -1): #Add palette options in params
  if (min > max):
    min = frame.min()
    max = frame.max()
  if (abs(max) > abs(min)):
    scalefactor = abs(max)
  else:
    scalefactor = abs(min)
  color=1/scalefactor*frame[:]             #color=delta_normalizado[:,delta]
  #color[0] = min/scalefactor
  #color[1] = max/scalefactor
  color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
  thisplot = plt
  thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color,cmap = cm.PuOr, vmin=min/scalefactor, vmax=max/scalefactor)
  thisplot.plot(boundary[:,0],boundary[:,1])
  thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])

def show_image_frame_onlypositive(frame, nodes, triangles, boundary, electrode_positions, min = 1, max = -1): #Add palette options in params
  if (min > max):
    min = frame.min()
    max = frame.max()
  if (abs(max) > abs(min)):
    scalefactor = abs(max)
  else:
    scalefactor = abs(min)
  color = []
  for k in range (0,frame.shape[0]):
    if (frame[k] > 0):
      color.append(frame[k]/scalefactor)
    else:
      color.append(0)
  color[0] = min/scalefactor
  color[1] = max/scalefactor
  color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
  thisplot = plt
  thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color,cmap = cm.PuOr)
  thisplot.plot(boundary[:,0],boundary[:,1])
  thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])
  
def show_image_frame_onlynegative(frame, nodes, triangles, boundary, electrode_positions, min = 1, max = -1): #Add palette options in params
  if (min > max):
    min = frame.min()
    max = frame.max()
  if (abs(max) > abs(min)):
    scalefactor = abs(max)
  else:
    scalefactor = abs(min)
  color = []
  for k in range (0,frame.shape[0]):
    if (frame[k] < 0):
      color.append(frame[k]/scalefactor)
    else:
      color.append(0)
  #color[0] = min/scalefactor
  #color[1] = max/scalefactor
  color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
  thisplot = plt
  thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color, cmap = cm.gray, vmin=min/scalefactor, vmax=max/scalefactor)  # cmap = cm.PuOr
  thisplot.plot(boundary[:,0],boundary[:,1])
  thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])

def show_image_frame_dif(frame1, frame2, nodes, triangles, boundary, electrode_positions):
  color=frame2[:]-frame1[:]             #color=delta_normalizado[:,delta]
  color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
  thisplot = plt
  thisplot.tripcolor(nodes[0], nodes[1], triangles, facecolors=color,cmap = cm.PuOr)
  thisplot.plot(boundary[:,0],boundary[:,1])
  thisplot.scatter(electrode_positions[:,0], electrode_positions[:,1])
  return
  
def get_wavefunction(raw_image):
  wf = []
  deriv = []
  dims = np.shape(raw_image.RawImage)
  for k in range(0,dims[1]):
    wf.append(sum(raw_image.RawImage[:,k]))
    if (k > 0):
      deriv.append(wf[k]-wf[k-1])
  return wf, deriv
  
def get_frame_pixel_from_coordinates(frame, x, y):
  pass
  
def get_frame_pixel_from_pseudopolar(frame, r, theta):
# r is a normalized radius (being 0 the origin and 1 the boundary for the given theta
# theta: angle in degrees respect to x axis
  pass

def get_sd_image(raw_image, initial_frnr, final_frnr, threshold):
  [sizeFrame, NrFrames]  = np.shape(raw_image)
  if (initial_frnr > final_frnr)|(final_frnr > NrFrames):
    print("*** ERROR in frame nr limits ***")
    return
  aux = []
  for k in range(0,sizeFrame):
    aux.append(np.std(raw_image[k,slice(initial_frnr, final_frnr+1)]))
  aux_array = np.array(aux)
  max_sd = aux_array.max()
  for k in range(0,sizeFrame):
    if (aux[k]/max_sd < threshold):
      aux_array[k] = 0
  return aux_array

def get_linreg_image(raw_image, initial_frnr, final_frnr, threshold):
  [sizeFrame, NrFrames]  = np.shape(raw_image)
  if (initial_frnr > final_frnr)|(final_frnr > NrFrames):
    print("*** ERROR in frame nr limits ***")
    return
  waveform = [np.mean(raw_image[:,k]) for k in range(initial_frnr, final_frnr+1)]
  aux = []
  for pix in range(0,np.shape(raw_image)[0]):
    pixwaveform = [raw_image[pix,k] for k in range(initial_frnr, final_frnr+1)]
    aux.append(np.polyfit(waveform, pixwaveform, 1)[0])
  aux_array = np.array(aux)
  max_sd = aux_array.max()
  for k in range(0,sizeFrame):
    if (aux[k]/max_sd < threshold):
      aux_array[k] = 0
  return aux_array
