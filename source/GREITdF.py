#############################################################################################
#
# GREITdF module: routines for EIT data retrieving and reconstruction
# (part of EITdF)
#
# Daniel Badagnani
# Instituto de Ciencias Polares, Ambiente y Recursos Naturales
# Universidad Nacional de Tierra del Fuego - Ushuaia - Argentina
#
# 2020
#
#############################################################################################
#
# Four classes are defined: EITdata, EnvironOpt, ReconstrictionModel and RawImage
#
############### EITdata ################# 
# Imports raw data and produces raw dif images.
# If measure_active_electrodes is 1 (default), a data frame is all 256 dif voltages.
# Else, a dataframe is 16*(16-3)=208 dif voltages only between passive electrodes.
#
# FIELDS:
#
# EITdata.inputfile
# EITdata.measure_active_electrodes
# EITdata.environment
# EITdata.eit_data_meas
# EITdata.eit_data_meas_passive_electrodes

# EITdata.eit_data_meas_active_electrodes
# EITdata.EIT_normalized_dif
# EITdata.ref_dataframe
#
# Methods:
#
# Constructor(environ,inputfile,measure_active_electrodes=1)
# separate_dataframe_data()
# get_normalized_dataframes(ref_dataframe)
# estimate_contact_resistances()
#
############### EnvironOpt #################
# creates an object with the names of the relevant folders for GREITdF.
# The constructor takes as argument the absolute path to GREITdF.
#
############### ReconstructionModel #################
# self.inputpath  (string)
# self.outputpath (string)
# self.binpath    (string)
# self.y          (numpy array)
# self.d          (numpy sparse)
# self.rm         (numpy array)
# self.hyperparam (double)
#
# self.set_hyperparameter(hyperparam_fraction)
# self.get_raw_image(meas_filename,image_filename,meas_current=1)
# self.show_raw_image()
# self.show_raw_frame(frame_nr)  (show electrode resistances if available)
# self.filter_meas()     (returns no_current_meas and electrode resistance estimates)
#
############### RawImage #################
# Produces and stores the raw image.
#
# FIELDS:
#
# RawImage.rec_model (a ReconstructionModel object)
# RawImage.eitdata   (an EITdata object)
#
# METHODS:
#
# RawImage.raw_image(rec_model,eitdata)
# RawImage.calc_raw_image(reffframe) #refframe is the number of frame taken as reference. It is called by the constructor.
# RawImage.get_ref_frame()           #finds the minimum conductivity frame to act as a reference. Called by constructor.
# RawImage.ShowFrame(frame_number) #plots frame
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import io

class EnvironOpt(object):
  #
  # The object EnvironObj is the one containing system's environmental variables
  # for making them non-global and guarrantee modularity.
  #
  def __init__(self,basedir):
    self.basedir   = basedir
    self.modelsdir = basedir+"/models"
    self.inputdir  = basedir+"/data"
    self.bindir    = basedir+"/source/python"
    
  def change_dir(self,basedir):
    self.basedir   = basedir
    self.modelsdir = basedir+"/models"
    self.inputdir  = basedir+"/data"
    self.bindir    = basedir+"/source/python"


class EITdata:
  def __init__(self,environ,inputfile):
    self.inputfile                 = inputfile
    self.environment               = environ
    # Load input
    self.eit_data_meas = io.loadmat(self.environment.inputdir+"/"+self.inputfile)['FRAMES']
    self.set_data()


  def set_data(self):
    [nr_meas, nr_frames] = np.shape(self.eit_data_meas)
    self.nr_frames = nr_frames
    self.nr_meas   = nr_meas
    if (nr_meas == 256):
      self.measure_active_electrodes = 1
    else:
      self.measure_active_electrodes = 0
    self.eit_data_meas_passive_electrodes,self.eit_data_meas_active_electrodes = self.separate_dataframe_data()
    #
    # Order ciclically active electrodes: 15,0,1; 0,1,2; ... ;13,14,15.
    # Dim is 3*16, each triplet is an injection frame.
    #
    if (self.measure_active_electrodes == 1):
      for k in range(nr_frames):
        aux = self.eit_data_meas_active_electrodes[3*k*16]
        self.eit_data_meas_active_electrodes[3*k*16] = self.eit_data_meas_active_electrodes[3*k*16+2]
        self.eit_data_meas_active_electrodes[3*k*16+2] = aux
        aux = self.eit_data_meas_active_electrodes[3*k*16+45]
        self.eit_data_meas_active_electrodes[3*k*16+45] = self.eit_data_meas_active_electrodes[3*k*16+47]
        self.eit_data_meas_active_electrodes[3*k*16+47] = aux
		    
    
  def separate_dataframe_data(self):
    if self.measure_active_electrodes == 0:
      return self.eit_data_meas, []
    active        = []
    passive_meas  = []
    electrode_pairs = [14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1]
    for frame in range(self.nr_frames):
      passive_frame = []
      for meas_in_frame in range(self.nr_meas):
        inj      = meas_in_frame//16      # Injection frame
        meas_in_inj = meas_in_frame-16*inj    # Measurement within injection frame between 0 and 15
        skipped  = [electrode_pairs[inj+m+1] for m in range(3)]
        if meas_in_inj in skipped:
          active.append(self.eit_data_meas[meas_in_frame,frame])
        else:
          passive_frame.append(self.eit_data_meas[meas_in_frame,frame])
      if (len(passive_frame) > 0):
        passive_meas.append(passive_frame)
    return np.array(passive_meas).transpose(), active    
  
  def get_normalized_dataframes(self,ref_dataframe):
    rows, cols = np.shape(self.eit_data_meas_passive_electrodes)
    if (ref_dataframe+1  > cols) or (ref_dataframe < 0):
      print("WARNING: reference frame out of range. Using 0 instead")
      ref_dataframe = 0
    self.ref_frame = ref_dataframe
    self.EIT_normalized_dif = (self.eit_data_meas_passive_electrodes.transpose() / self.eit_data_meas_passive_electrodes.transpose()[ref_dataframe] - 1).transpose()
    
  
  def estimate_contact_resistances(self,injected_currents):
    if (self.measure_active_electrodes != 1):
      print("ERROR: Contact resistances cannot be estimated. No data on active electrodes.")
      return
    contact_resistances = []
    for k in range(16*self.nr_frames):
      v = abs(self.eit_data_meas_active_electrodes[3*k])+abs(self.eit_data_meas_active_electrodes[3*k-1])/2
      contact_resistances.append(v/injected_currents[k])
    self.contact_resistances = contact_resistances
    
        
  def add_noise(self,SNR,noise_covariance):
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
    [Nmeas, Nframes] = np.shape(self.eit_data_meas)
    # Generate frames of noise
    uncorrelated_noise = [np.random.normal(0,1,Nframes) for k in range(Nmeas)]
    uncorrelated_noise_frames = np.vstack(uncorrelated_noise)
    noise_frames = T.dot(uncorrelated_noise_frames)
    # Get renormalization factor for noise to the specified SNR
    noise_level = np.linalg.norm(self.eit_data_meas)/(SNR * np.linalg.norm(noise_frames))
    # Add noise frames to data frames as returned data
    self.eit_data_meas = self.eit_data_meas + noise_frames * noise_level
    self.set_data()
    
    
class ReconstructionModel(object):
  def __init__(self,model_name,hyperparameter_fraction,environ_options):
    # environ_options: an EnvironOpt object
    # 
    # 
    self.modelname      = model_name
    self.hyperparameter_fraction = hyperparameter_fraction
    self.y_matrix       = io.loadmat(environ_options.modelsdir+"/"+model_name+"/y.mat")['Y']
    self.d_matrix       = io.loadmat(environ_options.modelsdir+"/"+model_name+"/d.mat")['D']
    try:
      self.Sn = io.loadmat(environ_options.modelsdir+"/"+model_name+"/Sn.mat")['Sn']
    except OSError:
      print("Noise covariances not found. Using identity.")
      self.Sn = np.identity(208)
    self.IsTijonov = 1;
    try:
      self.J  = io.loadmat(environ_options.modelsdir+"/"+model_name+"/J.mat")['J']
    except OSError:
      self.IsTijonov = 0;
    self.recaspectratio = 1
    try:
      self.J_ar = io.loadmat(environ_options.modelsdir+"/"+model_name+"/J_aspectratio.mat")['J_ar'][0]
      self.ar_eigenval = io.loadmat(environ_options.modelsdir+"/"+model_name+"/aspectratio.eigenval.mat")['ar_eigenval']
    except OSError:
      self.recaspectratio = 0
      print("OOPS")
    filename = environ_options.modelsdir+"/"+model_name+"/noiselev.mat"
    with open(filename) as f:
      lines = f.readlines()
    k = 0
    isfloat = 0
    while (isfloat == 0):
      try:
        isfloat = 1
        self.hyperparameter = float(lines[k])
      except ValueError:
        k += 1
        isfloat = 0
    filename = environ_options.modelsdir+"/"+model_name+"/recogrid"
    with open(filename) as f:
      lines = f.readlines()
    self.Xgridsize = int(lines[0])
    self.Ygridsize = int(lines[1])
    # Compute inverse matrix
    self.set_hyperparameter(hyperparameter_fraction)
    # Triangulation (tri and nodes)
    self.tri                  = io.loadmat(environ_options.modelsdir+"/"+model_name+"/elems.mat")['elems']-1
    self.nodes                = (io.loadmat(environ_options.modelsdir+"/"+model_name+"/nodes.mat")['nodes']).transpose()
    self.boundary             = io.loadmat(environ_options.modelsdir+"/"+model_name+"/boundary.mat")["boundary_coords"]
    self.electrode_positions  = io.loadmat(environ_options.modelsdir+"/"+model_name+"/electrode_positions.mat")["electrode_positions"]
    # Triangulation to Grid
    self.x0     = self.nodes[0,0]
    self.xf     = self.nodes[0,self.Xgridsize]
    grdszx = (self.xf - self.x0)/self.Xgridsize
    self.y0     = self.nodes[1,0]
    self.yf     = self.nodes[1,(self.Xgridsize+1)*(self.Ygridsize+1)-1]
    grdszy = (self.yf - self.y0)/self.Ygridsize
    # self_condmap is the map from the grid to a row_image frame
    self.condmap = np.empty((self.Xgridsize, self.Ygridsize))
    for i in range(0,self.Ygridsize):
      for j in range(0,self.Xgridsize):
        self.condmap[i,j] = -1
    dims = np.shape(self.rm)
    for k in range(0,dims[0]):
      ta = round(self.tri[2*k,0])
      i  = round((self.nodes[0,ta]-self.x0)/grdszx)
      j  = round((self.nodes[1,ta]-self.y0)/grdszy)
      self.condmap[31-j,i] = k
      
      
  def set_hyperparameter(self, hyperparameter_fraction, movement_hyperparam = 1):
    self.hyperparameter_fraction = hyperparameter_fraction
    # Compute inverse matrix
    w = self.hyperparameter * self.hyperparameter_fraction
    if (self.IsTijonov):
      m_matrix            = self.J.dot(self.J.transpose()) + self.Sn*(w**2)
      self.rm             = self.J.transpose().dot(np.linalg.inv(m_matrix))
      print("TIJONOV")
    else:
      m_matrix            = self.y_matrix.dot(self.y_matrix.transpose()) + self.Sn*(w**2)
      m_inv               = np.linalg.inv(m_matrix)
      pjt_matrix          = self.d_matrix.dot(self.y_matrix.transpose())
      self.rm             = pjt_matrix.dot(m_inv)
      if (self.recaspectratio):
        MJ_ar             = self.J_ar.reshape(len(self.J_ar),1)
        self.proj_mov     = (1/self.ar_eigenval)*MJ_ar.dot(MJ_ar.transpose())
        self.proj_cond    = np.identity(np.shape(self.Sn)[0]) - self.proj_mov
        inv_proj          = self.proj_cond.dot(m_inv)
        self.rm           = pjt_matrix.dot(inv_proj)
        self.aspectratio_rm = self.J_ar.dot(self.proj_mov)
        

class RawImage:
  def __init__(self,rec_model,eitdata):
    self.rec_model = rec_model
    self.eitdata   = eitdata
    self.calc_raw_image(0)
      
  def calc_raw_image(self,refframe):
    self.eitdata.get_normalized_dataframes(refframe)
    self.raw_image=np.matmul(self.rec_model.rm, self.eitdata.EIT_normalized_dif)
    if (self.rec_model.recaspectratio):
      self.aspectratio = np.matmul(self.rec_model.aspectratio_rm, self.eitdata.EIT_normalized_dif)
    self.min = self.raw_image.min()
    self.max = self.raw_image.max()
  
    
  def get_ref_frame(self):
    meast,frames = np.shape(self.raw_image)
    resmin = -10000000
    for frame in range(frames):
      res = 0
      for meas in range(meast):
        res += self.raw_image[meas,frame]
      if res > resmin:
        framin = frame
        resmin = res
    return framin
