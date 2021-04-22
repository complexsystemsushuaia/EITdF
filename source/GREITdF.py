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
# EITdata.EIT_data_meas
# EITdata.EIT_data_meas_passive_electrodes

# EITdata.EIT_data_meas_active_electrodes
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
# RawImage.RawImage(rec_model,eitdata)
# RawImage.CalcRawImage(reffframe) #refframe is the number of frame taken as reference. It is called by the constructor.
# RawImage.GetRefFrame()           #finds the minimum conductivity frame to act as a reference. Called by constructor.
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
        return


class EITdata:
    def __init__(self,environ,inputfile):
        self.inputfile                 = inputfile
        self.environment               = environ
        # Load input
        self.EIT_data_meas = io.loadmat(self.environment.inputdir+"/"+self.inputfile)['FRAMES']
        self.set_data()

     
    def set_data(self):
        [NrMeas, NrFrames] = np.shape(self.EIT_data_meas)
        self.NrFrames = NrFrames
        self.NrMeas   = NrMeas
        if (NrMeas == 256):
            self.measure_active_electrodes = 1
        else:
            self.measure_active_electrodes = 0
        self.EIT_data_meas_passive_electrodes,self.EIT_data_meas_active_electrodes = self.separate_dataframe_data()
        #
        # Order ciclically active electrodes: 15,0,1; 0,1,2; ... ;13,14,15.
        # Dim is 3*16, each triplet is an injection frame.
        #
        if (self.measure_active_electrodes == 1):
            for k in range(NrFrames):
                aux = self.EIT_data_meas_active_electrodes[3*k*16]
                self.EIT_data_meas_active_electrodes[3*k*16] = self.EIT_data_meas_active_electrodes[3*k*16+2]
                self.EIT_data_meas_active_electrodes[3*k*16+2] = aux
                aux = self.EIT_data_meas_active_electrodes[3*k*16+45]
                self.EIT_data_meas_active_electrodes[3*k*16+45] = self.EIT_data_meas_active_electrodes[3*k*16+47]
                self.EIT_data_meas_active_electrodes[3*k*16+47] = aux
        return
		    
    
    def separate_dataframe_data(self):
        if self.measure_active_electrodes == 0:
            return self.EIT_data_meas, []
        active        = []
        passive_meas  = []
        electrode_pairs = [14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1]
        for frame in range(self.NrFrames):
            passive_frame = []
            for meas_in_frame in range(self.NrMeas):
                inj      = meas_in_frame//16      # Injection frame
                meas_in_inj = meas_in_frame-16*inj    # Measurement within injection frame between 0 and 15
                skipped  = [electrode_pairs[inj+m+1] for m in range(3)]
                if meas_in_inj in skipped:
                    active.append(self.EIT_data_meas[meas_in_frame,frame])
                else:
                    passive_frame.append(self.EIT_data_meas[meas_in_frame,frame])
                    #print(meas_in_frame, frame)
            if (len(passive_frame) > 0):
                passive_meas.append(passive_frame)
        return np.array(passive_meas).transpose(), active    
  
    def get_normalized_dataframes(self,ref_dataframe):
        rows, cols = np.shape(self.EIT_data_meas_passive_electrodes)
        if (ref_dataframe+1  > cols) or (ref_dataframe < 0):
            print("WARNING: reference frame out of range. Using 0 instead")
            ref_dataframe = 0
        self.ref_frame = ref_dataframe
        self.EIT_normalized_dif = (self.EIT_data_meas_passive_electrodes.transpose() / self.EIT_data_meas_passive_electrodes.transpose()[ref_dataframe] - 1).transpose()
        return
    
    def estimate_contact_resistances(self,injected_currents):
        if (self.measure_active_electrodes != 1):
            print("ERROR: Contact resistances cannot be estimated. No data on active electrodes.")
            return
        contact_resistances = []
        for k in range(16*self.NrFrames):
            v = abs(self.EIT_data_meas_active_electrodes[3*k])+abs(self.EIT_data_meas_active_electrodes[3*k-1])/2
            contact_resistances.append(v/injected_currents[k])
        self.contact_resistances = contact_resistances
        return
        
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
        [Nmeas, Nframes] = np.shape(self.EIT_data_meas)
        # Generate frames of noise
        uncorrelated_noise = [np.random.normal(0,1,Nframes) for k in range(Nmeas)]
        uncorrelated_noise_frames = np.vstack(uncorrelated_noise)
        noise_frames = T.dot(uncorrelated_noise_frames)
        # Get renormalization factor for noise to the specified SNR
        noise_level = np.linalg.norm(self.EIT_data_meas)/(SNR * np.linalg.norm(noise_frames))
        # Add noise frames to data frames as returned data
        self.EIT_data_meas = self.EIT_data_meas + noise_frames * noise_level
        self.set_data()
        return

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
        # Compute inverse matrix
        self.set_hyperparameter(hyperparameter_fraction)
        # Triangulation (tri and nodes)
        self.tri                  = io.loadmat(environ_options.modelsdir+"/"+model_name+"/elems.mat")['elems']-1
        self.nodes                = (io.loadmat(environ_options.modelsdir+"/"+model_name+"/nodes.mat")['nodes']).transpose()
        self.boundary             = io.loadmat(environ_options.modelsdir+"/"+model_name+"/boundary.mat")["boundary_coords"]
        self.electrode_positions  = io.loadmat(environ_options.modelsdir+"/"+model_name+"/electrode_positions.mat")["electrode_positions"]
        
    def set_hyperparameter(self,hyperparameter_fraction):
        self.hyperparameter_fraction = hyperparameter_fraction
        # Compute inverse matrix
        w = self.hyperparameter * self.hyperparameter_fraction
        if (self.IsTijonov):
            m_matrix            = self.J.dot(self.J.transpose()) + self.Sn*(w**2)
            self.rm             = self.J.transpose().dot(np.linalg.inv(m_matrix))
            print("TIJONOV")
        else:
            m_matrix            = self.y_matrix.dot(self.y_matrix.transpose()) + self.Sn*(w**2)
            pjt_matrix          = self.d_matrix.dot(self.y_matrix.transpose())
            self.rm             = pjt_matrix.dot(np.linalg.inv(m_matrix))
            print("GREIT")
        return

class RawImage:
    def __init__(self,rec_model,eitdata):
        self.rec_model = rec_model
        self.eitdata   = eitdata
        self.CalcRawImage(0)
        
    def CalcRawImage(self,refframe):
        self.eitdata.get_normalized_dataframes(refframe)
        self.RawImage=np.matmul(self.rec_model.rm,self.eitdata.EIT_normalized_dif)
        self.min = self.RawImage.min()
        self.max = self.RawImage.max()
        return
    
    def GetRefFrame(self):
        meast,frames = np.shape(self.RawImage)
        resmin = -10000000
        for frame in range(frames):
            res = 0
            for meas in range(meast):
                res += self.RawImage[meas,frame]
            if res > resmin:
                framin = frame
                resmin = res
        return framin
    
 #   def ShowFrame(self,frame_number):
 #       color=self.RawImage[:,frame_number]             #color=delta_normalizado[:,delta]
 #       color=np.array([[c,c] for c in color]).ravel()  #color=np.array([[c,c] for c in color]).ravel()
 #       plt.tripcolor(self.rec_model.nodes[0],self.rec_model.nodes[1],self.rec_model.tri,facecolors=color,cmap = cm.PuOr)
 #       #plt.tripcolor(nod[0],nod[1],tri,facecolors=color)
 #       return
   
   
   
