import matplotlib.pyplot as plt
import numpy as npy

from scipy import io
import math
import sys
import os
from GREITdFv2 import EnvironOpt,RawImage,EITdata
basedir   = dir=os.getcwdb().decode("utf-8")+'..'
envi = EnvironOpt(basedir)
print(envi.bindir)
print(envi.modelsdir)
#datitos = EITdata(envi,"simuresp_allelec_snr500.mat")
#datitos = EITdata(envi,"simuresp_mistery_nomeas.mat",0)
#datitos = EITdata(envi,"simuresp_ztest.mat",1)
datitos = EITdata(envi,"simuresp_mistery_meas.mat",1)
datitos.get_normalized_dataframes(2)
#print(datitos.EIT_data_meas_active_electrodes)
plt.plot([abs(x) for x in datitos.EIT_data_meas_active_electrodes])
recmodel = GREITdF.ReconstructionModel("male_16el_01",0.02,envi)
recmodel.set_hyperparameter(0.005)
estaimagen = GREITdF.RawImage(recmodel,datitos)
estaimagen.CalcRawImage(estaimagen.GetRefFrame())
print(np.shape(estaimagen.RawImage))
print(estaimagen.eitdata.ref_frame)
estaimagen.ShowFrame(10)
datitos.estimate_contact_resistances([0.001 for k in range(16*datitos.NrFrames)])
print(datitos.contact_resistances)