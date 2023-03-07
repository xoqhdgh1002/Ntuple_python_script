import pandas as pd
import numpy as np
class Analysis:
	def __init__(self):
		pass

	def read_file(self,file,SB):
		data = pd.read_csv(file)
		if SB == "S":
			data = data[data['Ele_isMatchTrue']==1]
		elif SB == "B":
			data = data[data['Ele_isMatchTrue']==0]
		self.data = data.groupby("entry").head(1)

	def choose_Ecal_region(self,region):
		mask = {
		"Common":
		{'Ele_Et':lambda data: data[data['Ele_Et']>35],
		'Ele_ecalDriven':lambda data: data[data['Ele_ecalDriven']==1],
		'Ele_dPhiIn':lambda data: data[np.abs(data['Ele_dPhiIn'])<0.06],
		'Ele_dr03TkSumPt':lambda data: data[data['Ele_dr03TkSumPt']<5],
		'Ele_MissingInnerHits':lambda data: data[data['Ele_MissingInnerHits']<=1]},
		"EB":
		{'Ele_isEB':lambda data: data[data["Ele_isEB"]==1],
		'Ele_dEtaIn':lambda data: data[np.abs(data['Ele_dEtaIn'])<0.004],
		'Ele_hOverE':lambda data: data[data['Ele_hOverE']<(1/data['Ele_Esc']+0.05)],
		'Ele_E2x5OverE5x5':lambda data: data[(data['Ele_E2x5OverE5x5']>0.94)|(data['Ele_E1x5OverE5x5']>0.83)],
		'Ele_isoEmHadDepth1':lambda data:data[((data['Ele_Et']<50)&(data['Ele_isoEmHadDepth1']<(2.5+0.28*data['rho'])))|((data['Ele_pt']>=50)&(data['Ele_isoEmHadDepth1']<(2.5+0.03*(data['Ele_Et']-50)+0.28*data['rho'])))],
		'Ele_dxy':lambda data:data[np.abs(data['Ele_dxy'])<0.02]},
		"EE":
		{'Ele_isEE':lambda data: data[data["Ele_isEE"]==1],
		'Ele_dEtaIn':lambda data: data[np.abs(data['Ele_dEtaIn'])<0.006],
		'Ele_hOverE':lambda data: data[data['Ele_hOverE']<(5/data['Ele_Esc']+0.05)],
		'Ele_full5x5_sigmaIetaIeta':lambda data: data[data['Ele_full5x5_sigmaIetaIeta']<0.03],
		'Ele_isoEmHadDepth1':lambda data: data[data['Ele_isoEmHadDepth1']<(2+0.03*data['Ele_Et']+0.28*data['rho'])],
		'Ele_dxy':lambda data: data[np.abs(data['Ele_dxy'])<0.05]}}

		self.mask = {**mask['Common'],**mask[region]}

	def apply_mask(self,vars,target='Ele_pt'):
		data = self.data
		if len(vars) == 0:
			return data
		elif len(vars) > 0:
			for var in vars:
				data = self.mask[var](data)
			return data[target].values

	def get_none(self,target='Ele_pt'):
		none = self.apply_mask([],target)
		return none

	def get_N_minus_1(self,vars,target='Ele_pt'):
		masks = self.mask.keys()
		N_minus_1 = self.apply_mask([x for x in masks if x not in vars],target)
		return N_minus_1

	def get_all(self,target='Ele_pt'):
		masks = self.mask.keys()
		all = self.apply_mask(masks,target)
		return all