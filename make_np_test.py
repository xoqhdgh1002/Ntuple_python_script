from glob import glob
import uproot
import numpy as np
import awkward as ak
import pandas as pd
import math
class skimming:
	def __init__(self,files):
		if "*" in files:
			self.files = glob(files)
		else:
			self.files = [files]

	def test(self,var):
		tree = self.files[0] + ":ntuple/tree"
		data = uproot.open(tree)
		print(data[var].array())

	def make(self,vars,name,sb):
		trees = [x+":ntuple/tree" for x in self.files]
		arr = []
		cnt = 0
		for tree in uproot.iterate(trees,vars):
			cnt += 1
			print('\r', f"{cnt}/{len(trees)}", end='',flush=True)
			for var in vars:
				if len(np.shape(tree[var][0])) == 0:
					tree[var] = ak.broadcast_arrays(tree[var][:,np.newaxis],tree['Ele_pt'])[0]
				if "Ele_d" in var and "Ele_dr" not in var:
					tree[var] = np.abs(tree[var])
			arr.append(tree)
		self.ak_trees = ak.concatenate(arr)
		if sb == "S":
			self.ak_trees = self.ak_trees[self.ak_trees["Ele_isMatchTrue"]==1]
		elif sb == "B":
			self.ak_trees = self.ak_trees[self.ak_trees["Ele_isMatchTrue"]==0]
		self.ak_trees = self.ak_trees[self.ak_trees["isPVGood"]==1]
		self.ak_trees = self.ak_trees[ak.num(self.ak_trees['Ele_pt'],axis=1)>0][:,0]
		ak.to_dataframe(self.ak_trees).to_csv(f"{name}.csv",index=False)

class analysis:
	def __init__(self,csvf,region):
		self.data = pd.read_csv(csvf)
		self.region = region

	def choose_region(self):
		mask = {
		"Common":
		{'Ele_Et':self.data['Ele_Et']>35,
		'Ele_ecalDriven':self.data['Ele_ecalDriven']==1,
		'Ele_dPhiIn':self.data['Ele_dPhiIn']<0.06,
		'Ele_dr03TkSumPt':self.data['Ele_dr03TkSumPt']<5,
		'Ele_MissingInnerHits':self.data['Ele_MissingInnerHits']<=1},
		"EB":
		{'Ele_isEB':self.data["Ele_isEB"]==1,
		'Ele_dEtaIn':self.data['Ele_dEtaIn']<0.004,
		'Ele_hOverE':self.data['Ele_hOverE']<(1/self.data['Ele_Esc']+0.05),
		'Ele_E2x5OverE5x5':(self.data['Ele_E2x5OverE5x5']>0.94)|(self.data['Ele_E1x5OverE5x5']>0.83),
		'Ele_isoEmHadDepth1':((self.data['Ele_Et']<50)&(self.data['Ele_isoEmHadDepth1']<(2.5+0.28*self.data['rho'])))|((self.data['Ele_pt']>=50)&(self.data['Ele_isoEmHadDepth1']<(2.5+0.03*(self.data['Ele_Et']-50)+0.28*self.data['rho']))),
		'Ele_dxy':self.data['Ele_dxy']<0.02},
		"EE":
		{'Ele_isEE':self.data["Ele_isEE"]==1,
		'Ele_dEtaIn':self.data['Ele_dEtaIn']<0.006,
		'Ele_hOverE':self.data['Ele_hOverE']<(5/self.data['Ele_Esc']+0.05),
		'Ele_full5x5_sigmaIetaIeta':self.data['Ele_full5x5_sigmaIetaIeta']<0.03,
		'Ele_isoEmHadDepth1':self.data['Ele_isoEmHadDepth1']<(2+0.03*self.data['Ele_Et']+0.28*self.data['rho']),
		'Ele_dxy':self.data['Ele_dxy']<0.05}}

		self.mask = {**mask['Common'],**mask[self.region]}

	def apply_mask(self,masks,var):
		selected_mask = [self.mask[x] for x in masks]
		if len(masks) == 0:
			return self.data[var]
		if len(masks) > 0:
			mask = selected_mask[0]
			for i in selected_mask:
				mask = mask&i
			return self.data[var][mask]

import matplotlib.pyplot as plt
from matplotlib import colors
import mplhep as hep
def draw_histogram(datas1,datas2,title):
	total = np.concatenate(tuple(datas1.values()))
	xmin = np.min(total)
	xmax = np.max(total)
	size = 10**(int(math.log10(xmax))-2)*2
	hep.style.use("CMS")
	plt.figure(figsize=(10,10))
	for data in datas1.keys():
		plt.hist(datas1[data],bins=np.arange(0,xmax+size,size),histtype='step',linewidth=1.2,label=data+"_before")
		plt.hist(datas2[data],bins=np.arange(0,xmax+size,size),histtype='stepfilled',alpha=0.5,label=data+"_after")
	plt.xlabel(f"{title.split('_')[1]}")
	plt.ylabel(f"# of events/({size})")
	plt.legend()
	plt.grid()
	plt.yscale('log')
	plt.title(title)
	plt.savefig(f"hist/{title}_hist.png")
	plt.close()

def draw_ratio(datas1,datas2,title):
	total = np.concatenate(tuple(datas1.values()))
	xmin = np.min(total)
	xmax = np.max(total)
	hep.style.use("CMS")
	plt.figure(figsize=(10,10))
	for data in datas1.keys():
		hist1, _ = np.histogram(datas1[data],bins=np.arange(0,xmax+50,50))
		hist2, _ = np.histogram(datas2[data],bins=np.arange(0,xmax+50,50))
		ratio = hist2/hist1
		err = np.sqrt(hist2)/hist1
		plt.errorbar(np.arange(0,xmax+50,50)[:-1],ratio,yerr=[np.minimum(err,ratio),np.minimum(err,1-ratio)],fmt='.',markersize=7,label=data)
	plt.xlabel("pt[GeV]")
	plt.ylabel("ratio")
	plt.ylim((-0.1,1.1))
	plt.legend()
	plt.grid()
	plt.title(title)
	plt.savefig(f"ratio/{title}_ratio.png")
	plt.close()

def draw_2dhistogram(datas1,datas2,title1,title2,region,sample):
	total1 = datas1[title1]
	xmin1 = np.min(total1)
	xmax1 = np.max(total1)
	size1 = 10**(int(math.log10(xmax1))-2)*2
	total2 = datas2[title2]
	xmin2 = np.min(total2)
	xmax2 = np.max(total2)
	size2 = 10**(int(math.log10(xmax2))-2)*2
	hep.style.use("CMS")
	fig = plt.figure(figsize=(10,10))
	h = plt.hist2d(datas1[title1],datas2[title2],[np.arange(0,xmax1+size1,size1),np.arange(0,xmax2+size2,size2)],norm = colors.LogNorm())
	cur_ax = plt.gca()
	fig.colorbar(h[3],ax=cur_ax)
	plt.xlabel(f"{title1}")
	plt.ylabel(f"{title2}")
	plt.legend()
	plt.grid()
	plt.title(f"{title1} - {title2} _ {region}")
	plt.savefig(f"2dhist/{title1}-{title2}_{region}_{sample}_2dhist.png")
	plt.close()
'''
for i in ["EB","EE"]:
	arr = {}
	for j in ["zp_3000","zp_4000","tt"]:
		a = analysis(f'csv/{j}.csv',i)
		a.choose_region()
		mask_list = list(a.mask.keys())
		for k in [x for x in mask_list if "Ele_isE" not in x]:
			l = [x for x in mask_list if x != k]
			if k not in arr.keys():
				arr[k] = {}
			arr[k][j] = a.apply_mask(l,'Ele_pt').values
		if "none" not in arr.keys():
			arr["none"] = {}
		arr["none"][j] = a.apply_mask([f'Ele_is{i}'],'Ele_pt').values
		if "all" not in arr.keys():
			arr["all"] = {}
		arr["all"][j] = a.apply_mask(mask_list,'Ele_pt').values
	for x in [x for x in mask_list if "Ele_isE" not in x]:
		draw_histogram(arr[x],arr["all"],f"{x}_{i}")
		draw_ratio(arr[x],arr["all"],f"{x}_{i}")
	draw_histogram(arr["none"],arr["all"],f"all_{i}")
	draw_ratio(arr["none"],arr["all"],f"all_{i}")
'''
for i in ["EB","EE"]:
	for j in ["zp_3000","zp_4000","tt"]:
		a = analysis(f'csv/{j}.csv',i)
		a.choose_region()
		mask_list = list(a.mask.keys())
		arr = {}
		for k in ["Ele_pt","Ele_hOverE"]:
			l = [x for x in mask_list if x != "Ele_hOverE"]
			#l = mask_list
			arr[k] = a.apply_mask(l,k).values
		draw_2dhistogram(arr,arr,"Ele_pt","Ele_hOverE",i,j)
'''

for i in ["EB","EE"]:
	arr = {}
	all = {}
	for j in ["zp_3000","zp_4000","tt"]:
		a = analysis(f'csv/{j}.csv',i)
		a.choose_region()
		mask_list = list(a.mask.keys())
		for k in [x for x in mask_list if "Ele_isE" not in x]:
			l = [x for x in mask_list if x != k]
			if k not in arr.keys():
				arr[k] = {}
			arr[k][j] = a.apply_mask(l,k).values
			if j == "zp_3000":
				arr[k][j] = arr[k][j][a.apply_mask(l,'Ele_pt').values>1500]
			if j == "zp_4000":
				arr[k][j] = arr[k][j][a.apply_mask(l,'Ele_pt').values>2000]
			if k not in all.keys():
				all[k] = {}
			all[k][j] = a.apply_mask(mask_list,k).values
			if j == "zp_3000":
				all[k][j] = all[k][j][a.apply_mask(mask_list,'Ele_pt').values>1500]
			if j == "zp_4000":
				all[k][j] = all[k][j][a.apply_mask(mask_list,'Ele_pt').values>2000]
	for x in [x for x in mask_list if "Ele_isE" not in x]:
		draw_histogram(arr[x],all[x],f"{x}_{i}")

	
'''
'''
if __name__ == "__main__":
	import sys
	var_list = [
	'Ele_isMatchTrue',
	'Ele_isEE3',
	'Ele_isEB',
	'Ele_isEE',
	'Ele_Et',
	'Ele_ecalDriven',
	'Ele_dEtaIn',
	'Ele_dPhiIn',
	'Ele_hOverE',
	'Ele_full5x5_sigmaIetaIeta',
	'Ele_E1x5OverE5x5',
	'Ele_E2x5OverE5x5',
	'Ele_isoEmHadDepth1',
	'Ele_MissingInnerHits',
	'Ele_dxy',
	'Ele_Esc',
	'Ele_dr03TkSumPt',
	'rho',
	'Ele_pt',
	'Ele_passHEEPId',
	'isPVGood'
	]
	a = skimming(sys.argv[1]+"*")
#	a.make(var_list,sys.argv[2],sys.argv[3])
	a.test("isPVGood")
'''
