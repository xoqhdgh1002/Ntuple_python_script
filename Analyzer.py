import pandas as pd
import sys
def read_file(file):
    return pd.read_csv(file,index_col=False)

def choose_Ecal_region(region):
	mask = {
	"Common":
	{'Ele_Et':lambda data: data[data['Ele_Et']>35],
	'Ele_ecalDriven':lambda data: data[data['Ele_ecalDriven']==1],
	'Ele_dPhiIn':lambda data: data[data['Ele_dPhiIn']<0.06],
	'Ele_dr03TkSumPt':lambda data: data[data['Ele_dr03TkSumPt']<5],
	'Ele_MissingInnerHits':lambda data: data[data['Ele_MissingInnerHits']<=1]},
	"EB":
	{'Ele_isEB':lambda data: data[data["Ele_isEB"]==1],
	'Ele_dEtaIn':lambda data: data[data['Ele_dEtaIn']<0.004],
	'Ele_hOverE':lambda data: data[data['Ele_hOverE']<(1/data['Ele_Esc']+0.05)],
	'Ele_E2x5OverE5x5':lambda data: data[(data['Ele_E2x5OverE5x5']>0.94)|(data['Ele_E1x5OverE5x5']>0.83)],
	'Ele_isoEmHadDepth1':lambda data:data[((data['Ele_Et']<50)&(data['Ele_isoEmHadDepth1']<(2.5+0.28*data['rho'])))|((data['Ele_pt']>=50)&(data['Ele_isoEmHadDepth1']<(2.5+0.03*(data['Ele_Et']-50)+0.28*data['rho'])))],
	'Ele_dxy':lambda data:data[data['Ele_dxy']<0.02]},
	"EE":
	{'Ele_isEE':lambda data: data[data["Ele_isEE"]==1],
	'Ele_dEtaIn':lambda data: data[data['Ele_dEtaIn']<0.006],
	'Ele_hOverE':lambda data: data[data['Ele_hOverE']<(5/data['Ele_Esc']+0.05)],
	'Ele_full5x5_sigmaIetaIeta':lambda data: data[data['Ele_full5x5_sigmaIetaIeta']<0.03],
	'Ele_isoEmHadDepth1':lambda data: data[data['Ele_isoEmHadDepth1']<(2+0.03*data['Ele_Et']+0.28*data['rho'])],
	'Ele_dxy':lambda data: data[data['Ele_dxy']<0.05]}}

	return {**mask['Common'],**mask[region]}

def apply_mask(data,masks,vars):
	if len(masks) == 0:
		return data
	elif len(masks) > 0:
		for var in vars:
			data = masks[var](data)
		return data

if __name__ == "__main__":
	file_list = {"ttbar":read_file("tt.csv"),"Zprime_M-3000":read_file("zp3000.csv")}
	regions = ["EB","EE"]
	vars = [
		'Ele_Et'
	]

	for region in regions:
		masks = choose_Ecal_region(region)
		for file_name in file_list.keys():
			for mask in masks.keys():
				vars = [x for x in masks.keys() if "isE" not in x]
				print(apply_mask(file_list[file_name],vars))