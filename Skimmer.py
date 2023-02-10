from glob import glob
import uproot
import awkward as ak
import sys
def make_tree_list(files):
    tree_path = ":ntuple/tree"
    trees = [x+tree_path for x in glob(files+"*")]
    return trees

def skim_ntuples(trees,vars):
    arr = []
    cnt = 0
    for tree in uproot.iterate(trees,vars):
        print('\r', f"{cnt}/{len(trees)}", end='',flush=True)
        arr.append(tree)
        cnt+=1
    ak_trees = ak.concatenate(arr)
    return ak_trees

if __name__ == "__main__":
    vars = [
	'isPVGood',
	'rho',
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
	'Ele_pt',
	'Ele_passHEEPId',
	'Ele_etaSC'
    ]

    trees = make_tree_list(sys.argv[1])
    print(trees)
    skimed_data = skim_ntuples(trees,vars)
    ak.to_dataframe(skimed_data).to_csv(f'{sys.argv[2]}.csv')
#   ak.to_json(skimed_data,"test.json",num_indent_spaces=4)