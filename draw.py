import matplotlib.pyplot as plt
from matplotlib import colors
import mplhep as hep
import numpy as np
import math
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