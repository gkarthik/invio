import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 20, "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 12})

import shapely
import fiona
import matplotlib.colors as colors
from matplotlib.patches import Polygon
from matplotlib.colorbar import ColorbarBase
import matplotlib.cm as cmx
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
from matplotlib.collections import PatchCollection
from pyproj import Proj, transform

from matplotlib import animation
import matplotlib.gridspec as gridspec

import datetime

plt.rcParams['animation.ffmpeg_path'] = u'/usr/bin/ffmpeg'

##############
# Setup data #
##############

pop = pd.read_csv("../data/nst-est2017-01.csv")
pop = pop.set_index("State")
for i in pop.columns:
    pop[i] = pop[i].astype(str).apply(lambda x: x.replace(",", ""))
    pop[i] = pop[i].astype(int)
p = pop.mean(axis = 1)

df = pd.read_csv("../data/West-Nile-virus-disease-cases-reported-to-CDC-by-state_1999-2016_09292017.csv")
df = df.set_value(8, "State ", "District of Columbia")
df = df.set_index("State ")
df = df.drop("Total ").drop("Total ", axis = 1)
df = df.drop("Puerto Rico ")
df.index = df.index.str.rstrip()
for i in df.columns:
    df[i] = df[i].astype(str).apply(lambda x: x.replace(",", ""))
    df[i] = df[i].astype(int)

# Neuroinvasive disease
ndf = pd.read_csv("../data/West-Nile-virus-neuroinvasive-disease-cases-reported-to-CDC-by-state_1999-2016_09292017.csv")
ndf = ndf.set_value(8, "State ", "District of Columbia")
ndf = ndf.set_index("State ")
ndf = ndf.drop("Total ", axis = 1)
ndf = ndf.drop("Puerto Rico ")
ndf.index = ndf.index.str.rstrip()
for i in ndf.columns:
    ndf[i] = ndf[i].astype(str).apply(lambda x: x.replace(",", ""))
    ndf[i] = ndf[i].astype(int)

ndf = ndf.transpose()
ndf.index = ndf.index.astype(int)
# Get total case counts
tdf = df.transpose()
tdf.index = tdf.index.astype(int)

for i in df.index:
    df.ix[i] = (df.ix[i]/p[i.rstrip()]) * 100000

df = np.log10(df)

_min = float("inf")
for i in df.index:
    for c in df.columns:
        if df[c][i] == float("-inf"):
            continue
        if _min > df[c][i]:
            _min = df[c][i]

_max = df.max().max()
_cmap = plt.get_cmap('Reds')
cNorm  = colors.Normalize(vmin=_min, vmax=_max)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=_cmap)

# USA
f = plt.figure(figsize=(13,10))
gs = gridspec.GridSpec(3, 2,width_ratios=[1,0.02], height_ratios=[0.2,1,0.3])
yearax = plt.subplot(gs[0, 0])
mapax = plt.subplot(gs[1, 0])
scaleax = plt.subplot(gs[1, 1])
caseax = plt.subplot(gs[2,0])

in_projection = Proj(init = "epsg:4269")
out_projection = Proj(init = "epsg:3857")

yearax.axis("off")
mapax.axis("off")
ll_lng = -128.194915
ll_lat = 23.428521
ur_lat = 51.087372
ur_lng = -51.116799
ll = transform(in_projection, out_projection, ll_lng, ll_lat)
ur = transform(in_projection, out_projection, ur_lng, ur_lat)

##################
# Plot shapefile #
##################

poly_names = []
patches = []
shp = fiona.open("../shapefiles/cb_2016_us_state_20m/cb_2016_us_state_20m.shp")
for feat in shp:
    if feat["geometry"]["type"] == "Polygon":
        coords = feat["geometry"]["coordinates"][0]
        coords = transform(in_projection, out_projection, [i[0] for i in coords], [i[1] for i in coords])
        poly = Polygon(np.dstack(coords)[0])
        poly.set_facecolor("#FFFFFF")
        patches.append(poly)
        r = feat["properties"]["NAME"]
        poly_names.append(r)
    else:
        for k in feat["geometry"]["coordinates"]:
            coords = k[0]
            coords = transform(in_projection, out_projection, [i[0] for i in coords], [i[1] for i in coords])
            poly = Polygon(np.dstack(coords)[0])
            poly.set_facecolor("#FFFFFF")
            patches.append(poly)
            r = feat["properties"]["NAME"]
            poly_names.append(r)
    print(feat["properties"]["NAME"])

poly_map_collection = mapax.add_collection(PatchCollection(patches, facecolor = "#ECECEC", edgecolor="#333333"))

cb = ColorbarBase(scaleax, cmap=_cmap, norm = cNorm, orientation='vertical')
# scaleax.yaxis.set_major_locator(MaxNLocator(integer=True))
cb.formatter   = StrMethodFormatter(r"$10^{{{x:,.1f}}}$")
cb.update_ticks()
scaleax.yaxis.set_label_position("left")
scaleax.set_ylabel("Incidence rate per 100,000 humans")

mapax.set_xlim([ll[0], ur[0]])
mapax.set_ylim([ll[1], ur[1]])
# Year
ytext = yearax.annotate("1999", (0,0), size = 30, family="sans-serif")

# Case Count Polygon
case_polygon = Polygon([(0,0)], closed=True, facecolor='#4c4cff', fill=True, alpha = 0.2)
case_polygon = caseax.add_patch(case_polygon)
ncase_polygon = Polygon([(0,0)], closed=True, facecolor='#ff4c4c', fill=True, alpha = 0.2)
ncase_polygon = caseax.add_patch(ncase_polygon)

# Case Counts
tdf.sum(axis = 1).cumsum().plot(ax=caseax, label="Cumm Human", color="#4c4cff")
ndf.sum(axis = 1).cumsum().plot(ax=caseax, label="Cumm Neuroinvasive Human", color="#ff4c4c")
caseax.axvline(2003, lw=0.75, color="#008000", linestyle="--")
axvl = caseax.axvline(1999, lw=1, color="#000000")
caseax.spines['right'].set_visible(False)
caseax.spines['top'].set_visible(False)
caseax.set_ylim([0, 46050])
caseax.set_yticks([0, 20000, 40000])
caseax.set_xlim([1999,2016])
caseax.set_xticks(list(range(1999,2017)))
caseax.legend(loc='upper left')

_years = df.columns.astype(int).tolist()
t = tdf.sum(axis=1).cumsum()
n = ndf.sum(axis=1).cumsum()

#########################
# Animate with tweening #
#########################

def animate(i):
    index = i/(_frames/(len(_years) - 1))
    year = _years[int(index)]
    ytext.set_text(year)
    d = np.array(axvl.get_data()).astype(float)
    _x = _years[0]+((i/_frames)*(_years[-1]-_years[0]))
    d[0] = [_x, _x]
    axvl.set_data(d)
    x = [0]
    y = [0]
    y.extend(t[t.index<=year].tolist())
    x.extend(t[t.index<=year].index.tolist())
    if year +1 <= _years[-1]:
        y.append((t[year+1] - t[year])*(_x - int(_x)) + t[year])
    else:
        y.append(t[year])
    x.append(_x)
    y.append(0)
    x.append(_x)
    case_polygon.set_xy(list(zip(x,y)))
    y = [0]
    y.extend(n[n.index<=year].tolist())
    if year +1 <= _years[-1]:
        y.append((n[year+1] - n[year])*(_x - int(_x)) + n[year])
    else:
        y.append(n[year])
    y.append(0)
    ncase_polygon.set_xy(list(zip(x,y)))
    poly_colors = []
    for poly_name, i in zip(poly_names, range(0, len(patches))):
        if poly_name not in df.index.tolist():
            continue
        prevVal = df.loc[poly_name,str(year)]
        if prevVal == -float("inf"):
            prevVal = _min
        val = prevVal
        if year != _years[-1]:
            nextVal = df.loc[poly_name,str(year+1)]
            if nextVal == -float("inf"):
                nextVal = _min
            val = prevVal + ((nextVal - prevVal)*(index - int(index)))
        poly_colors.append(scalarMap.to_rgba(val))
    poly_map_collection.set_facecolor(poly_colors)
    return  (poly_map_collection, ytext,axvl,case_polygon,ncase_polygon,)

print("Generating video ... ")
plt.suptitle("Yearly Human Incidence Rates of West Nile Virus", fontsize = 30)
_frames = (len(_years) - 1) * 45
anim = animation.FuncAnimation(plt.gcf(), animate, frames=_frames + 1, interval = 20, blit=True, repeat=False)

gs.tight_layout(f)

FFwriter = animation.FFMpegWriter(fps=45)
anim.save('yearly_incidence_rate_usa.mp4',writer=FFwriter, dpi = 300)
plt.clf()
plt.close()

print("Finished generation of video")
