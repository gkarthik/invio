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
import matplotlib.dates as mdates
from dateutil.relativedelta import *
import datetime

# Add path to ffmpeg here
plt.rcParams['animation.ffmpeg_path'] = u'/usr/bin/ffmpeg'

##############
# Setup data #
##############

disease = "ZIKA"              # Chikungunya or ZIKA

df = pd.read_csv("../data/DENV_CHKV_ZIKV_cleaned_no_metadata_mapped.csv", index_col = 0)
disease_df = df[df["Suspected.Pathogen"] == disease]
disease_df["Date.of.Onset.of.Symptoms"] = pd.to_datetime(disease_df["Date.of.Onset.of.Symptoms"], format="%m/%d/%Y")

# Groupby Year - Month to make transitions slower
disease_df["Year-Month"] = pd.to_datetime(disease_df["Date.of.Onset.of.Symptoms"].dt.strftime("%Y-%m-15"), format="%Y-%m")

_times = [disease_df["Year-Month"].min()]
while _times[-1] < disease_df["Year-Month"].max():
    _times.append(_times[-1] + relativedelta(months=1))

# Group by individual days
# disease_df["Year-Month"] = pd.to_datetime(disease_df["Date.of.Onset.of.Symptoms"])
# _times = [disease_df["Year-Month"].min()]
# while _times[-1] < disease_df["Year-Month"].max():
#     _times.append(_times[-1] + relativedelta(days=1))

grouped_df = disease_df.groupby(["Municipality_map", "Year-Month"]).count()
grouped_df["count"] = grouped_df["Suspected.Pathogen"]
country_cases = grouped_df.reset_index().groupby("Year-Month").count()["count"]

grouped_df["count"] = np.log10(grouped_df["count"])

_max = grouped_df["count"].max()
_min = grouped_df["count"].min()
_cmap = plt.get_cmap('Reds')
cNorm  = colors.Normalize(vmin=_min, vmax=_max)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=_cmap)

# DMR
f = plt.figure(figsize=(13,10))
gs = gridspec.GridSpec(3, 2,width_ratios=[1,0.02], height_ratios=[0.2,1,0.3])
yearax = plt.subplot(gs[0, 0])
mapax = plt.subplot(gs[1, 0])
scaleax = plt.subplot(gs[1, 1])
caseax = plt.subplot(gs[2,0])

in_projection = Proj(init = "epsg:32619")
out_projection = Proj(init = "epsg:3857")

yearax.axis("off")
mapax.axis("off")

##################
# Plot shapefile #
##################

poly_names = []
patches = []
shp = fiona.open("../shapefiles/muncenso2010/MUNCenso2010.shp")

ll = transform(in_projection, out_projection, shp.bounds[0], shp.bounds[1])
ur = transform(in_projection, out_projection, shp.bounds[2], shp.bounds[3])

for feat in shp:
    if feat["geometry"]["type"] == "Polygon":
        coords = feat["geometry"]["coordinates"][0]
        coords = transform(in_projection, out_projection, [i[0] for i in coords], [i[1] for i in coords])
        poly = Polygon(np.dstack(coords)[0])
        poly.set_facecolor("#FFFFFF")
        patches.append(poly)
        r = feat["properties"]["TOPONIMIA"]
        poly_names.append(r)
    else:
        for k in feat["geometry"]["coordinates"]:
            coords = k[0]
            coords = transform(in_projection, out_projection, [i[0] for i in coords], [i[1] for i in coords])
            poly = Polygon(np.dstack(coords)[0])
            poly.set_facecolor("#FFFFFF")
            patches.append(poly)
            r = feat["properties"]["TOPONIMIA"]
            poly_names.append(r)
    print(feat["properties"]["TOPONIMIA"])

poly_map_collection = mapax.add_collection(PatchCollection(patches, facecolor = "#ECECEC", edgecolor="#333333"))

cb = ColorbarBase(scaleax, cmap=_cmap, norm = cNorm, orientation='vertical')
cb.formatter   = StrMethodFormatter(r"$10^{{{x:,.1f}}}$")
cb.update_ticks()
scaleax.yaxis.set_label_position("left")
scaleax.set_ylabel("Case counts")

mapax.set_xlim([ll[0], ur[0]])
mapax.set_ylim([ll[1], ur[1]])
mapax.axis('equal')
# Year
ytext = yearax.annotate("1999", (0,0), size = 30, family="sans-serif")

# Case Count Polygon
case_polygon = Polygon([(0,0)], closed=True, facecolor='#4c4cff', fill=True, alpha = 0.2)
case_polygon = caseax.add_patch(case_polygon)

# # Case Counts
country_cases.plot(ax=caseax, label="Cases", color="#4c4cff")
caseax.axvline(country_cases.index.min(), lw=0.75, color="#008000", linestyle="--")
axvl = caseax.axvline(country_cases.index.min(), lw=1, color="#000000")
caseax.spines['right'].set_visible(False)
caseax.spines['top'].set_visible(False)
caseax.legend(loc='upper left')

country_cases = country_cases.reindex(_times).fillna(0)

#########################
# Animate with tweening #
#########################

def animate(i):
    if i%100 == 0:
        print(i)
    index = i/(_frames/(len(_times) - 1))
    _time = _times[int(index)]
    ytext.set_text(_time.strftime("%Y-%m"))
    d = np.array(axvl.get_data())
    _x = _times[0]+((i/_frames)*(_times[-1]-_times[0]))
    d[0] = [_x, _x]
    axvl.set_data(d)
    x = [0]
    y = [0]
    y.extend(country_cases.ix[country_cases.index<=_time].tolist())
    x.extend([mdates.date2num(j) for j in country_cases.ix[country_cases.index<=_time].index.tolist()])
    # y.append(country_cases.ix[_time])
    if _time != _times[-1]:
        y.append((country_cases.ix[_times[int(index)+1]] - country_cases.ix[_time])*(mdates.date2num(_x) - int(mdates.date2num(_x))) + country_cases.ix[_time])
    else:
        y.append(country_cases.ix[_time])
    x.append(mdates.date2num(_x))
    y.append(0)
    x.append(mdates.date2num(_x))
    case_polygon.set_xy(list(zip(x,y)))
    poly_colors = []
    for poly_name, i in zip([i.lower() for i in poly_names], range(0, len(patches))):
        if poly_name not in grouped_df.index.get_level_values(0).tolist():
            poly_colors.append("#FFFFFF")
            continue
        prevVal = _min
        if _time in grouped_df.loc[poly_name].index:
            prevVal = grouped_df.loc[(poly_name, _time),"count"]
        if prevVal == -float("inf"):
            prevVal = _min
        val = prevVal
        if _time != _times[-1]:
            nextVal = 0
            if _times[int(index)+1] in grouped_df.loc[poly_name].index:
                nextVal = grouped_df.loc[(poly_name, _times[int(index)+1]),"count"]
            if nextVal == -float("inf"):
                nextVal = _min
            val = prevVal + ((nextVal - prevVal)*(index - int(index)))
        poly_colors.append(scalarMap.to_rgba(val))
    poly_map_collection.set_facecolor(poly_colors)
    return  (poly_map_collection, ytext,)

print("Generating video ... ")
plt.suptitle("{} case counts".format(disease.title()), fontsize = 30)
_frames = (len(_times) - 1) * 10
# _frames = 500
anim = animation.FuncAnimation(plt.gcf(), animate, frames=_frames + 1, interval = 20, blit=True, repeat=False)

gs.tight_layout(f)

FFwriter = animation.FFMpegWriter(fps=45)
anim.save('{}_DO.mp4'.format(disease),writer=FFwriter, dpi = 300)
plt.clf()
plt.close()

print("Finished generation of video")
