# Invio

Animate incidence rates on geographic projections and cummulative case counts over time with tweening.

The script broadly has three parts,
1. Compute incidence rates and case counts.
2. Read shapefile of country/region and create polygons
3. Animate incidence rates and case counts with tweening

## Install dependencies

* Python dependencies
```
pip install -r requirements.txt
```

* [ffmpeg](https://ffmpeg.org/)

Add the path to ffmpeg on line 23 of [invio.py](./scripts/invio.py)

## Output

Video output: [yearly_incidence_rate_usa.mp4](./scripts/yearly_incidence_rate_usa.mp4)

## TODO

* Generalize invio.py to accept shapefile, incidence rates and case counts
