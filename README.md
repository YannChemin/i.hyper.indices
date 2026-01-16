# i.indices.hyper - GRASS GIS Hyperspectral/Multispectral Indices Calculator

[![GRASS GIS](https://grass.osgeo.org/logos/grass-legend.png)](https://grass.osgeo.org/)

**Calculates 86+ spectral indices from hyperspectral and multispectral imagery.** Supports automatic wavelength-based band matching with thematic organization (vegetation, pigments, water, soil, urban, stress, biochemical, materials).

## Features

- **86+ Spectral Indices** across 9 themes:
  - Vegetation (15): NDVI, EVI, SAVI, GNDVI, NDRE, MTCI...
  - Pigments (16): ARI1, CRI550, CARI, MCARI/OSAVI...
  - Metabolism (13): WI, NDWI1240, LWCI, NDII...
  - Materials (24): PLASTIC, PDI, HI, THI (plastic/oil detection)...
- **Automatic Band Matching** - provides wavelengths, gets best spectral matches
- **Thematic Organization** - vegetation, water, pigments, stress, materials, etc.
- **GRASS Native** - uses `r.mapcalc`, color tables, full integration
- **Normalization Support** - scales to 0-1 range (-n flag)
- **Robust Error Handling** - skips unavailable bands with warnings

## Installation

```bash
# Clone/download to GRASS addons directory
cd $GISBASE/etc/python/scripts  # or ~/.grass8/addons/
git clone https://www.github.com/YannChemin/i.indices.hyper
cd i.indices.hyper

# Compile
make

# Install permanently
g.extension i.indices.hyper

## Quick Start

~~~bash
# 1. List all 86+ indices by theme
i.indices.hyper -l

# 2. Calculate NDVI (default)
i.indices.hyper input=b_red,b_nir wavelengths=665,850 output_prefix=my_ndvi

# 3. All vegetation indices (normalized)
i.indices.hyper -n input=b1,b2,b3,b4,b5 \
    wavelengths=480,560,665,705,842 \
    output_prefix=forest theme=vegetation

# 4. Plastic detection (marine debris)
i.indices.hyper input=coastal_b1,b4,b8,b12 \
    wavelengths=490,665,865,1600 \
    output_prefix=marine indices=PLASTIC,PDI,FPI
~~~

## Full Usage

# Key Parameters:
    input: comma-separated raster bands
    wavelengths: corresponding wavelengths in nm
    output_prefix: prefix for output rasters
    indices or theme: specific indices or theme name

# Flags:
    -l: list indices
    -i: detailed index info
    -n: normalize to 0-1 range

# Example Themes
Theme	Indices	Applications
vegetation	15	NDVI, EVI, SAVI, MTCI, MCARI
pigments	16	ARI1, CRI550, CARI (chlorophyll)
materials	24	PLASTIC, PDI, HI (debris, oil)
water	3	NDWI, MNDWI, NDMI

# Files:
~~~bash
├── i.hyper.indices.py     # Main module (55k LOC)
├── Makefile              # GRASS build system
├── i.hyper.indices.html  # HTML manual
├── i.hyper.indices.md    # Markdown manual
└── README.md            # This file
~~~

# References

Each index includes original scientific citation. View with:

~~~bash
i.indices.hyper -li
~~~
