# i.hyper.indices - GRASS GIS Hyperspectral/Multispectral Indices Calculator

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
git clone https://www.github.com/YannChemin/i.hyper.indices
cd i.hyper.indices
```

# Compile
```bash
make
```

# Install permanently
```bash
g.extension i.hyper.indices
```

## Quick Start

```bash
# 1. List all 86+ indices by theme
i.hyper.indices -l

# 2. Calculate NDVI (default)
i.hyper.indices input=b_red,b_nir wavelengths=665,850 output_prefix=my_ndvi

# 3. All vegetation indices (normalized)
i.hyper.indices -n input=b1,b2,b3,b4,b5 \
    wavelengths=480,560,665,705,842 \
    output_prefix=forest theme=vegetation

# 4. Plastic detection (marine debris)
i.hyper.indices input=coastal_b1,b4,b8,b12 \
    wavelengths=490,665,865,1600 \
    output_prefix=marine indices=PLASTIC,PDI,FPI
```

## Full Usage

# Key Parameters:
```bash
    input: comma-separated raster bands
    wavelengths: corresponding wavelengths in nm
    output_prefix: prefix for output rasters
    indices or theme: specific indices or theme name
```

# Flags:
```bash
    -l: list indices
    -i: detailed index info
    -n: normalize to 0-1 range
```

# Example Themes
|Theme	|Indices	|Applications|
|----------|----------|----------|
|vegetation	|15	|NDVI, EVI, SAVI, MTCI, MCARI|
|pigments	|16	|ARI1, CRI550, CARI (chlorophyll)|
|materials	|24	|PLASTIC, PDI, HI (debris, oil)|
|water	|3	|NDWI, MNDWI, NDMI|

# Files:
```bash
├── i.hyper.indices.py     # Main module (55k LOC)
├── Makefile              # GRASS build system
├── i.hyper.indices.html  # HTML manual
├── i.hyper.indices.md    # Markdown manual
└── README.md            # This file
```

# References

Each index includes original scientific citation. View with:

```bash
i.hyper.indices -li
```
# Exhaustive list to date
```bash
grass --exec i.hyper.indices -l
======================================================================
AVAILABLE SPECTRAL INDICES
======================================================================

BIOCHEMICAL INDICES (5):
----------------------------------------------------------------------
CAI - Cellulose Absorption Index
NDLI - Normalized Difference Lignin Index
PRI - Photochemical Reflectance Index
PSRI - Plant Senescence Reflectance Index
SIPI - Structure Insensitive Pigment Index

MATERIALS INDICES (19):
----------------------------------------------------------------------
ASPHALT - Asphalt/Bitumen Index
COAL - Coal/Carbon Index
CONCRETE - Concrete Detection Index
FERRIC - Ferric Iron Index
FERROUS - Ferrous Iron Index
FPI - Floating Plastic Index
GOSSAN - Gossan Index - Weathered sulfides
HI - Hydrocarbon Index - Oil & gas detection
LATERITE - Laterite Index - Iron-rich materials
MPDI - Marine Plastic Detection Index
NDPI - Normalized Difference Plastic Index
OHI - Oil and Hydrocarbon Index
PAINT - Paint Detection Index (synthetic coatings)
PDI - Plastic Debris Index
PLASTIC - Plastic Detection Index
RSWIR - Reversed SWIR Index - Marine debris
SINDEX - S-Index - Soil/sediment composition
THI - Tentative Hydrocarbon Index
TPI - Tar/Petroleum Index

METABOLISM INDICES (13):
----------------------------------------------------------------------
CARTER1 - Carter Index 1 - Stress
CARTER2 - Carter Index 2 - Stress
DATT - Datt Index - Leaf pigment
GMI - Gamon Index - Photosynthetic efficiency
LWCI - Leaf Water Content Index
NDII - Normalized Difference Infrared Index
NDWI1240 - Normalized Difference Water Index 1240
NDWI2130 - Normalized Difference Water Index 2130
NPCI - Normalized Pigment Chlorophyll Index
NPQI - Normalized Phaeophytinization Index
RGRI - Red-Green Ratio Index
SRWI - Simple Ratio Water Index
WI - Water Index - Leaf water content

PIGMENTS INDICES (15):
----------------------------------------------------------------------
ARI1 - Anthocyanin Reflectance Index 1
ARI2 - Anthocyanin Reflectance Index 2
CARI - Chlorophyll Absorption Ratio Index
CRI550 - Carotenoid Reflectance Index 550
CRI700 - Carotenoid Reflectance Index 700
DD - Double Difference Index - Chlorophyll
GITELSON - Gitelson Chlorophyll Index
GITELSON2 - Gitelson Chlorophyll Index 2
MARI - Modified Anthocyanin Reflectance Index
MCARI1 - Modified Chlorophyll Absorption Ratio Index 1
MCARI2 - Modified Chlorophyll Absorption Ratio Index 2
MCARIOSAVI - MCARI/OSAVI ratio - Chlorophyll content
RDVI - Renormalized Difference Vegetation Index
VREI1 - Vogelmann Red Edge Index 1
VREI2 - Vogelmann Red Edge Index 2

SOIL INDICES (3):
----------------------------------------------------------------------
BI - Brightness Index
CI - Coloration Index
RI - Redness Index

STRESS INDICES (3):
----------------------------------------------------------------------
MSI - Moisture Stress Index
NDNI - Normalized Difference Nitrogen Index
TCARI - Transformed Chlorophyll Absorption Ratio

URBAN INDICES (2):
----------------------------------------------------------------------
NDBI - Normalized Difference Built-up Index
UI - Urban Index

VEGETATION INDICES (14):
----------------------------------------------------------------------
ARVI - Atmospherically Resistant Vegetation Index
CIrededge - Chlorophyll Index Red Edge
DVI - Difference Vegetation Index
EVI - Enhanced Vegetation Index
GNDVI - Green Normalized Difference Vegetation Index
MCARI - Modified Chlorophyll Absorption Ratio Index
MSAVI - Modified Soil Adjusted Vegetation Index
MTCI - MERIS Terrestrial Chlorophyll Index
NDRE - Normalized Difference Red Edge
NDVI - Normalized Difference Vegetation Index
REIP - Red Edge Inflection Point
SAVI - Soil Adjusted Vegetation Index
TVI - Triangular Vegetation Index
VARI - Visible Atmospherically Resistant Index

WATER INDICES (3):
----------------------------------------------------------------------
MNDWI - Modified Normalized Difference Water Index
NDMI - Normalized Difference Moisture Index
NDWI - Normalized Difference Water Index

======================================================================
Total indices available: 77
======================================================================
```
