## i.hyper.indices

## NAME

i.hyper.indices - Calculates spectral indices from
hyperspectral/multispectral imagery

## KEYWORDS

imagery, indices, hyperspectral, vegetation, classification

## SYNOPSIS

**i.hyper.indices** **-l** **\--i** **-n**
**input**=*name\[,name,\...\]* **wavelengths**=*string*
**output_prefix**=*string* \[**indices**=*string*\]
\[**theme**=*string*\]

## DESCRIPTION

EM Calculates 86+ spectral indices from hyperspectral and multispectral
imagery. Supports wavelength-based automatic band matching with thematic
organization.

## OPTIONS

### Flags:

**-l**
:   List available indices and themes

**-i**
:   Print detailed information about selected indices

**-n**
:   Normalize indices to 0-1 range where applicable

### Parameters:

**input**=*name\[,name,\...\]* *(required)*
:   Input raster bands (comma-separated list)

**wavelengths**=*string* *(required)*
:   Wavelengths for input bands in nm (comma-separated, e.g.,
    450,550,670,800)

**output_prefix**=*string* *(required)*
:   Prefix for output index rasters

**indices**=*string* *(default: ndvi)*
:   Indices to calculate (comma-separated or \'all\' or theme name)

**theme**=*string*
:   Calculate all indices from a specific theme\
    *options:
    vegetation,water,soil,urban,stress,biochemical,pigments,metabolism,materials,all*

## EXAMPLES

::: code
    # List all available indices organized by theme
    i.hyper.indices -l

    # Show detailed information about indices
    i.hyper.indices -li

    # Calculate NDVI from multispectral image
    i.hyper.indices input=band_red,band_nir \
                    wavelengths=665,850 \
                    output_prefix=field \
                    indices=NDVI

    # Calculate all vegetation indices with normalization
    i.hyper.indices -n input=b1,b2,b3,b4,b5 \
                    wavelengths=480,560,665,705,842 \
                    output_prefix=forest \
                    theme=vegetation

    # Detect plastics in coastal imagery
    i.hyper.indices input=coastal_bands \
                    wavelengths=490,560,665,865,1600 \
                    output_prefix=marine \
                    indices=PLASTIC,PDI,MPDI
:::

## THEMES

-   **vegetation**: 15 indices (NDVI, EVI, SAVI, GNDVI, NDRE, MTCI,
    MCARI\...)
-   **pigments**: 16 indices (ARI1, ARI2, CRI550, CRI700, CARI\...)
-   **metabolism**: 13 indices (WI, NDWI1240, LWCI, NDII, SRWI\...)
-   **biochemical**: 6 indices (CAI, NDLI, PRI, SIPI, PSRI)
-   **water**: 3 indices (NDWI, MNDWI, NDMI)
-   **soil**: 3 indices (BI, CI, RI)
-   **urban**: 2 indices (NDBI, UI)
-   **stress**: 4 indices (MSI, NDNI, TCARI)
-   **materials**: 24 indices (HI, THI, PLASTIC, PDI, FPI\...)

## NOTES

-   Automatic color table application (ndvi, water, viridis)
-   Normalization to 0-1 range where applicable (-n flag)
-   Graceful handling of missing spectral bands
-   Dynamic r.mapcalc expression generation

## AUTHOR

GRASS Development Team (2026)

## REFERENCES

See individual index references via: i.hyper.indices -li
