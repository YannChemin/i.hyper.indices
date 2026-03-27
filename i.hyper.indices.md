## i.hyper.indices

## NAME

i.hyper.indices - Calculates spectral indices from
hyperspectral/multispectral imagery

## KEYWORDS

imagery, indices, hyperspectral, vegetation, classification

## SYNOPSIS

**i.hyper.indices** **-l** **\--i** **-n**
\[**input**=*name\[,name,\...\]*\] \[**wavelengths**=*string*\]
\[**input3d**=*name*\] \[**band_wavelengths**=*string*\]
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

**input**=*name\[,name,\...\]*
:   Input raster bands (comma-separated list). Required unless
    **input3d** is set.

**wavelengths**=*string*
:   Wavelengths for input bands in nm (comma-separated, e.g.,
    450,550,670,800). Required with **input**.

**input3d**=*name*
:   Input 3D raster whose Z-slices are the spectral bands. Each slice
    is extracted in bulk using `Rast3d_extract_z_slice()` (tile-based
    reads with `RASTER3D_NO_CACHE`), which is significantly faster than
    per-voxel access. Mutually exclusive with **input**.

**band_wavelengths**=*string*
:   Wavelengths in nm for each Z-slice of **input3d** (comma-separated,
    bottom slice first). Required with **input3d**.

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

    # Calculate vegetation indices from a 3D hyperspectral raster (bands as Z-slices)
    i.hyper.indices input3d=hyper_cube \
                    band_wavelengths=450,490,550,665,705,740,783,842,865,945 \
                    output_prefix=field \
                    theme=vegetation

    # Calculate all indices from a 3D raster with normalization
    i.hyper.indices -n input3d=hyper_cube \
                    band_wavelengths=450,490,550,665,705,740,783,842 \
                    output_prefix=crop \
                    indices=all
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
-   3D input: slices are extracted via `Rast3d_extract_z_slice()` from
    `libgrass_raster3d`, opened with `RASTER3D_NO_CACHE` so
    `Rast3d_get_block()` uses the tile-bulk path — each tile is read
    exactly once rather than once per voxel
-   Temporary 2D rasters created during 3D extraction are removed
    automatically on exit
-   Either **input**/**wavelengths** or **input3d**/**band_wavelengths**
    must be provided, but not both

## AUTHOR

GRASS Development Team (2026)

## REFERENCES

See individual index references via: i.hyper.indices -li
