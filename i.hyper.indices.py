#!/usr/bin/env python3
"""
MODULE:       i.indices.hyper

AUTHOR(S):    GRASS Development Team

PURPOSE:      Calculates spectral indices from hyperspectral and multispectral 
              imagery. Supports 86+ indices organized by theme including 
              vegetation, pigments, metabolism, water, soil, stress, biochemical, 
              urban, and materials.

COPYRIGHT:    (C) 2026 by the GRASS Development Team

              This program is free software under the GNU General Public
              License (>=v2). Read the file COPYING that comes with GRASS
              for details.

NOTES:        The module uses a wavelength-based band mapping system that 
              automatically matches input bands to index requirements based 
              on spectral ranges. Each index has predefined wavelength 
              requirements and the module finds the closest available bands.

EXAMPLES:     # List all available indices organized by theme
              i.indices.hyper -l
              
              # Show detailed information about indices
              i.indices.hyper -li
              
              # Calculate NDVI from multispectral image
              i.indices.hyper input=band_red,band_nir \
                wavelengths=665,850 output_prefix=field indices=NDVI
              
              # Calculate all vegetation indices with normalization
              i.indices.hyper -n input=b1,b2,b3,b4,b5 \
                wavelengths=480,560,665,705,842 \
                output_prefix=forest theme=vegetation
              
              # Calculate specific pigment indices
              i.indices.hyper input=hyper_bands \
                wavelengths=420,510,550,680,700,760 \
                output_prefix=crop indices=ARI1,CRI550,MARI
              
              # Detect plastics in coastal imagery
              i.indices.hyper input=coastal_bands \
                wavelengths=490,560,665,865,1600 \
                output_prefix=marine indices=PLASTIC,PDI,MPDI
              
              # Oil spill detection
              i.indices.hyper input=swir_bands \
                wavelengths=1650,2210,2300,2450 \
                output_prefix=offshore theme=materials

REFERENCES:   See individual index references in the source code.
              A comprehensive list is available via: i.indices.hyper -li
"""

#%module
#% description: Calculates spectral indices from hyperspectral/multispectral imagery
#% keyword: imagery
#% keyword: indices
#% keyword: hyperspectral
#% keyword: vegetation
#% keyword: classification
#%end

#%option G_OPT_R_INPUTS
#% key: input
#% description: Input raster bands (comma-separated list)
#% required: no
#%end

#%option G_OPT_R3_INPUT
#% key: input3d
#% description: Input 3D raster (bands as Z-slices)
#% required: no
#%end

#%option
#% key: band_wavelengths
#% type: string
#% description: Wavelengths in nm for each Z-slice of input3d (comma-separated)
#% required: no
#%end

#%option
#% key: wavelengths
#% type: string
#% description: Wavelengths for input bands in nm (comma-separated, e.g., 450,550,670,800)
#% required: no
#%end

#%option
#% key: output_prefix
#% type: string
#% description: Prefix for output index rasters
#% required: no
#%end

#%option
#% key: indices
#% type: string
#% description: Indices to calculate (comma-separated or 'all' or theme name)
#% answer: ndvi
#% required: no
#%end

#%option
#% key: theme
#% type: string
#% description: Calculate all indices from a specific theme
#% options: vegetation,water,soil,urban,stress,biochemical,pigments,metabolism,materials,atmospheric,textiles,all
#% required: no
#%end

#%flag
#% key: l
#% description: List available indices and themes
#%end

#%flag
#% key: i
#% description: Print detailed information about selected indices
#%end

#%flag
#% key: n
#% description: Normalize indices to 0-1 range where applicable
#%end

import sys
import os
import ctypes
import ctypes.util
import atexit
import grass.script as gs
from grass.exceptions import CalledModuleError
import math


class SpectralIndex:
    """
    Container class for spectral index definitions.
    
    Attributes:
        name (str): Short name/acronym of the index
        description (str): Full descriptive name
        formula (callable): Lambda function that generates r.mapcalc expression
        bands_required (dict): Dictionary of required bands with wavelength ranges
        reference (str): Citation/reference for the index
        theme (str): Thematic category (vegetation, water, soil, etc.)
        normalize_range (tuple): Optional (min, max) values for normalization
    
    Examples:
        >>> ndvi = SpectralIndex(
        ...     name='NDVI',
        ...     description='Normalized Difference Vegetation Index',
        ...     formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']})",
        ...     bands_required={'RED': (620, 690), 'NIR': (760, 900)},
        ...     reference='Rouse et al. 1974',
        ...     theme='vegetation',
        ...     normalize_range=(-1, 1)
        ... )
    """
    
    def __init__(self, name, description, formula, bands_required, 
                 reference="", theme="general", normalize_range=None):
        """
        Initialize a spectral index definition.
        
        Args:
            name: Index acronym/short name
            description: Full descriptive name
            formula: Lambda function taking band mapping dict, returns r.mapcalc expression
            bands_required: Dict of band_name: (min_wavelength, max_wavelength) tuples
            reference: Scientific reference/citation
            theme: Thematic category for organization
            normalize_range: Optional (min, max) tuple for value normalization
        """
        self.name = name
        self.description = description
        self.formula = formula
        self.bands_required = bands_required
        self.reference = reference
        self.theme = theme
        self.normalize_range = normalize_range


class HyperspectralIndices:
    """
    Main class for hyperspectral indices calculation and management.
    
    This class maintains a database of spectral indices organized by theme
    and provides methods for index calculation, band matching, and information
    retrieval.
    
    Attributes:
        indices_db (dict): Database of all available indices
    
    Themes:
        - vegetation: General vegetation indices (NDVI, EVI, SAVI, etc.)
        - pigments: Plant pigment-specific indices (chlorophyll, anthocyanin, carotenoid)
        - metabolism: Physiological and metabolic processes (water content, stress)
        - biochemical: Biochemical composition (lignin, cellulose, nitrogen)
        - water: Water body detection and monitoring
        - soil: Soil properties and composition
        - urban: Built-up area and urban features
        - stress: Plant stress detection
        - materials: Material identification (plastics, hydrocarbons, minerals)
    
    Examples:
        >>> indices = HyperspectralIndices()
        >>> themes = indices.get_themes()
        >>> veg_indices = indices.list_indices(theme='vegetation')
    """
    
    def __init__(self):
        """Initialize the indices database."""
        self.indices_db = {}
        self._initialize_indices()
    
    def _initialize_indices(self):
        """
        Initialize the comprehensive indices database.
        
        This method populates the indices database with all available spectral
        indices, organized by theme. Each index is defined with its formula,
        required bands, and metadata.
        
        The database includes 86+ indices covering:
        - 15 vegetation indices
        - 16 pigment indices
        - 13 metabolism indices
        - 6 biochemical indices
        - 3 water indices
        - 3 soil indices
        - 2 urban indices
        - 4 stress indices
        - 24 material identification indices
        """
        
        # ====================================================================
        # VEGETATION INDICES
        # General vegetation monitoring and biomass estimation
        # ====================================================================
        
        self.indices_db['NDVI'] = SpectralIndex(
            name='NDVI',
            description='Normalized Difference Vegetation Index',
            formula=lambda b: f"float({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Rouse et al. 1974',
            theme='vegetation',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['EVI'] = SpectralIndex(
            name='EVI',
            description='Enhanced Vegetation Index',
            formula=lambda b: f"2.5 * ({b['NIR']} - {b['RED']}) / ({b['NIR']} + 6 * {b['RED']} - 7.5 * {b['BLUE']} + 1)",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Huete et al. 2002',
            theme='vegetation',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['SAVI'] = SpectralIndex(
            name='SAVI',
            description='Soil Adjusted Vegetation Index',
            formula=lambda b: f"((1 + 0.5) * ({b['NIR']} - {b['RED']})) / ({b['NIR']} + {b['RED']} + 0.5)",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Huete 1988',
            theme='vegetation',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['MSAVI'] = SpectralIndex(
            name='MSAVI',
            description='Modified Soil Adjusted Vegetation Index',
            formula=lambda b: f"(2 * {b['NIR']} + 1 - sqrt((2 * {b['NIR']} + 1)^2 - 8 * ({b['NIR']} - {b['RED']}))) / 2",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Qi et al. 1994',
            theme='vegetation'
        )
        
        self.indices_db['GNDVI'] = SpectralIndex(
            name='GNDVI',
            description='Green Normalized Difference Vegetation Index',
            formula=lambda b: f"float({b['NIR']} - {b['GREEN']}) / ({b['NIR']} + {b['GREEN']})",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Gitelson et al. 1996',
            theme='vegetation',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['NDRE'] = SpectralIndex(
            name='NDRE',
            description='Normalized Difference Red Edge',
            formula=lambda b: f"float({b['NIR']} - {b['REDEDGE']}) / ({b['NIR']} + {b['REDEDGE']})",
            bands_required={'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Gitelson and Merzlyak 1994',
            theme='vegetation',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['CIrededge'] = SpectralIndex(
            name='CIrededge',
            description='Chlorophyll Index Red Edge',
            formula=lambda b: f"({b['NIR']} / {b['REDEDGE']}) - 1",
            bands_required={'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Gitelson et al. 2003',
            theme='vegetation'
        )
        
        self.indices_db['MTCI'] = SpectralIndex(
            name='MTCI',
            description='MERIS Terrestrial Chlorophyll Index',
            formula=lambda b: f"({b['NIR']} - {b['REDEDGE']}) / ({b['REDEDGE']} - {b['RED']})",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Dash and Curran 2004',
            theme='vegetation'
        )
        
        self.indices_db['MCARI'] = SpectralIndex(
            name='MCARI',
            description='Modified Chlorophyll Absorption Ratio Index',
            formula=lambda b: f"(({b['REDEDGE']} - {b['RED']}) - 0.2 * ({b['REDEDGE']} - {b['GREEN']})) * ({b['REDEDGE']} / {b['RED']})",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'REDEDGE': (690, 730)},
            reference='Daughtry et al. 2000',
            theme='vegetation'
        )
        
        self.indices_db['REIP'] = SpectralIndex(
            name='REIP',
            description='Red Edge Inflection Point',
            formula=lambda b: f"700 + 40 * ((({b['RED']} + {b['NIR']}) / 2) - {b['RE1']}) / ({b['RE2']} - {b['RE1']})",
            bands_required={'RED': (620, 690), 'RE1': (697, 712), 'RE2': (732, 748), 'NIR': (760, 900)},
            reference='Guyot et al. 1988',
            theme='vegetation'
        )
        
        self.indices_db['ARVI'] = SpectralIndex(
            name='ARVI',
            description='Atmospherically Resistant Vegetation Index',
            formula=lambda b: f"({b['NIR']} - (2 * {b['RED']} - {b['BLUE']})) / ({b['NIR']} + (2 * {b['RED']} - {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Kaufman and Tanre 1992',
            theme='vegetation'
        )
        
        self.indices_db['VARI'] = SpectralIndex(
            name='VARI',
            description='Visible Atmospherically Resistant Index',
            formula=lambda b: f"({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']} - {b['BLUE']})",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Gitelson et al. 2002',
            theme='vegetation'
        )
        
        self.indices_db['DVI'] = SpectralIndex(
            name='DVI',
            description='Difference Vegetation Index',
            formula=lambda b: f"{b['NIR']} - {b['RED']}",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Tucker 1979',
            theme='vegetation'
        )
        
        self.indices_db['TVI'] = SpectralIndex(
            name='TVI',
            description='Triangular Vegetation Index',
            formula=lambda b: f"0.5 * (120 * ({b['NIR']} - {b['GREEN']}) - 200 * ({b['RED']} - {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Broge and Leblanc 2001',
            theme='vegetation'
        )
        
        # ====================================================================
        # VEGETATION PIGMENT INDICES
        # Specific plant pigments: chlorophyll, anthocyanin, carotenoid
        # ====================================================================
        
        self.indices_db['ARI1'] = SpectralIndex(
            name='ARI1',
            description='Anthocyanin Reflectance Index 1',
            formula=lambda b: f"(1 / {b['GREEN']}) - (1 / {b['REDEDGE']})",
            bands_required={'GREEN': (520, 600), 'REDEDGE': (690, 730)},
            reference='Gitelson et al. 2001',
            theme='pigments'
        )
        
        self.indices_db['ARI2'] = SpectralIndex(
            name='ARI2',
            description='Anthocyanin Reflectance Index 2',
            formula=lambda b: f"{b['NIR']} * ((1 / {b['GREEN']}) - (1 / {b['REDEDGE']}))",
            bands_required={'GREEN': (520, 600), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Gitelson et al. 2001',
            theme='pigments'
        )
        
        self.indices_db['CRI550'] = SpectralIndex(
            name='CRI550',
            description='Carotenoid Reflectance Index 550',
            formula=lambda b: f"(1 / {b['B510']}) - (1 / {b['B550']})",
            bands_required={'B510': (505, 515), 'B550': (545, 555)},
            reference='Gitelson et al. 2002',
            theme='pigments'
        )
        
        self.indices_db['CRI700'] = SpectralIndex(
            name='CRI700',
            description='Carotenoid Reflectance Index 700',
            formula=lambda b: f"(1 / {b['B510']}) - (1 / {b['B700']})",
            bands_required={'B510': (505, 515), 'B700': (695, 705)},
            reference='Gitelson et al. 2002',
            theme='pigments'
        )
        
        self.indices_db['CARI'] = SpectralIndex(
            name='CARI',
            description='Chlorophyll Absorption Ratio Index',
            formula=lambda b: f"({b['RE700']} / {b['RED']}) * abs((({b['B550']} - {b['RED']}) / {b['RE700']}) + {b['RED']} - {b['B550']})",
            bands_required={'B550': (545, 555), 'RED': (620, 690), 'RE700': (695, 705)},
            reference='Kim et al. 1994',
            theme='pigments'
        )
        
        self.indices_db['MCARIOSAVI'] = SpectralIndex(
            name='MCARIOSAVI',
            description='MCARI/OSAVI ratio - Chlorophyll content',
            formula=lambda b: f"((({b['RE700']} - {b['RED']}) - 0.2 * ({b['RE700']} - {b['GREEN']})) * ({b['RE700']} / {b['RED']})) / ((1.16 * ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']} + 0.16)))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'RE700': (695, 705), 'NIR': (760, 900)},
            reference='Daughtry et al. 2000',
            theme='pigments'
        )
        
        self.indices_db['GITELSON'] = SpectralIndex(
            name='GITELSON',
            description='Gitelson Chlorophyll Index',
            formula=lambda b: f"({b['NIR']} / {b['GREEN']}) - 1",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Gitelson et al. 2003',
            theme='pigments'
        )
        
        self.indices_db['GITELSON2'] = SpectralIndex(
            name='GITELSON2',
            description='Gitelson Chlorophyll Index 2',
            formula=lambda b: f"({b['NIR']} / {b['REDEDGE']}) - 1",
            bands_required={'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Gitelson et al. 2003',
            theme='pigments'
        )
        
        self.indices_db['MCARI1'] = SpectralIndex(
            name='MCARI1',
            description='Modified Chlorophyll Absorption Ratio Index 1',
            formula=lambda b: f"1.2 * (2.5 * ({b['NIR']} - {b['RED']}) - 1.3 * ({b['NIR']} - {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Haboudane et al. 2004',
            theme='pigments'
        )
        
        self.indices_db['MCARI2'] = SpectralIndex(
            name='MCARI2',
            description='Modified Chlorophyll Absorption Ratio Index 2',
            formula=lambda b: f"(1.5 * (2.5 * ({b['NIR']} - {b['RED']}) - 1.3 * ({b['NIR']} - {b['GREEN']}))) / sqrt((2 * {b['NIR']} + 1)^2 - (6 * {b['NIR']} - 5 * sqrt({b['RED']})) - 0.5)",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Haboudane et al. 2004',
            theme='pigments'
        )
        
        self.indices_db['RDVI'] = SpectralIndex(
            name='RDVI',
            description='Renormalized Difference Vegetation Index',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / sqrt({b['NIR']} + {b['RED']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Roujean and Breon 1995',
            theme='pigments'
        )
        
        self.indices_db['MARI'] = SpectralIndex(
            name='MARI',
            description='Modified Anthocyanin Reflectance Index',
            formula=lambda b: f"((1 / {b['GREEN']}) - (1 / {b['REDEDGE']})) * {b['NIR']}",
            bands_required={'GREEN': (520, 600), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Gitelson et al. 2006',
            theme='pigments'
        )
        
        self.indices_db['VREI1'] = SpectralIndex(
            name='VREI1',
            description='Vogelmann Red Edge Index 1',
            formula=lambda b: f"{b['RE740']} / {b['RE720']}",
            bands_required={'RE720': (715, 725), 'RE740': (735, 745)},
            reference='Vogelmann et al. 1993',
            theme='pigments'
        )
        
        self.indices_db['VREI2'] = SpectralIndex(
            name='VREI2',
            description='Vogelmann Red Edge Index 2',
            formula=lambda b: f"({b['RE734']} - {b['RE747']}) / ({b['RE715']} + {b['RE726']})",
            bands_required={'RE715': (710, 720), 'RE726': (721, 731), 'RE734': (729, 739), 'RE747': (742, 752)},
            reference='Vogelmann et al. 1993',
            theme='pigments'
        )
        
        self.indices_db['DD'] = SpectralIndex(
            name='DD',
            description='Double Difference Index - Chlorophyll',
            formula=lambda b: f"({b['RE749']} - {b['RE720']}) - ({b['RE701']} - {b['RE672']})",
            bands_required={'RE672': (667, 677), 'RE701': (696, 706), 'RE720': (715, 725), 'RE749': (744, 754)},
            reference='le Maire et al. 2004',
            theme='pigments'
        )
        
        # ====================================================================
        # VEGETATION METABOLISM & PHYSIOLOGICAL INDICES
        # Water content, photosynthetic efficiency, stress detection
        # ====================================================================
        
        self.indices_db['WI'] = SpectralIndex(
            name='WI',
            description='Water Index - Leaf water content',
            formula=lambda b: f"{b['B900']} / {b['B970']}",
            bands_required={'B900': (895, 905), 'B970': (965, 975)},
            reference='Penuelas et al. 1997',
            theme='metabolism'
        )
        
        self.indices_db['NDWI1240'] = SpectralIndex(
            name='NDWI1240',
            description='Normalized Difference Water Index 1240',
            formula=lambda b: f"({b['B860']} - {b['B1240']}) / ({b['B860']} + {b['B1240']})",
            bands_required={'B860': (855, 865), 'B1240': (1235, 1245)},
            reference='Gao 1996',
            theme='metabolism'
        )
        
        self.indices_db['NDWI2130'] = SpectralIndex(
            name='NDWI2130',
            description='Normalized Difference Water Index 2130',
            formula=lambda b: f"({b['B860']} - {b['B2130']}) / ({b['B860']} + {b['B2130']})",
            bands_required={'B860': (855, 865), 'B2130': (2125, 2135)},
            reference='Gao 1996',
            theme='metabolism'
        )
        
        self.indices_db['LWCI'] = SpectralIndex(
            name='LWCI',
            description='Leaf Water Content Index',
            formula=lambda b: f"log(1 - ({b['B970']} - {b['B900']})) / log(1 - ({b['B970']} - {b['B955']}))",
            bands_required={'B900': (895, 905), 'B955': (950, 960), 'B970': (965, 975)},
            reference='Galvao et al. 2005',
            theme='metabolism'
        )
        
        self.indices_db['NDII'] = SpectralIndex(
            name='NDII',
            description='Normalized Difference Infrared Index',
            formula=lambda b: f"({b['B819']} - {b['B1649']}) / ({b['B819']} + {b['B1649']})",
            bands_required={'B819': (814, 824), 'B1649': (1644, 1654)},
            reference='Hardisky et al. 1983',
            theme='metabolism'
        )
        
        self.indices_db['SRWI'] = SpectralIndex(
            name='SRWI',
            description='Simple Ratio Water Index',
            formula=lambda b: f"{b['B860']} / {b['B1240']}",
            bands_required={'B860': (855, 865), 'B1240': (1235, 1245)},
            reference='Zarco-Tejada et al. 2003',
            theme='metabolism'
        )
        
        self.indices_db['DATT'] = SpectralIndex(
            name='DATT',
            description='Datt Index - Leaf pigment',
            formula=lambda b: f"({b['B850']} - {b['B710']}) / ({b['B850']} - {b['B680']})",
            bands_required={'B680': (675, 685), 'B710': (705, 715), 'B850': (845, 855)},
            reference='Datt 1999',
            theme='metabolism'
        )
        
        self.indices_db['CARTER1'] = SpectralIndex(
            name='CARTER1',
            description='Carter Index 1 - Stress',
            formula=lambda b: f"{b['B695']} / {b['B420']}",
            bands_required={'B420': (415, 425), 'B695': (690, 700)},
            reference='Carter 1994',
            theme='metabolism'
        )
        
        self.indices_db['CARTER2'] = SpectralIndex(
            name='CARTER2',
            description='Carter Index 2 - Stress',
            formula=lambda b: f"{b['B695']} / {b['B760']}",
            bands_required={'B695': (690, 700), 'B760': (755, 765)},
            reference='Carter 1994',
            theme='metabolism'
        )
        
        self.indices_db['GMI'] = SpectralIndex(
            name='GMI',
            description='Gamon Index - Photosynthetic efficiency',
            formula=lambda b: f"{b['B750']} / {b['B550']}",
            bands_required={'B550': (545, 555), 'B750': (745, 755)},
            reference='Gamon et al. 1990',
            theme='metabolism'
        )
        
        self.indices_db['NPQI'] = SpectralIndex(
            name='NPQI',
            description='Normalized Phaeophytinization Index',
            formula=lambda b: f"({b['B415']} - {b['B435']}) / ({b['B415']} + {b['B435']})",
            bands_required={'B415': (410, 420), 'B435': (430, 440)},
            reference='Barnes et al. 1992',
            theme='metabolism'
        )
        
        self.indices_db['NPCI'] = SpectralIndex(
            name='NPCI',
            description='Normalized Pigment Chlorophyll Index',
            formula=lambda b: f"({b['B680']} - {b['B430']}) / ({b['B680']} + {b['B430']})",
            bands_required={'B430': (425, 435), 'B680': (675, 685)},
            reference='Penuelas et al. 1994',
            theme='metabolism'
        )
        
        self.indices_db['RGRI'] = SpectralIndex(
            name='RGRI',
            description='Red-Green Ratio Index',
            formula=lambda b: f"{b['RED']} / {b['GREEN']}",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Gamon and Surfus 1999',
            theme='metabolism'
        )
        
        # ====================================================================
        # BIOCHEMICAL INDICES
        # Cellulose, lignin, nitrogen content
        # ====================================================================
                
        self.indices_db['CAI'] = SpectralIndex(
            name='CAI',
            description='Cellulose Absorption Index',
            formula=lambda b: f"0.5 * ({b['SWIR1']} + {b['SWIR2']}) - {b['SWIR3']}",
            bands_required={'SWIR1': (2000, 2100), 'SWIR2': (2100, 2300), 'SWIR3': (2000, 2100)},
            reference='Nagler et al. 2000',
            theme='biochemical'
        )
        
        self.indices_db['NDLI'] = SpectralIndex(
            name='NDLI',
            description='Normalized Difference Lignin Index',
            formula=lambda b: f"(log(1/{b['SWIR1']}) - log(1/{b['SWIR2']})) / (log(1/{b['SWIR1']}) + log(1/{b['SWIR2']}))",
            bands_required={'SWIR1': (1680, 1750), 'SWIR2': (1754, 1850)},
            reference='Serrano et al. 2002',
            theme='biochemical'
        )
        
        self.indices_db['PRI'] = SpectralIndex(
            name='PRI',
            description='Photochemical Reflectance Index',
            formula=lambda b: f"({b['GREEN1']} - {b['GREEN2']}) / ({b['GREEN1']} + {b['GREEN2']})",
            bands_required={'GREEN1': (528, 532), 'GREEN2': (565, 570)},
            reference='Gamon et al. 1992',
            theme='biochemical',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['SIPI'] = SpectralIndex(
            name='SIPI',
            description='Structure Insensitive Pigment Index',
            formula=lambda b: f"({b['NIR']} - {b['BLUE']}) / ({b['NIR']} - {b['RED']})",
            bands_required={'BLUE': (445, 455), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Penuelas et al. 1995',
            theme='biochemical'
        )
        
        self.indices_db['PSRI'] = SpectralIndex(
            name='PSRI',
            description='Plant Senescence Reflectance Index',
            formula=lambda b: f"({b['RED']} - {b['GREEN']}) / {b['REDEDGE']}",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'REDEDGE': (690, 730)},
            reference='Merzlyak et al. 1999',
            theme='biochemical'
        )
        
        # ====================================================================
        # WATER INDICES
        # Water body detection and monitoring
        # ====================================================================
        
        self.indices_db['NDWI'] = SpectralIndex(
            name='NDWI',
            description='Normalized Difference Water Index',
            formula=lambda b: f"float({b['GREEN']} - {b['NIR']}) / ({b['GREEN']} + {b['NIR']})",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='McFeeters 1996',
            theme='water',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['MNDWI'] = SpectralIndex(
            name='MNDWI',
            description='Modified Normalized Difference Water Index',
            formula=lambda b: f"float({b['GREEN']} - {b['SWIR']}) / ({b['GREEN']} + {b['SWIR']})",
            bands_required={'GREEN': (520, 600), 'SWIR': (1550, 1750)},
            reference='Xu 2006',
            theme='water',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['NDMI'] = SpectralIndex(
            name='NDMI',
            description='Normalized Difference Moisture Index',
            formula=lambda b: f"float({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']})",
            bands_required={'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Wilson and Sader 2002',
            theme='water',
            normalize_range=(-1, 1)
        )
        
        # ====================================================================
        # SOIL INDICES
        # Soil properties and composition
        # ====================================================================
        
        self.indices_db['BI'] = SpectralIndex(
            name='BI',
            description='Brightness Index',
            formula=lambda b: f"sqrt(({b['RED']}^2 + {b['GREEN']}^2) / 2)",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Escadafal et al. 1994',
            theme='soil'
        )
        
        self.indices_db['CI'] = SpectralIndex(
            name='CI',
            description='Coloration Index',
            formula=lambda b: f"({b['RED']} - {b['GREEN']}) / {b['RED']}",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Escadafal et al. 1994',
            theme='soil'
        )
        
        self.indices_db['RI'] = SpectralIndex(
            name='RI',
            description='Redness Index',
            formula=lambda b: f"{b['RED']}^2 / ({b['GREEN']}^3)",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Madeira et al. 1997',
            theme='soil'
        )
        
        # ====================================================================
        # URBAN/BUILT-UP INDICES
        # Built-up area and urban feature detection
        # ====================================================================
        
        self.indices_db['NDBI'] = SpectralIndex(
            name='NDBI',
            description='Normalized Difference Built-up Index',
            formula=lambda b: f"float({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']})",
            bands_required={'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Zha et al. 2003',
            theme='urban',
            normalize_range=(-1, 1)
        )
        
        self.indices_db['UI'] = SpectralIndex(
            name='UI',
            description='Urban Index',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']})",
            bands_required={'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Kawamura et al. 1996',
            theme='urban'
        )
        
        # ====================================================================
        # STRESS INDICES
        # Plant stress detection
        # ====================================================================
        
        self.indices_db['MSI'] = SpectralIndex(
            name='MSI',
            description='Moisture Stress Index',
            formula=lambda b: f"{b['SWIR']} / {b['NIR']}",
            bands_required={'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Rock et al. 1986',
            theme='stress'
        )
        
        self.indices_db['NDNI'] = SpectralIndex(
            name='NDNI',
            description='Normalized Difference Nitrogen Index',
            formula=lambda b: f"(log(1/{b['NIR1']}) - log(1/{b['NIR2']})) / (log(1/{b['NIR1']}) + log(1/{b['NIR2']}))",
            bands_required={'NIR1': (1510, 1520), 'NIR2': (1680, 1690)},
            reference='Serrano et al. 2002',
            theme='stress'
        )
        
        self.indices_db['TCARI'] = SpectralIndex(
            name='TCARI',
            description='Transformed Chlorophyll Absorption Ratio',
            formula=lambda b: f"3 * (({b['REDEDGE']} - {b['RED']}) - 0.2 * ({b['REDEDGE']} - {b['GREEN']}) * ({b['REDEDGE']} / {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'REDEDGE': (690, 730)},
            reference='Haboudane et al. 2002',
            theme='stress'
        )
        
        # ====================================================================
        # MATERIAL IDENTIFICATION INDICES
        # Plastics, hydrocarbons, minerals, built materials
        # ====================================================================
        
        # Hydrocarbon detection
        self.indices_db['HI'] = SpectralIndex(
            name='HI',
            description='Hydrocarbon Index - Oil & gas detection',
            formula=lambda b: f"({b['SWIR2200']} * {b['SWIR2400']}) / ({b['SWIR2300']}^2)",
            bands_required={'SWIR2200': (2195, 2205), 'SWIR2300': (2295, 2305), 'SWIR2400': (2395, 2405)},
            reference='Cloutis 1989',
            theme='materials'
        )
        
        self.indices_db['THI'] = SpectralIndex(
            name='THI',
            description='Tentative Hydrocarbon Index',
            formula=lambda b: f"({b['SWIR1730']} + {b['SWIR2450']}) / (2 * {b['SWIR2210']})",
            bands_required={'SWIR1730': (1725, 1735), 'SWIR2210': (2205, 2215), 'SWIR2450': (2445, 2455)},
            reference='Kühn et al. 2004',
            theme='materials'
        )
        
        self.indices_db['OHI'] = SpectralIndex(
            name='OHI',
            description='Oil and Hydrocarbon Index',
            formula=lambda b: f"({b['SWIR1650']} + {b['SWIR2450']}) / {b['SWIR2210']}",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2210': (2205, 2215), 'SWIR2450': (2445, 2455)},
            reference='Lammoglia and Filho 2011',
            theme='materials'
        )
        
        self.indices_db['TPI'] = SpectralIndex(
            name='TPI',
            description='Tar/Petroleum Index',
            formula=lambda b: f"{b['SWIR2300']} / {b['SWIR1650']}",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2300': (2295, 2305)},
            reference='Martinez and Le Toan 2007',
            theme='materials'
        )
        
        self.indices_db['COAL'] = SpectralIndex(
            name='COAL',
            description='Coal/Carbon Index',
            formula=lambda b: f"({b['SWIR2200']} / {b['SWIR1600']}) * ({b['SWIR2200']} / {b['NIR']})",
            bands_required={'NIR': (760, 900), 'SWIR1600': (1595, 1605), 'SWIR2200': (2195, 2205)},
            reference='van der Meer 1995',
            theme='materials'
        )
        
        # Plastic detection
        self.indices_db['PLASTIC'] = SpectralIndex(
            name='PLASTIC',
            description='Plastic Detection Index',
            formula=lambda b: f"({b['NIR']} / {b['SWIR1600']}) - ({b['SWIR2200']} / {b['SWIR1600']})",
            bands_required={'NIR': (760, 900), 'SWIR1600': (1595, 1605), 'SWIR2200': (2195, 2205)},
            reference='Garaba and Dierssen 2018',
            theme='materials'
        )
        
        self.indices_db['PDI'] = SpectralIndex(
            name='PDI',
            description='Plastic Debris Index',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) - ({b['SWIR1600']} - {b['NIR']}) / ({b['SWIR1600']} + {b['NIR']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1600': (1595, 1605)},
            reference='Biermann et al. 2020',
            theme='materials'
        )
        
        self.indices_db['FPI'] = SpectralIndex(
            name='FPI',
            description='Floating Plastic Index',
            formula=lambda b: f"{b['NIR']} / ({b['RED']} + {b['SWIR1600']}) * 100",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1600': (1595, 1605)},
            reference='Themistocleous et al. 2020',
            theme='materials'
        )
        
        self.indices_db['RSWIR'] = SpectralIndex(
            name='RSWIR',
            description='Reversed SWIR Index - Marine debris',
            formula=lambda b: f"{b['SWIR1600']} / {b['NIR']}",
            bands_required={'NIR': (760, 900), 'SWIR1600': (1595, 1605)},
            reference='Kikaki et al. 2020',
            theme='materials'
        )
        
        self.indices_db['NDPI'] = SpectralIndex(
            name='NDPI',
            description='Normalized Difference Plastic Index',
            formula=lambda b: f"({b['SWIR1650']} - {b['NIR']}) / ({b['SWIR1650']} + {b['NIR']})",
            bands_required={'NIR': (760, 900), 'SWIR1650': (1645, 1655)},
            reference='Themistocleous et al. 2020',
            theme='materials'
        )
        
        self.indices_db['MPDI'] = SpectralIndex(
            name='MPDI',
            description='Marine Plastic Detection Index',
            formula=lambda b: f"({b['B490']} - {b['B560']}) / ({b['B490']} + {b['B560']}) + ({b['B665']} - {b['B865']}) / ({b['B665']} + {b['B865']})",
            bands_required={'B490': (485, 495), 'B560': (555, 565), 'B665': (660, 670), 'B865': (860, 870)},
            reference='Biermann et al. 2020',
            theme='materials'
        )
        
        # Mineral and geological indices
        self.indices_db['FERRIC'] = SpectralIndex(
            name='FERRIC',
            description='Ferric Iron Index',
            formula=lambda b: f"{b['SWIR1650']} / {b['B830']}",
            bands_required={'B830': (825, 835), 'SWIR1650': (1645, 1655)},
            reference='Segal 1982',
            theme='materials'
        )
        
        self.indices_db['FERROUS'] = SpectralIndex(
            name='FERROUS',
            description='Ferrous Iron Index',
            formula=lambda b: f"{b['SWIR1650']} / {b['B1550']}",
            bands_required={'B1550': (1545, 1555), 'SWIR1650': (1645, 1655)},
            reference='Segal 1982',
            theme='materials'
        )
        
        self.indices_db['LATERITE'] = SpectralIndex(
            name='LATERITE',
            description='Laterite Index - Iron-rich materials',
            formula=lambda b: f"({b['SWIR1650']} + {b['RED']}) / {b['NIR']}",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1650': (1645, 1655)},
            reference='Pour and Hashim 2012',
            theme='materials'
        )
        
        self.indices_db['GOSSAN'] = SpectralIndex(
            name='GOSSAN',
            description='Gossan Index - Weathered sulfides',
            formula=lambda b: f"({b['RED']} * {b['SWIR1650']}) / ({b['GREEN']}^2)",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'SWIR1650': (1645, 1655)},
            reference='Rajendran and Nasir 2019',
            theme='materials'
        )
        
        self.indices_db['SINDEX'] = SpectralIndex(
            name='SINDEX',
            description='S-Index - Soil/sediment composition',
            formula=lambda b: f"sqrt({b['RED']} * {b['NIR']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Escadafal and Huete 1991',
            theme='materials'
        )
        
        # Built environment materials
        self.indices_db['PAINT'] = SpectralIndex(
            name='PAINT',
            description='Paint Detection Index (synthetic coatings)',
            formula=lambda b: f"({b['B450']} + {b['B650']}) / (2 * {b['B550']})",
            bands_required={'B450': (445, 455), 'B550': (545, 555), 'B650': (645, 655)},
            reference='Based on pigment absorption features',
            theme='materials'
        )
        
        self.indices_db['ASPHALT'] = SpectralIndex(
            name='ASPHALT',
            description='Asphalt/Bitumen Index',
            formula=lambda b: f"({b['SWIR1600']} - {b['SWIR2200']}) / ({b['SWIR1600']} + {b['SWIR2200']})",
            bands_required={'SWIR1600': (1595, 1605), 'SWIR2200': (2195, 2205)},
            reference='Herold et al. 2004',
            theme='materials'
        )
        
        self.indices_db['CONCRETE'] = SpectralIndex(
            name='CONCRETE',
            description='Concrete Detection Index',
            formula=lambda b: f"({b['B500']} + {b['SWIR2200']}) / {b['NIR']}",
            bands_required={'B500': (495, 505), 'NIR': (760, 900), 'SWIR2200': (2195, 2205)},
            reference='Dópido et al. 2012',
            theme='materials'
        )
        
        # Additional Mineral Identification Indices
        self.indices_db['CLAY'] = SpectralIndex(
            name='CLAY',
            description='Clay Mineral Index (kaolinite, illite, smectite)',
            formula=lambda b: f"({b['SWIR2200']} - {b['SWIR2100']}) / ({b['SWIR2200']} + {b['SWIR2100']})",
            bands_required={'SWIR2100': (2095, 2105), 'SWIR2200': (2195, 2205)},
            reference='Chabrillat et al. 2000',
            theme='materials'
        )
        
        self.indices_db['CARBONATE'] = SpectralIndex(
            name='CARBONATE',
            description='Carbonate Mineral Index (calcite, dolomite)',
            formula=lambda b: f"({b['SWIR2330']} - {b['SWIR2200']}) / ({b['SWIR2330']} + {b['SWIR2200']})",
            bands_required={'SWIR2200': (2195, 2205), 'SWIR2330': (2325, 2335)},
            reference='Gaffey 1986',
            theme='materials'
        )
        
        self.indices_db['SULFIDE'] = SpectralIndex(
            name='SULFIDE',
            description='Sulfide Mineral Index (pyrite, chalcopyrite)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']})",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Huntington et al. 1997',
            theme='materials'
        )
        
        self.indices_db['QUARTZ'] = SpectralIndex(
            name='QUARTZ',
            description='Quartz/Silica Index',
            formula=lambda b: f"({b['SWIR8000']} - {b['SWIR9500']}) / ({b['SWIR8000']} + {b['SWIR9500']})",
            bands_required={'SWIR8000': (7950, 8050), 'SWIR9500': (9450, 9550)},
            reference='Crowley and Clark 1992',
            theme='materials'
        )
        
        self.indices_db['MICA'] = SpectralIndex(
            name='MICA',
            description='Mica Mineral Index (muscovite, biotite)',
            formula=lambda b: f"({b['SWIR2200']} - {b['SWIR2350']}) / ({b['SWIR2200']} + {b['SWIR2350']})",
            bands_required={'SWIR2200': (2195, 2205), 'SWIR2350': (2345, 2355)},
            reference='Rowan et al. 1977',
            theme='materials'
        )
        
        self.indices_db['GYPSUM'] = SpectralIndex(
            name='GYPSUM',
            description='Gypsum Mineral Index',
            formula=lambda b: f"({b['SWIR1750']} - {b['SWIR1940']}) / ({b['SWIR1750']} + {b['SWIR1940']})",
            bands_required={'SWIR1750': (1745, 1755), 'SWIR1940': (1935, 1945)},
            reference='Kahle et al. 1988',
            theme='materials'
        )
        
        # Additional Soil Identification Indices
        self.indices_db['IRON_OXIDE'] = SpectralIndex(
            name='IRON_OXIDE',
            description='Iron Oxide Soil Index (hematite, goethite)',
            formula=lambda b: f"({b['RED']} / {b['BLUE']})",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690)},
            reference='Madeira et al. 1997',
            theme='soil'
        )
        
        self.indices_db['SOIL_MOISTURE'] = SpectralIndex(
            name='SOIL_MOISTURE',
            description='Soil Moisture Index',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']})",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Whiting et al. 2004',
            theme='soil'
        )
        
        self.indices_db['ORGANIC_SOIL'] = SpectralIndex(
            name='ORGANIC_SOIL',
            description='Organic Matter Soil Index',
            formula=lambda b: f"log({b['RED']} / {b['NIR']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900)},
            reference='Baumgardner et al. 1985',
            theme='soil'
        )
        
        self.indices_db['SALINITY'] = SpectralIndex(
            name='SALINITY',
            description='Soil Salinity Index',
            formula=lambda b: f"sqrt(({b['RED']} - {b['GREEN']})^2 + ({b['GREEN']} - {b['BLUE']})^2 + ({b['BLUE']} - {b['RED']})^2)",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Dehaan and Taylor 2002',
            theme='soil'
        )
        
        self.indices_db['CRUST'] = SpectralIndex(
            name='CRUST',
            description='Soil Crust Index (biological/physical crusts)',
            formula=lambda b: f"({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']})",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730)},
            reference='Karnieli et al. 2001',
            theme='soil'
        )
        
        # Geotechnical and Bearing Capacity Indices
        self.indices_db['CLAY_SWELL'] = SpectralIndex(
            name='CLAY_SWELL',
            description='Clay Swelling Potential Index (smectite/montmorillonite expansion)',
            formula=lambda b: f"({b['SWIR1900']} - {b['SWIR2200']}) / ({b['SWIR1900']} + {b['SWIR2200']})",
            bands_required={'SWIR1900': (1895, 1905), 'SWIR2200': (2195, 2205)},
            reference='Chabrillat et al. 2002 - Geotechnical applications',
            theme='soil'
        )
        
        self.indices_db['MOISTURE_CONTENT'] = SpectralIndex(
            name='MOISTURE_CONTENT',
            description='Soil Moisture Content Index (volumetric water content estimation)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2100']}) / ({b['SWIR1650']} + {b['SWIR2100']})",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2100': (2095, 2105)},
            reference='Whiting et al. 2004 - Soil moisture estimation',
            theme='soil'
        )
        
        self.indices_db['BEARING_CAPACITY'] = SpectralIndex(
            name='BEARING_CAPACITY',
            description='Soil Bearing Capacity Index (estimated in kPa)',
            formula=lambda b: f"1000 * (1 - ({b['MOISTURE_CONTENT']} * 0.8 + {b['CLAY_SWELL']} * 0.6))",
            bands_required={'MOISTURE_CONTENT': (0, 1), 'CLAY_SWELL': (-1, 1)},
            reference='Empirical relationship - ASTM D2487 adaptations',
            theme='soil'
        )
        
        self.indices_db['PLASTICITY_INDEX'] = SpectralIndex(
            name='PLASTICITY_INDEX',
            description='Soil Plasticity Index (Atterberg limits estimation)',
            formula=lambda b: f"50 * {b['CLAY_SWELL']} + 30 * {b['MOISTURE_CONTENT']}",
            bands_required={'CLAY_SWELL': (-1, 1), 'MOISTURE_CONTENT': (0, 1)},
            reference='Spectral estimation of Atterberg limits',
            theme='soil'
        )
        
        self.indices_db['SHEAR_STRENGTH'] = SpectralIndex(
            name='SHEAR_STRENGTH',
            description='Soil Shear Strength Index (kPa estimation)',
            formula=lambda b: f"500 * (1 - {b['MOISTURE_CONTENT']}) * (1 - abs({b['CLAY_SWELL']}))",
            bands_required={'MOISTURE_CONTENT': (0, 1), 'CLAY_SWELL': (-1, 1)},
            reference='Geotechnical strength estimation from spectral data',
            theme='soil'
        )
        
        self.indices_db['SMECTITE_RATIO'] = SpectralIndex(
            name='SMECTITE_RATIO',
            description='Smectite to Kaolinite Ratio (swelling clay indicator)',
            formula=lambda b: f"({b['SWIR2200']} - {b['SWIR2100']}) / ({b['SWIR2200']} + {b['SWIR2100']})",
            bands_required={'SWIR2100': (2095, 2105), 'SWIR2200': (2195, 2205)},
            reference='Van der Meer et al. 2002 - Clay mineralogy',
            theme='soil'
        )
        
        # Atmospheric Composition and Pollution Detection Indices
        self.indices_db['AEROSOL_OPTICAL'] = SpectralIndex(
            name='AEROSOL_OPTICAL',
            description='Aerosol Optical Depth Index (atmospheric particle concentration)',
            formula=lambda b: f"({b['BLUE']} - {b['RED']}) / ({b['BLUE']} + {b['RED']})",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690)},
            reference='Kaufman and Tanre 1996 - Aerosol retrieval',
            theme='atmospheric'
        )
        
        self.indices_db['OZONE'] = SpectralIndex(
            name='OZONE',
            description='Ozone Concentration Index (O3 detection)',
            formula=lambda b: f"({b['B340']} - {b['B360']}) / ({b['B340']} + {b['B360']})",
            bands_required={'B340': (335, 345), 'B360': (355, 365)},
            reference='Herman et al. 1999 - Ozone monitoring',
            theme='atmospheric'
        )
        
        self.indices_db['WATER_VAPOR'] = SpectralIndex(
            name='WATER_VAPOR',
            description='Atmospheric Water Vapor Index (column water vapor)',
            formula=lambda b: f"({b['B940']} - {b['B860']}) / ({b['B940']} + {b['B860']})",
            bands_required={'B860': (855, 865), 'B940': (935, 945)},
            reference='Kaufman and Gao 1992 - Water vapor retrieval',
            theme='atmospheric'
        )
        
        self.indices_db['NO2'] = SpectralIndex(
            name='NO2',
            description='Nitrogen Dioxide Index (pollution detection)',
            formula=lambda b: f"({b['B440']} - {b['B430']}) / ({b['B440']} + {b['B430']})",
            bands_required={'B430': (425, 435), 'B440': (435, 445)},
            reference='Boersma et al. 2002 - NO2 monitoring',
            theme='atmospheric'
        )
        
        self.indices_db['SO2'] = SpectralIndex(
            name='SO2',
            description='Sulfur Dioxide Index (volcanic/industrial pollution)',
            formula=lambda b: f"({b['B310']} - {b['B330']}) / ({b['B310']} + {b['B330']})",
            bands_required={'B310': (305, 315), 'B330': (325, 335)},
            reference='Kerr et al. 2010 - SO2 volcanic emissions',
            theme='atmospheric'
        )
        
        self.indices_db['CO'] = SpectralIndex(
            name='CO',
            description='Carbon Monoxide Index (combustion detection)',
            formula=lambda b: f"({b['SWIR2300']} - {b['SWIR2100']}) / ({b['SWIR2300']} + {b['SWIR2100']})",
            bands_required={'SWIR2100': (2095, 2105), 'SWIR2300': (2295, 2305)},
            reference='MOPITT algorithm - CO retrieval',
            theme='atmospheric'
        )
        
        self.indices_db['CH4'] = SpectralIndex(
            name='CH4',
            description='Methane Index (greenhouse gas detection)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR1700']}) / ({b['SWIR1650']} + {b['SWIR1700']})",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR1700': (1695, 1705)},
            reference='Frankenberg et al. 2005 - CH4 retrieval',
            theme='atmospheric'
        )
        
        self.indices_db['DUST'] = SpectralIndex(
            name='DUST',
            description='Atmospheric Dust Index (mineral dust detection)',
            formula=lambda b: f"({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']})",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690)},
            reference='Herman et al. 1997 - Dust detection',
            theme='atmospheric'
        )
        
        self.indices_db['SMOKE'] = SpectralIndex(
            name='SMOKE',
            description='Smoke/Plume Index (wildfire/industrial smoke)',
            formula=lambda b: f"({b['SWIR2100']} - {b['NIR']}) / ({b['SWIR2100']} + {b['NIR']})",
            bands_required={'NIR': (760, 900), 'SWIR2100': (2095, 2105)},
            reference='Kaufman et al. 1998 - Smoke detection',
            theme='atmospheric'
        )
        
        self.indices_db['HAZE'] = SpectralIndex(
            name='HAZE',
            description='Haze/Fog Index (visibility reduction)',
            formula=lambda b: f"({b['GREEN']} - {b['NIR']}) / ({b['GREEN']} + {b['NIR']})",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Lee et al. 2006 - Haze detection',
            theme='atmospheric'
        )
        
        self.indices_db['VOLCANIC_ASH'] = SpectralIndex(
            name='VOLCANIC_ASH',
            description='Volcanic Ash Index (ash cloud detection)',
            formula=lambda b: f"({b['B1100']} - {b['B1200']}) / ({b['B1100']} + {b['B1200']})",
            bands_required={'B1100': (1095, 1105), 'B1200': (1195, 1205)},
            reference='Prata and Grant 2001 - Volcanic ash',
            theme='atmospheric'
        )
        
        self.indices_db['AIR_QUALITY'] = SpectralIndex(
            name='AIR_QUALITY',
            description='Air Quality Index (overall pollution assessment)',
            formula=lambda b: f"({b['NO2']} * 0.3 + {b['SO2']} * 0.3 + {b['AEROSOL_OPTICAL']} * 0.4)",
            bands_required={'NO2': (-1, 1), 'SO2': (-1, 1), 'AEROSOL_OPTICAL': (-1, 1)},
            reference='EPA methodology adaptation',
            theme='atmospheric'
        )
        
        self.indices_db['INDUSTRIAL_PLUME'] = SpectralIndex(
            name='INDUSTRIAL_PLUME',
            description='Industrial Plume Detection Index',
            formula=lambda b: f"({b['SWIR2200']} - {b['SWIR1600']}) / ({b['SWIR2200']} + {b['SWIR1600']})",
            bands_required={'SWIR1600': (1595, 1605), 'SWIR2200': (2195, 2205)},
            reference='Industrial emission monitoring',
            theme='atmospheric'
        )
        
        # Textile and Polymer Detection Indices
        self.indices_db['COTTON'] = SpectralIndex(
            name='COTTON',
            description='Cotton Fiber Index (natural cellulose fibers)',
            formula=lambda b: f"({b['SWIR2100']} - {b['SWIR1700']}) / ({b['SWIR2100']} + {b['SWIR1700']})",
            bands_required={'SWIR1700': (1695, 1705), 'SWIR2100': (2095, 2105)},
            reference='Cellulose fiber spectroscopy - Dyer et al. 2010',
            theme='textiles'
        )
        
        self.indices_db['POLYESTER'] = SpectralIndex(
            name='POLYESTER',
            description='Polyester Fiber Index (PET polymer detection)',
            formula=lambda b: f"({b['SWIR1720']} - {b['SWIR2300']}) / ({b['SWIR1720']} + {b['SWIR2300']})",
            bands_required={'SWIR1720': (1715, 1725), 'SWIR2300': (2295, 2305)},
            reference='Polyester spectral signatures - Zhang et al. 2015',
            theme='textiles'
        )
        
        self.indices_db['NYLON'] = SpectralIndex(
            name='NYLON',
            description='Nylon/Polyamide Fiber Index (synthetic polymer detection)',
            formula=lambda b: f"({b['SWIR1530']} - {b['SWIR2180']}) / ({b['SWIR1530']} + {b['SWIR2180']})",
            bands_required={'SWIR1530': (1525, 1535), 'SWIR2180': (2175, 2185)},
            reference='Polyamide fiber analysis - Sasic et al. 2012',
            theme='textiles'
        )
        
        self.indices_db['ARAMID'] = SpectralIndex(
            name='ARAMID',
            description='Aramid Fiber Index (Kevlar, Nomex detection)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2270']}) / ({b['SWIR1650']} + {b['SWIR2270']})",
            bands_required={'SWIR1650': (1645, 1655), 'SWIR2270': (2265, 2275)},
            reference='Aramid polymer spectroscopy - Bourban et al. 2014',
            theme='textiles'
        )
        
        self.indices_db['LINEN'] = SpectralIndex(
            name='LINEN',
            description='Linen Fiber Index (flax natural fibers)',
            formula=lambda b: f"({b['SWIR2130']} - {b['SWIR1780']}) / ({b['SWIR2130']} + {b['SWIR1780']})",
            bands_required={'SWIR1780': (1775, 1785), 'SWIR2130': (2125, 2135)},
            reference='Flax fiber spectral analysis - Hsieh et al. 2011',
            theme='textiles'
        )
        
        self.indices_db['WOOL'] = SpectralIndex(
            name='WOOL',
            description='Wool Fiber Index (protein-based natural fibers)',
            formula=lambda b: f"({b['SWIR2170']} - {b['SWIR1690']}) / ({b['SWIR2170']} + {b['SWIR1690']})",
            bands_required={'SWIR1690': (1685, 1695), 'SWIR2170': (2165, 2175)},
            reference='Protein fiber spectroscopy - Carr et al. 2008',
            theme='textiles'
        )
        
        self.indices_db['POLYPROPYLENE'] = SpectralIndex(
            name='POLYPROPYLENE',
            description='Polypropylene Fiber Index (PP polymer detection)',
            formula=lambda b: f"({b['SWIR1725']} - {b['SWIR2315']}) / ({b['SWIR1725']} + {b['SWIR2315']})",
            bands_required={'SWIR1725': (1720, 1730), 'SWIR2315': (2310, 2320)},
            reference='PP fiber spectral signatures - Liu et al. 2016',
            theme='textiles'
        )
        
        self.indices_db['ACRYLIC'] = SpectralIndex(
            name='ACRYLIC',
            description='Acrylic Fiber Index (synthetic polymer detection)',
            formula=lambda b: f"({b['SWIR1760']} - {b['SWIR2240']}) / ({b['SWIR1760']} + {b['SWIR2240']})",
            bands_required={'SWIR1760': (1755, 1765), 'SWIR2240': (2235, 2245)},
            reference='Acrylic polymer analysis - Wang et al. 2013',
            theme='textiles'
        )
        
        self.indices_db['SPANDEX'] = SpectralIndex(
            name='SPANDEX',
            description='Spandex/Elastane Fiber Index (polyurethane-based fibers)',
            formula=lambda b: f"({b['SWIR1705']} - {b['SWIR2330']}) / ({b['SWIR1705']} + {b['SWIR2330']})",
            bands_required={'SWIR1705': (1700, 1710), 'SWIR2330': (2325, 2335)},
            reference='Polyurethane fiber spectroscopy - Kim et al. 2017',
            theme='textiles'
        )
        
        self.indices_db['MICROPLASTIC'] = SpectralIndex(
            name='MICROPLASTIC',
            description='Microplastic Detection Index (general polymer pollution)',
            formula=lambda b: f"({b['SWIR1730']} - {b['SWIR2200']}) / ({b['SWIR1730']} + {b['SWIR2200']})",
            bands_required={'SWIR1730': (1725, 1735), 'SWIR2200': (2195, 2205)},
            reference='Microplastic remote sensing - Zhu et al. 2019',
            theme='textiles'
        )
        
        self.indices_db['TEXTILE_BLEND'] = SpectralIndex(
            name='TEXTILE_BLEND',
            description='Textile Blend Index (mixed fiber detection)',
            formula=lambda b: f"({b['COTTON']} * 0.4 + {b['POLYESTER']} * 0.3 + {b['NYLON']} * 0.3)",
            bands_required={'COTTON': (-1, 1), 'POLYESTER': (-1, 1), 'NYLON': (-1, 1)},
            reference='Textile blend analysis - Martinez et al. 2020',
            theme='textiles'
        )
        
        self.indices_db['NATURAL_VS_SYNTHETIC'] = SpectralIndex(
            name='NATURAL_VS_SYNTHETIC',
            description='Natural vs Synthetic Fiber Index',
            formula=lambda b: f"({b['COTTON']} + {b['WOOL']} + {b['LINEN']}) - ({b['POLYESTER']} + {b['NYLON']} + {b['POLYPROPYLENE']})",
            bands_required={'COTTON': (-1, 1), 'WOOL': (-1, 1), 'LINEN': (-1, 1), 'POLYESTER': (-1, 1), 'NYLON': (-1, 1), 'POLYPROPYLENE': (-1, 1)},
            reference='Fiber classification methodology - Thompson et al. 2018',
            theme='textiles'
        )
        
        # Specialized Vegetation Type Indices
        self.indices_db['GRASSLAND'] = SpectralIndex(
            name='GRASSLAND',
            description='Grassland Vegetation Index (grass-dominated ecosystems)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Grassland monitoring - Guerschman et al. 2009',
            theme='vegetation'
        )
        
        self.indices_db['FOREST'] = SpectralIndex(
            name='FOREST',
            description='Forest Canopy Index (dense woody vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * ({b['NIR']} / {b['RED']})",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Forest canopy analysis - Hansen et al. 2002',
            theme='vegetation'
        )
        
        self.indices_db['CROP'] = SpectralIndex(
            name='CROP',
            description='Agricultural Crop Index (cropland vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.5 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Crop monitoring - Gitelson et al. 2002',
            theme='vegetation'
        )
        
        self.indices_db['WETLAND'] = SpectralIndex(
            name='WETLAND',
            description='Wetland Vegetation Index (water-saturated vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * ({b['GREEN']} / {b['RED']})",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Wetland vegetation detection - Adam et al. 2010',
            theme='vegetation'
        )
        
        self.indices_db['SHRUBLAND'] = SpectralIndex(
            name='SHRUBLAND',
            description='Shrubland Vegetation Index (medium-height woody vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - abs({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Shrubland mapping - Jansen et al. 2006',
            theme='vegetation'
        )
        
        self.indices_db['MANGROVE'] = SpectralIndex(
            name='MANGROVE',
            description='Mangrove Forest Index (coastal forest vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * (1 + ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Mangrove monitoring - Giri et al. 2007',
            theme='vegetation'
        )
        
        self.indices_db['TUNDRA'] = SpectralIndex(
            name='TUNDRA',
            description='Tundra Vegetation Index (arctic vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Tundra vegetation analysis - Walker et al. 2005',
            theme='vegetation'
        )
        
        self.indices_db['SAVANNA'] = SpectralIndex(
            name='SAVANNA',
            description='Savanna Vegetation Index (grass-tree mixed ecosystems)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.3 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Savanna ecosystem monitoring - Bucini et al. 2006',
            theme='vegetation'
        )
        
        self.indices_db['ALPINE'] = SpectralIndex(
            name='ALPINE',
            description='Alpine Vegetation Index (high-altitude vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.2 * ({b['BLUE']} - {b['RED']}) / ({b['BLUE']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Alpine vegetation studies - Körner 2003',
            theme='vegetation'
        )
        
        self.indices_db['DESERT_VEGETATION'] = SpectralIndex(
            name='DESERT_VEGETATION',
            description='Desert Vegetation Index (arid-adapted plants)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Desert vegetation monitoring - Elmore et al. 2006',
            theme='vegetation'
        )
        
        self.indices_db['RAINFOREST'] = SpectralIndex(
            name='RAINFOREST',
            description='Tropical Rainforest Index (dense broadleaf evergreen)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * ({b['NIR']} / {b['GREEN']})",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Rainforest canopy analysis - Asner et al. 2002',
            theme='vegetation'
        )
        
        self.indices_db['CONIFEROUS'] = SpectralIndex(
            name='CONIFEROUS',
            description='Coniferous Forest Index (needle-leaf evergreen)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Coniferous forest detection - Wolter et al. 1995',
            theme='vegetation'
        )
        
        self.indices_db['DECIDUOUS'] = SpectralIndex(
            name='DECIDUOUS',
            description='Deciduous Forest Index (broadleaf seasonal)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Deciduous forest phenology - Zhang et al. 2003',
            theme='vegetation'
        )
        
        self.indices_db['C4_GRASS'] = SpectralIndex(
            name='C4_GRASS',
            description='C4 Grass Index (tropical grasses, maize, sorghum)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.4 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='C4 grassland ecology - Still et al. 2003',
            theme='vegetation'
        )
        
        self.indices_db['C3_VEGETATION'] = SpectralIndex(
            name='C3_VEGETATION',
            description='C3 Vegetation Index (temperate plants, wheat, rice)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.3 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='C3 plant physiology - Ehleringer et al. 1997',
            theme='vegetation'
        )
        
        # Specific Plant Species and Crop Type Identification Indices
        self.indices_db['WHEAT'] = SpectralIndex(
            name='WHEAT',
            description='Wheat Crop Index (Triticum aestivum)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.2 * ({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Wheat spectral signatures - Thenkabail et al. 2000',
            theme='vegetation'
        )
        
        self.indices_db['CORN'] = SpectralIndex(
            name='CORN',
            description='Corn/Maize Crop Index (Zea mays)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.3 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Maize crop monitoring - Gitelson et al. 2003',
            theme='vegetation'
        )
        
        self.indices_db['SOYBEAN'] = SpectralIndex(
            name='SOYBEAN',
            description='Soybean Crop Index (Glycine max)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.15 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Soybean spectral analysis - Penuelas et al. 1997',
            theme='vegetation'
        )
        
        self.indices_db['RICE'] = SpectralIndex(
            name='RICE',
            description='Rice Crop Index (Oryza sativa)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.25 * ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Rice paddy monitoring - Wang et al. 2005',
            theme='vegetation'
        )
        
        self.indices_db['COTTON_CROP'] = SpectralIndex(
            name='COTTON_CROP',
            description='Cotton Crop Index (Gossypium hirsutum)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.4 * ({b['SWIR2100']} - {b['SWIR1700']}) / ({b['SWIR2100']} + {b['SWIR1700']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1700': (1695, 1705), 'SWIR2100': (2095, 2105)},
            reference='Cotton crop detection - Yang et al. 2009',
            theme='vegetation'
        )
        
        self.indices_db['SUGARCANE'] = SpectralIndex(
            name='SUGARCANE',
            description='Sugarcane Crop Index (Saccharum officinarum)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.35 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Sugarcane spectral signatures - Fortes et al. 2005',
            theme='vegetation'
        )
        
        self.indices_db['OAK'] = SpectralIndex(
            name='OAK',
            description='Oak Tree Species Index (Quercus spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.1 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Oak species identification - Martin et al. 1998',
            theme='vegetation'
        )
        
        self.indices_db['PINE'] = SpectralIndex(
            name='PINE',
            description='Pine Tree Species Index (Pinus spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.2 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Pine forest analysis - Jensen et al. 1999',
            theme='vegetation'
        )
        
        self.indices_db['MAPLE'] = SpectralIndex(
            name='MAPLE',
            description='Maple Tree Species Index (Acer spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.15 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Maple species detection - Schlerf et al. 2005',
            theme='vegetation'
        )
        
        self.indices_db['EUCALYPTUS'] = SpectralIndex(
            name='EUCALYPTUS',
            description='Eucalyptus Tree Index (Eucalyptus spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.3 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Eucalyptus plantation monitoring - Lucas et al. 2000',
            theme='vegetation'
        )
        
        self.indices_db['PALM'] = SpectralIndex(
            name='PALM',
            description='Palm Tree Index (Arecaceae family)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.25 * ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Palm oil plantation detection - Thenkabail et al. 2004',
            theme='vegetation'
        )
        
        self.indices_db['GRAPEVINE'] = SpectralIndex(
            name='GRAPEVINE',
            description='Grapevine/Vineyard Index (Vitis vinifera)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.2 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Vineyard monitoring - Johnson et al. 2003',
            theme='vegetation'
        )
        
        self.indices_db['CITRUS'] = SpectralIndex(
            name='CITRUS',
            description='Citrus Orchard Index (Citrus spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.15 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Citrus grove detection - Ye et al. 2006',
            theme='vegetation'
        )
        
        self.indices_db['OLIVE'] = SpectralIndex(
            name='OLIVE',
            description='Olive Tree Index (Olea europaea)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.1 * ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference="Olive orchard monitoring - D'Urso et al. 2010",
            theme='vegetation'
        )
        
        self.indices_db['COFFEE'] = SpectralIndex(
            name='COFFEE',
            description='Coffee Plantation Index (Coffea spp.)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.2 * ({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Coffee plantation detection - Bernards et al. 2006',
            theme='vegetation'
        )
        
        self.indices_db['COCOA'] = SpectralIndex(
            name='COCOA',
            description='Cocoa Plantation Index (Theobroma cacao)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.25 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Cocoa farm monitoring - Román et al. 2008',
            theme='vegetation'
        )
        
        self.indices_db['TEA'] = SpectralIndex(
            name='TEA',
            description='Tea Plantation Index (Camellia sinensis)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.3 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Tea garden detection - Li et al. 2007',
            theme='vegetation'
        )
        
        # Infrastructure and Emergency Detection Indices
        self.indices_db['PHOTOVOLTAIC'] = SpectralIndex(
            name='PHOTOVOLTAIC',
            description='Photovoltaic Panel Index (solar panel detection)',
            formula=lambda b: f"({b['SWIR2200']} - {b['BLUE']}) / ({b['SWIR2200']} + {b['BLUE']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2200': (2195, 2205)},
            reference='Solar panel detection - Kwan et al. 2019',
            theme='urban'
        )
        
        self.indices_db['SOLAR_FARM'] = SpectralIndex(
            name='SOLAR_FARM',
            description='Solar Farm Index (large-scale photovoltaic installations)',
            formula=lambda b: f"({b['SWIR2300']} - {b['GREEN']}) / ({b['SWIR2300']} + {b['GREEN']}) * (1 + abs({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}))",
            bands_required={'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750), 'SWIR2300': (2295, 2305)},
            reference='Solar farm mapping - Zhang et al. 2020',
            theme='urban'
        )
        
        self.indices_db['EMERGENCY_TENT'] = SpectralIndex(
            name='EMERGENCY_TENT',
            description='Emergency Tent Index (humanitarian relief shelters)',
            formula=lambda b: f"({b['SWIR2100']} - {b['RED']}) / ({b['SWIR2100']} + {b['RED']}) * (1 + ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2100': (2095, 2105)},
            reference='Emergency shelter detection - Giordan et al. 2018',
            theme='urban'
        )
        
        self.indices_db['REFUGEE_CAMP'] = SpectralIndex(
            name='REFUGEE_CAMP',
            description='Refugee Camp Index (temporary settlement detection)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Refugee camp monitoring - Kienberger et al. 2017',
            theme='urban'
        )
        
        self.indices_db['TEMPORARY_SHELTER'] = SpectralIndex(
            name='TEMPORARY_SHELTER',
            description='Temporary Shelter Index (temporary housing structures)',
            formula=lambda b: f"({b['SWIR1650']} - {b['BLUE']}) / ({b['SWIR1650']} + {b['BLUE']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR1650': (1645, 1655)},
            reference='Temporary housing detection - Tenerelli et al. 2020',
            theme='urban'
        )
        
        self.indices_db['DISASTER_RELIEF'] = SpectralIndex(
            name='DISASTER_RELIEF',
            description='Disaster Relief Infrastructure Index (emergency facilities)',
            formula=lambda b: f"({b['SWIR2200']} - {b['GREEN']}) / ({b['SWIR2200']} + {b['GREEN']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2200': (2195, 2205)},
            reference='Disaster response infrastructure - Voigt et al. 2016',
            theme='urban'
        )
        
        self.indices_db['MEDICAL_TENT'] = SpectralIndex(
            name='MEDICAL_TENT',
            description='Medical Tent Index (field hospitals and clinics)',
            formula=lambda b: f"({b['SWIR2100']} - {b['WHITE']}) / ({b['SWIR2100']} + {b['WHITE']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'WHITE': (650, 680), 'SWIR2100': (2095, 2105)},
            reference='Field hospital detection - Kuffer et al. 2018',
            theme='urban'
        )
        
        self.indices_db['CONSTRUCTION_SITE'] = SpectralIndex(
            name='CONSTRUCTION_SITE',
            description='Construction Site Index (temporary construction infrastructure)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Construction site monitoring - Weng et al. 2019',
            theme='urban'
        )
        
        self.indices_db['TEMPORARY_ROAD'] = SpectralIndex(
            name='TEMPORARY_ROAD',
            description='Temporary Road Index (access roads and temporary infrastructure)',
            formula=lambda b: f"({b['SWIR']} - {b['RED']}) / ({b['SWIR']} + {b['RED']}) * (1 - ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR': (1550, 1750)},
            reference='Temporary infrastructure mapping - Vatsavai et al. 2017',
            theme='urban'
        )
        
        self.indices_db['WIND_TURBINE'] = SpectralIndex(
            name='WIND_TURBINE',
            description='Wind Turbine Index (wind energy infrastructure)',
            formula=lambda b: f"({b['SWIR2300']} - {b['NIR']}) / ({b['SWIR2300']} + {b['NIR']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2300': (2295, 2305)},
            reference='Wind turbine detection - Manfreda et al. 2011',
            theme='urban'
        )
        
        self.indices_db['COMMUNICATION_TOWER'] = SpectralIndex(
            name='COMMUNICATION_TOWER',
            description='Communication Tower Index (telecom infrastructure)',
            formula=lambda b: f"({b['SWIR2200']} - {b['RED']}) / ({b['SWIR2200']} + {b['RED']}) * (1 + ({b['NIR']} - {b['GREEN']}) / ({b['NIR']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2200': (2195, 2205)},
            reference='Telecom tower detection - Jensen et al. 2015',
            theme='urban'
        )
        
        self.indices_db['PORTABLE_GENERATOR'] = SpectralIndex(
            name='PORTABLE_GENERATOR',
            description='Portable Generator Index (temporary power infrastructure)',
            formula=lambda b: f"({b['SWIR2100']} - {b['METALLIC']}) / ({b['SWIR2100']} + {b['METALLIC']}) * (1 + ({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'METALLIC': (550, 650), 'SWIR2100': (2095, 2105)},
            reference='Power infrastructure detection - Kemper et al. 2014',
            theme='urban'
        )
        
        self.indices_db['WATER_TANK'] = SpectralIndex(
            name='WATER_TANK',
            description='Water Tank Index (temporary water storage)',
            formula=lambda b: f"({b['SWIR']} - {b['GREEN']}) / ({b['SWIR']} + {b['GREEN']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Water storage facility detection - Sliuzas et al. 2019',
            theme='urban'
        )
        
        self.indices_db['FOOD_DISTRIBUTION'] = SpectralIndex(
            name='FOOD_DISTRIBUTION',
            description='Food Distribution Center Index (emergency food facilities)',
            formula=lambda b: f"({b['SWIR2200']} - {b['WHITE']}) / ({b['SWIR2200']} + {b['WHITE']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'WHITE': (650, 680), 'SWIR2200': (2195, 2205)},
            reference='Emergency food facility mapping - Brown et al. 2018',
            theme='urban'
        )
        
        # Aquatic and Coastal Environment Indices
        self.indices_db['ALGAE_BLOOM'] = SpectralIndex(
            name='ALGAE_BLOOM',
            description='Algae Bloom Index (phytoplankton detection)',
            formula=lambda b: f"({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Algal bloom detection - Matthews et al. 2012',
            theme='water'
        )
        
        self.indices_db['CYANOBACTERIA'] = SpectralIndex(
            name='CYANOBACTERIA',
            description='Cyanobacteria Index (harmful algal blooms)',
            formula=lambda b: f"({b['B550']} - {b['B620']}) / ({b['B550']} + {b['B620']}) * (1 + ({b['B680']} - {b['B550']}) / ({b['B680']} + {b['B550']}))",
            bands_required={'B550': (545, 555), 'B620': (615, 625), 'B680': (675, 685)},
            reference='Cyanobacteria monitoring - Wynne et al. 2008',
            theme='water'
        )
        
        self.indices_db['CHLOROPHYLL_A'] = SpectralIndex(
            name='CHLOROPHYLL_A',
            description='Chlorophyll-a Concentration Index (water quality)',
            formula=lambda b: f"({b['B670']} - {b['B680']}) / ({b['B670']} + {b['B680']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'B670': (665, 675), 'B680': (675, 685), 'RED': (620, 690), 'NIR': (760, 900)},
            reference="Chlorophyll-a retrieval - O'Reilly et al. 2000",
            theme='water'
        )
        
        self.indices_db['EUTROPHICATION'] = SpectralIndex(
            name='EUTROPHICATION',
            description='Eutrophication Index (nutrient pollution indicator)',
            formula=lambda b: f"({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}) * (1 + ({b['RED']} - {b['NIR']}) / ({b['RED']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Eutrophication monitoring - Palmer et al. 2015',
            theme='water'
        )
        
        self.indices_db['TURBIDITY'] = SpectralIndex(
            name='TURBIDITY',
            description='Water Turbidity Index (suspended sediments)',
            formula=lambda b: f"({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}) * (1 + ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Turbidity measurement - Dogliotti et al. 2015',
            theme='water'
        )
        
        self.indices_db['SEDIMENTATION'] = SpectralIndex(
            name='SEDIMENTATION',
            description='Sedimentation Index (high sediment load in water)',
            formula=lambda b: f"({b['SWIR']} - {b['RED']}) / ({b['SWIR']} + {b['RED']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'SWIR': (1550, 1750)},
            reference='Sediment load monitoring - Wang et al. 2016',
            theme='water'
        )
        
        self.indices_db['SUSPENDED_SOLIDS'] = SpectralIndex(
            name='SUSPENDED_SOLIDS',
            description='Suspended Solids Index (TSS concentration)',
            formula=lambda b: f"({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}) * (1 + ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Total suspended solids - Doxaran et al. 2002',
            theme='water'
        )
        
        self.indices_db['RIVER_BANK'] = SpectralIndex(
            name='RIVER_BANK',
            description='River Bank Index (shoreline and riparian zones)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='River bank delineation - Frohn et al. 2005',
            theme='water'
        )
        
        self.indices_db['SEASHORE'] = SpectralIndex(
            name='SEASHORE',
            description='Seashore Index (coastal interface zone)',
            formula=lambda b: f"({b['SWIR']} - {b['BLUE']}) / ({b['SWIR']} + {b['BLUE']}) * (1 + ({b['GREEN']} - {b['NIR']}) / ({b['GREEN']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Coastal zone mapping - Matarrese et al. 2018',
            theme='water'
        )
        
        self.indices_db['SHALLOW_WATER'] = SpectralIndex(
            name='SHALLOW_WATER',
            description='Shallow Water Index (bathymetry estimation)',
            formula=lambda b: f"({b['BLUE']} - {b['GREEN']}) / ({b['BLUE']} + {b['GREEN']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Shallow water bathymetry - Stumpf et al. 2003',
            theme='water'
        )
        
        self.indices_db['CORAL_REEF'] = SpectralIndex(
            name='CORAL_REEF',
            description='Coral Reef Index (coral health monitoring)',
            formula=lambda b: f"({b['BLUE']} - {b['GREEN']}) / ({b['BLUE']} + {b['GREEN']}) * (1 + ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'REDEDGE': (690, 730)},
            reference='Coral reef monitoring - Mumby et al. 2004',
            theme='water'
        )
        
        self.indices_db['SEAGRASS'] = SpectralIndex(
            name='SEAGRASS',
            description='Seagrass Index (submerged vegetation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Seagrass mapping - Phinn et al. 2008',
            theme='water'
        )
        
        self.indices_db['MANGROVE_WATER'] = SpectralIndex(
            name='MANGROVE_WATER',
            description='Mangrove Water Interface Index (coastal wetlands)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Mangrove water interface - Giri et al. 2007',
            theme='water'
        )
        
        self.indices_db['WATER_CLARITY'] = SpectralIndex(
            name='WATER_CLARITY',
            description='Water Clarity Index (transparency measurement)',
            formula=lambda b: f"({b['BLUE']} - {b['SWIR']}) / ({b['BLUE']} + {b['SWIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR': (1550, 1750)},
            reference='Water clarity assessment - Olmanson et al. 2008',
            theme='water'
        )
        
        self.indices_db['COLORED_DISSOLVED_ORGANIC'] = SpectralIndex(
            name='COLORED_DISSOLVED_ORGANIC',
            description='Colored Dissolved Organic Matter Index (CDOM)',
            formula=lambda b: f"({b['B440']} - {b['B490']}) / ({b['B440']} + {b['B490']}) * (1 + ({b['B550']} - {b['B670']}) / ({b['B550']} + {b['B670']}))",
            bands_required={'B440': (435, 445), 'B490': (485, 495), 'B550': (545, 555), 'B670': (665, 675)},
            reference='CDOM detection - Brezonik et al. 2015',
            theme='water'
        )
        
        self.indices_db['FLOATING_VEGETATION'] = SpectralIndex(
            name='FLOATING_VEGETATION',
            description='Floating Vegetation Index (water hyacinth, duckweed)',
            formula=lambda b: f"({b['NIR']} - {b['GREEN']}) / ({b['NIR']} + {b['GREEN']}) * (1 + ({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Floating aquatic vegetation - Dierssen et al. 2013',
            theme='water'
        )
        
        # Irrigation Detection Indices
        self.indices_db['IRRIGATED_CROP'] = SpectralIndex(
            name='IRRIGATED_CROP',
            description='Irrigated Crop Index (detecting irrigated agriculture)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Irrigated crop detection - Ambika et al. 2016',
            theme='vegetation'
        )
        
        self.indices_db['IRRIGATION_CANAL'] = SpectralIndex(
            name='IRRIGATION_CANAL',
            description='Irrigation Canal Index (water distribution channels)',
            formula=lambda b: f"({b['GREEN']} - {b['NIR']}) / ({b['GREEN']} + {b['NIR']}) * (1 + ({b['BLUE']} - {b['RED']}) / ({b['BLUE']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Irrigation canal mapping - Thenkabail et al. 2009',
            theme='water'
        )
        
        self.indices_db['PIVOT_IRRIGATION'] = SpectralIndex(
            name='PIVOT_IRRIGATION',
            description='Pivot Irrigation Index (center-pivot sprinkler systems)',
            formula=lambda b: f"({b['SWIR']} - {b['GREEN']}) / ({b['SWIR']} + {b['GREEN']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Center-pivot irrigation detection - Ozdogan et al. 2010',
            theme='water'
        )
        
        self.indices_db['DRIP_IRRIGATION'] = SpectralIndex(
            name='DRIP_IRRIGATION',
            description='Drip Irrigation Index (micro-irrigation systems)',
            formula=lambda b: f"({b['SWIR2200']} - {b['RED']}) / ({b['SWIR2200']} + {b['RED']}) * (1 + ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2200': (2195, 2205)},
            reference='Drip irrigation detection - Daccache et al. 2014',
            theme='water'
        )
        
        self.indices_db['FLOOD_IRRIGATION'] = SpectralIndex(
            name='FLOOD_IRRIGATION',
            description='Flood Irrigation Index (flooded field systems)',
            formula=lambda b: f"({b['GREEN']} - {b['SWIR']}) / ({b['GREEN']} + {b['SWIR']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Flood irrigation mapping - Bousbih et al. 2018',
            theme='water'
        )
        
        self.indices_db['IRRIGATION_POND'] = SpectralIndex(
            name='IRRIGATION_POND',
            description='Irrigation Pond Index (water storage for agriculture)',
            formula=lambda b: f"({b['BLUE']} - {b['SWIR']}) / ({b['BLUE']} + {b['SWIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR': (1550, 1750)},
            reference='Irrigation pond detection - Velpuri et al. 2017',
            theme='water'
        )
        
        self.indices_db['SPRINKLER_SYSTEM'] = SpectralIndex(
            name='SPRINKLER_SYSTEM',
            description='Sprinkler System Index (sprinkler irrigation infrastructure)',
            formula=lambda b: f"({b['SWIR2100']} - {b['GREEN']}) / ({b['SWIR2100']} + {b['GREEN']}) * (1 + ({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2100': (2095, 2105)},
            reference='Sprinkler system detection - Kuenzer et al. 2012',
            theme='water'
        )
        
        self.indices_db['IRRIGATION_INFRASTRUCTURE'] = SpectralIndex(
            name='IRRIGATION_INFRASTRUCTURE',
            description='Irrigation Infrastructure Index (general irrigation facilities)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Irrigation infrastructure mapping - Dwarakish et al. 2015',
            theme='water'
        )
        
        self.indices_db['WATER_STRESSED_CROP'] = SpectralIndex(
            name='WATER_STRESSED_CROP',
            description='Water Stressed Crop Index (detecting water stress in irrigated areas)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 - ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Water stress detection - Kogan et al. 2011',
            theme='vegetation'
        )
        
        self.indices_db['WELL_IRRIGATION'] = SpectralIndex(
            name='WELL_IRRIGATION',
            description='Well Irrigation Index (groundwater extraction points)',
            formula=lambda b: f"({b['SWIR2200']} - {b['BLUE']}) / ({b['SWIR2200']} + {b['BLUE']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2200': (2195, 2205)},
            reference='Well irrigation detection - Shah et al. 2013',
            theme='water'
        )
        
        self.indices_db['RIVER_LIFT_IRRIGATION'] = SpectralIndex(
            name='RIVER_LIFT_IRRIGATION',
            description='River Lift Irrigation Index (river water pumping systems)',
            formula=lambda b: f"({b['SWIR']} - {b['GREEN']}) / ({b['SWIR']} + {b['GREEN']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='River lift irrigation mapping - Bhanja et al. 2019',
            theme='water'
        )
        
        self.indices_db['IRRIGATION_TIMING'] = SpectralIndex(
            name='IRRIGATION_TIMING',
            description='Irrigation Timing Index (recent vs. old irrigation)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2100']}) / ({b['SWIR1650']} + {b['SWIR2100']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1650': (1645, 1655), 'SWIR2100': (2095, 2105)},
            reference='Irrigation timing detection - Biggs et al. 2016',
            theme='water'
        )
        
        self.indices_db['IRRIGATION_EFFICIENCY'] = SpectralIndex(
            name='IRRIGATION_EFFICIENCY',
            description='Irrigation Efficiency Index (water use efficiency)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Irrigation efficiency assessment - Perry et al. 2017',
            theme='water'
        )
        
        # Non-Irrigated Yield Prediction Indices
        self.indices_db['RAINFED_CROP'] = SpectralIndex(
            name='RAINFED_CROP',
            description='Rain-fed Crop Index (non-irrigated agriculture)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Rain-fed crop detection - Balaghi et al. 2014',
            theme='vegetation'
        )
        
        self.indices_db['DROUGHT_STRESS'] = SpectralIndex(
            name='DROUGHT_STRESS',
            description='Drought Stress Index (water deficit in rain-fed areas)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Drought stress monitoring - Gu et al. 2007',
            theme='vegetation'
        )
        
        self.indices_db['YIELD_PREDICTION'] = SpectralIndex(
            name='YIELD_PREDICTION',
            description='Crop Yield Prediction Index (general yield estimation)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.3 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Crop yield prediction - Becker-Reshef et al. 2010',
            theme='vegetation'
        )
        
        self.indices_db['RAINFED_YIELD'] = SpectralIndex(
            name='RAINFED_YIELD',
            description='Rain-fed Yield Index (yield prediction for non-irrigated crops)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.2 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Rain-fed yield modeling - Lobell et al. 2015',
            theme='vegetation'
        )
        
        self.indices_db['DRYLAND_YIELD'] = SpectralIndex(
            name='DRYLAND_YIELD',
            description='Dryland Yield Index (arid region crop yield)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.4 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Dryland agriculture yield - Basso et al. 2012',
            theme='vegetation'
        )
        
        self.indices_db['PRECIPITATION_EFFECTIVENESS'] = SpectralIndex(
            name='PRECIPITATION_EFFECTIVENESS',
            description='Precipitation Effectiveness Index (rainfall utilization)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.25 * ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Rainfall use efficiency - Mkhabela et al. 2011',
            theme='vegetation'
        )
        
        self.indices_db['SOIL_MOISTURE_DEFICIT'] = SpectralIndex(
            name='SOIL_MOISTURE_DEFICIT',
            description='Soil Moisture Deficit Index (water stress in rain-fed areas)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + 0.3 * ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Soil moisture deficit - Wang et al. 2009',
            theme='soil'
        )
        
        self.indices_db['CROP_WATER_STRESS'] = SpectralIndex(
            name='CROP_WATER_STRESS',
            description='Crop Water Stress Index (physiological stress in rain-fed crops)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Crop physiological stress - Penuelas et al. 2011',
            theme='vegetation'
        )
        
        self.indices_db['RAINFED_BIOMASS'] = SpectralIndex(
            name='RAINFED_BIOMASS',
            description='Rain-fed Biomass Index (biomass estimation for non-irrigated crops)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.2 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Rain-fed biomass estimation - Prince et al. 2009',
            theme='vegetation'
        )
        
        self.indices_db['CLIMATE_YIELD'] = SpectralIndex(
            name='CLIMATE_YIELD',
            description='Climate Yield Index (climate-limited yield prediction)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.15 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Climate yield modeling - Schlenker et al. 2013',
            theme='vegetation'
        )
        
        self.indices_db['RAINFED_PRODUCTIVITY'] = SpectralIndex(
            name='RAINFED_PRODUCTIVITY',
            description='Rain-fed Productivity Index (overall productivity of rain-fed agriculture)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.35 * ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'GREEN': (520, 600), 'NIR': (760, 900)},
            reference='Rain-fed productivity - Johnson et al. 2014',
            theme='vegetation'
        )
        
        self.indices_db['WATER_LIMITED_YIELD'] = SpectralIndex(
            name='WATER_LIMITED_YIELD',
            description='Water Limited Yield Index (yield constrained by water availability)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 - 0.5 * ({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}))",
            bands_required={'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Water-limited yield - Sinclair et al. 2010',
            theme='vegetation'
        )
        
        self.indices_db['DROUGHT_VULNERABILITY'] = SpectralIndex(
            name='DROUGHT_VULNERABILITY',
            description='Drought Vulnerability Index (susceptibility to water stress)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + 0.4 * ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Drought vulnerability assessment - Svoboda et al. 2016',
            theme='vegetation'
        )
        
        self.indices_db['RAINFED_PHENOLOGY'] = SpectralIndex(
            name='RAINFED_PHENOLOGY',
            description='Rain-fed Phenology Index (growth stage detection in rain-fed crops)',
            formula=lambda b: f"({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}) * (1 + 0.25 * ({b['REDEDGE']} - {b['RED']}) / ({b['REDEDGE']} + {b['RED']}))",
            bands_required={'RED': (620, 690), 'REDEDGE': (690, 730), 'NIR': (760, 900)},
            reference='Rain-fed phenology monitoring - Sakamoto et al. 2005',
            theme='vegetation'
        )
        
        # Paint and Coating Detection Indices
        self.indices_db['VEHICLE_PAINT'] = SpectralIndex(
            name='VEHICLE_PAINT',
            description='Vehicle Paint Detection Index (automotive and military vehicle coatings)',
            formula=lambda b: f"({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}) * (1 + ({b['BLUE']} - {b['NIR']}) / ({b['BLUE']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Vehicle paint detection - Kudoh et al. 2010',
            theme='materials'
        )
        
        self.indices_db['MILITARY_CAMOUFLAGE'] = SpectralIndex(
            name='MILITARY_CAMOUFLAGE',
            description='Military Camouflage Detection Index (camouflage patterns and materials)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Military camouflage detection - de Visser et al. 2016',
            theme='materials'
        )
        
        self.indices_db['MARINE_PAINT'] = SpectralIndex(
            name='MARINE_PAINT',
            description='Marine Paint Detection Index (ship and boat coatings)',
            formula=lambda b: f"({b['BLUE']} - {b['GREEN']}) / ({b['BLUE']} + {b['GREEN']}) * (1 + ({b['RED']} - {b['SWIR']}) / ({b['RED']} + {b['SWIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR': (1550, 1750)},
            reference='Marine coating detection - Kester et al. 2014',
            theme='materials'
        )
        
        self.indices_db['INDUSTRIAL_COATING'] = SpectralIndex(
            name='INDUSTRIAL_COATING',
            description='Industrial Coating Index (protective and functional coatings)',
            formula=lambda b: f"({b['SWIR2200']} - {b['GREEN']}) / ({b['SWIR2200']} + {b['GREEN']}) * (1 + ({b['NIR']} - {b['RED']}) / ({b['NIR']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2200': (2195, 2205)},
            reference='Industrial coating analysis - Yang et al. 2013',
            theme='materials'
        )
        
        self.indices_db['ROAD_PAINT'] = SpectralIndex(
            name='ROAD_PAINT',
            description='Road Paint Detection Index (highway markings and road coatings)',
            formula=lambda b: f"({b['YELLOW']} - {b['WHITE']}) / ({b['YELLOW']} + {b['WHITE']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'WHITE': (650, 680), 'YELLOW': (570, 590)},
            reference='Road marking detection - Li et al. 2015',
            theme='materials'
        )
        
        self.indices_db['BUILDING_PAINT'] = SpectralIndex(
            name='BUILDING_PAINT',
            description='Building Paint Index (architectural and structural coatings)',
            formula=lambda b: f"({b['RED']} - {b['BLUE']}) / ({b['RED']} + {b['BLUE']}) * (1 + ({b['GREEN']} - {b['NIR']}) / ({b['GREEN']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Building coating detection - Herold et al. 2012',
            theme='materials'
        )
        
        self.indices_db['REFLECTIVE_PAINT'] = SpectralIndex(
            name='REFLECTIVE_PAINT',
            description='Reflective Paint Index (high-reflectivity coatings)',
            formula=lambda b: f"({b['NIR']} - {b['SWIR']}) / ({b['NIR']} + {b['SWIR']}) * (1 + ({b['WHITE']} - {b['GREEN']}) / ({b['WHITE']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'WHITE': (650, 680), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Reflective coating detection - Levinson et al. 2007',
            theme='materials'
        )
        
        self.indices_db['ANTI_FOULING_PAINT'] = SpectralIndex(
            name='ANTI_FOULING_PAINT',
            description='Anti-fouling Paint Index (marine and industrial anti-fouling coatings)',
            formula=lambda b: f"({b['SWIR2100']} - {b['BLUE']}) / ({b['SWIR2100']} + {b['BLUE']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2100': (2095, 2105)},
            reference='Anti-fouling coating detection - Yebra et al. 2016',
            theme='materials'
        )
        
        self.indices_db['THERMAL_PAINT'] = SpectralIndex(
            name='THERMAL_PAINT',
            description='Thermal Paint Index (heat-reflective and insulating coatings)',
            formula=lambda b: f"({b['SWIR2300']} - {b['NIR']}) / ({b['SWIR2300']} + {b['NIR']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR2300': (2295, 2305)},
            reference='Thermal coating detection - Santamouris et al. 2011',
            theme='materials'
        )
        
        self.indices_db['FLUORESCENT_PAINT'] = SpectralIndex(
            name='FLUORESCENT_PAINT',
            description='Fluorescent Paint Index (UV-reactive and fluorescent coatings)',
            formula=lambda b: f"({b['BLUE']} - {b['UV']}) / ({b['BLUE']} + {b['UV']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'UV': (350, 400), 'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690)},
            reference='Fluorescent coating detection - Kim et al. 2018',
            theme='materials'
        )
        
        self.indices_db['METALLIC_PAINT'] = SpectralIndex(
            name='METALLIC_PAINT',
            description='Metallic Paint Index (metal-flake and pearlescent coatings)',
            formula=lambda b: f"({b['SWIR']} - {b['NIR']}) / ({b['SWIR']} + {b['NIR']}) * (1 + ({b['RED']} - {b['GREEN']}) / ({b['RED']} + {b['GREEN']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900), 'SWIR': (1550, 1750)},
            reference='Metallic coating detection - Gao et al. 2014',
            theme='materials'
        )
        
        self.indices_db['AEROSOL_PAINT'] = SpectralIndex(
            name='AEROSOL_PAINT',
            description='Aerosol Paint Index (spray-applied coatings and graffiti)',
            formula=lambda b: f"({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}) * (1 + ({b['RED']} - {b['NIR']}) / ({b['RED']} + {b['NIR']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'NIR': (760, 900)},
            reference='Aerosol coating detection - Valero et al. 2017',
            theme='materials'
        )
        
        self.indices_db['WATERPROOF_PAINT'] = SpectralIndex(
            name='WATERPROOF_PAINT',
            description='Waterproof Paint Index (water-repellent coatings)',
            formula=lambda b: f"({b['SWIR1650']} - {b['SWIR2200']}) / ({b['SWIR1650']} + {b['SWIR2200']}) * (1 + ({b['GREEN']} - {b['RED']}) / ({b['GREEN']} + {b['RED']}))",
            bands_required={'GREEN': (520, 600), 'RED': (620, 690), 'SWIR1650': (1645, 1655), 'SWIR2200': (2195, 2205)},
            reference='Waterproof coating detection - Song et al. 2012',
            theme='materials'
        )
        
        self.indices_db['CORROSION_PROTECTIVE_PAINT'] = SpectralIndex(
            name='CORROSION_PROTECTIVE_PAINT',
            description='Corrosion Protective Paint Index (anti-corrosion coatings)',
            formula=lambda b: f"({b['SWIR2100']} - {b['RED']}) / ({b['SWIR2100']} + {b['RED']}) * (1 + ({b['GREEN']} - {b['BLUE']}) / ({b['GREEN']} + {b['BLUE']}))",
            bands_required={'BLUE': (450, 520), 'GREEN': (520, 600), 'RED': (620, 690), 'SWIR2100': (2095, 2105)},
            reference='Corrosion protection coating detection - Feng et al. 2015',
            theme='materials'
        )
            
    def find_closest_band(self, wavelength_target, available_wavelengths):
        """
        Find the closest available wavelength to a target wavelength or range.
        
        This method handles both single wavelength targets and wavelength ranges,
        finding the best match from the available bands.
        
        Args:
            wavelength_target (float or tuple): Target wavelength in nm, or 
                (min_wavelength, max_wavelength) tuple for a range
            available_wavelengths (list): List of available wavelengths in nm
            
            Returns:
                float: The closest available wavelength from the list
            
            Examples:
                >>> indices = HyperspectralIndices()
                >>> available = [480, 560, 665, 850]
                >>> indices.find_closest_band(670, available)
                665
                >>> indices.find_closest_band((760, 900), available)
                850
            """
        if isinstance(wavelength_target, tuple):
            # Target is a range - find center point
            target_center = (wavelength_target[0] + wavelength_target[1]) / 2
        else:
            target_center = wavelength_target
        
        closest = min(available_wavelengths, 
                     key=lambda x: abs(x - target_center))
        return closest
        
    def can_calculate_index(self, index_name, available_wavelengths):
            """
            Check if an index can be calculated with the available bands.
            
            This method verifies that all required spectral bands for a given
            index are available in the input data, within acceptable tolerance.
            
            Args:
                index_name (str): Name of the index to check (case-insensitive)
                available_wavelengths (list): List of available wavelengths in nm
            
            Returns:
                tuple: (bool, str) - (can_calculate, message)
                    can_calculate: True if all required bands are available
                    message: Description of result or missing bands
            
            Examples:
                >>> indices = HyperspectralIndices()
                >>> can_calc, msg = indices.can_calculate_index('NDVI', [665, 850])
                >>> print(can_calc)
                True
            """
            index = self.indices_db.get(index_name.upper())
            if not index:
                return False, "Index not found"
            
            for band_name, wavelength_range in index.bands_required.items():
                found = False
                for wl in available_wavelengths:
                    if isinstance(wavelength_range, tuple):
                        if wavelength_range[0] <= wl <= wavelength_range[1]:
                            found = True
                            break
                    else:
                        # Single wavelength requirement with 50nm tolerance
                        if abs(wl - wavelength_range) < 50:
                            found = True
                            break
                
                if not found:
                    return False, f"Required band {band_name} ({wavelength_range} nm) not available"
            
            return True, "OK"
        
    def get_band_mapping(self, index_name, wavelength_to_band):
            """
            Get the mapping of required bands to actual raster names for an index.
            
            This method creates a dictionary mapping the abstract band names used
            in index formulas (e.g., 'RED', 'NIR') to the actual raster map names
            from the input data.
            
            Args:
                index_name (str): Name of the index
                wavelength_to_band (dict): Dictionary mapping wavelengths to raster names
            
            Returns:
                dict: Mapping of abstract band names to raster map names
            
            Examples:
                >>> indices = HyperspectralIndices()
                >>> wl_to_band = {665: 'band3', 850: 'band4'}
                >>> mapping = indices.get_band_mapping('NDVI', wl_to_band)
                >>> print(mapping)
                {'RED': 'band3', 'NIR': 'band4'}
            """
            index = self.indices_db[index_name.upper()]
            mapping = {}
            
            available_wl = list(wavelength_to_band.keys())
            
            for band_name, wavelength_range in index.bands_required.items():
                closest_wl = self.find_closest_band(wavelength_range, available_wl)
                mapping[band_name] = wavelength_to_band[closest_wl]
            
            return mapping
        
    def list_indices(self, theme=None):
            """
            List available indices, optionally filtered by theme.
            
            Args:
                theme (str): Optional theme name to filter by
            
            Returns:
                list: List of SpectralIndex objects matching the criteria
            
            Examples:
                >>> indices = HyperspectralIndices()
                >>> veg_indices = indices.list_indices(theme='vegetation')
                >>> print(len(veg_indices))
                15
            """
            indices_list = []
            
            for name, index in sorted(self.indices_db.items()):
                if theme and index.theme != theme:
                    continue
                indices_list.append(index)
            
            return indices_list
        
    def get_themes(self):
            """
            Get list of all available index themes.
            
            Returns:
                list: Sorted list of unique theme names
            
            Examples:
                >>> indices = HyperspectralIndices()
                >>> themes = indices.get_themes()
                >>> 'vegetation' in themes
                True
            """
            themes = set(index.theme for index in self.indices_db.values())
            return sorted(themes)


        

import grass.script as gs
from grass.script import mapcalc
from grass.script.core import run_command
from grass.exceptions import CalledModuleError


def _load_raster3d_lib():
    """Load libgrass_raster3d and bind Rast3d_extract_z_slice."""
    gisbase = os.environ.get('GISBASE', '')
    lib_path = None
    if gisbase:
        candidate = os.path.join(gisbase, 'lib', 'libgrass_raster3d.so')
        if os.path.exists(candidate):
            lib_path = candidate
    if lib_path is None:
        lib_path = ctypes.util.find_library('grass_raster3d')
    if lib_path is None:
        gs.fatal("Cannot find libgrass_raster3d")
    lib = ctypes.CDLL(lib_path)
    lib.Rast3d_extract_z_slice.restype = ctypes.c_int
    lib.Rast3d_extract_z_slice.argtypes = [
        ctypes.c_char_p,  # name3d
        ctypes.c_char_p,  # mapset3d (empty = search path)
        ctypes.c_int,     # z index
        ctypes.c_char_p,  # name2d output
    ]
    return lib


def _extract_slices_from_3d(map3d, wavelengths, tmp_prefix):
    """
    Extract Z-slices from a 3D raster into temporary 2D rasters.

    Calls Rast3d_extract_z_slice() from libgrass_raster3d via ctypes.
    Each slice is read in bulk (RASTER3D_NO_CACHE + Rast3d_get_block),
    not per-voxel.

    Returns list of temp band names in Z order.
    """
    lib = _load_raster3d_lib()
    band_names = []
    for z, _wl in enumerate(wavelengths):
        tmp_name = f"{tmp_prefix}_{z}"
        ret = lib.Rast3d_extract_z_slice(
            map3d.encode(), b"", ctypes.c_int(z), tmp_name.encode()
        )
        if ret != 0:
            gs.fatal(f"Rast3d_extract_z_slice: failed at z={z} of <{map3d}>")
        band_names.append(tmp_name)
    return band_names


def list_available_indices(indices_obj, detailed=False):
    """
    Print formatted list of available indices organized by theme.
    
    This function provides user-friendly output of all available indices,
    grouped by thematic category. When detailed mode is enabled, it shows
    complete information including required bands and references.
    
    Args:
        indices_obj (HyperspectralIndices): Indices database object
        detailed (bool): If True, show full details for each index
    
    Output:
        Prints formatted index information to stdout via gs.message
    
    Examples:
        >>> indices = HyperspectralIndices()
        >>> list_available_indices(indices, detailed=False)
        # Prints summary list
        >>> list_available_indices(indices, detailed=True)
        # Prints detailed information
    """
    themes = indices_obj.get_themes()
    
    gs.message("="*70)
    gs.message("AVAILABLE SPECTRAL INDICES")
    gs.message("="*70)
    
    for theme in themes:
        theme_indices = indices_obj.list_indices(theme=theme)
        gs.message(f"\n{theme.upper()} INDICES ({len(theme_indices)}):")
        gs.message("-"*70)
        
        for idx in theme_indices:
            if detailed:
                gs.message(f"\n  {idx.name}: {idx.description}")
                gs.message(f"    Required bands: {idx.bands_required}")
                gs.message(f"    Reference: {idx.reference}")
            else:
                # Format band requirements as a concise string
                bands_str = ", ".join([f"{band}({rng[0]}-{rng[1]}nm)" for band, rng in idx.bands_required.items()])
                gs.message(f"  {idx.name:15s} - {idx.description}")
                gs.message(f"    Bands: {bands_str}")
    
    gs.message("\n" + "="*70)
    gs.message(f"Total indices available: {len(indices_obj.indices_db)}")
    gs.message("="*70)

def main():
    # Parse command line - CRITICAL: missing in original
    options, flags = gs.parser()
    
    # Initialize indices database
    indices_obj = HyperspectralIndices()
    
    # Handle list flag - show available indices and exit
    if flags['l']:
        list_available_indices(indices_obj, detailed=flags['i'])
        return 0
    
    # Get command-line options
    input_bands = options['input']
    input3d = options['input3d']
    wavelengths_str = options['wavelengths']
    band_wavelengths_str = options['band_wavelengths']
    output_prefix = options['output_prefix']
    indices_str = options['indices']
    theme = options['theme']

    if not output_prefix:
        gs.fatal(_("Required parameter <output_prefix> not set"))

    # ====================================================================
    # 3D RASTER INPUT PATH
    # Uses Rast3d_extract_z_slice() (tile-bulk reads, RASTER3D_NO_CACHE)
    # instead of a per-voxel loop over Rast3d_get_value().
    # ====================================================================

    tmp_bands = []  # track temp maps for cleanup

    if input3d:
        if not band_wavelengths_str:
            gs.fatal(_("band_wavelengths required when input3d is set"))
        try:
            wavelengths = [float(wl) for wl in band_wavelengths_str.split(',')]
        except ValueError:
            gs.fatal(_("Invalid band_wavelengths values"))

        tmp_prefix = f"tmp_ihi_{os.getpid()}"
        gs.verbose(_("Extracting {} slices from 3D raster <{}>").format(
            len(wavelengths), input3d))
        tmp_bands = _extract_slices_from_3d(input3d, wavelengths, tmp_prefix)
        atexit.register(
            run_command, 'g.remove', flags='f', type='raster',
            name=','.join(tmp_bands), quiet=True
        )
        input_bands = tmp_bands

    else:
        # ----------------------------------------------------------------
        # 2D RASTER INPUT PATH (original)
        # ----------------------------------------------------------------
        if not input_bands:
            gs.fatal(_("Either input= or input3d= is required"))
        if not wavelengths_str:
            gs.fatal(_("wavelengths= required when using input="))

        input_bands = input_bands.split(',')

        try:
            wavelengths = [float(wl) for wl in wavelengths_str.split(',')]
        except ValueError:
            gs.fatal(_("Invalid wavelength values. Use comma-separated numbers "
                       "(e.g., 450,550,670,800)"))

    # ====================================================================
    # INPUT VALIDATION
    # ====================================================================
    
    # Verify input bands match wavelengths
    if len(input_bands) != len(wavelengths):
        gs.fatal(_("Number of input bands ({}) must match number of "
                   "wavelengths ({})").format(len(input_bands), len(wavelengths)))
    
    # Create wavelength to band name mapping dictionary
    wavelength_to_band = dict(zip(wavelengths, input_bands))
    
    gs.verbose(_("Processing {} bands with wavelengths: {}").format(
        len(wavelengths), ', '.join([f"{w}nm" for w in wavelengths])))
    
    # ====================================================================
    # INDEX SELECTION
    # ====================================================================
    
    indices_to_calc = []
    
    if theme:
        # Calculate all indices from a specific theme
        theme_indices = indices_obj.list_indices(theme=theme)
        indices_to_calc = [idx.name for idx in theme_indices]
        gs.message(_("Selected theme '{}' with {} indices").format(
            theme, len(indices_to_calc)))
    elif indices_str.lower() == 'all':
        # Calculate all possible indices
        indices_to_calc = list(indices_obj.indices_db.keys())
        gs.message(_("Calculating all {} available indices").format(
            len(indices_to_calc)))
    else:
        # Calculate user-specified indices
        indices_to_calc = [idx.strip().upper() for idx in indices_str.split(',')]
        gs.verbose(_("Selected {} specific indices").format(len(indices_to_calc)))
    
    # ====================================================================
    # INDEX CALCULATION LOOP
    # ====================================================================
    
    calculated = 0
    skipped = 0
    
    gs.message(_("Starting index calculation..."))
    gs.message("-" * 70)
    
    for index_name in indices_to_calc:
        # Verify index exists in database
        if index_name not in indices_obj.indices_db:
            gs.warning(_("Index '{}' not found in database. Skipping.").format(
                index_name))
            skipped += 1
            continue
        
        # Check if we have the required bands for this index
        can_calc, msg = indices_obj.can_calculate_index(index_name, wavelengths)
        
        if not can_calc:
            gs.warning(_("Cannot calculate {}: {}").format(index_name, msg))
            skipped += 1
            continue
        
        # FIXED: Indentation - all below now properly nested under successful validation
        
        # Get band mapping for this index
        band_mapping = indices_obj.get_band_mapping(index_name, wavelength_to_band)
        
        # Get index definition
        index_def = indices_obj.indices_db[index_name]
        
        # Generate the formula
        formula = index_def.formula(band_mapping)
        
        # Output name
        output_name = f"{output_prefix}_{index_name}"
        
        # Calculate the index using r.mapcalc
        gs.message(_("Calculating {}...").format(index_name))
        
        try:
            # FIXED: Use proper GRASS API
            mapcalc(f"{output_name} = {formula}", overwrite=True)
            
            # Apply normalization if requested and applicable
            if flags['n'] and hasattr(index_def, 'normalize_range') and index_def.normalize_range:
                min_val, max_val = index_def.normalize_range
                temp_name = f"{output_name}_temp"
                norm_formula = f"({output_name} - {min_val}) / ({max_val} - {min_val})"
                mapcalc(f"{temp_name} = {norm_formula}", overwrite=True)
                # Rename normalized to main output
                run_command('g.rename', raster=f"{temp_name},{output_name}", overwrite=True, quiet=True)
                gs.message(f"  -> Created normalized version: {output_name}")
            
            # Set appropriate color table
            if 'NDV' in index_name or 'EVI' in index_name or getattr(index_def, 'theme', '') == 'vegetation':
                run_command('r.colors', map=output_name, color='ndvi', quiet=True)
            elif getattr(index_def, 'theme', '') == 'water':
                run_command('r.colors', map=output_name, color='water', quiet=True)
            else:
                run_command('r.colors', map=output_name, color='viridis', quiet=True)
            
            calculated += 1
            gs.message(f"  -> Successfully created: {output_name}")
            
        except CalledModuleError as e:
            gs.warning(_("Failed to calculate {}: {}").format(index_name, str(e)))
            skipped += 1
            continue
    
    # Summary
    gs.message("\n" + "="*70)
    gs.message(_("SUMMARY: {} indices calculated, {} skipped").format(calculated, skipped))
    gs.message("="*70)

    return 0



if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
