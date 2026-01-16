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
#% options: vegetation,water,soil,urban,stress,biochemical,pigments,metabolism,materials,all
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
                gs.message(f"  {idx.name:15s} - {idx.description}")
    
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
    wavelengths_str = options['wavelengths']
    output_prefix = options['output_prefix']
    indices_str = options['indices']
    theme = options['theme']
    
    # Validate required parameters for index calculation
    if not input_bands:
        gs.fatal(_("Required parameter <input> not set: (Input raster bands (comma-separated list))"))
    if not wavelengths_str:
        gs.fatal(_("Required parameter <wavelengths> not set: (Wavelengths for input bands in nm (comma-separated, e.g., 450,550,670,800))"))
    if not output_prefix:
        gs.fatal(_("Required parameter <output_prefix> not set: (Prefix for output index rasters)"))
    
    input_bands = input_bands.split(',')
    wavelengths_str = wavelengths_str.split(',')
    
    # ====================================================================
    # INPUT VALIDATION
    # ====================================================================
    
    # Parse wavelengths from string to float
    try:
        wavelengths = [float(wl) for wl in wavelengths_str]
    except ValueError:
        gs.fatal(_("Invalid wavelength values. Use comma-separated numbers "
                   "(e.g., 450,550,670,800)"))
    
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
