"""Tests that libras3d is functional for i.hyper.indices standalone mode."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, TANAGER_PATH,
    skip_without_ras3d, skip_without_wyvern, skip_without_tanager,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)

@skip_without_ras3d
def test_shim_installs():
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'raster3d_info') and hasattr(gs, 'parser')

@skip_without_ras3d
@skip_without_wyvern
def test_open_wyvern_geotiff():
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    assert r['depths'] == 23 and r['rows'] == 7825 and r['cols'] == 6003
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_tanager
def test_open_tanager_hdf5():
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    assert r['depths'] == 426 and r['rows'] == 732 and r['cols'] == 607
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_wyvern
def test_read_all_bands_wyvern():
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    cube = ras3d.read_all_bands(h)
    assert cube.shape == (r['depths'], r['rows'], r['cols'])
    assert cube.dtype == np.float32
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_tanager
def test_read_all_bands_tanager():
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    cube = ras3d.read_all_bands(h)
    assert_band_valid(cube[0], 'Tanager band 0')
    assert_band_valid(cube[425], 'Tanager band 425')
    ras3d.close_cube(h)

@skip_without_ras3d
@skip_without_wyvern
def test_extract_slices_ras3d(tmp_path):
    install_ras3d_shim()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)
    sys.path.insert(0, '/home/yann/dev/i.hyper.indices')
    import i_hyper_indices
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    wl = list(range(r['depths']))
    ras3d.close_cube(h)
    names = i_hyper_indices._extract_slices_from_3d(WYVERN_PATH, wl, '_tst')
    assert len(names) == r['depths']
    from ras3d_grass_shim import get_band_cache
    for n in names[:3]:
        assert n in get_band_cache()
        assert_band_valid(get_band_cache()[n], f'indices slice {n}')

@skip_without_ras3d
@skip_without_wyvern
def test_get_region_from_shim():
    install_ras3d_shim()
    import grass.script as gs
    info = gs.raster3d_info(WYVERN_PATH)
    assert info['depths'] == 23 and info['rows'] == 7825 and info['cols'] == 6003
