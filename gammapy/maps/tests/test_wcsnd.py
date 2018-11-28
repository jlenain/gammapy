# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from ...utils.testing import requires_dependency
from ..utils import fill_poisson
from ..geom import MapAxis, MapCoord, coordsys_to_frame
from ..base import Map
from ..wcs import WcsGeom
from ..hpx import HpxGeom
from ..wcsnd import WcsNDMap

pytest.importorskip('scipy')
pytest.importorskip('reproject')

axes1 = [MapAxis(np.logspace(0., 3., 3), interp='log', name='spam')]
axes2 = [MapAxis(np.logspace(0., 3., 3), interp='log'),
         MapAxis(np.logspace(1., 3., 4), interp='lin')]
skydir = SkyCoord(110., 75.0, unit='deg', frame='icrs')

wcs_allsky_test_geoms = [
    (None, 10.0, 'GAL', 'AIT', skydir, None),
    (None, 10.0, 'GAL', 'AIT', skydir, axes1),
    (None, [10.0, 20.0], 'GAL', 'AIT', skydir, axes1),
    (None, 10.0, 'GAL', 'AIT', skydir, axes2),
    (None, [[10.0, 20.0, 30.0], [10.0, 20.0, 30.0]],
     'GAL', 'AIT', skydir, axes2),
]

wcs_partialsky_test_geoms = [
    (10, 1.0, 'GAL', 'AIT', skydir, None),
    (10, 1.0, 'GAL', 'AIT', skydir, axes1),
    (10, [1.0, 2.0], 'GAL', 'AIT', skydir, axes1),
    (10, 1.0, 'GAL', 'AIT', skydir, axes2),
    (10, [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]],
     'GAL', 'AIT', skydir, axes2),
]

wcs_test_geoms = wcs_allsky_test_geoms + wcs_partialsky_test_geoms


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_init(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m0 = WcsNDMap(geom)
    coords = m0.geom.get_coord()
    m0.set_by_coord(coords, coords[1])
    m1 = WcsNDMap(geom, m0.data)
    assert_allclose(m0.data, m1.data)


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_read_write(tmpdir, npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    filename = str(tmpdir / 'map.fits')

    m0 = WcsNDMap(geom)
    fill_poisson(m0, mu=0.5)
    m0.write(filename, overwrite=True)
    m1 = WcsNDMap.read(filename)
    m2 = Map.read(filename)
    m3 = Map.read(filename, map_type='wcs')
    assert_allclose(m0.data, m1.data)
    assert_allclose(m0.data, m2.data)
    assert_allclose(m0.data, m3.data)

    m0.write(filename, sparse=True, overwrite=True)
    m1 = WcsNDMap.read(filename)
    m2 = Map.read(filename)
    m3 = Map.read(filename, map_type='wcs')
    assert_allclose(m0.data, m1.data)
    assert_allclose(m0.data, m2.data)
    assert_allclose(m0.data, m3.data)

    # Specify alternate HDU name for IMAGE and BANDS table
    m0.write(filename, hdu='IMAGE', hdu_bands='TEST', overwrite=True)
    m1 = WcsNDMap.read(filename)
    m2 = Map.read(filename)
    m3 = Map.read(filename, map_type='wcs')


def test_wcsndmap_read_write_fgst(tmpdir):
    filename = str(tmpdir / 'map.fits')

    axis = MapAxis.from_bounds(100., 1000., 4, name='energy', unit='MeV')
    geom = WcsGeom.create(npix=10, binsz=1.0,
                          proj='AIT', coordsys='GAL', axes=[axis])

    # Test Counts Cube
    m = WcsNDMap(geom)
    m.write(filename, conv='fgst-ccube', overwrite=True)
    with fits.open(filename) as h:
        assert 'EBOUNDS' in h

    m2 = Map.read(filename)
    assert m2.geom.conv == 'fgst-ccube'

    # Test Model Cube
    m.write(filename, conv='fgst-template', overwrite=True)
    with fits.open(filename) as h:
        assert 'ENERGIES' in h

    m2 = Map.read(filename)
    assert m2.geom.conv == 'fgst-template'


def test_wcs_nd_map_data_transpose_issue(tmpdir):
    # Regression test for https://github.com/gammapy/gammapy/issues/1346

    # Our test case: a little map with WCS shape (3, 2), i.e. numpy array shape (2, 3)
    data = np.array([[0, 1, 2], [np.nan, np.inf, -np.inf]])
    geom = WcsGeom.create(npix=(3, 2))

    # Data should be unmodified after init
    m = WcsNDMap(data=data, geom=geom)
    assert_equal(m.data, data)

    # Data should be unmodified if initialised like this
    m = WcsNDMap(geom=geom)
    # and then filled via an in-place Numpy array operation
    m.data += data
    assert_equal(m.data, data)
    # This is done e.g. in `m.interp_image` or probably also other operations,
    # sometimes they operate on `m.data` in-place.

    # Data should be unmodified after write / read to normal image format
    filename = str(tmpdir / 'normal.fits.gz')
    m.write(filename)
    assert_equal(Map.read(filename).data, data)

    # Data should be unmodified after write / read to sparse image format
    filename = str(tmpdir / 'sparse.fits.gz')
    m.write(filename)
    assert_equal(Map.read(filename).data, data)


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_set_get_by_pix(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, skydir=skydir,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    pix = m.geom.get_idx()
    m.set_by_pix(pix, coords[0])
    assert_allclose(coords[0], m.get_by_pix(pix))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_set_get_by_coord(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, skydir=skydir,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    m.set_by_coord(coords, coords[0])
    assert_allclose(coords[0], m.get_by_coord(coords))

    if not geom.is_allsky:
        coords[1][...] = 0.0
        assert_allclose(
            np.nan * np.ones(coords[0].shape), m.get_by_coord(coords))

    # Test with SkyCoords
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    skydir = SkyCoord(coords[0], coords[1], unit='deg',
                      frame=coordsys_to_frame(geom.coordsys))
    skydir_cel = skydir.transform_to('icrs')
    skydir_gal = skydir.transform_to('galactic')

    m.set_by_coord((skydir_gal,) + tuple(coords[2:]), coords[0])
    assert_allclose(coords[0], m.get_by_coord(coords))
    assert_allclose(m.get_by_coord((skydir_cel,) + tuple(coords[2:])),
                    m.get_by_coord((skydir_gal,) + tuple(coords[2:])))

    # Test with MapCoord
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    coords_dict = dict(lon=coords[0], lat=coords[1])
    if axes:
        for i, ax in enumerate(axes):
            coords_dict[ax.name] = coords[i + 2]
    map_coords = MapCoord.create(coords_dict, coordsys=coordsys)
    m.set_by_coord(map_coords, coords[0])
    assert_allclose(coords[0], m.get_by_coord(map_coords))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_fill_by_coord(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, skydir=skydir,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    fill_coords = tuple([np.concatenate((t, t)) for t in coords])
    fill_vals = fill_coords[1]
    m.fill_by_coord(fill_coords, fill_vals)
    assert_allclose(m.get_by_coord(coords), 2.0 * coords[1])

    # Test with SkyCoords
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    skydir = SkyCoord(coords[0], coords[1], unit='deg',
                      frame=coordsys_to_frame(geom.coordsys))
    skydir_cel = skydir.transform_to('icrs')
    skydir_gal = skydir.transform_to('galactic')
    fill_coords_cel = (skydir_cel,) + tuple(coords[2:])
    fill_coords_gal = (skydir_gal,) + tuple(coords[2:])
    m.fill_by_coord(fill_coords_cel, coords[1])
    m.fill_by_coord(fill_coords_gal, coords[1])
    assert_allclose(m.get_by_coord(coords), 2.0 * coords[1])


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_coadd(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, skydir=skydir,
                          proj=proj, coordsys=coordsys, axes=axes)
    m0 = WcsNDMap(geom)
    m1 = WcsNDMap(geom.upsample(2))
    coords = m0.geom.get_coord()
    m1.fill_by_coord(tuple([np.concatenate((t, t)) for t in coords]),
                     np.concatenate((coords[1], coords[1])))
    m0.coadd(m1)
    assert_allclose(np.nansum(m0.data), np.nansum(m1.data), rtol=1E-4)


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_interp_by_coord(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, skydir=skydir,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord(flat=True)
    m.set_by_coord(coords, coords[1])
    assert_allclose(coords[1], m.interp_by_coord(coords, interp='nearest'))
    assert_allclose(coords[1], m.interp_by_coord(coords, interp='linear'))
    assert_allclose(coords[1], m.interp_by_coord(coords, interp=1))
    if geom.is_regular and not geom.is_allsky:
        assert_allclose(coords[1], m.interp_by_coord(coords, interp='cubic'))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_iter(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    m.fill_by_coord(coords, coords[0])
    for vals, pix in m.iter_by_pix(buffersize=100):
        assert_allclose(vals, m.get_by_pix(pix))
    for vals, coords in m.iter_by_coord(buffersize=100):
        assert_allclose(vals, m.get_by_coord(coords))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_sum_over_axes(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    m.fill_by_coord(coords, coords[0])
    msum = m.sum_over_axes()

    if m.geom.is_regular:
        assert_allclose(np.nansum(m.data), np.nansum(msum.data))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_reproject(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz, proj=proj,
                          skydir=skydir, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)

    if geom.projection == 'AIT' and geom.is_allsky:
        pytest.xfail('Bug in reproject version <= 0.3.1')

    if geom.ndim > 3 or geom.npix[0].size > 1:
        pytest.xfail(
            "> 3 dimensions or multi-resolution geometries not supported")

    geom0 = WcsGeom.create(npix=npix, binsz=binsz, proj=proj,
                           skydir=skydir, coordsys=coordsys, axes=axes)
    m0 = m.reproject(geom0, order=1)

    assert_allclose(m.data, m0.data)

    # TODO : Reproject to a different spatial geometry


def test_wcsndmap_reproject_allsky_car():
    geom = WcsGeom.create(binsz=10.0, proj='CAR', coordsys='CEL')
    m = WcsNDMap(geom)
    coords = m.geom.get_coord()
    m.set_by_coord(coords, coords[0])

    geom0 = WcsGeom.create(binsz=1.0, proj='CAR', coordsys='CEL',
                           skydir=(180.0, 0.0), width=30.0)
    m0 = m.reproject(geom0, order=1)
    coords0 = m0.geom.get_coord()
    assert_allclose(m0.get_by_coord(coords0), coords0[0])

    geom1 = HpxGeom.create(binsz=5.0, coordsys='CEL')
    m1 = m.reproject(geom1, order=1)
    coords1 = m1.geom.get_coord()

    m = (coords1[0] > 10.0) & (coords1[0] < 350.0)
    assert_allclose(m1.get_by_coord((coords1[0][m], coords1[1][m])),
                    coords1[0][m])


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_pad(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    m2 = m.pad(1, mode='constant', cval=2.2)
    if not geom.is_allsky:
        coords = m2.geom.get_coord()
        msk = m2.geom.contains(coords)
        coords = tuple([c[~msk] for c in coords])
        assert_allclose(m2.get_by_coord(coords), 2.2)
    m.pad(1, mode='edge')
    m.pad(1, mode='interp')


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_crop(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    m.crop(1)


@requires_dependency('skimage')
@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_downsample(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    # Check whether we can downsample
    if (np.all(np.mod(geom.npix[0], 2) == 0) and
            np.all(np.mod(geom.npix[1], 2) == 0)):
        m2 = m.downsample(2, preserve_counts=True)
        assert_allclose(np.nansum(m.data), np.nansum(m2.data))


@pytest.mark.parametrize(('npix', 'binsz', 'coordsys', 'proj', 'skydir', 'axes'),
                         wcs_test_geoms)
def test_wcsndmap_upsample(npix, binsz, coordsys, proj, skydir, axes):
    geom = WcsGeom.create(npix=npix, binsz=binsz,
                          proj=proj, coordsys=coordsys, axes=axes)
    m = WcsNDMap(geom)
    m2 = m.upsample(2, order=0, preserve_counts=True)
    assert_allclose(np.nansum(m.data), np.nansum(m2.data))

def test_coadd_unit():
    geom = WcsGeom.create(npix=(10,10), binsz=1,
                          proj='CAR', coordsys='GAL')
    m1 = WcsNDMap(geom, data=np.ones((10,10)), unit='m2')
    m2 = WcsNDMap(geom, data=np.ones((10,10)), unit='cm2')

    m1.coadd(m2)

    assert_allclose(m1.data, 1.0001)

def test_make_region_mask():
    from regions import CircleSkyRegion
    geom = WcsGeom.create(npix=(3,3), binsz=2,
                          proj='CAR', coordsys='GAL')
    m = WcsNDMap(geom)
    region = CircleSkyRegion(SkyCoord(0, 0, unit='deg', frame='galactic'), 1.0*u.deg)
    maskmap = m.make_region_mask(region)

    assert maskmap.data.dtype == bool
    assert np.sum(maskmap.data) == 1

    maskmap = m.make_region_mask(region, inside=False)
    assert np.sum(maskmap.data) == 8


@requires_dependency('scipy')
@pytest.mark.parametrize('kernel', ['gauss', 'box', 'disk'])
def test_smooth(kernel):
    geom = WcsGeom.create(npix=(10, 10), binsz=1,
                          proj='CAR', coordsys='GAL')
    m = WcsNDMap(geom, data=np.ones((10, 10)), unit='m2')

    desired = m.data.sum()
    smoothed = m.smooth(0.2 * u.deg, kernel)
    actual = smoothed.data.sum()
    assert_allclose(actual, desired)
