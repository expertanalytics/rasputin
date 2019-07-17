import pytest
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from PIL import Image, TiffImagePlugin
from shapely.geometry import Polygon

from rasputin.reader import read_raster_file, crop_image_to_polygon, get_image_extents
from rasputin.reader import GeoTiffTags

@contextmanager
def image_from_array(array: np.ndarray):
    # Provides a tiff image file interface to a numpy array
    image = Image.fromarray(array)
    m, n = array.shape

    # Store image extent data in tiff tags
    tags = TiffImagePlugin.ImageFileDirectory_v2()
    scale_tag = GeoTiffTags.ModelPixelScaleTag.value
    tiept_tag = GeoTiffTags.ModelTiePointTag.value

    tags[scale_tag] = (1, 1, 0)
    tags[tiept_tag] = (0, 0, 0, 0, n-1, 0)

    # Store to temporary file
    _, temp_name = tempfile.mkstemp(suffix=".tif")
    temp_path = Path(temp_name)
    try:
        fp = temp_path.open(mode="w+b")
        image.save(fp, format="tiff", tiffinfo=tags)

        yield Image.open(fp)

        fp.close()
    finally:
        # Delete temporary file
        temp_path.unlink()


@pytest.fixture
def two_level_array():
    M, N = 24, 16

    x = np.linspace(0, M-1, M)
    y = np.linspace(N-1, 0, N)

    A = np.zeros((M, N))
    for (i,j) in np.ndindex((M, N)):
        A[i,j] = 2.0 if (y[j] > max(y)/2)  else  1.0

    return A


@pytest.fixture
def random_array():
    M, N = 24, 16

    return np.random.randn(M, N)


@pytest.fixture
def upper_half_rectangle():
    return Polygon.from_bounds(0, 8, 23, 15)


def test_image_extents(two_level_array):
    with image_from_array(two_level_array) as image:
        extents = get_image_extents(image=image)
        m, n = extents.shape

        assert np.array(image).shape == (m, n)
        assert (extents.x_min, extents.y_min) == (0, 0)
        assert (extents.x_max, extents.y_max) == (m-1, n-1)


def test_image_crop(two_level_array, upper_half_rectangle):
    with image_from_array(two_level_array) as image:
        M, N = two_level_array.shape
        sub_image, extents = crop_image_to_polygon(image=image, polygon=upper_half_rectangle)
        sub_array = np.array(sub_image)

        # Check extents
        assert extents.shape == (M, np.ceil(N/2))
        assert extents.delta_x == extents.delta_y == 1
        assert extents.x_min == 0
        assert extents.x_max == M - 1
        assert extents.y_min == np.floor(N / 2)
        assert extents.y_max == N - 1

        # Check contents
        assert not (sub_array - 2.0).any()
