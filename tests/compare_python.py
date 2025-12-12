
import sys
import os
import geopandas as gpd
import shapely.geometry
import warnings
import pandas as pd

# Add the cloned repos to the path so we can import them
sys.path.insert(0, os.path.abspath("repos/neatnet_python"))
sys.path.insert(0, os.path.abspath("repos/parenx/src"))

from neatnet import neatify
from parenx.skeletonize import skeletonize_frame

# Create a simple test network (parallel lines)
# Same coordinates as in R test
l1 = shapely.geometry.LineString([(0, 0), (100, 0)])
l2 = shapely.geometry.LineString([(0, 5), (100, 5)])

gdf = gpd.GeoDataFrame(
    {"id": [1, 2], "geometry": [l1, l2]},
    crs="EPSG:27700"
)

# Define output paths
out_neatnet = "tests/output_neatnet_python.geojson"
out_parenx = "tests/output_parenx.geojson"

print("Running neatnet (python)...")
try:
    # neatify signature: neatify(streets, ...)
    # It expects projected CRS (which we have)
    # Using defaults similar to what R might use
    # Note: neatnet needs topologically fixed input usually, but we have simple lines
    simplified_neatnet = neatify(gdf, max_segment_length=5, simplification_factor=1)
    simplified_neatnet.to_file(out_neatnet, driver="GeoJSON")
    print(f"Saved neatnet output to {out_neatnet}")
except Exception as e:
    print(f"Neatnet (python) failed: {e}")

print("Running parenx...")
try:
    # skeletonize_frame signature: skeletonize_frame(this_gs, parameter)
    # parameter is a dict with: simplify, buffer, scale, knot, segment
    # buffer=5 corresponds to 10m width approx?
    # The R test used dist=5 (radius)
    params = {
        "simplify": 0.0,
        "buffer": 5.0,
        "scale": 1.0,
        "knot": False,
        "segment": False
    }
    # It returns a GeoDataFrame (or GeoSeries?)
    simplified_parenx = skeletonize_frame(gdf.geometry, params)
    
    # Check if it's a GeoDataFrame or GeoSeries
    if isinstance(simplified_parenx, gpd.GeoSeries):
        simplified_parenx = gpd.GeoDataFrame(geometry=simplified_parenx)
        
    simplified_parenx.to_file(out_parenx, driver="GeoJSON")
    print(f"Saved parenx output to {out_parenx}")
except Exception as e:
    print(f"Parenx failed: {e}")
