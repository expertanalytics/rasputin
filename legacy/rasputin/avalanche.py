from typing import List, Dict, Any
import requests
import json
import pyproj
import datetime

base_avalanche_url = "https://api01.nve.no/hydrology/forecast/avalanche/v4.0.2/api"

varsom_angles = [( -22.5,  22.5),
                 (  22.5,  67.5),
                 (  67.5, 112.5),
                 ( 112.5, 157.5),
                 ( 157.5,-157.5),
                 (-157.5,-112.5),
                 (-112.5, -67.5),
                 ( -67.5, -22.5)]


def get_forecasts(*, x: float, y: float, crs: pyproj.CRS) -> Dict[str, Any]:
    target_crs = pyproj.CRS.from_epsg(4326)
    proj = pyproj.Transformer.from_crs(crs, target_crs)
    lat, lon = proj.transform(x, y)
    lang = 1
    now = datetime.datetime.now()
    start = now.strftime("%Y-%m-%d")
    end = (now + datetime.timedelta(days=1)).strftime("%Y-%m-%d")

    query = f"/AvalancheWarningByCoordinates/Detail/{lat}/{lon}/{lang}/{start}/{end}"
    full = base_avalanche_url + query
    response = requests.get(full)
    if response.ok:
        res = json.loads(response.text)
        return res
    return {}
