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


def get_forecasts(*, x: float, y: float, proj: str) -> Dict[str, Any]:
    p_source = pyproj.Proj(proj)
    p_target = pyproj.Proj("+init=epsg:4326")
    lat, lon = pyproj.transform(p_source, p_target, x, y)
    lang = 1
    now = datetime.datetime.now()
    start = now.strftime("%Y-%m-%d")
    end = (now + datetime.timedelta(days=1)).strftime("%Y-%m-%d")

    query = f"/AvalancheWarningByCoordinates/Detail/{lon}/{lat}/{lang}/{start}/{end}"
    full = base_avalanche_url + query
    response = requests.get(full)
    if response.ok:
        return json.loads(response.text)
    return {}
