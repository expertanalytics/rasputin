""" Specifications of terrains."""

lake_material = {
    "material_type": "THREE.MeshPhongMaterial",
    "material_params": {
        "specular": "0xffffff",
        "shininess": 25,
        "color": "0x006994",
        "reflectivity": 0.3
        }
    }


avalanche_material = {
    "material_type": "THREE.MeshPhongMaterial",
    "material_params": {
        "specular": "0xffffff",
        "shininess": 25,
        "vertexColors": "THREE.FaceColors",
        "reflectivity": 0.3
        }
    }


terrain_material = {
    "material_type": "THREE.MeshPhysicalMaterial",
    "material_params": {
        "metalness": 0.0, 
        "roughness": 0.5, 
        "reflectivity": 0.5,
        "clearcoat": 0.0, 
        "vertexColors": "THREE.FaceColors",
        }
    }
