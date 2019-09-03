lake_material = """\
THREE.MeshPhongMaterial( {
    specular: 0xffffff,
    shininess: 25,
    color: 0x006994,
    reflectivity: 0.3,
} )
"""

avalanche_material = """\
THREE.MeshPhongMaterial( {
    specular: 0xffffff,
    shininess: 75,
    reflectivity: 0.3,
    vertexColors: THREE.FaceColors
} )
"""

terrain_material = """\
THREE.MeshPhysicalMaterial( {
    metalness: 0.0,
    roughness: 0.5,
    reflectivity: 0.7,
    clearCoat: 0.0,
    vertexColors: THREE.FaceColors
} )
"""

