lake_material = """\
function create_material(geometry){
    var material = new THREE.MeshPhongMaterial( {
        specular: 0xffffff,
        shininess: 25,
        color: 0x006994,
        reflectivity: 0.3
    });
    return [new THREE.Mesh(geometry, material), null];
}
"""

avalanche_material = """\
function create_material(geometry){
    var material = new THREE.MeshPhongMaterial( {
        specular: 0xffffff,
        shininess: 75,
        reflectivity: 0.3,
        vertexColors: THREE.FaceColors
} );
    return [new THREE.Mesh(geometry, material), null];
}   
"""

terrain_material = """\
function create_material(geometry){
    var material = new THREE.MeshPhysicalMaterial({
        metalness: 0.0, 
        roughness: 0.5, 
        reflectivity: 0.5,
        clearCoat: 0.0, 
        vertexColors: THREE.FaceColors 
    });
    return [new THREE.Mesh(geometry, material), null];
}   
"""

