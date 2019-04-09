
function init({geometries}) {

    const scene = new THREE.Scene();
    scene.background = new THREE.Color( 0xdddddd );

    const ambientLight = new THREE.AmbientLight( 0x222222 );
    scene.add( ambientLight );

    const light1 = new THREE.DirectionalLight( 0xffffff, 0.5 );
    light1.position.set( 1, 1, 1 );
    scene.add( light1 );

    const light2 = new THREE.DirectionalLight( 0xffffff, 1 );
    light2.position.set( 0, - 1, 0 );
    scene.add( light2 );

    const meshes = [];
    for ( var i = 0; i < geometries.length; i ++ ) {
        const geom = new THREE.BufferGeometry();
        geom.addAttribute( 'position', new THREE.Float32BufferAttribute( geometries[i].vertices, 3) );
        geom.addAttribute( 'normal', new THREE.Float32BufferAttribute( geometries[i].normals, 3 ) );
        geom.addAttribute( 'color', new THREE.Float32BufferAttribute( geometries[i].colors, 3 ) );
        geom.computeVertexNormals();
        var mesh = new THREE.Mesh(geom, geometries[i].material);
        scene.add(mesh);
        meshes.push(mesh);
    }
    const test_geom = new THREE.BufferGeometry();
    test_geom.addAttribute( 'position', new THREE.Float32BufferAttribute( geometries[2].vertices, 3) );
    test_geom.addAttribute( 'normal', new THREE.Float32BufferAttribute( geometries[2].normals, 3 ) );
    test_geom.addAttribute( 'color', new THREE.Float32BufferAttribute( geometries[2].colors, 3 ) );

    var test_mesh = new THREE.Mesh(test_geom, geometries[2].material);
    scene.add(test_mesh);

    const center = getCenterPoint(meshes);
    const {z_min, z_max} = getMinMaxZ(meshes);
    var cam_z_pos = z_min + (z_max - z_min)*3;
    const camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, cam_z_pos*4);
    camera.position.set(center.x, center.y, cam_z_pos);
    camera.up.set(0, 0, 1);

    const renderer = new THREE.WebGLRenderer( { antialias: true } );
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize( window.innerWidth, window.innerHeight );

    const controls = new THREE.OrbitControls( camera, renderer.domElement );
    controls.enableDamping = true; // an animation loop is required when either damping or auto-rotation are enabled
    controls.dampingFactor = 0.25;
    controls.target.set(center.x, center.y, z_min);

    const state = { meshes, scene, camera, renderer, controls, test_mesh};
    return state
}

function assign_uvs(geometry) {
    geometry.computeBoundingBox();

    var max = geometry.boundingBox.max,
        min = geometry.boundingBox.min;
    var offset = new THREE.Vector2(0 - min.x, 0 - min.y);
    var range = new THREE.Vector2(max.x - min.x, max.y - min.y);
    var vertices = geometry.getAttribute("position").array;
    console.log(vertices.length);

    var uvs = new Float32Array(2*vertices.length);

    for (var i = 0; i < vertices.count/3 ; i++) {

        var v1 = new THREE.Vector2(vertices[3*i], vertices[3*i + 1]);
        var v2 = new THREE.Vector2(vertices[3*i + 3], vertices[3*i + 4]);
        var v3 = new THREE.Vector2(vertices[3*i + 6], vertices[3*i + 7]);

        uvs[3*i    ] = (v1.x + offset.x)/range.x;
        uvs[3*i + 1] = (v1.y + offset.y)/range.y;
        uvs[3*i + 2] = (v2.x + offset.x)/range.x;
        uvs[3*i + 3] = (v2.y + offset.y)/range.y;
        uvs[3*i + 4] = (v3.x + offset.x)/range.x;
        uvs[3*i + 5] = (v3.y + offset.y)/range.y;
    }
    geometry.addAttribute("uv", new THREE.BufferAttribute(uvs, 2));
    geometry.attributes.uv.needsUpdate = true;
    return geometry;
}

function animate() {
    const {renderer, camera, scene, controls} = state;
    requestAnimationFrame( animate );

    controls.update();
	renderer.render( scene, camera );
}

function getCenterPoint(meshes) {
    var max_area = 0;
    var center = new THREE.Vector3();
    for (var i = 0; i < meshes.length; i ++ ) {
        var geometry = meshes[i].geometry;
        geometry.computeBoundingBox();
        var size = new THREE.Vector3();
        geometry.boundingBox.getSize(size);
        var area = size.x*size.y;
        if (area > max_area) {
            max_area = area;
            geometry.boundingBox.getCenter(center);
            meshes[i].localToWorld(center);
        }
    }
    return center;
}

function getMinMaxZ(meshes) {
    var min_z = 10000;
    var max_z = -10000;
    for (var i = 0; i < meshes.length; i ++ ) {
        var geometry = meshes[i].geometry;
        geometry.computeBoundingBox();
        if (geometry.boundingBox.min.z < min_z) {
            min_z = geometry.boundingBox.min.z;
        }
        if (geometry.boundingBox.max.z > max_z) {
            max_z = geometry.boundingBox.max.z;
        }
    }
    return  {min_z, max_z};
}


function onWindowResize({renderer, camera}) {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
}

const state = init(data);
document.getElementById("thecanvas").appendChild( state.renderer.domElement );
window.addEventListener( 'resize', () => onWindowResize(state), false );
animate();
