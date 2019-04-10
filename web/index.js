
function init({vertices, normals, face_field, vertex_field, features}) {

    THREE.ImageUtils.crossOrigin = "";
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

    const geometry = new THREE.BufferGeometry();
    //geometry.setIndex(indices);
    geometry.addAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );
    geometry.addAttribute( 'normal', new THREE.Float32BufferAttribute( normals, 3 ) );
    geometry.addAttribute( 'color', new THREE.Float32BufferAttribute( face_field, 3 ) );

    //var sea_texture = new THREE.TextureLoader().load("textures/sea.jpg");
    //sea_texture.encoding = THREE.sRGBEncoding;
    //sea_texture.anisotropy = 16;
    //sea_texture.magFilter = THREE.LinearFilter;
    //sea_texture.minFilter = THREE.LinearFilter;
    //sea_texture.mapping = THREE.CubeReflectionMapping;
    //sea_texture.repeat.set(4096, 4096);
    //sea_texture.wrapS = THREE.RepeatWrapping;
    //sea_texture.wrapT = THREE.RepeatWrapping;

    const feature_material = new THREE.MeshPhongMaterial( {
    //const feature_material = new THREE.MeshBasicMaterial( {
        //map: sea_texture,
        //bumpScale: 50,
        specular: 0xffffff,
        shininess: 25,
        //side: THREE.DoubleSide,
        color: 0x006994,
        reflectivity: 0.3,
        //vertexColors: THREE.FaceColors
    } );

    const phys_material = new THREE.MeshPhysicalMaterial( {
        metalness: 0.0,
        roughness: 0.5,
        reflectivity: 0.7,
        clearCoat: 0.0,
        //side: THREE.DoubleSide,
        vertexColors: THREE.FaceColors
    } );

    fs = [];

    for ( var i = 0; i < features.length; i ++ ) {
        const geom = new THREE.BufferGeometry();
        geom.addAttribute( 'position', new THREE.Float32BufferAttribute(features[i], 3));
        geom.computeVertexNormals();
        var uv_geom = assign_uvs(geom);
        var f = new THREE.Mesh(uv_geom, feature_material);
        scene.add(f);
        fs.push(f);
    }

    const mesh = new THREE.Mesh( geometry, phys_material );
    scene.add(mesh);

    const center = getCenterPoint(mesh);

    var z_max = mesh.geometry.boundingBox.max.z;
    var z_min = mesh.geometry.boundingBox.min.z;
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

    const state = { mesh, scene, camera, renderer, controls, fs};
    return state
}

function assign_uvs(geometry) {
    geometry.computeBoundingBox();

    var max = geometry.boundingBox.max,
        min = geometry.boundingBox.min;
    var offset = new THREE.Vector2(0 - min.x, 0 - min.y);
    var range = new THREE.Vector2(max.x - min.x, max.y - min.y);
    //var vertices = geometry.attributes.position.array;
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

function getCenterPoint(mesh) {
    var geometry = mesh.geometry;
    geometry.computeBoundingBox();   
    const center = new THREE.Vector3()
    geometry.boundingBox.getCenter(center);
    mesh.localToWorld(center);
    return center;
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
