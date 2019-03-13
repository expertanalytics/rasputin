
function init({vertices, normals, face_field, vertex_field}) {
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

    //geometry.computeVertexNormals();
    //geometry.normalizeNormals();

    //for ( var i = 0; i < geometry.faces.length; i ++ ) {

    //    var face = geometry.faces[ i ];
    //    face.color.setRGB( vertex_colors[3*i], vertex_colors[3*i + 1], vertex_colors[3*i] + 2);

    //}

    const phong_material = new THREE.MeshPhongMaterial( {
        specular: 0x111111,
        shininess: 250,
        //side: THREE.DoubleSide,
        vertexColors: THREE.FaceColors
    } );

    const phys_material = new THREE.MeshPhysicalMaterial( {
        metalness: 0.0,
        roughness: 0.5,
        reflectivity: 0.7,
        clearCoat: 0.0,
        //side: THREE.DoubleSide,
        vertexColors: THREE.FaceColors
    } );

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

    const state = { mesh, scene, camera, renderer, controls }
    return state
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
