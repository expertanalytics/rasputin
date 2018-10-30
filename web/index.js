
function init({indices, vertices, colors}) {
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
    geometry.setIndex(indices);
    geometry.addAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );
    // geometry.addAttribute( 'normal', new THREE.Float32BufferAttribute( normals, 3 ) );
    geometry.addAttribute( 'color', new THREE.Float32BufferAttribute( colors, 3 ) );

    const material = new THREE.MeshPhongMaterial( {
        specular: 0x111111,
        shininess: 250,
        side: THREE.DoubleSide,
        vertexColors: THREE.VertexColors,
    } );

    const mesh = new THREE.Mesh( geometry, material );
    scene.add( mesh );

    const camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 1000 );
    camera.position.z = 5;

    const renderer = new THREE.WebGLRenderer( { antialias: true } );
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize( window.innerWidth, window.innerHeight );

    const controls = new THREE.OrbitControls( camera, renderer.domElement );
    controls.enableDamping = true; // an animation loop is required when either damping or auto-rotation are enabled
    controls.dampingFactor = 0.25;
    controls.screenSpacePanning = false;
    controls.minDistance = 1;
    controls.maxDistance = 5;
    controls.maxPolarAngle = Math.PI / 2;

    const state = {mesh, scene, camera, renderer, controls}
    return state
}

function animate() {
    const {renderer, camera, scene, controls} = state;
    requestAnimationFrame( animate );

    controls.update();

	renderer.render( scene, camera );
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
