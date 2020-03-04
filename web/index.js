import {GUI} from 'https://threejsfundamentals.org/threejs/../3rdparty/dat.gui.module.js';

function init({geometries}) {

    const wrapModes = {
    'RepeatWrapping': THREE.RepeatWrapping,
    'ClampToEdgeWrapping': THREE.ClampToEdgeWrapping,
    'MirroredRepeatWrapping': THREE.MirroredRepeatWrapping,
    }
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

    //const textureLoader = new THREE.TextureLoader();
    //var texture = textureLoader.load("textures/sea.jpg");
    //texture.wrapS = THREE.RepeatWrapping;
    //texture.wrapT = THREE.RepeatWrapping;

    const meshes = [];
    for ( var i = 0; i < geometries.length; i ++ ) {
        const geom = new THREE.BufferGeometry();
        geom.addAttribute( 'position', new THREE.Float32BufferAttribute( geometries[i].vertices, 3) );
        geom.addAttribute( 'normal'  , new THREE.Float32BufferAttribute( geometries[i].normals, 3 ) );
        geom.addAttribute( 'color'   , new THREE.Float32BufferAttribute( geometries[i].colors, 3 ) );
        geom.addAttribute( 'uv'      , new THREE.Float32BufferAttribute( geometries[i].uvs, 2));
        geom.computeVertexNormals();
        const result = geometries[i].material_constructor(geom);
        scene.add(result[0]);
        meshes.push(result[0]);
    }

    //const gui = new GUI({width: 280});
    //gui.add(new StringToNumberHelper(texture, 'wrapS'), 'value', wrapModes).name('Wrap Horizontal').onChange(updateTexture);
    //gui.add(new StringToNumberHelper(texture, 'wrapT'), 'value', wrapModes).name('Wrap Vertical').onChange(updateTexture);
    //gui.add(texture.repeat, 'x', 0, 5, .01).name('Repeat X');
    //gui.add(texture.repeat, 'y', 0, 5, .01).name('Repeat Y');
    //gui.add(texture.offset, 'x', -2, 2, .01).name('Offset X');
    //gui.add(texture.offset, 'y', -2, 2, .01).name('Offset Y');
    //gui.add(texture.center, 'x', -.5, 1.5, .01).name('Center X');
    //gui.add(texture.center, 'y', -.5, 1.5, .01).name('Center Y');

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

    const state = { meshes, scene, camera, renderer, controls };
    return state
}

function create_material(geom, texture) {
    var geomMaterial = new THREE.MeshPhysicalMaterial({map: texture, bumpMap: texture, vertexColors: THREE.FaceColors});
    return new THREE.Mesh(geom, geomMaterial);
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
    var z_min = 10000;
    var z_max = -10000;
    for (var i = 0; i < meshes.length; i ++ ) {
        var geometry = meshes[i].geometry;
        geometry.computeBoundingBox();
        if (geometry.boundingBox.min.z < z_min) {
            z_min = geometry.boundingBox.min.z;
        }
        if (geometry.boundingBox.max.z > z_max) {
            z_max = geometry.boundingBox.max.z;
        }
    }
    return {z_min, z_max};
}

function updateTexture() {
    texture.needsUpdate = true;
}

class StringToNumberHelper {
    constructor(obj, prop) {
        this.obj = obj;
        this.prop = prop;
    }
    get value() {
        return this.obj[this.prop];
    }
    set value(v) {
        this.obj[this.prop] = parseFloat(v);
    }
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
