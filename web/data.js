const indices = [
    0, 1, 2,
    3, 4, 5,
];

const vertices = new Float32Array( [
    -1.0, -1.0,  1.0,
     1.0, -1.0,  2.0,
     1.0,  1.0,  1.5,

     1.0,  1.0,  1.3,
    -1.0,  1.0,  1.7,
    -1.0, -1.0,  1.1,
] );

const colors = new Float32Array( [
     0.0,  1.0,  0.0,
     1.0,  0.0,  0.0,
     0.0,  0.0,  1.0,

     1.0,  1.0,  0.0,
     0.0,  1.0,  1.0,
     1.0,  0.0,  1.0,
] );

const data = {indices, vertices, colors};