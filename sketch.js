/// <reference types="p5" />
// Parameters
const n = 6; // Subdivision level (increase for more dimples)

let R = 200;        // Radius of the sphere (e.g., golf ball radius)
let Rd = 10;        // Dimple radius
let dimplePositions = [];

let rotationMatrix;
let lastMouseX, lastMouseY;
let isDragging = false;
let angularVelocity = [0, 0, 0]; // Angular velocity of the sphere
let friction = 0.98; // Friction to slow down the spin
let vibrationInterval = null; // Interval for vibration

function setup() {
    createCanvas(windowWidth, windowHeight, WEBGL);
    updateDimensions();

    // Initialize rotation matrix to identity
    rotationMatrix = createIdentityMatrix();

    // Generate dimple positions
    dimplePositions = geodesicPolyhedronSubdivision(n, R);
    print(dimplePositions.length);
}

function draw() {
    background(30);

    // Apply the angular velocity to the rotation matrix
    if (!isDragging) {
        let angle = createVector(...angularVelocity).mag();
        if (angle > 0.001) { // Only rotate if there's significant angular velocity
            let axis = createVector(...angularVelocity).normalize();
            let rotMatrix = rotationMatrixFromAxisAngle(axis, angle);
            rotationMatrix = multiplyMatrices(rotMatrix, rotationMatrix);
            angularVelocity = angularVelocity.map(v => v * friction); // Apply friction

            // Trigger vibration proportional to angular velocity
            if (typeof navigator !== 'undefined' && navigator.vibrate) {
                if (vibrationInterval === null) {
                    vibrationInterval = setInterval(() => {
                        let vibrationDuration = map(angle, 0.001, 1, 100, 5); // Start with higher vibration and decrease
                        navigator.vibrate(vibrationDuration);
                    }, 100);
                }
            }
        } else {
            angularVelocity = [0, 0, 0]; // Stop rotation when velocity is minimal
            if (vibrationInterval !== null) {
                clearInterval(vibrationInterval);
                vibrationInterval = null;
            }
        }
    }

    push();

    // Apply the accumulated rotation
    applyMatrix(...rotationMatrix);

    // Define the light direction in world space (fixed relative to the viewer)
    directionalLight(253, 251, 211, -0.35, 0.25, -1);
    directionalLight(0, 0, 50, 1, -1, -1);
    directionalLight(55, 0, 0, -1, -0.75, -1);
    ambientLight(50);

    drawBall(false);

    pop();
}

function windowResized() {
    resizeCanvas(windowWidth, windowHeight);
    updateDimensions();
    dimplePositions = geodesicPolyhedronSubdivision(n, R);
}

function updateDimensions() {
    R = min(windowWidth, windowHeight) * 0.7 / 2;
    Rd = R / n / 3.5;
}

function drawBall(skeletonOnly) {
    strokeWeight(1);
    // Draw the sphere
    if (skeletonOnly) {
        stroke(255, 100);
        noFill();
    } else {
        noStroke();
        fill(220);
    }
    sphere(R);

    // Draw the dimples
    if (!skeletonOnly) {
        fill(150);
        noStroke();
        for (let i = 0; i < dimplePositions.length; i++) {
            let pos = dimplePositions[i];
            push();
            translate(pos.x, pos.y, pos.z);
            fill(lerpColor(color(255, 0, 0), color(255, 255, 0), map(pos.z, -R, R, 0, 1)));
            sphere(Rd);
            pop();
        }
    }
}

function mousePressed() {
    isDragging = true;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
}

function mouseReleased() {
    isDragging = false;
}

function multiplyMatricesWithVector(rotationMatrix, vector) {
    let result = [];
    for (let i = 0; i < 4; i++) {
        let sum = 0;
        for (let j = 0; j < 4; j++) {
            sum += rotationMatrix[i * 4 + j] * vector[j];
        }
        result.push(sum);
    }
    return result;
}

function mouseDragged() {
    if (isDragging) {
        let va = getArcballVector(lastMouseX, lastMouseY); // Vector from the center of the sphere to the last mouse position
        let vb = getArcballVector(mouseX, mouseY); // Vector from the center of the sphere to the current mouse position

        let rotateVa = multiplyMatricesWithVector(rotationMatrix, [va.x, va.y, va.z, 1]);
        let rotateVb = multiplyMatricesWithVector(rotationMatrix, [vb.x, vb.y, vb.z, 1]);
        va = createVector(rotateVa[0], rotateVa[1], rotateVa[2]);
        vb = createVector(rotateVb[0], rotateVb[1], rotateVb[2]);

        let angle = acos(min(1.0, va.dot(vb)));
        let axis = va.cross(vb);
        if (axis.magSq() > 0.0001) {
            axis.normalize();
            let rotMatrix = rotationMatrixFromAxisAngle(axis, angle);
            rotationMatrix = multiplyMatrices(rotMatrix, rotationMatrix);
            angularVelocity = [axis.x * angle, axis.y * angle, axis.z * angle]; // Set initial angular velocity based on mouse drag
        }

        lastMouseX = mouseX;
        lastMouseY = mouseY;
    }
}

function touchStarted() {
    isDragging = true;
    lastMouseX = touches[0].x;
    lastMouseY = touches[0].y;
    return false;
}

function touchEnded() {
    isDragging = false;
    return false;
}

function touchMoved() {
    if (isDragging) {
        let va = getArcballVector(lastMouseX, lastMouseY); // Vector from the center of the sphere to the last touch position
        let vb = getArcballVector(touches[0].x, touches[0].y); // Vector from the center of the sphere to the current touch position

        let rotateVa = multiplyMatricesWithVector(rotationMatrix, [va.x, va.y, va.z, 1]);
        let rotateVb = multiplyMatricesWithVector(rotationMatrix, [vb.x, vb.y, vb.z, 1]);
        va = createVector(rotateVa[0], rotateVa[1], rotateVa[2]);
        vb = createVector(rotateVb[0], rotateVb[1], rotateVb[2]);

        let angle = acos(min(1.0, va.dot(vb)));
        let axis = va.cross(vb);
        if (axis.magSq() > 0.0001) {
            axis.normalize();
            let rotMatrix = rotationMatrixFromAxisAngle(axis, angle);
            rotationMatrix = multiplyMatrices(rotMatrix, rotationMatrix);
            angularVelocity = [axis.x * angle, axis.y * angle, axis.z * angle]; // Set initial angular velocity based on touch drag
        }

        lastMouseX = touches[0].x;
        lastMouseY = touches[0].y;
    }
    return false;
}

function getArcballVector(x, y) {
    let vec = createVector();
    vec.x = -(2.0 * x - width) / width;
    vec.y = (height - 2.0 * y) / height;
    vec.z = 0.0;
    let lengthSquared = vec.x * vec.x + vec.y * vec.y;
    if (lengthSquared <= 1.0) {
        vec.z = sqrt(1.0 - lengthSquared);
    } else {
        vec.normalize();
    }
    return vec;
}

function rotationMatrixFromAxisAngle(axis, angle) {
    let c = cos(angle);
    let s = sin(angle);
    let t = 1 - c;
    let x = axis.x, y = axis.y, z = axis.z;

    return [
        t * x * x + c, t * x * y - s * z, t * x * z + s * y, 0,
        t * x * y + s * z, t * y * y + c, t * y * z - s * x, 0,
        t * x * z - s * y, t * y * z + s * x, t * z * z + c, 0,
        0, 0, 0, 1
    ];
}

function createIdentityMatrix() {
    return [
        1, 0, 0, 0, // Column 1
        0, 1, 0, 0, // Column 2
        0, 0, 1, 0, // Column 3
        0, 0, 0, 1  // Column 4
    ];
}

function multiplyMatrices(a, b) {
    let result = [];
    for (let i = 0; i < 16; i++) result[i] = 0;
    for (let row = 0; row < 4; row++) {
        for (let col = 0; col < 4; col++) {
            for (let k = 0; k < 4; k++) {
                result[row * 4 + col] += a[row * 4 + k] * b[k * 4 + col];
            }
        }
    }
    return result;
}

function geodesicPolyhedronSubdivision(n, R) {
    let vertices;
    let faces;

    // Generate the icosahedron
    [vertices, faces] = generateIcosahedron(R);

    let vertexSet = new Set();

    // Subdivide faces
    for (let face of faces) {
        let subdividedFaces = subdivideFace(vertices, face, n);
        for (let smallFace of subdividedFaces) {
            for (let P of smallFace) {
                let P_normalized = normalizeToSphere(P, R);
                // Use a string key to ensure uniqueness of vertices
                let key = `${P_normalized.x.toFixed(5)},${P_normalized.y.toFixed(5)},${P_normalized.z.toFixed(5)}`;
                vertexSet.add(key);
            }
        }
    }

    // Collect unique vertices
    let dimplePositions = [];
    vertexSet.forEach(function (value) {
        let coords = value.split(',');
        let x = parseFloat(coords[0]);
        let y = parseFloat(coords[1]);
        let z = parseFloat(coords[2]);
        dimplePositions.push(createVector(x, y, z));
    });

    return dimplePositions;
}

function generateIcosahedron(R) {
    let t = (1 + sqrt(5)) / 2;

    // Create 12 vertices of an icosahedron
    let vertices = [
        createVector(-1, t, 0),
        createVector(1, t, 0),
        createVector(-1, -t, 0),
        createVector(1, -t, 0),

        createVector(0, -1, t),
        createVector(0, 1, t),
        createVector(0, -1, -t),
        createVector(0, 1, -t),

        createVector(t, 0, -1),
        createVector(t, 0, 1),
        createVector(-t, 0, -1),
        createVector(-t, 0, 1)
    ];

    // Scale vertices to the sphere's radius
    for (let i = 0; i < vertices.length; i++) {
        vertices[i].normalize();
        vertices[i].mult(R);
    }

    // Create 20 faces of the icosahedron (indices of vertices)
    let faces = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],

        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],

        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],

        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1]
    ];

    return [vertices, faces];
}

function subdivideFace(vertices, face, n) {
    let A = vertices[face[0]].copy();
    let B = vertices[face[1]].copy();
    let C = vertices[face[2]].copy();

    let points = [];

    // Generate points on the face using barycentric coordinates
    for (let i = 0; i <= n; i++) {
        for (let j = 0; j <= n - i; j++) {
            let k = n - i - j;
            let u = i / n;
            let v = j / n;
            let w = k / n;

            // Interpolate to find the point inside the triangle
            let P = p5.Vector.mult(A, u)
                .add(p5.Vector.mult(B, v))
                .add(p5.Vector.mult(C, w));
            points.push({i: i, j: j, P: P});
        }
    }

    // Create small faces from the grid points
    let smallFaces = [];
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n - i; j++) {
            let P1 = getPoint(points, i, j);
            let P2 = getPoint(points, i + 1, j);
            let P3 = getPoint(points, i, j + 1);
            smallFaces.push([P1.P.copy(), P2.P.copy(), P3.P.copy()]);

            if (i + j < n - 1) {
                let P4 = getPoint(points, i + 1, j);
                let P5 = getPoint(points, i + 1, j + 1);

                let P6 = getPoint(points, i, j + 1);
                smallFaces.push([P4.P.copy(), P5.P.copy(), P6.P.copy()]);
            }
        }
    }
    return smallFaces;
}

function getPoint(points, i, j) {
    // Retrieve the point with indices (i, j)
    for (let pt of points) {
        if (pt.i === i && pt.j === j) {
            return pt;
        }
    }
    return null;
}

function normalizeToSphere(P, R) {
    // Normalize vector P to lie on the sphere's surface
    let P_normalized = P.copy();
    P_normalized.normalize();
    P_normalized.mult(R);
    return P_normalized;
}
