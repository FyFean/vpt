// #part /glsl/shaders/renderers/DELTA/integrate/vertex

#version 300 es

const vec2 vertices[] = vec2[](
    vec2(-1, -1),
    vec2( 3, -1),
    vec2(-1,  3)
);

out vec2 vPosition;
out vec4 oDebug; // Add this line to declare oDebug


void main() {
    vec2 position = vertices[gl_VertexID];
    vPosition = position;
    gl_Position = vec4(position, 0, 1);
}

// ta dela s fotoni

// #part /glsl/shaders/renderers/DELTA/integrate/fragment

#version 300 es
precision mediump float;
precision mediump sampler2D;
precision mediump sampler3D;

#define EPS 1e-5

// #link /glsl/mixins/Photon
@Photon
// #link /glsl/mixins/intersectCube
@intersectCube

@constants
@random/hash/pcg
@random/hash/squashlinear
@random/distribution/uniformdivision
@random/distribution/square
@random/distribution/disk
@random/distribution/sphere
@random/distribution/exponential

@unprojectRand

uniform sampler2D uPosition;
uniform sampler2D uDirection;
uniform sampler2D uTransmittance;
uniform sampler2D uRadiance;

uniform sampler3D uVolume;
uniform sampler2D uTransferFunction;
uniform sampler2D uEnvironment;

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

uniform float uExtinction;
uniform float uAnisotropy;
uniform uint uMaxBounces;
uniform uint uSteps;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

void resetPhoton(inout uint state, inout Photon photon) {
    vec3 from, to;
    unprojectRand(state, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from); //v vsakmu pikslu shranjujemo ne sam barve ampak tud smer fotona, bounces, 
    photon.bounces = 0u; //kokrat se je foton odbiu znotraj volumna
    vec2 tbounds = max(intersectCube(from, photon.direction), 0.0);
    photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1); //prepustnost materiala do te tocke
}

vec4 sampleEnvironmentMap(vec3 d) {
    vec2 texCoord = vec2(atan(d.x, -d.z), asin(-d.y) * 2.0) * INVPI * 0.5 + 0.5;
    return texture(uEnvironment, texCoord);
}

// returm color and opacity of a sampled voxel based on density, uTransferFunction -> translates the raw scalar values from your volume data into visually meaningful colors and opacities for rendering.
vec4 sampleVolumeColor(vec3 position) {
    vec2 volumeSample = texture(uVolume, position).rg; // raw values, sample a 2D vector from the volume texture (uVolume) at the specified position
    vec4 transferSample = texture(uTransferFunction, volumeSample); //samples uTransferFunction, mapping between input values (often representing density or other scalar information) and output color and opacity values.
    return transferSample;
}

float sampleHenyeyGreensteinAngleCosine(inout uint state, float g) {
    float g2 = g * g;
    float c = (1.0 - g2) / (1.0 - g + 2.0 * g * random_uniform(state));
    return (1.0 + g2 - c * c) / (2.0 * g);
}

vec3 sampleHenyeyGreenstein(inout uint state, float g, vec3 direction) {
    // generate random direction and adjust it so that the angle is HG-sampled
    vec3 u = random_sphere(state);
    if (abs(g) < EPS) {
        return u;
    }
    float hgcos = sampleHenyeyGreensteinAngleCosine(state, g);
    vec3 circle = normalize(u - dot(u, direction) * direction);
    return sqrt(1.0 - hgcos * hgcos) * circle + hgcos * direction;
}

float max3(vec3 v) {
    return max(max(v.x, v.y), v.z);
}

float mean3(vec3 v) {
    return dot(v, vec3(1.0 / 3.0));
}

void main() {
    Photon photon;
    vec2 mappedPosition = vPosition * 0.5 + 0.5;
    photon.position = texture(uPosition, mappedPosition).xyz; // initial position for the photon based on the current pixel's location within the volume.
    vec4 directionAndBounces = texture(uDirection, mappedPosition); 
    photon.direction = directionAndBounces.xyz;
    photon.bounces = uint(directionAndBounces.w + 0.5);
    photon.transmittance = texture(uTransmittance, mappedPosition).rgb;
    vec4 radianceAndSamples = texture(uRadiance, mappedPosition);
    photon.radiance = radianceAndSamples.rgb;
    photon.samples = uint(radianceAndSamples.w + 0.5);
    float w = 1.0;

    uint state = hash(uvec3(floatBitsToUint(mappedPosition.x), floatBitsToUint(mappedPosition.y), floatBitsToUint(uRandSeed)));
    for (uint i = 0u; i < uSteps; i++) { // simulira rekurzijo

        float dist = random_exponential(state, uExtinction); //sampling distance to the first collision t = - ln (1-random) / uExtinction, samplamo u = [0,1], da ga plugamo v P^-1, majorant = mu_t + mu_null => constant
        photon.position += dist * photon.direction; // x = x + tw, (t-dist, w-direction)
        vec4 volumeSample = sampleVolumeColor(photon.position);

        
  
        float majorant = 1.0; // extinction = 1
        float mu_n = majorant - volumeSample.a; 
        float mu_s = volumeSample.a * max3(volumeSample.rgb);
        float mu_a = majorant - mu_n - mu_s;
        float mu_t = mu_s + mu_a;

        float Pa = mu_a / (mu_t +  max(mu_n, -mu_n)); // mu_a / (mu_t + abs(mu_n))
        float Ps = (photon.bounces >= uMaxBounces) ? 0.0 :  mu_s / (mu_t +  max(mu_n, -mu_n)); // mu_s / (mu_t + abs(mu_n))
        float Pn = max(mu_n, -mu_n) / (mu_t +  max(mu_n, -mu_n)); // abs(mu_n) / (mu_t + abs(mu_n))




        float totalP = Pa + Ps + Pn;
        // float fortuneWheel = random_uniform(state) * totalP; // [0,majorant]
        float fortuneWheel = random_uniform(state); // [0,majorant]

        // float fortuneWheel = random_uniform(state) * majorant; // [0,majorant]
        // debugTexture[gl_FragCoord.xy] = vec4(Pa, Ps, Pn, 1.0);

         vec3 debugValues = vec3(Pa, Ps, Pn);



       




        // float PNull = 1.0 - volumeSample.a; // probability of null collision is 1 - opacity
        // float PScattering;
        // if (photon.bounces >= uMaxBounces) {
        //     PScattering = 0.0;   //ce je ze max bounces potem je scattering = 0
        // } else {
        //     PScattering = volumeSample.a * max3(volumeSample.rgb);  // based on the volumes opacity and max color channel value 
        // }
        // float PAbsorption = 1.0 - PNull - PScattering; // Pa + Ps + Pn = 1
        // float fortuneWheel = random_uniform(state); // [0,1] need to fix this line

        if (any(greaterThan(photon.position, vec3(1))) || any(lessThan(photon.position, vec3(0)))) { // we hit a boundary
            // out of bounds
            vec4 envSample = sampleEnvironmentMap(photon.direction);  //kao zadi za volumom (recimo sky)
            vec3 radiance = photon.transmittance * envSample.rgb; //sevalnost L = prosojnost * environment * w, L je trenutna intenziteta
            photon.samples++;
            photon.radiance += (radiance - photon.radiance) / float(photon.samples); //incremental averaging, curr + (new - curr) / sample count 
            resetPhoton(state, photon);
            w = 1.0;
        } else if (fortuneWheel < Pa) { 
            // absorption -> photon radiance = 0, reset the photon
            vec3 radiance = vec3(0); //emisije ni, Le
            photon.samples++;
            photon.radiance += (radiance - photon.radiance) / float(photon.samples);
            resetPhoton(state, photon);
            // Apply weight for absorption, Le = 0, w = mu_a * Le/... -> no weight bc no absorbtion
            w = 1.0;
        } else if (fortuneWheel < Pa + Ps) {
            // scattering -> update transmittance and direction of a photon, mu_s * Ls, Ls je izracunan kot integral phase functiona * L
            photon.transmittance *= volumeSample.rgb;
            photon.direction = sampleHenyeyGreenstein(state, uAnisotropy, photon.direction); //phase function km se odbija
            photon.bounces++;
            float weightS = mu_s / (majorant * Ps);
            // float weightS = Ps * (mu_s / 1.0);

            w *= weightS; // Apply weight for scattering
        } else {
            // null collision, Ln je ubistvu L, seva samo v eno smer
            // float weightN = max(mu_n, -mu_n) / (majorant * Pn);
            float weightN = mu_n / (majorant * Pn);

            // float weightN = Pn * (mu_n / 1.0);

            w *= weightN;
        }

        // oDebug = vec4(Pa, Ps, Pn, 1.0); // Output debug values
    }
    photon.radiance *= w;
    // photon.radiance *= 0.0001;


    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}

// #part /glsl/shaders/renderers/DELTA/render/vertex

#version 300 es

const vec2 vertices[] = vec2[](
    vec2(-1, -1),
    vec2( 3, -1),
    vec2(-1,  3)
);

out vec2 vPosition;

void main() {
    vec2 position = vertices[gl_VertexID];
    vPosition = position * 0.5 + 0.5;
    gl_Position = vec4(position, 0, 1);
}

// #part /glsl/shaders/renderers/DELTA/render/fragment

#version 300 es
precision mediump float;
precision mediump sampler2D;

uniform sampler2D uColor;

in vec2 vPosition;

out vec4 oColor;

void main() {
    oColor = vec4(texture(uColor, vPosition).rgb, 1);
}

// #part /glsl/shaders/renderers/DELTA/reset/vertex

#version 300 es

const vec2 vertices[] = vec2[](
    vec2(-1, -1),
    vec2( 3, -1),
    vec2(-1,  3)
);

out vec2 vPosition;

void main() {
    vec2 position = vertices[gl_VertexID];
    vPosition = position;
    gl_Position = vec4(position, 0, 1);
}

// #part /glsl/shaders/renderers/DELTA/reset/fragment

#version 300 es
precision mediump float;

// #link /glsl/mixins/Photon
@Photon
// #link /glsl/mixins/intersectCube
@intersectCube

@constants
@random/hash/pcg
@random/hash/squashlinear
@random/distribution/uniformdivision
@random/distribution/square
@random/distribution/disk
@random/distribution/sphere
@random/distribution/exponential

@unprojectRand

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

void main() {
    Photon photon;
    vec3 from, to;
    uint state = hash(uvec3(floatBitsToUint(vPosition.x), floatBitsToUint(vPosition.y), floatBitsToUint(uRandSeed)));
    unprojectRand(state, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    vec2 tbounds = max(intersectCube(from, photon.direction), 0.0);
    photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
    photon.radiance = vec3(1);
    photon.bounces = 0u;
    photon.samples = 0u;
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}


// // returm color and opacity of a sampled voxel based on density, uTransferFunction -> translates the raw scalar values from your volume data into visually meaningful colors and opacities for rendering.
// vec4 sampleVolumeColor(vec3 position) {
//     vec2 volumeSample = texture(uVolume, position).rg; // raw values, sample a 2D vector from the volume texture (uVolume) at the specified position
//     vec4 transferSample = texture(uTransferFunction, volumeSample); //samples uTransferFunction, mapping between input values (often representing density or other scalar information) and output color and opacity values.
//     return transferSample;
// }


// void main() { //gremo cez volumen, kjer se zark seka skozi kocko, izracunas fiksno kaksna je dolzina koraka, prinesemo bias s tem notr -> to das lhko v monte carlo obliko kjer zgeneriras sample in integriras spodi v to kar je ze prej bilo na sliki
//     vec3 rayDirection = vRayTo - vRayFrom; //  Calculates the direction of the ray by subtracting the origin (vRayFrom) from the target (vRayTo).
//     vec2 tbounds = max(intersectCube(vRayFrom, rayDirection), 0.0); // calculates the intersection points of the ray with the cube. 
//     // It returns a vector vec2 where x represents the distance along the ray where it enters the cube, and y represents the distance where it exits


//     if (tbounds.x >= tbounds.y) { //no intersection
//         oColor = vec4(0, 0, 0, 1);
//     } else {        
//         const float epsilon = 0.001;
//         // 1. samplamo u = [0,1], da ga plugamo v P^-1, majorant = mu_t + mu_null => constant
//         const u = Math.random();
//         const majorant = uExtinction + 10 // should it be zgornji bound moje kocke?
//         // 2. sampling distance to the first collision t, P^-1(u) = - ln(1 - u) / mu
//         const t = - Math.log(1 - u) / majorant;
//         // 3. current ray position
//         vec3 position = mix(from, to, t); 
//         // 4. random sample a je real ali null collision, scale [0,majorant]
//         const random_sample_majorant = Math.random() * majorant; 
//         const mu_t = texture(uVolume, position).r; // a je okay ce sam red channel vzamem? should i scale it
//         const p_real = mu_t / majorant;
//         if (random_sample_majorant < p_real){

//         }
        
        







//     // if (tbounds.x >= tbounds.y) { //no intersection
//     //     oColor = vec4(0, 0, 0, 1);
//     // } else {         
//     //     vec3 from = mix(vRayFrom, vRayTo, tbounds.x);
//     //     vec3 to = mix(vRayFrom, vRayTo, tbounds.y);
//     //     float rayStepLength = distance(from, to) * uStepSize;

//     //     float t = uStepSize * uOffset;
//     //     vec4 accumulator = vec4(0); // to store accumulated color

//     //     while (t < 1.0 && accumulator.a < 0.99) {
//     //         vec3 position = mix(from, to, t); // current position
//     //         vec4 colorSample = sampleVolumeColor(position); //samples the current color of the volume
//     //         colorSample.a *= rayStepLength * uExtinction; //uExtinction = user defined mu coeff kok absorbira, proportional to step
//     //         colorSample.rgb *= colorSample.a; //color sample pomnozis z alfo to account for transparency
//     //         accumulator += (1.0 - accumulator.a) * colorSample; // 1-transparency das v accumulated value k ga nakonc returnas
//     //         t += uStepSize;
//     //     }

//     //     if (accumulator.a > 1.0) {
//     //         accumulator.rgb /= accumulator.a; //overflow of alpha ce je vec kot 1
//     //     }

//     //     oColor = vec4(accumulator.rgb, 1); // alpha to 1

//     }
// }
