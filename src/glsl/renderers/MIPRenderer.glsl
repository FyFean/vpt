// #part /glsl/shaders/renderers/MIP/generate/vertex

#version 300 es

uniform mat4 uMvpInverseMatrix;

out vec3 vRayFrom;
out vec3 vRayTo;

// vsi shaderji k so za MIP renderer, za vsak shader je vertex in fragment shader
// shaderji se zapakirajo nakonc v en json format, v shaders.json


// ta @unproject je mixin kar je neka funckija k jo lhko klices skoz in se pr buildu notr poturi, unproject je ime mixina
// v webgl.js je buildPrograms funkcija ki najde to afno in jo zamena z mixinom k builda, src/glsl/mixins

// #link /glsl/mixins/unproject
@unproject

const vec2 vertices[] = vec2[](
    vec2(-1, -1),
    vec2( 3, -1),
    vec2(-1,  3)
);

void main() {
    vec2 position = vertices[gl_VertexID];
    unproject(position, uMvpInverseMatrix, vRayFrom, vRayTo); 
    gl_Position = vec4(position, 0, 1);
}

// #part /glsl/shaders/renderers/MIP/generate/fragment

#version 300 es
precision mediump float;
precision mediump sampler2D;
precision mediump sampler3D;

uniform sampler3D uVolume;
uniform sampler2D uTransferFunction;
uniform float uStepSize;
uniform float uOffset;

in vec3 vRayFrom;
in vec3 vRayTo;

out float oColor;

// #link /glsl/mixins/intersectCube
@intersectCube

vec4 sampleVolumeColor(vec3 position) {
    vec2 volumeSample = texture(uVolume, position).rg;
    vec4 transferSample = texture(uTransferFunction, volumeSample);
    return transferSample;
}

void main() {
    vec3 rayDirection = vRayTo - vRayFrom; //generiramo zarek v unproject zgoraj
    vec2 tbounds = max(intersectCube(vRayFrom, rayDirection), 0.0); // kje se zarek seka s kocko
    if (tbounds.x >= tbounds.y) {
        oColor = 0.0; // ce se ne seka je output 0, ce se samplamo ta volumn
    } else {
        vec3 from = mix(vRayFrom, vRayTo, tbounds.x);
        vec3 to = mix(vRayFrom, vRayTo, tbounds.y);

        float t = 0.0;
        float val = 0.0;
        float offset = uOffset;
        vec3 pos;
        do {
            pos = mix(from, to, offset);
            val = max(sampleVolumeColor(pos).a, val); //samplamo ta max in ga v korakih posodabljamo
            t += uStepSize;
            offset = mod(offset + uStepSize, 1.0);
        } while (t < 1.0);
        oColor = val;
    }
}

// #part /glsl/shaders/renderers/MIP/integrate/vertex

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

// #part /glsl/shaders/renderers/MIP/integrate/fragment

#version 300 es
precision mediump float;
precision mediump sampler2D;

uniform sampler2D uAccumulator;
uniform sampler2D uFrame;

in vec2 vPosition;

out float oColor;

void main() {
    float acc = texture(uAccumulator, vPosition).r; // uAccumulator je to kar je blo ze prej na sliki
    float frame = texture(uFrame, vPosition).r; // uFrame je kar smo lihkar zgeneriral, to je max sample k smo ga zgori izracunal
    oColor = max(acc, frame); // vzame max teh dveh, nato gre to v render korak
}

// #part /glsl/shaders/renderers/MIP/render/vertex

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

// #part /glsl/shaders/renderers/MIP/render/fragment -> tuki dejansko izrisemo na piksl rgb barvo

#version 300 es
precision mediump float;
precision mediump sampler2D;

uniform sampler2D uAccumulator;

in vec2 vPosition;

out vec4 oColor;

void main() {
    float acc = texture(uAccumulator, vPosition).r;
    oColor = vec4(acc, acc, acc, 1); // dobimo akumulirano vrednosti in jo zapisemo na 3 rbg kanale, lhko to tud poljubno spremenimo
}

// #part /glsl/shaders/renderers/MIP/reset/vertex -> ko se pogled kamere spremeni se poenostavijo vrendnosti 

#version 300 es

const vec2 vertices[] = vec2[](
    vec2(-1, -1),
    vec2( 3, -1),
    vec2(-1,  3)
);

void main() {
    vec2 position = vertices[gl_VertexID];
    gl_Position = vec4(position, 0, 1);
}

// #part /glsl/shaders/renderers/MIP/reset/fragment

#version 300 es
precision mediump float;

out float oColor;

void main() {
    oColor = 0.0;
}
