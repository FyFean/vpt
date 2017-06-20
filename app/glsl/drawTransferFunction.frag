#version 300 es
precision mediump float;

#define NBUMPS @NBUMPS

struct Bump {
    vec2 position;
    float size;
    vec4 color;
};

uniform Bump uBumps[NBUMPS];

in vec2 vPosition;
out vec4 oColor;

void main() {
    vec4 color = vec4(0.0);
    for (int i = 0; i < NBUMPS; i++) {
        Bump b = uBumps[i];
        float dist = distance(b.position, vPosition);
        float x = dist / b.size;
        b.color.a *= exp(-x * x);
        b.color.rgb *= b.color.a;
        color += b.color * (1.0 - color.a);
    }
    oColor = color;
}