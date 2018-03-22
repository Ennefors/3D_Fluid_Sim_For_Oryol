@vs vs
uniform params {
    vec3 pos;
    mat4 model;
    mat4 view;
    mat4 projection;
};

in vec4 position;
in vec4 color0;
out vec4 color;

void main() {
    gl_Position = projection * view * model * position;
    color = color0;
}
@end

@fs fs 
in vec4 color;
out vec4 fragColor;
void main() {
    fragColor = color;
}
@end

@program Shader vs fs
