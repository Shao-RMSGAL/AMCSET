#version 440
// Particle paths and target outlines (AMCSET trajectory views).
layout(location = 0) in vec3 position;   // Angstrom, scene axes
layout(location = 1) in vec4 color;      // straight alpha
layout(location = 2) in float generation; // -1 = outline, 0 = incident particle, >0 recoils

layout(location = 0) out vec4 v_color;

layout(std140, binding = 0) uniform buf {
  mat4 mvp;
  vec4 params; // x: incident paths shown, y: recoil paths shown, z: path opacity, w: deepest generation shown
};

void main() {
  bool outline = generation < 0.0;
  bool shown = outline || (generation < 0.5 ? params.x > 0.5
                                            : (params.y > 0.5 && generation <= params.w + 0.5));
  float a = outline ? color.a : color.a * params.z;
  v_color = vec4(color.rgb * a, a); // premultiplied
  gl_Position = shown ? mvp * vec4(position, 1.0) : vec4(2.0, 2.0, 2.0, 1.0);
}
