uniform vec3 uLight, uLight2, uColor;

varying vec3 vNormal;
varying vec3 vPosition;

void main() {
  vec3 tolight = normalize(uLight - vPosition); // direction to the point light
  vec3 tolight2 = normalize(uLight2 - vPosition);
  vec3 normal = normalize(vNormal);

  float diffuse = max(0.0, dot(normal, tolight)); // negative clamping of the cosine law
  diffuse += max(0.0, dot(normal, tolight2)); // add light energy
  vec3 intensity = uColor * diffuse; // only apply diffuse reflectance function of BRDF

  gl_FragColor = vec4(intensity, 1.0); // 4th element is opacity.
}
