#version 130

uniform vec2 iResolution;
uniform float iTime;

void draw(in vec2 uv, in vec3 col, out vec3 dst);

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.yy-.5;

    // Time varying pixel color
    vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));
	draw(uv, col, col);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
