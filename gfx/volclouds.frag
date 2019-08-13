#version 130

uniform float iTime;
uniform vec2 iResolution;
uniform float iFader0;
uniform float iFader1;
uniform float iFader2;
uniform float iFader3;
uniform float iFader4;
uniform float iFader5;
uniform float iFader6;
uniform float iFader7;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);
float a = 1.0;

float nbeats, iScale;

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW
void hash13(in vec3 p3, out float d)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    d = fract((p3.x + p3.y) * p3.z);
}

// Arbitrary-frequency 2D noise
void lfnoise3(in vec3 t, out float num)
{
    t -= vec3(11.,13.,5.);
    vec3 i = floor(t);
    t = fract(t);
    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise
    t = smoothstep(c.yyy, c.xxx, t); // TODO: add this for faster value noise
    vec2 v1, v2, v3, v4;
    hash13(i, v1.x);
    hash13(i+c.xyy, v1.y);
    hash13(i+c.yxy, v2.x);
    hash13(i+c.xxy, v2.y);
    hash13(i+c.yyx, v3.x);
    hash13(i+c.xyx, v3.y);
    hash13(i+c.yxx, v4.x);
    hash13(i+c.xxx, v4.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    v3 = c.zz+2.*mix(v3, v4, t.y);
    v2.x = mix(v1.x, v1.y, t.x);
    v2.y = mix(v3.x, v3.y, t.x);
    num = mix(v2.x, v2.y, t.z);
}

void mfnoise3(in vec3 x, in float d, in float b, in float e, out float n)
{
    n = 0.;
    float a = 1., nf = 0., buf;
    for(float f = d; f<b; f *= 2.)
    {
        lfnoise3(f*x-vec3(11.,13.,5.), buf);
        n += a*buf;
        a *= e;
        nf += 1.;
    }
    n *= (1.-e)/(1.-pow(e, nf));
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

mat3 R;
vec3 ind;
void scene(in vec3 x, out vec2 sdf)
{
    //x = R * x;
    
    float n;
    mfnoise3(x-.1*iTime*c.yyx,10.,400.,.25,n);
    n = .5+.5*n;
    sdf = vec2(-n, 2.);
}

void normal(in vec3 x, out vec3 n, in float dx)
{
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

void palette1(in float scale, out vec3 col)
{
    const int N = 5;
   
    /*
    const vec3 colors[N] = vec3[N](
            vec3(0.82,0.27,0.13),
            vec3(0.85,0.77,0.68),
            vec3(0.65,0.59,0.55),
            vec3(0.45,0.29,0.24),
            vec3(0.85,0.27,0.15)
        );
    //*/
    
    /*
	const vec3 colors[N] = vec3[N](
       	vec3(0.86,0.21,0.13),
        vec3(0.85,0.80,0.62),
        vec3(0.22,0.25,0.25),
        vec3(0.16,0.17,0.17),
        vec3(0.12,0.12,0.13)
    );
    //*/
    
	//*
    const vec3 colors[N] = vec3[N](
       	vec3(0.99,0.33,0.05),
        vec3(0.94,0.94,0.94),
        vec3(0.75,0.82,0.88),
        vec3(0.25,0.34,0.39),
        vec3(0.17,0.22,0.27)
    );
    //*/
	float index = floor(scale*float(N)), 
        remainder = scale*float(N)-index;
    col = mix(colors[int(index)],colors[int(index)+1], remainder);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    rot3(.3*vec3(1.1,1.3,1.5)*iTime, R);
    
    float a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;
    
    float d = 0.;
    vec2 s;
    vec3 o, t, dir, x, n;
    
	o = c.yyx;
    t = c.yyy;
    int N = 50,
        i;
    dir = normalize(vec3(uv,-1.));//normalize(t-o);
    
    for(i = 0; i<N; ++i)
    {
        d += .5/float(N);
        x = o + d * dir;
        scene(x,s);
        normal(x,n,5.e-4);
        vec3 l = normalize(x+.1*n);
        vec3 c1;
        palette1(-s.x, c1);
        c1 = .1*c1
                            + .1*c1 * abs(dot(l,n))
                            + mix(1.,10.,iFader0) * c1 * abs(pow(dot(reflect(-l,n),dir),2.));
    	col = mix(col, c1, d*d);
    }

    col *= col;
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
