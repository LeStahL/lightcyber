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

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW
void hash13(in vec3 p3, out float d)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    d = fract((p3.x + p3.y) * p3.z);
}

// Add sdfs
void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = sda.x<sdb.x?sda:sdb;
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void dvoronoi3(in vec3 x, out float d, out vec3 z)
{
    vec3 y = floor(x);
    float ret = 1.;
    vec3 pf=c.yyy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                float pa;
                hash13(p, pa);
                p += pa;

                d = length(x-p);

                if(d < df)
                {
                    df = d;
                    pf = p;
                }
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                float pa;
                hash13(p, pa);
                p += pa;

                vec3 o = p - pf;
                d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
                ret = min(ret, d);
            }
        }
    
    d = ret;
    z = pf;
}

vec3 ind;
void scene(in vec3 x, out vec2 sdf)
{
    mat3 R;
    rot3(.1*vec3(1.1,1.3,1.5)*iTime, R);
    x = R * x;
    x.z -= .1*iTime;
    
    sdf = c.xy;
    
    float bsize = 2.;
    float d, da;
    
    dvoronoi3(bsize*x, d, ind);
    vec3 y = x-ind/bsize;
    
    float n;
    hash13(ind,n);
    
    //add(sdf, vec2(abs(length(y)-.3)-.001,2.), sdf);
	add(sdf, vec2(abs(length(y)-.3)-.001,2.), sdf);
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    
    vec3 col = c.yyy;
    
    float d = 0.;
    vec2 s;
    vec3 o, t, dir, x, n;
    
    //mat3 Ra;
    //rot3(mix(c.yyy,vec3(-3.*pi/4.,3.*pi/4.,-7.*pi/4.),clamp((iTime-6.)/1.5,0.,1.)), Ra);
    //Ra *= mix(1.,-1.,clamp((iTime-6.)/1.5,0.,1.));
       
    
    //o = Ra * mix(mix(mix(c.yyy-.1*c.yxy,c.yyx,clamp(iTime/2.,0.,1.)),10.*c.yyx,clamp((iTime-2.)/2.,0.,1.)), 100.*c.yyx, clamp((iTime-4.)/2.,0.,1.));
	o = c.yyx;
    t = c.yyy;
    int N = 250,
        i;
    dir = normalize(vec3(uv,-1.));//normalize(t-o);
    
    for(i = 0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        d += s.x;
        //d += min(s.x, .1);
    }
        
    if(s.x < 1.e-4)
    {
        normal(x,n, 2.e-5);
        vec3 l = normalize(x+.1*n);
        
        if(s.y == 2.)
        {
            col = vec3(0.86,0.21,0.13);
            col = .1*col
                            + .1*col * abs(dot(l,n))
                            + .5 * col * abs(pow(dot(reflect(-l,n),dir),2.));
            
            vec3 c1 = c.yyy;
            for(float fraction = 0.; fraction <= 4.; fraction += 1.)
    		{
                o = x;
                dir = refract(dir,n,.9);
                //dir = reflect(dir,n);
                d = 1.e-2;

                for(i = 0; i<N; ++i)
                {
                    x = o + d * dir;
                    scene(x,s);
                    if(s.x < 1.e-4)break;
                    d += s.x;
                    //d += min(s.x, .01);
                }

                if(s.x < 1.e-4)
                {
                    normal(x,n, 2.e-4);
                    vec3 l = normalize(x+.1*n);

                    if(s.y == 2.)
                    {
                        //c1 = .7*c.xxx;
                        c1 = (fraction == 0.)?vec3(0.86,0.21,0.13):
                        	(fraction == 1.)?vec3(0.85,0.80,0.62):
                        	(fraction == 2.)?vec3(0.22,0.25,0.25):
                        	(fraction == 3.)?vec3(0.16,0.17,0.17):
                        vec3(0.12,0.12,0.13);
                        c1 = .1*c1
                            + .4*c1 * abs(dot(l,n))
                            + 5.8 * c1 * abs(pow(dot(reflect(-l,n),dir),2.));
                    }
					//col = clamp(col, 0., 1.);
                    col = mix(col, c1, .15);
                }

                col = clamp(col, 0., 1.);
            }
        }
        
    }
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}


void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
