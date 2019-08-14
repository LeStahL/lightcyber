/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19
* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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

float nbeats;
float iScale;

// Global constants
const vec3 c = vec3(1.0, 0.0, -1.0);
const float pi = acos(-1.);

void scale(out float s);
void dsmoothvoronoi(in vec2 x, out float d, out vec2 z);
void rand(in vec2 x, out float n);
void hash31(in float p, out vec3 d);
void lfnoise(in vec2 t, out float n);
void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);
void dbox(in vec2 x, in vec2 b, out float d);
void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
void smoothmin(in float a, in float b, in float k, out float dst);
void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);
void dvoronoi(in vec2 x, out float d, out vec2 z);

vec2 vind,vind2;
float v, fn, r1, fb;
void scene(in vec3 x, out vec2 sdf)
{
    x.y += .2*iTime;
    
    dvoronoi(1.5*x.xy, v, vind);
    
    vec3 y = vec3(vind/1.5-x.xy,x.z);
    
    float n, n2;
    
    lfnoise(c.xx-.3*iTime+vind*3., n);
    lfnoise(5.*x.z*c.xx-iTime-vind*4., n2);
    n2 *= .2;
    
    mat2 RR = mat2(cos(n2), sin(n2), -sin(n2), cos(n2));
    vec2 a = x.xy;
    x.xy = RR * x.xy;
    rand(vind, r1);
    float r2;
    rand(vind -1336., r2);
    
    float phi = atan(y.y, y.x),
        dp = pi/24.,
        phii = mod(phi, dp)-.5*dp,
        pa = phi - phii, 
        R1 = .05,
        R2 = mix(.4,.25,1.-r2);
    
    R2 = mix(R1, R2, .5+.5*n2);
    
    float r0;
    rand(pa*c.xx, r0);
    r0 = mix(r0,.5+.5*n,.5);
    
    dspline3(y, vec3(1.4*R1*cos(pa), 1.4*R1*sin(pa), -.5), vec3(R1*cos(pa), R1*sin(pa), .1*r1), vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), sdf.x);
    float da;
    dspline3(y, vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0), da);
    sdf.x = min(sdf.x, da);
    stroke(sdf.x, .25*mix(.02,.05, .5+.5*n2), sdf.x);
    sdf.y = 2.;
    
    add(sdf, vec2(length(y-vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0))-.01, 3.), sdf);
    
    float fa;
    lfnoise(4.*a,  fa);
    dvoronoi(a,fn, vind2); 
    fa = x.z+.4+.1*mix((v+fn),fa,.5);
    add(sdf, vec2(fa,4.), sdf);
    smoothmin(sdf.x, fa, .1, sdf.x);
}

void normal(in vec3 x, out vec3 n, in float dx);

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void vs(in vec3 x, out vec2 sdf)
{
    vec2 vi;
    dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx), sdf.x, vi);
    sdf.x = x.z-.1-.2*sdf.x;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, 
        s;
    
    scale(iScale);
    
    uv *= 2.;
    
    vec3 col = c.yyy, 
        o = c.yzx,
        r = c.xyy, 
        u = normalize(c.yxx), 
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 400,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.3)/dir.z;
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        vs(x,s);
        if(s.x < 1.e-4)break;
        d += s.x;
    }
    float v1, rar, dx = 1.e-3;
    vec2 vi1, na;
    vec3 cv, l;
    x = o + d * dir;
    dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx), v1, vi1);

    rand(vi1, rar);
    cv = mix(c.yyy,vec3(.23,.23,.23), rar);
    v1 = abs(v1)-.01;
    cv = mix(cv, c.yyy, sm(v1));
    v1 = abs(v1-.01)-.005;
    cv = mix(cv, c.xxx, sm(v1));
    vs(x,s);
    vs(x+dx*c.xyy, na);
    n.x = na.x;
    vs(x+dx*c.yxy, na);
    n.y = na.x;
    vs(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
    l = normalize(x+c.yyx);
    cv = .2*cv
        +.2*cv*abs(dot(l,n))
        +.4*cv*pow(abs(dot(reflect(-l,n),dir)),3.);
    cv = mix(cv, 1.5*vec3(0.76,0.20,0.13), smoothstep(0.858, 1.02, dot(n, c.yyx)));
    dir = refract(dir, n, .98);
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;

        d += min(s.x,3.e-2);
    }
    
    if(i < N)
    {
        normal(x,n, 1.e-4);
       
        if(s.y == 3.)
        {
            l = normalize(x+c.yyx);
            float r;
            
            col = mix(c.xxx, vec3(0.76,0.20,0.13), .8);
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .8*col * pow(abs(dot(reflect(-l,n),dir)),2.);
        }
        else if(s.y == 4.)
        {
            l = normalize(x+c.yyx);
            float r;
            rand(vind+vind2,r);
            
            col = mix(.023*c.xxx, vec3(0.76,0.40,0.23), r);
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
            stroke(v, .01, v);
            stroke(fn, .01, fn);
            col = mix(col, c.yyy, sm(v));
            col = mix(col, col*col, sm(fn));
        }
		if(s.y == 2.)
        {
            l = normalize(x+c.yyx);
            float r;
            rand(vind, r);
            
            col = mix(mix(vec3(0.76,0.20,0.23),vec3(.18,.32,.13), r), vec3(0.23,0.23,0.23),  clamp((x.z)/r1/.1,0.,1.));
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
            col = mix(col, 5.*col, .25*n.x*n.x);
        }
        
    }
    
    col *=3.6;
    col *= col;

    if(s.y != 3.)
    {
        col = mix(length(col)/sqrt(3.)*c.xxx, col,.3);
    }
    
    col = mix(col,cv,.8);
    col = mix(col, c.yyy, smoothstep(1.,5.,d));
    
    col *= mix(1.,15.,mix(.28,.88, 0.*iScale));
    col *= col;
    
    col = mix(col, .01*col, smoothstep(-.6,-1.,uv.y));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}	

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
