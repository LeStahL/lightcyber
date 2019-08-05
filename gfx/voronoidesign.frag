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

/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
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

// Update1: Changes implementing FabriceNeyret2's comments.

// Global constants
const vec3 c = vec3(1.0, 0.0, -1.0);
const float pi = acos(-1.);

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW
void hash31(in float p, out vec3 d)
{
   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
   p3 += dot(p3, p3.yzx+33.33);
   d = fract((p3.xxy+p3.yzz)*p3.zyx); 
}
// End of Hoskins Hash

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)
{
    n = 0.;
    float a = 1., nf = 0., buf;
    for(float f = d; f<b; f *= 2.)
    {
        lfnoise(f*x, buf);
        n += a*buf;
        a *= e;
        nf += 1.;
    }
    n *= (1.-e)/(1.-pow(e, nf));
}

// Box sdf
void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)
{
    vec3 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
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

//distance to spline with parameter t
float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum dist3ance to spline
void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,
    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        ds = dist3(p0,p1,p2,x,ui.x+ui.y-tau);
        return;
    }
    
    //three dist3inct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    ds = min(
        dist3(p0,p1,p2,x, t.x),
        min(
            dist3(p0,p1,p2,x,t.y),
            dist3(p0,p1,p2,x,t.z)
        )
    );
}

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    z = pf;
}

vec2 vind,vind2;
float v, fn;
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
    
    float phi = atan(y.y, y.x),
        dp = pi/24.,
        phii = mod(phi, dp)-.5*dp,
        pa = phi - phii, 
        R1 = .05,
        R2 = .4;
    
    R2 = mix(R1, R2, .5+.5*n2);
    
    float r0, r1;
    rand(pa*c.xx, r0);
    r0 = mix(r0,.5+.5*n,.5);
    
    
    rand(vind, r1);
    
    dspline3(y, vec3(R1*cos(pa), R1*sin(pa), -.5), vec3(R1*cos(pa), R1*sin(pa), .1*r1), vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), sdf.x);
    float da;
    dspline3(y, vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0), da);
    sdf.x = min(sdf.x, da);
    stroke(sdf.x, .25*mix(.02,.05, .5+.5*n2), sdf.x);
    sdf.y = 2.;
    
    add(sdf, vec2(length(y-vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0))-.01, 3.), sdf);
    
    float fa;
    lfnoise(4.*a,  fa);
    dvoronoi(a,fn, vind2); 
    add(sdf, vec2(x.z+.4+.1*mix((v+fn),fa,.5),4.), sdf);
    //x.xy = x.yx;
    /*
    float t, ta;

    x.y += .3*iTime;

    
    lfnoise(5.*x.y*c.xx, t);
    lfnoise(5.*(x.y*c.xx-iTime), ta);
    
    
    mat2 R = mat2(cos(ta), sin(ta), -sin(ta), cos(ta));
    vec2 z = R * x.xz;
    x.x = z.x;
    x.z = z.y;

    float size = .4,
        r = (.5-.3*t)*size,
        phi = atan(x.z,x.x),
        modsize = pi/32.,
        phii = (phi-mod(phi, modsize)+.5*modsize);
    vec3 dx = vec3(r*cos(phii), 0., r*sin(phii));
    
    x -= dx;
    
    dlinesegment3(x, 1.e3*c.yxy, -1.e3*c.yxy, sdf.x);
    stroke(sdf.x, .01, sdf.x);
    sdf.y = 2.;*/
    
/*    
    vec3 y = vec3(x.x, mod(x.y,size)-.5*size, x.z),
        dxi, 
        dxip1, 
        dxim1;
    float yi = (y.y-x.y)/size,
        d;
    
    
    
    hash31(yi, dxi);
    hash31(yi-1., dxim1);
    hash31(yi+1., dxip1);
    
    
    dlinesegment3(y, size*(.5*dxi-.25), size*(.5*dxip1-.25-c.yxy), sdf.x);
    dlinesegment3(y, size*(.5*dxi-.25), size*(.5*dxim1-.25+c.yxy), d);
    d = min(d, sdf.x);
    stroke(d, .01, sdf.x);
    sdf.y = 2.;
*/
}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 1.e-4;
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

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, 
        s;
    
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

    float d = -(o.z-.1)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        /*
        if(x.z<-.05)
        {
            col = .2*c.xxx;
            i = N;
            break;
        }*/
        d += min(s.x,3.e-2);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n);
       
        if(s.y == 3.)
        {
            vec3 l = normalize(x+c.yyx);
            float r;
            
            col = mix(c.xxx, vec3(0.76,0.20,0.13), .2);
            //col = mix(vec3(0.76,0.20,0.13), vec3(0.23,0.23,0.23),  clamp((x.z-.2*length(x.xy))/.1,0.,1.));
            col = .4*col
                + .4*col * abs(dot(l,n))
                + .2*col * pow(abs(dot(reflect(-l,n),dir)),3.);
        }
        else if(s.y == 4.)
        {
            vec3 l = normalize(x+c.yyx);
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
            vec3 l = normalize(x+c.yyx);
            float r;
            rand(vind, r);
            
            col = mix(mix(vec3(0.76,0.20,0.23),vec3(.18,.32,.13), r), vec3(0.23,0.23,0.23),  clamp((x.z)/.05,0.,1.));
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
        }

    }
    
    col = mix(col, 0.*.23*c.xxx*vec3(0.76,0.20,0.13),smoothstep(1.,5.,d));
    col *=3.6;
    //col = mix(col, 2.*col, iScale);
    
    //col = atan(col);
    col *= col;
    fragColor = vec4(clamp(col,0.,1.),1.0);
}	

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
