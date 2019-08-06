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

// Global constants
const vec3 c = vec3(1.0, 0.0, -1.0);
const float pi = acos(-1.);

float iScale;

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

void dbox3(in vec3 x, in vec3 b, out float d)
{
  	vec3 da = abs(x) - b;
  	d = length(max(da,0.0)) + min(max(da.x,max(da.y,da.z)),0.0);
}

// Compute distance to regular polygon
void dpolygon(in vec2 x, in float N, out float d)
{
    d = 2.0*pi/N;
    float t = mod(acos(x.x/length(x)), d)-0.5*d;
    d = -0.5+length(x)*cos(t)/cos(0.5*d);
}

// Distance to circle
void dcircle(in vec2 x, out float d)
{
    d = length(x)-1.0;
}

// Distance to line segment
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// 2D rotational matrix
void rot(in float phi, out mat2 m)
{
    vec2 cs = vec2(cos(phi), sin(phi));
    m = mat2(cs.x, -cs.y, cs.y, cs.x);
}

// Distance to pig ear
void dear(in vec2 x, out float d)
{
    d = abs(2.*x.y)
        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))
        -.5*min(-abs(x.x),.01);
}

// Distance to a triangle
void dtriangle(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float d)
{
    vec2 d1 = c.xz*(p1-p0).yx, d2 = c.xz*(p2-p1).yx, d3 = c.xz*(p0-p2).yx;
    d = -min(
        dot(p0-x,d1)/length(d1),
        min(
            dot(p1-x,d2)/length(d2),
            dot(p2-x,d3)/length(d3)
        )
    );
}

// Distance to schnappsgirls logo in hexagon
void dschnappsgirls(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Dress
    dtriangle(x, vec2(-.1,-.3), vec2(.5,-.3), vec2(.2, .6), d0);
    dlinesegment(x, vec2(-.1,.325), vec2(.5,.325), da);
    stroke(da,.06,da);
    d0 = max(d0,-da);
    
    // Head
    dcircle(7.*(x-vec2(.2,.5)), da);
    d0 = max(d0, -da+.5);
    d0 = min(d0, da/7.);
    
    // Legs
    dlinesegment(x, vec2(.125,-.3), vec2(.125,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.275,-.3), vec2(.275,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    
    // Shoulders
    dlinesegment(x, vec2(0.05,.25), vec2(.35,.25), da);
    stroke(da, .085, da);
    d0 = min(d0, da);
    
    // Arms
    dlinesegment(x, vec2(.385,.25), vec2(.5, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.017,.25), vec2(-.1, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    
    // Glass
    dtriangle(x, vec2(-.6,.3), vec2(-.4,.1), vec2(-.2,.3), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(-.4,.15), vec2(-.4,-.1), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dtriangle(x, vec2(-.5,-.15), vec2(-.3,-.15), vec2(-.4,-.1), da);
    d0 = min(d0, da);
    
    // Liquid
    dtriangle(x, vec2(-.55,.25), vec2(-.4,.1), vec2(-.25,.25), da);
    d0 = min(d0, da);
    
    // Salad
    dlinesegment(x, vec2(-.4,.1), vec2(-.2,.5), da);
    stroke(da, .01, da);
    d0 = min(d0, da);
    dcircle(24.*(x-vec2(-.3,.3)), da);
    d0 = min(d0, da/24.);
    dcircle(24.*(x-vec2(-.25,.4)), da);
    d0 = min(d0, da/24.);
    
    d = max(d, -d0);
}

// Distance to spacepigs logo in hexagon
void dspacepigs(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Head
    dcircle(2.5*x,d0);
    d0 /= 2.5;
    
    // Ears
    dear(vec2(2.,5.)*x-vec2(.8,1.3), da);
    d0 = min(d0,da/10.);
    dear(vec2(2.,5.)*x+vec2(.8,-1.3), da);
    d0 = min(d0,da/10.);
    
    // Nose
    dcircle(6.*x-vec2(0.,-.5),da);
    d0 = max(d0,-da/6.);
    dcircle(24.*x-vec2(-1.5,-2.),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(1.5,-2.),da);
    d0 = min(d0,da/24.);
    
    // Eyes
    dcircle(16.*x-vec2(-3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(16.*x-vec2(3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(24.*x-vec2(-5.,3.5),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(5.,3.5),da);
    d0 = min(d0,da/24.);
    
    d = max(d, -d0);
}

// Distance to kewlers logo in hexagon
void dkewlers(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    x *= 1.2;
    
    dbox(x-vec2(0.,-.3),vec2(.6,.1),d0);
    dbox(x-vec2(-.5,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    dbox(x-vec2(-.5+1./3.,.25),vec2(.1,.5),da);
    d0 = min(d0,da);
    dbox(x-vec2(-.5+2./3.,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    dbox(x-vec2(.5,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    
    d = max(d, -d0);
}

// Distance to farbrausch logo in hexagon
void dfarbrausch(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    x += vec2(.1,0.);
    x *= 1.2;
    
    dlinesegment(x,vec2(-.65,.05),vec2(-.5,.05),d0);
    dlinesegment(x,vec2(-.5,.05),vec2(-.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.2,-.49),vec2(-.0,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.0,-.49),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.07,0.),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.2,-.49),vec2(-.07,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.6,-.15),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.53,0.),vec2(.6,-.15),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.53,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.68,-.04),da);
    d0 = min(d0, da);
    dpolygon(5.*(x+vec2(.3,.65)),6.,da);
    d0 = min(d0, da/5.);
    dpolygon(5.*(x+vec2(-.5,.65)),6.,da);
    d0 = min(d0, da/5.);
    
    stroke(d0,.035, d0);
    d = max(d, -d0);
}

// Distance to haujobb logo in hexagon
void dhaujobb(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    mat2 m;
	rot(.3,m);
    x = 1.1*x*m;
    x.x *= 1.1;
        
    x += vec2(-.05,.2);
    
    // Left leg
    dbox(x+.35*c.xx,vec2(.1,.05),d0);
    dbox(x+vec2(.3,.25),vec2(.05,.15),da);
    d0 = min(d0,da);
    dbox(x+vec2(.2,.15),vec2(.1,.05),da);
    d0 = min(d0,da);
    dbox(x+vec2(.15,.05),vec2(.05,.15),da);
    d0 = min(d0,da);
    
    // Right leg
    dbox(x-vec2(.65,.35),vec2(.05,.15),da);
    d0 = min(d0,da);

    // Torso
    rot(.2, m);
    dbox(m*(x-vec2(.25,.15)),vec2(.45,.05),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(-.15,.35)),vec2(.45,.05),da);
    d0 = min(d0,da);
    rot(pi/8.,m);
    dbox(m*(x-vec2(.0,.25)),vec2(.1,.15),da);
    d0 = min(d0,da);
    
    // Penis
    dbox(m*(x-vec2(.1,-.0)),vec2(.025,.1),da);
    d0 = min(d0,da);
    
    // Left hand
    rot(.3,m);
    dbox(m*(x-vec2(.235,.535)),vec2(.035,.15),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(.225,.7)),vec2(.075,.025),da);
    d0 = min(d0,da);
    
    // Right hand
    rot(-.2,m);
    dbox(m*(x+vec2(.585,-.2)),vec2(.0375,.1),da);
    d0 = min(d0,da);
    
    // Head
    dcircle(6.*(x-vec2(-.15,.58)),da);
    d0 = min(d0,da/6.);
    
    d0 -= .05*(abs(x.x)+abs(x.y)-.2);
    d = max(d,-d0);
}

// Distance to mercury logo in hexagon
void dmercury(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da;

    x += .1*c.yx;

    // Upper part
    dbox(x-.35*c.yx,vec2(.4,.35), da);
    d = max(d, -da);
    dbox(x-.7*c.yx, vec2(.2,.2), da);
    d = min(d,da);
    dbox(x-.25*c.yx,vec2(.2,.05),da);
    d = min(d,da);
    
    // Lower part
    dbox(x+.2*c.yx,vec2(.1,.4),da);
    d = max(d, -da);
    dbox(x+.2*c.yx, vec2(.4,.1),da);
    d = max(d, -da);
}

vec2 ind, indc;
void scene(in vec3 x, out vec2 sdf)
{
    x.z -= .05*iTime;
    
    float d;
    
    // Corridor
    dbox3(x, vec3(.1,.1,1.e3), d);
    sdf = vec2(-d, 2.);
    
    // Wall tiles
    float distortion;
	lfnoise(5.2e2*x.yz, distortion);
    float tsize = .005,
    	dy = mod(x.y, tsize)-.5*tsize,
        yi = (x.y-dy)/tsize,
        zpar = x.z+mix(0., .5*tsize, mod(yi,2.)),
        dz = mod(zpar, tsize)-.5*tsize,
        zi = (zpar-dz)/tsize;
    dbox3(vec3(abs(x.x)-.1, dy, dz), vec3(.0005+.00001*distortion, .39*tsize*c.xx), d);
    add(sdf, vec2(d, 3.), sdf);
    smoothmin(sdf.x, d, .001, sdf.x);
    ind = vec2(yi, zi);
    
    // Ceiling 
    tsize = .025;
    dz = mod(x.z, tsize)-.5*tsize;
    float dx = mod(x.x, tsize)-.5*tsize;
    zi = (x.z-dz)/tsize;
    float xi = (x.x-dx)/tsize;
    dbox3(vec3(dx, abs(x.y)-.1, dz), vec3(.48*tsize, .0005, .48*tsize), d);
    add(sdf, vec2(d, 4.), sdf);
    smoothmin(sdf.x, d, .001, sdf.x);
    indc = vec2(xi, zi);
    
    // Tentacles
    /*
    // Generate Block coordinates
    tsize = .02;
    vec2 xz = mod(x.xz, tsize)-.5*tsize,
        xzi = (x.xz-xz)/tsize;
    
    // Determine lower and upper rotation
    vec2 rotation;
    rand(xzi, rotation.x);
    rand(xzi-1337., rotation.y);
    rotation = pi*floor(4.*rotation)/2.;
    
    // Determine graph
    float graph;
    rand(xzi-2337., graph);
    
    vec2 lower = vec2(-.5*tsize, 0.),
        upper = vec2(-.5*tsize, 0.);
    mat2 Rlower = mat2(cos(rotation.x), sin(rotation.x), -sin(rotation.x), cos(rotation.x)),
        Rupper = mat2(cos(rotation.y), sin(rotation.y), -sin(rotation.y), cos(rotation.y));
    
    lower = Rlower * lower;
    upper = Rupper * upper;
    
    if(graph < 1./3.)
    {
        
    }
*/
    
    // Logos
    tsize = .25;
    float tw = .0005;
    dz = mod(x.z-.5*tsize, tsize)-.5*tsize;
    zi = floor((x.z-dz)/tsize);
    zi = mod(zi, 6.);
    if(zi < .5)
    {
        dmercury(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }
    else if(zi < 1.5)
    {
        dhaujobb(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }
    else if(zi < 2.5)
    {
        dfarbrausch(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }
    else if(zi < 3.5)
    {
        dkewlers(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }
	else if(zi < 4.5)
    {
        dspacepigs(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }
	else if(zi < 5.5)
    {
        dschnappsgirls(20.*x.xy, d);
        stroke(d/20.,tw, d);
        zextrude(dz, -d, .005, d);
        add(sdf, vec2(d, 5.), sdf);
    }

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
    iScale = 1.;
    
    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, 
        s;
    
    //uv *= 2.;
    
    vec3 col = c.yyy, 
        o = c.yyx,
        r = c.xyy, 
        u = c.yxy, 
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 400,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = 0.;//-(o.z-.1)/dir.z;
    
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
        vec3 l = mix(1.5,.2,abs(pow(sin(2.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);
       
		if(s.y == 2.)
        {
            col = .23*c.xxx;
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
        }
        else if(s.y == 3.)
        {
            float r;
            rand(ind, r);
            
            col = .2*vec3(0.02,0.11,0.24)
                + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))
                + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);
            
            // Reflections
            d = 1.e-2;
            o = x;
            dir = reflect(dir, n);
            
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
                d += min(s.x,3.e-1);
                //d += s.x;
            }
            
            normal(x,n);
        	vec3 l = mix(1.5,.2,abs(pow(sin(4.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);
       
            vec3 c1;
            if(s.y == 2.)
            {
                c1 = .23*c.xxx;
                c1 = .2*c1
                    + .2*c1 * abs(dot(l,n))
                    + .6*c1 * pow(abs(dot(reflect(-l,n),dir)),3.);
            }
            else if(s.y == 3.)
            {
                float r;
                rand(ind, r);

                c1 = .2*vec3(0.02,0.11,0.24)
                    + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))
                    + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);
            }
            else if(s.y == 4.)
            {
                float r;
                rand(indc, r);

                c1 = .2*.2*c.xxx
                    + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))
                    + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);
            }
            else if(s.y == 5.)
            {
                c1 = .2*.2*c.xyy
                    + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))
                    + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);
                c1 = mix(c1, c.xxx, .1);
            }
            
            col = mix(col, c1, .5);
        }
        else if(s.y == 4.)
        {
            float r;
            rand(indc, r);
            
            col = .2*.2*c.xxx
                + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))
                + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);
        }
        else if(s.y == 5.)
        {
            col = .2*.2*c.xyy
                + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))
                + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);
            col = mix(col, c.xxx, .1);
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
