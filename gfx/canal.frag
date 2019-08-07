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

void dsmoothvoronoi(in vec2 x, out float d, out vec2 z)
{
    float n;
//     lfnoise(x-iTime*c.xy, n);
    
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
            smoothmin(ret, d, .05, ret);
        }
    
    d = ret;
    z = pf;
}

vec2 ind;
void scene(in vec3 x, out vec2 sdf)
{
    x.z -= 1.3*iTime;
    
    float dx,
        d, v;
    
    lfnoise(.5*x.z*c.xx, dx);
    x.xy-=.4*dx*c.xy;
    
    // Voronoi
    float phi = atan(x.y,x.x);
    dsmoothvoronoi(2.*vec2(mod(phi+pi/4., 2.*pi), x.z), v, ind);
    stroke(v, .01, v);
    d = length(x.xy) - mix(1.,1.1, smoothstep(.0,.2,v));
    
    zextrude(length(x.xy)-1.0,d,.05, d);
    d -= .05;
    sdf = vec2(d,1.);
    
    // Smaller voronoi
    dsmoothvoronoi(8.*vec2(mod(phi+pi/4., 2.*pi), x.z), v, ind);
    stroke(v, .02, v);
    d = length(x.xy) - mix(1.1,1.2, smoothstep(.0,.2,v));
    
    zextrude(length(x.xy)-1.11,d,.01, d);
    d -= .1;
    add(sdf, vec2(d,1.), sdf);
    smoothmin(d, sdf.x , .1, sdf.x);
    
    float dy = 0.;
    //lfnoise(33.*x.x*c.xx-iTime,dy); 
    
    for(float i=1.; i<=15.; i+=1.)
    {
        float f, a;
        vec2 dir;
        
        float n;
        
        lfnoise(i*c.xx-n, f);
        f = .5+.5*f;
        lfnoise(i*c.xx+1337.-n, a);
        a = .5+.5*a;
        lfnoise(i*c.xx+2337.-2.*n, dir.x);
        lfnoise(i*c.xx+3337.-3.*n, dir.y);
        dir = mix(c.yx-.2*c.xy, x.yx+.2*c.xy, 2.*dir);
        dir = normalize(dir);
        
        
        float dya = pow(1.01,f)* a * sin(-2.e-3*2.*pi*pow(1.95,abs(f+.01*a))*(1.*f-.01*a)*iTime-2.e-4*2.*pi*pow(1.99,abs(i-.1*a))*dot(dir,vec2(.5,4.)*(2.*(x.xz+1.3*(iTime))*c.yx)));
    	dy += 2.*pow((dya+1.)/2.,4.)-1.;
    }
    dy = .4*dy;
    
    add(sdf, vec2(x.y+.4+.001*dy, 2.), sdf);
    //smoothmin(x.y+.4-.02*dy, sdf.x , .5, sdf.x);
    
    

}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 5.e-4;
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
    
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    
    float dx, dx2, d0;
    lfnoise(-.5*1.3*iTime*c.xx, dx);
    lfnoise(-.5*1.3*(iTime+1.e-3)*c.xx, dx2);

    vec3 col = c.yyy, 
        o = c.yyx+.1*c.yxy+.4*dx*c.xyy,
        r = c.xyy, 
        u = c.yxy, 
        t = c.yyy+.4*dx2*c.xyy, 
        dir,
        n,
        x;
    int N = 400,
        i, a = 0;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = .5/length(dir.xy);// -(o.z-.2)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(length(x.xy-.4*dx*c.xy)>1.5)
        {
            col = c.yyy;
            i = N;
            break;
        }
        d += min(s.x,2.e-2);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n);
        vec3 l = normalize(x+.5*n);
       
		if(s.y == 1.)
        {
            col = vec3(0.76,0.20,0.23);
            
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
            
            col = mix(col, 2.*vec3(0.76,0.20,0.13), smoothstep(0.658, 1.02, 1.-abs(dot(n, c.yyz))));
        	col = mix(col, vec3(0.96,0.7,0.423), smoothstep(0.658, 1.02, abs(dot(n, c.yyz))));
        }
        else if(s.y == 2.) // Mirror material
        {
            col = .3*c.xxx;
            col = .2*col
                + .2*col * abs(dot(l,n))
                + .3*col * pow(abs(dot(reflect(-l,n),dir)),2.);
            
            
            N = 50;
            o = x;
            dir = reflect(dir, n);
            d0 = d;
            d = 1.e-2;
            vec3 c1 = c.yyy;
            
            for(i = 0; i<N; ++i)
            {
                x = o + d * dir;
                scene(x,s);
                if(s.x < 1.e-4)break;
                if(length(x.xy)>1.5)
                {
                    c1 = c.yyy;
                    i = N;
                    break;
                }
                //d += min(s.x,5.e-3);
                d += s.x;
            }
            
            if(i < N)
            {
                normal(x,n);
                vec3 l = normalize(x+.5*n);

                if(s.y == 1.)
                {
                    c1 = vec3(0.76,0.20,0.23);
                    c1 = .2*c1
                        + .2*c1 * abs(dot(l,n))
                        + .6*c1 * pow(abs(dot(reflect(-l,n),dir)),3.);
                    c1 = mix(c1, 2.*vec3(0.76,0.20,0.13), smoothstep(0.658, 1.02, clamp(1.-abs(dot(n, c.yyz)),0.,1.)));
        			c1 = mix(c1, vec3(0.96,0.7,0.423), smoothstep(0.658, 1.02, abs(dot(n, c.yyz))));
                }
            }
            //c1 = mix(c1, c.yyy, smoothstep(3.,6.,d));
            
            col = mix(col, c1, .35);
            col = mix(col, .1*c.yxx, .3);
            //a = 1;
        }

    }
    
    col = mix(col, c.yyy, smoothstep(2.,22.,d+d0));
    //col = mix(col, 2.*col, iScale);
    
    float nn;
    lfnoise(12.*(x.z-1.3*iTime)*c.xx, nn);
    
    col *=mix(1.1,2.6,mix(.5+.5*nn,1.,iScale));
    col *= col;
    col = clamp(col, 0., 1.);
    if(col == c.xxx) col = c.yyy;
    
//     col = mix(col,c.yyy, .1);
    
    col = mix(col, mix(col, length(col)/sqrt(3.)*c.xxx, .7), iScale);
    col = mix(col, 1.7*1.7*col*col, iScale);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
