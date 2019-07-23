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

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float iScale, nbeats;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

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

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
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

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

// iq's smooth minimum
void smoothmin(in float a, in float b, in float k, out float dst)
{
    float h = max( k-abs(a-b), 0.0 )/k;
    dst = min( a, b ) - h*h*h*k*(1.0/6.0);
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = mix(sda, sdb, step(sdb.x, sda.x));
}  

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    float n;
    lfnoise(x-iTime*c.xy, n);
    
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

void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

vec2 ind=c.yy;
void scene(in vec3 x, out vec2 sdf)
{    
    
    float d,
        s = .1;
    vec2 x2 = mod(x.xz,s)-.5*s;
	
    ind = (x.xz - x2)/s;
    dbox(x2, .5*s*c.xx, d);
    zextrude(x.y+.4, -d-.005, .05, d);
    d = max(x.y,d);
    d = abs(d);
    sdf = vec2(d,2.);
    
    //sdf = c.xy;
    //sdf = vec2(x.y+.3, -2.);
    
    /**
    mat3 RR;
    rot3(-pi/4.*c.yxy, RR);
    x = RR * x;
    //*/
    
    for(float size = .1; size < .35; size += .1)
    {
        dbox3(x, size*c.xxx, d);
        stroke(d, .001, d);
        vec2 sda = vec2(d,3.+size);

        float n;
        vec3 y = mod(x,.25*size)-.5*.25*size,
            yi = (x-y)/size;
		ind = yi.xy+yi.yz+yi.xz;
        lfnoise(1.6*ind+5.*size-1.1*iTime, n);
        
        if(n>-.3)
        {
            
            dbox3(y, .25*size*c.xxx, d);
            sda.x = max(sda.x, -d);
        }
        
        add(sdf, sda, sdf);
        
        
    }   
    
    //add(sdf, vec2(x.x+.5,1.), sdf);
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

void colorize(in vec2 x, out vec3 col)
{
    col = .5*c.xxx;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    mat3 RR;
    rot3(-pi/4.*c.yxy, RR);
    
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 5./29., iScale));
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    vec3 col = c.yyy, 
        o = RR * c.yyx,
        r = normalize(c.xyy), 
        u = normalize(c.yxy),
        t = c.yyy, 
        dir,
        n,
        x,
        size = .35*c.xxx;
        
    int N = 950,
        i = 0;
    t = uv.x * r + uv.y * u;
    t = RR * t;
    dir = normalize(t-o);
    
    vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes
    vec2 abxlo = abs(o.yz + tlo.x*dir.yz),
        abylo = abs(o.xz + tlo.y*dir.xz),
        abzlo = abs(o.xy + tlo.z*dir.xy);

    vec4 dn = 100.*c.xyyy;
    float d = 0., ra;
    
    dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));
    dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));
    dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));
    
    d = dn.r;
    /*
    if(d < 10.) 
    {
        fragColor = c.xyyx;
        return;
    }*/
    
    int inside = 0;
    if(d < 2.)
    {
        inside = 1;
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4)break;
            if(x.y<-.5)
            {
                col = .2*c.xxx;
                i = N;
                break;
            }
            //d += s.x;
            //d += exp(-1.5e0*d)*s.x;
            d += s.x<1.e-1?min(s.x,1.e-3):s.x;
        }
    }
    if(d > 2. || i==N) 
    {
        d = -(o.y+.37)/dir.y;
        x = o + d * dir;
        scene(x, s);
        
        for(i = 0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x,s);
            if(s.x < 1.e-4)break;
            if(x.y<-.5)
            {
                col = .2*c.xxx;
                i = N;
                break;
            }
            d += s.x;
            //d += exp(-1.5e0*d)*s.x;
            //d += s.x<1.e-1?min(s.x,1.e-3):s.x;
        }
    }
    
    if(i < N)
    {
        normal(x,n, 5.e-4);
        
        lfnoise(x.xy*vec2(3.,8.)-iTime, ra);
        
        vec3 f;
        float dd = 5.e-1;
        
        r = RR*c.xyy;
        f = RR*c.yzy;
        u = RR*c.yyx;
        
        vec3 dp = abs(vec3(dot(n,r), dot(n,f), dot(n,u)));
        
        if(dp.y < dd && dp.z < dd) n = r;
        else if(dp.x < dd && dp.z < dd) n = f;
        else if(dp.x < dd && dp.y < dd) n = u;
        else s.y = -1.;
        /*
        if(dp.x+dp.y < .1*dd) ;
        else if(dp.y+dp.z < .1*dd) ;
        else if(dp.z+dp.x < .1*dd) ;
        else s.y = -1.;
        */
        vec3 l = normalize(x+.2*mix(c.yxx,c.yxz,iScale));
        
        
        if(s.y == 1.)
        {
            colorize(x.xy, col);
            col = .1*col
                + 1.*col * abs(dot(l,n))
                + 1.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));
            
        }
        else if(s.y == 2.)
        {
            
            col = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));
            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*ra*x));
            col = mix(col, c1, .5+.5*ra);
            col = .3*c.xxx;
            col = .1*col
                + .8*col * abs(dot(l,n))
                + 1.*col * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
        else if(s.y >= 3.)
        {
            lfnoise(.1*ind+3.*s.y*c.xx-iTime, col.x);
            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, col.y);
            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, col.z);
            float na;
            ind = x.xy + x.yz + x.zx;
            ind = mod(ind, .01)-.005;
            ind = x.xy + x.yz + x.zx-ind;
            rand(ind, na);
            na = 0.;
            col = .8+.2*col;
            col *= na;
            col = mix(.1,mix(.2,.4,iScale),step(na,.05))*col
                + mix(.1,.2,step(na,.95))*col * abs(dot(l,n))
                + .5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));
        }
        else if(s.y == -1.)
        {
            lfnoise(.1*ind+3.*s.y*c.xx-iTime, col.x);
            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, col.y);
            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, col.z);
            col = .8+.2*col;
            col = .1*col
                + 1.*col * abs(dot(l,n))
                + 1.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));
        }
        else if(s.y == -2.)
        {
            float s = .05;
            col = .5*c.xxx;
            vec2 dd = mod(x.xz, s)-.5*s;
            stroke(dd.x, .005, dd.x);
            stroke(dd.y, .005, dd.y);
            col = mix(col, c.xxx, sm(min(dd.x, dd.y)));
        }

    }
    
    if(inside == 1)
    {
        vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));
    	col = mix(col, c1, .3);
    }
    
    col = mix(col, c.yyy, smoothstep(2., 3., d));
    
    col *= 1.2;
    //col = .3*sqrt(col);
    //col *= col;
    col *= col*col;
    //col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    
    //col *= mix(col, length(col)/sqrt(3.)*c.xxx, iScale);

    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
