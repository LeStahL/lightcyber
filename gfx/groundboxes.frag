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
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float iScale, nbeats;

void rand(in vec2 x, out float n);
void lfnoise(in vec2 t, out float n);
void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);
void dbox3(in vec3 x, in vec3 b, out float d);
void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void scale(out float s);

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void smoothmin(in float a, in float b, in float k, out float dst);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
void dvoronoi(in vec2 x, out float d, out vec2 z);
void dbox(in vec2 x, in vec2 b, out float d);
void rot3(in vec3 p, out mat3 rot);

vec2 ind=c.yy;

void scene(in vec3 x, out vec2 sdf)
{    
    
    float d,
        s = .1;

    sdf = c.xy;
    
    for(float size = .0; size <= .5; size += .025)
    {
        dbox3(x, size*c.xxx, d);
        stroke(d, .001, d);
        vec2 sda = vec2(d,3.+size);

        float n;
        vec3 y = mod(x,.125*size)-.5*.125*size,
            yi = (x-y)/size;
		ind = yi.xy+yi.yz+yi.xz;
        lfnoise(3.6*ind+15.*size-1.1*iTime, n); // TODO: use shifts to 13.6 from 3.6
        
        if(n>-.3)
        {
            
            dbox3(y, .25*size*c.xxx, d);
            sda.x = max(sda.x, -d);
        }
        
        add(sdf, sda, sdf);
    }   
}  

void normal(in vec3 x, out vec3 n, in float dx);

void colorize(in vec2 x, out vec3 col)
{
    col = .5*c.xxx;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    mat3 RR;
    rot3(.2*iTime*vec3(1.1,1.4,1.6), RR);
    
    scale(iScale);
    
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
        size = .301*c.xxx;
        
    int N = 150,
        i = 0;
    t = uv.x * r + uv.y * u;
    t = RR * t;
    dir = normalize(t-o);
    float d = 0., ra, ss, inside = 0., dd;
    
    for(ss = .5; ss >= .0; ss -= .025)
    {
        size = ss*c.xxx+1.e-4;

        vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes
        vec2 abxlo = abs(o.yz + tlo.x*dir.yz),
            abylo = abs(o.xz + tlo.y*dir.xz),
            abzlo = abs(o.xy + tlo.z*dir.xy);

        vec4 dn = 100.*c.xyyy;
        
        dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));
        dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));
        dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));

        inside += .05;
        
        if(ss == 3.)dd = dn.r;
        d = dn.r;
        
        x = o + d * dir;
        scene(x, s);
        if(s.x < 1.e-4)break;
        
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
        s.y = -1.;

        vec3 l = normalize(x+.03*normalize(x-o));
        
        vec3 c1 = col;
        
        if(s.y == 1.)
        {
            colorize(x.xy, c1);
            c1 = .1*c1
                + 1.*c1 * abs(dot(l,n))
                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));
            
        }
        else if(s.y == 2.)
        {
            
            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));
            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*ra*x));
            c1 = mix(c1, c1, .5+.5*ra);
            c1 = .3*c.xxx;
            c1 = .1*c1
                + .4*c1 * abs(dot(l,n))
                + .3*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
        else if(s.y >= 3.)
        {
            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);
            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);
            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);
            float na;
            ind = x.xy + x.yz + x.zx;
            ind = mod(ind, .01)-.005;
            ind = x.xy + x.yz + x.zx-ind;
            rand(ind, na);
            na = 0.;
            c1 = .8+.2*c1;
            c1 *= na;
            c1 = mix(.1,mix(.2,.4,iScale),step(na,.05))*c1
                + mix(.1,.2,step(na,.95))*c1 * abs(dot(l,n))
                + .5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));
        }
        else if(s.y == -1.)
        {
            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);
            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);
            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);
            c1 = .8+.2*c1;
            c1 = .1*c1
                + 1.*c1 * abs(dot(l,n))
                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),3.));
                
            vec3 dc;

            vec3 zz = mod(x, .025)-.5*.025, zi = x-zz;
            rand(zi.xy+zi.yz+zi.zx, dc.x);
            rand(zi.xy+zi.yz+zi.zx+1337., dc.y);
            rand(zi.xy+zi.yz+zi.zx+2337., dc.z);
            float da;
            dbox3(zz, .01*c.xxx, da);
            stroke(da, .001, da);
            c1 = mix(c1, 1.2*c1*c1+.2*dc, sm(da));
            stroke(da-.002, .001, da);
            c1 = mix(c1, 1.6*c1*c1, sm(da));
        }
        else if(s.y == -2.)
        {
            float s = .05;
            c1 = .5*c.xxx;
            vec2 dd = mod(x.xz, s)-.5*s;
            stroke(dd.x, .005, dd.x);
            stroke(dd.y, .005, dd.y);
            c1 = mix(c1, c.xxx, sm(min(dd.x, dd.y)));
        }
        col = mix(col, c1, mix(.2,.6,iScale));
        col = .9*col;
    }
    
    if(s.x > 1.e-4 && uv.y<-2.e-4)
    {
        d = -(o.y+.375)/dir.y;
        x = o + d * dir;
        scene(x, s);
        i=N;
    }
    else x = o + d * dir;
    
    if(s.y == -1.)
    {
        vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));
    	
        col = mix(col, c1, .3+.5*inside);
    }
    
    col = mix(col, c.yyy, smoothstep(2., 3., d));
    
    col *= 1.2;
    col *= col;
    
    col *= col*col;
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
