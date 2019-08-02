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

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
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

void colorize(in vec2 x, out vec3 col)
{
    vec3 c1;
    vec2 ind,
        xv,
        xi;
    float d,
        vs = 16.,
        n,
        size = .1,
        xix = mod(x.x, size)-.5*size,
        xixj = (x.x - xix),
        ri,
        rim1,
        rip1,
        lines = 8.,
        da,
        op,
        s;
    
    // Background blending
    s = smoothstep(0.,.5,.5-abs(x.y));
    col = mix(c.yyy, vec3(0.02,0.16,0.25), s);
    
    // Background circles
    dvoronoi(vs*x, d, ind);
    xv = ind/vs-x;
    lfnoise(vec2(3.,33.)*ind/vs-3.*iTime*c.xy,n);
    n = .5+.5*n;
    d = length(xv)-mix(.0,.35,n)/vs;
    col = mix(col, n*vec3(0.07,0.29,0.42), sm(d));
    d = abs(d-.005) -.002;
    col = mix(col, (1.-n)*vec3(0.49,0.71,0.78), sm(d));
    
    for(float i = 1.; i < 9.; i += 1.)
    {
        rand(i*c.xx, op);
        op = .5+.5*round(16.*op)/16.;
        x += -.1+.2*op;
        
        xix = mod(x.x, size)-.5*size;
        xixj = (x.x - xix);
        
        // Edges
        lfnoise(2.e0*xixj*c.xx+14.*i, ri);
        lfnoise(2.e0*(xixj+size)*c.xx+14.*i, rip1);
        lfnoise(2.e0*(xixj-size)*c.xx+14.*i, rim1);

        float h = .2;
        
        ri = h*round(lines*ri)/lines;
        rip1 = h*round(lines*rip1)/lines;
        rim1 = h*round(lines*rim1)/lines;

        //if(ri < 0.)
        {
            dlinesegment(vec2(xix, x.y), vec2(-.5*size, mix(ri,rim1,.5)), vec2(-.25*size, ri), d);
            dlinesegment(vec2(xix, x.y), vec2(-.25*size, ri), vec2(.25*size, ri), da);
            d = min(d, da);
            dlinesegment(vec2(xix, x.y), vec2(.25*size, ri), vec2(.5*size, mix(ri,rip1,.5)), da);
            d = min(d, da);
            stroke(d, .002+.002*op, d);
            col = mix(col, op*(1.-n)*vec3(0.44,1.00,1.00), sm(d));

            // Dots
            lfnoise(8.*xixj*c.xx-3.*iTime*c.xy+14.*i, n);
            n = .5+.5*n;
            d = length(vec2(xix, x.y-ri))-mix(.0,.35,n)/vs;
            c1 = mix(vec3(0.27,0.63,0.87), vec3(0.49,0.97,1.00), n);
            col = mix(col, op*(1.-n)*c1, sm(d));
            stroke(d - .009, (1.-n)*.005, d);
            c1 *= 2.4;
            col = mix(col, op*(1.-n)*c1, sm(d));
        }
        
        x -= -.1+.2*op;
    }
    
    //mix to blackish
    lfnoise(3.*x.xy-vec2(1.,.1)*iTime, n);
    stroke(n, .3, n);
    col = mix(col, c.yyy, n);
    col = mix(col, .1*col, 1.-s);
    
    col = mix(col, mix(col, vec3(0.27,0.63,0.87), mix(.4,.8,.5+.5*x.y/.1)), sm(abs(x.y)-.1));
    col = mix(col, c.xxx, sm(abs(abs(x.y)-.11)-.001));
    
    col = mix(col, col*col, clamp(-x.y/.1,0.,1.));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;
    
    colorize(1.5*uv-12.*c.xy/*-.3*iTime*c.xy*/, col);

    fragColor = vec4(clamp(col,0.,1.),1.0);
}
