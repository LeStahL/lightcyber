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

float nbeats;
float iScale;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void hsv2rgb(in vec3 hsv, out vec3 rgb);
void rgb2hsv(in vec3 rgb, out vec3 hsv);
void rand(in vec2 x, out float n);
void lfnoise(in vec2 t, out float n);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t);
void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}
void smoothmin(in float a, in float b, in float k, out float dst);
float ind;
void graf(in vec2 x, in vec2 size, in float w, out float dst)
{
    x.x += .3*iTime;
    x.y -= .1;
    
    vec2 y = vec2(mod(x.x, size.x)-.5*size.x, x.y);
    ind = (y.x-x.x)/size.x;
    vec3 r,
        d;
    
    size.x -= 2.*w; 
    
    rand(vec2(ind,0.), r.x);
    rand(vec2(ind,1337.), r.y);
    rand(vec2(ind,2337.), r.z);
    
    vec2 pc = vec2(.5*size.x*(-1.+2.*r.y), 0.);
    
    // x component of r selects first actor
    if(r.x < .5) // Lin
    {
        vec2 p1 = vec2(.5*size.x*(-1.+4.*r.x),-.5*size.y),
            p2 = pc;
        dlinesegment(y, p1, p2, dst);
    }
    else // Quad
    {
        vec2 p1 = vec2(.5*size.x*(-1.+4.*(r.x-.5)),-.5*size.y),
            p2 = vec2(.5*size.x*(-1.+2.*r.y), -.5*size.y),
            p3 = pc;
        dspline2(y, p1, p2, p3, dst);
    }
    
    // y component of r selects connection type
    
    // z component of r selects second actor
    if(r.z < .5) // Lin
    {
        vec2 p1 = vec2(.5*size.x*(-1.+4.*r.z),.5*size.y),
            p2 = pc;
        dlinesegment(y, p1, p2, d.x);
    }
    else // Quad
    {
        vec2 p1 = vec2(.5*size.x*(-1.+4.*(r.z-.5)),.5*size.y),
            p2 = vec2(.5*size.x*(-1.+2.*r.y), .5*size.y),
            p3 = pc;
        dspline2(y, p1, p2, p3, d.x);
    }
    
    dst = min(dst, d.x);
    
    // Generate displacement
    lfnoise(12.*x, d.y);
    lfnoise(22.*x, d.z);
    d.y += .3*d.z;
    d.y = floor(.2+d.y);
    
    stroke(dst, w+.08*d.y, dst);
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
void scene(in vec3 x, out vec2 sdf)
{
    float n;
    lfnoise(5.*x.xy+iTime*c.xy, n);
    x.z += .03*n;
    
    float d = 1., da;
    for(float i=0.; i<2.; i+=1.)
    {
    	graf(x.xy-1337.333*i*c.xy, vec2(.4,.6), .05, da);
        //d = min(d, da);
        float rr;
        rand(i*c.xx*1.e2,rr); 
        zextrude(x.z, -d+.5*abs(x.z), .07+.06*rr, d);
        stroke(d,mix(.02,.04,iScale), d);
        smoothmin(d, da, .2, d);
    }
    d = mix(1., d, smoothstep(19., 23., iTime));
    sdf = vec2(d, 2.);
    
    add(sdf, vec2(x.z+.05,1.), sdf);
}   

void normal(in vec3 x, out vec3 n, in float dx);
void colorize(in vec2 x, out vec3 col)
{
    x.x += .3*iTime;
    x.y -= .05;
    
    col = .5*c.xxx;
    
    float s = .1;
    vec2 dd = mod(x, s)-.5*s;
    stroke(dd.x, .005, dd.x);
    stroke(dd.y, .005, dd.y);
    col = mix(col, c.xxx, sm(min(dd.x, dd.y)));
    
    float d = 1., da;
    for(float i=0.; i<2.; i+=1.)
    {
    	graf(x.xy-1337.333*i*c.xy-.3*iTime*c.xy, vec2(.4,.6), .05, da);
        //d = min(d, da);
        float rr;
        rand(i*c.xx*1.e2,rr); 
        stroke(da,.18, da);
        smoothmin(d, da, .2, d);
    }
    d = mix(1., d, smoothstep(15., 19., iTime));
    vec3 c1 = vec3(.78,.61*abs(2.*x.y),.15);
    
//     vec3 hsv;
//     rgb2hsv(c1, hsv);
//     float na;
//     lfnoise(x.xy-iTime+4.*hsv.x, na);
//     hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);
//     hsv2rgb(hsv, c1);
    
    col = mix(col, c1, sm(d));
    
    if(d != 1.)
    {
        stroke(d-.03, .03, d);
        col = mix(col, c.yyy, sm(d));
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    vec3 col = c.yyy, 
        o = c.yzx,
        r = c.xyy, 
        u = normalize(c.yxx),
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 150,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.08)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(x.z<-.1)
        {
            col = .2*c.xxx;
            i = N;
            break;
        }
        d += min(s.x,5.e-3);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n, 5.e-3);
        
        if(s.y == 1.)
        {
            vec3 l = normalize(x+.5*c.yzx);
            colorize(x.xy, col);
            col = .1*col
                + 1.*col * abs(dot(l,n))
                + 1.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));
            
        }
        else if(s.y == 2.)
        {
            vec3 l = normalize(x+c.xzx);
            float r;
            lfnoise(x.xy, r);
            col = vec3(0.99,0.43,0.15);
            
            vec3 hsv;
//             rgb2hsv(col, hsv);
//             float na;
//             lfnoise(x.xy-iTime+4.*hsv.x, na);
//             hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);
//             hsv2rgb(hsv, col);
            
            vec3 c2 = vec3(0.44,0.07,0.66);
//             rgb2hsv(c2, hsv);
//             lfnoise(x.xy+iTime+4.*hsv.x, na);
//             hsv.x = mod(1.*hsv.x+.2*na-iTime, 2.*pi);
//             hsv2rgb(hsv, c2);
            
            col = mix(col,c2,sin(2.*iScale*r*x));
            col = .1*col
                + .8*col * abs(dot(l,n))
                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
    }
    
    //col += col;
    
    col *= col;
    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    
    col = mix(c.yyy, col, smoothstep(0., 1., iTime));
    col = mix(col, c.yyy, smoothstep(48.655, 49.655, iTime));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
