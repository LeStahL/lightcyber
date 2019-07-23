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


float nbeats;
float iScale;
float smoothdis;

out vec4 gl_FragColor;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;


void hsv2rgb(in vec3 hsv, out vec3 rgb);
void rgb2hsv(in vec3 rgb, out vec3 hsv);
void rand(in vec2 x, out float n);
void lfnoise(in vec2 t, out float n);
void rot3(in vec3 p, out mat3 rot);
void dbox(in vec2 x, in vec2 b, out float d);
void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);

vec2 ind;
void scene(in vec3 x, out vec2 sdf)
{
    x.y += mix(1.,1.01,iScale)*iTime;
    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));
    x.xy = R*x.xy;
    
    float d,
        size = .1;
    vec2 x2 = mod(x.xy,size)-.5*size;
	
    ind = (x.xy - x2)/size;
    dbox(x2, .5*size*c.xx, d);
    zextrude(x.z, -d-.005, .05, d);
    d = max(x.z,d);
    d = abs(d);
    sdf = vec2(d,2.);
    
    float r;
    rand(1.e0*ind-1.e2*floor(iTime), r);
    //lfnoise(12.*ind-1.*iTime, r);
    //r = .5+.5*r;
    if(r > .7)
    {
        float r2;
        rand(ind-1337., r2);
        r2 = .5+.5*iScale;
        dbox(x2, .5*size*c.xx, d);
        zextrude(x.z, -d-.02, .3*(r-.7)/.3*r2, d);
        stroke(d, .001, d);
        add(sdf, vec2(d,1.), sdf);
    }
}

void normal(in vec3 x, out vec3 n, in float dx);
float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void colorize(in vec2 x, out vec3 col)
{
    x.y += mix(1.,1.01,iScale)*iTime;
    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));
    x.xy = R*x.xy;
    
    float d,
        size = .1,
        r;
    vec2 x2 = mod(x.xy,size)-.5*size;
    
    rand(1.e0*ind-1.e2*floor(iTime), r);
    //lfnoise(12.*ind-1.*iTime, r);
    //r = .5+.5*r;
    col = mix(.14*c.xxx, .33*c.xxx, r);
    dbox(x2, .35*size*c.xx, d);
    if(r > .9)
    {
        col = mix(col, mix(c.xxy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(c.xyy,c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, c.xyy, sm(d));
    }
	else if(r > .8)
    {
        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(.7*c.xxy,c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, .7*c.xxy, sm(d));
    }
    else if(r > .7)
    {
        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(mix(c.xxy, c.xyy, .5),c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, mix(c.xxy, c.xyy, .5), sm(d));
    }
    
    // Truchet
    /*
    float da = floor(4.*r)*pi/2.;
    R = mat2(cos(da), sin(da),-sin(da),cos(da));
    x2 = R * x2;
    if(r > .3)
    {
    	dspline2(x2,-.5*size*c.xy, c.yy, -.5*size*c.yx, d);
        dspline2(x2,.5*size*c.xy, c.yy, .5*size*c.yx, da);
        
    }
    else
    {
        dlinesegment(x2,-.5*size*c.xy, .5*size*c.xy, d);
        dlinesegment(x2,-.5*size*c.yx, .5*size*c.yx, da);
    }
    d = min(d, da);
    stroke(d, .001, d);
    col = mix(col,  .9*c.xyy, sm(d));
    stroke(d-.004, .002, d);
    col = mix(col, .0*c.xyy, sm(d));
*/
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    nbeats = mod(iTime-30./29./4., 60./29./4.);
    iScale = nbeats-30./29./4.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 5./29., iScale));
    
    vec3 ra = vec3(nbeats-60./29./4., nbeats, nbeats+60./29./4.);
    rand(ra.x*c.xx, ra.x);
    rand(ra.y*c.xx, ra.y);
    rand(ra.z*c.xx, ra.z);
    smoothdis = mix(
                    mix(ra.x,ra.y, smoothstep(0.,.5,(iTime-nbeats)/(60./29./4.))),
                    ra.z,
                    smoothstep(.5,1.,(iTime-nbeats)/(60./29./4.)));
    
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
    int N = 100,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.15)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(x.z<-.05)
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
        normal(x,n, 5.e-4);
        
        if(s.y == 1.)
        {
            vec3 l = normalize(x+c.xzx);
            vec3 c1;
            
            float r,si = 1.;
		    rand(ind-1.e2*floor(iTime), r);
            if(r > .9)
                col = c.xyy;
            else if(r > .8)
            {
                col = .7*c.xxy;
                si = -1.;
            }
            else if(r > .7)
                col = mix(c.xxy, c.xyy, .5);
            
//             vec3 hsv;
//             rgb2hsv(col, hsv);
//             float na;
//             lfnoise(x.xy-si*iTime+4.*hsv.x, na);
//             hsv.x = mod(1.*hsv.x+.2*na+si*iTime, 2.*pi);
//             hsv2rgb(hsv, col);
            
            float sc = clamp((r-.7)/.3,0.,1.);
            col = mix(mix(col, c.xxx, .1*sc), .4*c.xyy, sc);
            col = .3*col
                + .9*col * abs(dot(l,n))
                + 1.3*col * abs(pow(dot(reflect(-l,n),dir),3.));
            col = mix(col, c.xxx, .4);
            col *= col;
            
            d = -(o.z)/dir.z;
            x = o + d * dir;
            scene(x,s);
            l = normalize(x+c.xzx);
            colorize(x.xy, c1);
            n = c.yyx;
            
            c1 = .1*c1
                + .8*c1 * abs(dot(l,n))
                + c1 * abs(pow(dot(reflect(-l,n),dir),3.));
            col = mix(col, c1, .3);
        }
        else if(s.y == 2.)
        {
            vec3 l = normalize(x+c.xzx);
            float r;
            
            colorize(x.xy, col);
            col = .1*col
                + .8*col * abs(dot(l,n))
                + col * abs(pow(dot(reflect(-l,n),dir),3.));
                
//             vec3 hsv;
//             rgb2hsv(col, hsv);
//             float na;
//             lfnoise(x.xy-iTime+4.*hsv.x, na);
//             hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);
//             hsv2rgb(hsv, col);
        }
    }
    col += col;
    col *= col;
    col *= mix(c.xxx, col, iScale);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
