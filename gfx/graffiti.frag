/* Lightcyber by Team210 - 64k intro by Team210 at Solskogen 2k19
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

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float iScale, nbeats;

void rand(in vec2 x, out float n);
void lfnoise(in vec2 t, out float n);
void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);
void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);
void dbox(in vec2 x, in vec2 b, out float d);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void stroke(in float d0, in float s, out float d);
void dvoronoi(in vec2 x, out float d, out vec2 z);

void graf(in vec2 x, out float d)
{
    x.y *= .7;
    float size = .4,
        n,
        da;
    vec2 y = vec2(mod(x.x, size)-.5*size, x.y),
        yi = (x-y)/size,
        x1,
        x2,
        x3;
    
    dbox(y,vec2(.75,.75)*size, d);
    
    // lines
    rand(yi, n);
    x1 = vec2(-.5+.02+n*.96, .75)*size,
    x2 = vec2(.5-.02-n*.96, -.75)*size;
    x1.x = floor(5.*x1.x)/5.;
    x2.x = floor(5.*x2.x)/5.;
    x1.x = max(x1.x,-.4*size);
    x1.x = min(x1.x,.4*size);
    x2.x = max(x2.x,-.4*size);
    x2.x = min(x2.x,.4*size);
    dlinesegment(y, x1, x2, da);
    stroke(da, .02, da);
    d = max(d, -da);
    
    // upper triangles
    rand(yi+1337., n);
	x1 = vec2(-.55+n,.75)*size*1.05,
    x2 = vec2(.5-n,.75)*size*1.05,
    x3 = .75*(.8-n)*size*1.05*c.yx-.01*c.yx;
    x1 = round(15.*x1)/15.;
    x2 = round(15.*x2)/15.;
    x3 = round(15.*x3)/15.;
    dtriangle(y, x1, x2, x3, da);
    d = max(d, -da);
    
    // lower triangles
    rand(yi+2337., n);
    x1 = vec2(-.5+n,-.75)*size*1.05,
    x2 = vec2(.55-n,-.75)*size*1.05,
    x3 = -.75*(.8-n)*size*1.05*c.yx+.01*c.yx;
    x1 = round(15.*x1)/15.;
    x2 = round(15.*x2)/15.;
    x3 = round(15.*x3)/15.;
    dtriangle(y, x1, x2, x3, da);
    d = max(d, -da);
}

void zextrude(in float z, in float d2d, in float h, out float d);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);

void scene(in vec3 x, out vec2 sdf)
{
    x.x += .3*iTime;
    
    vec2 n;
    lfnoise(x.x*c.xx-iTime, n.x);
    lfnoise(x.x*c.xx-iTime-1337., n.y);
    x.yz += .3*vec2(cos(x.x), sin(x.x))*n;
    
    float d, da;
    
    graf(x.xy, d);
    stroke(d+mix(.01,.04, iScale), mix(.01,.04, iScale), da);
    //stroke(d,.01,da);
    
    float v;
    vec2 ind;
    dvoronoi(12.*x.xy, v, ind);
    
    zextrude(x.z, -d, .1-.1*v, d);
    
	sdf = vec2(d,1.);
    float modsize = .05,
		y = mod(d-.3-.02*iTime,modsize)-.5*modsize,
        yi = (d-y)/modsize;
    
    float na;
    lfnoise(2.*yi*c.xx-.3*iTime, na);

    zextrude(x.z-.05*na, -y, mix(0.,.05+.05*na,iScale), d);
    //stroke(d,mix(0.,.035,iScale),d);
    
    
    
    add(sdf, vec2(d, 1.), sdf);
    
    zextrude(x.z, -da, .25, da);
	add(sdf, vec2(da, 1.), sdf);
    
    add(sdf, vec2(x.z+.25,1.), sdf);
}

void normal(in vec3 x, out vec3 n, in float dx);

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void colorize(in vec2 x, out vec3 col)
{
    x.x += .3*iTime;

    float n;
    lfnoise(x.x*c.xx-iTime, n);
    x.y += .3*cos(x.x)*n;
    
    float d;
    graf(x, d);
    col = mix(col, vec3(.9,.3,.7), sm(d-.2));
    col = mix(col, vec3(.3,.7,.9), sm(d));
    float da = d;
    stroke(d+mix(.01,.03, iScale), mix(.01,.04,iScale), d);
    //stroke(d,.01, d);
    col = mix(col, 1.4*col, sm(d));
    stroke(d, .001, d);
    col = mix(col, 1.3*col, sm(d));
    
    if(da < .02 && da > -.02)
    {
        lfnoise(5.*x, da);
	    mfnoise(x, 32., 422., .45, d);
        d = .5*(d+da);
		col = mix(col, vec3(.4,.2,.9), sm(d));
        stroke(d, .1, d);
        col = mix(col, 1.5*col, sm(d));
    }
    
    col *= mix(1., 1.6, iScale);
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
    int N = 250,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.35)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(x.z<-.35)
        {
            col = .2*c.xxx;
            i = N;
            break;
        }
        d += min(s.x,5.e-2);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n, 1.e-2);
        
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
            col = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),sin(2.*iScale*r*x));
            col = .1*col
                + .8*col * abs(dot(l,n))
                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
    }
    
    //col += col;
    
    col *= col*col;
    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}