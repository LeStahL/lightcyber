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
void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);
void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);
void dbox(in vec2 x, in vec2 b, out float d);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void stroke(in float d0, in float s, out float d);
void dvoronoi(in vec2 x, out float d, out vec2 z);
void rot3(in vec3 phi, out mat3 R); 
void scale(out float s);

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

mat3 R;
void scene(in vec3 x, out vec2 sdf)
{
    x = R * x;
    x.x += .3*iTime;
    x *= 2.;
    
    vec3 n;
    lfnoise(x.x*c.xx-iTime, n.x);
    lfnoise(2.*x.x*c.xx-iTime-1337., n.y);
    lfnoise(x.x*c.xx+2.*iTime-2337., n.z);

    x.yz += .1*vec2(cos(x.x), sin(x.x))*n.xy;
    
    mat3 RR;
    rot3(1.3*mix(.2,1.5, .5+.5*n.x)*n.z * c.xyy, RR);
    x = RR * x;
    
//     x = R * x;
    
    x.z = abs(x.z);
    
    float d, da, db;
    
    graf(x.xy, d);
    stroke(d+mix(.01,.04, iScale), mix(.01,.04, iScale), da);
    //stroke(d,.01,da);
    
    float v;
    vec2 ind;
    dvoronoi(12.*x.xy, v, ind); // 12.
    
    zextrude(x.z, -d, .1-.1*v, d);
    
	sdf = vec2(d,1.);
    float modsize = .025,
		y = mod(d-.3-.02*iTime,modsize)-.5*modsize,
        yi = (d-y)/modsize;
    
    float na;
    lfnoise(2.*yi*c.xx-.3*iTime, na);

    zextrude(x.z-.05*na, -y, mix(0.,.05+.05*na,iScale), d);
    stroke(d,.035,d);
        
    
    
    zextrude(x.z, -da, .25, da);
    
	add(sdf, vec2(da, 1.), sdf);
	
	lfnoise(5.*x.xy, da);
	    mfnoise(x.xy, 32., 422., .45, db);
        da = .5*(db+da);
		sdf.x -= .001*da; // .001
        stroke(da, .1, da);
        sdf.x -=.005*da; // .005
    add(sdf, vec2(d, 1.), sdf);
    add(sdf, vec2(x.z+.25,1.), sdf);
    
    // xa(t)
    float xa = mix(x.x+3.*a,x.x-3.*a,clamp(iTime/3.,0.,1.));
    xa = mix(xa,-xa, clamp((iTime-3.)/6., 0.,1.));
    sdf.x += mix(0., 2., smoothstep(-.5, .5, xa));
    
}

void normal(in vec3 x, out vec3 n, in float dx);

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void colorize(in vec2 x, out vec3 col)
{
    x.x += .3*iTime;
    x *= 2.;

    float n;
    lfnoise(x.x*c.xx-iTime, n);
    x.y += .3*cos(x.x)*n;
 
    float d;
    graf(x, d);
    col = mix(col, mix(mix(vec3(0.85,0.87,0.89), c.xxx, step(50., iTime)), mix(vec3(0.04,0.18,0.24),vec3(0.00,0.20,0.36),step(50.,iTime)), clamp(abs(x.y/2.),0.,1.)), sm(d-.2));
    col = mix(col, mix(vec3(1.00,0.40,0.39), vec3(0.00,0.67,0.91), step(50., iTime)), sm(d));
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
		col = mix(col, vec3(0.27,0.27,0.27), sm(d));
        stroke(d, .1, d);
        col = mix(col, 1.5*col, sm(d));
    }
    
    col *= mix(1., 1.6, iScale);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    scale(iScale);
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    
    if(iTime > 71.)
        rot3(mix(pi/4.*c.xyy, 7.*pi/4.*c.xxy,smoothstep(71.,86.,iTime)), R);
    else R = mat3(1.);
    
    float sc2 = 0.,
        sc3 = 0.;
        
    vec3 col = c.yyy, 
        o = mix(1.,.5,smoothstep(0.,5.,clamp(iTime-71.,0.,5.)))*c.yzx,
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

    vec3 c1;
    float d = -(o.z-.35)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(x.z<-.15)
        {
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
            colorize(x.xy, c1);
            c1 = .1*c1
                + 1.*c1 * abs(dot(l,n))
                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));
        }
        else if(s.y == 2.)
        {
            vec3 l = normalize(x+c.xzx);
            float r;
            lfnoise(x.xy, r);
            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),sin(2.*iScale*r*x));
            c1 = .1*c1
                + .8*c1 * abs(dot(l,n))
                + 6.5*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
        col = c1;
    }
    
    
    //col += col;
    
    col *= col*col;
    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    
    // Background geraffel
    if(length(col) < .001)
    {
        float v, ra, v2;
        vec2 ind, ind2;
        lfnoise(iTime*c.xx, ind2.x);
        lfnoise(iTime*c.xx-1337., ind2.y);
        dvoronoi(12.*(uv-.03*ind2), v, ind);
//         dvoronoi(43.*uv, v2, ind2);
//         ind += .1*ind2;
        rand(ind, ra);
        stroke(-v, .05, v);
        v = -v;
        col = mix(col, .3*ra*mix( .5*vec3(1.00,0.40,0.39), .05*c.xxx, clamp(tanh(1.5*length(uv)),0.,1.)), sm(v));
        col *= mix(13.,1.,smoothstep(0.,.5,clamp((iTime-6.),0.,1.)));
    }
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
