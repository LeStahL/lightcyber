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

// void dpolygon(in vec2 x, in float N, out float d);
// void rot(in float phi, out mat2 m);
// void dcircle(in vec2 x, out float d);
// void dbox(in vec2 x, in vec2 b, out float d);
// void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
// void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);

void rand(in vec2 x, out float n);
void hash31(in float p, out vec3 d);
void lfnoise(in vec2 t, out float n);
void dbox(in vec2 x, in vec2 b, out float d);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
void smoothmin(in float a, in float b, in float k, out float dst);
void dbox3(in vec3 x, in vec3 b, out float d);
void dschnappsgirls(in vec2 x, out float d);
void dspacepigs(in vec2 x, out float d);
void dkewlers(in vec2 x, out float d);
void dfarbrausch(in vec2 x, out float d);
void dhaujobb(in vec2 x, out float d);
void dmercury(in vec2 x, out float d);

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
    smoothmin(sdf.x, d, .002, sdf.x);
    indc = vec2(xi, zi);

    // Logos
    tsize = .25;
    float tw = .0005;
    dz = mod(x.z-.5*tsize, tsize)-.5*tsize;
    zi = round((x.z-dz)/tsize);
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

void normal(in vec3 x, out vec3 n, in float dx);

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
        normal(x,n, 5.e-4);
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
            
            normal(x,n, 5.e-4);
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
