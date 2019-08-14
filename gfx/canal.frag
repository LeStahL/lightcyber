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

void scale(out float s);
void rand(in vec2 x, out float n);
void hash31(in float p, out vec3 d);
void lfnoise(in vec2 t, out float n);
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
void smoothmin(in float a, in float b, in float k, out float dst);
void dsmoothvoronoi(in vec2 x, out float d, out vec2 z);

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

void normal(in vec3 x, out vec3 n, in float dx);

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, 
        s;
    
    scale(iScale);
    
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
        normal(x,n, 5.e-4);
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
                normal(x,n, 5.e-4);
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
    
    col *=mix(1.1,2.6,mix(.5+.5*nn,1.,0.*iScale));
    col *= col;
    col = clamp(col, 0., 1.);
    if(col == c.xxx) col = c.yyy;
    
//     col = mix(col,c.yyy, .1);
    
//     col = mix(col, mix(col, length(col)/sqrt(3.)*c.xxx, .7), iScale);
//     col = mix(col, 1.7*1.7*col*col, iScale);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
