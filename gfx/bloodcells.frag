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
void stroke(in float d0, in float s, out float d);
void zextrude(in float z, in float d2d, in float h, out float d);
float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}
void smoothmin(in float a, in float b, in float k, out float dst);
void dsmoothvoronoi(in vec2 x, out float d, out vec2 z)
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
            smoothmin(ret, d, .4+.38*n, ret);
        }
    
    d = ret;
    z = pf;
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf);
vec2 ind;
void scene(in vec3 x, out vec2 sdf)
{    
    x.y += .3*iTime;
    float d;
    dsmoothvoronoi(mix(2.,3.,smoothstep(10.,12.,iTime))*x.xy-1337.,d,ind);
    stroke(d, .1, d);
    float modsize = .04,
		y = mod(d-.02*iTime,modsize)-.5*modsize,
        yi = (d-y)/modsize;
    
    float n;
    lfnoise(2.*yi*c.xx-.3*iTime, n);
    
    zextrude(x.z-.05*n, -y, mix(0.,.05+.05*n,iScale), d);
    
    stroke(d,mix(0.,.02,iScale),d);
    
    sdf = vec2(d, 2.);
    
    add(sdf, vec2(x.z+.05,1.), sdf);
}   

void normal(in vec3 x, out vec3 n, in float dx);

void colorize(in vec2 x, out vec3 col)
{
    col = .5*c.xxx;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    nbeats = mod(iTime, 60./29.);
    iScale = nbeats-30./29.;
    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));
    
    iScale *= (1.-smoothstep(51., 52., iTime));
    
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
    int N = 200,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.1)/dir.z;
    
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
        d += s.x<1.e-2?min(s.x,5.e-4):s.x;
        //d += s.x<5.e-2?min(s.x,2.e-3):s.x;
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
            lfnoise(x.xy-iTime, r);
            col = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*r*x));
            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*r*x));
            col = mix(col, c1, .5+.5*r);
            col = .1*col
                + .8*col * abs(dot(l,n))
                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));
        }
    }
    
    col *= col*col;
    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    
    col *= mix(col, length(col)/sqrt(3.)*c.xxx, iScale);

    col = mix(c.yyy, col, smoothstep(0., 1., iTime));
    col = mix(col, c.yyy, smoothstep(53.69, 54.69, iTime));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
