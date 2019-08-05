/* Fuer Elite - 64k Intro by Team210 at Underground Conference 9
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#version 130

uniform float iTime, iProgress;
uniform vec2 iResolution;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

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

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

float sm(float d)
{
    return step(d,0.);
//    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.0*pi < philo.x && p+2.0*pi > philo.y) || (p-2.0*pi < philo.x && p-2.0*pi > philo.y))
        d = abs(length(x)-R);
    else d = min(
        length(x-vec2(cos(p0), sin(p0))),
        length(x-vec2(cos(p1), sin(p1)))
        );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    
    // background
    vec3 col = c.yyy,
        o = c.yzx,
        r = c.xyy, 
        u = cross(r,normalize(o)),
        t = uv.x * r + uv.y * u,
        dir = normalize(t-o),
        size,
        x,
        n;
    
    // Floor
    float d = -(o.z-.35)/dir.z,
        da;
    vec2 y = (o + d * dir).xy;
    d = length(y)-.3;
    stroke(d, .04, d);
    float p = atan(y.y, y.x);
    p = mod(p, pi/8.)-.5*pi/8.;
    d = mix(d, 1., sm(p));
    col = mix(col, vec3(0.81,0.15,0.18), sm(d));
    stroke(d, .005, da);
    col = mix(col, vec3(0.04,0.18,0.24), sm(da));
    stroke(d-.01, .0025, d);
    col = mix(col, vec3(0.81,0.15,0.18), sm(d));
    
    // Layer above
    d = -(o.z-.3)/dir.z;
    y = (o + d * dir).xy;
    d = length(y)-.3;
    stroke(d, .04, d);
    p = atan(y.y, y.x);
    p = mod(p, pi/8.)-.5*pi/8.;
    d = mix(d, 1., sm(p));
    stroke(d, .0035, d);
    col = mix(col, vec3(0.85,0.87,0.89), sm(d));
    
    // Progress
    d = -(o.z-.3)/dir.z;
    y = (o + d * dir).xy;
    dcirclesegment(y, .3, -pi, mix(-pi, pi, iProgress), d);
    stroke(d, .03, d);
    col = mix(col, mix(col, vec3(0.85,0.87,0.89), .8), sm(d));
    stroke(d, .0035, d);
    col = mix(col, mix(col, vec3(0.04,0.18,0.24), .5), sm(d));
    
    // dots above
    for(float i=-.1; i>-.8; i-= .01)
    {
        d = -(o.z-.25-i)/dir.z;
        y = (o + d * dir).xy;
        float n;
        p = atan(y.y, y.x);
        lfnoise(vec2(i,p), n);
        d = length(y)-.3+n;
        stroke(d, .04, d);
        
        float pp = pi/(32.*(-i));
        p = mod(p, pp)-.5*pp;
        d = mix(d, 1., sm(p));
        stroke(d, .017, d);
        col = mix(col, mix(col, mix(col,vec3(0.81,0.15,0.18),step(atan(y.y,y.x)-p,mix(-pi,pi,iProgress))) , .2), sm(d));
        stroke(d, .005, d);
        col = mix(col, mix(col, mix(col,vec3(0.04,0.18,0.24),step(atan(y.y,y.x)-p,mix(-pi,pi,iProgress))) , .2), sm(d));
        
    }
    
    vec2 ff = abs(mod(uv,.1)-.05)-.0005;
    col = mix(col, vec3(.7,.11,.42), sm(min(ff.x,ff.y)));
    ff = abs(mod(uv-.0005,.1)-.05)-.0005;
    col = mix(col, c.xxx, sm(min(ff.x,ff.y)));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
