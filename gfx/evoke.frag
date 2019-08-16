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
void dbox(in vec2 x, in vec2 b, out float d);
void stroke(in float d0, in float s, out float d);
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);
void lfnoise(in vec2 t, out float n);
void dvoronoi(in vec2 x, out float d, out vec2 z);
void scale(out float s);

void devoke(in vec2 x, out float d)
{
    x.x += .225;
    x *= 1.1;
    
    // o
    d = length(x+.35*c.xy)-.1;
    stroke(d,.06,d);
    
    // I
    float da;
    dbox(x+.1*c.xy, vec2(.05, .25), da);
    d = min(d, da);
    
    x = 2.*x - vec2(.4,-.2);
    // Mercury
    // Upper part
    dbox(x-.35*c.yx,vec2(.4,.35), da);
    d = min(d,da);
    dbox(x-.7*c.yx, vec2(.2,.2), da);
    d = max(d,-da);
    dbox(x-.25*c.yx,vec2(.2,.05),da);
    d = max(d,-da);
    
    // Lower part
    dbox(x+.1*c.yx,vec2(.1,.2),da);
    d = min(d,da);
    dbox(x+.2*c.yx, vec2(.4,.1),da);
    d = min(d,da);
    
    x = .5*(x + vec2(.4,-.2));
    
    // E
    // Right
    dbox(x-.9*c.xy, vec2(.05, .25), da);
    d = min(d,da);
    
    // Top/bot
    dbox(vec2(x.x-.7, abs(x.y)-.2), vec2(.2, .05), da);
    d = min(d,da);
    
    // Middle
    dbox(x-.7*c.xy, vec2(.2, .05), da);
    d = min(d,da);
    
    // Appendix
    dbox(vec2(x.x-.95,x.y+.2), vec2(.05,.05), da);
    d = min(d,da);
    
    stroke(d,.001,d);
}

void dstripe(in vec2 x, out float d)
{
    dlinesegment(x-a*mix(-.4*c.xy, .4*c.xy, clamp(iTime/6.,0.,1.)), -.5*c.yx, .5*c.yx, d);
    d -= .005;
    float dd;
    vec2 vi;
    dvoronoi(5.*x, dd, vi);
    vi = x-vi/5.;
    dd = abs(length(vi)-.002)-.001;
    d = min(d,dd);
    stroke(d,.001,d);
    d = mix(1.,d, clamp(iTime, 0., 1.));
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    scale(iScale);
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    vec3 col = vec3(0.20,0.01,0.14), 
        o = c.yyx,
        r = c.xyy, 
        u = c.yxy,
        t = c.yyy, 
        dir,
        n,
        x;
    float d, i, ra;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    for(i=1.4; i>=0.; i -= .01)
    {
        lfnoise(102.*i*c.xx-mix(102.,0.,smoothstep(0.,1.,clamp(iTime-6.,0.,1.)))*iTime, ra);
        ra = .5+.5*ra;
        
		d = -(o.z-.2+i)/dir.z;
        x = o + d * dir;
        
        float da;
        dstripe(x.xy, da);
        devoke(x.xy, s.x);
        s.x = mix(da, s.x, smoothstep(0.,1., clamp(iTime-5.5,0.,1.)));
        s.x -= .01*iScale;
        
        if(ra < .5)
        {
            vec3 c1 = mix(mix(vec3(0.75,0.24,0.31), vec3(1.00,0.87,0.57), smoothstep(1.25,1.4,1.4-i)),vec3(0.20,0.01,0.14),i/1.4);
	        if(iTime > 6.)col = mix(col, c1, sm(s.x));
            col = mix(col, mix(1.1,1.,clamp(-iScale+smoothstep(0.,1.,clamp(iTime-6.,0.,1.)),0.,1.))*mix(col,vec3(.7,.45,.3), mix(.02,.1,iScale)), sm(s.x/64.));
        }
    }
    
    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));
    
    // Blend to voronoi background
    vec3 c1 = c.yyy;
    float v, v2;
    vec2 ind, ind2;
    lfnoise((iTime-12.)*c.xx, ind2.x);
    lfnoise((iTime-12.)*c.xx-1337., ind2.y);
    dvoronoi(12.*(uv-.03*ind2), v, ind);
    rand(ind, ra);
    stroke(-v, .05, v);
    v = -v;
    c1 = mix(c1, .3*ra*mix( .5*vec3(1.00,0.40,0.39), .05*c.xxx, clamp(tanh(1.5*length(uv)),0.,1.)), sm(v));
    
    c1 *= mix(1.,13., smoothstep(0.,1., clamp((iTime-11.), 0., 1.)));
    
    col = mix(col, c1, smoothstep(0.,1., clamp((iTime-11.),0.,1.)));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
