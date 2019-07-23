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
 
 uniform float iTime;
 uniform float iMaxTime;
 uniform vec2 iMouse;
 uniform float iPlaying;
 uniform vec2 iResolution;
 uniform sampler2D iChannel0;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

// compute distance to regular triangle
void dtriangle2(in vec2 uv, in float r, out float d)
{
    float dp = 2.*pi/3.;
    vec2 p0 = r*vec2(cos(pi/2.), -sin(pi/2.)),
        p1 = r*vec2(cos(pi/2.+dp), -sin(pi/2.+dp)),
        p2 = r*vec2(cos(pi/2.+2.*dp), -sin(pi/2.+2.*dp)), 
        pd = p2-p1;
    
    d = min(dot(uv-p0,c.xz*(p1-p0).yx),dot(uv-p1, pd.yx*c.xz));
	d = min(d, dot(uv-p2, (p0-p2).yx*c.xz))/length(pd);
}

// 2D box
void dbox(in vec2 x, in vec2 b, out float d)
{
	vec2 da = abs(x)-b;
	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

void dpause2(in vec2 uv, in float r, out float d)
{
    dbox(vec2(abs(uv.x)-r/4.,uv.y), vec2(.15*r,.4*r), d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = texture(iChannel0, fragCoord/iResolution.xy).rgb, ccol;
    float d = abs(uv.y+.45)-.05;
    
    // Play/pause controls
    ccol = mix(col,c.yyy,.5*sm(d));
    vec3 c1 = 2.*c.xxx;
    if(iPlaying == 0.) // Check if mouse is in play sign
    {
        dtriangle2(((iMouse.xy/iResolution.yy-.5*vec2(a,1.)-vec2(-.44*a,-.455))).yx*c.xz, .03, d);
    	c1 = mix(ccol, mix(c.xxx, c1, sm(-d)), .5);
        dtriangle2(((uv-vec2(-.44*a,-.455))).yx*c.xz, .03, d);
	    ccol = mix(ccol, c1, sm(-d));
    }
    else // Check if mouse is in pause sign
    {
        dpause2(((iMouse.xy/iResolution.yy-.5*vec2(a,1.)-vec2(-.44*a,-.455))), .05, d);
    	c1 = mix(ccol, mix(c.xxx, c1, sm(d)), .5);
        dpause2(((uv-vec2(-.44*a,-.455))), .05, d);
	    ccol = mix(ccol, c1, sm(d));
    }
    
    // Progress bar
    vec2 xpr = mix(vec2(-.3*a,-.45), vec2(.3*a,-.45), clamp(iTime/iMaxTime, 0., 1.));
    dlinesegment(iMouse.xy/iResolution.yy-.5*vec2(a,1.), vec2(-.3*a,-.45), vec2(.3*a, -.45), d);
    stroke(d, .0025, d);
    d = min(d, length(iMouse.xy/iResolution.yy-.5*vec2(a,1.)-xpr)-.015);
    c1 = mix(ccol, mix(c.xxx, 2.*c.xxx, sm(d)), .5);
    dlinesegment(uv, vec2(-.3*a,-.45), vec2(.3*a, -.45), d);
    stroke(d, .0025, d);
    d = min(d, length(uv-xpr)-.015);
    ccol = mix(ccol, c1, sm(d));

    
    col = mix(ccol, col, smoothstep(.1,.2,iMouse.y/iResolution.y));
    
    dlinesegment(uv, vec2(-.9*a,-.9), vec2(.9*a, -.9), d);
    stroke(d, .01, d);
	
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
