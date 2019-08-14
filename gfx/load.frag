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

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

void palette1(in float scale, out vec3 col)
{
    const int N = 5;
   
    //*
    const vec3 colors[N] = vec3[N](
            vec3(0.82,0.27,0.13),
            vec3(0.85,0.77,0.68),
            vec3(0.65,0.59,0.55),
            vec3(0.45,0.29,0.24),
            vec3(0.85,0.27,0.15)
        );
    //*/
    
    /*
	const vec3 colors[N] = vec3[N](
       	vec3(0.86,0.21,0.13),
        vec3(0.85,0.80,0.62),
        vec3(0.22,0.25,0.25),
        vec3(0.16,0.17,0.17),
        vec3(0.12,0.12,0.13)
    );
    //*/
    
	/*
    const vec3 colors[N] = vec3[N](
       	vec3(0.00,0.00,0.00),
        vec3(0.64,0.05,0.05),
        vec3(0.91,0.06,0.05),
        vec3(0.96,0.82,0.65),
        vec3(0.65,0.49,0.36)
    );
    //*/
	float index = floor(scale*float(N)), 
        remainder = scale*float(N)-index;
    col = mix(colors[int(index)],colors[int(index)+1], remainder);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = .1*c.xxx;
    col = mix(col, .3*c.xxx, sm(-min(abs(mod(uv.x,.01)-.005)-.003,abs(mod(uv.y,.01)-.005)-.003)));
    
    float d, d0;
    dlinesegment(uv, -.25*a*c.xy, .25*a*c.xy, d);
 	d = abs(d)-.05;
    col = mix(col, .3*c.xxx, sm(d));
    d = abs(d)-.001;
    col = mix(col, c.xxx, sm(d));
    dlinesegment(uv, -.25*a*c.xy, -.25*a*c.xy+.5*a*iProgress*c.xy, d);
    d = abs(d)-.045;
    
    vec3 c1;
    palette1(clamp((uv.y+.045)/.15,0.,1.), c1);
    
    col = mix(col, c1, sm(d));
    
    //col += col;
    //col *= col;
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
