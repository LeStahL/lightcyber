/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
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

uniform float iFSAA;
uniform vec2 iResolution;
uniform sampler2D iChannel0;
uniform float iTime;

out vec4 gl_FragColor;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

float nbeats;
float iScale;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord_ )
{
    vec2 fragCoord = fragCoord_;
    float a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    
    vec3 col = vec3(0.);
    float delta = 0.;
    vec2 n = c.yy;
    
    if(iTime > 66.27 && iTime < 82.76) // blend in to caleidoscope
    {
        n = vec2(7.,1.);
        float phi = abs(mod(atan(uv.y, uv.x),pi/n.x)-.5*pi/n.x);
        uv = length(uv)*vec2(cos(phi+2.*pi), sin(phi+2.*pi));
        fragCoord = mix(fragCoord, (uv + .5*vec2(a,1.))*iResolution.yy, smoothstep(66.27,67.27,iTime)*(1.-smoothstep(81.76,82.76,iTime)));
    }
    else if(iTime > 120.0 && iTime < 136.0) // blend in to caleidoscope
    {
        n = vec2(7.,1.);
        float phi = abs(mod(atan(uv.y, uv.x),pi/n.x)-.5*pi/n.x);
        uv = length(uv)*vec2(cos(phi+2.*pi), sin(phi+2.*pi));
        fragCoord = mix(fragCoord, (uv + .5*vec2(a,1.))*iResolution.yy, smoothstep(120.,121.,iTime)*(1.-smoothstep(135.,136.,iTime)));
    }
    
    float bound = sqrt(iFSAA)-1.;
//     bound = mix(bound, 4., smoothstep(66.27,67.27,iTime)*(1.-smoothstep(81.76,82.76,iTime)));
    
    if((iTime > 66.27 && iTime < 82.76) || (iTime > 120.0 && iTime < 136.0))
    {
        for(float i = -.5*bound; i<=.5*bound; i+=1.)
            for(float j=-.5*bound; j<=.5*bound; j+=1.)
            {
                col += texture(iChannel0, fragCoord/iResolution.xy+vec2(i,j)*3./max(bound, 1.)/iResolution.xy).xyz;
            }
        col /= iFSAA;
    }
    else
    {
        for(float i = -.5*bound; i<=.5*bound; i+=1.)
            for(float j=-.5*bound; j<=.5*bound; j+=1.)
            {
                col += texture(iChannel0, fragCoord/iResolution.xy+vec2(i,j)*mix(3.,15.,2.*abs(fragCoord.y/iResolution.y-.5))*exp(-abs(1.e-2*length(fragCoord.xy)/iResolution.y-.5))/max(bound, 1.)/iResolution.xy).xyz;
            }
        col /= iFSAA;
    }
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
