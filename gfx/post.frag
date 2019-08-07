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

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = mix(sda, sdb, step(sdb.x, sda.x));
}

mat3 R;
void scene(in vec3 x, out vec2 sdf)
{
    float d;
    
	// Big red box    
    dbox3(x, .2*c.xxx, sdf.x);
    sdf.y = 1.;
    
    // Holes
    
    // 2 upper bar
    dbox3(x-.1*c.xyy, vec3(.02,.3,.12), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 2 right bar
    dbox3(x-.05*c.xyy-.1*c.yyx, vec3(.07,.3,.02), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 2 mid bar
    dbox3(x, vec3(.02,.3,.1), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 2 left bar
    dbox3(x+.05*c.xyy+.1*c.yyx, vec3(.07,.3,.02), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 2 dot
    dbox3(x+.1*c.xyy-.1*c.yyx, vec3(.02,.3,.02), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 1 bar
    dbox3(x+.04*c.yyx, vec3(.3,.02,.08), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 1 dot
    dbox3(x-.1*c.yyx, vec3(.3,.02,.02), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
    
    // 0 big stripes
    vec3 y = vec3(x.x, abs(x.y), x.z);
    dbox3(y-.05*c.yxy, vec3(.1,.03,.3), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));

	// 0 small stripes
    dbox3(y-.1*c.yxy-.06*c.xyy, vec3(.08,.021,.3), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));

    // 0 upper/lower stripes
    vec3 z = vec3(abs(x.x), x.yz);
	dbox3(z-.119*c.xyy, vec3(.021,.08,.3), d);
    sdf.x = max(-d, sdf.x);
    sdf.y = mix(sdf.y, 2., step(d, sdf.x));
}

void normal(in vec3 x, out vec3 n, in float dx)
{
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord_ )
{
    vec2 fragCoord = fragCoord_;
    float a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    
    vec3 col = vec3(0.);
    float delta = 0.;
//     vec2 n = c.yy;
    
    // Box
    rot3(vec3(-2.*pi/8.,2.*pi/8.,2.*pi/4.)-iTime*vec3(1.1,1.3,1.5), R);
    
    float d;
    vec2 s;
    vec3 o, r, u, t, size, dir, x, n;
    vec2 uv2 = 10.*(uv-vec2(-.45*a,.45));
    o = R * c.yyx;
	r = c.xyy; 
	u = c.yxy;
	t = c.yyy; 
    int N = 250,
        i;
    t = uv2.x * r + uv2.y * u;
    t = R * t;
    dir = normalize(t-o);

    size = .2*c.xxx;

	vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes
    vec2 abxlo = abs(o.yz + tlo.x*dir.yz),
        abylo = abs(o.xz + tlo.y*dir.xz),
        abzlo = abs(o.xy + tlo.z*dir.xy);

    vec4 dn = 100.*c.xyyy;

    dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));
    dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));
    dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));
    
    d = dn.r;
    if(d<=2.)
    {
        x = o + d * dir;
        scene(x,s);
        
        if(s.x > 1.e-4)
        {
            for(i = 0; i<N; ++i)
            {
                x = o + d * dir;
                scene(x,s);
                if(s.x < 1.e-4)break;
                d += s.x;
            }
        }
        
        if(i<N)
        {
            normal(x,n, 5.e-4);
            
            if(s.y == 1.)
            {
                vec3 l = normalize(x+c.zzx*vec3(1.3,.9,1.2));
                col = vec3(0.81,0.15,0.18);
                col = .3*col
                    + .4*col * abs(dot(l,n))
                    + .6 * col * abs(pow(dot(reflect(-l,n),dir),2.));
            }
            else if(s.y == 2.)
            {
                vec3 l = normalize(x+c.zzx*vec3(1.3,.9,1.2));
                col = .7*c.xxx;
                col = .5*col
                    + .4*col * abs(dot(l,n))
                    + .8 * col * abs(pow(dot(reflect(-l,n),dir),2.));
            }
        }
    }
    else
    {
    /*
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
    }*/
    
        float bound = sqrt(iFSAA)-1.;
    //     bound = mix(bound, 4., smoothstep(66.27,67.27,iTime)*(1.-smoothstep(81.76,82.76,iTime)));
        
        //if((iTime > 66.27 && iTime < 82.76) || (iTime > 120.0 && iTime < 136.0))
        if(iTime > 0.)
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
    }
    
    // Scan lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
