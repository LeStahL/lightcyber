//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

extern float progress;int hsv2rgb_handle, rgb2hsv_handle, rand_handle, lfnoise_handle, mfnoise_handle, stroke_handle, add_handle, dvoronoi_handle, normal_handle, rot3_handle, dbox_handle, dlinesegment3_handle, zextrude_handle, dlinesegment_handle, dspline2_handle, smoothmin_handle;
const int nsymbols = 16;
const char *hsv2rgb_source = "#version 130\n\n"
"const float pi = acos(-1.);\n"
"void hsv2rgb(in vec3 hsv, out vec3 rgb)\n"
"{\n"
"    float C = hsv.y * hsv.z,\n"
"        Hprime = hsv.x / pi * 3.,\n"
"        X = C * (1.-abs(mod(Hprime,2.)-1.));\n"
"    \n"
"    if(0. <= Hprime && Hprime <= 1.) rgb = vec3(C, X, 0.);\n"
"    else if( 1. < Hprime && Hprime <= 2.) rgb = vec3(X, C, 0.);\n"
"    else if( 2. < Hprime && Hprime <= 3.) rgb = vec3(0., C, X);\n"
"    else if( 3. < Hprime && Hprime <= 4.) rgb = vec3(0., X, C);\n"
"    else if( 4. < Hprime && Hprime <= 5.) rgb = vec3(X, 0., C);\n"
"    else if( 5. < Hprime && Hprime <= 6.) rgb = vec3(C, 0., X);\n"
"        \n"
"    float m = hsv.z - C;\n"
"    rgb += m;\n"
"}\n"
"\0";
const char *rgb2hsv_source = "#version 130\n\n"
"const float pi = acos(-1.);\n"
"void rgb2hsv(in vec3 rgb, out vec3 hsv)\n"
"{\n"
"    float MAX = max(rgb.r, max(rgb.g, rgb.b)),\n"
"        MIN = min(rgb.r, min(rgb.g, rgb.b)),\n"
"        C = MAX-MIN;\n"
"    \n"
"    if(MAX == MIN) hsv.x = 0.;\n"
"    else if(MAX == rgb.r) hsv.x = pi/3.*(rgb.g-rgb.b)/C;\n"
"    else if(MAX == rgb.g) hsv.x = pi/3.*(2.+(rgb.b-rgb.r)/C);\n"
"    else if(MAX == rgb.b) hsv.x = pi/3.*(4.+(rgb.r-rgb.g)/C);\n"
"    hsv.x = mod(hsv.x, 2.*pi);\n"
"        \n"
"    if(MAX == 0.) hsv.y = 0.;\n"
"    else hsv.y = (MAX-MIN)/MAX;\n"
"        \n"
"    hsv.z = MAX;\n"
"}\n"
"\0";
const char *rand_source = "#version 130\n\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\0";
const char *lfnoise_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
"void lfnoise(in vec2 t, out float n)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    t = smoothstep(c.yy, c.xx, t);\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    n = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\0";
const char *mfnoise_source = "#version 130\n\n"
"// const vec3 c = vec3(1.,0.,-1.);\n"
"void lfnoise(in vec2 x, out float d);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)\n"
"{\n"
"    n = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = d; f<b; f *= 2.)\n"
"    {\n"
"        lfnoise(f*x, buf);\n"
"        n += a*buf;\n"
"        a *= e;\n"
"        nf += 1.;\n"
"    }\n"
"    n *= (1.-e)/(1.-pow(e, nf));\n"
"}\n"
"\0";
const char *stroke_source = "// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\0";
const char *add_source = "void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = mix(sda, sdb, step(sdb.x, sda.x));\n"
"}\n"
"\0";
const char *dvoronoi_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z)\n"
"{\n"
"    vec2 y = floor(x);\n"
"       float ret = 1.;\n"
"    vec2 pf=c.yy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            if(d < df)\n"
"            {\n"
"                df = d;\n"
"                pf = p;\n"
"            }\n"
"        }\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            vec2 o = p - pf;\n"
"            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"            ret = min(ret, d);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\0";
const char *normal_source = "const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"void scene(in vec3 x, out vec2 s);\n"
"void normal(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\0";
const char *rot3_source = "const vec3 c = vec3(1.,0.,-1.);\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\0";
const char *dbox_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\0";
const char *dlinesegment3_source = "#version 130\n\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *zextrude_source = "// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\0";
const char *dlinesegment_source = "#version 130\n\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *dspline2_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"\n"
"//distance to spline with parameter t\n"
"float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum dist3ance to spline\n"
"void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;\n"
"    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist2(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three dist3inct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist2(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist2(p0,p1,p2,x,t.y),\n"
"            dist2(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\0";
const char *smoothmin_source = "// iq's smooth minimum\n"
"void smoothmin(in float a, in float b, in float k, out float dst)\n"
"{\n"
"    float h = max( k-abs(a-b), 0.0 )/k;\n"
"    dst = min( a, b ) - h*h*h*k*(1.0/6.0);\n"
"}\n"
"\0";
const char *voronoidesign_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void hsv2rgb(in vec3 hsv, out vec3 rgb);\n"
"void rgb2hsv(in vec3 rgb, out vec3 hsv);\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.y += .3*iTime;\n"
"    \n"
"    float d, \n"
"        d2,\n"
"        dz,\n"
"        da,\n"
"        r;\n"
"    vec2 p, \n"
"        p2,\n"
"        y;\n"
"        \n"
"    // Net\n"
"    dvoronoi(6.*x.xy, d, p);\n"
"    dvoronoi(12.*x.xy, d2, p2);\n"
"    ind = p;\n"
"    \n"
"    x.z -= -.1+.4*p.x*p2.x/6./12.;\n"
"    \n"
"    // Displacement\n"
"    mfnoise(x.xy-iTime*c.yx, 2., 12., .25, dz);\n"
"    d = x.z-.1-d/12.+.1*dz-.2*d2/12.;\n"
"    //d = mix(d, length(x)-1., .01);\n"
"	sdf = vec2(d, 1.);\n"
"    \n"
"    // Discrete displacement\n"
"    mfnoise(p/6.-iTime*c.yx, 2., 12., .25, dz);\n"
"	add(sdf, vec2(length(x-vec3(p/6.,.14-.1*dz))-.02, 2.), sdf);\n"
"    \n"
"}   \n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    col = vec3(0.09,0.10,0.07);\n"
"    \n"
"    x.y += .3*iTime;\n"
"    \n"
"    float d, \n"
"        d2,\n"
"        dz,\n"
"        da,\n"
"        r;\n"
"    vec2 p, \n"
"        p2,\n"
"        y;\n"
"        \n"
"    // Net\n"
"    dvoronoi(6.*x.xy, d, p);\n"
"    dvoronoi(12.*x.xy, d2, p2);\n"
"    \n"
"    stroke(d, .02, d);\n"
"    col = mix(col, vec3(0.79,0.44,0.36), sm(d));\n"
"    \n"
"    stroke(d2, .02, d2);\n"
"    col = mix(col, vec3(0.89,0.44,0.26), sm(d2));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 5./29., iScale));\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 200,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.2)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.1)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += s.x<5.e-2?min(s.x,2.e-3):s.x;\n"
"        //d += min(s.x,3.e-3);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n,1.e-3);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+.5*c.yzx);\n"
"            colorize(x.xy, col);\n"
"            col = .1*col\n"
"                + 1.8*col * abs(dot(l,n))\n"
"                + 2.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"//             vec3 hsv;\n"
"//             rgb2hsv(col, hsv);\n"
"//             float na;\n"
"//             lfnoise(x.xy-iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, col);\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            rand(ind, r);\n"
"            col = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),iScale*r);\n"
"            col = .1*col\n"
"                + .8*col * abs(dot(l,n))\n"
"                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"                        vec3 hsv;\n"
"//             rgb2hsv(col, hsv);\n"
"//             float na;\n"
"//             lfnoise(x.xy+iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na-iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, col);\n"
"        }\n"
"    }\n"
"    \n"
"    //col += col;\n"
"    \n"
"    col *= col;\n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"//    col *= mix(c.xxx, 2.*c.xxx, iScale);\n"
"    col = mix(col, length(col)/1.732*c.xxx, .5*iScale);\n"
"\n"
"    col = mix(c.yyy, col, smoothstep(0., 1., iTime));\n"
"    col = mix(col, c.yyy, smoothstep(15.55, 16.55, iTime));\n"
"    \n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}	\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *groundboxes_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"float smoothdis;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"\n"
"void hsv2rgb(in vec3 hsv, out vec3 rgb);\n"
"void rgb2hsv(in vec3 rgb, out vec3 hsv);\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"\n"
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.y += mix(1.,1.01,iScale)*iTime;\n"
"    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));\n"
"    x.xy = R*x.xy;\n"
"    \n"
"    float d,\n"
"        size = .1;\n"
"    vec2 x2 = mod(x.xy,size)-.5*size;\n"
"	\n"
"    ind = (x.xy - x2)/size;\n"
"    dbox(x2, .5*size*c.xx, d);\n"
"    zextrude(x.z, -d-.005, .05, d);\n"
"    d = max(x.z,d);\n"
"    d = abs(d);\n"
"    sdf = vec2(d,2.);\n"
"    \n"
"    float r;\n"
"    rand(1.e0*ind-1.e2*floor(iTime), r);\n"
"    //lfnoise(12.*ind-1.*iTime, r);\n"
"    //r = .5+.5*r;\n"
"    if(r > .7)\n"
"    {\n"
"        float r2;\n"
"        rand(ind-1337., r2);\n"
"        r2 = .5+.5*iScale;\n"
"        dbox(x2, .5*size*c.xx, d);\n"
"        zextrude(x.z, -d-.02, .3*(r-.7)/.3*r2, d);\n"
"        stroke(d, .001, d);\n"
"        add(sdf, vec2(d,1.), sdf);\n"
"    }\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    x.y += mix(1.,1.01,iScale)*iTime;\n"
"    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));\n"
"    x.xy = R*x.xy;\n"
"    \n"
"    float d,\n"
"        size = .1,\n"
"        r;\n"
"    vec2 x2 = mod(x.xy,size)-.5*size;\n"
"    \n"
"    rand(1.e0*ind-1.e2*floor(iTime), r);\n"
"    //lfnoise(12.*ind-1.*iTime, r);\n"
"    //r = .5+.5*r;\n"
"    col = mix(.14*c.xxx, .33*c.xxx, r);\n"
"    dbox(x2, .35*size*c.xx, d);\n"
"    if(r > .9)\n"
"    {\n"
"        col = mix(col, mix(c.xxy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(c.xyy,c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, c.xyy, sm(d));\n"
"    }\n"
"	else if(r > .8)\n"
"    {\n"
"        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(.7*c.xxy,c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, .7*c.xxy, sm(d));\n"
"    }\n"
"    else if(r > .7)\n"
"    {\n"
"        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(mix(c.xxy, c.xyy, .5),c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, mix(c.xxy, c.xyy, .5), sm(d));\n"
"    }\n"
"    \n"
"    // Truchet\n"
"    /*\n"
"    float da = floor(4.*r)*pi/2.;\n"
"    R = mat2(cos(da), sin(da),-sin(da),cos(da));\n"
"    x2 = R * x2;\n"
"    if(r > .3)\n"
"    {\n"
"    	dspline2(x2,-.5*size*c.xy, c.yy, -.5*size*c.yx, d);\n"
"        dspline2(x2,.5*size*c.xy, c.yy, .5*size*c.yx, da);\n"
"        \n"
"    }\n"
"    else\n"
"    {\n"
"        dlinesegment(x2,-.5*size*c.xy, .5*size*c.xy, d);\n"
"        dlinesegment(x2,-.5*size*c.yx, .5*size*c.yx, da);\n"
"    }\n"
"    d = min(d, da);\n"
"    stroke(d, .001, d);\n"
"    col = mix(col,  .9*c.xyy, sm(d));\n"
"    stroke(d-.004, .002, d);\n"
"    col = mix(col, .0*c.xyy, sm(d));\n"
"*/\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    nbeats = mod(iTime-30./29./4., 60./29./4.);\n"
"    iScale = nbeats-30./29./4.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 5./29., iScale));\n"
"    \n"
"    vec3 ra = vec3(nbeats-60./29./4., nbeats, nbeats+60./29./4.);\n"
"    rand(ra.x*c.xx, ra.x);\n"
"    rand(ra.y*c.xx, ra.y);\n"
"    rand(ra.z*c.xx, ra.z);\n"
"    smoothdis = mix(\n"
"                    mix(ra.x,ra.y, smoothstep(0.,.5,(iTime-nbeats)/(60./29./4.))),\n"
"                    ra.z,\n"
"                    smoothstep(.5,1.,(iTime-nbeats)/(60./29./4.)));\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx), \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 100,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.15)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.05)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,5.e-3);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            vec3 c1;\n"
"            \n"
"            float r,si = 1.;\n"
"		    rand(ind-1.e2*floor(iTime), r);\n"
"            if(r > .9)\n"
"                col = c.xyy;\n"
"            else if(r > .8)\n"
"            {\n"
"                col = .7*c.xxy;\n"
"                si = -1.;\n"
"            }\n"
"            else if(r > .7)\n"
"                col = mix(c.xxy, c.xyy, .5);\n"
"            \n"
"//             vec3 hsv;\n"
"//             rgb2hsv(col, hsv);\n"
"//             float na;\n"
"//             lfnoise(x.xy-si*iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na+si*iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, col);\n"
"            \n"
"            float sc = clamp((r-.7)/.3,0.,1.);\n"
"            col = mix(mix(col, c.xxx, .1*sc), .4*c.xyy, sc);\n"
"            col = .3*col\n"
"                + .9*col * abs(dot(l,n))\n"
"                + 1.3*col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            col = mix(col, c.xxx, .4);\n"
"            col *= col;\n"
"            \n"
"            d = -(o.z)/dir.z;\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            l = normalize(x+c.xzx);\n"
"            colorize(x.xy, c1);\n"
"            n = c.yyx;\n"
"            \n"
"            c1 = .1*c1\n"
"                + .8*c1 * abs(dot(l,n))\n"
"                + c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            col = mix(col, c1, .3);\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            \n"
"            colorize(x.xy, col);\n"
"            col = .1*col\n"
"                + .8*col * abs(dot(l,n))\n"
"                + col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                \n"
"//             vec3 hsv;\n"
"//             rgb2hsv(col, hsv);\n"
"//             float na;\n"
"//             lfnoise(x.xy-iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, col);\n"
"        }\n"
"    }\n"
"    col += col;\n"
"    col *= col;\n"
"    col *= mix(c.xxx, col, iScale);\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *graffiti_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void hsv2rgb(in vec3 hsv, out vec3 rgb);\n"
"void rgb2hsv(in vec3 rgb, out vec3 hsv);\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t);\n"
"void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"float ind;\n"
"void graf(in vec2 x, in vec2 size, in float w, out float dst)\n"
"{\n"
"    x.x += .3*iTime;\n"
"    x.y -= .1;\n"
"    \n"
"    vec2 y = vec2(mod(x.x, size.x)-.5*size.x, x.y);\n"
"    ind = (y.x-x.x)/size.x;\n"
"    vec3 r,\n"
"        d;\n"
"    \n"
"    size.x -= 2.*w; \n"
"    \n"
"    rand(vec2(ind,0.), r.x);\n"
"    rand(vec2(ind,1337.), r.y);\n"
"    rand(vec2(ind,2337.), r.z);\n"
"    \n"
"    vec2 pc = vec2(.5*size.x*(-1.+2.*r.y), 0.);\n"
"    \n"
"    // x component of r selects first actor\n"
"    if(r.x < .5) // Lin\n"
"    {\n"
"        vec2 p1 = vec2(.5*size.x*(-1.+4.*r.x),-.5*size.y),\n"
"            p2 = pc;\n"
"        dlinesegment(y, p1, p2, dst);\n"
"    }\n"
"    else // Quad\n"
"    {\n"
"        vec2 p1 = vec2(.5*size.x*(-1.+4.*(r.x-.5)),-.5*size.y),\n"
"            p2 = vec2(.5*size.x*(-1.+2.*r.y), -.5*size.y),\n"
"            p3 = pc;\n"
"        dspline2(y, p1, p2, p3, dst);\n"
"    }\n"
"    \n"
"    // y component of r selects connection type\n"
"    \n"
"    // z component of r selects second actor\n"
"    if(r.z < .5) // Lin\n"
"    {\n"
"        vec2 p1 = vec2(.5*size.x*(-1.+4.*r.z),.5*size.y),\n"
"            p2 = pc;\n"
"        dlinesegment(y, p1, p2, d.x);\n"
"    }\n"
"    else // Quad\n"
"    {\n"
"        vec2 p1 = vec2(.5*size.x*(-1.+4.*(r.z-.5)),.5*size.y),\n"
"            p2 = vec2(.5*size.x*(-1.+2.*r.y), .5*size.y),\n"
"            p3 = pc;\n"
"        dspline2(y, p1, p2, p3, d.x);\n"
"    }\n"
"    \n"
"    dst = min(dst, d.x);\n"
"    \n"
"    // Generate displacement\n"
"    lfnoise(12.*x, d.y);\n"
"    lfnoise(22.*x, d.z);\n"
"    d.y += .3*d.z;\n"
"    d.y = floor(.2+d.y);\n"
"    \n"
"    stroke(dst, w+.08*d.y, dst);\n"
"}\n"
"\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float n;\n"
"    lfnoise(5.*x.xy+iTime*c.xy, n);\n"
"    x.z += .03*n;\n"
"    \n"
"    float d = 1., da;\n"
"    for(float i=0.; i<2.; i+=1.)\n"
"    {\n"
"    	graf(x.xy-1337.333*i*c.xy, vec2(.4,.6), .05, da);\n"
"        //d = min(d, da);\n"
"        float rr;\n"
"        rand(i*c.xx*1.e2,rr); \n"
"        zextrude(x.z, -d+.5*abs(x.z), .07+.06*rr, d);\n"
"        stroke(d,mix(.02,.04,iScale), d);\n"
"        smoothmin(d, da, .2, d);\n"
"    }\n"
"    d = mix(1., d, smoothstep(19., 23., iTime));\n"
"    sdf = vec2(d, 2.);\n"
"    \n"
"    add(sdf, vec2(x.z+.05,1.), sdf);\n"
"}   \n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    x.x += .3*iTime;\n"
"    x.y -= .05;\n"
"    \n"
"    col = .5*c.xxx;\n"
"    \n"
"    float s = .1;\n"
"    vec2 dd = mod(x, s)-.5*s;\n"
"    stroke(dd.x, .005, dd.x);\n"
"    stroke(dd.y, .005, dd.y);\n"
"    col = mix(col, c.xxx, sm(min(dd.x, dd.y)));\n"
"    \n"
"    float d = 1., da;\n"
"    for(float i=0.; i<2.; i+=1.)\n"
"    {\n"
"    	graf(x.xy-1337.333*i*c.xy-.3*iTime*c.xy, vec2(.4,.6), .05, da);\n"
"        //d = min(d, da);\n"
"        float rr;\n"
"        rand(i*c.xx*1.e2,rr); \n"
"        stroke(da,.18, da);\n"
"        smoothmin(d, da, .2, d);\n"
"    }\n"
"    d = mix(1., d, smoothstep(15., 19., iTime));\n"
"    vec3 c1 = vec3(.78,.61*abs(2.*x.y),.15);\n"
"    \n"
"//     vec3 hsv;\n"
"//     rgb2hsv(c1, hsv);\n"
"//     float na;\n"
"//     lfnoise(x.xy-iTime+4.*hsv.x, na);\n"
"//     hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);\n"
"//     hsv2rgb(hsv, c1);\n"
"    \n"
"    col = mix(col, c1, sm(d));\n"
"    \n"
"    if(d != 1.)\n"
"    {\n"
"        stroke(d-.03, .03, d);\n"
"        col = mix(col, c.yyy, sm(d));\n"
"    }\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 150,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.08)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.1)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,5.e-3);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-3);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+.5*c.yzx);\n"
"            colorize(x.xy, col);\n"
"            col = .1*col\n"
"                + 1.*col * abs(dot(l,n))\n"
"                + 1.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"            \n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            lfnoise(x.xy, r);\n"
"            col = vec3(0.99,0.43,0.15);\n"
"            \n"
"            vec3 hsv;\n"
"//             rgb2hsv(col, hsv);\n"
"//             float na;\n"
"//             lfnoise(x.xy-iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na+iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, col);\n"
"            \n"
"            vec3 c2 = vec3(0.44,0.07,0.66);\n"
"//             rgb2hsv(c2, hsv);\n"
"//             lfnoise(x.xy+iTime+4.*hsv.x, na);\n"
"//             hsv.x = mod(1.*hsv.x+.2*na-iTime, 2.*pi);\n"
"//             hsv2rgb(hsv, c2);\n"
"            \n"
"            col = mix(col,c2,sin(2.*iScale*r*x));\n"
"            col = .1*col\n"
"                + .8*col * abs(dot(l,n))\n"
"                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"    }\n"
"    \n"
"    //col += col;\n"
"    \n"
"    col *= col;\n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"    col = mix(c.yyy, col, smoothstep(0., 1., iTime));\n"
"    col = mix(col, c.yyy, smoothstep(48.655, 49.655, iTime));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *bloodcells_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void hsv2rgb(in vec3 hsv, out vec3 rgb);\n"
"void rgb2hsv(in vec3 rgb, out vec3 hsv);\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dsmoothvoronoi(in vec2 x, out float d, out vec2 z)\n"
"{\n"
"    float n;\n"
"    lfnoise(x-iTime*c.xy, n);\n"
"    \n"
"    vec2 y = floor(x);\n"
"       float ret = 1.;\n"
"    vec2 pf=c.yy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            if(d < df)\n"
"            {\n"
"                df = d;\n"
"                pf = p;\n"
"            }\n"
"        }\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            vec2 o = p - pf;\n"
"            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"            smoothmin(ret, d, .4+.38*n, ret);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{    \n"
"    x.y += .3*iTime;\n"
"    float d;\n"
"    dsmoothvoronoi(mix(2.,3.,smoothstep(10.,12.,iTime))*x.xy-1337.,d,ind);\n"
"    stroke(d, .1, d);\n"
"    float modsize = .04,\n"
"		y = mod(d-.02*iTime,modsize)-.5*modsize,\n"
"        yi = (d-y)/modsize;\n"
"    \n"
"    float n;\n"
"    lfnoise(2.*yi*c.xx-.3*iTime, n);\n"
"    \n"
"    zextrude(x.z-.05*n, -y, mix(0.,.05+.05*n,iScale), d);\n"
"    \n"
"    stroke(d,mix(0.,.02,iScale),d);\n"
"    \n"
"    sdf = vec2(d, 2.);\n"
"    \n"
"    add(sdf, vec2(x.z+.05,1.), sdf);\n"
"}   \n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    col = .5*c.xxx;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));\n"
"    \n"
"    iScale *= (1.-smoothstep(51., 52., iTime));\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 200,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.1)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.1)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += s.x<1.e-2?min(s.x,5.e-4):s.x;\n"
"        //d += s.x<5.e-2?min(s.x,2.e-3):s.x;\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-3);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+.5*c.yzx);\n"
"            colorize(x.xy, col);\n"
"            col = .1*col\n"
"                + 1.*col * abs(dot(l,n))\n"
"                + 1.5 * col * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            lfnoise(x.xy-iTime, r);\n"
"            col = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*r*x));\n"
"            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*r*x));\n"
"            col = mix(col, c1, .5+.5*r);\n"
"            col = .1*col\n"
"                + .8*col * abs(dot(l,n))\n"
"                + 6.5*col * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"    }\n"
"    \n"
"    col *= col*col;\n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"    col *= mix(col, length(col)/sqrt(3.)*c.xxx, iScale);\n"
"\n"
"    col = mix(c.yyy, col, smoothstep(0., 1., iTime));\n"
"    col = mix(col, c.yyy, smoothstep(53.69, 54.69, iTime));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
void Loadhsv2rgb()
{
    int hsv2rgb_size = strlen(hsv2rgb_source);
    hsv2rgb_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(hsv2rgb_handle, 1, (GLchar **)&hsv2rgb_source, &hsv2rgb_size);
    glCompileShader(hsv2rgb_handle);
#ifdef DEBUG
    printf("---> hsv2rgb Shader:\n");
    debug(hsv2rgb_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrgb2hsv()
{
    int rgb2hsv_size = strlen(rgb2hsv_source);
    rgb2hsv_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rgb2hsv_handle, 1, (GLchar **)&rgb2hsv_source, &rgb2hsv_size);
    glCompileShader(rgb2hsv_handle);
#ifdef DEBUG
    printf("---> rgb2hsv Shader:\n");
    debug(rgb2hsv_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrand()
{
    int rand_size = strlen(rand_source);
    rand_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rand_handle, 1, (GLchar **)&rand_source, &rand_size);
    glCompileShader(rand_handle);
#ifdef DEBUG
    printf("---> rand Shader:\n");
    debug(rand_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadlfnoise()
{
    int lfnoise_size = strlen(lfnoise_source);
    lfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(lfnoise_handle, 1, (GLchar **)&lfnoise_source, &lfnoise_size);
    glCompileShader(lfnoise_handle);
#ifdef DEBUG
    printf("---> lfnoise Shader:\n");
    debug(lfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadmfnoise()
{
    int mfnoise_size = strlen(mfnoise_source);
    mfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(mfnoise_handle, 1, (GLchar **)&mfnoise_source, &mfnoise_size);
    glCompileShader(mfnoise_handle);
#ifdef DEBUG
    printf("---> mfnoise Shader:\n");
    debug(mfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadstroke()
{
    int stroke_size = strlen(stroke_source);
    stroke_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(stroke_handle, 1, (GLchar **)&stroke_source, &stroke_size);
    glCompileShader(stroke_handle);
#ifdef DEBUG
    printf("---> stroke Shader:\n");
    debug(stroke_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadadd()
{
    int add_size = strlen(add_source);
    add_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(add_handle, 1, (GLchar **)&add_source, &add_size);
    glCompileShader(add_handle);
#ifdef DEBUG
    printf("---> add Shader:\n");
    debug(add_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddvoronoi()
{
    int dvoronoi_size = strlen(dvoronoi_source);
    dvoronoi_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dvoronoi_handle, 1, (GLchar **)&dvoronoi_source, &dvoronoi_size);
    glCompileShader(dvoronoi_handle);
#ifdef DEBUG
    printf("---> dvoronoi Shader:\n");
    debug(dvoronoi_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadnormal()
{
    int normal_size = strlen(normal_source);
    normal_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(normal_handle, 1, (GLchar **)&normal_source, &normal_size);
    glCompileShader(normal_handle);
#ifdef DEBUG
    printf("---> normal Shader:\n");
    debug(normal_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrot3()
{
    int rot3_size = strlen(rot3_source);
    rot3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rot3_handle, 1, (GLchar **)&rot3_source, &rot3_size);
    glCompileShader(rot3_handle);
#ifdef DEBUG
    printf("---> rot3 Shader:\n");
    debug(rot3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox()
{
    int dbox_size = strlen(dbox_source);
    dbox_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox_handle, 1, (GLchar **)&dbox_source, &dbox_size);
    glCompileShader(dbox_handle);
#ifdef DEBUG
    printf("---> dbox Shader:\n");
    debug(dbox_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlinesegment3()
{
    int dlinesegment3_size = strlen(dlinesegment3_source);
    dlinesegment3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlinesegment3_handle, 1, (GLchar **)&dlinesegment3_source, &dlinesegment3_size);
    glCompileShader(dlinesegment3_handle);
#ifdef DEBUG
    printf("---> dlinesegment3 Shader:\n");
    debug(dlinesegment3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadzextrude()
{
    int zextrude_size = strlen(zextrude_source);
    zextrude_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(zextrude_handle, 1, (GLchar **)&zextrude_source, &zextrude_size);
    glCompileShader(zextrude_handle);
#ifdef DEBUG
    printf("---> zextrude Shader:\n");
    debug(zextrude_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlinesegment()
{
    int dlinesegment_size = strlen(dlinesegment_source);
    dlinesegment_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlinesegment_handle, 1, (GLchar **)&dlinesegment_source, &dlinesegment_size);
    glCompileShader(dlinesegment_handle);
#ifdef DEBUG
    printf("---> dlinesegment Shader:\n");
    debug(dlinesegment_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddspline2()
{
    int dspline2_size = strlen(dspline2_source);
    dspline2_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspline2_handle, 1, (GLchar **)&dspline2_source, &dspline2_size);
    glCompileShader(dspline2_handle);
#ifdef DEBUG
    printf("---> dspline2 Shader:\n");
    debug(dspline2_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadsmoothmin()
{
    int smoothmin_size = strlen(smoothmin_source);
    smoothmin_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(smoothmin_handle, 1, (GLchar **)&smoothmin_source, &smoothmin_size);
    glCompileShader(smoothmin_handle);
#ifdef DEBUG
    printf("---> smoothmin Shader:\n");
    debug(smoothmin_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}

void LoadSymbols()
{
    Loadhsv2rgb();
    updateBar();
    Loadrgb2hsv();
    updateBar();
    Loadrand();
    updateBar();
    Loadlfnoise();
    updateBar();
    Loadmfnoise();
    updateBar();
    Loadstroke();
    updateBar();
    Loadadd();
    updateBar();
    Loaddvoronoi();
    updateBar();
    Loadnormal();
    updateBar();
    Loadrot3();
    updateBar();
    Loaddbox();
    updateBar();
    Loaddlinesegment3();
    updateBar();
    Loadzextrude();
    updateBar();
    Loaddlinesegment();
    updateBar();
    Loaddspline2();
    updateBar();
    Loadsmoothmin();
    updateBar();
}
int voronoidesign_program, voronoidesign_handle, groundboxes_program, groundboxes_handle, graffiti_program, graffiti_handle, bloodcells_program, bloodcells_handle;
int voronoidesign_iTime_location;
voronoidesign_iResolution_location;
int groundboxes_iTime_location;
groundboxes_iResolution_location;
int graffiti_iTime_location;
graffiti_iResolution_location;
int bloodcells_iTime_location;
bloodcells_iResolution_location;
const int nprograms = 4;

void Loadvoronoidesign()
{
    int voronoidesign_size = strlen(voronoidesign_source);
    voronoidesign_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(voronoidesign_handle, 1, (GLchar **)&voronoidesign_source, &voronoidesign_size);
    glCompileShader(voronoidesign_handle);
#ifdef DEBUG
    printf("---> voronoidesign Shader:\n");
    debug(voronoidesign_handle);
    printf(">>>>\n");
#endif
    voronoidesign_program = glCreateProgram();
    glAttachShader(voronoidesign_program,voronoidesign_handle);
    glAttachShader(voronoidesign_program,hsv2rgb_handle);
    glAttachShader(voronoidesign_program,rgb2hsv_handle);
    glAttachShader(voronoidesign_program,rand_handle);
    glAttachShader(voronoidesign_program,lfnoise_handle);
    glAttachShader(voronoidesign_program,mfnoise_handle);
    glAttachShader(voronoidesign_program,stroke_handle);
    glAttachShader(voronoidesign_program,add_handle);
    glAttachShader(voronoidesign_program,dvoronoi_handle);
    glAttachShader(voronoidesign_program,normal_handle);
    glLinkProgram(voronoidesign_program);
#ifdef DEBUG
    printf("---> voronoidesign Program:\n");
    debugp(voronoidesign_program);
    printf(">>>>\n");
#endif
    glUseProgram(voronoidesign_program);
    voronoidesign_iTime_location = glGetUniformLocation(voronoidesign_program, "iTime");
    voronoidesign_iResolution_location = glGetUniformLocation(voronoidesign_program, "iResolution");
    progress += .2/(float)nprograms;
}

void Loadgroundboxes()
{
    int groundboxes_size = strlen(groundboxes_source);
    groundboxes_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(groundboxes_handle, 1, (GLchar **)&groundboxes_source, &groundboxes_size);
    glCompileShader(groundboxes_handle);
#ifdef DEBUG
    printf("---> groundboxes Shader:\n");
    debug(groundboxes_handle);
    printf(">>>>\n");
#endif
    groundboxes_program = glCreateProgram();
    glAttachShader(groundboxes_program,groundboxes_handle);
    glAttachShader(groundboxes_program,hsv2rgb_handle);
    glAttachShader(groundboxes_program,rgb2hsv_handle);
    glAttachShader(groundboxes_program,rand_handle);
    glAttachShader(groundboxes_program,lfnoise_handle);
    glAttachShader(groundboxes_program,rot3_handle);
    glAttachShader(groundboxes_program,dbox_handle);
    glAttachShader(groundboxes_program,dlinesegment3_handle);
    glAttachShader(groundboxes_program,stroke_handle);
    glAttachShader(groundboxes_program,zextrude_handle);
    glAttachShader(groundboxes_program,add_handle);
    glAttachShader(groundboxes_program,normal_handle);
    glLinkProgram(groundboxes_program);
#ifdef DEBUG
    printf("---> groundboxes Program:\n");
    debugp(groundboxes_program);
    printf(">>>>\n");
#endif
    glUseProgram(groundboxes_program);
    groundboxes_iTime_location = glGetUniformLocation(groundboxes_program, "iTime");
    groundboxes_iResolution_location = glGetUniformLocation(groundboxes_program, "iResolution");
    progress += .2/(float)nprograms;
}

void Loadgraffiti()
{
    int graffiti_size = strlen(graffiti_source);
    graffiti_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(graffiti_handle, 1, (GLchar **)&graffiti_source, &graffiti_size);
    glCompileShader(graffiti_handle);
#ifdef DEBUG
    printf("---> graffiti Shader:\n");
    debug(graffiti_handle);
    printf(">>>>\n");
#endif
    graffiti_program = glCreateProgram();
    glAttachShader(graffiti_program,graffiti_handle);
    glAttachShader(graffiti_program,hsv2rgb_handle);
    glAttachShader(graffiti_program,rgb2hsv_handle);
    glAttachShader(graffiti_program,rand_handle);
    glAttachShader(graffiti_program,lfnoise_handle);
    glAttachShader(graffiti_program,dlinesegment_handle);
    glAttachShader(graffiti_program,dspline2_handle);
    glAttachShader(graffiti_program,stroke_handle);
    glAttachShader(graffiti_program,zextrude_handle);
    glAttachShader(graffiti_program,smoothmin_handle);
    glAttachShader(graffiti_program,add_handle);
    glAttachShader(graffiti_program,normal_handle);
    glLinkProgram(graffiti_program);
#ifdef DEBUG
    printf("---> graffiti Program:\n");
    debugp(graffiti_program);
    printf(">>>>\n");
#endif
    glUseProgram(graffiti_program);
    graffiti_iTime_location = glGetUniformLocation(graffiti_program, "iTime");
    graffiti_iResolution_location = glGetUniformLocation(graffiti_program, "iResolution");
    progress += .2/(float)nprograms;
}

void Loadbloodcells()
{
    int bloodcells_size = strlen(bloodcells_source);
    bloodcells_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(bloodcells_handle, 1, (GLchar **)&bloodcells_source, &bloodcells_size);
    glCompileShader(bloodcells_handle);
#ifdef DEBUG
    printf("---> bloodcells Shader:\n");
    debug(bloodcells_handle);
    printf(">>>>\n");
#endif
    bloodcells_program = glCreateProgram();
    glAttachShader(bloodcells_program,bloodcells_handle);
    glAttachShader(bloodcells_program,hsv2rgb_handle);
    glAttachShader(bloodcells_program,rgb2hsv_handle);
    glAttachShader(bloodcells_program,rand_handle);
    glAttachShader(bloodcells_program,lfnoise_handle);
    glAttachShader(bloodcells_program,stroke_handle);
    glAttachShader(bloodcells_program,zextrude_handle);
    glAttachShader(bloodcells_program,smoothmin_handle);
    glAttachShader(bloodcells_program,add_handle);
    glAttachShader(bloodcells_program,normal_handle);
    glLinkProgram(bloodcells_program);
#ifdef DEBUG
    printf("---> bloodcells Program:\n");
    debugp(bloodcells_program);
    printf(">>>>\n");
#endif
    glUseProgram(bloodcells_program);
    bloodcells_iTime_location = glGetUniformLocation(bloodcells_program, "iTime");
    bloodcells_iResolution_location = glGetUniformLocation(bloodcells_program, "iResolution");
    progress += .2/(float)nprograms;
}

void LoadPrograms()
{
    Loadvoronoidesign();
    updateBar();
    Loadgroundboxes();
    updateBar();
    Loadgraffiti();
    updateBar();
    Loadbloodcells();
    updateBar();
}
#endif
