//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

extern float progress;int rand_handle, lfnoise_handle, dspline3_handle, dbox3_handle, dlinesegment3_handle, stroke_handle, zextrude_handle, smoothmin_handle, add_handle, dvoronoi_handle, dbox_handle, rot3_handle, normal_handle, mfnoise_handle, dtriangle_handle, dlinesegment_handle, dpolygon_handle, rot_handle, dcircle_handle, hash31_handle, dschnappsgirls_handle, dspacepigs_handle, dkewlers_handle, dfarbrausch_handle, dhaujobb_handle, dmercury_handle;
const int nsymbols = 26;
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
const char *dspline3_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"\n"
"//distance to spline with parameter t\n"
"float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum dist3ance to spline\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,\n"
"    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist3(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three dist3inct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist3(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist3(p0,p1,p2,x,t.y),\n"
"            dist3(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\0";
const char *dbox3_source = "#version 130\n\n"
"void dbox3(in vec3 x, in vec3 b, out float d)\n"
"{\n"
"  vec3 da = abs(x) - b;\n"
"  d = length(max(da,0.0))\n"
"         + min(max(da.x,max(da.y,da.z)),0.0);\n"
"}\n"
"\0";
const char *dlinesegment3_source = "#version 130\n\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *stroke_source = "// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\0";
const char *zextrude_source = "// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\0";
const char *smoothmin_source = "// iq's smooth minimum\n"
"void smoothmin(in float a, in float b, in float k, out float dst)\n"
"{\n"
"    float h = max( k-abs(a-b), 0.0 )/k;\n"
"    dst = min( a, b ) - h*h*h*k*(1.0/6.0);\n"
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
const char *dbox_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
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
const char *dtriangle_source = "// Adapted from iq, https://www.shadertoy.com/view/XsXSz4\n"
"void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst)\n"
"{\n"
"	vec2 e0 = p1 - p0;\n"
"	vec2 e1 = p2 - p1;\n"
"	vec2 e2 = p0 - p2;\n"
"\n"
"	vec2 v0 = p - p0;\n"
"	vec2 v1 = p - p1;\n"
"	vec2 v2 = p - p2;\n"
"\n"
"	vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );\n"
"	vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );\n"
"	vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );\n"
"    \n"
"    float s = sign( e0.x*e2.y - e0.y*e2.x );\n"
"    vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),\n"
"                       vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),\n"
"                       vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));\n"
"\n"
"	dst = -sqrt(d.x)*sign(d.y);\n"
"}\n"
"\0";
const char *dlinesegment_source = "#version 130\n\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *dpolygon_source = "#version 130\n\n"
"// compute distance to regular polygon\n"
"const float pi = acos(-1.);\n"
"void dpolygon(in vec2 x, in float N, out float d)\n"
"{\n"
"    d = 2.0*pi/N;\n"
"    float t = mod(acos(x.x/length(x)), d)-0.5*d;\n"
"    d = -0.5+length(x)*cos(t)/cos(0.5*d);\n"
"}\n"
"\0";
const char *rot_source = "#version 130\n\n"
"void rot(in float phi, out mat2 m)\n"
"{\n"
"    vec2 cs = vec2(cos(phi), sin(phi));\n"
"    m = mat2(cs.x, -cs.y, cs.y, cs.x);\n"
"}\n"
"\0";
const char *dcircle_source = "#version 130\n\n"
"\n"
"void dcircle(in vec2 x, out float d)\n"
"{\n"
"    d = length(x)-1.0;\n"
"}\n"
"\0";
const char *hash31_source = "// Creative Commons Attribution-ShareAlike 4.0 International Public License\n"
"// Created by David Hoskins.\n"
"// See https://www.shadertoy.com/view/4djSRW\n"
"void hash31(in float p, out vec3 d)\n"
"{\n"
"   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));\n"
"   p3 += dot(p3, p3.yzx+33.33);\n"
"   d = fract((p3.xxy+p3.yzz)*p3.zyx); \n"
"}\n"
"\0";
const char *dschnappsgirls_source = "#version 130\n\n"
"\n"
"void dtriangle(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dschnappsgirls(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    // Dress\n"
"    dtriangle(x, vec2(-.1,-.3), vec2(.5,-.3), vec2(.2, .6), d0);\n"
"    dlinesegment(x, vec2(-.1,.325), vec2(.5,.325), da);\n"
"    stroke(da,.06,da);\n"
"    d0 = max(d0,-da);\n"
"    \n"
"    // Head\n"
"    dcircle(7.*(x-vec2(.2,.5)), da);\n"
"    d0 = max(d0, -da+.5);\n"
"    d0 = min(d0, da/7.);\n"
"    \n"
"    // Legs\n"
"    dlinesegment(x, vec2(.125,-.3), vec2(.125,-.6), da);\n"
"    stroke(da, .06, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(.275,-.3), vec2(.275,-.6), da);\n"
"    stroke(da, .06, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Shoulders\n"
"    dlinesegment(x, vec2(0.05,.25), vec2(.35,.25), da);\n"
"    stroke(da, .085, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Arms\n"
"    dlinesegment(x, vec2(.385,.25), vec2(.5, -.1), da);\n"
"    stroke(da, .055, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(.017,.25), vec2(-.1, -.1), da);\n"
"    stroke(da, .055, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Glass\n"
"    dtriangle(x, vec2(-.6,.3), vec2(-.4,.1), vec2(-.2,.3), da);\n"
"    stroke(da, .0125, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(-.4,.15), vec2(-.4,-.1), da);\n"
"    stroke(da, .0125, da);\n"
"    d0 = min(d0, da);\n"
"    dtriangle(x, vec2(-.5,-.15), vec2(-.3,-.15), vec2(-.4,-.1), da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Liquid\n"
"    dtriangle(x, vec2(-.55,.25), vec2(-.4,.1), vec2(-.25,.25), da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Salad\n"
"    dlinesegment(x, vec2(-.4,.1), vec2(-.2,.5), da);\n"
"    stroke(da, .01, da);\n"
"    d0 = min(d0, da);\n"
"    dcircle(24.*(x-vec2(-.3,.3)), da);\n"
"    d0 = min(d0, da/24.);\n"
"    dcircle(24.*(x-vec2(-.25,.4)), da);\n"
"    d0 = min(d0, da/24.);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dspacepigs_source = "#version 130\n\n"
"\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void dcircle(in vec2 x, out float d);\n"
"\n"
"void dear(in vec2 x, out float d)\n"
"{\n"
"    d = abs(2.*x.y)\n"
"        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))\n"
"        -.5*min(-abs(x.x),.01);\n"
"}\n"
"\n"
"void dspacepigs(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    // Head\n"
"    dcircle(2.5*x,d0);\n"
"    d0 /= 2.5;\n"
"    \n"
"    // Ears\n"
"    dear(vec2(2.,5.)*x-vec2(.8,1.3), da);\n"
"    d0 = min(d0,da/10.);\n"
"    dear(vec2(2.,5.)*x+vec2(.8,-1.3), da);\n"
"    d0 = min(d0,da/10.);\n"
"    \n"
"    // Nose\n"
"    dcircle(6.*x-vec2(0.,-.5),da);\n"
"    d0 = max(d0,-da/6.);\n"
"    dcircle(24.*x-vec2(-1.5,-2.),da);\n"
"    d0 = min(d0,da/24.);\n"
"    dcircle(24.*x-vec2(1.5,-2.),da);\n"
"    d0 = min(d0,da/24.);\n"
"    \n"
"    // Eyes\n"
"    dcircle(16.*x-vec2(-3.5,2.5),da);\n"
"    d0 = max(d0,-da/16.);\n"
"    dcircle(16.*x-vec2(3.5,2.5),da);\n"
"    d0 = max(d0,-da/16.);\n"
"    dcircle(24.*x-vec2(-5.,3.5),da);\n"
"    d0 = min(d0,da/24.);\n"
"    dcircle(24.*x-vec2(5.,3.5),da);\n"
"    d0 = min(d0,da/24.);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dkewlers_source = "#version 130\n\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dkewlers(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    x *= 1.2;\n"
"    \n"
"    dbox(x-vec2(0.,-.3),vec2(.6,.1),d0);\n"
"    dbox(x-vec2(-.5,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(-.5+1./3.,.25),vec2(.1,.5),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(-.5+2./3.,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(.5,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dfarbrausch_source = "#version 130\n\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"\n"
"void dfarbrausch(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    x += vec2(.1,0.);\n"
"    x *= 1.2;\n"
"    \n"
"    dlinesegment(x,vec2(-.65,.05),vec2(-.5,.05),d0);\n"
"    dlinesegment(x,vec2(-.5,.05),vec2(-.2,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.2,-.49),vec2(-.0,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.0,-.49),vec2(-.27,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.07,0.),vec2(-.27,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.2,-.49),vec2(-.07,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.4,-.49),vec2(.13,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.4,-.49),vec2(.2,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.33,0.),vec2(.13,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.33,0.),vec2(.51,-.33),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.6,-.15),vec2(.51,-.33),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.53,0.),vec2(.6,-.15),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.7,0.),vec2(.53,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.7,0.),vec2(.68,-.04),da);\n"
"    d0 = min(d0, da);\n"
"    dpolygon(5.*(x+vec2(.3,.65)),6.,da);\n"
"    d0 = min(d0, da/5.);\n"
"    dpolygon(5.*(x+vec2(-.5,.65)),6.,da);\n"
"    d0 = min(d0, da/5.);\n"
"    \n"
"    stroke(d0,.035, d0);\n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dhaujobb_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void rot(in float phi, out mat2 m);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"\n"
"void dhaujobb(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    mat2 m;\n"
"	rot(.3,m);\n"
"    x = 1.1*x*m;\n"
"    x.x *= 1.1;\n"
"        \n"
"    x += vec2(-.05,.2);\n"
"    \n"
"    // Left leg\n"
"    dbox(x+.35*c.xx,vec2(.1,.05),d0);\n"
"    dbox(x+vec2(.3,.25),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x+vec2(.2,.15),vec2(.1,.05),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x+vec2(.15,.05),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Right leg\n"
"    dbox(x-vec2(.65,.35),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"\n"
"    // Torso\n"
"    rot(.2, m);\n"
"    dbox(m*(x-vec2(.25,.15)),vec2(.45,.05),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(m*(x-vec2(-.15,.35)),vec2(.45,.05),da);\n"
"    d0 = min(d0,da);\n"
"    rot(pi/8.,m);\n"
"    dbox(m*(x-vec2(.0,.25)),vec2(.1,.15),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Penis\n"
"    dbox(m*(x-vec2(.1,-.0)),vec2(.025,.1),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Left hand\n"
"    rot(.3,m);\n"
"    dbox(m*(x-vec2(.235,.535)),vec2(.035,.15),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(m*(x-vec2(.225,.7)),vec2(.075,.025),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Right hand\n"
"    rot(-.2,m);\n"
"    dbox(m*(x+vec2(.585,-.2)),vec2(.0375,.1),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Head\n"
"    dcircle(6.*(x-vec2(-.15,.58)),da);\n"
"    d0 = min(d0,da/6.);\n"
"    \n"
"    d0 -= .05*(abs(x.x)+abs(x.y)-.2);\n"
"    d = max(d,-d0);\n"
"}\n"
"\0";
const char *dmercury_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dmercury(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da;\n"
"\n"
"    x += .1*c.yx;\n"
"\n"
"    // Upper part\n"
"    dbox(x-.35*c.yx,vec2(.4,.35), da);\n"
"    d = max(d, -da);\n"
"    dbox(x-.7*c.yx, vec2(.2,.2), da);\n"
"    d = min(d,da);\n"
"    dbox(x-.25*c.yx,vec2(.2,.05),da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Lower part\n"
"    dbox(x+.2*c.yx,vec2(.1,.4),da);\n"
"    d = max(d, -da);\n"
"    dbox(x+.2*c.yx, vec2(.4,.1),da);\n"
"    d = max(d, -da);\n"
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
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"// Creative Commons Attribution-ShareAlike 4.0 International Public License\n"
"// Created by David Hoskins.\n"
"// See https://www.shadertoy.com/view/4djSRW\n"
"void hash31(in float p, out vec3 d)\n"
"{\n"
"   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));\n"
"   p3 += dot(p3, p3.yzx+33.33);\n"
"   d = fract((p3.xxy+p3.yzz)*p3.zyx); \n"
"}\n"
"// End of Hoskins Hash\n"
"\n"
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
"\n"
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
"\n"
"// Box sdf\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\n"
"// Add sdfs\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = sda.x<sdb.x?sda:sdb;\n"
"}\n"
"\n"
"// iq's smooth minimum\n"
"void smoothmin(in float a, in float b, in float k, out float dst)\n"
"{\n"
"    float h = max( k-abs(a-b), 0.0 )/k;\n"
"    dst = min( a, b ) - h*h*h*k*(1.0/6.0);\n"
"}\n"
"\n"
"//distance to spline with parameter t\n"
"float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum dist3ance to spline\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,\n"
"    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist3(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three dist3inct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist3(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist3(p0,p1,p2,x,t.y),\n"
"            dist3(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\n"
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
"\n"
"vec2 vind,vind2;\n"
"float v, fn, r1, fb;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.y += .2*iTime;\n"
"    \n"
"    dvoronoi(1.5*x.xy, v, vind);\n"
"    \n"
"    vec3 y = vec3(vind/1.5-x.xy,x.z);\n"
"    \n"
"    float n, n2;\n"
"    \n"
"    lfnoise(c.xx-.3*iTime+vind*3., n);\n"
"    lfnoise(5.*x.z*c.xx-iTime-vind*4., n2);\n"
"    n2 *= .2;\n"
"    \n"
"    mat2 RR = mat2(cos(n2), sin(n2), -sin(n2), cos(n2));\n"
"    vec2 a = x.xy;\n"
"    x.xy = RR * x.xy;\n"
"    rand(vind, r1);\n"
"    float r2;\n"
"    rand(vind -1336., r2);\n"
"    \n"
"    float phi = atan(y.y, y.x),\n"
"        dp = pi/24.,\n"
"        phii = mod(phi, dp)-.5*dp,\n"
"        pa = phi - phii, \n"
"        R1 = .05,\n"
"        R2 = mix(.4,.25,1.-r2);\n"
"    \n"
"    R2 = mix(R1, R2, .5+.5*n2);\n"
"    \n"
"    float r0;\n"
"    rand(pa*c.xx, r0);\n"
"    r0 = mix(r0,.5+.5*n,.5);\n"
"    \n"
"    dspline3(y, vec3(1.4*R1*cos(pa), 1.4*R1*sin(pa), -.5), vec3(R1*cos(pa), R1*sin(pa), .1*r1), vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), sdf.x);\n"
"    float da;\n"
"    dspline3(y, vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0), da);\n"
"    sdf.x = min(sdf.x, da);\n"
"    stroke(sdf.x, .25*mix(.02,.05, .5+.5*n2), sdf.x);\n"
"    sdf.y = 2.;\n"
"    \n"
"    add(sdf, vec2(length(y-vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0))-.01, 3.), sdf);\n"
"    \n"
"    float fa;\n"
"    lfnoise(4.*a,  fa);\n"
"    dvoronoi(a,fn, vind2); \n"
"    fa = x.z+.4+.1*mix((v+fn),fa,.5);\n"
"    add(sdf, vec2(fa,4.), sdf);\n"
"    smoothmin(sdf.x, fa, .1, sdf.x);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 1.e-4;\n"
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
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void vs(in vec3 x, out vec2 sdf)\n"
"{\n"
"    vec2 vi;\n"
"    dvoronoi(3.*(x.xy+.2*iTime*c.yx), sdf.x, vi);\n"
"    sdf.x = x.z-.1-.2*sdf.x;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, \n"
"        s;\n"
"    \n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));\n"
"    \n"
"    uv *= 2.;\n"
"    \n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx), \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.3)/dir.z;\n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        vs(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        d += s.x;\n"
"    }\n"
"    float v1, rar, dx = 1.e-3;\n"
"    vec2 vi1, na;\n"
"    vec3 cv, l;\n"
"    x = o + d * dir;\n"
"    dvoronoi(3.*(x.xy+.2*iTime*c.yx), v1, vi1);\n"
"\n"
"    rand(vi1, rar);\n"
"    cv = mix(c.yyy,vec3(.23,.23,.23), rar);\n"
"    v1 = abs(v1)-.01;\n"
"    cv = mix(cv, c.yyy, sm(v1));\n"
"    v1 = abs(v1-.01)-.005;\n"
"    cv = mix(cv, /*vec3(0.76,0.20,0.13)*/c.xxx, sm(v1));\n"
"    vs(x,s);\n"
"    vs(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    vs(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    vs(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"    l = normalize(x+c.yyx);\n"
"    cv = .2*cv\n"
"        +.2*cv*abs(dot(l,n))\n"
"        +.6*cv*pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"    dir = refract(dir, n, .99);\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"\n"
"        d += min(s.x,3.e-2);\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n);\n"
"       \n"
"        if(s.y == 3.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            \n"
"            col = mix(c.xxx, vec3(0.76,0.20,0.13), .8);\n"
"            col = .4*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .8*col * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"        }\n"
"        else if(s.y == 4.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            rand(vind+vind2,r);\n"
"            \n"
"            col = mix(.023*c.xxx, vec3(0.76,0.40,0.23), r);\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            stroke(v, .01, v);\n"
"            stroke(fn, .01, fn);\n"
"            col = mix(col, c.yyy, sm(v));\n"
"            col = mix(col, col*col, sm(fn));\n"
"        }\n"
"		if(s.y == 2.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            rand(vind, r);\n"
"            \n"
"            col = mix(mix(vec3(0.76,0.20,0.23),vec3(.18,.32,.13), r), vec3(0.23,0.23,0.23),  clamp((x.z)/r1/.1,0.,1.));\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            col = mix(col, 5.*col, .25*n.x*n.x);\n"
"        }\n"
"        \n"
"    }\n"
"    \n"
"    col *=3.6;\n"
"    col *= col;\n"
"\n"
"    if(s.y != 3.)\n"
"    {\n"
"        col = mix(length(col)/sqrt(3.)*c.xxx, col, mix(.3,1.,iScale));\n"
"    }\n"
"    \n"
"    col = mix(col,cv,.8);\n"
"    col = mix(col, c.yyy, smoothstep(1.,5.,d));\n"
"\n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}	\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *groundboxes_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" *\n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"vec2 ind=c.yy;\n"
"\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{    \n"
"    \n"
"    float d,\n"
"        s = .1;\n"
"\n"
"    sdf = c.xy;\n"
"    \n"
"    for(float size = .0; size <= .5; size += .025)\n"
"    {\n"
"        dbox3(x, size*c.xxx, d);\n"
"        stroke(d, .001, d);\n"
"        vec2 sda = vec2(d,3.+size);\n"
"\n"
"        float n;\n"
"        vec3 y = mod(x,.125*size)-.5*.125*size,\n"
"            yi = (x-y)/size;\n"
"		ind = yi.xy+yi.yz+yi.xz;\n"
"        lfnoise(3.6*ind+15.*size-1.1*iTime, n); // TODO: use shifts to 13.6 from 3.6\n"
"        \n"
"        if(n>-.3)\n"
"        {\n"
"            \n"
"            dbox3(y, .25*size*c.xxx, d);\n"
"            sda.x = max(sda.x, -d);\n"
"        }\n"
"        \n"
"        add(sdf, sda, sdf);\n"
"    }   \n"
"}  \n"
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
"    mat3 RR;\n"
"    rot3(/*-pi/4.*c.yxy*/.2*iTime*vec3(1.1,1.4,1.6), RR);\n"
"    \n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = RR * c.yyx,\n"
"        r = normalize(c.xyy), \n"
"        u = normalize(c.yxy),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x,\n"
"        size = .301*c.xxx;\n"
"        \n"
"    int N = 150,\n"
"        i = 0;\n"
"    t = uv.x * r + uv.y * u;\n"
"    t = RR * t;\n"
"    dir = normalize(t-o);\n"
"    float d = 0., ra, ss, inside = 0., dd;\n"
"    \n"
"    for(ss = .5; ss >= .0; ss -= .025)\n"
"    {\n"
"        size = ss*c.xxx+1.e-4;\n"
"\n"
"        vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes\n"
"        vec2 abxlo = abs(o.yz + tlo.x*dir.yz),\n"
"            abylo = abs(o.xz + tlo.y*dir.xz),\n"
"            abzlo = abs(o.xy + tlo.z*dir.xy);\n"
"\n"
"        vec4 dn = 100.*c.xyyy;\n"
"        \n"
"        dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));\n"
"        dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));\n"
"        dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));\n"
"\n"
"        inside += .05;\n"
"        \n"
"        if(ss == 3.)dd = dn.r;\n"
"        d = dn.r;\n"
"        \n"
"        x = o + d * dir;\n"
"        scene(x, s);\n"
"        if(s.x < 1.e-4)break;\n"
"        \n"
"                normal(x,n, 5.e-4);\n"
"        \n"
"        lfnoise(x.xy*vec2(3.,8.)-iTime, ra);\n"
"        \n"
"        vec3 f;\n"
"        float dd = 5.e-1;\n"
"        \n"
"        r = RR*c.xyy;\n"
"        f = RR*c.yzy;\n"
"        u = RR*c.yyx;\n"
"        \n"
"        vec3 dp = abs(vec3(dot(n,r), dot(n,f), dot(n,u)));\n"
"        \n"
"        if(dp.y < dd && dp.z < dd) n = r;\n"
"        else if(dp.x < dd && dp.z < dd) n = f;\n"
"        else if(dp.x < dd && dp.y < dd) n = u;\n"
"        s.y = -1.;\n"
"\n"
"        vec3 l = normalize(x+.03*normalize(x-o));\n"
"        \n"
"        vec3 c1 = col;\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            colorize(x.xy, c1);\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"            \n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            \n"
"            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));\n"
"            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*ra*x));\n"
"            c1 = mix(c1, c1, .5+.5*ra);\n"
"            c1 = .3*c.xxx;\n"
"            c1 = .1*c1\n"
"                + .4*c1 * abs(dot(l,n))\n"
"                + .3*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"        else if(s.y >= 3.)\n"
"        {\n"
"            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);\n"
"            float na;\n"
"            ind = x.xy + x.yz + x.zx;\n"
"            ind = mod(ind, .01)-.005;\n"
"            ind = x.xy + x.yz + x.zx-ind;\n"
"            rand(ind, na);\n"
"            na = 0.;\n"
"            c1 = .8+.2*c1;\n"
"            c1 *= na;\n"
"            c1 = mix(.1,mix(.2,.4,iScale),step(na,.05))*c1\n"
"                + mix(.1,.2,step(na,.95))*c1 * abs(dot(l,n))\n"
"                + .5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"        }\n"
"        else if(s.y == -1.)\n"
"        {\n"
"            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);\n"
"            c1 = .8+.2*c1;\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"                \n"
"            vec3 dc;\n"
"\n"
"            vec3 zz = mod(x, .025)-.5*.025, zi = x-zz;\n"
"            rand(zi.xy+zi.yz+zi.zx, dc.x);\n"
"            rand(zi.xy+zi.yz+zi.zx+1337., dc.y);\n"
"            rand(zi.xy+zi.yz+zi.zx+2337., dc.z);\n"
"            float da;\n"
"            dbox3(zz, .01*c.xxx, da);\n"
"            stroke(da, .001, da);\n"
"            c1 = mix(c1, 1.2*c1*c1+.2*dc, sm(da));\n"
"            stroke(da-.002, .001, da);\n"
"            c1 = mix(c1, 1.6*c1*c1, sm(da));\n"
"        }\n"
"        else if(s.y == -2.)\n"
"        {\n"
"            float s = .05;\n"
"            c1 = .5*c.xxx;\n"
"            vec2 dd = mod(x.xz, s)-.5*s;\n"
"            stroke(dd.x, .005, dd.x);\n"
"            stroke(dd.y, .005, dd.y);\n"
"            c1 = mix(c1, c.xxx, sm(min(dd.x, dd.y)));\n"
"        }\n"
"        col = mix(col, c1, mix(.2,.6,iScale));\n"
"        col = .9*col;\n"
"    }\n"
"    \n"
"    if(s.x > 1.e-4 && uv.y<-2.e-4)\n"
"    {\n"
"        d = -(o.y+.375)/dir.y;\n"
"        x = o + d * dir;\n"
"        scene(x, s);\n"
"        i=N;\n"
"    }\n"
"    else x = o + d * dir;\n"
"    \n"
"    if(s.y == -1.)\n"
"    {\n"
"        vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));\n"
"    	\n"
"        col = mix(col, c1, .3+.5*inside);\n"
"    }\n"
"    \n"
"    col = mix(col, c.yyy, smoothstep(2., 3., d));\n"
"    \n"
"    col *= 1.2;\n"
"    col *= col;\n"
"    \n"
"    col *= col*col;\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *graffiti_source = "/* Lightcyber by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
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
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);\n"
"void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void rot3(in vec3 phi, out mat3 R); \n"
"\n"
"void graf(in vec2 x, out float d)\n"
"{\n"
"    x.y *= .7;\n"
"    float size = .4,\n"
"        n,\n"
"        da;\n"
"    vec2 y = vec2(mod(x.x, size)-.5*size, x.y),\n"
"        yi = (x-y)/size,\n"
"        x1,\n"
"        x2,\n"
"        x3;\n"
"    \n"
"    dbox(y,vec2(.75,.75)*size, d);\n"
"    \n"
"    // lines\n"
"    rand(yi, n);\n"
"    x1 = vec2(-.5+.02+n*.96, .75)*size,\n"
"    x2 = vec2(.5-.02-n*.96, -.75)*size;\n"
"    x1.x = floor(5.*x1.x)/5.;\n"
"    x2.x = floor(5.*x2.x)/5.;\n"
"    x1.x = max(x1.x,-.4*size);\n"
"    x1.x = min(x1.x,.4*size);\n"
"    x2.x = max(x2.x,-.4*size);\n"
"    x2.x = min(x2.x,.4*size);\n"
"    dlinesegment(y, x1, x2, da);\n"
"    stroke(da, .02, da);\n"
"    d = max(d, -da);\n"
"    \n"
"    // upper triangles\n"
"    rand(yi+1337., n);\n"
"	x1 = vec2(-.55+n,.75)*size*1.05,\n"
"    x2 = vec2(.5-n,.75)*size*1.05,\n"
"    x3 = .75*(.8-n)*size*1.05*c.yx-.01*c.yx;\n"
"    x1 = round(15.*x1)/15.;\n"
"    x2 = round(15.*x2)/15.;\n"
"    x3 = round(15.*x3)/15.;\n"
"    dtriangle(y, x1, x2, x3, da);\n"
"    d = max(d, -da);\n"
"    \n"
"    // lower triangles\n"
"    rand(yi+2337., n);\n"
"    x1 = vec2(-.5+n,-.75)*size*1.05,\n"
"    x2 = vec2(.55-n,-.75)*size*1.05,\n"
"    x3 = -.75*(.8-n)*size*1.05*c.yx+.01*c.yx;\n"
"    x1 = round(15.*x1)/15.;\n"
"    x2 = round(15.*x2)/15.;\n"
"    x3 = round(15.*x3)/15.;\n"
"    dtriangle(y, x1, x2, x3, da);\n"
"    d = max(d, -da);\n"
"}\n"
"\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.x += .3*iTime;\n"
"    x *= 2.;\n"
"    \n"
"    vec3 n;\n"
"    lfnoise(x.x*c.xx-iTime, n.x);\n"
"    lfnoise(2.*x.x*c.xx-iTime-1337., n.y);\n"
"    lfnoise(x.x*c.xx+2.*iTime-2337., n.z);\n"
"\n"
"    x.yz += .1*vec2(cos(x.x), sin(x.x))*n.xy;\n"
"    \n"
"    mat3 RR;\n"
"    rot3(1.3*mix(.2,1.5, .5+.5*n.x)*n.z * c.xyy, RR);\n"
"    x = RR * x;\n"
"    \n"
"    float d, da, db;\n"
"    \n"
"    graf(x.xy, d);\n"
"    stroke(d+mix(.01,.04, iScale), mix(.01,.04, iScale), da);\n"
"    //stroke(d,.01,da);\n"
"    \n"
"    float v;\n"
"    vec2 ind;\n"
"    dvoronoi(12.*x.xy, v, ind); // 12.\n"
"    \n"
"    zextrude(x.z, -d, .1-.1*v, d);\n"
"    \n"
"	sdf = vec2(d,1.);\n"
"    float modsize = .025,\n"
"		y = mod(d-.3-.02*iTime,modsize)-.5*modsize,\n"
"        yi = (d-y)/modsize;\n"
"    \n"
"    float na;\n"
"    lfnoise(2.*yi*c.xx-.3*iTime, na);\n"
"\n"
"    zextrude(x.z-.05*na, -y, mix(0.,.05+.05*na,iScale), d);\n"
"    stroke(d,.035,d);\n"
"        \n"
"    \n"
"    \n"
"    zextrude(x.z, -da, .25, da);\n"
"    \n"
"	add(sdf, vec2(da, 1.), sdf);\n"
"	\n"
"	lfnoise(5.*x.xy, da);\n"
"	    mfnoise(x.xy, 32., 422., .45, db);\n"
"        da = .5*(db+da);\n"
"		sdf.x -= .001*da; // .001\n"
"        stroke(da, .1, da);\n"
"        sdf.x -=.005*da; // .005\n"
"    add(sdf, vec2(d, 1.), sdf);\n"
"    add(sdf, vec2(x.z+.25,1.), sdf);\n"
"    \n"
"    // xa(t)\n"
"    float xa = mix(x.x+3.*a,x.x-3.*a,clamp(iTime/3.,0.,1.));\n"
"    xa = mix(xa,-xa, clamp((iTime-3.)/6., 0.,1.));\n"
"    sdf.x += mix(0., 2., smoothstep(-.5, .5, xa));\n"
"    \n"
"}\n"
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
"    x.x += .3*iTime;\n"
"    x *= 2.;\n"
"\n"
"    float n;\n"
"    lfnoise(x.x*c.xx-iTime, n);\n"
"    x.y += .3*cos(x.x)*n;\n"
" \n"
"    float d;\n"
"    graf(x, d);\n"
"    col = mix(col, mix(vec3(0.85,0.87,0.89), vec3(0.04,0.18,0.24), clamp(abs(x.y/2.),0.,1.)), sm(d-.2));\n"
"    col = mix(col, vec3(1.00,0.40,0.39), sm(d));\n"
"    float da = d;\n"
"    stroke(d+mix(.01,.03, iScale), mix(.01,.04,iScale), d);\n"
"    //stroke(d,.01, d);\n"
"    col = mix(col, 1.4*col, sm(d));\n"
"    stroke(d, .001, d);\n"
"    col = mix(col, 1.3*col, sm(d));\n"
"    \n"
"    if(da < .02 && da > -.02)\n"
"    {\n"
"        lfnoise(5.*x, da);\n"
"	    mfnoise(x, 32., 422., .45, d);\n"
"        d = .5*(d+da);\n"
"		col = mix(col, vec3(0.27,0.27,0.27), sm(d));\n"
"        stroke(d, .1, d);\n"
"        col = mix(col, 1.5*col, sm(d));\n"
"    }\n"
"    \n"
"    col *= mix(1., 1.6, iScale);\n"
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
"    int N = 250,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    vec3 c1;\n"
"    float d = -(o.z-.35)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.15)\n"
"        {\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,5.e-2);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 1.e-2);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+.5*c.yzx);\n"
"            colorize(x.xy, c1);\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            lfnoise(x.xy, r);\n"
"            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),sin(2.*iScale*r*x));\n"
"            c1 = .1*c1\n"
"                + .8*c1 * abs(dot(l,n))\n"
"                + 6.5*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"        col = c1;\n"
"    }\n"
"    \n"
"    \n"
"    //col += col;\n"
"    \n"
"    col *= col*col;\n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"    // Background geraffel\n"
"    if(length(col) < .001)\n"
"    {\n"
"        float v, ra, v2;\n"
"        vec2 ind, ind2;\n"
"        lfnoise(iTime*c.xx, ind2.x);\n"
"        lfnoise(iTime*c.xx-1337., ind2.y);\n"
"        dvoronoi(12.*(uv-.03*ind2), v, ind);\n"
"//         dvoronoi(43.*uv, v2, ind2);\n"
"//         ind += .1*ind2;\n"
"        rand(ind, ra);\n"
"        stroke(-v, .05, v);\n"
"        v = -v;\n"
"        col = mix(col, .3*ra*mix( .5*vec3(1.00,0.40,0.39), .05*c.xxx, clamp(tanh(1.5*length(uv)),0.,1.)), sm(v));\n"
"    }\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *greet_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" *\n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"float iScale;\n"
"\n"
"// void dpolygon(in vec2 x, in float N, out float d);\n"
"// void rot(in float phi, out mat2 m);\n"
"// void dcircle(in vec2 x, out float d);\n"
"// void dbox(in vec2 x, in vec2 b, out float d);\n"
"// void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"// void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void hash31(in float p, out vec3 d);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void dschnappsgirls(in vec2 x, out float d);\n"
"void dspacepigs(in vec2 x, out float d);\n"
"void dkewlers(in vec2 x, out float d);\n"
"void dfarbrausch(in vec2 x, out float d);\n"
"void dhaujobb(in vec2 x, out float d);\n"
"void dmercury(in vec2 x, out float d);\n"
"\n"
"vec2 ind, indc;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.z -= .05*iTime;\n"
"    \n"
"    float d;\n"
"    \n"
"    // Corridor\n"
"    dbox3(x, vec3(.1,.1,1.e3), d);\n"
"    sdf = vec2(-d, 2.);\n"
"    \n"
"    // Wall tiles\n"
"    float distortion;\n"
"	lfnoise(5.2e2*x.yz, distortion);\n"
"    float tsize = .005,\n"
"    	dy = mod(x.y, tsize)-.5*tsize,\n"
"        yi = (x.y-dy)/tsize,\n"
"        zpar = x.z+mix(0., .5*tsize, mod(yi,2.)),\n"
"        dz = mod(zpar, tsize)-.5*tsize,\n"
"        zi = (zpar-dz)/tsize;\n"
"    dbox3(vec3(abs(x.x)-.1, dy, dz), vec3(.0005+.00001*distortion, .39*tsize*c.xx), d);\n"
"    add(sdf, vec2(d, 3.), sdf);\n"
"    smoothmin(sdf.x, d, .001, sdf.x);\n"
"    ind = vec2(yi, zi);\n"
"    \n"
"    // Ceiling \n"
"    tsize = .025;\n"
"    dz = mod(x.z, tsize)-.5*tsize;\n"
"    float dx = mod(x.x, tsize)-.5*tsize;\n"
"    zi = (x.z-dz)/tsize;\n"
"    float xi = (x.x-dx)/tsize;\n"
"    dbox3(vec3(dx, abs(x.y)-.1, dz), vec3(.48*tsize, .0005, .48*tsize), d);\n"
"    add(sdf, vec2(d, 4.), sdf);\n"
"    smoothmin(sdf.x, d, .002, sdf.x);\n"
"    indc = vec2(xi, zi);\n"
"\n"
"    // Logos\n"
"    tsize = .25;\n"
"    float tw = .0005;\n"
"    dz = mod(x.z-.5*tsize, tsize)-.5*tsize;\n"
"    zi = round((x.z-dz)/tsize);\n"
"    zi = mod(zi, 6.);\n"
"    if(zi < .5)\n"
"    {\n"
"        dmercury(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"    else if(zi < 1.5)\n"
"    {\n"
"        dhaujobb(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"    else if(zi < 2.5)\n"
"    {\n"
"        dfarbrausch(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"    else if(zi < 3.5)\n"
"    {\n"
"        dkewlers(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"	else if(zi < 4.5)\n"
"    {\n"
"        dspacepigs(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"	else if(zi < 5.5)\n"
"    {\n"
"        dschnappsgirls(20.*x.xy, d);\n"
"        stroke(d/20.,tw, d);\n"
"        zextrude(dz, -d, .005, d);\n"
"        add(sdf, vec2(d, 5.), sdf);\n"
"    }\n"
"\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    iScale = 1.;\n"
"    \n"
"    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, \n"
"        s;\n"
"    \n"
"    //uv *= 2.;\n"
"    \n"
"    vec3 col = c.yyy, \n"
"        o = c.yyx,\n"
"        r = c.xyy, \n"
"        u = c.yxy, \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = 0.;//-(o.z-.1)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        /*\n"
"        if(x.z<-.05)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }*/\n"
"        d += min(s.x,3.e-2);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        vec3 l = mix(1.5,.2,abs(pow(sin(2.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);\n"
"       \n"
"		if(s.y == 2.)\n"
"        {\n"
"            col = .23*c.xxx;\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"        }\n"
"        else if(s.y == 3.)\n"
"        {\n"
"            float r;\n"
"            rand(ind, r);\n"
"            \n"
"            col = .2*vec3(0.02,0.11,0.24)\n"
"                + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))\n"
"                + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            \n"
"            // Reflections\n"
"            d = 1.e-2;\n"
"            o = x;\n"
"            dir = reflect(dir, n);\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                /*\n"
"                if(x.z<-.05)\n"
"                {\n"
"                    col = .2*c.xxx;\n"
"                    i = N;\n"
"                    break;\n"
"                }*/\n"
"                d += min(s.x,3.e-1);\n"
"                //d += s.x;\n"
"            }\n"
"            \n"
"            normal(x,n, 5.e-4);\n"
"        	vec3 l = mix(1.5,.2,abs(pow(sin(4.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);\n"
"       \n"
"            vec3 c1;\n"
"            if(s.y == 2.)\n"
"            {\n"
"                c1 = .23*c.xxx;\n"
"                c1 = .2*c1\n"
"                    + .2*c1 * abs(dot(l,n))\n"
"                    + .6*c1 * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            }\n"
"            else if(s.y == 3.)\n"
"            {\n"
"                float r;\n"
"                rand(ind, r);\n"
"\n"
"                c1 = .2*vec3(0.02,0.11,0.24)\n"
"                    + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))\n"
"                    + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            }\n"
"            else if(s.y == 4.)\n"
"            {\n"
"                float r;\n"
"                rand(indc, r);\n"
"\n"
"                c1 = .2*.2*c.xxx\n"
"                    + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))\n"
"                    + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            }\n"
"            else if(s.y == 5.)\n"
"            {\n"
"                c1 = .2*.2*c.xyy\n"
"                    + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))\n"
"                    + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"                c1 = mix(c1, c.xxx, .1);\n"
"            }\n"
"            \n"
"            col = mix(col, c1, .5);\n"
"        }\n"
"        else if(s.y == 4.)\n"
"        {\n"
"            float r;\n"
"            rand(indc, r);\n"
"            \n"
"            col = .2*.2*c.xxx\n"
"                + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))\n"
"                + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"        }\n"
"        else if(s.y == 5.)\n"
"        {\n"
"            col = .2*.2*c.xyy\n"
"                + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))\n"
"                + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            col = mix(col, c.xxx, .1);\n"
"        }\n"
"\n"
"    }\n"
"    \n"
"    col = mix(col, 0.*.23*c.xxx*vec3(0.76,0.20,0.13),smoothstep(1.,5.,d));\n"
"    col *=3.6;\n"
"    //col = mix(col, 2.*col, iScale);\n"
"    \n"
"    //col = atan(col);\n"
"    col *= col;\n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}	\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *evoke_source = "/* Lightcyber by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
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
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void devoke(in vec2 x, out float d)\n"
"{\n"
"    x.x += .225;\n"
"    x *= 1.1;\n"
"    \n"
"    // o\n"
"    d = length(x+.35*c.xy)-.1;\n"
"    stroke(d,.06,d);\n"
"    \n"
"    // I\n"
"    float da;\n"
"    dbox(x+.1*c.xy, vec2(.05, .25), da);\n"
"    d = min(d, da);\n"
"    \n"
"    x = 2.*x - vec2(.4,-.2);\n"
"    // Mercury\n"
"    // Upper part\n"
"    dbox(x-.35*c.yx,vec2(.4,.35), da);\n"
"    d = min(d,da);\n"
"    dbox(x-.7*c.yx, vec2(.2,.2), da);\n"
"    d = max(d,-da);\n"
"    dbox(x-.25*c.yx,vec2(.2,.05),da);\n"
"    d = max(d,-da);\n"
"    \n"
"    // Lower part\n"
"    dbox(x+.1*c.yx,vec2(.1,.2),da);\n"
"    d = min(d,da);\n"
"    dbox(x+.2*c.yx, vec2(.4,.1),da);\n"
"    d = min(d,da);\n"
"    \n"
"    x = .5*(x + vec2(.4,-.2));\n"
"    \n"
"    // E\n"
"    // Right\n"
"    dbox(x-.9*c.xy, vec2(.05, .25), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Top/bot\n"
"    dbox(vec2(x.x-.7, abs(x.y)-.2), vec2(.2, .05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Middle\n"
"    dbox(x-.7*c.xy, vec2(.2, .05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Appendix\n"
"    dbox(vec2(x.x-.95,x.y+.2), vec2(.05,.05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    stroke(d,.001,d);\n"
"}\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
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
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = vec3(0.20,0.01,0.14), \n"
"        o = c.yyx,\n"
"        r = c.xyy, \n"
"        u = c.yxy,\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    float d, i, ra;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    for(i=mix(0., 1.4, smoothstep(0., 1., iTime)); i>=0.; i -= .01)\n"
"    {\n"
"        rand(i*c.xx, ra);\n"
"        \n"
"		d = -(o.z-.2+i)/dir.z;\n"
"        x = o + d * dir;\n"
"        \n"
"        devoke(x.xy, s.x);\n"
"        s.x -= .01*iScale;\n"
"        \n"
"        if(ra < .5)\n"
"        {\n"
"            vec3 c1 = mix(mix(vec3(0.75,0.24,0.31), vec3(1.00,0.87,0.57), smoothstep(1.25,1.4,1.4-i)),vec3(0.20,0.01,0.14),i/1.4);\n"
"	        col = mix(col, c1, sm(s.x));\n"
"            col = mix(col, mix(col,vec3(.7,.45,.3), mix(.02,.1,iScale)), sm(s.x/64.));\n"
"        }\n"
"    }\n"
"    \n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
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
void Loaddspline3()
{
    int dspline3_size = strlen(dspline3_source);
    dspline3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspline3_handle, 1, (GLchar **)&dspline3_source, &dspline3_size);
    glCompileShader(dspline3_handle);
#ifdef DEBUG
    printf("---> dspline3 Shader:\n");
    debug(dspline3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox3()
{
    int dbox3_size = strlen(dbox3_source);
    dbox3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox3_handle, 1, (GLchar **)&dbox3_source, &dbox3_size);
    glCompileShader(dbox3_handle);
#ifdef DEBUG
    printf("---> dbox3 Shader:\n");
    debug(dbox3_handle);
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
void Loaddtriangle()
{
    int dtriangle_size = strlen(dtriangle_source);
    dtriangle_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dtriangle_handle, 1, (GLchar **)&dtriangle_source, &dtriangle_size);
    glCompileShader(dtriangle_handle);
#ifdef DEBUG
    printf("---> dtriangle Shader:\n");
    debug(dtriangle_handle);
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
void Loaddpolygon()
{
    int dpolygon_size = strlen(dpolygon_source);
    dpolygon_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dpolygon_handle, 1, (GLchar **)&dpolygon_source, &dpolygon_size);
    glCompileShader(dpolygon_handle);
#ifdef DEBUG
    printf("---> dpolygon Shader:\n");
    debug(dpolygon_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrot()
{
    int rot_size = strlen(rot_source);
    rot_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rot_handle, 1, (GLchar **)&rot_source, &rot_size);
    glCompileShader(rot_handle);
#ifdef DEBUG
    printf("---> rot Shader:\n");
    debug(rot_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddcircle()
{
    int dcircle_size = strlen(dcircle_source);
    dcircle_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dcircle_handle, 1, (GLchar **)&dcircle_source, &dcircle_size);
    glCompileShader(dcircle_handle);
#ifdef DEBUG
    printf("---> dcircle Shader:\n");
    debug(dcircle_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadhash31()
{
    int hash31_size = strlen(hash31_source);
    hash31_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(hash31_handle, 1, (GLchar **)&hash31_source, &hash31_size);
    glCompileShader(hash31_handle);
#ifdef DEBUG
    printf("---> hash31 Shader:\n");
    debug(hash31_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddschnappsgirls()
{
    int dschnappsgirls_size = strlen(dschnappsgirls_source);
    dschnappsgirls_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dschnappsgirls_handle, 1, (GLchar **)&dschnappsgirls_source, &dschnappsgirls_size);
    glCompileShader(dschnappsgirls_handle);
#ifdef DEBUG
    printf("---> dschnappsgirls Shader:\n");
    debug(dschnappsgirls_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddspacepigs()
{
    int dspacepigs_size = strlen(dspacepigs_source);
    dspacepigs_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspacepigs_handle, 1, (GLchar **)&dspacepigs_source, &dspacepigs_size);
    glCompileShader(dspacepigs_handle);
#ifdef DEBUG
    printf("---> dspacepigs Shader:\n");
    debug(dspacepigs_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddkewlers()
{
    int dkewlers_size = strlen(dkewlers_source);
    dkewlers_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dkewlers_handle, 1, (GLchar **)&dkewlers_source, &dkewlers_size);
    glCompileShader(dkewlers_handle);
#ifdef DEBUG
    printf("---> dkewlers Shader:\n");
    debug(dkewlers_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddfarbrausch()
{
    int dfarbrausch_size = strlen(dfarbrausch_source);
    dfarbrausch_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dfarbrausch_handle, 1, (GLchar **)&dfarbrausch_source, &dfarbrausch_size);
    glCompileShader(dfarbrausch_handle);
#ifdef DEBUG
    printf("---> dfarbrausch Shader:\n");
    debug(dfarbrausch_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddhaujobb()
{
    int dhaujobb_size = strlen(dhaujobb_source);
    dhaujobb_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dhaujobb_handle, 1, (GLchar **)&dhaujobb_source, &dhaujobb_size);
    glCompileShader(dhaujobb_handle);
#ifdef DEBUG
    printf("---> dhaujobb Shader:\n");
    debug(dhaujobb_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddmercury()
{
    int dmercury_size = strlen(dmercury_source);
    dmercury_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dmercury_handle, 1, (GLchar **)&dmercury_source, &dmercury_size);
    glCompileShader(dmercury_handle);
#ifdef DEBUG
    printf("---> dmercury Shader:\n");
    debug(dmercury_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}

void LoadSymbols()
{
    Loadrand();
    updateBar();
    Loadlfnoise();
    updateBar();
    Loaddspline3();
    updateBar();
    Loaddbox3();
    updateBar();
    Loaddlinesegment3();
    updateBar();
    Loadstroke();
    updateBar();
    Loadzextrude();
    updateBar();
    Loadsmoothmin();
    updateBar();
    Loadadd();
    updateBar();
    Loaddvoronoi();
    updateBar();
    Loaddbox();
    updateBar();
    Loadrot3();
    updateBar();
    Loadnormal();
    updateBar();
    Loadmfnoise();
    updateBar();
    Loaddtriangle();
    updateBar();
    Loaddlinesegment();
    updateBar();
    Loaddpolygon();
    updateBar();
    Loadrot();
    updateBar();
    Loaddcircle();
    updateBar();
    Loadhash31();
    updateBar();
    Loaddschnappsgirls();
    updateBar();
    Loaddspacepigs();
    updateBar();
    Loaddkewlers();
    updateBar();
    Loaddfarbrausch();
    updateBar();
    Loaddhaujobb();
    updateBar();
    Loaddmercury();
    updateBar();
}
int voronoidesign_program, voronoidesign_handle, groundboxes_program, groundboxes_handle, graffiti_program, graffiti_handle, greet_program, greet_handle, evoke_program, evoke_handle;
int voronoidesign_iTime_location;
voronoidesign_iResolution_location;
voronoidesign_iFader0_location;
voronoidesign_iFader1_location;
voronoidesign_iFader2_location;
voronoidesign_iFader3_location;
voronoidesign_iFader4_location;
voronoidesign_iFader5_location;
voronoidesign_iFader6_location;
voronoidesign_iFader7_location;
int groundboxes_iTime_location;
groundboxes_iResolution_location;
groundboxes_iFader0_location;
groundboxes_iFader1_location;
groundboxes_iFader2_location;
groundboxes_iFader3_location;
groundboxes_iFader4_location;
groundboxes_iFader5_location;
groundboxes_iFader6_location;
groundboxes_iFader7_location;
int graffiti_iTime_location;
graffiti_iResolution_location;
graffiti_iFader0_location;
graffiti_iFader1_location;
graffiti_iFader2_location;
graffiti_iFader3_location;
graffiti_iFader4_location;
graffiti_iFader5_location;
graffiti_iFader6_location;
graffiti_iFader7_location;
int greet_iTime_location;
greet_iResolution_location;
greet_iFader0_location;
greet_iFader1_location;
greet_iFader2_location;
greet_iFader3_location;
greet_iFader4_location;
greet_iFader5_location;
greet_iFader6_location;
greet_iFader7_location;
int evoke_iTime_location;
evoke_iResolution_location;
evoke_iFader0_location;
evoke_iFader1_location;
evoke_iFader2_location;
evoke_iFader3_location;
evoke_iFader4_location;
evoke_iFader5_location;
evoke_iFader6_location;
evoke_iFader7_location;
const int nprograms = 5;

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
    glLinkProgram(voronoidesign_program);
#ifdef DEBUG
    printf("---> voronoidesign Program:\n");
    debugp(voronoidesign_program);
    printf(">>>>\n");
#endif
    glUseProgram(voronoidesign_program);
    voronoidesign_iTime_location = glGetUniformLocation(voronoidesign_program, "iTime");
    voronoidesign_iResolution_location = glGetUniformLocation(voronoidesign_program, "iResolution");
    voronoidesign_iFader0_location = glGetUniformLocation(voronoidesign_program, "iFader0");
    voronoidesign_iFader1_location = glGetUniformLocation(voronoidesign_program, "iFader1");
    voronoidesign_iFader2_location = glGetUniformLocation(voronoidesign_program, "iFader2");
    voronoidesign_iFader3_location = glGetUniformLocation(voronoidesign_program, "iFader3");
    voronoidesign_iFader4_location = glGetUniformLocation(voronoidesign_program, "iFader4");
    voronoidesign_iFader5_location = glGetUniformLocation(voronoidesign_program, "iFader5");
    voronoidesign_iFader6_location = glGetUniformLocation(voronoidesign_program, "iFader6");
    voronoidesign_iFader7_location = glGetUniformLocation(voronoidesign_program, "iFader7");
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
    glAttachShader(groundboxes_program,rand_handle);
    glAttachShader(groundboxes_program,lfnoise_handle);
    glAttachShader(groundboxes_program,dspline3_handle);
    glAttachShader(groundboxes_program,dbox3_handle);
    glAttachShader(groundboxes_program,dlinesegment3_handle);
    glAttachShader(groundboxes_program,stroke_handle);
    glAttachShader(groundboxes_program,zextrude_handle);
    glAttachShader(groundboxes_program,smoothmin_handle);
    glAttachShader(groundboxes_program,add_handle);
    glAttachShader(groundboxes_program,dvoronoi_handle);
    glAttachShader(groundboxes_program,dbox_handle);
    glAttachShader(groundboxes_program,rot3_handle);
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
    groundboxes_iFader0_location = glGetUniformLocation(groundboxes_program, "iFader0");
    groundboxes_iFader1_location = glGetUniformLocation(groundboxes_program, "iFader1");
    groundboxes_iFader2_location = glGetUniformLocation(groundboxes_program, "iFader2");
    groundboxes_iFader3_location = glGetUniformLocation(groundboxes_program, "iFader3");
    groundboxes_iFader4_location = glGetUniformLocation(groundboxes_program, "iFader4");
    groundboxes_iFader5_location = glGetUniformLocation(groundboxes_program, "iFader5");
    groundboxes_iFader6_location = glGetUniformLocation(groundboxes_program, "iFader6");
    groundboxes_iFader7_location = glGetUniformLocation(groundboxes_program, "iFader7");
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
    glAttachShader(graffiti_program,rand_handle);
    glAttachShader(graffiti_program,lfnoise_handle);
    glAttachShader(graffiti_program,mfnoise_handle);
    glAttachShader(graffiti_program,dtriangle_handle);
    glAttachShader(graffiti_program,dbox_handle);
    glAttachShader(graffiti_program,dlinesegment_handle);
    glAttachShader(graffiti_program,stroke_handle);
    glAttachShader(graffiti_program,dvoronoi_handle);
    glAttachShader(graffiti_program,rot3_handle);
    glAttachShader(graffiti_program,zextrude_handle);
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
    graffiti_iFader0_location = glGetUniformLocation(graffiti_program, "iFader0");
    graffiti_iFader1_location = glGetUniformLocation(graffiti_program, "iFader1");
    graffiti_iFader2_location = glGetUniformLocation(graffiti_program, "iFader2");
    graffiti_iFader3_location = glGetUniformLocation(graffiti_program, "iFader3");
    graffiti_iFader4_location = glGetUniformLocation(graffiti_program, "iFader4");
    graffiti_iFader5_location = glGetUniformLocation(graffiti_program, "iFader5");
    graffiti_iFader6_location = glGetUniformLocation(graffiti_program, "iFader6");
    graffiti_iFader7_location = glGetUniformLocation(graffiti_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadgreet()
{
    int greet_size = strlen(greet_source);
    greet_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(greet_handle, 1, (GLchar **)&greet_source, &greet_size);
    glCompileShader(greet_handle);
#ifdef DEBUG
    printf("---> greet Shader:\n");
    debug(greet_handle);
    printf(">>>>\n");
#endif
    greet_program = glCreateProgram();
    glAttachShader(greet_program,greet_handle);
    glAttachShader(greet_program,dpolygon_handle);
    glAttachShader(greet_program,rot_handle);
    glAttachShader(greet_program,dcircle_handle);
    glAttachShader(greet_program,dbox_handle);
    glAttachShader(greet_program,dlinesegment_handle);
    glAttachShader(greet_program,dtriangle_handle);
    glAttachShader(greet_program,rand_handle);
    glAttachShader(greet_program,hash31_handle);
    glAttachShader(greet_program,lfnoise_handle);
    glAttachShader(greet_program,dbox_handle);
    glAttachShader(greet_program,stroke_handle);
    glAttachShader(greet_program,zextrude_handle);
    glAttachShader(greet_program,add_handle);
    glAttachShader(greet_program,smoothmin_handle);
    glAttachShader(greet_program,dbox3_handle);
    glAttachShader(greet_program,dschnappsgirls_handle);
    glAttachShader(greet_program,dspacepigs_handle);
    glAttachShader(greet_program,dkewlers_handle);
    glAttachShader(greet_program,dfarbrausch_handle);
    glAttachShader(greet_program,dhaujobb_handle);
    glAttachShader(greet_program,dmercury_handle);
    glAttachShader(greet_program,normal_handle);
    glLinkProgram(greet_program);
#ifdef DEBUG
    printf("---> greet Program:\n");
    debugp(greet_program);
    printf(">>>>\n");
#endif
    glUseProgram(greet_program);
    greet_iTime_location = glGetUniformLocation(greet_program, "iTime");
    greet_iResolution_location = glGetUniformLocation(greet_program, "iResolution");
    greet_iFader0_location = glGetUniformLocation(greet_program, "iFader0");
    greet_iFader1_location = glGetUniformLocation(greet_program, "iFader1");
    greet_iFader2_location = glGetUniformLocation(greet_program, "iFader2");
    greet_iFader3_location = glGetUniformLocation(greet_program, "iFader3");
    greet_iFader4_location = glGetUniformLocation(greet_program, "iFader4");
    greet_iFader5_location = glGetUniformLocation(greet_program, "iFader5");
    greet_iFader6_location = glGetUniformLocation(greet_program, "iFader6");
    greet_iFader7_location = glGetUniformLocation(greet_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadevoke()
{
    int evoke_size = strlen(evoke_source);
    evoke_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(evoke_handle, 1, (GLchar **)&evoke_source, &evoke_size);
    glCompileShader(evoke_handle);
#ifdef DEBUG
    printf("---> evoke Shader:\n");
    debug(evoke_handle);
    printf(">>>>\n");
#endif
    evoke_program = glCreateProgram();
    glAttachShader(evoke_program,evoke_handle);
    glAttachShader(evoke_program,rand_handle);
    glAttachShader(evoke_program,dbox_handle);
    glAttachShader(evoke_program,stroke_handle);
    glLinkProgram(evoke_program);
#ifdef DEBUG
    printf("---> evoke Program:\n");
    debugp(evoke_program);
    printf(">>>>\n");
#endif
    glUseProgram(evoke_program);
    evoke_iTime_location = glGetUniformLocation(evoke_program, "iTime");
    evoke_iResolution_location = glGetUniformLocation(evoke_program, "iResolution");
    evoke_iFader0_location = glGetUniformLocation(evoke_program, "iFader0");
    evoke_iFader1_location = glGetUniformLocation(evoke_program, "iFader1");
    evoke_iFader2_location = glGetUniformLocation(evoke_program, "iFader2");
    evoke_iFader3_location = glGetUniformLocation(evoke_program, "iFader3");
    evoke_iFader4_location = glGetUniformLocation(evoke_program, "iFader4");
    evoke_iFader5_location = glGetUniformLocation(evoke_program, "iFader5");
    evoke_iFader6_location = glGetUniformLocation(evoke_program, "iFader6");
    evoke_iFader7_location = glGetUniformLocation(evoke_program, "iFader7");
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
    Loadgreet();
    updateBar();
    Loadevoke();
    updateBar();
}
#endif
