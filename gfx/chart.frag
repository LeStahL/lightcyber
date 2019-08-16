#version 130 
 
uniform float iTime;
uniform vec2 iResolution;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);
float a = 1.0;

float nbeats, iScale;

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

void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)
{
    n = 0.;
    float a = 1., nf = 0., buf;
    for(float f = d; f<b; f *= 2.)
    {
        lfnoise(f*x, buf);
        n += a*buf;
        a *= e;
        nf += 1.;
    }
    n *= (1.-e)/(1.-pow(e, nf));
}

void dbox3(in vec3 x, in vec3 b, out float d)
{
  vec3 da = abs(x) - b;
  d = length(max(da,0.0))
         + min(max(da.x,max(da.y,da.z)),0.0);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Add sdfs
void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = sda.x<sdb.x?sda:sdb;
}

void rotAB(in vec3 a, in vec3 b, out mat3 R)
{
    a = normalize(a);
    b = normalize(b);
    vec3 v = cross(a,b);
    float co = dot(a,b);
    R = mat3(0.,v.z,-v.y,-v.z,0.,v.x,v.y,-v.x,0.);
    R = R*R/(1.+co) + R;
    R += mat3(1.);
}

mat3 R;
float ind;
void scene(in vec3 x, out vec2 sdf)
{
    sdf = 3.*c.yx;
    
    float d, da;
    dbox3(x, 2.*c.xxx, sdf.x);
    sdf.x *= -1.;
	sdf.x = min(sdf.x, x.z);
    
    vec3 y = x;
    x = vec3(x.x, mod(x.y,.2)-.1, x.z);
    vec3 xi = (y-x)/.2;
    ind = xi.y;
    
    float dx;
    lfnoise(ind*c.xx, dx);
    
    float n;    
    mfnoise((x.x+.3*iTime)*c.xx-xi.y,12.,120.,.45, n);
    n *= clamp(.2-.01*xi.y,0.,1.);
    
    vec2 sda;
    dbox3(x-1.4*dx*c.xyy, vec3(.5,.01,.3*(.5+.5*dx)-n), sda.x);
    sda.y = 1.;
    add(sdf,sda,sdf);
    
    dbox3(x-(.3*(.5+.5*dx)-n)*c.yyx-1.4*dx*c.xyy, vec3(.5,.03,.01), sda.x);
    sda.y = 4.;
    add(sdf, sda, sdf);
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

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void palette1(in float scale, out vec3 col)
{
    const int N = 8;

    const vec3 colors[N] = vec3[N](
       	vec3(0.18,0.30,0.65),
        vec3(0.00,0.69,0.80),
        vec3(0.45,0.72,0.50),
        vec3(0.93,0.36,0.44),
        vec3(1.00,0.51,0.45),
        vec3(0.78,0.57,0.78),
        vec3(0.54,0.33,0.60),
        vec3(0.98,0.78,0.38)
    );

	float index = floor(scale*float(N)), 
        remainder = scale*float(N)-index;
    col = mix(colors[int(index)],colors[int(index)+1], remainder);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);

    vec3 col = c.yyy, 
        o = c.yzx,
        r = c.xyy, 
        u = normalize(c.yxx), 
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 400,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);
    
    float d = 0.;
    vec2 s;
    
    for(i = 0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        //d += s.x;
        d += min(s.x, 7.e-3);
    }
    
    if(s.x < 1.e-4)
    {
        normal(x,n, 5.e-4);
        vec3 l = normalize(x+.1*n);
        
        if(s.y == 1.)
        {
            ind = .5+.5*clamp(ind/10.,-1.,1.);
            palette1(ind, col);
            vec3 c1;
            palette1(ind+.1, c1);
            col = mix(col,c1,clamp(1.-x.z/.3,0.,1.));
            col = .3*col
                + .4*col * abs(dot(l,n))
                + .9 * col * abs(pow(dot(reflect(-l,n),dir),3.));
        }
        else if(s.y == 4.)
        {
        	col = c.xxx;
            col = .5*col
            	+ .4*col * abs(dot(l,n))
                + .8 * col * abs(pow(dot(reflect(-l,n),dir),3.));
        }
        else if(s.y == 2. || s.y == 3.)
        {
            col = c.xxx;
            vec2 ma = abs(mod(x.xy,.05)-.025)-.0015;
            col = mix(col,.5*c.xxx, sm(min(ma.x,ma.y)));
            col = .5*col
                + .4*col * abs(dot(l,n))
                + .8 * col * abs(pow(dot(reflect(-l,n),dir),3.));
            
            vec3 c1 = c.yyy;
            
            o = x;
            dir = reflect(dir,n);
            d = 1.e-2;
            
            N = 100;
            
            for(i = 0; i<N; ++i)
            {
                x = o + d * dir;
                scene(x,s);
                if(s.x < 1.e-4)break;
                d += s.x;
            }
                
            if(s.x < 1.e-4)
            {
                normal(x,n, 5.e-4);
                vec3 l = normalize(x+.1*n);

                if(s.y == 1.)
                {
                    ind = .5+.5*clamp(ind/10.,-1.,1.);
                    palette1(ind, c1);
                    vec3 c2;
                    palette1(ind+.1, c2);
                    c1 = mix(c1,c1,clamp(1.-x.z/.3,0.,1.));
                    c1 = .3*c1
                        + .4*c1 * abs(dot(l,n))
                        + .9 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));
                }
                else if(s.y == 2.)
                {
                    c1 = .7*c.xxx;
                    c1 = .5*c1
                        + .4*c1 * abs(dot(l,n))
                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));
                    
                }
                else if(s.y == 3.)
                {
                    c1 = .7*c.xxx;
                    c1 = .5*c1
                        + .4*c1 * abs(dot(l,n))
                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));
                }
                else if(s.y == 4.)
                {
                    c1 = c.xxx;
                    c1 = .5*c1
                        + .4*c1 * abs(dot(l,n))
                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));
                }
                
                c1 = clamp(c1, 0., 1.);
                
                col = mix(col, c1, .2);
            }
            
            col = clamp(col, 0., 1.);
        }
     	
    }
    //col = mix(col,vec3(0.20,0.01,0.14),smoothstep(0.,1.,iTime-10.));
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
